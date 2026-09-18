#!/usr/bin/env python3
"""Diagnose projected noise, angular/radial transfer, and interpolation on the saved baseline.

The previously fixed 16 common split sites and all 24 covariance policies are
development diagnostics. Train/validation stencils read disjoint native half-planes;
normalized fits condition on the existing shared radial profile. No new reductions,
detection thresholds, production policies, or uncertainty corrections are introduced.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.linalg import cho_solve, eigh

import compare_p4_step5_radial_pooling as radial
from compare_p4_step5_variance_floor import FLOORS, compare_saved, policy_name
from run_p4_step5_full_injections import fingerprint, write_json


def score_stats(scores: np.ndarray) -> dict:
    """Separate centered dispersion and mean offset in standardized projected residuals."""
    if len(scores) < 2:
        return {'valid': False, 'patches': len(scores)}
    mean = float(np.mean(scores))
    variance = float(np.var(scores, ddof=1))
    mse = float(np.mean(scores * scores))
    radial.require(np.isclose(mse, (len(scores) - 1) / len(scores) * variance + mean * mean,
                              rtol=1e-12, atol=1e-12), 'projected MSE decomposition failed')
    return {'valid': True, 'patches': len(scores), 'mean_over_sigma': mean,
            'variance_over_predicted': variance, 'mse_over_predicted': mse,
            'mean_squared_fraction_of_mse': mean * mean / mse if mse > 0 else 0,
            'within_one_sigma': int(np.sum(np.abs(scores) <= 1)),
            'within_one_sigma_after_validation_demeaning': int(np.sum(np.abs(scores - mean) <= 1)),
            'scores': scores.tolist()}


def nearest_samples(image: np.ndarray, ring: dict) -> np.ndarray:
    """Use the nearest native neighbor at each accepted bilinear output coordinate."""
    chosen = np.argmax(ring['weights'], axis=2)[..., None]
    indices = np.take_along_axis(ring['indices'], chosen, axis=2)[..., 0]
    return image.ravel()[indices]


def native_support(ring: dict, selected: np.ndarray) -> np.ndarray:
    """Return exact native pixels read by the selected patch stencils."""
    return np.unique(ring['indices'][selected][ring['weights'][selected] != 0])


def white_noise_gains(ring: dict, selected: np.ndarray, weight: np.ndarray) -> dict:
    """Calculate exact projected variance gains for unit independent native pixels."""
    bilinear, nearest = [], []
    for indices, coefficients in zip(ring['indices'][selected], ring['weights'][selected]):
        unique, inverse = np.unique(indices.ravel(), return_inverse=True)
        projected = np.bincount(inverse, weights=(coefficients * weight[:, None]).ravel(), minlength=len(unique))
        selected_indices = indices[np.arange(len(weight)), np.argmax(coefficients, axis=1)]
        _, nearest_inverse = np.unique(selected_indices, return_inverse=True)
        projected_nearest = np.bincount(nearest_inverse, weights=weight)
        native_variance = float(weight @ weight)
        bilinear.append(float(projected @ projected / native_variance))
        nearest.append(float(projected_nearest @ projected_nearest / native_variance))
    return {'bilinear_over_native': bilinear, 'nearest_over_native': nearest}


def variance_budget(training: np.ndarray, model: dict, template: np.ndarray,
                    weight: np.ndarray, sigma: float) -> dict:
    """Resolve the training projection variance into retained and discarded eigenspaces."""
    centered = training - model['mean']
    covariance = centered.T @ centered / (len(training) - 1)
    eigenvalues, eigenvectors = eigh(covariance, subset_by_index=(118, 120), check_finite=False)
    retained = eigenvalues > model['floor']
    modes = eigenvectors[:, retained]
    coefficients = modes.T @ weight
    complement = weight - modes @ coefficients
    empirical_total = float(weight @ covariance @ weight)
    empirical_retained = float(np.sum(eigenvalues[retained] * coefficients**2))
    empirical_complement = empirical_total - empirical_retained
    modeled_complement = float(model['floor'] * (complement @ complement))
    radial.require(empirical_complement >= -1e-12 and modeled_complement > 0, 'invalid variance partition')
    radial.require(np.isclose(empirical_retained + modeled_complement, sigma*sigma, rtol=1e-10), 'modeled variance partition mismatch')
    outside_template = template - modes @ (modes.T @ template)
    return {'retained_modes': int(retained.sum()),
        'retained_empirical_over_predicted': empirical_retained / (sigma*sigma),
        'complement_empirical_over_predicted': empirical_complement / (sigma*sigma),
        'complement_modeled_over_predicted': modeled_complement / (sigma*sigma),
        'complement_empirical_over_modeled': empirical_complement / modeled_complement,
        'complement_fraction_of_empirical_variance': empirical_complement / empirical_total,
        'template_energy_in_complement_fraction': float(outside_template @ outside_template / (template @ template))}


def summarize(records: list, nulls: list, common: list) -> dict:
    """Summarize paired directional fits without treating reused patches as independent trials."""
    summary = {'common_sites': common, 'directional_fits_per_policy': 2 * len(common), 'policies': {}}
    for fraction in FLOORS:
        for base, _, _ in radial.POLICIES:
            name = policy_name(base, fraction)
            selected = [r for r in records if r['policy'] == name]
            radial.require(len(selected) == 2 * len(common), 'incomplete paired-fit comparison')
            one = {}
            for key in ('training', 'same_radius', 'same_radius_nearest', 'matched_band'):
                one[key] = {'median_' + metric: float(np.median([r[key][metric] for r in selected])) for metric in
                    ('variance_over_predicted', 'mse_over_predicted', 'mean_squared_fraction_of_mse')}
                one[key]['median_absolute_mean_over_sigma'] = float(np.median([abs(r[key]['mean_over_sigma']) for r in selected]))
                one[key]['median_within_one_sigma_fraction'] = float(np.median([r[key]['within_one_sigma'] / r[key]['patches'] for r in selected]))
            one['median_paired_bilinear_over_nearest_variance'] = float(np.median([
                r['same_radius']['variance_over_predicted'] / r['same_radius_nearest']['variance_over_predicted'] for r in selected]))
            one['median_per_fit_bilinear_white_noise_gain'] = float(np.median([
                np.median(r['white_noise']['bilinear_over_native']) for r in selected]))
            one['median_per_fit_nearest_white_noise_gain'] = float(np.median([
                np.median(r['white_noise']['nearest_over_native']) for r in selected]))
            one['variance_budget'] = {'median_' + key: float(np.median([r['variance_budget'][key] for r in selected]))
                for key in ('complement_empirical_over_modeled', 'complement_fraction_of_empirical_variance',
                            'complement_modeled_over_predicted', 'template_energy_in_complement_fraction')}
            one['radial_transfer'] = {}
            for offset in range(-20, 21, 5):
                usable = [r for r in selected if r['rings'].get(str(offset), {}).get('patches', 0) >= 8]
                one['radial_transfer'][str(offset)] = {'directional_fits': len(usable)}
                if usable:
                    one['radial_transfer'][str(offset)].update({
                        'median_variance_over_predicted': float(np.median([r['rings'][str(offset)]['variance_over_predicted'] for r in usable])),
                        'median_paired_ring_over_same_radius_variance': float(np.median([
                            r['rings'][str(offset)]['variance_over_predicted'] / r['same_radius']['variance_over_predicted'] for r in usable]))})
            native = [r['policies'][name]['center']['score'] for r in nulls
                      if r['role'] == 'evaluation' and r['original_common_eligible']]
            one['original_native_evaluation_centers'] = score_stats(np.array(native))
            summary['policies'][name] = one
    return summary


def plots(summary: dict, output: Path) -> None:
    """Plot the fixed-sample projection diagnostic and its interpolation control."""
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout='constrained')
    panels = [
        ('Training projected variance / model', lambda r: r['training']['median_variance_over_predicted'], True),
        ('Held-out same-radius variance / model', lambda r: r['same_radius']['median_variance_over_predicted'], True),
        ('Held-out mean² / mean squared error', lambda r: r['same_radius']['median_mean_squared_fraction_of_mse'], False),
        ('Held-out nearest-neighbor variance / model', lambda r: r['same_radius_nearest']['median_variance_over_predicted'], True),
        ('Paired bilinear / nearest variance', lambda r: r['median_paired_bilinear_over_nearest_variance'], False),
        ('Bilinear white-noise gain vs native grid', lambda r: r['median_per_fit_bilinear_white_noise_gain'], False)]
    for axis, (title, getter, logarithmic) in zip(axes.flat, panels):
        values = np.array([[getter(summary['policies'][policy_name(base, f)]) for f in FLOORS]
                           for base, _, _ in radial.POLICIES])
        rendered = np.log10(values) if logarithmic else values
        lower, upper = (min(-1, rendered.min()), max(2, rendered.max())) if logarithmic else (0, max(1, rendered.max()))
        mesh = axis.imshow(rendered, vmin=lower, vmax=upper, cmap='Blues', aspect='auto')
        axis.set(title=title, xticks=range(3), xticklabels=FLOORS, xlabel='Variance-floor fraction',
                 yticks=range(8), yticklabels=[f'{"Norm" if n else "Raw"} ±{w}' for _, w, n in radial.POLICIES])
        for y in range(8):
            for x in range(3):
                axis.text(x, y, f'{values[y,x]:.2g}', ha='center', va='center',
                          color='white' if rendered[y,x] > lower + 0.6 * (upper - lower) else 'black')
        fig.colorbar(mesh, ax=axis, shrink=0.7, label='log10 ratio' if logarithmic else 'ratio')
    fig.suptitle('Projected-noise development diagnostic: 16 sites × two angular directions\n'
                 'Same held-out ring for every training width; normalized fits share a fixed radial profile')
    fig.savefig(output / 'projection.png', dpi=170)
    plt.close(fig)


def main() -> None:
    """Freeze the diagnostic, verify archived controls, and project untouched opposite-half patches."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--floor-comparison', type=Path, required=True)
    parser.add_argument('--radial-comparison', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    floor_root, radial_root = args.floor_comparison.resolve(), args.radial_comparison.resolve()
    prior = radial.read_json(radial_root / 'manifest.json')
    floor_manifest = radial.read_json(floor_root / 'manifest.json')
    nulls = radial.read_json(floor_root / 'nulls.json')
    common = radial.read_json(floor_root / 'summary.json')['common_split_stability_trials']
    trials = [row for row in nulls if row['name'] in common]
    radial.require(len(trials) == len(common) == 16, 'changed common-site protocol')
    science_path = Path(prior['inputs'][0]['path'])
    field = Path(prior['inputs'][1]['path']).parent
    paths = [science_path, field / 'p4PSF_coordinates.fits', field / 'p4PSF_model_0000.fits', field / 'p4PSF_validity_0000.fits']
    known = {r['path']: r for r in prior['inputs']}
    for path in paths:
        radial.require(fingerprint(path) == known[str(path)], 'changed original FITS input')
    radial.require(fingerprint(radial_root / 'manifest.json') in floor_manifest['inputs'], 'changed original sampling manifest')
    inputs = [fingerprint(p) for p in paths + [radial_root / 'manifest.json', floor_root / 'manifest.json',
              floor_root / 'nulls.json', floor_root / 'summary.json', floor_root / 'variance_profiles.json']]
    scripts = [fingerprint(Path(__file__)), fingerprint(Path(radial.__file__)),
               fingerprint(Path(__file__).with_name('compare_p4_step5_variance_floor.py')),
               fingerprint(Path(__file__).with_name('run_p4_step5_full_injections.py'))]
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'manifest.json', {'purpose': 'development-only projected-noise diagnostic',
        'common_sites': common, 'floor_fractions': FLOORS, 'sampling_policies': radial.POLICIES,
        'training': 'one native y half-plane, straddlers removed; both directional fits',
        'primary_validation': 'opposite half-plane, same-radius ring, at least eight patches; fixed across training widths',
        'radial_transfer': 'each opposite-half ring separately; summary requires at least eight patches',
        'normalization': 'condition on existing baseline radial profile estimated outside all original source/holdout circles; shared by halves',
        'interpolation': 'paired nearest-neighbor at identical accepted coordinates; no additional pixels read; analytic unit-white-noise propagation includes reused native pixels',
        'inputs': inputs, 'scripts': scripts})
    science = fits.getdata(science_path).squeeze().astype(float)
    yy, xx = np.indices(science.shape)
    forbidden = np.zeros(science.shape, bool)
    for x, y, radius in prior['exclusions']:
        forbidden |= np.hypot(xx - x, yy - y) <= radius
    scale_map, profile = radial.variance_profile(science, forbidden)
    controls = compare_saved(profile, radial.read_json(floor_root / 'variance_profiles.json')['baseline'], 'profile')
    coordinates = fits.getdata(paths[1]).T
    responses = fits.getdata(paths[2]).astype(float)
    validity = fits.getdata(paths[3]).ravel()
    templates = {(int(c[0]), int(c[1])): t.ravel() for c, t, v in zip(coordinates, responses, validity)
                 if v == 1 and np.isfinite(t).all()}
    records, geometry_checks, unit_checks = [], 0, 0
    for trial in trials:
        position = (trial['row'], trial['column'])
        rings = radial.geometry(science, position, forbidden)
        x, y = position
        candidate_scale = scale_map[y-5:y+6, x-5:x+6].ravel()
        template = templates[position]
        for offset, ring in rings.items():
            first, second = [native_support(ring, ring['halves'] == half) for half in (0, 1)]
            radial.require(not np.intersect1d(first, second).size, 'angular stencils overlap')
            radial.require(not np.any(forbidden.ravel()[native_support(ring, np.ones(len(ring['halves']), bool))]), 'forbidden pixel used')
            radial.require(np.allclose(ring['weights'].sum(axis=2), 1, rtol=0, atol=1e-14), 'interpolation does not preserve a constant')
            geometry_checks += 1
        # Cross-ring disjointness matters when wider bands reuse neighboring rings.
        support = [np.unique(np.concatenate([native_support(r, r['halves'] == half) for r in rings.values()])) for half in (0, 1)]
        radial.require(not np.intersect1d(*support).size, 'cross-ring angular stencils overlap')
        for fraction in FLOORS:
            full, _ = radial.analyze_position(science, scale_map, position, template, rings, split=True, floor_fraction=fraction)
            for base, _, _ in radial.POLICIES:
                controls += compare_saved(full[base], trial['policies'][policy_name(base, fraction)]['center'], trial['name'])
            for base, width, normalized in radial.POLICIES:
                image = science / scale_map if normalized else science
                matrices = {offset: radial.extract(image, ring) for offset, ring in rings.items()}
                nearest = nearest_samples(image, rings[0])
                t = template / candidate_scale if normalized else template
                selected = [offset for offset in rings if abs(offset) <= width]
                for training_half in (0, 1):
                    validation_half = 1 - training_half
                    training = np.vstack([matrices[o][rings[o]['halves'] == training_half] for o in selected])
                    model = radial.fit(training, fraction)
                    radial.require(model is not None, 'common-site training unexpectedly invalid')
                    precision_template = cho_solve(model['factorization'], t, check_finite=False)
                    energy = float(t @ precision_template)
                    weight, sigma = precision_template / energy, 1 / math.sqrt(energy)
                    radial.require(np.isclose(weight @ t, 1, rtol=1e-12) and
                                   np.isclose(weight @ model['covariance'] @ weight, sigma*sigma, rtol=1e-10), 'filter normalization mismatch')
                    projected = {str(o): score_stats((matrices[o][rings[o]['halves'] == validation_half] - model['mean']) @ weight / sigma) for o in rings}
                    radial.require(projected['0']['patches'] >= 8, 'changed common validation ring')
                    heldout_band = np.vstack([matrices[o][rings[o]['halves'] == validation_half] for o in selected])
                    same_nearest = nearest[rings[0]['halves'] == validation_half]
                    white = white_noise_gains(rings[0], rings[0]['halves'] == validation_half, weight)
                    row = {'site': trial['name'], 'policy': policy_name(base, fraction), 'training_half': training_half,
                        'sigma': sigma, 'training': score_stats((training - model['mean']) @ weight / sigma),
                        'same_radius': projected['0'], 'matched_band': score_stats((heldout_band - model['mean']) @ weight / sigma),
                        'same_radius_nearest': score_stats((same_nearest - model['mean']) @ weight / sigma),
                        'rings': projected, 'white_noise': white,
                        'variance_budget': variance_budget(training, model, t, weight, sigma)}
                    budget = row['variance_budget']
                    radial.require(np.isclose(budget['retained_empirical_over_predicted'] + budget['complement_empirical_over_predicted'],
                                             row['training']['variance_over_predicted'], rtol=1e-10), 'empirical variance partition mismatch')
                    # The direct sample-covariance projection must equal the score variance.
                    centered = matrices[0][rings[0]['halves'] == validation_half]
                    centered = centered - centered.mean(axis=0)
                    covariance = centered.T @ centered / (len(centered) - 1)
                    radial.require(np.isclose(weight @ covariance @ weight / (sigma*sigma),
                                             row['same_radius']['variance_over_predicted'], rtol=1e-10), 'projected covariance mismatch')
                    unit_checks += 1
                    records.append(row)
        print(f'completed {trial["name"]}', flush=True)
    summary = summarize(records, nulls, common)
    write_json(args.output / 'projections.json', records)
    write_json(args.output / 'summary.json', summary)
    for record in [*inputs, *scripts]:
        radial.require(fingerprint(Path(record['path'])) == record, 'input or script changed during diagnostic')
    verification = {'archived_control_scalar_values': controls, 'ring_geometry_checks': geometry_checks,
        'directional_projection_checks': unit_checks, 'train_validation_native_support_disjoint_including_cross_ring': True,
        'original_forbidden_pixels_unused': True, 'filter_unit_response_and_model_variance_agree': True,
        'retained_and_complement_variance_budgets_agree': True,
        'direct_covariance_projection_and_MSE_decomposition_agree': True, 'source_and_script_fingerprints_unchanged': True}
    write_json(args.output / 'verification.json', verification)
    plots(summary, args.output)
    write_json(args.output / 'complete.json', {'baseline_images': 1, 'new_reductions': 0, 'new_positive_images': 0,
        'policies': 24, 'common_sites': 16, 'directional_fits': len(records), 'verification': verification})


if __name__ == '__main__':
    main()
