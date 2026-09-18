#!/usr/bin/env python3
"""Test trace-preserving covariance shrinkage on fixed Step-5 noise projections.

Use the previously declared strengths, sampling grid, 16 common sites, and two
angular directions. Compare conditional calibration, physical amplitude variance,
and weight stability with the verified three-mode floor-1 control. No injection
recovery, new reductions, or production configuration changes are performed.
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
from scipy.linalg import cho_factor, cho_solve, eigvalsh

import compare_p4_step5_radial_pooling as radial
import diagnose_p4_step5_projected_noise as projection
from compare_p4_step5_variance_floor import compare_saved
from run_p4_step5_full_injections import fingerprint, write_json

STRENGTHS = (0.1, 0.3, 1.0)


def fit_shrinkage(samples: np.ndarray, strength: float) -> dict | None:
    """Shrink all empirical covariance modes toward their mean variance, preserving the trace."""
    radial.require(np.isfinite(strength) and 0 < strength <= 1, 'shrinkage strength must be in (0, 1]')
    if len(samples) < 8:
        return None
    radial.require(np.isfinite(samples).all(), 'nonfinite covariance training sample')
    mean = samples.mean(axis=0)
    centered = samples - mean
    covariance = centered.T @ centered / (len(samples) - 1)
    target = float(np.trace(covariance) / covariance.shape[0])
    if not np.isfinite(target) or target <= 0:
        return None
    model = (1 - strength) * covariance + strength * target * np.eye(covariance.shape[0])
    radial.require(np.isclose(np.trace(model), np.trace(covariance), rtol=1e-12, atol=0), 'shrinkage changes sample trace')
    eigenvalues = eigvalsh(model, check_finite=False)
    radial.require(eigenvalues[0] >= strength * target * (1 - 1e-10), 'shrinkage loses positive variance bound')
    return {'mean': mean, 'covariance': model, 'factorization': cho_factor(model, lower=True, check_finite=False),
            'target_variance': target, 'samples': len(samples), 'strength': strength,
            'condition_number': float(eigenvalues[-1] / eigenvalues[0])}


def measure(training: np.ndarray, matrices: dict, rings: dict, selected: list,
            model: dict, template: np.ndarray, scale: np.ndarray, half: int) -> tuple[dict, dict]:
    """Freeze a fitted mean and unit-response weights before projecting the opposite half-plane."""
    t = template / scale
    precision_template = cho_solve(model['factorization'], t, check_finite=False)
    energy = float(t @ precision_template)
    radial.require(energy > 0 and np.isfinite(energy), 'invalid amplitude normalization')
    weight, sigma = precision_template / energy, 1 / math.sqrt(energy)
    radial.require(np.isclose(weight @ t, 1, rtol=1e-12) and
                   np.isclose(weight @ model['covariance'] @ weight, sigma*sigma, rtol=1e-10, atol=0), 'conditional variance identity failed')
    opposite = 1 - half
    projected = {str(o): projection.score_stats((matrices[o][rings[o]['halves'] == opposite] - model['mean']) @ weight / sigma)
                 for o in rings}
    heldout_band = np.vstack([matrices[o][rings[o]['halves'] == opposite] for o in selected])
    row = {'sigma': sigma, 'training': projection.score_stats((training - model['mean']) @ weight / sigma),
           'same_radius': projected['0'], 'matched_band': projection.score_stats((heldout_band - model['mean']) @ weight / sigma),
           'rings': projected}
    radial.require(row['same_radius']['patches'] >= 8, 'changed primary validation support')
    # Verify the standardized projection using a separate original-unit covariance solve.
    physical_covariance = scale[:, None] * model['covariance'] * scale[None, :]
    direct = cho_solve(cho_factor(physical_covariance, lower=True, check_finite=False), template, check_finite=False)
    direct_energy = float(template @ direct)
    physical_weight = weight / scale
    radial.require(np.allclose(direct / direct_energy, physical_weight, rtol=1e-9, atol=1e-12) and
                   np.isclose(1 / math.sqrt(direct_energy), sigma, rtol=1e-10, atol=0), 'physical-unit weights disagree')
    return row, {'weight': physical_weight, 'covariance': physical_covariance, 'fit_weight': weight}


def split_stability(parts: list) -> dict:
    """Compare physical-unit weights and covariances from the same site's two training halves."""
    a, b = [p['weight'] for p in parts]
    ca, cb = [p['covariance'] for p in parts]
    return {'weight_cosine': float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b))),
            'relative_covariance_difference': float(2 * np.linalg.norm(ca-cb) / (np.linalg.norm(ca)+np.linalg.norm(cb)))}


def aggregate(rows: list, stability: list) -> dict:
    """Summarize fixed directional projections and expose dispersion as well as their median."""
    result = {'directional_fits': len(rows), 'split_sites': len(stability)}
    for key in ('training', 'same_radius', 'matched_band'):
        values = np.array([r[key]['variance_over_predicted'] for r in rows])
        result[key] = {'median_variance_over_predicted': float(np.median(values)),
            'variance_ratio_10_90_percentiles': np.percentile(values, [10, 90]).tolist(),
            'median_mse_over_predicted': float(np.median([r[key]['mse_over_predicted'] for r in rows])),
            'median_mean_squared_fraction_of_mse': float(np.median([r[key]['mean_squared_fraction_of_mse'] for r in rows])),
            'median_within_one_sigma_fraction': float(np.median([r[key]['within_one_sigma']/r[key]['patches'] for r in rows]))}
    result['median_split_weight_cosine'] = float(np.median([r['weight_cosine'] for r in stability]))
    result['median_split_relative_covariance_difference'] = float(np.median([r['relative_covariance_difference'] for r in stability]))
    result['radial_transfer'] = {}
    for offset in range(-20, 21, 5):
        usable = [r for r in rows if r['rings'].get(str(offset), {}).get('patches', 0) >= 8]
        result['radial_transfer'][str(offset)] = {'directional_fits': len(usable)}
        if usable:
            result['radial_transfer'][str(offset)].update({
                'median_variance_over_predicted': float(np.median([r['rings'][str(offset)]['variance_over_predicted'] for r in usable])),
                'median_paired_ring_over_same_radius_variance': float(np.median([
                    r['rings'][str(offset)]['variance_over_predicted']/r['same_radius']['variance_over_predicted'] for r in usable]))})
    return result


def plots(summary: dict, output: Path) -> None:
    """Show conditional calibration, actual projected noise, and split-weight stability together."""
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5), layout='constrained')
    panels = [
        ('Training variance / prediction', lambda r: r['training']['median_variance_over_predicted'], True),
        ('Held-out variance / prediction', lambda r: r['same_radius']['median_variance_over_predicted'], True),
        ('Held-out amplitude variance / PCA control', lambda r: r['median_paired_variance_over_pca1'], True),
        ('Split-weight cosine: same 16 sites', lambda r: r['median_split_weight_cosine'], False),
        ('Held-out mean² / mean squared error', lambda r: r['same_radius']['median_mean_squared_fraction_of_mse'], False),
        ('Median fraction within ±1 sigma', lambda r: r['same_radius']['median_within_one_sigma_fraction'], False)]
    for axis, (title, getter, logarithmic) in zip(axes.flat, panels):
        values = np.array([[getter(summary['controls'][base])] +
                           [getter(summary['policies'][f'{base}_g{g:g}']) for g in STRENGTHS]
                           for base, _, _ in radial.POLICIES])
        rendered = np.log10(values) if logarithmic else values
        lower, upper = (min(-2, math.floor(rendered.min())), max(2, math.ceil(rendered.max()))) if logarithmic else (0, 1)
        mesh = axis.imshow(rendered, vmin=lower, vmax=upper, cmap='RdBu_r' if logarithmic else 'Blues', aspect='auto')
        axis.set(title=title, xticks=range(4), xticklabels=('PCA f=1', 'γ=0.1', 'γ=0.3', 'γ=1'),
                 yticks=range(8), yticklabels=[f'{"Norm" if n else "Raw"} ±{w}' for _, w, n in radial.POLICIES])
        axis.tick_params(axis='x', labelsize=9)
        for y in range(8):
            for x in range(4):
                dark = abs(rendered[y, x]) > 1.2 if logarithmic else rendered[y, x] > 0.6
                axis.text(x, y, f'{values[y,x]:.2g}', ha='center', va='center', color='white' if dark else 'black')
        fig.colorbar(mesh, ax=axis, shrink=0.7, label='log10 ratio' if logarithmic else 'fraction / cosine')
    fig.suptitle('Step 5 shrinkage: all empirical modes, fixed held-out same-radius patches\n'
                 '16 reused sites × two directions; shared radial profile; development diagnostic')
    fig.savefig(output / 'comparison.png', dpi=170)
    plt.close(fig)


def main() -> None:
    """Freeze the shrinkage grid, reproduce three-mode controls, and evaluate fixed noise projections."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--projected-noise', type=Path, required=True)
    parser.add_argument('--floor-comparison', type=Path, required=True)
    parser.add_argument('--radial-comparison', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prior_root = args.projected_noise.resolve()
    prior = radial.read_json(prior_root / 'manifest.json')
    previous = {(r['site'], r['policy'], r['training_half']): r
                for r in radial.read_json(prior_root / 'projections.json')}
    common = prior['common_sites']
    radial.require(len(common) == 16, 'changed common-site experiment')
    sampling = radial.read_json(args.radial_comparison / 'manifest.json')
    nulls = radial.read_json(args.floor_comparison / 'nulls.json')
    trials = [r for r in nulls if r['name'] in common]
    radial.require(len(trials) == 16, 'missing common trial')
    inputs = [*prior['inputs'], *[fingerprint(prior_root / n) for n in ('manifest.json', 'projections.json', 'summary.json')]]
    for root, name in ((args.radial_comparison, 'manifest.json'), (args.floor_comparison, 'nulls.json')):
        radial.require(fingerprint(root / name) in inputs, 'changed control input')
    for record in inputs:
        radial.require(fingerprint(Path(record['path'])) == record, 'changed original input')
    scripts = [fingerprint(Path(__file__)), *[fingerprint(Path(__file__).with_name(n)) for n in
        ('compare_p4_step5_radial_pooling.py', 'diagnose_p4_step5_projected_noise.py', 'compare_p4_step5_variance_floor.py',
         'analyze_p4_step5_development.py', 'run_p4_step5_full_injections.py')]]
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'manifest.json', {'purpose': 'development shrinkage test on fixed noise projections',
        'formula': '(1-gamma)*sample_covariance + gamma*trace(sample_covariance)/p*I', 'strengths': STRENGTHS,
        'policies': radial.POLICIES, 'common_sites': common, 'directions_per_site': 2,
        'training': prior['training'], 'primary_validation': prior['primary_validation'],
        'normalization': prior['normalization'], 'control': 'three modes, floor fraction 1.0; same mean and native support',
        'comparison': 'conditional variance ratios plus paired absolute amplitude variance and MSE; physical-unit weight stability',
        'sample_support': 'SVD rank cutoff max(n,p)*machine_epsilon*largest singular value; squared weight norm share in sample nullspace is not a variance share',
        'inputs': inputs, 'scripts': scripts})
    science = fits.getdata(sampling['inputs'][0]['path']).squeeze().astype(float)
    field = Path(sampling['inputs'][1]['path']).parent
    yy, xx = np.indices(science.shape)
    forbidden = np.zeros(science.shape, bool)
    for x, y, radius in sampling['exclusions']:
        forbidden |= np.hypot(xx-x, yy-y) <= radius
    scale_map, profile = radial.variance_profile(science, forbidden)
    checks = compare_saved(profile, radial.read_json(args.floor_comparison / 'variance_profiles.json')['baseline'], 'profile')
    coordinates = fits.getdata(field / 'p4PSF_coordinates.fits').T
    responses = fits.getdata(field / 'p4PSF_model_0000.fits').astype(float)
    validity = fits.getdata(field / 'p4PSF_validity_0000.fits').ravel()
    templates = {(int(c[0]), int(c[1])): t.ravel() for c, t, v in zip(coordinates, responses, validity)
                 if v == 1 and np.isfinite(t).all()}
    records, controls, stability, control_stability = [], [], [], []
    endpoint_variances = {}
    for trial in trials:
        position = (trial['row'], trial['column'])
        rings = radial.geometry(science, position, forbidden)
        supports = [np.unique(np.concatenate([projection.native_support(r, r['halves']==h) for r in rings.values()])) for h in (0,1)]
        radial.require(not np.intersect1d(*supports).size, 'train/validation native support overlaps')
        radial.require(not np.any(forbidden.ravel()[np.concatenate(supports)]), 'forbidden native input used')
        x, y = position
        template = templates[position]
        for base, width, normalized in radial.POLICIES:
            image = science / scale_map if normalized else science
            scale = scale_map[y-5:y+6, x-5:x+6].ravel() if normalized else np.ones(121)
            matrices = {o: radial.extract(image, r) for o, r in rings.items()}
            selected = [o for o in rings if abs(o) <= width]
            models_by_policy = {'pca1': []} | {f'{g:g}': [] for g in STRENGTHS}
            for half in (0,1):
                training = np.vstack([matrices[o][rings[o]['halves']==half] for o in selected])
                control_model = radial.fit(training, 1.0)
                radial.require(control_model is not None, 'invalid fixed control')
                control, detail = measure(training, matrices, rings, selected, control_model, template, scale, half)
                saved = previous[(trial['name'], base+'_f1', half)]
                checks += compare_saved(control, {k: saved[k] for k in control}, trial['name'])
                control.update({'site': trial['name'], 'policy': base, 'training_half': half})
                controls.append(control)
                models_by_policy['pca1'].append(detail)
                centered = training - training.mean(axis=0)
                _, singular, vectors = np.linalg.svd(centered, full_matrices=False)
                measured = singular > max(centered.shape) * np.finfo(float).eps * singular[0]
                measured_basis = vectors[measured]
                for gamma in STRENGTHS:
                    model = fit_shrinkage(training, gamma)
                    radial.require(model is not None, 'invalid shrinkage fit on common support')
                    row, detail = measure(training, matrices, rings, selected, model, template, scale, half)
                    weight = detail['fit_weight']
                    sample_variance = row['training']['variance_over_predicted'] * row['sigma']**2
                    isotropic_variance = model['target_variance'] * float(weight @ weight)
                    radial.require(np.isclose((1-gamma)*sample_variance + gamma*isotropic_variance,
                                             row['sigma']**2, rtol=1e-10, atol=0), 'projected shrinkage variance identity failed')
                    row.update({'site': trial['name'], 'policy': f'{base}_g{gamma:g}', 'training_half': half,
                        'gamma': gamma, 'target_variance': model['target_variance'], 'condition_number': model['condition_number'],
                        'paired_sigma_over_pca1': row['sigma']/control['sigma'], 'training_sample_rank': int(measured.sum())})
                    null_weight = weight - measured_basis.T @ (measured_basis @ weight)
                    row['weight_squared_norm_fraction_in_sample_nullspace'] = float(null_weight @ null_weight / (weight @ weight))
                    for metric, label in (('variance_over_predicted', 'variance'), ('mse_over_predicted', 'mse')):
                        row[f'paired_{label}_over_pca1'] = row['same_radius'][metric]*row['sigma']**2 / (
                            control['same_radius'][metric]*control['sigma']**2)
                    if gamma == 1:
                        t = template / scale
                        radial.require(np.allclose(weight, t/(t@t), rtol=1e-10, atol=1e-12), 'isotropic endpoint weight mismatch')
                        variance = row['same_radius']['variance_over_predicted']*row['sigma']**2
                        key = (trial['name'], normalized, half)
                        if key in endpoint_variances:
                            radial.require(np.isclose(variance, endpoint_variances[key], rtol=1e-10, atol=0), 'isotropic variance changes with training band')
                        endpoint_variances[key] = variance
                    records.append(row)
                    models_by_policy[f'{gamma:g}'].append(detail)
            control_pair = split_stability(models_by_policy['pca1'])
            checks += compare_saved(control_pair, {k: trial['policies'][base+'_f1']['center']['split_stability'][k] for k in control_pair}, 'PCA split control')
            control_stability.append({'site': trial['name'], 'policy': base, **control_pair})
            for gamma in STRENGTHS:
                pair = split_stability(models_by_policy[f'{gamma:g}'])
                if gamma == 1:
                    radial.require(np.isclose(pair['weight_cosine'], 1, rtol=1e-12), 'isotropic endpoint halves disagree')
                stability.append({'site': trial['name'], 'policy': f'{base}_g{gamma:g}', **pair})
        print(f'completed {trial["name"]}', flush=True)
    lookup = {(r['site'], r['policy'], r['training_half']): r for r in records}
    for row in records:
        isotropic = lookup[(row['site'], row['policy'].split('_g')[0]+'_g1', row['training_half'])]
        row['paired_variance_over_isotropic'] = row['same_radius']['variance_over_predicted']*row['sigma']**2 / (
            isotropic['same_radius']['variance_over_predicted']*isotropic['sigma']**2)
    summary = {'common_sites': common, 'directional_fits_per_policy': 32, 'controls': {}, 'policies': {}}
    for base, _, _ in radial.POLICIES:
        summary['controls'][base] = aggregate([r for r in controls if r['policy']==base], [r for r in control_stability if r['policy']==base])
        summary['controls'][base]['median_paired_variance_over_pca1'] = 1.0
        for gamma in STRENGTHS:
            name = f'{base}_g{gamma:g}'
            rows = [r for r in records if r['policy']==name]
            one = aggregate(rows, [r for r in stability if r['policy']==name])
            for metric in ('paired_variance_over_pca1', 'paired_mse_over_pca1', 'paired_variance_over_isotropic', 'paired_sigma_over_pca1'):
                one['median_'+metric] = float(np.median([r[metric] for r in rows]))
            one['median_condition_number'] = float(np.median([r['condition_number'] for r in rows]))
            one['training_patch_count_range'] = [min(r['training']['patches'] for r in rows), max(r['training']['patches'] for r in rows)]
            one['training_sample_rank_range'] = [min(r['training_sample_rank'] for r in rows), max(r['training_sample_rank'] for r in rows)]
            one['median_weight_squared_norm_fraction_in_sample_nullspace'] = float(np.median([
                r['weight_squared_norm_fraction_in_sample_nullspace'] for r in rows]))
            summary['policies'][name] = one
    write_json(args.output / 'projections.json', records)
    write_json(args.output / 'stability.json', stability)
    write_json(args.output / 'summary.json', summary)
    for record in [*inputs, *scripts]:
        radial.require(fingerprint(Path(record['path'])) == record, 'input or script changed during experiment')
    verification = {'archived_control_scalar_values': checks, 'shrinkage_directional_fits': len(records),
        'pca_control_directional_fits': len(controls), 'split_comparisons': len(stability),
        'native_train_validation_support_disjoint': True, 'trace_and_positive_variance_bound_preserved': True,
        'projected_shrinkage_variance_identity_agrees': True, 'physical_unit_and_standardized_weights_agree': True,
        'isotropic_weights_and_band_invariant_absolute_variance_agree': True, 'input_and_script_fingerprints_unchanged': True}
    write_json(args.output / 'verification.json', verification)
    plots(summary, args.output)
    write_json(args.output / 'complete.json', {'baseline_images': 1, 'new_reductions': 0, 'new_positive_images': 0,
        'policies': 24, 'control_policies': 8, 'common_sites': 16, 'directional_fits': len(records), 'verification': verification})


if __name__ == '__main__':
    main()
