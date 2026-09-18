#!/usr/bin/env python3
"""Compare fixed variance-floor fractions on the saved Step-5 development images.

Reuse the radial-pooling geometry, profiles, and three-mode covariance model.
Reproduce all archived floor-0.1 measurements, freeze calibration thresholds before
positive filtering, and preserve the original common comparison support.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import compare_p4_step5_radial_pooling as radial
from run_p4_step5_full_injections import fingerprint, write_json

FLOORS = (0.1, 0.3, 1.0)


def policy_name(base: str, fraction: float) -> str:
    """Name each sampling/normalization/floor combination without changing the original base names."""
    return f'{base}_f{fraction:g}'


def compare_saved(actual, expected, path: str = '') -> int:
    """Check every archived measurement, decision, count, and diagnostic at the original floor."""
    if isinstance(expected, dict):
        return sum(compare_saved(actual[key], value, f'{path}/{key}') for key, value in expected.items())
    if isinstance(expected, list):
        radial.require(len(actual) == len(expected), f'control list length mismatch: {path}')
        return sum(compare_saved(a, e, f'{path}/{i}') for i, (a, e) in enumerate(zip(actual, expected)))
    if isinstance(expected, float):
        radial.require(np.isclose(actual, expected, rtol=1e-10, atol=1e-12), f'control numeric mismatch: {path}')
    else:
        radial.require(actual == expected, f'control value mismatch: {path}')
    return 1


def analyze_grid(science: np.ndarray, scale: np.ndarray, position: tuple, template: np.ndarray,
                 rings: dict, split: bool = False) -> tuple[dict, dict]:
    """Apply every floor and verify nondecreasing conditional sigma on identical training data."""
    rows, details, previous = {}, {}, {}
    for fraction in FLOORS:
        measured, models = radial.analyze_position(science, scale, position, template, rings,
                                                   split=split, floor_fraction=fraction)
        for base, value in measured.items():
            if base in previous:
                radial.require(value['valid'] == previous[base]['valid'] and
                               value['samples'] == previous[base]['samples'], 'floor changes eligibility')
                if value['valid']:
                    radial.require(value['sigma'] >= previous[base]['sigma'] * (1 - 1e-12),
                                   'conditional sigma decreases as covariance floor increases')
            rows[policy_name(base, fraction)] = value
        details.update({policy_name(base, fraction): model for base, model in models.items()})
        previous = measured
    return rows, details


def at_floor(records: list, fraction: float) -> list:
    """Project the combined experiment onto the original eight-policy summary interface."""
    return [{**row, 'policies': {base: row['policies'][policy_name(base, fraction)]
                               for base, _, _ in radial.POLICIES}} for row in records]


def plots(summary: dict, output: Path) -> None:
    """Plot recovery, uncertainty, and stability across the prespecified floor grid."""
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout='constrained')
    panels = [('Faint recoveries / 6', lambda p: p['groups'][0]['detections'], 0, 6, 'd'),
              ('Null exceedances / 28', lambda p: p['evaluation_exceedances'], 0, 4, 'd'),
              ('Positive ±1-sigma coverage / 18', lambda p: p['positive_one_sigma'], 0, 18, 'd'),
              ('Null-center ±1-sigma coverage / 28', lambda p: p['null_center_one_sigma'], 0, 28, 'd'),
              ('Weight cosine: same 16 split sites', lambda p: p['common_median_split_weight_cosine'], 0.5, 1, '.2f'),
              ('Modeled / sample covariance trace', lambda p: p['median_modeled_sample_trace_fraction'], 0, 1.5, '.2f')]
    for axis, (title, getter, lower, upper, number_format) in zip(axes.flat, panels):
        values = np.array([[getter(summary['policies'][policy_name(base, f)]) for f in FLOORS]
                           for base, _, _ in radial.POLICIES])
        mesh = axis.imshow(values, vmin=lower, vmax=upper, cmap='Blues', aspect='auto')
        axis.set(title=title, xticks=range(3), xticklabels=FLOORS, xlabel='Variance-floor fraction',
                 yticks=range(8), yticklabels=[f'{"Norm" if n else "Raw"} ±{w}' for _, w, n in radial.POLICIES])
        for y in range(8):
            for x in range(3):
                value = int(values[y, x]) if number_format == 'd' else values[y, x]
                axis.text(x, y, format(value, number_format), ha='center', va='center',
                          color='white' if values[y, x] > lower + 0.6 * (upper - lower) else 'black')
        fig.colorbar(mesh, ax=axis, shrink=0.7)
    fig.suptitle('Step 5 development: fixed rank 3, saved images reused\n'
                 'Identity / Gaussian-SNR references: 3/6 faint recoveries, 2/28 null exceedances')
    fig.savefig(output / 'comparison.png', dpi=170)
    plt.close(fig)


def main() -> None:
    """Run the complete saved-image experiment and check the archived floor-0.1 control."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--radial-comparison', type=Path, required=True)
    parser.add_argument('--evaluation', type=Path, required=True)
    parser.add_argument('--gaussian', type=Path, required=True)
    parser.add_argument('--development', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prior = args.radial_comparison.resolve()
    previous_manifest = radial.read_json(prior / 'manifest.json')
    for record in previous_manifest['inputs']:
        radial.require(fingerprint(Path(record['path'])) == record, 'original experiment input changed')
    evaluation = args.evaluation.resolve()
    protocol = radial.read_json(evaluation / 'protocol.json')
    original_nulls = radial.read_json(evaluation / 'null_results.json')
    original = radial.read_json(evaluation / 'results.json')
    gaussian = radial.read_json(args.gaussian / 'results.json')
    saved_nulls = {t['name']: t for t in radial.read_json(prior / 'nulls.json')}
    saved_positives = {t['trial']['name']: t for t in radial.read_json(prior / 'injections.json')}
    saved_summary = radial.read_json(prior / 'summary.json')
    saved_profiles = radial.read_json(prior / 'variance_profiles.json')
    science = fits.getdata(protocol['science']['path']).squeeze().astype(float)
    yy, xx = np.indices(science.shape)
    forbidden = np.zeros(science.shape, bool)
    for x, y, radius in previous_manifest['exclusions']:
        forbidden |= np.hypot(xx - x, yy - y) <= radius
    baseline_scale, baseline_profile = radial.variance_profile(science, forbidden)
    control_checks = compare_saved(baseline_profile, saved_profiles['baseline'], 'baseline profile')
    field = Path(protocol['source_response_inputs'][0]['path']).parent
    coordinates = fits.getdata(field / 'p4PSF_coordinates.fits').T
    responses = fits.getdata(field / 'p4PSF_model_0000.fits').astype(float)
    validity = fits.getdata(field / 'p4PSF_validity_0000.fits').ravel()
    templates = {(int(c[0]), int(c[1])): t.ravel() for c, t, v in zip(coordinates, responses, validity)
                 if v == 1 and np.isfinite(t).all()}
    jobs = [{'trial': j, 'role': 'original_evaluation', 'image': evaluation / 'reductions' / j['name'] / 'finim.fits'}
            for j in radial.read_json(evaluation / 'jobs.json')]
    development_manifest = radial.read_json(args.development / 'manifest.json')
    jobs += [{'trial': j, 'role': 'original_development', 'image': args.development.resolve() / j['name'] / 'finim.fits'}
             for j in development_manifest['trials'] if j['contrast'] > 0]
    centers = {(t['row'], t['column']) for t in protocol['trials']} | {(j['trial']['row'], j['trial']['column']) for j in jobs}
    positions = sorted({(x + dx, y + dy) for x, y in centers for dx, dy in radial.OFFSETS})
    positive_positions = {(j['trial']['row'] + dx, j['trial']['column'] + dy) for j in jobs for dx, dy in radial.OFFSETS}
    policies = [policy_name(base, f) for f in FLOORS for base, _, _ in radial.POLICIES]
    inputs = [*previous_manifest['inputs'], *[fingerprint(prior / name) for name in
              ('manifest.json', 'nulls.json', 'injections.json', 'summary.json', 'thresholds.json', 'variance_profiles.json')],
              fingerprint(evaluation / 'jobs.json')]
    scripts = [fingerprint(Path(__file__)), fingerprint(Path(radial.__file__)),
               fingerprint(Path(radial.__file__).with_name('analyze_p4_step5_development.py')),
               fingerprint(Path(radial.__file__).with_name('run_p4_step5_full_injections.py'))]
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output / 'manifest.json', {'purpose': 'development reuse; fixed variance-floor comparison',
        'floor_fractions': FLOORS, 'maximum_modes': 3, 'policies': policies,
        'sampling_and_normalization': previous_manifest['policies'], 'original_geometry': str(prior / 'manifest.json'),
        'comparison_support': previous_manifest['comparison_support'],
        'threshold_rule': 'maximum of the original 28 calibration five-pixel searches; strict exceedance',
        'inputs': inputs, 'scripts': scripts})
    baseline_pixels, baseline_details, cached_geometry = {}, {}, {}
    holdout = Path(protocol['software'][0]['path']).parent.parent
    production_maps = {role: fits.getdata(holdout / f'analysis/pca/science_{role}.fits').squeeze()
                       for role in ('psf_amplitude', 'psf_sigma', 'psf_score', 'noise_samples', 'noise_status')}
    original_positions = {(t['row'] + dx, t['column'] + dy) for t in protocol['trials'] for dx, dy in radial.OFFSETS}
    production_checks, production_errors = 0, []
    for index, position in enumerate(positions):
        rings = radial.geometry(science, position, forbidden)
        rows, details = analyze_grid(science, baseline_scale, position, templates[position], rings, split=position in centers)
        baseline_pixels[position] = rows
        if position in positive_positions:
            cached_geometry[position] = rings
        if position in centers:
            baseline_details[position] = details
        if position in original_positions:
            x, y = position
            control = rows[policy_name('raw_b0', 0.1)]
            radial.require(control['samples'] == int(production_maps['noise_samples'][y, x]) and
                           control['valid'] == (production_maps['noise_status'][y, x] == 0), 'production geometry mismatch')
            if control['valid']:
                for key, role in (('amplitude', 'psf_amplitude'), ('sigma', 'psf_sigma'), ('score', 'psf_score')):
                    expected = float(production_maps[role][y, x])
                    radial.require(np.isclose(control[key], expected, rtol=2e-5, atol=1e-10), 'production numeric mismatch')
                    production_errors.append(abs(control[key] - expected) / max(abs(expected), 1e-12))
            production_checks += 1
        if index % 50 == 0:
            print(f'baseline candidate {index + 1}/{len(positions)}', flush=True)
    nulls = [{**{k: t[k] for k in ('name', 'role', 'nominal_radius', 'angular_block', 'row', 'column')},
              'original_common_eligible': t['common_eligible'],
              'policies': {name: radial.trial_measurement(t, baseline_pixels, name) for name in policies}}
             for t in original_nulls['trials']]
    for row in at_floor(nulls, 0.1):
        control_checks += compare_saved(row, saved_nulls[row['name']], row['name'])
    thresholds = {}
    for name in policies:
        calibration = [t['policies'][name] for t in nulls if t['role'] == 'calibration' and t['original_common_eligible']]
        radial.require(len(calibration) == 28 and all(c['valid'] for c in calibration), 'invalid calibration support')
        thresholds[name] = max(c['search_score'] for c in calibration)
    write_json(args.output / 'thresholds.json', thresholds)
    threshold_record = fingerprint(args.output / 'thresholds.json')
    write_json(args.output / 'nulls.json', nulls)
    positives, profiles = [], {'baseline': baseline_profile}
    for index, job in enumerate(jobs):
        trial = job['trial']
        positive = fits.getdata(job['image']).squeeze().astype(float)
        radial.require(np.array_equal(np.isfinite(science), np.isfinite(positive)), 'changed finite support')
        scale, profiles[trial['name']] = radial.variance_profile(positive, forbidden)
        control_checks += compare_saved(profiles[trial['name']], saved_profiles[trial['name']], 'positive profile')
        pixels, central_details = {}, {}
        position = (trial['row'], trial['column'])
        for dx, dy in radial.OFFSETS:
            p = (position[0] + dx, position[1] + dy)
            pixels[p], details = analyze_grid(positive, scale, p, templates[p], cached_geometry[p])
            if dx == dy == 0:
                central_details = details
        measured = {}
        for policy in policies:
            value = radial.trial_measurement(trial, pixels, policy)
            value['detected'] = value['valid'] and value['search_score'] > thresholds[policy]
            center = value['center']
            value['inside_one_conditional_sigma'] = center['valid'] and abs(center['amplitude'] - trial['contrast']) <= center['sigma']
            base = baseline_pixels[position][policy]
            if center['valid']:
                value['raw_contrast_error'] = center['amplitude'] / trial['contrast'] - 1
            if center['valid'] and base['valid']:
                value['paired_increment_error'] = (center['amplitude'] - base['amplitude']) / trial['contrast'] - 1
                before, after = baseline_details[position][policy], central_details[policy]
                centered = before['samples'] - before['samples'].mean(axis=0)
                value['training_change_over_centered_baseline_norm'] = float(np.linalg.norm(after['samples'] - before['samples']) / np.linalg.norm(centered))
                value['physical_covariance_relative_change'] = float(np.linalg.norm(after['physical_covariance'] - before['physical_covariance']) / np.linalg.norm(before['physical_covariance']))
            measured[policy] = value
        row = {'trial': trial, 'role': job['role'], 'policies': measured}
        control_checks += compare_saved(at_floor([row], 0.1)[0], saved_positives[trial['name']], trial['name'])
        positives.append(row)
        print(f'positive image {index + 1}/{len(jobs)}: {trial["name"]}', flush=True)
    summary = {'policies': {}}
    for fraction in FLOORS:
        one = radial.summarize(at_floor(nulls, fraction), at_floor(positives, fraction),
                              {base: thresholds[policy_name(base, fraction)] for base, _, _ in radial.POLICIES}, original, gaussian)
        if fraction == 0.1:
            control_checks += compare_saved(one, saved_summary, 'summary')
            summary.update({key: value for key, value in one.items() if key != 'policies'})
        radial.require(one['common_split_stability_trials'] == summary['common_split_stability_trials'], 'floor changes stability support')
        summary['policies'].update({policy_name(base, fraction): value for base, value in one['policies'].items()})
    write_json(args.output / 'summary.json', summary)
    write_json(args.output / 'injections.json', positives)
    write_json(args.output / 'variance_profiles.json', profiles)
    for record in [*inputs, *scripts, threshold_record]:
        radial.require(fingerprint(Path(record['path'])) == record, 'frozen input, script, or threshold changed')
    verification = {'archived_floor_0p1_scalar_controls': control_checks, 'production_baseline_pixels_checked': production_checks,
        'maximum_production_control_relative_difference': max(production_errors), 'original_positive_images_checked': len(jobs),
        'sigma_nondecreasing_with_floor_at_every_valid_fit': True, 'eligibility_unchanged_with_floor': True,
        'normalized_and_original_unit_solves_agree': True, 'common_split_stability_centers': len(summary['common_split_stability_trials']),
        'source_script_threshold_fingerprints_unchanged': True}
    write_json(args.output / 'verification.json', verification)
    plots(summary, args.output)
    write_json(args.output / 'complete.json', {'baseline_images': 1, 'positive_images': len(jobs), 'policies': len(policies),
        'new_reductions': 0, 'purpose': 'development reuse, not blind validation', 'verification': verification})
    for name, value in summary['policies'].items():
        print(name, json.dumps(value), flush=True)


if __name__ == '__main__':
    main()
