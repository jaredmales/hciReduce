#!/usr/bin/env python3
"""Evaluate the fixed Welch-style PSD grid on saved Step-5 null and injection images.

Freeze calibration-only thresholds before positive filtering. Retain all original
comparison sites, reproduce PCA controls and Gaussian/identity references, and
report native-center photometry separately from five-pixel search recovery.
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
from scipy.linalg import cho_factor, cho_solve

import compare_p4_step5_radial_pooling as radial
import compare_p4_step5_shrinkage as shrinkage
import compare_p4_step5_welch_psd as psd
from compare_p4_step5_variance_floor import compare_saved
from diagnose_p4_step5_projected_noise import score_stats
from run_p4_step5_full_injections import fingerprint, write_json

VARIANTS = [('pca', None), ('isotropic', None), *psd.VARIANTS]


def policy_name(base: str, kind: str, mixing: float | None) -> str:
    """Give the controls and spectral variants stable, distinct names."""
    return base + '_pca1' if kind == 'pca' else base + '_isotropic' if kind == 'isotropic' else psd.policy_name(base, kind, mixing)


def analyze_grid(science: np.ndarray, scale_map: np.ndarray, position: tuple,
                 template: np.ndarray, rings: dict) -> tuple[dict, dict]:
    """Fit all accepted patches outside the fixed exclusions, then filter the native candidate."""
    x, y = position
    data = science[y-5:y+6, x-5:x+6].ravel()
    scale = scale_map[y-5:y+6, x-5:x+6].ravel()
    radial.require(np.isfinite(data).all() and np.isfinite(template).all() and np.isfinite(scale).all(),
                   'invalid native candidate support')
    matrices = {normalized: {o: radial.extract(science/scale_map if normalized else science, ring)
                            for o, ring in rings.items()} for normalized in (False, True)}
    rows, details = {}, {}
    for base, width, normalized in radial.POLICIES:
        selected = [o for o in rings if abs(o) <= width]
        training = np.vstack([matrices[normalized][o] for o in selected])
        factors = scale if normalized else np.ones(121)
        validity = []
        for kind, mixing in VARIANTS:
            model = (radial.fit(training, 1) if kind == 'pca' else
                     shrinkage.fit_shrinkage(training, 1) if kind == 'isotropic' else psd.fit_psd(training, kind, mixing))
            name = policy_name(base, kind, mixing)
            rows[name] = {'valid': model is not None, 'samples': len(training)}
            validity.append(model is not None)
            if model is None:
                continue
            measured = radial.filter_stamp(data, template, model, factors)
            covariance = factors[:, None] * model['covariance'] * factors[None, :]
            mean = factors * model['mean']
            weight = measured['physical_weight']
            radial.require(np.isclose(weight @ template, 1, rtol=1e-12, atol=0) and
                           np.isclose(weight @ covariance @ weight, measured['sigma']**2, rtol=1e-10, atol=0),
                           'unit-response or conditional-variance identity failed')
            direct = cho_solve(cho_factor(covariance, lower=True, check_finite=False), template, check_finite=False)
            energy = float(template @ direct)
            radial.require(np.allclose(direct/energy, weight, rtol=1e-9, atol=1e-12) and
                           np.isclose((direct @ (data-mean))/energy, measured['amplitude'], rtol=1e-9, atol=1e-12) and
                           np.isclose(1/math.sqrt(energy), measured['sigma'], rtol=1e-10, atol=0),
                           'physical-unit solve disagrees')
            rows[name].update({k: measured[k] for k in ('amplitude', 'sigma', 'score')})
            details[name] = {'samples': training, 'physical_covariance': covariance, 'physical_weight': weight}
        radial.require(len(set(validity)) == 1, 'estimator changes common sample eligibility')
    return rows, details


def summarize(nulls: list, positives: list, thresholds: dict, policies: list) -> dict:
    """Keep original support, raw photometry, paired increments, and native-null scores distinct."""
    evaluation = [n for n in nulls if n['role'] == 'evaluation' and n['original_common_eligible']]
    injections = [p for p in positives if p['role'] == 'original_evaluation']
    radial.require(len(evaluation) == 28 and len(injections) == 18, 'changed primary comparison support')
    summary = {}
    for name in policies:
        centers = [n['policies'][name]['center'] for n in evaluation]
        radial.require(all(n['policies'][name]['valid'] for n in evaluation), 'invalid primary null')
        result = {'threshold': thresholds[name], 'evaluation_nulls': 28,
                  'evaluation_exceedances': sum(n['policies'][name]['search_score'] > thresholds[name] for n in evaluation),
                  'null_center_one_sigma': sum(abs(c['amplitude']) <= c['sigma'] for c in centers),
                  'native_null_center_scores': score_stats(np.array([c['score'] for c in centers])),
                  'positive_one_sigma': sum(p['policies'][name]['inside_one_conditional_sigma'] for p in injections),
                  'new_radius20_valid_searches': sum(n['nominal_radius'] == 20 and n['policies'][name]['valid'] for n in nulls),
                  'groups': []}
        for level in (.5, 1, 2):
            selected = [p for p in injections if p['trial']['brightness_multiplier'] == level]
            values = [p['policies'][name] for p in selected]
            radial.require(len(values) == 6 and all(v['valid'] for v in values), 'changed injection support')
            errors = np.array([v['raw_contrast_error'] for v in values])
            result['groups'].append({'brightness_multiplier': level, 'trials': 6,
                'detections': sum(v['detected'] for v in values), 'median_raw_contrast_error': float(np.median(errors)),
                'raw_contrast_error_range': [float(errors.min()), float(errors.max())],
                'median_absolute_raw_contrast_error': float(np.median(np.abs(errors))),
                'rms_raw_contrast_error': float(np.sqrt(np.mean(errors**2))),
                'inside_one_conditional_sigma': sum(v['inside_one_conditional_sigma'] for v in values)})
        for field in ('raw_contrast_error', 'paired_increment_error', 'frozen_weight_increment_error',
                      'training_change_over_centered_baseline_norm', 'physical_covariance_relative_change'):
            values = np.array([p['policies'][name][field] for p in injections])
            result['median_'+field] = float(np.median(values))
            result['median_absolute_'+field] = float(np.median(np.abs(values)))
            result[field+'_range'] = [float(values.min()), float(values.max())]
        summary[name] = result
    return summary


def references(original: dict, gaussian: dict, nulls: dict) -> dict:
    """Recompute saved Gaussian and identity decisions with their unchanged calibration thresholds."""
    result = {}
    for model in ('identity', 'gaussian_snr', 'gaussian_raw', 'identity_snr'):
        calibration = [r['models'][model]['search_score'] for r in gaussian['null_trials']
                       if r['common_eligible'] and r['role'] == 'calibration']
        evaluation = [r for r in gaussian['null_trials'] if r['common_eligible'] and r['role'] == 'evaluation']
        radial.require(len(calibration) == len(evaluation) == 28, 'changed reference support')
        threshold = max(calibration)
        rows = [r for r in gaussian['injections'] if r['model'] == model]
        radial.require(len(rows) == 18 and all(r['detected'] == (r['search_score'] > threshold) for r in rows),
                       'reference decisions changed')
        counts = [sum(r['detected'] for r in rows if r['trial']['brightness_multiplier'] == level) for level in (.5, 1, 2)]
        one = {'threshold': threshold, 'evaluation_exceedances': sum(r['models'][model]['search_score'] > threshold for r in evaluation),
               'detections_by_brightness': counts}
        compare_saved(one, {k: gaussian['models'][model][k] for k in one}, 'reference '+model)
        if model == 'identity':
            measured = [r for r in original['measurements'] if r['model'] == model]
            one['positive_one_sigma'] = sum(r['inside_one_conditional_sigma'] for r in measured)
            one['null_center_one_sigma'] = sum(abs(r['models'][model]['center_amplitude']) <= r['models'][model]['center_conditional_sigma']
                                               for r in nulls['trials'] if r['common_eligible'] and r['role'] == 'evaluation')
            one['groups'] = [r for r in original['groups'] if r['model'] == model]
            one['sigma_interpretation'] = 'identity uses C=I; its original conditional sigma is an algebraic normalization, not fitted physical noise'
        result[model] = one
    return result


def plots(summary: dict, output: Path) -> None:
    """Show every fixed PSD setting and both covariance controls beside the detection references."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
    panels = [('Faint recovery / 6', lambda p: p['groups'][0]['detections'], 0, 6, 'd'),
              ('Middle-brightness recovery / 6', lambda p: p['groups'][1]['detections'], 0, 6, 'd'),
              ('Evaluation null exceedances / 28', lambda p: p['evaluation_exceedances'], 0, 6, 'd'),
              ('Positive ±1-sigma count / 18', lambda p: p['positive_one_sigma'], 0, 18, 'd'),
              ('Native null-center ±1-sigma count / 28', lambda p: p['null_center_one_sigma'], 0, 28, 'd'),
              ('Median absolute raw contrast error', lambda p: p['median_absolute_raw_contrast_error'], 0, 1, '.2f')]
    for axis, (title, getter, lo, hi, fmt) in zip(axes.flat, panels):
        values = np.array([[getter(summary['policies'][policy_name(base, k, m)]) for k, m in VARIANTS]
                           for base, _, _ in radial.POLICIES])
        hi = max(hi, float(values.max()))
        mesh = axis.imshow(values, vmin=lo, vmax=hi, cmap='Blues', aspect='auto')
        axis.set(title=title, xticks=range(6), xticklabels=('PCA f=1', 'Iso', 'Rect .1', 'Rect .3', 'Hann .1', 'Hann .3'),
                 yticks=range(8), yticklabels=[f'{"Norm" if n else "Raw"} ±{w}' for _, w, n in radial.POLICIES])
        axis.tick_params(axis='x', labelsize=9)
        for y in range(8):
            for x in range(6):
                value = int(values[y, x]) if fmt == 'd' else values[y, x]
                axis.text(x, y, format(value, fmt), ha='center', va='center', fontsize=9,
                          color='white' if values[y, x] > .6*hi else 'black')
        fig.colorbar(mesh, ax=axis, shrink=.7)
    refs = summary['references']
    fig.suptitle('Step 5 PSD: original 28+28 null searches and 18 saved injections; development reuse\n'
                 f'Identity / Gaussian-SNR: {refs["identity"]["detections_by_brightness"]} / '
                 f'{refs["gaussian_snr"]["detections_by_brightness"]} detections; '
                 f'{refs["identity"]["evaluation_exceedances"]} / {refs["gaussian_snr"]["evaluation_exceedances"]} null exceedances')
    fig.savefig(output/'comparison.png', dpi=170)
    plt.close(fig)


def main() -> None:
    """Freeze the full PSD grid, reproduce controls, set null thresholds, and process saved positives."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('psd-comparison', 'floor-comparison', 'radial-comparison', 'evaluation', 'gaussian', 'development', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    previous = radial.read_json(args.psd_comparison/'manifest.json')
    radial.require(previous['variants'] == [{'window': w, 'spectral_isotropic_mixing': m} for w, m in psd.VARIANTS], 'changed PSD grid')
    floor_manifest = radial.read_json(args.floor_comparison/'manifest.json')
    sampling = radial.read_json(args.radial_comparison/'manifest.json')
    protocol = radial.read_json(args.evaluation/'protocol.json')
    original_nulls = radial.read_json(args.evaluation/'null_results.json')
    original = radial.read_json(args.evaluation/'results.json')
    gaussian = radial.read_json(args.gaussian/'results.json')
    saved_nulls = {r['name']: r for r in radial.read_json(args.floor_comparison/'nulls.json')}
    saved_positives = {r['trial']['name']: r for r in radial.read_json(args.floor_comparison/'injections.json')}
    saved_summary = radial.read_json(args.floor_comparison/'summary.json')
    saved_profiles = radial.read_json(args.floor_comparison/'variance_profiles.json')
    jobs = [{'trial': j, 'role': 'original_evaluation', 'image': args.evaluation.resolve()/'reductions'/j['name']/'finim.fits'}
            for j in radial.read_json(args.evaluation/'jobs.json')]
    development = radial.read_json(args.development/'manifest.json')
    jobs += [{'trial': j, 'role': 'original_development', 'image': args.development.resolve()/j['name']/'finim.fits'}
             for j in development['trials'] if j['contrast'] > 0]
    inputs = [*floor_manifest['inputs'], *[fingerprint(args.floor_comparison/n) for n in
              ('manifest.json', 'nulls.json', 'injections.json', 'summary.json', 'variance_profiles.json')],
              fingerprint(args.psd_comparison/'manifest.json'), fingerprint(args.psd_comparison/'summary.json')]
    for root, names in ((args.evaluation, ('protocol.json', 'null_results.json', 'results.json', 'jobs.json')),
                        (args.gaussian, ('results.json',)), (args.development, ('manifest.json',)),
                        (args.radial_comparison, ('manifest.json',))):
        for name in names:
            radial.require(fingerprint(root/name) in inputs, 'unfrozen experiment source')
    for job in jobs:
        radial.require(fingerprint(job['image']) in inputs, 'unfrozen positive image')
    for item in [*inputs, *previous['scripts']]:
        radial.require(fingerprint(Path(item['path'])) == item, 'changed source or PSD estimator')
    scripts = [fingerprint(Path(__file__)), *previous['scripts']]
    policies = [policy_name(base, k, m) for base, _, _ in radial.POLICIES for k, m in VARIANTS]
    args.output.mkdir(parents=True, exist_ok=False)
    write_json(args.output/'manifest.json', {'purpose': 'development-only PSD recovery comparison on saved images',
        'policies': policies, 'psd_variants': previous['variants'], 'sampling_policies': radial.POLICIES,
        'estimator': 'unchanged zero-padded Welch-style PSD; window only inside covariance estimation',
        'training': 'all accepted annular/band patches, including split straddlers; original native source/holdout union excluded',
        'comparison_support': sampling['comparison_support'], 'positive_jobs': [{**j, 'image': str(j['image'])} for j in jobs],
        'profile': sampling['profile'], 'exclusions': sampling['exclusions'],
        'threshold_rule': 'maximum of original 28 calibration five-pixel searches; strict exceedance; frozen before positive filtering',
        'photometry': 'exact native center; raw error and conditional sigma; paired/adaptive and frozen-weight increments separate',
        'controls': 'eight three-mode floor-1 controls, eight fitted-mean isotropic controls; unchanged Gaussian and identity detection references',
        'inputs': inputs, 'scripts': scripts})
    science = fits.getdata(protocol['science']['path']).squeeze().astype(float)
    yy, xx = np.indices(science.shape)
    forbidden = np.zeros(science.shape, bool)
    for x, y, radius in sampling['exclusions']:
        forbidden |= np.hypot(xx-x, yy-y) <= radius
    baseline_scale, baseline_profile = radial.variance_profile(science, forbidden)
    checks = compare_saved(baseline_profile, saved_profiles['baseline'], 'baseline profile')
    field = Path(protocol['source_response_inputs'][0]['path']).parent
    coordinates = fits.getdata(field/'p4PSF_coordinates.fits').T
    responses = fits.getdata(field/'p4PSF_model_0000.fits').astype(float)
    validity = fits.getdata(field/'p4PSF_validity_0000.fits').ravel()
    templates = {(int(c[0]), int(c[1])): t.ravel() for c, t, v in zip(coordinates, responses, validity) if v == 1 and np.isfinite(t).all()}
    positive_centers = {(j['trial']['row'], j['trial']['column']) for j in jobs}
    centers = {(t['row'], t['column']) for t in protocol['trials']} | positive_centers
    positions = sorted({(x+dx, y+dy) for x, y in centers for dx, dy in radial.OFFSETS})
    positive_positions = {(x+dx, y+dy) for x, y in positive_centers for dx, dy in radial.OFFSETS}
    baseline_pixels, baseline_details, cached = {}, {}, {}
    fits_count = 0
    for index, position in enumerate(positions):
        rings = radial.geometry(science, position, forbidden)
        for ring in rings.values():
            used = ring['indices'][ring['weights'] != 0]
            radial.require(not np.any(forbidden.ravel()[used]), 'training uses an excluded native pixel')
            x, y = position
            radial.require(not np.any((np.abs(used % science.shape[1]-x) <= 5) &
                                      (np.abs(used // science.shape[1]-y) <= 5)), 'candidate stamp enters training')
        rows, details = analyze_grid(science, baseline_scale, position, templates[position], rings)
        baseline_pixels[position] = rows
        fits_count += sum(r['valid'] for r in rows.values())
        if position in positive_centers:
            baseline_details[position] = details
        if position in positive_positions:
            cached[position] = rings
        if index % 50 == 0:
            print(f'baseline candidate {index+1}/{len(positions)}', flush=True)
    nulls = [{**{k: t[k] for k in ('name', 'role', 'nominal_radius', 'angular_block', 'row', 'column')},
              'original_common_eligible': t['common_eligible'],
              'policies': {name: radial.trial_measurement(t, baseline_pixels, name) for name in policies}}
             for t in original_nulls['trials']]
    for trial in nulls:
        for base, _, _ in radial.POLICIES:
            row = trial['policies'][base+'_pca1']
            expected = saved_nulls[trial['name']]['policies'][base+'_f1']
            checks += compare_saved(row, {**{k: expected[k] for k in row if k != 'center'},
                                         'center': {k: expected['center'][k] for k in row['center']}}, 'PCA null')
    thresholds = {}
    for name in policies:
        calibration = [n['policies'][name] for n in nulls if n['role'] == 'calibration' and n['original_common_eligible']]
        radial.require(len(calibration) == 28 and all(c['valid'] for c in calibration), 'invalid calibration support')
        thresholds[name] = max(c['search_score'] for c in calibration)
    write_json(args.output/'thresholds.json', thresholds)
    threshold_record = fingerprint(args.output/'thresholds.json')
    write_json(args.output/'nulls.json', nulls)
    write_json(args.output/'baseline_pixels.json', [{'position': p, 'policies': rows} for p, rows in baseline_pixels.items()])
    positives, profiles = [], {'baseline': baseline_profile}
    for index, job in enumerate(jobs):
        trial = job['trial']
        positive = fits.getdata(job['image']).squeeze().astype(float)
        radial.require(np.array_equal(np.isfinite(science), np.isfinite(positive)), 'positive changes finite support')
        scale, profiles[trial['name']] = radial.variance_profile(positive, forbidden)
        checks += compare_saved(profiles[trial['name']], saved_profiles[trial['name']], 'positive profile')
        x, y = position = (trial['row'], trial['column'])
        pixels, central_details = {}, {}
        for dx, dy in radial.OFFSETS:
            p = (x+dx, y+dy)
            pixels[p], details = analyze_grid(positive, scale, p, templates[p], cached[p])
            fits_count += sum(r['valid'] for r in pixels[p].values())
            if dx == dy == 0:
                central_details = details
        delta = (positive[y-5:y+6, x-5:x+6]-science[y-5:y+6, x-5:x+6]).ravel()
        measured = {}
        for name in policies:
            value = radial.trial_measurement(trial, pixels, name)
            value['detected'] = value['valid'] and value['search_score'] > thresholds[name]
            center = value['center']
            value['inside_one_conditional_sigma'] = center['valid'] and abs(center['amplitude']-trial['contrast']) <= center['sigma']
            before = baseline_pixels[position][name]
            if center['valid']:
                value['raw_contrast_error'] = center['amplitude']/trial['contrast']-1
            if center['valid'] and before['valid']:
                value['paired_increment_error'] = (center['amplitude']-before['amplitude'])/trial['contrast']-1
                baseline = baseline_details[position][name]
                after = central_details[name]
                value['frozen_weight_increment_error'] = float(baseline['physical_weight'] @ delta/trial['contrast']-1)
                centered = baseline['samples']-baseline['samples'].mean(axis=0)
                value['training_change_over_centered_baseline_norm'] = float(np.linalg.norm(after['samples']-baseline['samples'])/np.linalg.norm(centered))
                value['physical_covariance_relative_change'] = float(np.linalg.norm(after['physical_covariance']-baseline['physical_covariance'])/np.linalg.norm(baseline['physical_covariance']))
            measured[name] = value
        for base, _, _ in radial.POLICIES:
            row = measured[base+'_pca1']
            expected = saved_positives[trial['name']]['policies'][base+'_f1']
            checks += compare_saved(row, {**{k: expected[k] for k in row if k not in ('center', 'frozen_weight_increment_error')},
                                         'center': {k: expected['center'][k] for k in row['center']}}, 'PCA positive')
        positives.append({'trial': trial, 'role': job['role'], 'policies': measured,
                          'search_pixels': [{'position': p, 'policies': rows} for p, rows in pixels.items()]})
        print(f'positive image {index+1}/{len(jobs)}: {trial["name"]}', flush=True)
    summary = {'policies': summarize(nulls, positives, thresholds, policies),
               'references': references(original, gaussian, original_nulls)}
    for base, _, _ in radial.POLICIES:
        actual = summary['policies'][base+'_pca1']
        expected = saved_summary['policies'][base+'_f1']
        for key in ('threshold', 'evaluation_nulls', 'evaluation_exceedances', 'null_center_one_sigma',
                    'positive_one_sigma', 'new_radius20_valid_searches'):
            checks += compare_saved(actual[key], expected[key], 'PCA summary '+key)
        for a, e in zip(actual['groups'], expected['groups']):
            checks += compare_saved(a, {k: v for k, v in e.items() if k != 'invalid'}, 'PCA group')
    write_json(args.output/'summary.json', summary)
    write_json(args.output/'injections.json', positives)
    write_json(args.output/'variance_profiles.json', profiles)
    for item in [*inputs, *scripts, threshold_record]:
        radial.require(fingerprint(Path(item['path'])) == item, 'frozen input, estimator, or threshold changed')
    verification = {'archived_pca_control_and_profile_scalars': checks, 'valid_candidate_fits': fits_count,
        'baseline_search_pixels': len(positions), 'positive_images': len(jobs), 'source_and_candidate_native_exclusions_checked': True,
        'physical_unit_solves_unit_response_and_variance_identities': True, 'eligibility_same_for_all_estimators_at_fixed_sampling': True,
        'references_recomputed_from_archived_trials': True, 'thresholds_frozen_before_positive_filtering': True,
        'source_script_threshold_fingerprints_unchanged': True}
    write_json(args.output/'verification.json', verification)
    plots(summary, args.output)
    write_json(args.output/'complete.json', {'new_reductions': 0, 'positive_images': len(jobs),
        'psd_policies': 32, 'control_policies': 16, 'purpose': 'development reuse', 'verification': verification})
    print('completed PSD recovery comparison', flush=True)


if __name__ == '__main__':
    main()
