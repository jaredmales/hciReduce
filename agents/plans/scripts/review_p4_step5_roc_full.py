#!/usr/bin/env python3
"""Verify and summarize the completed full ROC study without changing its design.

Reconstruct native measurements with generic covariance solves, check every
calibration threshold and decision, and separate raw photometry from increments.
The covariance estimators and interpolation geometry reuse the tested study code;
this audit independently checks the masks, solve, reference pixels, and summaries.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import shutil

import numpy as np
from astropy.io import fits

import run_p4_step5_roc_full as full
from run_p4_step5_full_injections import fingerprint, write_json


def require(condition: bool, message: str) -> None:
    """Reject an inconsistent result even when Python assertions are disabled."""
    if not condition:
        raise RuntimeError(message)


def review(root: Path, output: Path) -> None:
    """Audit all saved searches and archive descriptive results with unchanged thresholds."""
    read = full.read
    protocol, results = read(root/'protocol.json'), read(root/'results.json')
    jobs, thresholds = read(root/'jobs.json'), read(root/'thresholds.json')
    calibration, nulls = read(root/'calibration.json'), read(root/'baseline_nulls.json')
    require(read(root/'state.json')['status'] == 'complete', 'study is unfinished')
    require(len(jobs) == len(results['measurements']) == 90 and len(nulls) == 30, 'incomplete sample')
    require(results['nulls'] == nulls and results['thresholds'] == thresholds, 'summary inputs changed')
    records = [*read(root/'complete.json')['products'], *read(root/'calibration_complete.json')['products']]
    for directory in ('reductions', 'references'):
        receipts = list((root/directory).glob('*/complete.json'))
        require(len(receipts) == 91, 'incomplete '+directory)
        for receipt in receipts:
            records += read(receipt)['products']
    unique = {record['path']: record for record in records}
    full.verify(list(unique.values()))
    templates = full.load_templates(root/'payload/response')
    baseline = fits.getdata(root/'reductions/baseline/finim.fits').squeeze().astype(float)
    baseline_maps = full.gaussian.maps(root/'references/baseline')
    counts = {'searches': 0, 'covariance_pixel_solves': 0, 'identity_pixels': 0,
              'reference_searches': 0, 'thresholds': 0, 'raw_photometry': 0}
    maximum = {'amplitude_absolute_difference': 0., 'sigma_relative_difference': 0., 'score_absolute_difference': 0.}
    masks = {}
    yy, xx = np.indices(baseline.shape)
    sx, sy, radius = protocol['known_source_circle']
    # A 15-square kernel shifted over an axial search occupies a 17-square
    # footprint with only its four extreme corner pixels absent.
    for site in protocol['sites']:
        mask = (xx-sx)**2+(yy-sy)**2 <= radius**2
        for trial in [*protocol['calibration_trials'], site]:
            x, y = trial['row'], trial['column']
            mask[y-7:y+8, x-8:x+9] = True
            mask[y-8:y+9, x-7:x+8] = True
        require(np.array_equal(mask, full.holdout_mask(baseline.shape, protocol, site)), 'holdout mask mismatch')
        masks[site['name']] = mask

    def close(actual: float, expected: float, label: str) -> None:
        """Check reproducibility while allowing rounding differences between CPU libraries."""
        require(np.isclose(actual, expected, rtol=2e-8, atol=1e-12), label)

    def audit(image: np.ndarray, trial: dict, site_name: str, saved: dict, maps: dict) -> None:
        """Recalculate one complete native search without the runner's measurement routine."""
        x, y = trial['row'], trial['column']
        observed = {name: [] for name in (*full.COVARIANCE_MODELS, 'identity')}
        for index, (ox, oy) in enumerate(full.radial.OFFSETS):
            position = (x+ox, y+oy)
            data = image[y+oy-5:y+oy+6, x+ox-5:x+ox+6].ravel()
            template = templates[position]
            models = full.candidate_models(image, position, masks[site_name])
            for name, model in models.items():
                require(model is not None, 'unexpected invalid covariance')
                inverse_template = np.linalg.solve(model['covariance'], template)
                energy = float(template @ inverse_template)
                weight = inverse_template/energy
                amplitude = float(weight @ (data-model['mean']))
                sigma = 1/np.sqrt(energy)
                score = amplitude/sigma
                old = saved[name]['pixels'][index]
                require(old['valid'] and model['samples'] == old['training_samples'], 'changed training support')
                for key, value in (('amplitude', amplitude), ('sigma', sigma), ('score', score)):
                    close(value, old[key], trial['name']+'/'+name+'/'+key)
                close(float(weight @ template), 1., 'unit response')
                close(float(weight @ model['covariance'] @ weight)/sigma**2, 1., 'conditional variance')
                maximum['amplitude_absolute_difference'] = max(maximum['amplitude_absolute_difference'], abs(amplitude-old['amplitude']))
                maximum['sigma_relative_difference'] = max(maximum['sigma_relative_difference'], abs(sigma/old['sigma']-1))
                maximum['score_absolute_difference'] = max(maximum['score_absolute_difference'], abs(score-old['score']))
                observed[name].append(score)
                counts['covariance_pixel_solves'] += 1
            energy = float(template @ template)
            amplitude, sigma = float(template @ data/energy), 1/np.sqrt(energy)
            for key, value in (('amplitude', amplitude), ('sigma', sigma), ('score', amplitude/sigma)):
                close(value, saved['identity']['pixels'][index][key], 'identity '+key)
            observed['identity'].append(amplitude/sigma)
            counts['identity_pixels'] += 1
        for name, scores in observed.items():
            require(saved[name]['valid'] and saved[name]['center'] == saved[name]['pixels'][0], 'invalid center/search')
            close(max(scores), saved[name]['search_score'], 'search maximum')
        for name, image_map in maps.items():
            values = [float(image_map[y+oy, x+ox]) for ox, oy in full.radial.OFFSETS]
            close(max(values), saved[name]['search_score'], 'reference search')
            close(values[0], saved[name]['center_value'], 'reference center')
            peak = full.radial.OFFSETS[int(np.argmax(values))]
            require(saved[name]['peak_row_column'] == [x+peak[0], y+peak[1]], 'reference peak position')
            counts['reference_searches'] += 1
        counts['searches'] += 1

    for site_record in calibration:
        site_name = site_record['site']
        rows = site_record['calibration_trials']
        require([r['trial'] for r in rows] == protocol['calibration_trials'], 'calibration set changed')
        for row in rows:
            audit(baseline, row['trial'], site_name, row['models'], baseline_maps)
        for name in full.MODELS:
            require(max(r['models'][name]['search_score'] for r in rows) == thresholds[site_name][name], 'threshold changed')
            counts['thresholds'] += 1
        site = next(s for s in protocol['sites'] if s['name'] == site_name)
        model = full.candidate_models(baseline, (site['row'], site['column']), masks[site_name])['isotropic_b5']
        template = templates[(site['row'], site['column'])]
        sigma = 1/np.sqrt(template @ np.linalg.solve(model['covariance'], template))
        scale = thresholds[site_name]['isotropic_b5']*sigma
        close(scale, site_record['reference_contrast_scale'], 'contrast reference')
        require(not site_record['target_values_used_for_contrast'], 'source brightness uses target values')
        for job in (j for j in jobs if j['site'] == site_name):
            close(job['contrast'], job['brightness_multiplier']*scale, 'source contrast')
        print('verified calibration '+site_name, flush=True)
    for row in nulls:
        audit(baseline, row['site'], row['site']['name'], row['models'], baseline_maps)
    for row, job in zip(results['measurements'], jobs):
        require(row['trial'] == job == read(root/'reductions'/job['name']/'complete.json')['trial'], 'injection job mismatch')
        image, header = fits.getdata(root/'reductions'/job['name']/'finim.fits', header=True)
        image = image.squeeze().astype(float)
        require(image.shape == (256, 256) and np.array_equal(np.isfinite(image), np.isfinite(baseline)), 'image support changed')
        require(header['P4 LOCAL STAMP SIZE'] == 0 and header['COMBINATION METHOD'].strip() == 'mean', 'reduction contract changed')
        audit(image, job, job['site'], row['models'], full.gaussian.maps(root/'references'/job['name']))
        before = next(n for n in nulls if n['site']['name'] == job['site'])
        for name in (*full.COVARIANCE_MODELS, 'identity'):
            value = row['models'][name]
            center = value['center']
            close(center['amplitude']/job['contrast']-1, value['raw_contrast_error'], 'raw error')
            close((center['amplitude']-before['models'][name]['center']['amplitude'])/job['contrast']-1,
                  value['paired_increment_error'], 'paired increment')
            require(value['inside_one_conditional_sigma'] == (abs(center['amplitude']-job['contrast']) <= center['sigma']), 'conditional inclusion')
            counts['raw_photometry'] += 1

    summary = {}
    for name in full.MODELS:
        for row in nulls:
            value = row['models'][name]
            require(value['detected'] == (value['search_score'] > thresholds[row['site']['name']][name]), 'null decision')
        for row in results['measurements']:
            value = row['models'][name]
            require(value['detected'] == (value['search_score'] > thresholds[row['trial']['site']][name]), 'positive decision')
        entry = {'null_exceedances': sum(n['models'][name]['detected'] for n in nulls), 'by_level': []}
        require(entry['null_exceedances'] == results['models'][name]['evaluation_exceedances'], 'null summary')
        values = [r['models'][name] for r in results['measurements']]
        for level in full.LEVELS:
            selected = [r for r in results['measurements'] if r['trial']['brightness_multiplier'] == level]
            clean = [r for r in selected if not next(n['models'][name]['detected'] for n in nulls if n['site']['name'] == r['trial']['site'])]
            group = next(g for g in results['groups'] if g['model'] == name and g['brightness_multiplier'] == level)
            detections = sum(r['models'][name]['detected'] for r in selected)
            require(detections == group['detections'] and group['trials'] == 30 and group['invalid_searches'] == 0, 'recovery summary')
            entry['by_level'].append({'level': level, 'detections': detections,
                'baseline_below_threshold': len(clean), 'new_crossings': sum(r['models'][name]['detected'] for r in clean),
                'by_radius': {str(radius): sum(r['models'][name]['detected'] for r in selected if r['trial']['nominal_radius'] == radius)
                              for radius in protocol['nominal_radii']}})
        if name in (*full.COVARIANCE_MODELS, 'identity'):
            errors = np.array([v['raw_contrast_error'] for v in values])
            increments = [v['paired_increment_error'] for v in values]
            null_scores = np.array([n['models'][name]['center']['score'] for n in nulls])
            entry['photometry'] = {'median_error': float(np.median(errors)), 'median_absolute_error': float(np.median(abs(errors))),
                'rms_error': float(np.sqrt(np.mean(errors**2))), 'within_one_conditional_sigma': sum(v['inside_one_conditional_sigma'] for v in values),
                'median_paired_increment_error': float(np.median(increments)), 'null_centers_within_one_sigma': int(np.sum(abs(null_scores) <= 1)),
                'null_center_score_mean': float(np.mean(null_scores)), 'null_center_score_sample_variance': float(np.var(null_scores, ddof=1))}
        summary[name] = entry
    paired = {}
    for name in ('gaussian_snr', 'identity', 'identity_snr', 'pca_b5_f1'):
        paired[name] = []
        for level in full.LEVELS:
            selected = [r for r in results['measurements'] if r['trial']['brightness_multiplier'] == level]
            paired[name].append({'level': level,
                'hann_only': [r['trial']['site'] for r in selected if r['models'][full.MODELS[0]]['detected'] and not r['models'][name]['detected']],
                'reference_only': [r['trial']['site'] for r in selected if r['models'][name]['detected'] and not r['models'][full.MODELS[0]]['detected']]})
    output.mkdir(parents=True, exist_ok=True)
    write_json(output/'review.json', {'models': summary, 'hann_paired_decisions': paired,
        'verification': {'all_passed': True, 'verified_local_product_hashes': len(unique), 'counts': counts, 'maximum_differences': maximum,
                         'independent_masks': len(masks), 'covariance_estimator_and_geometry_reused': True},
        'inputs': [fingerprint(root/name) for name in ('results.json', 'protocol.json', 'jobs.json', 'thresholds.json', 'calibration.json')]})
    for name in ('results.json', 'results.md', 'comparison.png', 'baseline_nulls.json', 'calibration.json',
                 'complete.json', 'completion_verification.json', 'state.json'):
        require(not (output/name).exists(), 'refuse to overwrite archived '+name)
        shutil.copy2(root/name, output/name)
    print(counts, maximum, flush=True)


def main() -> None:
    """Review an existing completed run and archive its scientific results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    review(args.root.resolve(), args.output.resolve())


if __name__ == '__main__':
    main()
