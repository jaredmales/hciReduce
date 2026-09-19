#!/usr/bin/env python3
"""Verify post-mean patch normalization, archived controls, maps, and calibrated decisions."""
from pathlib import Path
import sys
import time

import numpy as np
from astropy.io import fits

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO/'agents/plans/scripts'))
from run_p4_step5_roc_full import read, verify
from run_p4_step5_full_injections import fingerprint, write_json
from compare_p4_step5_radial_pooling import OFFSETS
import compare_p4_step5_patch_rms as test
from setup_check import direct_covariance

ROOT = REPO/'working/roc/p4_patch_rms_step5_20260919'
REPORT = Path(__file__).resolve().parent
protocol, result = read(ROOT/'protocol.json'), read(ROOT/'results.json')
methods = protocol['methods']
# Preserve the parent search ordering from its production-independent archived driver.
offsets = OFFSETS
start = time.monotonic()
verified = []
verify(read(ROOT/'complete.json')['products'])
verify(read(ROOT/'calibration_complete.json')['products'])
manifest = read(ROOT/'manifest.json')
missing = []
for record in [*manifest['inputs'], manifest['protocol']]:
    if Path(record['path']).exists():
        verify([record])
        verified.append(record)
    else:
        missing.append(record['path'])

settings = read(Path(protocol['study'])/'protocol.json')
yy, xx = np.indices((256, 256))
radius64 = np.hypot(xx-127.5, yy-127.5)
radius = radius64.astype('f4')
sx, sy, sr = settings['known_source_circle']
noise_mask = np.hypot(xx-sx, yy-sy).astype('f4') > sr+.5
counts = {'tasks': 0, 'searches': 0, 'snr_pixels': 0, 'profiles': 0, 'thresholds': 0, 'decisions': 0, 'mean_toggle_pixels': 0}
maximum_snr_difference = 0.
maximum_generic_difference = 0.
counts.update(raw_control_planes=0, generic_patch_covariance_solves=0)
templates = test.full.load_templates(Path(protocol['study'])/'payload/response')
prior_root = Path(protocol['previous_ablation'])
records = {'baseline': [], 'positive': []}
for phase in records:
    for file in sorted((ROOT/phase).glob('*/measurements.json')):
        directory = file.parent
        receipt = read(directory/'complete.json')
        verify(receipt['products'])
        verified.extend(receipt['products'])
        record = read(file)
        records[phase].append(record)
        amplitude = fits.getdata(directory/'amplitudes.fits')
        snr = fits.getdata(directory/'amplitudes_snr.fits')
        assert amplitude.shape == snr.shape == (10, 256, 256)
        for filename, actual in (('amplitudes.fits', amplitude), ('amplitudes_snr.fits', snr)):
            original = fits.getdata(prior_root/phase/directory.name/filename)
            for name in methods:
                if name in test.previous.METHODS:
                    assert np.array_equal(actual[methods.index(name)], original[test.previous.METHODS.index(name)], equal_nan=True)
                    counts['raw_control_planes'] += 1
        saved_profiles = read(directory/'profiles.json')
        for m, name in enumerate(methods):
            profile = saved_profiles[name]
            for p in profile:
                lo = p['radius']-.5
                values = amplitude[m][noise_mask & np.isfinite(amplitude[m]) & (radius > lo) & (radius <= lo+1)].astype(float)
                assert p['pixels'] == len(values)
                assert np.isclose(p['mean'], np.mean(values), rtol=1e-12, atol=1e-15)
                assert np.isclose(p['stddev'], np.std(values, ddof=1), rtol=1e-12, atol=1e-15)
                counts['profiles'] += 1
            for row in record['rows']:
                trial, value = row['trial'], row['models'][name]
                x = np.array([trial['row']+dx for dx, dy in offsets])
                y = np.array([trial['column']+dy for dx, dy in offsets])
                assert value['valid'] and np.isfinite(amplitude[m, y, x]).all()
                assert value['snr_pixels'] == snr[m, y, x].tolist()
                assert value['amplitude_pixels'] == amplitude[m, y, x].tolist()
                assert value['search_score'] == float(np.max(snr[m, y, x]))
                assert value['center_amplitude'] == float(amplitude[m, trial['column'], trial['row']])
                avg = np.interp(radius[y, x], [p['radius'] for p in profile], [p['mean'] for p in profile])
                std = np.interp(radius[y, x], [p['radius'] for p in profile], [p['stddev'] for p in profile]).astype('f4')
                n = 2*np.pi*radius64[y, x]/float(np.float32(3.6))-1
                expected = ((amplitude[m, y, x]-avg)/std).astype('f4') * (1/np.sqrt(1+1/n)).astype('f4')
                difference = float(np.max(abs(expected-snr[m, y, x])))
                maximum_snr_difference = max(maximum_snr_difference, difference)
                assert np.allclose(expected, snr[m, y, x], rtol=3e-6, atol=3e-6)
                counts['snr_pixels'] += len(x)
                counts['searches'] += 1
        for row in record['rows']:
            for pixel in row['fit_pixels']:
                for prefix, diagnostic in (('hann_psd', 'hann'), ('rect_psd', 'rect'),
                                           ('hann_patch_rms', 'hann_patch_rms'), ('rect_patch_rms', 'rect_patch_rms')):
                    a, d = pixel['amplitudes'], pixel['diagnostic'][diagnostic]
                    assert np.isclose(a[prefix+'_no_mean']-a[prefix+'_mean'], d['mean_projection'], rtol=1e-9, atol=1e-13)
                    counts['mean_toggle_pixels'] += 1
        site = record['task']['site']
        x, y = site['row'], site['column']
        science = fits.getdata(Path(protocol['study'])/'reductions'/record['task']['image']/'finim.fits').squeeze().astype(float)
        matrices = test.previous.narrow_samples(science, (x, y), test.full.holdout_mask(science.shape, settings, site))
        template = templates[(x, y)]
        for prefix, window, mixing, width, _ in test.FAMILIES:
            samples = matrices[0] if width == 0 else np.vstack(list(matrices.values()))
            covariance = direct_covariance(samples, window, mixing)
            q = np.linalg.solve(covariance, template)
            weight = q/(template @ q)
            data = science[y-5:y+6, x-5:x+6].ravel()
            on = float(weight @ (data-samples.mean(axis=0)))
            off = float(weight @ data)
            reference = record['rows'][-1]['fit_pixels'][0]
            for kind, expected in (('mean', on), ('no_mean', off)):
                actual = reference['amplitudes'][prefix+'_patch_rms_'+kind]
                assert np.isclose(actual, expected, rtol=1e-10, atol=1e-13)
                maximum_generic_difference = max(maximum_generic_difference, abs(actual-expected))
            assert np.isclose(reference['diagnostic'][prefix+'_patch_rms']['sigma'], 1/np.sqrt(template @ q), rtol=1e-10, atol=1e-13)
            counts['generic_patch_covariance_solves'] += 1
        counts['tasks'] += 1
assert [len(records[p]) for p in records] == [30, 90]
thresholds = read(ROOT/'thresholds.json')
for record in records['baseline']:
    site = record['task']['site']['name']
    assert len(record['rows']) == 29
    for name in methods:
        assert thresholds[site][name] == max(r['models'][name]['search_score'] for r in record['rows'][:-1])
        counts['thresholds'] += 1
assert result['thresholds'] == thresholds
for phase, key in (('baseline', 'nulls'), ('positive', 'injections')):
    saved = {r['job']['name'] if phase == 'positive' else r['trial']['name']: r for r in result[key]}
    for record in records[phase]:
        row = record['rows'][-1]
        target = saved[record['task']['name']]
        for name in methods:
            actual = row['models'][name]
            assert target['models'][name] == {**actual, 'detected': actual['search_score'] > thresholds[row['trial']['name']][name]}
            counts['decisions'] += 1

for name in methods:
    summary = result['models'][name]
    assert summary['null_exceedances'] == sum(r['models'][name]['detected'] for r in result['nulls'])
    assert summary['invalid_nulls'] == 0
    for group in summary['levels']:
        rows = [r for r in result['injections'] if r['job']['brightness_multiplier'] == group['level']]
        assert group['trials'] == len(rows) == 30 and group['invalid'] == 0
        assert group['detections'] == sum(r['models'][name]['detected'] for r in rows)
        if name != 'gaussian':
            errors = np.array([r['models'][name]['center_amplitude']/r['job']['contrast']-1 for r in rows])
            assert np.isclose(group['median_raw_error'], np.median(errors), rtol=1e-12)
            assert np.isclose(group['rms_raw_error'], np.sqrt(np.mean(errors**2)), rtol=1e-12)
for key, levels in result['paired_decisions'].items():
    tested, reference = key.split('/')
    for group in levels:
        rows = [r for r in result['injections'] if r['job']['brightness_multiplier'] == group['level']]
        assert group['test_only'] == [r['trial']['name'] for r in rows if r['models'][tested]['detected'] and not r['models'][reference]['detected']]
        assert group['reference_only'] == [r['trial']['name'] for r in rows if r['models'][reference]['detected'] and not r['models'][tested]['detected']]
derived = {'photometry': {}, 'mean_toggles': {}, 'by_radius': {}, 'normalization_effect': {}}
for name in methods:
    if name != 'gaussian':
        errors = np.array([r['models'][name]['center_amplitude']/r['job']['contrast']-1 for r in result['injections']])
        derived['photometry'][name] = {'median_raw_error': float(np.median(errors)),
                                     'rms_raw_error': float(np.sqrt(np.mean(errors**2)))}
    derived['by_radius'][name] = [{'radius': radius, 'trials': 18,
        'detections': sum(r['models'][name]['detected'] for r in result['injections'] if r['trial']['nominal_radius'] == radius)}
        for radius in (26, 32, 38, 44, 50)]
for family in ('hann', 'rect'):
    off, on = family+'_patch_rms_no_mean', family+'_patch_rms_mean'
    derived['mean_toggles'][family] = {}
    for key in ('nulls', 'injections'):
        rows = result[key]
        difference = np.array([r['models'][on]['search_score']-r['models'][off]['search_score'] for r in rows])
        derived['mean_toggles'][family][key] = {'mean_absolute_score_change': float(np.mean(abs(difference))),
            'maximum_absolute_score_change': float(np.max(abs(difference))),
            'flips': [r['trial']['name'] for r in rows if r['models'][on]['detected'] != r['models'][off]['detected']]}
    for subtract in ('mean', 'no_mean'):
        tested, raw = family+'_patch_rms_'+subtract, family+'_psd_'+subtract
        difference = np.array([r['models'][tested]['search_score']-r['models'][raw]['search_score'] for r in result['injections']])
        ratios = [p['diagnostic'][family+'_patch_rms']['rms_max']/p['diagnostic'][family+'_patch_rms']['rms_min']
                  for row in result['injections'] for p in row['fit_pixels']]
        derived['normalization_effect'][tested] = {'mean_absolute_positive_score_change': float(np.mean(abs(difference))),
            'maximum_absolute_positive_score_change': float(np.max(abs(difference))),
            'median_training_max_min_rms_ratio': float(np.median(ratios)), 'maximum_training_max_min_rms_ratio': float(np.max(ratios)),
            'null_decision_flips': [r['trial']['name'] for r in result['nulls']
                                    if r['models'][tested]['detected'] != r['models'][raw]['detected']]}
write_json(REPORT/'derived_comparisons.json', derived)
write_json(REPORT/'review.json', {'all_passed': True, 'checks': counts,
    'verified_distinct_local_fingerprints': len({r['path'] for r in verified}),
    'remote_only_dependencies': missing, 'maximum_search_snr_difference': maximum_snr_difference,
    'maximum_generic_amplitude_difference': maximum_generic_difference,
    'elapsed_seconds': time.monotonic()-start, 'driver': fingerprint(REPO/'agents/plans/scripts/compare_p4_step5_patch_rms.py'),
    'results': fingerprint(ROOT/'results.json'), 'new_reductions': 0})
print(counts)
