#!/usr/bin/env python3
"""Test equal-RMS training patches after ensemble mean subtraction on saved injections.

Retain the fixed raw PSD controls and mean-on/off pairs, sampling, source masks,
five-pixel searches, and production annular SNR. Freeze baseline thresholds before
positives. Candidate data and response retain their original contrast units.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
import math
import multiprocessing as mp
import os
from pathlib import Path
import subprocess
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.linalg import cho_factor, cho_solve, eigvalsh

import run_p4_step5_roc_full as full
import compare_p4_step5_mean_snr as previous
from measure_p4_step5_planet import annular_oracle
from run_p4_step5_full_injections import fingerprint, write_json

METHODS = ('gaussian', 'identity', 'hann_psd_no_mean', 'hann_psd_mean',
           'hann_patch_rms_no_mean', 'hann_patch_rms_mean', 'rect_psd_no_mean', 'rect_psd_mean',
           'rect_patch_rms_no_mean', 'rect_patch_rms_mean')
FAMILIES = previous.FAMILIES
LABELS = ('Gaussian 3.6', 'Identity', 'Hann raw, no mean', 'Hann raw + mean',
          'Hann patch RMS, no mean', 'Hann patch RMS + mean', 'Rect. raw, no mean', 'Rect. raw + mean',
          'Rect. patch RMS, no mean', 'Rect. patch RMS + mean')
_CONTEXT = None


def fit_patch_rms(samples: np.ndarray, window_name: str, mixing: float) -> dict | None:
    """Average equal-RMS residual-patch powers without a second ensemble centering."""
    full.radial.require(samples.ndim == 2 and samples.shape[1] == 121 and np.isfinite(samples).all(),
                        'expected finite 11x11 training patches')
    full.radial.require(window_name in ('rectangular', 'hann') and 0 < mixing <= 1, 'invalid fixed PSD settings')
    if len(samples) < 8:
        return None
    mean = samples.mean(axis=0)
    centered = samples-mean
    rms = np.sqrt(np.mean(centered**2, axis=1))
    roundoff = 16*np.finfo(float).eps*float(np.max(abs(samples)))
    if not np.isfinite(rms).all() or np.any(rms <= roundoff):
        return None
    normalized = centered/rms[:, None]
    full.radial.require(np.allclose(np.mean(normalized**2, axis=1), 1, rtol=1e-12), 'patch RMS normalization failed')
    target = float(np.sum(centered**2)/((len(samples)-1)*121))
    window = np.ones((11, 11)) if window_name == 'rectangular' else np.outer(np.hanning(11), np.hanning(11))
    power = np.sum(abs(np.fft.fft2(normalized.reshape((-1, 11, 11))*window,
                                  s=(21, 21), axes=(-2, -1)))**2, axis=0)/((len(samples)-1)*np.sum(window**2))
    # Keep the old physical variance scale; only the relative patch contribution changes.
    power *= target/float(power.mean())
    power = (1-mixing)*power+mixing*target
    lag = np.fft.ifft2(power)
    full.radial.require(np.max(abs(lag.imag)) <= 1e-12*target, 'complex lag covariance')
    yy, xx = np.indices((11, 11))
    covariance = lag.real[(yy.ravel()[:, None]-yy.ravel()[None, :]) % 21,
                          (xx.ravel()[:, None]-xx.ravel()[None, :]) % 21]
    covariance = (covariance+covariance.T)/2
    eigenvalues = eigvalsh(covariance, check_finite=False)
    full.radial.require(np.isclose(np.trace(covariance), 121*target, rtol=1e-12) and
                        eigenvalues[0] >= mixing*target*(1-1e-10), 'variance or positive-floor identity failed')
    return {'mean': mean, 'covariance': covariance, 'factorization': cho_factor(covariance, lower=True, check_finite=False),
        'target_variance': target, 'power': power, 'rms': rms, 'normalized_mean': normalized.mean(axis=0)}


def filter_pixel(science: np.ndarray, position: tuple, template: np.ndarray, forbidden: np.ndarray) -> tuple:
    """Change only the training residual-patch power weighting within each fixed PSD family."""
    x, y = position
    data = science[y-5:y+6, x-5:x+6].ravel()
    matrices = previous.narrow_samples(science, position, forbidden)
    amplitudes = {'identity': float(template @ data/(template @ template))}
    diagnostics = {}
    for prefix, window, mixing, width, _ in FAMILIES:
        samples = matrices[0] if width == 0 else np.vstack(list(matrices.values()))
        models = {'psd': full.psd.fit_psd(samples, window, mixing), 'patch_rms': fit_patch_rms(samples, window, mixing)}
        if all(m is not None for m in models.values()):
            full.radial.require(np.array_equal(models['psd']['mean'], models['patch_rms']['mean']) and
                                np.isclose(models['psd']['target_variance'], models['patch_rms']['target_variance'], rtol=1e-12),
                                'normalization changed raw fitted mean or variance scale')
        for kind, model in models.items():
            key = prefix if kind == 'psd' else prefix+'_patch_rms'
            diagnostics[key] = {'valid': model is not None, 'samples': len(samples)}
            if model is None:
                continue
            q = cho_solve(model['factorization'], template, check_finite=False)
            energy = float(template @ q)
            weight = q/energy
            off, on = float(weight @ data), float(weight @ (data-model['mean']))
            full.radial.require(np.isclose(weight @ template, 1, rtol=1e-12) and
                                np.isclose(off-on, weight @ model['mean'], rtol=1e-9, atol=1e-13), 'filter normalization failed')
            amplitudes.update({prefix+'_'+kind+'_no_mean': off, prefix+'_'+kind+'_mean': on})
            diagnostics[key].update(amplitude=on, sigma=1/math.sqrt(energy), score=on*math.sqrt(energy),
                                    mean_projection=float(weight @ model['mean']))
            if kind == 'patch_rms':
                diagnostics[key].update(rms_min=float(model['rms'].min()), rms_median=float(np.median(model['rms'])),
                    rms_max=float(model['rms'].max()), normalized_mean_rms=float(np.sqrt(np.mean(model['normalized_mean']**2))))
    return amplitudes, diagnostics


def required_bins(trials: list, radius: np.ndarray) -> tuple:
    """Select complete one-pixel annuli bracketing every native search-pixel radius."""
    bins, positions = set(), set()
    for trial in trials:
        for ox, oy in full.radial.OFFSETS:
            x, y = trial['row']+ox, trial['column']+oy
            positions.add((x, y))
            lower = math.floor(float(radius[y, x])-.5)
            bins.update((lower, lower+1))
    selected = np.zeros(radius.shape, bool)
    for lower in bins:
        selected |= (radius > lower) & (radius <= lower+1)
    return sorted(bins), positions, selected


def prepare(root: Path, study: Path) -> None:
    """Freeze the post-mean patch normalization and its raw controls before new measurements."""
    full.radial.require(not (root/'manifest.json').exists(), 'comparison already prepared')
    parent = full.read(study/'protocol.json')
    full.radial.require(full.read(study/'state.json')['status'] == 'complete', 'parent study incomplete')
    records = [fingerprint(study/name) for name in ('protocol.json', 'jobs.json', 'results.json', 'calibration.json',
                'baseline_nulls.json', 'thresholds.json', 'complete.json', 'manifest.json')]
    records += [fingerprint(p) for p in sorted((study/'payload/response').glob('*.fits'))]
    for directory in ('reductions', 'references'):
        receipts = sorted((study/directory).glob('*/complete.json'))
        full.radial.require(len(receipts) == 91, 'incomplete saved image set')
        for path in receipts:
            records += [fingerprint(path), *full.read(path)['products']]
    native = [r for r in full.read(study/'manifest.json')['frozen_records']
              if Path(r['path']).parent == study/'software' and Path(r['path']).suffix != '.py']
    ablation = study.parent/'p4_mean_snr_step5_20260919'
    full.radial.require(full.read(ablation/'state.json')['status'] == 'complete', 'previous annular comparison incomplete')
    records += [fingerprint(ablation/name) for name in ('results.json', 'thresholds.json', 'complete.json', 'manifest.json', 'protocol.json')]
    for phase in ('baseline', 'positive'):
        for receipt in sorted((ablation/phase).glob('*/complete.json')):
            records += [fingerprint(receipt), *full.read(receipt)['products']]
    scripts = [fingerprint(p) for p in sorted(Path(__file__).parent.glob('*.py'))]
    records = list({r['path']: r for r in [*records, *native, *scripts]}.values())
    full.verify(records)
    protocol = {'schema': 1, 'purpose': 'development comparison: individual residual-patch RMS after ensemble mean subtraction',
        'study': str(study), 'previous_ablation': str(ablation), 'methods': list(METHODS), 'families': FAMILIES, 'new_reductions': 0,
        'sites': parent['sites'], 'positive_images': 90, 'brightnesses': list(full.LEVELS),
        'search': 'unchanged native center plus four axial one-pixel neighbors; all five required; invalid search is nondetection',
        'training': 'unchanged sitewise parent holdout and complete finite 11x11 raw patches; offsets -5,0,+5; step5; minimum8',
        'mean_toggle': 'within each raw/RMS condition hold covariance and weights fixed; toggle only candidate subtraction of the raw fitted mean',
        'patch_normalization': 'subtract raw ensemble mean patch, divide each training residual by sqrt(mean(residual**2)); no second centering; no spatial scalar-mean subtraction',
        'variance_scale': 'rescale PSD mean to original raw across-patch mean pixel variance before the same spectral mixing',
        'candidate_units': 'candidate and template remain raw; candidate RMS is never fitted; raw training mean is identical in each raw/RMS pair',
        'invalid_normalization': 'any nonfinite RMS or RMS <= 16*float64 epsilon*max(abs(raw samples)) invalidates the normalized fit; no sample removal or tuning',
        'annular_snr': 'production hciAnalyze mean/stddev profile in one-pixel bins, interpolated to radius; small-sample correction; lambdaD3.6; range0..60',
        'annular_exclusion': 'unchanged legacy known source only: planet.R7.3 plus production0.5 pixel buffer; trial neighborhoods remain in radial normalization',
        'threshold': 'each site/method maximum of original28 calibration search scores; strict exceedance; all thresholds frozen before positive analysis',
        'radial_computation': 'complete radial bins required by search interpolation; baseline covers calibration+site; positives cover site; sparse maps elsewhere',
        'validation_aperture': 'CLI apertureR60 only avoids invalid known-planet aperture in sparse maps; does not change SNR maps; experiment searches remain five pixels',
        'workers': 12, 'cpu_ids': list(range(12)), 'threads_per_worker': 1,
        'dependence': 'same correlated residual field and already inspected injections; this is not fresh validation'}
    write_json(root/'protocol.json', protocol)
    write_json(root/'manifest.json', {'inputs': records, 'protocol': fingerprint(root/'protocol.json'), 'host': os.uname().nodename})
    write_json(root/'state.json', {'status': 'prepared'})


def initialize(root: str, queue) -> None:
    """Assign one physical CPU per process and load immutable shared study context."""
    global _CONTEXT
    os.sched_setaffinity(0, {queue.get()})
    path = Path(root)
    contract = full.read(path/'protocol.json')
    study = Path(contract['study'])
    parent = full.read(study/'protocol.json')
    yy, xx = np.indices((256, 256))
    radius = np.hypot(xx-127.5, yy-127.5).astype('f4')
    calibration = {r['site']: r['calibration_trials'] for r in full.read(study/'calibration.json')}
    nulls = {r['site']['name']: r['models'] for r in full.read(study/'baseline_nulls.json')}
    positives = {r['trial']['name']: r['models'] for r in full.read(study/'results.json')['measurements']}
    _CONTEXT = (path, study, parent, radius, full.load_templates(study/'payload/response'), calibration, nulls, positives, Path(contract['previous_ablation']))


def measure(task: dict) -> str:
    """Produce one site's baseline or positive amplitude/SNR cubes and audited search records."""
    root, study, parent, radius, templates, old_calibration, old_nulls, old_positives, ablation = _CONTEXT
    site, image_name, phase = task['site'], task['image'], task['phase']
    directory = root/phase/task['name']
    if (directory/'complete.json').exists():
        full.verify(full.read(directory/'complete.json')['products'])
        return str(directory/'measurements.json')
    directory.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    source = study/'reductions'/image_name/'finim.fits'
    full.verify(full.read(source.parent/'complete.json')['products'])
    science, header = fits.getdata(source, header=True)
    science = science.squeeze().astype(float)
    baseline_image = fits.getdata(study/'reductions/baseline/finim.fits').squeeze()
    full.radial.require(np.array_equal(np.isfinite(science), np.isfinite(baseline_image)), 'changed science support')
    trials = [*parent['calibration_trials'], site] if phase == 'baseline' else [site]
    bins, positions, selected = required_bins(trials, radius)
    forbidden = full.holdout_mask(science.shape, parent, site)
    maps = np.full((len(METHODS), 256, 256), np.nan, dtype='f4')
    gaussian = fits.getdata(study/'references'/image_name/'gaussian.fits').squeeze()
    maps[0][selected] = gaussian[selected]
    details = {}
    fitted_pixels, fitted_models = 0, 0
    for (x, y), template in templates.items():
        if not selected[y, x]:
            continue
        data = science[y-5:y+6, x-5:x+6]
        if data.size != 121 or not np.isfinite(data).all():
            continue
        amplitudes, diagnostic = filter_pixel(science, (x, y), template, forbidden)
        for name, amplitude in amplitudes.items():
            maps[METHODS.index(name), y, x] = amplitude
        fitted_pixels += 1
        fitted_models += sum(d['valid'] for d in diagnostic.values())
        if (x, y) in positions:
            details[(x, y)] = {'amplitudes': amplitudes, 'diagnostic': diagnostic}
    full.radial.require(all(p in details for p in positions), 'missing preassigned search support')
    header['HCI FILTER LABELS'] = ','.join(METHODS)
    header['HCI RADIAL BINS'] = ','.join(map(str, bins))
    fits.writeto(directory/'amplitudes.fits', maps, header)
    command = [str(study/'software/hciAnalyze'), '--file='+str(directory/'amplitudes.fits'), '--lambdaD=3.6',
        '--planet.sep=11.782', '--planet.PA=262.051', '--planet.R=7.3', '--snr.apertureR=60',
        '--snr.minRad=0', '--snr.maxRad=60', '--filter.psfResponse=', '--filter.lpfGaussFW=0', '--filter.hpfGaussFW=0',
        '--noise.model=identity', '--noise.only=false', '--noise.outputDiagnostics=false']
    write_json(directory/'command.json', command)
    environment = full.binary_environment(study, reduction=False)
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (directory/'analysis.log').open('w') as log:
        subprocess.run(command, env=environment, cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)
    snr, snr_header = fits.getdata(directory/'amplitudes_snr.fits', header=True)
    full.radial.require(snr.shape == maps.shape and snr_header['SNRMINR'] == 0 and snr_header['SNRMAXR'] == 60 and
                        snr_header['SNRMEAN'] == snr_header['SNRSMALL'] == 1, 'changed production normalization')
    settings = {'source_x': parent['known_source_circle'][0], 'source_y': parent['known_source_circle'][1],
        'source_radius': parent['known_source_circle'][2], 'min_radius': 0, 'max_radius': 60, 'lambda_d': float(np.float32(3.6))}
    profiles, oracle_errors = {}, {}
    for index, name in enumerate(METHODS):
        expected, profile = annular_oracle(maps[index], settings)
        valid = np.isfinite(expected)
        full.radial.require(np.allclose(snr[index][valid], expected[valid], rtol=3e-6, atol=3e-6), 'annular oracle mismatch')
        oracle_errors[name] = float(np.max(abs(snr[index][valid]-expected[valid]))) if valid.any() else None
        profiles[name] = [p for p in profile if int(p['radius']-.5) in bins]
    earlier_maps = {label: fits.getdata(ablation/phase/task['name']/filename) for label, filename in
                    (('amplitude', 'amplitudes.fits'), ('snr', 'amplitudes_snr.fits'))}
    for name in METHODS:
        if name in previous.METHODS:
            for label, current in (('amplitude', maps), ('snr', snr)):
                full.radial.require(np.allclose(current[METHODS.index(name)], earlier_maps[label][previous.METHODS.index(name)],
                                    rtol=1e-6, atol=1e-6 if label == 'snr' else 1e-12, equal_nan=True), 'raw control map changed')
    old_rows = {r['trial']['name']: r['models'] for r in old_calibration[site['name']]}
    old_rows[site['name']] = old_nulls[site['name']]
    old_snr_maps = {name: fits.getdata(study/'references'/image_name/(name+'_snr.fits')).squeeze()
                    for name in ('gaussian', 'identity')}
    scalar_checks, maximum_difference = 0, 0.
    rows = []
    for trial in trials:
        old = old_rows[trial['name']] if phase == 'baseline' else old_positives[image_name]
        pixels = [details[(trial['row']+dx, trial['column']+dy)] for dx, dy in full.radial.OFFSETS]
        for pixel_index, pixel in enumerate(pixels):
            comparisons = [(pixel['amplitudes']['identity'], old['identity']['pixels'][pixel_index]['amplitude'])]
            for prefix, _, _, _, old_name in FAMILIES:
                diagnostic = pixel['diagnostic'][prefix]
                reference = old[old_name]['pixels'][pixel_index]
                full.radial.require(diagnostic['valid'] == reference['valid'], 'old PSD validity changed')
                if diagnostic['valid']:
                    full.radial.require(diagnostic['samples'] == reference['training_samples'], 'training patch count changed')
                    comparisons += [(diagnostic[k], reference[k]) for k in ('amplitude', 'sigma', 'score')]
            for actual, expected in comparisons:
                full.radial.require(np.isclose(actual, expected, rtol=1e-8, atol=1e-12), 'archived candidate fit mismatch')
                scalar_checks += 1
                maximum_difference = max(maximum_difference, abs(actual-expected))
        models = {}
        for index, name in enumerate(METHODS):
            values = [float(snr[index, trial['column']+dy, trial['row']+dx]) for dx, dy in full.radial.OFFSETS]
            valid = all(np.isfinite(maps[index, trial['column']+dy, trial['row']+dx]) for dx, dy in full.radial.OFFSETS)
            full.radial.require(all(np.isfinite(v) for v in values), 'nonfinite SNR search')
            amplitudes = [float(maps[index, trial['column']+dy, trial['row']+dx]) for dx, dy in full.radial.OFFSETS]
            models[name] = {'valid': bool(valid), 'search_score': max(values) if valid else None,
                'center_amplitude': amplitudes[0] if valid else None, 'snr_pixels': values if valid else [],
                'amplitude_pixels': amplitudes if valid else []}
            if name in ('gaussian', 'identity'):
                old_name = 'gaussian_snr' if name == 'gaussian' else 'identity_snr'
                full.radial.require(np.isclose(max(values), old[old_name]['search_score'], rtol=1e-6, atol=1e-6), 'reference search changed')
                original = [float(old_snr_maps[name][trial['column']+dy, trial['row']+dx]) for dx, dy in full.radial.OFFSETS]
                full.radial.require(np.allclose(values, original, rtol=1e-6, atol=1e-6), 'reference search pixels changed')
        rows.append({'trial': trial, 'models': models, 'fit_pixels': pixels})
    write_json(directory/'profiles.json', profiles)
    write_json(directory/'measurements.json', {'task': task, 'rows': rows, 'radial_bins': bins,
        'elapsed_seconds': time.monotonic()-start, 'fitted_pixels': fitted_pixels, 'fitted_psd_models': fitted_models,
        'archived_scalar_checks': scalar_checks, 'maximum_archived_scalar_difference': maximum_difference, 'oracle_errors': oracle_errors})
    products = [fingerprint(p) for p in sorted(directory.iterdir()) if p.is_file()]
    write_json(directory/'complete.json', {'products': products})
    return str(directory/'measurements.json')


def pool_phase(root: Path, tasks: list, phase: str) -> list:
    """Run one phase on separate physical cores and preserve atomic progress for restart."""
    protocol = full.read(root/'protocol.json')
    context = mp.get_context('spawn')
    queue = context.Queue()
    for cpu in protocol['cpu_ids']:
        queue.put(cpu)
    completed = []
    with ProcessPoolExecutor(max_workers=protocol['workers'], mp_context=context, initializer=initialize,
                             initargs=(str(root), queue)) as executor:
        futures = {executor.submit(measure, task): task['name'] for task in tasks}
        for future in as_completed(futures):
            path = future.result()
            completed.append(path)
            write_json(root/'state.json', {'status': phase, 'pid': os.getpid(), 'finished': len(completed), 'total': len(tasks)})
            print(f'{phase}: {len(completed)}/{len(tasks)} {futures[future]}', flush=True)
    return sorted(completed)


def summarize(root: Path, baseline: list, positives: list, thresholds: dict) -> None:
    """Compare calibrated recovery, null exceedances, and paired gains for every fixed ablation."""
    nulls, rows, summary = [], [], {}
    for record in baseline:
        row = record['rows'][-1]
        nulls.append(row)
    for record in positives:
        row = record['rows'][0]
        row['job'] = record['task']['job']
        rows.append(row)
    for row in [*nulls, *rows]:
        for name, value in row['models'].items():
            value['detected'] = value['valid'] and value['search_score'] > thresholds[row['trial']['name']][name]
    parent = full.read(Path(full.read(root/'protocol.json')['study'])/'results.json')
    for name, old_name in (('gaussian', 'gaussian_snr'), ('identity', 'identity_snr')):
        for row in nulls:
            old = next(n for n in parent['nulls'] if n['site']['name'] == row['trial']['name'])
            full.radial.require(row['models'][name]['detected'] == old['models'][old_name]['detected'], 'reference null decision changed')
        for row in rows:
            old = next(r for r in parent['measurements'] if r['trial']['name'] == row['job']['name'])
            full.radial.require(row['models'][name]['detected'] == old['models'][old_name]['detected'], 'reference injection decision changed')
    for name in METHODS:
        groups = []
        for level in full.LEVELS:
            selected = [r for r in rows if r['job']['brightness_multiplier'] == level]
            errors = [r['models'][name]['center_amplitude']/r['job']['contrast']-1 for r in selected
                      if name != 'gaussian' and r['models'][name]['valid']]
            groups.append({'level': level, 'trials': len(selected), 'detections': sum(r['models'][name]['detected'] for r in selected),
                'invalid': sum(not r['models'][name]['valid'] for r in selected),
                'median_raw_error': float(np.median(errors)) if errors else None,
                'rms_raw_error': float(np.sqrt(np.mean(np.array(errors)**2))) if errors else None})
        summary[name] = {'null_exceedances': sum(r['models'][name]['detected'] for r in nulls),
                         'invalid_nulls': sum(not r['models'][name]['valid'] for r in nulls), 'levels': groups}
    pairs = {}
    for prefix, *_ in FAMILIES:
        for subtract in ('no_mean', 'mean'):
            tested = prefix+'_patch_rms_'+subtract
            for reference in (prefix+'_psd_'+subtract, 'identity', 'gaussian'):
                pairs[tested+'/'+reference] = [{'level': level,
                    'test_only': [r['trial']['name'] for r in rows if r['job']['brightness_multiplier'] == level and
                                  r['models'][tested]['detected'] and not r['models'][reference]['detected']],
                    'reference_only': [r['trial']['name'] for r in rows if r['job']['brightness_multiplier'] == level and
                                       r['models'][reference]['detected'] and not r['models'][tested]['detected']]}
                    for level in full.LEVELS]
    prior = full.read(Path(full.read(root/'protocol.json')['previous_ablation'])/'results.json')
    for name in METHODS:
        if name not in previous.METHODS:
            continue
        full.radial.require(summary[name] == prior['models'][name], 'raw control summary changed')
        for row in [*nulls, *rows]:
            key = 'injections' if 'job' in row else 'nulls'
            matching = next(r for r in prior[key] if (r['job']['name'] if key == 'injections' else r['trial']['name']) ==
                            (row['job']['name'] if key == 'injections' else row['trial']['name']))
            full.radial.require(row['models'][name]['detected'] == matching['models'][name]['detected'], 'raw control decision changed')
    verification = {'archived_scalar_checks': sum(r['archived_scalar_checks'] for r in [*baseline, *positives]),
        'fitted_psd_models': sum(r['fitted_psd_models'] for r in [*baseline, *positives]),
        'maximum_archived_scalar_difference': max(r['maximum_archived_scalar_difference'] for r in [*baseline, *positives]),
        'maximum_annular_oracle_difference': max(v for r in [*baseline, *positives] for v in r['oracle_errors'].values() if v is not None)}
    write_json(root/'results.json', {'models': summary, 'paired_decisions': pairs, 'nulls': nulls, 'injections': rows,
        'thresholds': thresholds, 'verification': verification, 'baseline_records': baseline, 'new_reductions': 0,
        'previous_ablation_summary': prior['models'], 'reference_decisions_unchanged': True})
    lines = ['# Per-patch normalization after mean subtraction', '', '| Method | 0.5× /30 | 0.75× /30 | 1× /30 | Nulls /30 |',
             '| --- | --- | --- | --- | --- |']
    for name in METHODS:
        entry = summary[name]
        lines.append('| '+name+' | '+' | '.join(str(g['detections']) for g in entry['levels'])+f' | {entry["null_exceedances"]} |')
    (root/'results.md').write_text('\n'.join(lines)+'\n')
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    for name, label in zip(METHODS, LABELS):
        axes[0].plot(full.LEVELS, [g['detections']/30 for g in summary[name]['levels']], 'o-', label=label)
    axes[0].set(xlabel='Brightness / fixed reference contrast', ylabel='Recovery fraction', ylim=(0, 1))
    axes[0].legend(fontsize=7)
    axes[1].barh(LABELS, [summary[n]['null_exceedances'] for n in METHODS])
    axes[1].set(xlabel='Baseline exceedances out of 30')
    fig.suptitle('Saved 90 injections: equal-RMS training patches after mean subtraction')
    fig.savefig(root/'comparison.png', dpi=170)
    plt.close(fig)


def run(root: Path) -> None:
    """Complete baseline calibration before the positive phase, then summarize without retuning."""
    manifest = full.read(root/'manifest.json')
    protocol = full.read(root/'protocol.json')
    full.verify([*manifest['inputs'], manifest['protocol']])
    study = Path(protocol['study'])
    tasks = [{'name': s['name'], 'site': s, 'phase': 'baseline', 'image': 'baseline'} for s in protocol['sites']]
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (root/'complete.json').exists():
            full.verify(full.read(root/'complete.json')['products'])
            return
        try:
            baseline = [full.read(Path(p)) for p in pool_phase(root, tasks, 'baseline')]
            full.radial.require(all(row['models'][name]['valid'] for r in baseline for row in r['rows'][:-1]
                                for name in METHODS), 'calibration search became invalid')
            thresholds = {r['task']['site']['name']: {name: max(row['models'][name]['search_score']
                for row in r['rows'][:-1]) for name in METHODS} for r in baseline}
            if (root/'thresholds.json').exists():
                full.radial.require(full.read(root/'thresholds.json') == thresholds, 'changed frozen thresholds')
                full.verify(full.read(root/'calibration_complete.json')['products'])
            else:
                write_json(root/'thresholds.json', thresholds)
                write_json(root/'calibration_complete.json', {'all_thresholds_frozen_before_positive_analysis': True,
                    'products': [fingerprint(root/'thresholds.json')]})
            jobs = full.read(study/'jobs.json')
            tasks = [{'name': j['name'], 'site': next(s for s in protocol['sites'] if s['name'] == j['site']),
                      'image': j['name'], 'phase': 'positive', 'job': j} for j in jobs]
            positives = [full.read(Path(p)) for p in pool_phase(root, tasks, 'positive')]
            summarize(root, baseline, positives, thresholds)
            full.verify([*manifest['inputs'], manifest['protocol'], *full.read(root/'calibration_complete.json')['products']])
            write_json(root/'complete.json', {'baseline_maps': 30, 'positive_images': 90, 'new_reductions': 0,
                'inputs_unchanged': True, 'products': [fingerprint(root/name) for name in ('results.json', 'results.md', 'comparison.png', 'thresholds.json')]})
            write_json(root/'state.json', {'status': 'complete', 'new_reductions': 0, 'positive_images': 90})
        except Exception as error:
            write_json(root/'state.json', {'status': 'failed', 'error': str(error)})
            raise


def main() -> None:
    """Prepare or execute the immutable analysis-only comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--study', type=Path, default=Path('working/roc/p4_psd_full_20260918'))
    args = parser.parse_args()
    if args.action == 'prepare':
        args.root.mkdir(parents=True, exist_ok=True)
        prepare(args.root.resolve(), args.study.resolve())
    else:
        run(args.root.resolve())


if __name__ == '__main__':
    main()
