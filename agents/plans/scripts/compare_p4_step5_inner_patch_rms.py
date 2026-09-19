#!/usr/bin/env python3
"""Compare raw and post-mean patch-RMS PSDs on saved inner-radius injections.

Reuse the completed inner-radius baseline and 108 positive reductions.  Freeze
all baseline thresholds before positive analysis, retain the original per-trial
holdouts, and require complete five-pixel annular normalization.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import fcntl
import math
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import subprocess
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import compare_p4_step5_patch_rms as patch_rms
from measure_p4_step5_planet import annular_oracle
import run_p4_step5_inner_rectangular as inner
from run_p4_step5_full_injections import fingerprint, write_json

WIDTHS = inner.WIDTHS
RADII = inner.RADII
LEVELS = inner.LEVELS
RAW = {width: inner.RECTANGULAR[width] for width in WIDTHS}
RMS = {width: RAW[width]+'_patch_rms' for width in WIDTHS}
METHODS = tuple(name for width in WIDTHS for name in (RAW[width], RMS[width])) + ('identity', 'gaussian')
LABELS = {**{RAW[width]: f'Raw ±{width}' for width in WIDTHS},
          **{RMS[width]: f'Patch RMS ±{width}' for width in WIDTHS},
          'identity': 'Identity', 'gaussian': 'Gaussian FWHM 3.6'}
PARENT_INDEX = {name: inner.METHODS.index(name) for name in inner.METHODS}
CALIBRATION_SEARCHES = inner.CALIBRATION_SEARCHES
_CONTEXT = None


def deduplicate(records: list[dict]) -> list[dict]:
    """Return one identical fingerprint for each resolved input path."""
    result = {}
    for record in records:
        path = record['path']
        if path in result:
            inner.full.radial.require(result[path] == record, 'conflicting input fingerprint: '+path)
        result[path] = record
    return list(result.values())


def parent_inputs(study: Path) -> list[dict]:
    """Collect only completed parent products read by this comparison."""
    records = [fingerprint(study/name) for name in ('manifest.json', 'state.json', 'complete.json', 'calibration_complete.json',
        'protocol.json', 'jobs.json', 'results.json', 'baseline.json', 'thresholds.json',
        'effective_calibration_pools.json')]
    records += [fingerprint(path) for path in sorted((study/'payload/response').glob('*.fits'))]
    records.append(fingerprint(study/'payload/baseline.fits'))
    manifest = inner.full.read(study/'manifest.json')
    records += [record for record in manifest['frozen_records']
                if Path(record['path']).parent == study/'software' and Path(record['path']).suffix != '.py']
    analysis_receipts = sorted((study/'analysis').glob('*/complete.json'))
    inner.full.radial.require(len(analysis_receipts) == 272, 'expected 272 completed parent analyses')
    for receipt in analysis_receipts:
        complete = inner.full.read(receipt)
        products = {Path(record['path']).name: record for record in complete['products']}
        inner.full.radial.require('amplitudes.fits' in products and 'amplitudes_snr.fits' in products and
                                  'measurements.json' in products, 'parent analysis receipt is incomplete')
        records += [fingerprint(receipt), products['amplitudes.fits'], products['amplitudes_snr.fits'],
                    products['measurements.json']]
    jobs = inner.full.read(study/'jobs.json')
    inner.full.radial.require(len(jobs) == 108, 'expected 108 completed parent reductions')
    for job in jobs:
        receipt = study/'reductions'/job['name']/'complete.json'
        record = inner.full.read(study/'measurements'/(job['name']+'.json'))
        inner.full.radial.require(record['trial']['name'] == job['name'] and len(record['products']) == 2,
                                  'parent positive measurement does not match its frozen job')
        records += [fingerprint(receipt), record['products'][0], fingerprint(study/'measurements'/(job['name']+'.json'))]
    records += [fingerprint(path) for path in sorted(Path(__file__).parent.glob('*.py'))]
    return deduplicate(records)


def prepare(args: argparse.Namespace) -> None:
    """Freeze the saved-image comparison and all parent products it reads."""
    root, study = args.root.resolve(), args.study.resolve()
    inner.full.radial.require(not root.exists(), 'comparison root already exists')
    inner.full.radial.require(inner.full.read(study/'state.json')['status'] == 'complete' and
                              (study/'complete.json').exists(), 'inner-radius parent study is incomplete')
    inner.full.radial.require(bool(args.cpus) and len(args.cpus) == len(set(args.cpus)) and
                              set(args.cpus).issubset(os.sched_getaffinity(0)), 'requested CPUs are unavailable')
    records = parent_inputs(study)
    inner.full.verify(records)
    parent = inner.full.read(study/'protocol.json')
    jobs = inner.full.read(study/'jobs.json')
    completion = inner.full.read(study/'complete.json')
    inner.full.radial.require(len(parent['calibration_trials']) == 128 and len(parent['sites']) == 36 and
                              len(jobs) == 108 and completion['positive_reductions'] == 108,
                              'changed or incomplete inner-radius parent design')
    trial_names = [trial['name'] for trial in [*parent['calibration_trials'], *parent['sites']]]
    job_names = [job['name'] for job in jobs]
    inner.full.radial.require(len(set(trial_names)) == 164 and len(set(job_names)) == 108,
                              'duplicate parent trial or positive names')
    expected_analyses = {'baseline__'+name for name in trial_names} | set(job_names)
    actual_analyses = {receipt.parent.name for receipt in (study/'analysis').glob('*/complete.json')}
    inner.full.radial.require(actual_analyses == expected_analyses,
                              'parent analysis receipts do not match the frozen trial and job names')
    baseline_tasks = [{'name': 'baseline__'+trial['name'], 'phase': 'baseline',
                       'trial_name': trial['name'], 'source_name': 'baseline'}
                      for trial in [*parent['calibration_trials'], *parent['sites']]]
    positive_tasks = [{'name': job['name'], 'phase': 'positive', 'trial_name': job['site'],
                       'source_name': job['name'], 'job': job} for job in jobs]
    protocol = {'schema': 1,
        'purpose': 'paired raw versus post-mean patch-RMS rectangular PSD at inner radii',
        'parent_study': str(study), 'methods': list(METHODS), 'labels': LABELS,
        'widths': list(WIDTHS), 'radii': list(RADII), 'brightnesses': list(LEVELS),
        'baseline_tasks': baseline_tasks, 'positive_tasks': positive_tasks,
        'new_reductions': 0, 'saved_positive_images': len(positive_tasks),
        'training': 'raw aligned 11x11 patches; five-pixel angular/radial center step; per-trial five-search holdout and buffered known planet excluded; minimum8',
        'raw_estimator': 'ensemble mean subtraction; rectangular 21x21 padded PSD; mixture0.3',
        'patch_normalization': 'subtract the raw ensemble mean patch, divide each residual patch by sqrt(mean(residual**2)), average rectangular powers without second centering, then restore raw mean pixel variance before mixture0.3',
        'candidate_units': 'candidate stamp, fitted raw mean patch, and response remain in raw contrast units',
        'invalid_normalization': 'nonfinite or roundoff-scale patch RMS invalidates the normalized fit; no samples are removed and no fallback is used',
        'search': 'native center plus four axial one-pixel neighbors; all five covariance fits and annular normalizations required; invalid search is a nondetection',
        'snr_summary': 'use the frozen production hciAnalyze SNR maps; at each radius, brightness, and method report the arithmetic mean of the valid five-pixel maximum search SNRs; invalid searches are omitted and counted; also retain mean center-pixel SNR',
        'thresholds': '20 valid geometry-selected baseline searches in the parent radial band; score-free validity replacements only; all thresholds frozen before positives',
        'controls': 'raw PSD, identity, and Gaussian amplitude/SNR maps and decisions must reproduce the completed parent study',
        'workers': len(args.cpus), 'cpu_ids': args.cpus, 'threads_per_worker': 1,
        'dependence': 'reuses the same correlated residual field and inspected injections; paired development comparison, not independent validation'}
    root.mkdir(parents=True)
    write_json(root/'protocol.json', protocol)
    write_json(root/'manifest.json', {'schema': 1, 'inputs': records,
        'protocol': fingerprint(root/'protocol.json'), 'host': os.uname().nodename})
    write_json(root/'state.json', {'status': 'prepared', 'baseline_tasks': len(baseline_tasks),
        'positive_tasks': len(positive_tasks), 'new_reductions': 0})
    print(f'prepared {len(baseline_tasks)} baseline and {len(positive_tasks)} saved-positive analyses', flush=True)


def initialize(root: str, queue) -> None:
    """Pin one worker to one CPU and load its immutable parent context."""
    global _CONTEXT
    os.sched_setaffinity(0, {queue.get()})
    output = Path(root)
    contract = inner.full.read(output/'protocol.json')
    study = Path(contract['parent_study'])
    parent = inner.full.read(study/'protocol.json')
    trials = {trial['name']: trial for trial in [*parent['calibration_trials'], *parent['sites']]}
    _CONTEXT = (output, study, parent, trials, inner.full.load_templates(study/'payload/response'))


def source_path(study: Path, task: dict) -> Path:
    """Return the frozen baseline or positive final image for one task."""
    if task['phase'] == 'baseline':
        return study/'payload/baseline.fits'
    return study/'reductions'/task['source_name']/'finim.fits'


def archive_incomplete(root: Path, directory: Path) -> Path:
    """Preserve one unreceipted task directory before recomputing it."""
    archive = root/'interrupted'/directory.parent.name/directory.name
    archive.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while (archive/f'attempt_{attempt:04d}').exists():
        attempt += 1
    destination = archive/f'attempt_{attempt:04d}'
    shutil.move(directory, destination)
    return destination


def fit_pixel(science: np.ndarray, position: tuple[int, int], template: np.ndarray,
              forbidden: np.ndarray, validate_raw: bool) -> tuple[dict, dict]:
    """Fit all three post-mean patch-RMS models and optional raw controls."""
    x, y = position
    data = science[y-5:y+6, x-5:x+6].ravel()
    rings = inner.full.radial.geometry(science, position, forbidden)
    matrices = {offset: inner.full.radial.extract(science, ring) for offset, ring in rings.items()}
    amplitudes, diagnostics = {}, {}
    for width in WIDTHS:
        samples = np.vstack([matrix for offset, matrix in matrices.items() if abs(offset) <= width])
        model = patch_rms.fit_patch_rms(samples, 'rectangular', .3)
        diagnostic = {'valid': model is not None, 'samples': len(samples)}
        if model is not None:
            result = inner.full.radial.filter_stamp(data, template, model, np.ones(121))
            amplitudes[RMS[width]] = result['amplitude']
            diagnostic.update({key: float(result[key]) for key in ('amplitude', 'sigma', 'score')})
            diagnostic.update(rms_min=float(model['rms'].min()), rms_median=float(np.median(model['rms'])),
                              rms_max=float(model['rms'].max()))
        if validate_raw:
            raw = inner.full.psd.fit_psd(samples, 'rectangular', .3)
            diagnostic['raw_valid'] = raw is not None
            if raw is not None:
                raw_result = inner.full.radial.filter_stamp(data, template, raw, np.ones(121))
                diagnostic['raw_amplitude'] = float(raw_result['amplitude'])
                if model is not None:
                    inner.full.radial.require(np.array_equal(raw['mean'], model['mean']) and
                                              np.isclose(raw['target_variance'], model['target_variance'],
                                                         rtol=1e-12, atol=0),
                                              'patch normalization changed fitted mean or raw variance scale')
        diagnostics[RMS[width]] = diagnostic
    return amplitudes, diagnostics


def select_methods(maps: np.ndarray, trial: dict, settings: dict) -> tuple[list, dict, dict, dict]:
    """Select planes with finite amplitudes and annular SNR at all five search pixels."""
    positions = [(trial['row']+dx, trial['column']+dy) for dx, dy in inner.SEARCH_OFFSETS]
    enabled, expected, profiles, diagnostics = [], {}, {}, {}
    for index, name in enumerate(METHODS):
        expected[name], profiles[name] = annular_oracle(maps[index], settings)
        amplitude_pixels = [bool(np.isfinite(maps[index, y, x])) for x, y in positions]
        annular_pixels = [bool(np.isfinite(expected[name][y, x])) for x, y in positions]
        valid = all(amplitude_pixels) and all(annular_pixels)
        diagnostics[name] = {'amplitude_pixels': amplitude_pixels, 'annular_pixels': annular_pixels,
                             'valid': valid}
        if valid:
            enabled.append(name)
    return enabled, expected, profiles, diagnostics


def score_trial(trial: dict, snr: np.ndarray, amplitudes: np.ndarray) -> dict:
    """Measure the fixed five-pixel search, retaining nulls for unavailable pixels."""
    result = {}
    for index, name in enumerate(METHODS):
        values = [float(snr[index, trial['column']+dy, trial['row']+dx]) for dx, dy in inner.SEARCH_OFFSETS]
        raw = [float(amplitudes[index, trial['column']+dy, trial['row']+dx])
               for dx, dy in inner.SEARCH_OFFSETS]
        valid = all(np.isfinite(values)) and all(np.isfinite(raw))
        result[name] = {'valid': valid, 'search_score': max(values) if valid else None,
            'snr_pixels': [value if np.isfinite(value) else None for value in values],
            'amplitude_pixels': [value if np.isfinite(value) else None for value in raw],
            'center_amplitude': raw[0] if np.isfinite(raw[0]) else None}
    return result


def analyze(task: dict) -> str:
    """Build normalized maps for one frozen image/trial and reproduce all raw controls."""
    root, study, parent, trials, templates = _CONTEXT
    trial = trials[task['trial_name']]
    directory = root/task['phase']/task['name']
    complete = directory/'complete.json'
    if complete.exists():
        inner.full.verify(inner.full.read(complete)['products'])
        return str(directory/'measurements.json')
    if directory.exists():
        destination = archive_incomplete(root, directory)
        print(f'{task["name"]}: archived incomplete task as {destination}', flush=True)
    directory.mkdir(parents=True)
    start = time.monotonic()
    source = source_path(study, task)
    science, header = fits.getdata(source, header=True)
    science = science.squeeze().astype(float)
    parent_directory = study/'analysis'/task['name']
    parent_amplitudes = fits.getdata(parent_directory/'amplitudes.fits').reshape((len(inner.METHODS), *science.shape))
    parent_snr = fits.getdata(parent_directory/'amplitudes_snr.fits').reshape(parent_amplitudes.shape)
    parent_measurement = inner.full.read(parent_directory/'measurements.json')
    yy, xx = np.indices(science.shape)
    radius = np.hypot(xx-127.5, yy-127.5).astype('f4')
    bins, positions, selected = inner.required_bins([trial], radius)
    source_mask, clean_centers = inner.planet_masks(science.shape, parent['known_source_circle'])
    forbidden = source_mask.copy()
    inner.full.mark_footprint(forbidden, trial['row'], trial['column'])
    maps = np.full((len(METHODS), *science.shape), np.nan, dtype='f4')
    for name in (*RAW.values(), 'identity', 'gaussian'):
        maps[METHODS.index(name)] = parent_amplitudes[PARENT_INDEX[name]]
    details = {}
    fitted_pixels = 0
    for (x, y), template in sorted(templates.items()):
        if not selected[y, x] or not clean_centers[y, x]:
            continue
        data = science[y-5:y+6, x-5:x+6]
        if data.size != 121 or not np.isfinite(data).all():
            continue
        amplitudes, diagnostics = fit_pixel(science, (x, y), template, forbidden, (x, y) in positions)
        for name, amplitude in amplitudes.items():
            maps[METHODS.index(name), y, x] = amplitude
        if (x, y) in positions:
            for width in WIDTHS:
                expected_raw = parent_amplitudes[PARENT_INDEX[RAW[width]], y, x]
                diagnostic = diagnostics[RMS[width]]
                inner.full.radial.require(diagnostic['raw_valid'] == bool(np.isfinite(expected_raw)),
                                          'raw PSD validity changed')
                if diagnostic['raw_valid']:
                    inner.full.radial.require(np.isclose(diagnostic['raw_amplitude'], expected_raw,
                                              rtol=2e-6, atol=1e-12), 'raw PSD amplitude changed')
            details[f'{x},{y}'] = diagnostics
        fitted_pixels += 1
    inner.full.radial.require(set(details) == {f'{x},{y}' for x, y in positions},
                              'missing fixed search-pixel covariance diagnostics')
    settings = {'source_x': parent['known_source_circle'][0], 'source_y': parent['known_source_circle'][1],
        'source_radius': parent['known_source_circle'][2], 'lambda_d': 3.6, 'min_radius': 0, 'max_radius': 60}
    enabled, expected, profiles, annular_support = select_methods(maps, trial, settings)
    inner.full.radial.require('identity' in enabled and 'gaussian' in enabled,
                              'reference method lacks complete five-pixel normalization')
    parent_active = set(parent_measurement['active_methods'])
    for name in (*RAW.values(), 'identity', 'gaussian'):
        inner.full.radial.require((name in enabled) == (name in parent_active),
                                  'raw/reference annular validity changed for '+name)
    enabled_indices = [METHODS.index(name) for name in enabled]
    header['HCI FILTER LABELS'] = ','.join(METHODS)
    header['HCI RADIAL BINS'] = ','.join(map(str, bins))
    fits.writeto(directory/'amplitudes.fits', maps, header)
    active_header = header.copy()
    active_header['HCI FILTER LABELS'] = ','.join(enabled)
    fits.writeto(directory/'active_amplitudes.fits', maps[enabled_indices], active_header)
    command = [str(study/'software/hciAnalyze'), '--file='+str(directory/'active_amplitudes.fits'), '--lambdaD=3.6',
        '--planet.sep=11.782', '--planet.PA=262.051', '--planet.R=7.3', '--snr.apertureR=60',
        '--snr.minRad=0', '--snr.maxRad=60', '--filter.psfResponse=', '--filter.lpfGaussFW=0',
        '--filter.hpfGaussFW=0', '--noise.model=identity', '--noise.only=false', '--noise.outputDiagnostics=false']
    write_json(directory/'command.json', command)
    environment = inner.full.binary_environment(study, reduction=False)
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (directory/'analysis.log').open('w') as log:
        subprocess.run(command, cwd=directory, env=environment, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    active_snr, snr_header = fits.getdata(directory/'active_amplitudes_snr.fits', header=True)
    active_snr = active_snr.reshape((len(enabled), *science.shape))
    inner.full.radial.require(active_snr.shape == maps[enabled_indices].shape and
                              snr_header['SNRMEAN'] == snr_header['SNRSMALL'] == 1,
                              'changed production annular SNR contract')
    snr = np.full(maps.shape, np.nan, dtype='f4')
    snr[enabled_indices] = active_snr
    fits.writeto(directory/'amplitudes_snr.fits', snr, snr_header)
    oracle_errors = {}
    for index, name in enumerate(METHODS):
        if name not in enabled:
            oracle_errors[name] = None
            continue
        valid = np.isfinite(maps[index]) & np.isfinite(expected[name])
        difference = np.abs(snr[index][valid]-expected[name][valid])
        oracle_errors[name] = float(np.max(difference))
        inner.full.radial.require(np.allclose(snr[index][valid], expected[name][valid],
                                  rtol=2e-6, atol=2e-6), 'annular oracle mismatch for '+name)
    control_errors = {}
    for name in (*RAW.values(), 'identity', 'gaussian'):
        current = snr[METHODS.index(name)]
        prior = parent_snr[PARENT_INDEX[name]]
        control_errors[name] = float(np.nanmax(np.abs(current-prior))) if np.isfinite(prior).any() else None
        inner.full.radial.require(np.allclose(current, prior, rtol=2e-6, atol=2e-6, equal_nan=True),
                                  'parent SNR map changed for '+name)
    measurement = {'task': task, 'trial': trial, 'source': fingerprint(source),
        'models': score_trial(trial, snr, maps), 'active_methods': enabled,
        'annular_support': annular_support, 'required_bins': bins, 'fit_details': details,
        'fitted_pixels': fitted_pixels, 'elapsed_seconds': time.monotonic()-start,
        'annular_oracle_max_errors': oracle_errors, 'parent_control_max_errors': control_errors}
    write_json(directory/'profiles.json', {name: [row for row in profile if int(row['radius']-.5) in bins]
                                           for name, profile in profiles.items()})
    write_json(directory/'measurements.json', measurement)
    products = [fingerprint(directory/name) for name in ('amplitudes.fits', 'active_amplitudes.fits',
        'active_amplitudes_snr.fits', 'amplitudes_snr.fits', 'measurements.json', 'profiles.json', 'command.json')]
    write_json(complete, {'products': products})
    return str(directory/'measurements.json')


def run_phase(root: Path, tasks: list, phase: str) -> list[dict]:
    """Run one resumable analysis phase on pinned physical CPU workers."""
    protocol = inner.full.read(root/'protocol.json')
    context = mp.get_context('spawn')
    queue = context.Queue()
    for cpu in protocol['cpu_ids']:
        queue.put(cpu)
    paths = []
    with ProcessPoolExecutor(max_workers=protocol['workers'], mp_context=context,
                             initializer=initialize, initargs=(str(root), queue)) as executor:
        futures = {executor.submit(analyze, task): task['name'] for task in tasks}
        for future in as_completed(futures):
            path = future.result()
            paths.append(path)
            write_json(root/'state.json', {'status': phase, 'pid': os.getpid(),
                'finished': len(paths), 'total': len(tasks), 'last': futures[future]})
            print(f'{phase}: {len(paths)}/{len(tasks)} {futures[future]}', flush=True)
    return [inner.full.read(Path(path)) for path in sorted(paths)]


def base_method(method: str) -> str:
    """Map a normalized method to the corresponding parent calibration policy."""
    for width in WIDTHS:
        if method in (RAW[width], RMS[width]):
            return RAW[width]
    return method


def effective_pools(parent: dict, baseline: list[dict]) -> dict:
    """Choose fixed-size valid null pools by geometry without reading score values."""
    rows = {record['trial']['name']: record for record in baseline}
    trials = {trial['name']: trial for trial in parent['calibration_trials']}
    result = {}
    for nominal in RADII:
        result[str(nominal)] = {}
        for method in METHODS:
            pool = parent['calibration_pools'][str(nominal)][base_method(method)]
            if pool is None:
                result[str(nominal)][method] = None
                continue
            original = list(pool['trials'])
            chosen = [name for name in original if rows[name]['models'][method]['valid']]
            invalid = [name for name in original if name not in chosen]
            lower, upper = pool['radius_range']
            available = [trial for trial in trials.values() if trial['name'] not in chosen and
                         lower <= trial['nominal_radius'] <= upper and rows[trial['name']]['models'][method]['valid']]
            replacements = []
            while len(chosen) < CALIBRATION_SEARCHES:
                inner.full.radial.require(bool(available),
                    f'insufficient valid calibration searches for radius {nominal} {method}')
                prior = [trials[name] for name in chosen]
                choice = max(available, key=lambda trial: (
                    min(math.dist((trial['row'], trial['column']), (row['row'], row['column'])) for row in prior),
                    -trial['azimuth_degrees'], -trial['row'], -trial['column'])) if prior else min(
                        available, key=lambda trial: (trial['azimuth_degrees'], trial['row'], trial['column']))
                chosen.append(choice['name'])
                replacements.append(choice['name'])
                available.remove(choice)
            result[str(nominal)][method] = {**pool, 'trials': chosen, 'original_trials': original,
                'invalid_original_trials': invalid, 'replacement_trials': replacements}
    for nominal in RADII:
        for width in WIDTHS:
            raw, normalized = result[str(nominal)][RAW[width]], result[str(nominal)][RMS[width]]
            if raw is None or normalized is None:
                inner.full.radial.require(raw is normalized, 'paired method availability changed')
            else:
                inner.full.radial.require(raw['trials'] == normalized['trials'],
                                          f'patch normalization changed calibration support at radius {nominal} width {width}')
    return result


def calibrate(root: Path, baseline: list[dict]) -> tuple[dict, dict]:
    """Freeze paired effective pools and thresholds before reading positive maps."""
    record = root/'calibration_complete.json'
    if record.exists():
        inner.full.verify(inner.full.read(record)['products'])
        return inner.full.read(root/'thresholds.json'), inner.full.read(root/'effective_calibration_pools.json')['pools']
    study = Path(inner.full.read(root/'protocol.json')['parent_study'])
    parent = inner.full.read(study/'protocol.json')
    pools = effective_pools(parent, baseline)
    rows = {record['trial']['name']: record for record in baseline}
    thresholds = {}
    for nominal in RADII:
        thresholds[str(nominal)] = {}
        for method in METHODS:
            pool = pools[str(nominal)][method]
            if pool is None:
                thresholds[str(nominal)][method] = None
                continue
            values = [rows[name]['models'][method] for name in pool['trials']]
            inner.full.radial.require(len(values) == CALIBRATION_SEARCHES and all(value['valid'] for value in values),
                                      f'invalid effective pool for radius {nominal} {method}')
            thresholds[str(nominal)][method] = max(value['search_score'] for value in values)
    parent_thresholds = inner.full.read(study/'thresholds.json')
    parent_pools = inner.full.read(study/'effective_calibration_pools.json')['pools']
    for nominal in RADII:
        for method in (*RAW.values(), 'identity', 'gaussian'):
            inner.full.radial.require(pools[str(nominal)][method] == parent_pools[str(nominal)][method] and
                                      thresholds[str(nominal)][method] == parent_thresholds[str(nominal)][method],
                                      'raw calibration control changed for '+method)
    write_json(root/'baseline.json', {'rows': baseline, 'parent_baseline_reused': True})
    write_json(root/'effective_calibration_pools.json', {'schema': 1, 'score_values_used_for_selection': False,
        'paired_raw_normalized_locations': True, 'pools': pools})
    write_json(root/'thresholds.json', thresholds)
    products = [fingerprint(root/name) for name in
                ('baseline.json', 'effective_calibration_pools.json', 'thresholds.json')]
    write_json(record, {'thresholds_frozen_before_positive_analysis': True, 'products': products})
    return thresholds, pools


def summarize(root: Path, baseline: list[dict], positives: list[dict], thresholds: dict, pools: dict) -> None:
    """Record per-radius recovery and every paired raw/normalized decision change."""
    study = Path(inner.full.read(root/'protocol.json')['parent_study'])
    parent_results = inner.full.read(study/'results.json')
    baseline_rows = {record['trial']['name']: record for record in baseline}
    jobs = {record['task']['job']['name']: record for record in positives}
    parent_jobs = {record['trial']['name']: record for record in parent_results['measurements']}
    nulls, measurements = [], []
    for site in inner.full.read(study/'protocol.json')['sites']:
        record = baseline_rows[site['name']]
        for method, value in record['models'].items():
            threshold = thresholds[str(site['nominal_radius'])][method]
            value['detected'] = value['valid'] and threshold is not None and value['search_score'] > threshold
        nulls.append(record)
    for name, record in jobs.items():
        job = record['task']['job']
        for method, value in record['models'].items():
            threshold = thresholds[str(job['nominal_radius'])][method]
            value['detected'] = value['valid'] and threshold is not None and value['search_score'] > threshold
            if value['center_amplitude'] is not None and method != 'gaussian':
                value['raw_contrast_error'] = value['center_amplitude']/job['contrast']-1
        measurements.append(record)
        for method in (*RAW.values(), 'identity', 'gaussian'):
            inner.full.radial.require(record['models'][method]['detected'] ==
                                      parent_jobs[name]['models'][method]['detected'],
                                      'parent positive decision changed for '+method)
    parent_nulls = {record['site']['name']: record for record in parent_results['nulls']}
    for record in nulls:
        name = record['trial']['name']
        for method in (*RAW.values(), 'identity', 'gaussian'):
            inner.full.radial.require(record['models'][method]['detected'] ==
                                      parent_nulls[name]['models'][method]['detected'],
                                      'parent held-out null decision changed for '+method)
    groups = []
    for nominal in RADII:
        for method in METHODS:
            site_rows = [record for record in nulls if record['trial']['nominal_radius'] == nominal]
            one = {'radius': nominal, 'method': method,
                'null_exceedances': sum(record['models'][method]['detected'] for record in site_rows),
                'invalid_nulls': sum(not record['models'][method]['valid'] for record in site_rows), 'levels': []}
            for level in LEVELS:
                selected = [record for record in measurements if record['task']['job']['nominal_radius'] == nominal and
                            record['task']['job']['brightness_multiplier'] == level]
                values = [record['models'][method] for record in selected]
                errors = [value['raw_contrast_error'] for value in values if 'raw_contrast_error' in value]
                one['levels'].append({'brightness_multiplier': level, 'trials': len(values),
                    'detections': sum(value['detected'] for value in values),
                    'invalid_searches': sum(not value['valid'] for value in values),
                    'median_raw_contrast_error': float(np.median(errors)) if errors else None})
            groups.append(one)
            if method in (*RAW.values(), 'identity', 'gaussian'):
                parent = next(group for group in parent_results['groups'] if
                              group['radius'] == nominal and group['method'] == method)
                inner.full.radial.require(one == parent, 'parent summary changed for '+method)
            for level_row in one['levels']:
                selected = [record for record in measurements if
                            record['task']['job']['nominal_radius'] == nominal and
                            record['task']['job']['brightness_multiplier'] ==
                            level_row['brightness_multiplier']]
                values = [record['models'][method] for record in selected]
                valid = [value for value in values if value['valid']]
                level_row.update(valid_searches=len(valid),
                    mean_search_snr=float(np.mean([value['search_score'] for value in valid])) if valid else None,
                    mean_center_snr=float(np.mean([value['snr_pixels'][0] for value in valid])) if valid else None)
    paired = []
    for nominal in RADII:
        for width in WIDTHS:
            raw, normalized = RAW[width], RMS[width]
            for level in LEVELS:
                selected = [record for record in measurements if record['task']['job']['nominal_radius'] == nominal and
                            record['task']['job']['brightness_multiplier'] == level]
                paired.append({'radius': nominal, 'width': width, 'brightness_multiplier': level,
                    'normalized_only': [record['task']['job']['name'] for record in selected
                                        if record['models'][normalized]['detected'] and not record['models'][raw]['detected']],
                    'raw_only': [record['task']['job']['name'] for record in selected
                                 if record['models'][raw]['detected'] and not record['models'][normalized]['detected']]})
    aggregate = []
    for method in METHODS:
        method_groups = [group for group in groups if group['method'] == method]
        invalid_by_level = [sum(group['levels'][index]['invalid_searches'] for group in method_groups)
                            for index in range(len(LEVELS))]
        inner.full.radial.require(len(set(invalid_by_level)) == 1,
                                  'method support changed with injected brightness for '+method)
        aggregate_levels = []
        for level in LEVELS:
            selected = [record['models'][method] for record in measurements if
                        record['task']['job']['brightness_multiplier'] == level]
            valid = [value for value in selected if value['valid']]
            aggregate_levels.append({'brightness_multiplier': level,
                'detections': sum(value['detected'] for value in selected), 'trials': len(selected),
                'valid_searches': len(valid),
                'mean_search_snr': float(np.mean([value['search_score'] for value in valid])) if valid else None,
                'mean_center_snr': float(np.mean([value['snr_pixels'][0] for value in valid])) if valid else None})
        aggregate.append({'method': method, 'valid_sites': len(RADII)*6-invalid_by_level[0],
            'null_exceedances': sum(group['null_exceedances'] for group in method_groups),
            'invalid_nulls': sum(group['invalid_nulls'] for group in method_groups),
            'levels': aggregate_levels})
    verification = {'new_reductions': 0,
        'maximum_parent_control_snr_difference': max(value for record in [*baseline, *measurements]
            for value in record['parent_control_max_errors'].values() if value is not None),
        'maximum_annular_oracle_difference': max(value for record in [*baseline, *measurements]
            for value in record['annular_oracle_max_errors'].values() if value is not None),
        'raw_reference_decisions_reproduced': True, 'paired_calibration_locations': True}
    write_json(root/'results.json', {'groups': groups, 'aggregate': aggregate,
        'paired_decisions': paired, 'nulls': nulls,
        'measurements': measurements, 'thresholds': thresholds, 'calibration_pools': pools,
        'verification': verification, 'caveats': [inner.full.read(study/'protocol.json')['dependence'],
        'This reuses inspected images and does not provide independent validation.']})
    lines = ['# Inner-radius post-mean patch-RMS comparison', '',
        '| Radius | Method | 0.5× | 0.75× | 1× | Nulls | Invalid nulls |',
        '| ---: | --- | ---: | ---: | ---: | ---: | ---: |']
    for group in groups:
        lines.append(f'| {group["radius"]} | {LABELS[group["method"]]} | '+
            ' | '.join(f'{row["detections"]}/{row["trials"]}' for row in group['levels'])+
            f' | {group["null_exceedances"]}/6 | {group["invalid_nulls"]}/6 |')
    lines += ['', 'All thresholds were frozen before positive analysis. Invalid searches are nondetections. '
              'Raw controls reproduce the parent study; no P4 reductions were run.', '',
        '## Mean five-pixel search SNR by radius', '',
        'Each per-injection value comes from the frozen production hciAnalyze SNR map. The table gives the arithmetic '
        'mean over valid injections at that radius and shows the valid count. The reported SNR is the maximum over '
        'the same fixed five-pixel search used for detection.', '',
        '| Radius | Method | Valid sites | 0.5× mean SNR | 0.75× mean SNR | 1× mean SNR |',
        '| ---: | --- | ---: | ---: | ---: | ---: |']
    for group in groups:
        valid_counts = [row['valid_searches'] for row in group['levels']]
        inner.full.radial.require(len(set(valid_counts)) == 1,
                                  'method support changed with brightness within one radius')
        formatted = [('—' if row['mean_search_snr'] is None else f'{row["mean_search_snr"]:.4f}')
                     for row in group['levels']]
        lines.append(f'| {group["radius"]} | {LABELS[group["method"]]} | {valid_counts[0]}/6 | '+
                     ' | '.join(formatted)+' |')
    lines += ['',
        '## Aggregate across radii', '',
        '| Method | Valid sites | 0.5× | 0.75× | 1× | Nulls | Invalid nulls |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for row in aggregate:
        lines.append(f'| {LABELS[row["method"]]} | {row["valid_sites"]}/36 | '+
            ' | '.join(str(level['detections']) for level in row['levels'])+
            f' | {row["null_exceedances"]}/36 | {row["invalid_nulls"]}/36 |')
    lines += ['', '## Paired decision changes', '']
    for width in WIDTHS:
        for level in LEVELS:
            selected = [row for row in paired if row['width'] == width and row['brightness_multiplier'] == level]
            normalized_only = sum(len(row['normalized_only']) for row in selected)
            raw_only = sum(len(row['raw_only']) for row in selected)
            lines.append(f'- ±{width}, {level:g}×: patch RMS only {normalized_only}; raw only {raw_only}.')
    lines.append('')
    (root/'results.md').write_text('\n'.join(lines))
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True, layout='constrained')
    for axis, nominal in zip(axes.flat, RADII):
        for method in METHODS:
            group = next(row for row in groups if row['radius'] == nominal and row['method'] == method)
            axis.plot(LEVELS, [row['detections']/6 for row in group['levels']], 'o-', label=LABELS[method])
        axis.set(title=f'r = {nominal} px', xlabel='Brightness / identity threshold scale',
                 ylabel='Recovery fraction', ylim=(-.03, 1.03))
    axes[0, 0].legend(fontsize=6)
    fig.suptitle('Inner-radius P4: post-mean patch-RMS PSD weighting')
    fig.savefig(root/'comparison.png', dpi=170)
    plt.close(fig)


def run(root: Path) -> None:
    """Run baseline calibration, then analyze all saved positives without new reductions."""
    manifest = inner.full.read(root/'manifest.json')
    inner.full.verify([*manifest['inputs'], manifest['protocol']])
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (root/'complete.json').exists():
            inner.full.verify(inner.full.read(root/'complete.json')['products'])
            print('Inner-radius patch-RMS comparison already complete.', flush=True)
            return
        protocol = inner.full.read(root/'protocol.json')
        try:
            write_json(root/'state.json', {'status': 'baseline', 'pid': os.getpid(), 'new_reductions': 0})
            baseline = run_phase(root, protocol['baseline_tasks'], 'baseline')
            thresholds, pools = calibrate(root, baseline)
            calibration = inner.full.read(root/'calibration_complete.json')
            write_json(root/'state.json', {'status': 'positive', 'pid': os.getpid(), 'new_reductions': 0})
            positives = run_phase(root, protocol['positive_tasks'], 'positive')
            summarize(root, baseline, positives, thresholds, pools)
            inner.full.verify([*manifest['inputs'], manifest['protocol'], *calibration['products']])
            products = [fingerprint(root/name) for name in ('results.json', 'results.md', 'comparison.png',
                'thresholds.json', 'effective_calibration_pools.json', 'baseline.json')]
            write_json(root/'complete.json', {'baseline_analyses': len(baseline),
                'positive_analyses': len(positives), 'new_reductions': 0,
                'frozen_inputs_unchanged': True, 'products': products})
            write_json(root/'state.json', {'status': 'complete', 'new_reductions': 0,
                'results': str(root/'results.json')})
        except Exception as error:
            write_json(root/'state.json', {'status': 'failed', 'pid': os.getpid(),
                'new_reductions': 0, 'error': str(error)})
            raise


def main() -> None:
    """Expose immutable preparation and resumable paired comparison actions."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--study', type=Path, default=Path('working/roc/p4_inner_rectangular_20260919'))
    parser.add_argument('--cpus', type=int, nargs='+', default=list(range(12)))
    args = parser.parse_args()
    if args.action == 'prepare':
        prepare(args)
    else:
        run(args.root.resolve())


if __name__ == '__main__':
    main()
