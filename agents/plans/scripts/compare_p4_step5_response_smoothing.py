#!/usr/bin/env python3
"""Compare smoothed measured-response templates with identity and PSD weights.

Reuse the completed inner SNR-3/5/7 images and lambda/D-masked annular
analysis. Smooth only each 11x11 measured response, apply the same template
grid to identity and raw rectangular-PSD weights, and freeze all null
thresholds before reading positive measurements.
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
from scipy.ndimage import gaussian_filter

import compare_p4_step5_inner_patch_rms as previous
from measure_p4_step5_planet import annular_oracle
import run_p4_step5_inner_rectangular as inner
from run_p4_step5_full_injections import fingerprint, write_json

LAMBDA_D = 3.6
FWHMS = (0.0, 0.9, 1.8, 2.7, 3.6)
RADII = inner.RADII
RECTANGULAR_WIDTH = 20
PARENT_RECTANGULAR = inner.RECTANGULAR[RECTANGULAR_WIDTH]
CALIBRATION_SEARCHES = inner.CALIBRATION_SEARCHES


def width_tag(fwhm: float) -> str:
    """Return a stable method-name suffix for one smoothing width."""
    return ('%.1f' % fwhm).replace('.', 'p')


IDENTITY = {fwhm: 'identity_t'+width_tag(fwhm) for fwhm in FWHMS}
RECTANGULAR = {fwhm: 'rect20_t'+width_tag(fwhm) for fwhm in FWHMS}
METHODS = tuple(IDENTITY.values()) + tuple(RECTANGULAR.values()) + ('gaussian',)
LABELS = {**{name: f'Identity, response LPF {fwhm:.1f} px' for fwhm, name in IDENTITY.items()},
          **{name: f'Rectangular ±20, response LPF {fwhm:.1f} px'
             for fwhm, name in RECTANGULAR.items()},
          'gaussian': 'Gaussian FWHM 3.6'}
FAMILY = {**{name: 'identity' for name in IDENTITY.values()},
          **{name: 'rectangular20' for name in RECTANGULAR.values()}, 'gaussian': 'gaussian'}
BASE_IDENTITY = IDENTITY[0.0]
BASE_RECTANGULAR = RECTANGULAR[0.0]
CONTROL_PARENT = {BASE_IDENTITY: 'identity', BASE_RECTANGULAR: PARENT_RECTANGULAR,
                  'gaussian': 'gaussian'}
PARENT_INDEX = {name: previous.METHODS.index(name) for name in previous.METHODS}
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


def comparison_inputs(comparison: Path, protocol: dict) -> list[dict]:
    """Collect completed comparison products and inherited frozen inputs read here."""
    manifest = inner.full.read(comparison/'manifest.json')
    records = [*manifest['inputs'], manifest['protocol'], *manifest.get('repair_records', [])]
    records += [fingerprint(comparison/name) for name in ('manifest.json', 'state.json', 'complete.json',
        'calibration_complete.json', 'protocol.json', 'results.json', 'baseline.json', 'thresholds.json',
        'effective_calibration_pools.json')]
    tasks = [*protocol['baseline_tasks'], *protocol['positive_tasks']]
    inner.full.radial.require(len(tasks) == 272, 'expected 272 completed parent analyses')
    for task in tasks:
        directory = comparison/task['phase']/task['name']
        complete = inner.full.read(directory/'complete.json')
        products = {Path(record['path']).name: record for record in complete['products']}
        inner.full.radial.require('amplitudes.fits' in products and 'amplitudes_snr.fits' in products,
                                  'parent comparison task is incomplete')
        records += [fingerprint(directory/'complete.json'), products['amplitudes.fits'],
                    products['amplitudes_snr.fits']]
    records += [fingerprint(Path(__file__)), fingerprint(Path(previous.__file__)),
                fingerprint(Path(inner.__file__))]
    return deduplicate(records)


def prepare(args: argparse.Namespace) -> None:
    """Freeze the paired saved-image template-smoothing comparison."""
    root, comparison = args.root.resolve(), args.comparison.resolve()
    inner.full.radial.require(not root.exists(), 'comparison root already exists')
    inner.full.radial.require((comparison/'complete.json').exists() and
                              inner.full.read(comparison/'state.json')['status'] == 'complete',
                              'lambda/D-masked parent comparison is incomplete')
    inner.full.radial.require(bool(args.cpus) and len(args.cpus) == len(set(args.cpus)) and
                              set(args.cpus).issubset(os.sched_getaffinity(0)),
                              'requested CPUs are unavailable')
    prior = inner.full.read(comparison/'protocol.json')
    study = Path(prior['parent_study'])
    parent = inner.full.read(study/'protocol.json')
    inner.full.radial.require(prior.get('annular_trial_exclusion_radius') == 3.1 and
                              prior.get('parent_snr_replay') is False,
                              'parent must be the completed one-lambda/D masked comparison')
    inner.full.radial.require(parent.get('target_source_snrs') == [3.0, 5.0, 7.0] and
                              len(parent['sites']) == 36 and len(parent['calibration_trials']) == 128,
                              'changed SNR-3/5/7 parent design')
    records = comparison_inputs(comparison, prior)
    inner.full.verify(records)
    protocol = {'schema': 1,
        'purpose': 'measured-response low-pass grid under identity and raw rectangular ±20 PSD weights',
        'parent_comparison': str(comparison), 'parent_study': str(study),
        'methods': list(METHODS), 'labels': LABELS, 'families': FAMILY,
        'response_smoothing_fwhm_pixels': list(FWHMS),
        'response_smoothing_fwhm_lambda_d': [value/LAMBDA_D for value in FWHMS],
        'response_smoothing': 'scipy.ndimage.gaussian_filter on each native 11x11 measured response; constant-zero exterior; truncate4; no amplitude renormalization',
        'identity': 'template dot candidate divided by template energy; no fitted-mean subtraction',
        'rectangular_psd': 'raw patches, ±20-pixel training centers, ensemble-mean subtraction, rectangular 21x21 padded PSD, mixture0.3',
        'template_only_change': 'candidate images, covariance samples, fitted means, exclusions, annular SNR, search geometry, and source injections are unchanged',
        'radii': list(RADII), 'brightnesses': parent['brightness_multipliers'],
        'target_source_snrs': parent['target_source_snrs'], 'level_labels': parent['level_labels'],
        'baseline_tasks': prior['baseline_tasks'], 'positive_tasks': prior['positive_tasks'],
        'calibration_pools': 'reuse the frozen lambda/D-masked identity or raw ±20 locations for every smoothing width in that family; each method receives its own maximum-null threshold',
        'annular_trial_exclusion_radius': prior['annular_trial_exclusion_radius'],
        'annular_trial_effective_radius': prior['annular_trial_effective_radius'],
        'annular_exclusions': prior['annular_exclusions'],
        'search': 'native center plus four axial one-pixel neighbors; all five amplitude and annular-SNR pixels required',
        'new_reductions': 0, 'saved_positive_images': len(prior['positive_tasks']),
        'workers': len(args.cpus), 'cpu_ids': args.cpus, 'threads_per_worker': 1,
        'dependence': 'reuses the same correlated residual field and inspected injections; paired development comparison, not independent validation'}
    root.mkdir(parents=True)
    write_json(root/'protocol.json', protocol)
    write_json(root/'manifest.json', {'schema': 1, 'inputs': records,
        'protocol': fingerprint(root/'protocol.json'), 'host': os.uname().nodename})
    write_json(root/'state.json', {'status': 'prepared',
        'baseline_tasks': len(prior['baseline_tasks']), 'positive_tasks': len(prior['positive_tasks']),
        'new_reductions': 0})
    print(f'prepared {len(prior["baseline_tasks"])} baseline and '
          f'{len(prior["positive_tasks"])} saved-positive analyses', flush=True)


def initialize(root: str, queue) -> None:
    """Pin one worker to one CPU and load its immutable response field."""
    global _CONTEXT
    os.sched_setaffinity(0, {queue.get()})
    output = Path(root)
    protocol = inner.full.read(output/'protocol.json')
    comparison = Path(protocol['parent_comparison'])
    study = Path(protocol['parent_study'])
    parent = inner.full.read(study/'protocol.json')
    trials = {trial['name']: trial for trial in [*parent['calibration_trials'], *parent['sites']]}
    _CONTEXT = (output, comparison, study, parent, trials,
                inner.full.load_templates(study/'payload/response'))


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


def smooth_template(template: np.ndarray, fwhm: float) -> np.ndarray:
    """Low-pass one native response stamp without changing the zero-width control."""
    stamp = np.asarray(template, dtype=float).reshape((11, 11))
    if fwhm == 0:
        return stamp.ravel().copy()
    sigma = fwhm/math.sqrt(8*math.log(2))
    return gaussian_filter(stamp, sigma=sigma, mode='constant', cval=0.0,
                           truncate=4.0).ravel()


def template_diagnostic(filtered: np.ndarray, original: np.ndarray) -> dict:
    """Measure energy and alignment of one smoothed response with its original."""
    original_energy = float(original @ original)
    filtered_energy = float(filtered @ filtered)
    cross = float(filtered @ original)
    inner.full.radial.require(original_energy > 0 and filtered_energy > 0,
                              'response template has zero energy')
    return {'template_energy_fraction': filtered_energy/original_energy,
            'template_cosine': cross/math.sqrt(filtered_energy*original_energy)}


def fit_pixel(science: np.ndarray, position: tuple[int, int], template: np.ndarray,
              forbidden: np.ndarray) -> tuple[dict, dict]:
    """Evaluate both weight families across the fixed response-smoothing grid."""
    x, y = position
    data = science[y-5:y+6, x-5:x+6].ravel()
    rings = inner.full.radial.geometry(science, position, forbidden)
    matrices = {offset: inner.full.radial.extract(science, ring) for offset, ring in rings.items()
                if abs(offset) <= RECTANGULAR_WIDTH}
    samples = np.vstack(list(matrices.values()))
    model = inner.full.psd.fit_psd(samples, 'rectangular', .3)
    original_psd = (inner.full.radial.filter_stamp(data, template, model, np.ones(121))
                    if model is not None else None)
    amplitudes, diagnostics = {}, {}
    for fwhm in FWHMS:
        filtered = smooth_template(template, fwhm)
        common = template_diagnostic(filtered, template)
        identity_name = IDENTITY[fwhm]
        energy = float(filtered @ filtered)
        weight = filtered/energy
        amplitudes[identity_name] = float(weight @ data)
        diagnostics[identity_name] = {**common, 'valid': True,
            'amplitude': amplitudes[identity_name],
            'expected_original_response': float(weight @ template),
            'matched_snr_efficiency_if_original_exact': common['template_cosine']}
        rectangular_name = RECTANGULAR[fwhm]
        diagnostics[rectangular_name] = {**common, 'valid': model is not None,
                                         'samples': len(samples)}
        if model is not None:
            result = (original_psd if fwhm == 0 else
                      inner.full.radial.filter_stamp(data, filtered, model, np.ones(121)))
            amplitudes[rectangular_name] = float(result['amplitude'])
            weight = result['physical_weight']
            expected = float(weight @ template)
            efficiency = expected*float(original_psd['sigma'])/float(result['sigma'])
            inner.full.radial.require(abs(efficiency) <= 1+1e-12,
                                      'PSD template-match efficiency is outside its algebraic bound')
            diagnostics[rectangular_name].update(
                amplitude=float(result['amplitude']), sigma=float(result['sigma']),
                score=float(result['score']), expected_original_response=expected,
                matched_snr_efficiency_if_original_exact=efficiency)
    return amplitudes, diagnostics


def select_methods(maps: np.ndarray, trial: dict, settings: dict) -> tuple[list, dict, dict, dict]:
    """Select planes with finite amplitudes and annular SNR at all search pixels."""
    positions = [(trial['row']+dx, trial['column']+dy) for dx, dy in inner.SEARCH_OFFSETS]
    enabled, expected, profiles, diagnostics = [], {}, {}, {}
    for index, name in enumerate(METHODS):
        expected[name], profiles[name] = annular_oracle(maps[index], settings)
        amplitude_pixels = [bool(np.isfinite(maps[index, y, x])) for x, y in positions]
        annular_pixels = [bool(np.isfinite(expected[name][y, x])) for x, y in positions]
        valid = all(amplitude_pixels) and all(annular_pixels)
        diagnostics[name] = {'amplitude_pixels': amplitude_pixels,
                             'annular_pixels': annular_pixels, 'valid': valid}
        if valid:
            enabled.append(name)
    return enabled, expected, profiles, diagnostics


def score_trial(trial: dict, snr: np.ndarray, amplitudes: np.ndarray) -> dict:
    """Measure the fixed five-pixel search, retaining nulls for unavailable pixels."""
    result = {}
    for index, name in enumerate(METHODS):
        values = [float(snr[index, trial['column']+dy, trial['row']+dx])
                  for dx, dy in inner.SEARCH_OFFSETS]
        raw = [float(amplitudes[index, trial['column']+dy, trial['row']+dx])
               for dx, dy in inner.SEARCH_OFFSETS]
        valid = all(np.isfinite(values)) and all(np.isfinite(raw))
        result[name] = {'valid': valid, 'search_score': max(values) if valid else None,
            'snr_pixels': [value if np.isfinite(value) else None for value in values],
            'amplitude_pixels': [value if np.isfinite(value) else None for value in raw],
            'center_amplitude': raw[0] if np.isfinite(raw[0]) else None}
    return result


def maximum_finite_difference(current: np.ndarray, prior: np.ndarray) -> float | None:
    """Return the largest difference on common finite support, or null without overlap."""
    overlap = np.isfinite(current) & np.isfinite(prior)
    return float(np.max(np.abs(current[overlap]-prior[overlap]))) if overlap.any() else None


def analyze(task: dict) -> str:
    """Build every smoothed-template map and verify unsmoothed parent controls."""
    root, comparison, study, parent, trials, templates = _CONTEXT
    protocol = inner.full.read(root/'protocol.json')
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
    parent_directory = comparison/task['phase']/task['name']
    parent_amplitudes = fits.getdata(parent_directory/'amplitudes.fits').reshape(
        (len(previous.METHODS), *science.shape))
    parent_snr = fits.getdata(parent_directory/'amplitudes_snr.fits').reshape(parent_amplitudes.shape)
    yy, xx = np.indices(science.shape)
    radius = np.hypot(xx-127.5, yy-127.5).astype('f4')
    bins, positions, selected = inner.required_bins([trial], radius)
    source_mask, clean_centers = inner.planet_masks(science.shape, parent['known_source_circle'])
    forbidden = source_mask.copy()
    inner.full.mark_footprint(forbidden, trial['row'], trial['column'])
    maps = np.full((len(METHODS), *science.shape), np.nan, dtype='f4')
    maps[METHODS.index('gaussian')] = parent_amplitudes[PARENT_INDEX['gaussian']]
    details = {}
    fitted_pixels = 0
    for (x, y), template in sorted(templates.items()):
        if not selected[y, x] or not clean_centers[y, x]:
            continue
        data = science[y-5:y+6, x-5:x+6]
        if data.size != 121 or not np.isfinite(data).all():
            continue
        amplitudes, diagnostics = fit_pixel(science, (x, y), template, forbidden)
        for name, amplitude in amplitudes.items():
            maps[METHODS.index(name), y, x] = amplitude
        if (x, y) in positions:
            details[f'{x},{y}'] = diagnostics
        fitted_pixels += 1
    inner.full.radial.require(set(details) == {f'{x},{y}' for x, y in positions},
                              'missing fixed search-pixel diagnostics')
    amplitude_control_errors = {}
    for method, parent_method in CONTROL_PARENT.items():
        current = maps[METHODS.index(method)]
        prior = parent_amplitudes[PARENT_INDEX[parent_method]]
        amplitude_control_errors[method] = maximum_finite_difference(current, prior)
        inner.full.radial.require(np.allclose(current, prior, rtol=2e-6, atol=1e-12, equal_nan=True),
                                  'unsmoothed amplitude control changed for '+method)
    settings = {'source_x': parent['known_source_circle'][0],
        'source_y': parent['known_source_circle'][1], 'source_radius': parent['known_source_circle'][2],
        'lambda_d': LAMBDA_D, 'min_radius': 0, 'max_radius': 60,
        'source_exclusions': [parent['known_source_circle'],
            [trial['row'], trial['column'], protocol['annular_trial_exclusion_radius']]]}
    oracle_enabled, expected, profiles, annular_support = select_methods(maps, trial, settings)
    search_positions = [(trial['row']+dx, trial['column']+dy) for dx, dy in inner.SEARCH_OFFSETS]
    parent_active = {name: all(np.isfinite(parent_snr[PARENT_INDEX[parent_name], y, x])
                               for x, y in search_positions)
                     for name, parent_name in CONTROL_PARENT.items()}
    enabled = []
    for name in METHODS:
        base = (BASE_IDENTITY if FAMILY[name] == 'identity' else
                BASE_RECTANGULAR if FAMILY[name] == 'rectangular20' else 'gaussian')
        if name in oracle_enabled and parent_active[base]:
            enabled.append(name)
    for name in METHODS:
        base = (BASE_IDENTITY if FAMILY[name] == 'identity' else
                BASE_RECTANGULAR if FAMILY[name] == 'rectangular20' else 'gaussian')
        annular_support[name]['parent_family_support'] = parent_active[base]
    if task['phase'] == 'positive' or trial['name'] in {site['name'] for site in parent['sites']}:
        inner.full.radial.require(set(enabled) == set(METHODS),
                                  'an injection site lacks complete smoothing-grid support')
    enabled_indices = [METHODS.index(name) for name in enabled]
    header['HCI FILTER LABELS'] = ','.join(METHODS)
    header['HCI RADIAL BINS'] = ','.join(map(str, bins))
    fits.writeto(directory/'amplitudes.fits', maps, header)
    active_header = header.copy()
    active_header['HCI FILTER LABELS'] = ','.join(enabled)
    fits.writeto(directory/'active_amplitudes.fits', maps[enabled_indices], active_header)
    command = []
    if enabled:
        delta_x, delta_y = trial['row']-127.5, trial['column']-127.5
        command = [str(study/'software/hciAnalyze'),
            '--file='+str(directory/'active_amplitudes.fits'), '--lambdaD=3.6',
            '--planet.sep=11.782,'+str(math.hypot(delta_x, delta_y)),
            '--planet.PA=262.051,'+str(math.degrees(-math.atan2(delta_x, delta_y)) % 360),
            '--planet.R=7.3,'+str(protocol['annular_trial_exclusion_radius']),
            '--snr.apertureR=60', '--snr.minRad=0', '--snr.maxRad=60',
            '--filter.psfResponse=', '--filter.lpfGaussFW=0', '--filter.hpfGaussFW=0',
            '--noise.model=identity', '--noise.only=false', '--noise.outputDiagnostics=false']
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
    else:
        (directory/'analysis.log').write_text(
            'No method has complete five-pixel support; production annular analysis skipped.\n')
        active_snr = np.empty((0, *science.shape), dtype='f4')
        snr_header = active_header.copy()
        snr_header['SNRMEAN'] = snr_header['SNRSMALL'] = 1
        fits.writeto(directory/'active_amplitudes_snr.fits', active_snr, snr_header)
    write_json(directory/'command.json', command)
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
    snr_control_errors = {}
    for method, parent_method in CONTROL_PARENT.items():
        current = snr[METHODS.index(method)]
        prior = parent_snr[PARENT_INDEX[parent_method]]
        snr_control_errors[method] = maximum_finite_difference(current, prior)
        inner.full.radial.require(np.allclose(current, prior, rtol=2e-6, atol=2e-6, equal_nan=True),
                                  'unsmoothed SNR control changed for '+method)
    measurement = {'task': task, 'trial': trial, 'source': fingerprint(source),
        'models': score_trial(trial, snr, maps), 'active_methods': enabled,
        'annular_support': annular_support, 'required_bins': bins, 'fit_details': details,
        'annular_profile_min_pixels': {name: min(row['pixels'] for row in profile
            if int(row['radius']-.5) in bins) for name, profile in profiles.items()},
        'fitted_pixels': fitted_pixels, 'elapsed_seconds': time.monotonic()-start,
        'annular_oracle_max_errors': oracle_errors,
        'amplitude_control_max_errors': amplitude_control_errors,
        'snr_control_max_errors': snr_control_errors}
    write_json(directory/'profiles.json', {name: [row for row in profile
        if int(row['radius']-.5) in bins] for name, profile in profiles.items()})
    write_json(directory/'measurements.json', measurement)
    products = [fingerprint(directory/name) for name in ('amplitudes.fits', 'active_amplitudes.fits',
        'active_amplitudes_snr.fits', 'amplitudes_snr.fits', 'measurements.json', 'profiles.json',
        'command.json')]
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
                'finished': len(paths), 'total': len(tasks), 'last': futures[future],
                'new_reductions': 0})
            print(f'{phase}: {len(paths)}/{len(tasks)} {futures[future]}', flush=True)
    return [inner.full.read(Path(path)) for path in sorted(paths)]


def parent_pool_method(method: str) -> str:
    """Map each smoothing arm to its frozen parent calibration family."""
    if FAMILY[method] == 'identity':
        return 'identity'
    if FAMILY[method] == 'rectangular20':
        return PARENT_RECTANGULAR
    return 'gaussian'


def calibrate(root: Path, baseline: list[dict]) -> tuple[dict, dict]:
    """Freeze thresholds on identical locations within each weighting family."""
    record = root/'calibration_complete.json'
    if record.exists():
        inner.full.verify(inner.full.read(record)['products'])
        return (inner.full.read(root/'thresholds.json'),
                inner.full.read(root/'effective_calibration_pools.json')['pools'])
    protocol = inner.full.read(root/'protocol.json')
    comparison = Path(protocol['parent_comparison'])
    parent_pools = inner.full.read(comparison/'effective_calibration_pools.json')['pools']
    parent_thresholds = inner.full.read(comparison/'thresholds.json')
    rows = {row['trial']['name']: row for row in baseline}
    pools, thresholds = {}, {}
    for nominal in RADII:
        pools[str(nominal)], thresholds[str(nominal)] = {}, {}
        for method in METHODS:
            parent_method = parent_pool_method(method)
            pool = parent_pools[str(nominal)][parent_method]
            inner.full.radial.require(pool is not None and len(pool['trials']) == CALIBRATION_SEARCHES,
                                      f'missing frozen calibration pool for radius {nominal} {method}')
            pools[str(nominal)][method] = {**pool, 'source_parent_method': parent_method}
            values = [rows[name]['models'][method] for name in pool['trials']]
            inner.full.radial.require(all(value['valid'] for value in values),
                                      f'invalid frozen calibration pool for radius {nominal} {method}')
            thresholds[str(nominal)][method] = max(value['search_score'] for value in values)
        inner.full.radial.require(np.isclose(thresholds[str(nominal)][BASE_IDENTITY],
                                  parent_thresholds[str(nominal)]['identity'], rtol=2e-6, atol=2e-6),
                                  'identity threshold control changed')
        inner.full.radial.require(np.isclose(thresholds[str(nominal)][BASE_RECTANGULAR],
                                  parent_thresholds[str(nominal)][PARENT_RECTANGULAR], rtol=2e-6, atol=2e-6),
                                  'rectangular threshold control changed')
        inner.full.radial.require(np.isclose(thresholds[str(nominal)]['gaussian'],
                                  parent_thresholds[str(nominal)]['gaussian'], rtol=2e-6, atol=2e-6),
                                  'Gaussian threshold control changed')
    write_json(root/'baseline.json', {'rows': baseline, 'parent_baseline_reused': True})
    write_json(root/'effective_calibration_pools.json', {'schema': 1,
        'score_values_used_for_selection': False, 'locations_reused_from_parent': True,
        'paired_within_weight_family': True, 'pools': pools})
    write_json(root/'thresholds.json', thresholds)
    products = [fingerprint(root/name) for name in
                ('baseline.json', 'effective_calibration_pools.json', 'thresholds.json')]
    write_json(record, {'thresholds_frozen_before_positive_analysis': True, 'products': products})
    return thresholds, pools


def center_detail(record: dict, method: str) -> dict | None:
    """Return one trial center's saved template diagnostic when applicable."""
    key = f'{record["trial"]["row"]},{record["trial"]["column"]}'
    return record['fit_details'].get(key, {}).get(method)


def summarize(root: Path, baseline: list[dict], positives: list[dict],
              thresholds: dict, pools: dict) -> None:
    """Record recovery, SNR, throughput, and paired smoothing changes."""
    protocol = inner.full.read(root/'protocol.json')
    comparison = Path(protocol['parent_comparison'])
    study = Path(protocol['parent_study'])
    parent = inner.full.read(study/'protocol.json')
    parent_results = inner.full.read(comparison/'results.json')
    levels = tuple(protocol['brightnesses'])
    labels = protocol['level_labels']
    baseline_rows = {record['trial']['name']: record for record in baseline}
    nulls, measurements = [], []
    for site in parent['sites']:
        record = baseline_rows[site['name']]
        for method, value in record['models'].items():
            value['detected'] = value['valid'] and value['search_score'] > thresholds[
                str(site['nominal_radius'])][method]
        nulls.append(record)
    for record in positives:
        job = record['task']['job']
        baseline_record = baseline_rows[job['site']]
        for method, value in record['models'].items():
            value['detected'] = value['valid'] and value['search_score'] > thresholds[
                str(job['nominal_radius'])][method]
            if method != 'gaussian' and value['center_amplitude'] is not None:
                baseline_amplitude = baseline_record['models'][method]['center_amplitude']
                value['paired_amplitude_increment'] = value['center_amplitude']-baseline_amplitude
                value['paired_throughput'] = value['paired_amplitude_increment']/job['contrast']
                detail = center_detail(record, method)
                value['expected_original_response'] = (detail or {}).get('expected_original_response')
                value['template_cosine'] = (detail or {}).get('template_cosine')
                value['matched_snr_efficiency_if_original_exact'] = (detail or {}).get(
                    'matched_snr_efficiency_if_original_exact')
        measurements.append(record)
    parent_nulls = {record['trial']['name']: record for record in parent_results['nulls']}
    parent_measurements = {record['task']['job']['name']: record
                           for record in parent_results['measurements']}
    for record in nulls:
        for method, parent_method in CONTROL_PARENT.items():
            inner.full.radial.require(record['models'][method]['detected'] ==
                                      parent_nulls[record['trial']['name']]['models'][parent_method]['detected'],
                                      'parent null decision changed for '+method)
    for record in measurements:
        for method, parent_method in CONTROL_PARENT.items():
            inner.full.radial.require(record['models'][method]['detected'] ==
                                      parent_measurements[record['task']['job']['name']]['models'][parent_method]['detected'],
                                      'parent positive decision changed for '+method)
    groups = []
    for nominal in RADII:
        for method in METHODS:
            site_rows = [record for record in nulls if record['trial']['nominal_radius'] == nominal]
            one = {'radius': nominal, 'method': method,
                'null_exceedances': sum(record['models'][method]['detected'] for record in site_rows),
                'invalid_nulls': sum(not record['models'][method]['valid'] for record in site_rows),
                'levels': []}
            for level in levels:
                selected = [record for record in measurements
                    if record['task']['job']['nominal_radius'] == nominal and
                    record['task']['job']['brightness_multiplier'] == level]
                values = [record['models'][method] for record in selected]
                valid = [value for value in values if value['valid']]
                throughputs = [value['paired_throughput'] for value in valid
                               if 'paired_throughput' in value]
                predicted = [value['expected_original_response'] for value in valid
                             if value.get('expected_original_response') is not None]
                support = [record['annular_profile_min_pixels'][method] for record in selected
                           if record['models'][method]['valid']]
                one['levels'].append({'brightness_multiplier': level, 'trials': len(values),
                    'detections': sum(value['detected'] for value in values),
                    'invalid_searches': sum(not value['valid'] for value in values),
                    'valid_searches': len(valid),
                    'mean_search_snr': float(np.mean([value['search_score'] for value in valid])),
                    'mean_center_snr': float(np.mean([value['snr_pixels'][0] for value in valid])),
                    'mean_paired_throughput': float(np.mean(throughputs)) if throughputs else None,
                    'mean_expected_original_response': float(np.mean(predicted)) if predicted else None,
                    'minimum_annular_pixels': min(support) if support else None,
                    'median_minimum_annular_pixels': float(np.median(support)) if support else None})
            groups.append(one)
    aggregate = []
    for method in METHODS:
        method_groups = [group for group in groups if group['method'] == method]
        aggregate_levels = []
        for level in levels:
            values = [record['models'][method] for record in measurements
                      if record['task']['job']['brightness_multiplier'] == level]
            valid = [value for value in values if value['valid']]
            throughputs = [value['paired_throughput'] for value in valid
                           if 'paired_throughput' in value]
            efficiencies = [value['matched_snr_efficiency_if_original_exact'] for value in valid
                            if value.get('matched_snr_efficiency_if_original_exact') is not None]
            aggregate_levels.append({'brightness_multiplier': level, 'trials': len(values),
                'valid_searches': len(valid), 'detections': sum(value['detected'] for value in values),
                'mean_search_snr': float(np.mean([value['search_score'] for value in valid])),
                'mean_center_snr': float(np.mean([value['snr_pixels'][0] for value in valid])),
                'mean_paired_throughput': float(np.mean(throughputs)) if throughputs else None,
                'mean_matched_snr_efficiency_if_original_exact':
                    float(np.mean(efficiencies)) if efficiencies else None})
        aggregate.append({'method': method,
            'valid_sites': min(level['valid_searches'] for level in aggregate_levels),
            'null_exceedances': sum(group['null_exceedances'] for group in method_groups),
            'invalid_nulls': sum(group['invalid_nulls'] for group in method_groups),
            'levels': aggregate_levels})
    paired = []
    for family, variants in (('identity', IDENTITY), ('rectangular20', RECTANGULAR)):
        base = variants[0.0]
        for fwhm in FWHMS[1:]:
            tested = variants[fwhm]
            for level in levels:
                selected = [record for record in measurements
                            if record['task']['job']['brightness_multiplier'] == level]
                paired.append({'family': family, 'fwhm': fwhm, 'brightness_multiplier': level,
                    'smoothed_only': [record['task']['job']['name'] for record in selected
                        if record['models'][tested]['detected'] and not record['models'][base]['detected']],
                    'unsmoothed_only': [record['task']['job']['name'] for record in selected
                        if record['models'][base]['detected'] and not record['models'][tested]['detected']]})
    verification = {'new_reductions': 0, 'unsmoothed_controls_reproduced': True,
        'calibration_locations_paired_within_family': True,
        'maximum_control_amplitude_difference': max(value for record in [*baseline, *measurements]
            for value in record['amplitude_control_max_errors'].values() if value is not None),
        'maximum_control_snr_difference': max(value for record in [*baseline, *measurements]
            for value in record['snr_control_max_errors'].values() if value is not None),
        'maximum_annular_oracle_difference': max(value for record in [*baseline, *measurements]
            for value in record['annular_oracle_max_errors'].values() if value is not None)}
    write_json(root/'results.json', {'groups': groups, 'aggregate': aggregate,
        'paired_decisions': paired, 'nulls': nulls, 'measurements': measurements,
        'thresholds': thresholds, 'calibration_pools': pools, 'verification': verification,
        'caveats': [protocol['dependence'],
            'The smoothing grid and inspected injections make this a development comparison.',
            'Gaussian filtering assumes the measured 11x11 response is zero outside its stored support.']})
    lines = ['# Measured-response smoothing under identity and rectangular-PSD weights', '',
        'Only the measured 11x11 response template changes. The completed SNR-3/5/7 images, one-lambda/D '
        'trial exclusion, covariance samples, candidate data, searches, and production annular SNR are reused. '
        'Every method has its own threshold on fixed locations shared within its weighting family.', '',
        '## Aggregate recovery and measured SNR', '',
        '| Method | Valid sites | '+' | '.join(labels)+' detections | Nulls | '
        + ' | '.join(labels)+' mean max SNR |',
        '| --- | ---: | '+' | '.join('---:' for _ in levels)+' | ---: | '
        + ' | '.join('---:' for _ in levels)+' |']
    for row in aggregate:
        lines.append(f'| {LABELS[row["method"]]} | {row["valid_sites"]}/36 | '
            + ' | '.join(str(level['detections']) for level in row['levels'])
            + f' | {row["null_exceedances"]}/36 | '
            + ' | '.join(f'{level["mean_search_snr"]:.4f}' for level in row['levels'])+' |')
    lines += ['', '## Mean paired amplitude throughput and exact-template efficiency', '',
        'Each value is the positive-minus-baseline center amplitude divided by injected contrast. It measures '
        'response mismatch separately from annular SNR. The parenthesized value is the predicted matched-SNR '
        'ratio if the unsmoothed response were exact. The Gaussian image is not a contrast estimator.', '',
        '| Method | '+' | '.join(labels)+' throughput (exact-template efficiency) |',
        '| --- | '+' | '.join('---:' for _ in levels)+' |']
    for row in aggregate:
        if row['method'] == 'gaussian':
            continue
        lines.append(f'| {LABELS[row["method"]]} | '+' | '.join(
            f'{level["mean_paired_throughput"]:.4f} '
            f'({level["mean_matched_snr_efficiency_if_original_exact"]:.4f})'
            for level in row['levels'])+' |')
    lines += ['', '## Paired decision changes from the unsmoothed response', '']
    for row in paired:
        label = labels[levels.index(row['brightness_multiplier'])]
        lines.append(f'- {row["family"]}, FWHM {row["fwhm"]:.1f} px, {label}: '
            f'smoothed only {len(row["smoothed_only"])}; unsmoothed only {len(row["unsmoothed_only"])}.')
    lines += ['', 'All thresholds were frozen before positive analysis. No P4 reduction was run.', '']
    (root/'results.md').write_text('\n'.join(lines))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    x_values = protocol['target_source_snrs']
    for column, (family, variants) in enumerate((('Identity', IDENTITY),
                                                  ('Rectangular ±20', RECTANGULAR))):
        for fwhm, method in variants.items():
            row = next(item for item in aggregate if item['method'] == method)
            label = f'{fwhm:.1f} px'
            axes[0, column].plot(x_values, [level['mean_search_snr'] for level in row['levels']],
                                 'o-', label=label)
            axes[1, column].plot(x_values, [level['detections']/36 for level in row['levels']],
                                 'o-', label=label)
        axes[0, column].set(title=family, xlabel='Nominal identity source SNR',
                            ylabel='Mean five-pixel maximum SNR')
        axes[1, column].set(xlabel='Nominal identity source SNR', ylabel='Recovery fraction',
                            ylim=(-.03, 1.03))
        axes[0, column].legend(title='Response LPF FWHM')
    fig.suptitle('Measured-response smoothing with fixed identity and PSD weighting')
    fig.savefig(root/'comparison.png', dpi=170)
    plt.close(fig)


def self_check() -> None:
    """Verify smoothing, identity normalization, and PSD unit-response algebra."""
    rng = np.random.default_rng(20260920)
    template = rng.normal(size=121)
    inner.full.radial.require(np.array_equal(smooth_template(template, 0), template),
                              'zero-width smoothing is not exact')
    for fwhm in FWHMS[1:]:
        filtered = smooth_template(template, fwhm)
        diagnostic = template_diagnostic(filtered, template)
        inner.full.radial.require(np.isfinite(filtered).all() and
                                  0 < diagnostic['template_energy_fraction'] < 1 and
                                  0 < diagnostic['template_cosine'] <= 1,
                                  'invalid smoothed-template diagnostic')
    data = rng.normal(size=121)
    energy = float(template @ template)
    inner.full.radial.require(np.isclose((template/energy) @ template, 1,
                              rtol=1e-14, atol=1e-14), 'identity unit response failed')
    samples = rng.normal(size=(48, 121))
    model = inner.full.psd.fit_psd(samples, 'rectangular', .3)
    inner.full.radial.require(model is not None, 'synthetic PSD fit failed')
    result = inner.full.radial.filter_stamp(data, template, model, np.ones(121))
    inner.full.radial.require(np.isclose(result['physical_weight'] @ template, 1,
                              rtol=1e-12, atol=1e-12), 'PSD unit response failed')
    print('response-smoothing checks passed', flush=True)


def run(root: Path) -> None:
    """Run calibration, then analyze all saved positives without new reductions."""
    manifest = inner.full.read(root/'manifest.json')
    inner.full.verify([*manifest['inputs'], manifest['protocol']])
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (root/'complete.json').exists():
            inner.full.verify(inner.full.read(root/'complete.json')['products'])
            print('Response-smoothing comparison already complete.', flush=True)
            return
        protocol = inner.full.read(root/'protocol.json')
        try:
            write_json(root/'state.json', {'status': 'baseline', 'pid': os.getpid(),
                'new_reductions': 0})
            baseline = run_phase(root, protocol['baseline_tasks'], 'baseline')
            thresholds, pools = calibrate(root, baseline)
            calibration = inner.full.read(root/'calibration_complete.json')
            write_json(root/'state.json', {'status': 'positive', 'pid': os.getpid(),
                'new_reductions': 0})
            positives = run_phase(root, protocol['positive_tasks'], 'positive')
            summarize(root, baseline, positives, thresholds, pools)
            inner.full.verify([*manifest['inputs'], manifest['protocol'], *calibration['products']])
            products = [fingerprint(root/name) for name in ('results.json', 'results.md',
                'comparison.png', 'thresholds.json', 'effective_calibration_pools.json', 'baseline.json')]
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
    """Expose synthetic validation, immutable preparation, and resumable execution."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('check', 'prepare', 'run'))
    parser.add_argument('--root', type=Path)
    parser.add_argument('--comparison', type=Path,
                        default=Path('working/roc/p4_inner_snr357_patch_rms_lambdad_20260919'))
    parser.add_argument('--cpus', type=int, nargs='+', default=list(range(12)))
    args = parser.parse_args()
    if args.action == 'check':
        self_check()
        return
    inner.full.radial.require(args.root is not None, '--root is required for prepare and run')
    if args.action == 'prepare':
        prepare(args)
    else:
        run(args.root.resolve())


if __name__ == '__main__':
    main()
