#!/usr/bin/env python3
"""Compare truncated and clipped inverses of the rectangular PSD covariance.

Reuse the completed inner SNR-3/5/7 images and lambda/D-masked annular
analysis. Keep the measured response fixed, regularize only the eigenspectrum
of the finite PSD covariance, and freeze all null thresholds before reading
positive measurements.
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
from scipy.linalg import eigh

import compare_p4_step5_response_smoothing as previous
from measure_p4_step5_planet import annular_oracle
import run_p4_step5_inner_rectangular as inner
from run_p4_step5_full_injections import fingerprint, write_json

CUTOFFS = (0.5, 0.75, 1.0)
LAMBDA_D = 3.6
RADII = inner.RADII
RECTANGULAR_WIDTH = 20
CALIBRATION_SEARCHES = inner.CALIBRATION_SEARCHES


def cutoff_tag(cutoff: float) -> str:
    """Return a stable method-name suffix for one variance cutoff."""
    return ('%.2f' % cutoff).rstrip('0').rstrip('.').replace('.', 'p')


FULL = 'rect20_full'
CLIPPED = {cutoff: 'rect20_clip_'+cutoff_tag(cutoff) for cutoff in CUTOFFS}
TRUNCATED = {cutoff: 'rect20_tsvd_'+cutoff_tag(cutoff) for cutoff in CUTOFFS}
REGULARIZED = (*CLIPPED.values(), *TRUNCATED.values())
PRECISION_METHODS = (FULL, *REGULARIZED)
REFERENCES = ('identity', 'identity_response_lpf1p8', 'rect20_response_lpf2p7', 'gaussian')
METHODS = (FULL, *CLIPPED.values(), *TRUNCATED.values(), *REFERENCES)
LABELS = {FULL: 'Rectangular ±20, full inverse',
          **{name: f'Rectangular ±20, clipped at {cutoff:g} mean variance'
             for cutoff, name in CLIPPED.items()},
          **{name: f'Rectangular ±20, truncated at {cutoff:g} mean variance'
             for cutoff, name in TRUNCATED.items()},
          'identity': 'Identity', 'identity_response_lpf1p8': 'Identity, response LPF 1.8 px',
          'rect20_response_lpf2p7': 'Rectangular ±20, response LPF 2.7 px',
          'gaussian': 'Gaussian FWHM 3.6'}
FAMILY = {FULL: 'full', **{name: 'clipped' for name in CLIPPED.values()},
          **{name: 'truncated' for name in TRUNCATED.values()},
          **{name: 'reference' for name in REFERENCES}}
CONTROL_PARENT = {FULL: previous.RECTANGULAR[0.0], 'identity': previous.IDENTITY[0.0],
                  'identity_response_lpf1p8': previous.IDENTITY[1.8],
                  'rect20_response_lpf2p7': previous.RECTANGULAR[2.7], 'gaussian': 'gaussian'}
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
    """Freeze the paired saved-image PSD-precision comparison."""
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
                              prior.get('response_smoothing_fwhm_pixels') == [0.0, 0.9, 1.8, 2.7, 3.6],
                              'parent must be the completed response-smoothing comparison')
    inner.full.radial.require(parent.get('target_source_snrs') == [3.0, 5.0, 7.0] and
                              len(parent['sites']) == 36 and len(parent['calibration_trials']) == 128,
                              'changed SNR-3/5/7 parent design')
    records = comparison_inputs(comparison, prior)
    inner.full.verify(records)
    protocol = {'schema': 1,
        'purpose': 'hard-truncated and clipped inverses of raw rectangular ±20 PSD covariance',
        'parent_comparison': str(comparison), 'parent_study': str(study),
        'methods': list(METHODS), 'labels': LABELS, 'families': FAMILY,
        'variance_cutoffs_over_mean': list(CUTOFFS),
        'full_inverse': 'all 121 eigenmodes weighted by reciprocal covariance eigenvalue',
        'clipped_inverse': 'replace each covariance eigenvalue below cutoff times mean variance by that cutoff before inversion; retain all modes',
        'truncated_inverse': 'set inverse weight to zero for covariance eigenvalues below cutoff times mean variance; invert every retained mode exactly',
        'rectangular_psd': 'raw patches, ±20-pixel training centers, ensemble-mean subtraction, rectangular 21x21 padded PSD, mixture0.3',
        'fixed_response': 'use the unchanged unsmoothed measured 11x11 response for every new precision policy',
        'precision_only_change': 'candidate images, response, covariance estimate, fitted mean, exclusions, annular SNR, search geometry, and source injections are unchanged',
        'references': 'copy and exactly reproduce full rectangular inverse, identity, identity response-LPF1.8, rectangular response-LPF2.7, and production Gaussian from the parent comparison',
        'radii': list(RADII), 'brightnesses': parent['brightness_multipliers'],
        'target_source_snrs': parent['target_source_snrs'], 'level_labels': parent['level_labels'],
        'baseline_tasks': prior['baseline_tasks'], 'positive_tasks': prior['positive_tasks'],
        'calibration_pools': 'reuse frozen parent locations; every new precision method uses the raw rectangular ±20 pool and receives its own maximum-null threshold',
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


def precision_result(data: np.ndarray, template: np.ndarray, mean: np.ndarray,
                     covariance: np.ndarray, eigenvectors: np.ndarray, precision: np.ndarray,
                     full_variance: float) -> dict:
    """Apply one regularized covariance precision with unit response."""
    coefficients = eigenvectors.T @ template
    precision_template = eigenvectors @ (precision * coefficients)
    energy = float(template @ precision_template)
    inner.full.radial.require(np.isfinite(energy) and energy > 0,
                              'regularized precision has no template support')
    weight = precision_template/energy
    amplitude = float(weight @ (data-mean))
    modeled_variance = float(weight @ covariance @ weight)
    efficiency = math.sqrt(full_variance/modeled_variance)
    inner.full.radial.require(np.isclose(weight @ template, 1, rtol=1e-11, atol=1e-12) and
                              0 < efficiency <= 1+1e-10,
                              'invalid regularized precision normalization')
    return {'amplitude': amplitude, 'conditional_sigma': 1/math.sqrt(energy),
            'psd_sigma': math.sqrt(modeled_variance),
            'efficiency_if_psd_exact': efficiency,
            'inverse_energy': energy, 'squared_weight_norm': float(weight @ weight)}


def fit_pixel(science: np.ndarray, position: tuple[int, int], template: np.ndarray,
              forbidden: np.ndarray) -> tuple[dict, dict]:
    """Evaluate clipped and hard-truncated inverses of one fitted PSD covariance."""
    x, y = position
    data = science[y-5:y+6, x-5:x+6].ravel()
    rings = inner.full.radial.geometry(science, position, forbidden)
    matrices = {offset: inner.full.radial.extract(science, ring) for offset, ring in rings.items()
                if abs(offset) <= RECTANGULAR_WIDTH}
    samples = np.vstack(list(matrices.values()))
    model = inner.full.psd.fit_psd(samples, 'rectangular', .3)
    amplitudes, diagnostics = {}, {}
    for name in (FULL, *CLIPPED.values(), *TRUNCATED.values()):
        diagnostics[name] = {'valid': model is not None, 'samples': len(samples)}
    if model is None:
        return amplitudes, diagnostics
    eigenvalues, eigenvectors = eigh(model['covariance'], check_finite=False)
    inner.full.radial.require(np.all(eigenvalues > 0) and
                              np.isclose(eigenvalues.mean(), model['target_variance'],
                                         rtol=1e-12, atol=0),
                              'PSD eigensystem changed covariance trace or positivity')
    full = inner.full.radial.filter_stamp(data, template, model, np.ones(121))
    full_variance = float(full['sigma'])**2
    eigen_full = precision_result(data, template, model['mean'], model['covariance'],
                                  eigenvectors, 1/eigenvalues, full_variance)
    inner.full.radial.require(np.isclose(eigen_full['amplitude'], full['amplitude'],
                              rtol=1e-10, atol=1e-12) and
                              np.isclose(eigen_full['psd_sigma'], full['sigma'], rtol=1e-10, atol=0),
                              'full eigen-inverse differs from Cholesky control')
    diagnostics[FULL].update(eigen_full, condition_number=float(eigenvalues[-1]/eigenvalues[0]),
                             minimum_over_mean=float(eigenvalues[0]/model['target_variance']),
                             maximum_over_mean=float(eigenvalues[-1]/model['target_variance']),
                             retained_modes=121, modified_modes=0,
                             full_matched_energy_fraction=1.0)
    for cutoff in CUTOFFS:
        threshold = cutoff*model['target_variance']
        clipped = eigenvalues < threshold
        clipped_precision = 1/np.maximum(eigenvalues, threshold)
        clipped_result = precision_result(data, template, model['mean'], model['covariance'],
                                          eigenvectors, clipped_precision, full_variance)
        clipped_result.update(cutoff_over_mean=cutoff, retained_modes=121,
                              modified_modes=int(clipped.sum()),
                              full_matched_energy_fraction=float(clipped_result['inverse_energy']/
                                                                 eigen_full['inverse_energy']))
        name = CLIPPED[cutoff]
        amplitudes[name] = clipped_result['amplitude']
        diagnostics[name].update(clipped_result)
        retained = ~clipped
        inner.full.radial.require(retained.any(), 'hard truncation discards every covariance mode')
        truncated_precision = np.where(retained, 1/eigenvalues, 0)
        truncated_result = precision_result(data, template, model['mean'], model['covariance'],
                                            eigenvectors, truncated_precision, full_variance)
        truncated_result.update(cutoff_over_mean=cutoff, retained_modes=int(retained.sum()),
                                modified_modes=int(clipped.sum()),
                                full_matched_energy_fraction=float(truncated_result['inverse_energy']/
                                                                   eigen_full['inverse_energy']))
        name = TRUNCATED[cutoff]
        amplitudes[name] = truncated_result['amplitude']
        diagnostics[name].update(truncated_result)
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


def support_control(method: str) -> str:
    """Return the copied parent method that defines one method's usable support."""
    return FULL if method in PRECISION_METHODS else method


def analyze(task: dict) -> str:
    """Build regularized-precision maps and reproduce every copied parent control."""
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
    for method, parent_method in CONTROL_PARENT.items():
        maps[METHODS.index(method)] = parent_amplitudes[PARENT_INDEX[parent_method]]
    details = {}
    fitted_pixels = 0
    full_solver_errors = []
    for (x, y), template in sorted(templates.items()):
        if not selected[y, x] or not clean_centers[y, x]:
            continue
        data = science[y-5:y+6, x-5:x+6]
        if data.size != 121 or not np.isfinite(data).all():
            continue
        amplitudes, diagnostics = fit_pixel(science, (x, y), template, forbidden)
        for name, amplitude in amplitudes.items():
            maps[METHODS.index(name), y, x] = amplitude
        recomputed = diagnostics[FULL].get('amplitude')
        copied = float(maps[METHODS.index(FULL), y, x])
        inner.full.radial.require((recomputed is not None) == np.isfinite(copied),
                                  'full-inverse validity differs from copied parent control')
        if recomputed is not None:
            full_solver_errors.append(abs(recomputed-copied))
            inner.full.radial.require(np.isclose(recomputed, copied, rtol=2e-6, atol=1e-12),
                                      'eigen full inverse differs from copied parent control')
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
                                  'copied amplitude control changed for '+method)
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
        base = support_control(name)
        if name in oracle_enabled and parent_active[base]:
            enabled.append(name)
    for name in METHODS:
        base = support_control(name)
        annular_support[name]['parent_family_support'] = parent_active[base]
    if task['phase'] == 'positive' or trial['name'] in {site['name'] for site in parent['sites']}:
        inner.full.radial.require(set(enabled) == set(METHODS),
                                  'an injection site lacks complete precision-grid support')
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
                                  'copied SNR control changed for '+method)
    measurement = {'task': task, 'trial': trial, 'source': fingerprint(source),
        'models': score_trial(trial, snr, maps), 'active_methods': enabled,
        'annular_support': annular_support, 'required_bins': bins, 'fit_details': details,
        'annular_profile_min_pixels': {name: min(row['pixels'] for row in profile
            if int(row['radius']-.5) in bins) for name, profile in profiles.items()},
        'fitted_pixels': fitted_pixels, 'elapsed_seconds': time.monotonic()-start,
        'annular_oracle_max_errors': oracle_errors,
        'amplitude_control_max_errors': amplitude_control_errors,
        'snr_control_max_errors': snr_control_errors,
        'maximum_full_solver_amplitude_difference': max(full_solver_errors, default=None)}
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
    """Map each precision arm or reference to its frozen parent calibration pool."""
    if method in PRECISION_METHODS or method == 'rect20_response_lpf2p7':
        return previous.RECTANGULAR[0.0]
    if method in ('identity', 'identity_response_lpf1p8'):
        return previous.IDENTITY[0.0]
    return 'gaussian'


def calibrate(root: Path, baseline: list[dict]) -> tuple[dict, dict]:
    """Freeze thresholds on the raw-PSD pool and reproduce copied controls."""
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
        for method, parent_method in CONTROL_PARENT.items():
            inner.full.radial.require(np.isclose(thresholds[str(nominal)][method],
                                      parent_thresholds[str(nominal)][parent_method],
                                      rtol=2e-6, atol=2e-6),
                                      'copied threshold control changed for '+method)
    write_json(root/'baseline.json', {'rows': baseline, 'parent_baseline_reused': True})
    write_json(root/'effective_calibration_pools.json', {'schema': 1,
        'score_values_used_for_selection': False, 'locations_reused_from_parent': True,
        'new_precision_methods_paired_on_raw_rectangular_pool': True, 'pools': pools})
    write_json(root/'thresholds.json', thresholds)
    products = [fingerprint(root/name) for name in
                ('baseline.json', 'effective_calibration_pools.json', 'thresholds.json')]
    write_json(record, {'thresholds_frozen_before_positive_analysis': True, 'products': products})
    return thresholds, pools


def center_detail(record: dict, method: str) -> dict | None:
    """Return one trial center's saved precision diagnostic when applicable."""
    key = f'{record["trial"]["row"]},{record["trial"]["column"]}'
    return record['fit_details'].get(key, {}).get(method)


def summarize(root: Path, baseline: list[dict], positives: list[dict],
              thresholds: dict, pools: dict) -> None:
    """Record recovery, SNR, throughput, and paired precision-policy changes."""
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
                if detail is not None:
                    value['efficiency_if_psd_exact'] = detail['efficiency_if_psd_exact']
                    value['retained_modes'] = detail['retained_modes']
                    value['modified_modes'] = detail['modified_modes']
                    value['full_matched_energy_fraction'] = detail[
                        'full_matched_energy_fraction']
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
                efficiencies = [value['efficiency_if_psd_exact'] for value in valid
                                if value.get('efficiency_if_psd_exact') is not None]
                retained = [value['retained_modes'] for value in valid
                            if value.get('retained_modes') is not None]
                modified = [value['modified_modes'] for value in valid
                            if value.get('modified_modes') is not None]
                support = [record['annular_profile_min_pixels'][method] for record in selected
                           if record['models'][method]['valid']]
                one['levels'].append({'brightness_multiplier': level, 'trials': len(values),
                    'detections': sum(value['detected'] for value in values),
                    'invalid_searches': sum(not value['valid'] for value in values),
                    'valid_searches': len(valid),
                    'mean_search_snr': float(np.mean([value['search_score'] for value in valid])),
                    'mean_center_snr': float(np.mean([value['snr_pixels'][0] for value in valid])),
                    'mean_paired_throughput': float(np.mean(throughputs)) if throughputs else None,
                    'mean_efficiency_if_psd_exact': float(np.mean(efficiencies))
                        if efficiencies else None,
                    'mean_retained_modes': float(np.mean(retained)) if retained else None,
                    'mean_modified_modes': float(np.mean(modified)) if modified else None,
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
            efficiencies = [value['efficiency_if_psd_exact'] for value in valid
                            if value.get('efficiency_if_psd_exact') is not None]
            retained = [value['retained_modes'] for value in valid
                        if value.get('retained_modes') is not None]
            modified = [value['modified_modes'] for value in valid
                        if value.get('modified_modes') is not None]
            aggregate_levels.append({'brightness_multiplier': level, 'trials': len(values),
                'valid_searches': len(valid), 'detections': sum(value['detected'] for value in values),
                'mean_search_snr': float(np.mean([value['search_score'] for value in valid])),
                'mean_center_snr': float(np.mean([value['snr_pixels'][0] for value in valid])),
                'mean_paired_throughput': float(np.mean(throughputs)) if throughputs else None,
                'mean_efficiency_if_psd_exact': float(np.mean(efficiencies))
                    if efficiencies else None,
                'mean_retained_modes': float(np.mean(retained)) if retained else None,
                'mean_modified_modes': float(np.mean(modified)) if modified else None})
        aggregate.append({'method': method,
            'valid_sites': min(level['valid_searches'] for level in aggregate_levels),
            'null_exceedances': sum(group['null_exceedances'] for group in method_groups),
            'invalid_nulls': sum(group['invalid_nulls'] for group in method_groups),
            'levels': aggregate_levels})
    paired = []
    for policy, variants in (('clipped', CLIPPED), ('hard_truncated', TRUNCATED)):
        for cutoff, tested in variants.items():
            for level in levels:
                selected = [record for record in measurements
                            if record['task']['job']['brightness_multiplier'] == level]
                paired.append({'policy': policy, 'cutoff_over_mean': cutoff,
                    'brightness_multiplier': level,
                    'regularized_only': [record['task']['job']['name'] for record in selected
                        if record['models'][tested]['detected'] and
                        not record['models'][FULL]['detected']],
                    'full_inverse_only': [record['task']['job']['name'] for record in selected
                        if record['models'][FULL]['detected'] and
                        not record['models'][tested]['detected']]})
    solver_errors = [record['maximum_full_solver_amplitude_difference']
                     for record in [*baseline, *measurements]
                     if record['maximum_full_solver_amplitude_difference'] is not None]
    verification = {'new_reductions': 0, 'copied_controls_reproduced': True,
        'calibration_locations_paired_on_raw_rectangular_pool': True,
        'maximum_control_amplitude_difference': max(value for record in [*baseline, *measurements]
            for value in record['amplitude_control_max_errors'].values() if value is not None),
        'maximum_control_snr_difference': max(value for record in [*baseline, *measurements]
            for value in record['snr_control_max_errors'].values() if value is not None),
        'maximum_annular_oracle_difference': max(value for record in [*baseline, *measurements]
            for value in record['annular_oracle_max_errors'].values() if value is not None),
        'maximum_full_eigen_vs_parent_amplitude_difference': max(solver_errors)}
    write_json(root/'results.json', {'groups': groups, 'aggregate': aggregate,
        'paired_decisions': paired, 'nulls': nulls, 'measurements': measurements,
        'thresholds': thresholds, 'calibration_pools': pools, 'verification': verification,
        'caveats': [protocol['dependence'],
            'The fixed precision grid and inspected injections make this a development comparison.',
            'The exact-PSD efficiency treats the fitted PSD covariance as the true covariance.']})
    lines = ['# Rectangular-PSD precision regularization', '',
        'The raw rectangular ±20 PSD covariance and unsmoothed measured response stay fixed. '
        'Only the covariance inverse eigenvalue policy changes. The completed SNR-3/5/7 images, '
        'one-lambda/D trial exclusion, searches, and production annular SNR are reused. Each method '
        'has its own maximum-null threshold on the same raw-PSD calibration locations.', '',
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
    lines += ['', '## Mean paired amplitude throughput', '',
        'Each value is the positive-minus-baseline center amplitude divided by injected contrast. '
        'The Gaussian image is not a contrast estimator.', '',
        '| Method | '+' | '.join(labels)+' throughput |',
        '| --- | '+' | '.join('---:' for _ in levels)+' |']
    for row in aggregate:
        if row['method'] == 'gaussian':
            continue
        lines.append(f'| {LABELS[row["method"]]} | '+' | '.join(
            f'{level["mean_paired_throughput"]:.4f}' for level in row['levels'])+' |')
    lines += ['', '## Precision diagnostics if the fitted PSD were exact', '',
        'Efficiency is the regularized-filter SNR divided by the full-inverse SNR under the fitted '
        'covariance. Mode counts are the mean retained/modified eigenmodes out of 121.', '',
        '| Method | '+' | '.join(labels)+' efficiency (retained/modified) |',
        '| --- | '+' | '.join('---:' for _ in levels)+' |']
    for row in aggregate:
        if row['method'] not in PRECISION_METHODS:
            continue
        lines.append(f'| {LABELS[row["method"]]} | '+' | '.join(
            f'{level["mean_efficiency_if_psd_exact"]:.4f} '
            f'({level["mean_retained_modes"]:.1f}/{level["mean_modified_modes"]:.1f})'
            for level in row['levels'])+' |')
    lines += ['', '## Paired decision changes from the full inverse', '']
    for row in paired:
        label = labels[levels.index(row['brightness_multiplier'])]
        lines.append(f'- {row["policy"]}, cutoff {row["cutoff_over_mean"]:g} × mean variance, '
            f'{label}: regularized only {len(row["regularized_only"])}; '
            f'full inverse only {len(row["full_inverse_only"])}.')
    lines += ['', 'All thresholds were frozen before positive analysis. No P4 reduction was run.', '']
    (root/'results.md').write_text('\n'.join(lines))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    x_values = protocol['target_source_snrs']
    for column, (policy, variants) in enumerate((('Clipped inverse', CLIPPED),
                                                  ('Hard-truncated inverse', TRUNCATED))):
        plotted = [(None, FULL), *variants.items()]
        for cutoff, method in plotted:
            row = next(item for item in aggregate if item['method'] == method)
            label = 'Full inverse' if cutoff is None else f'{cutoff:g} × mean variance'
            axes[0, column].plot(x_values, [level['mean_search_snr'] for level in row['levels']],
                                 'o-', label=label)
            axes[1, column].plot(x_values, [level['detections']/36 for level in row['levels']],
                                 'o-', label=label)
        axes[0, column].set(title=policy, xlabel='Nominal identity source SNR',
                            ylabel='Mean five-pixel maximum SNR')
        axes[1, column].set(xlabel='Nominal identity source SNR', ylabel='Recovery fraction',
                            ylim=(-.03, 1.03))
        axes[0, column].legend(title='Eigenvalue cutoff')
    fig.suptitle('Raw rectangular ±20 PSD precision regularization')
    fig.savefig(root/'comparison.png', dpi=170)
    plt.close(fig)


def self_check() -> None:
    """Verify full, clipped, and hard-truncated precision algebra."""
    rng = np.random.default_rng(20260920)
    template = rng.normal(size=121)
    data = rng.normal(size=121)
    samples = rng.normal(size=(64, 121))
    model = inner.full.psd.fit_psd(samples, 'rectangular', .3)
    inner.full.radial.require(model is not None, 'synthetic PSD fit failed')
    eigenvalues, eigenvectors = eigh(model['covariance'], check_finite=False)
    full = inner.full.radial.filter_stamp(data, template, model, np.ones(121))
    reproduced = precision_result(data, template, model['mean'], model['covariance'],
                                  eigenvectors, 1/eigenvalues, float(full['sigma'])**2)
    inner.full.radial.require(np.isclose(reproduced['amplitude'], full['amplitude'],
                              rtol=1e-10, atol=1e-12) and
                              np.isclose(reproduced['psd_sigma'], full['sigma'], rtol=1e-10),
                              'full precision does not reproduce the production solve')
    for cutoff in CUTOFFS:
        threshold = cutoff*model['target_variance']
        for precision in (1/np.maximum(eigenvalues, threshold),
                          np.where(eigenvalues >= threshold, 1/eigenvalues, 0)):
            result = precision_result(data, template, model['mean'], model['covariance'],
                                      eigenvectors, precision, float(full['sigma'])**2)
            inner.full.radial.require(np.isfinite(list(result.values())).all() and
                                      result['efficiency_if_psd_exact'] <= 1+1e-10,
                                      'invalid regularized-precision result')
    print('PSD-precision checks passed', flush=True)


def run(root: Path) -> None:
    """Run calibration, then analyze all saved positives without new reductions."""
    manifest = inner.full.read(root/'manifest.json')
    inner.full.verify([*manifest['inputs'], manifest['protocol']])
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (root/'complete.json').exists():
            inner.full.verify(inner.full.read(root/'complete.json')['products'])
            print('PSD-precision comparison already complete.', flush=True)
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
                        default=Path('working/roc/p4_response_smoothing_20260920'))
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
