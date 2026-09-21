#!/usr/bin/env python3
"""Close the P4 precision study with frozen weights, null-tail checks, and planet SNR.

Reuse the completed precision comparison and every saved positive image. Fit
selected covariance weights only on the uninjected baseline, replay them on all
positives, audit the sensitivity of maximum-null thresholds on common support,
and measure the known planet with the same precision policies. No P4 reduction
is run.
"""
from __future__ import annotations

import argparse
import configparser
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

import compare_p4_step5_psd_precision as precision
import compare_p4_step5_response_smoothing as smoothing
from measure_p4_step5_planet import annular_oracle
from run_p4_step5_full_injections import fingerprint, write_json

inner = precision.inner
RADII = precision.RADII
RESAMPLES = 10000
RESAMPLE_SEED = 20260920

POLICIES = {'full': precision.FULL, 'clip1': precision.CLIPPED[1.0],
            'tsvd0p5': precision.TRUNCATED[0.5],
            'tsvd0p75': precision.TRUNCATED[0.75]}
REFIT = {key: 'refit_'+key for key in POLICIES}
FROZEN = {key: 'frozen_'+key for key in POLICIES}
REFERENCES = ('identity_response_lpf1p8', 'gaussian')
METHODS = (*REFIT.values(), *FROZEN.values(), *REFERENCES)
PARENT_METHOD = {**{REFIT[key]: parent for key, parent in POLICIES.items()},
                 **{FROZEN[key]: parent for key, parent in POLICIES.items()},
                 'identity_response_lpf1p8': 'identity_response_lpf1p8',
                 'gaussian': 'gaussian'}
LABELS = {REFIT['full']: 'Full inverse, refit on each image',
          FROZEN['full']: 'Full inverse, baseline-frozen',
          REFIT['clip1']: 'Clipped 1.0, refit on each image',
          FROZEN['clip1']: 'Clipped 1.0, baseline-frozen',
          REFIT['tsvd0p5']: 'Hard truncated 0.5, refit on each image',
          FROZEN['tsvd0p5']: 'Hard truncated 0.5, baseline-frozen',
          REFIT['tsvd0p75']: 'Hard truncated 0.75, refit on each image',
          FROZEN['tsvd0p75']: 'Hard truncated 0.75, baseline-frozen',
          'identity_response_lpf1p8': 'Identity, response LPF 1.8 px',
          'gaussian': 'Gaussian FWHM 3.6'}
PLANET_METHODS = (*POLICIES, 'identity_response_lpf1p8', 'gaussian')
PLANET_LABELS = {**{key: LABELS[FROZEN[key]].replace(', baseline-frozen', '')
                    for key in POLICIES},
                 'identity_response_lpf1p8': LABELS['identity_response_lpf1p8'],
                 'gaussian': LABELS['gaussian']}
PARENT_INDEX = {name: precision.METHODS.index(name) for name in precision.METHODS}
_CONTEXT = None


def deduplicate(records: list[dict]) -> list[dict]:
    """Return one identical fingerprint for each resolved input path."""
    result = {}
    for record in records:
        path = record['path']
        if path in result:
            inner.full.radial.require(result[path] == record,
                                      'conflicting input fingerprint: '+path)
        result[path] = record
    return list(result.values())


def parent_inputs(parent_root: Path, protocol: dict, site_names: set[str],
                  config: Path) -> list[dict]:
    """Collect every completed parent product read by the closure study."""
    manifest = inner.full.read(parent_root/'manifest.json')
    records = [*manifest['inputs'], manifest['protocol']]
    records += [fingerprint(parent_root/name) for name in ('manifest.json', 'protocol.json',
        'state.json', 'complete.json', 'results.json', 'baseline.json', 'thresholds.json',
        'effective_calibration_pools.json')]
    tasks = [task for task in [*protocol['baseline_tasks'], *protocol['positive_tasks']]
             if task['phase'] == 'positive' or task['trial_name'] in site_names]
    inner.full.radial.require(len(tasks) == 144, 'expected 36 site baselines and 108 positives')
    for task in tasks:
        directory = parent_root/task['phase']/task['name']
        complete = inner.full.read(directory/'complete.json')
        products = {Path(record['path']).name: record for record in complete['products']}
        inner.full.radial.require('amplitudes.fits' in products and
                                  'amplitudes_snr.fits' in products,
                                  'parent task lacks amplitude products')
        records += [fingerprint(directory/'complete.json'), products['amplitudes.fits'],
                    products['amplitudes_snr.fits']]
    study = Path(protocol['parent_study'])
    records += [fingerprint(config), fingerprint(study/'software/hciAnalyze'),
                fingerprint(study/'references/baseline/gaussian.fits'),
                fingerprint(Path(__file__)), fingerprint(Path(precision.__file__)),
                fingerprint(Path(smoothing.__file__)), fingerprint(Path(inner.__file__))]
    return deduplicate(records)


def read_planet_config(path: Path, study_protocol: dict, header) -> dict:
    """Read and validate the user's fixed planet-analysis configuration."""
    config = configparser.ConfigParser()
    config.read_string('[root]\n'+path.read_text())
    inner.full.radial.require(float(config['planet']['sep']) == float(header['PLANETSEP']) and
                              float(config['planet']['PA']) == float(header['PLANETPA']),
                              'planet coordinates changed')
    settings = {'lambda_d': float(config['root']['lambdaD']),
        'source_x': study_protocol['known_source_circle'][0],
        'source_y': study_protocol['known_source_circle'][1],
        'source_radius': float(config['planet']['R']),
        'aperture_radius': float(config['snr']['apertureR']),
        'min_radius': float(config['snr']['minRad']), 'max_radius': 60.0}
    inner.full.radial.require(settings['lambda_d'] == 3.6 and
                              settings['source_radius'] == 7 and
                              settings['aperture_radius'] == 3 and
                              settings['min_radius'] == 6 and
                              'maxRad' not in config['snr'],
                              'changed working/analyze.conf guidance')
    return settings


def prepare(args: argparse.Namespace) -> None:
    """Freeze the saved-image covariance-adaptation closure comparison."""
    root, parent_root = args.root.resolve(), args.parent.resolve()
    config = args.config.resolve()
    inner.full.radial.require(not root.exists(), 'closure root already exists')
    inner.full.radial.require((parent_root/'complete.json').exists() and
                              inner.full.read(parent_root/'state.json')['status'] == 'complete',
                              'precision parent is incomplete')
    inner.full.radial.require(bool(args.cpus) and len(args.cpus) == len(set(args.cpus)) and
                              set(args.cpus).issubset(os.sched_getaffinity(0)),
                              'requested CPUs are unavailable')
    prior = inner.full.read(parent_root/'protocol.json')
    study = Path(prior['parent_study'])
    study_protocol = inner.full.read(study/'protocol.json')
    sites = {site['name']: site for site in study_protocol['sites']}
    inner.full.radial.require(set(POLICIES.values()).issubset(prior['methods']) and
                              len(sites) == 36 and len(prior['positive_tasks']) == 108,
                              'changed completed precision design')
    science, header = fits.getdata(study/'payload/baseline.fits', header=True)
    inner.full.radial.require(science.squeeze().shape == (256, 256), 'changed baseline shape')
    planet_settings = read_planet_config(config, study_protocol, header)
    tasks_by_site = {name: [] for name in sites}
    for task in prior['positive_tasks']:
        tasks_by_site[task['trial_name']].append(task)
    baseline_task_by_site = {task['trial_name']: task['name'] for task in prior['baseline_tasks']
                             if task['trial_name'] in sites}
    inner.full.radial.require(all(len(tasks) == 3 for tasks in tasks_by_site.values()),
                              'expected three positives per site')
    inner.full.radial.require(set(baseline_task_by_site) == set(sites),
                              'missing parent baseline site task')
    records = parent_inputs(parent_root, prior, set(sites), config)
    inner.full.verify(records)
    protocol = {'schema': 1,
        'purpose': 'baseline-frozen PSD precision replay, common-null threshold audit, and known-planet SNR',
        'parent_comparison': str(parent_root), 'parent_study': str(study),
        'methods': list(METHODS), 'labels': LABELS,
        'policies': POLICIES, 'refit_methods': REFIT, 'frozen_methods': FROZEN,
        'references': list(REFERENCES), 'sites': list(sites), 'tasks_by_site': tasks_by_site,
        'baseline_task_by_site': baseline_task_by_site,
        'radii': list(RADII), 'brightnesses': prior['brightnesses'],
        'target_source_snrs': prior['target_source_snrs'], 'level_labels': prior['level_labels'],
        'frozen_covariance': 'fit raw rectangular ±20 PSD covariance and mean once on the uninjected baseline for each trial-specific exclusion and candidate pixel; apply the resulting unit-response weights and baseline mean unchanged to all three positive images',
        'refit_control': 'copy the completed per-image covariance-refit maps and reproduce their production annular SNR exactly',
        'thresholds': 'frozen and refit pairs share the completed parent baseline maximum-null threshold because their baseline maps are identical',
        'threshold_audit': {'common_support': 'all valid calibration trials inside each completed parent radial band',
            'draws': RESAMPLES, 'draw_size': 20, 'seed': RESAMPLE_SEED,
            'sampling': 'paired uniform subsets without replacement; bands with exactly 20 candidates have one fixed subset'},
        'planet': {**planet_settings,
            'methods': list(PLANET_METHODS), 'labels': PLANET_LABELS,
            'role': 'descriptive single-source endpoint; not used for policy selection'},
        'annular_trial_exclusion_radius': prior['annular_trial_exclusion_radius'],
        'workers': len(args.cpus), 'cpu_ids': args.cpus, 'threads_per_worker': 1,
        'new_reductions': 0,
        'dependence': 'reuses the same correlated residual field and repeatedly inspected injections; development closure, not independent validation'}
    root.mkdir(parents=True)
    shutil.copy2(config, root/'analyze.conf')
    write_json(root/'protocol.json', protocol)
    write_json(root/'manifest.json', {'schema': 1, 'inputs': records,
        'protocol': fingerprint(root/'protocol.json'),
        'config': fingerprint(root/'analyze.conf'), 'host': os.uname().nodename})
    write_json(root/'state.json', {'status': 'prepared', 'sites': len(sites),
        'saved_positive_images': len(prior['positive_tasks']), 'new_reductions': 0})
    print(f'prepared {len(sites)} frozen-weight site fits and '
          f'{len(prior["positive_tasks"])} saved-positive replays', flush=True)


def initialize(root: str, queue) -> None:
    """Pin one worker to one CPU and load the immutable baseline and response."""
    global _CONTEXT
    os.sched_setaffinity(0, {queue.get()})
    output = Path(root)
    protocol = inner.full.read(output/'protocol.json')
    parent_root = Path(protocol['parent_comparison'])
    study = Path(protocol['parent_study'])
    study_protocol = inner.full.read(study/'protocol.json')
    sites = {site['name']: site for site in study_protocol['sites']}
    baseline, header = fits.getdata(study/'payload/baseline.fits', header=True)
    _CONTEXT = (output, parent_root, study, study_protocol, sites,
                baseline.squeeze().astype(float), header,
                inner.full.load_templates(study/'payload/response'))


def archive_incomplete(root: Path, directory: Path) -> Path:
    """Preserve one unreceipted site directory before recomputing it."""
    archive = root/'interrupted'/directory.name
    archive.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while (archive/f'attempt_{attempt:04d}').exists():
        attempt += 1
    destination = archive/f'attempt_{attempt:04d}'
    shutil.move(directory, destination)
    return destination


def model_weights(science: np.ndarray, position: tuple[int, int], template: np.ndarray,
                  forbidden: np.ndarray) -> tuple[dict, np.ndarray, dict]:
    """Fit the baseline PSD model and return the four fixed unit-response weights."""
    rings = inner.full.radial.geometry(science, position, forbidden)
    matrices = {offset: inner.full.radial.extract(science, ring)
                for offset, ring in rings.items() if abs(offset) <= 20}
    samples = np.vstack(list(matrices.values()))
    model = inner.full.psd.fit_psd(samples, 'rectangular', .3)
    inner.full.radial.require(model is not None, 'baseline PSD model is unavailable')
    eigenvalues, eigenvectors = eigh(model['covariance'], check_finite=False)
    mean_variance = model['target_variance']
    precisions = {'full': 1/eigenvalues,
        'clip1': 1/np.maximum(eigenvalues, mean_variance),
        'tsvd0p5': np.where(eigenvalues >= .5*mean_variance, 1/eigenvalues, 0),
        'tsvd0p75': np.where(eigenvalues >= .75*mean_variance, 1/eigenvalues, 0)}
    coefficients = eigenvectors.T @ template
    weights, diagnostics = {}, {}
    for key, values in precisions.items():
        numerator = eigenvectors @ (values*coefficients)
        energy = float(template @ numerator)
        inner.full.radial.require(np.isfinite(energy) and energy > 0,
                                  'frozen precision lacks template support')
        weight = numerator/energy
        inner.full.radial.require(np.isclose(weight @ template, 1, rtol=1e-11, atol=1e-12),
                                  'frozen weight lacks unit response')
        weights[key] = weight
        cutoff = 0 if key == 'full' else 1 if key == 'clip1' else .5 if key == 'tsvd0p5' else .75
        modified = int(np.sum(eigenvalues < cutoff*mean_variance)) if cutoff else 0
        diagnostics[key] = {'retained_modes': 121 if key.startswith('clip') or key == 'full'
            else 121-modified, 'modified_modes': modified}
    production = inner.full.radial.filter_stamp(
        science[position[1]-5:position[1]+6, position[0]-5:position[0]+6].ravel(),
        template, model, np.ones(121))
    reproduced = float(weights['full'] @ (science[
        position[1]-5:position[1]+6, position[0]-5:position[0]+6].ravel()-model['mean']))
    inner.full.radial.require(np.isclose(reproduced, production['amplitude'],
                              rtol=1e-10, atol=1e-12),
                              'frozen full eigensolve differs from production solve')
    return weights, model['mean'], {'samples': len(samples), 'policies': diagnostics}


def select_methods(maps: np.ndarray, trial: dict, settings: dict) -> tuple[dict, dict]:
    """Evaluate the independent annular oracle for every output plane."""
    expected, profiles = {}, {}
    positions = [(trial['row']+dx, trial['column']+dy) for dx, dy in inner.SEARCH_OFFSETS]
    for index, name in enumerate(METHODS):
        expected[name], profiles[name] = annular_oracle(maps[index], settings)
        inner.full.radial.require(all(np.isfinite(maps[index, y, x]) and
                                      np.isfinite(expected[name][y, x])
                                      for x, y in positions),
                                  'injection site lacks complete support for '+name)
    return expected, profiles


def score_trial(trial: dict, snr: np.ndarray, amplitudes: np.ndarray) -> dict:
    """Measure the fixed center plus four axial neighboring pixels."""
    result = {}
    for index, name in enumerate(METHODS):
        values = [float(snr[index, trial['column']+dy, trial['row']+dx])
                  for dx, dy in inner.SEARCH_OFFSETS]
        raw = [float(amplitudes[index, trial['column']+dy, trial['row']+dx])
               for dx, dy in inner.SEARCH_OFFSETS]
        inner.full.radial.require(np.isfinite(values).all() and np.isfinite(raw).all(),
                                  'nonfinite fixed search for '+name)
        result[name] = {'valid': True, 'search_score': max(values),
            'snr_pixels': values, 'amplitude_pixels': raw, 'center_amplitude': raw[0]}
    return result


def analyze_site(site_name: str) -> str:
    """Fit baseline weights once and apply them to all positives at one site."""
    root, parent_root, study, study_protocol, sites, baseline, baseline_header, templates = _CONTEXT
    protocol = inner.full.read(root/'protocol.json')
    trial = sites[site_name]
    tasks = protocol['tasks_by_site'][site_name]
    directory = root/'sites'/site_name
    complete = directory/'complete.json'
    if complete.exists():
        inner.full.verify(inner.full.read(complete)['products'])
        return str(directory/'site_results.json')
    if directory.exists():
        destination = archive_incomplete(root, directory)
        print(f'{site_name}: archived incomplete task as {destination}', flush=True)
    directory.mkdir(parents=True)
    start = time.monotonic()
    positive = {}
    for task in tasks:
        image = fits.getdata(study/'reductions'/task['source_name']/'finim.fits')
        positive[task['name']] = image.squeeze().astype(float)
    yy, xx = np.indices(baseline.shape)
    radius = np.hypot(xx-127.5, yy-127.5).astype('f4')
    bins, positions, selected = inner.required_bins([trial], radius)
    source_mask, clean_centers = inner.planet_masks(baseline.shape,
                                                    study_protocol['known_source_circle'])
    forbidden = source_mask.copy()
    inner.full.mark_footprint(forbidden, trial['row'], trial['column'])
    parent_baseline_dir = parent_root/'baseline'/protocol['baseline_task_by_site'][site_name]
    parent_baseline = fits.getdata(parent_baseline_dir/'amplitudes.fits').reshape(
        (len(precision.METHODS), *baseline.shape))
    maps_by_task = {task['name']: np.full((len(METHODS), *baseline.shape), np.nan, dtype='f4')
                    for task in tasks}
    parent_maps, parent_snr = {}, {}
    for task in tasks:
        parent_dir = parent_root/'positive'/task['name']
        parent_maps[task['name']] = fits.getdata(parent_dir/'amplitudes.fits').reshape(
            (len(precision.METHODS), *baseline.shape))
        parent_snr[task['name']] = fits.getdata(parent_dir/'amplitudes_snr.fits').reshape(
            (len(precision.METHODS), *baseline.shape))
        for method in (*REFIT.values(), *REFERENCES):
            parent_method = PARENT_METHOD[method]
            maps_by_task[task['name']][METHODS.index(method)] = parent_maps[task['name']][
                PARENT_INDEX[parent_method]]
    control_errors = {key: [] for key in POLICIES}
    search_diagnostics = {}
    fitted_pixels = 0
    for (x, y), template in sorted(templates.items()):
        if not selected[y, x] or not clean_centers[y, x]:
            continue
        baseline_data = baseline[y-5:y+6, x-5:x+6].ravel()
        if baseline_data.size != 121 or not np.isfinite(baseline_data).all():
            continue
        weights, baseline_mean, diagnostic = model_weights(baseline, (x, y), template, forbidden)
        for key, parent_method in POLICIES.items():
            baseline_amplitude = float(weights[key] @ (baseline_data-baseline_mean))
            saved = float(parent_baseline[PARENT_INDEX[parent_method], y, x])
            inner.full.radial.require(np.isclose(baseline_amplitude, saved,
                                      rtol=2e-6, atol=1e-12),
                                      'frozen baseline differs from parent '+key)
            control_errors[key].append(abs(baseline_amplitude-saved))
            for task in tasks:
                data = positive[task['name']][y-5:y+6, x-5:x+6].ravel()
                maps_by_task[task['name']][METHODS.index(FROZEN[key]), y, x] = float(
                    weights[key] @ (data-baseline_mean))
        if (x, y) in positions:
            search_diagnostics[f'{x},{y}'] = diagnostic
        fitted_pixels += 1
    inner.full.radial.require(set(search_diagnostics) == {f'{x},{y}' for x, y in positions},
                              'missing frozen search-pixel diagnostics')
    baseline_control = {'site': site_name, 'trial': trial, 'required_bins': bins,
        'fitted_pixels': fitted_pixels, 'search_diagnostics': search_diagnostics,
        'maximum_parent_amplitude_difference': {key: max(values)
            for key, values in control_errors.items()}}
    write_json(directory/'baseline_control.json', baseline_control)
    measurements, products = [], [fingerprint(directory/'baseline_control.json')]
    for task in tasks:
        task_dir = directory/task['name']
        task_dir.mkdir()
        maps = maps_by_task[task['name']]
        header = baseline_header.copy()
        header['HCI FILTER LABELS'] = ','.join(METHODS)
        header['HCI RADIAL BINS'] = ','.join(map(str, bins))
        fits.writeto(task_dir/'amplitudes.fits', maps, header)
        delta_x, delta_y = trial['row']-127.5, trial['column']-127.5
        command = [str(study/'software/hciAnalyze'), '--file='+str(task_dir/'amplitudes.fits'),
            '--lambdaD=3.6', '--planet.sep=11.782,'+str(math.hypot(delta_x, delta_y)),
            '--planet.PA=262.051,'+str(math.degrees(-math.atan2(delta_x, delta_y)) % 360),
            '--planet.R=7.3,'+str(protocol['annular_trial_exclusion_radius']),
            '--snr.apertureR=60', '--snr.minRad=0', '--snr.maxRad=60',
            '--filter.psfResponse=', '--filter.lpfGaussFW=0', '--filter.hpfGaussFW=0',
            '--noise.model=identity', '--noise.only=false', '--noise.outputDiagnostics=false']
        environment = inner.full.binary_environment(study, reduction=False)
        environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        with (task_dir/'analysis.log').open('w') as log:
            subprocess.run(command, cwd=task_dir, env=environment, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        snr, snr_header = fits.getdata(task_dir/'amplitudes_snr.fits', header=True)
        snr = snr.reshape(maps.shape)
        inner.full.radial.require(snr_header['SNRMEAN'] == snr_header['SNRSMALL'] == 1,
                                  'changed production annular SNR contract')
        settings = {'source_x': study_protocol['known_source_circle'][0],
            'source_y': study_protocol['known_source_circle'][1],
            'source_radius': study_protocol['known_source_circle'][2],
            'lambda_d': 3.6, 'min_radius': 0, 'max_radius': 60,
            'source_exclusions': [study_protocol['known_source_circle'],
                [trial['row'], trial['column'], protocol['annular_trial_exclusion_radius']]]}
        expected, profiles = select_methods(maps, trial, settings)
        oracle_errors = {}
        for index, name in enumerate(METHODS):
            valid = np.isfinite(maps[index]) & np.isfinite(expected[name])
            oracle_errors[name] = float(np.max(np.abs(snr[index][valid]-expected[name][valid])))
            inner.full.radial.require(np.allclose(snr[index][valid], expected[name][valid],
                                      rtol=2e-6, atol=2e-6),
                                      'annular oracle mismatch for '+name)
        copied_snr_errors = {}
        for name in (*REFIT.values(), *REFERENCES):
            prior = parent_snr[task['name']][PARENT_INDEX[PARENT_METHOD[name]]]
            current = snr[METHODS.index(name)]
            copied_snr_errors[name] = precision.maximum_finite_difference(current, prior)
            inner.full.radial.require(np.allclose(current, prior, rtol=2e-6, atol=2e-6,
                                      equal_nan=True), 'copied parent SNR changed for '+name)
        measurement = {'task': task, 'trial': trial,
            'source': fingerprint(study/'reductions'/task['source_name']/'finim.fits'),
            'models': score_trial(trial, snr, maps),
            'annular_profile_min_pixels': {name: min(row['pixels'] for row in profile
                if int(row['radius']-.5) in bins) for name, profile in profiles.items()},
            'annular_oracle_max_errors': oracle_errors,
            'copied_parent_snr_max_errors': copied_snr_errors}
        write_json(task_dir/'measurements.json', measurement)
        write_json(task_dir/'command.json', command)
        task_products = [fingerprint(task_dir/name) for name in ('amplitudes.fits',
            'amplitudes_snr.fits', 'measurements.json', 'command.json')]
        write_json(task_dir/'complete.json', {'products': task_products})
        products += [*task_products, fingerprint(task_dir/'complete.json')]
        measurements.append(measurement)
    write_json(directory/'site_results.json', {'site': site_name, 'measurements': measurements,
        'elapsed_seconds': time.monotonic()-start})
    products.append(fingerprint(directory/'site_results.json'))
    write_json(complete, {'products': products})
    return str(directory/'site_results.json')


def run_sites(root: Path) -> list[dict]:
    """Run one resumable baseline-frozen fit for each injection site."""
    protocol = inner.full.read(root/'protocol.json')
    context = mp.get_context('spawn')
    queue = context.Queue()
    for cpu in protocol['cpu_ids']:
        queue.put(cpu)
    paths = []
    with ProcessPoolExecutor(max_workers=protocol['workers'], mp_context=context,
                             initializer=initialize, initargs=(str(root), queue)) as executor:
        futures = {executor.submit(analyze_site, name): name for name in protocol['sites']}
        for future in as_completed(futures):
            path = future.result()
            paths.append(path)
            write_json(root/'state.json', {'status': 'positive_replay', 'pid': os.getpid(),
                'finished_sites': len(paths), 'total_sites': len(futures),
                'last': futures[future], 'new_reductions': 0})
            print(f'site: {len(paths)}/{len(futures)} {futures[future]}', flush=True)
    measurements = []
    for path in sorted(paths):
        measurements += inner.full.read(Path(path))['measurements']
    return measurements


def parent_baseline_models(parent_results: dict, thresholds: dict) -> list[dict]:
    """Duplicate parent null measurements for refit and frozen baseline-identical arms."""
    result = []
    for parent_record in parent_results['nulls']:
        trial = parent_record['trial']
        models = {}
        for method in METHODS:
            value = dict(parent_record['models'][PARENT_METHOD[method]])
            value['detected'] = value['valid'] and value['search_score'] > thresholds[
                str(trial['nominal_radius'])][method]
            models[method] = value
        result.append({'trial': trial, 'models': models})
    return result


def output_thresholds(parent_thresholds: dict) -> dict:
    """Copy parent thresholds so each frozen/refit pair shares one baseline threshold."""
    return {radius: {method: parent_thresholds[radius][PARENT_METHOD[method]]
                     for method in METHODS} for radius in parent_thresholds}


def distribution(values: np.ndarray) -> dict:
    """Summarize one finite scalar distribution with fixed quantiles."""
    inner.full.radial.require(values.size and np.isfinite(values).all(),
                              'cannot summarize an empty or nonfinite distribution')
    return {'minimum': float(np.min(values)), 'p05': float(np.quantile(values, .05)),
        'median': float(np.median(values)), 'p95': float(np.quantile(values, .95)),
        'maximum': float(np.max(values))}


def threshold_audit(root: Path, measurements: list[dict], nulls: list[dict]) -> dict:
    """Audit common-support null tails and paired 20-location maximum thresholds."""
    protocol = inner.full.read(root/'protocol.json')
    parent_root = Path(protocol['parent_comparison'])
    parent_results = inner.full.read(parent_root/'results.json')
    baseline_rows = inner.full.read(parent_root/'baseline.json')['rows']
    study = inner.full.read(Path(protocol['parent_study'])/'protocol.json')
    calibration_names = {trial['name'] for trial in study['calibration_trials']}
    rows = {record['trial']['name']: record for record in baseline_rows}
    rng = np.random.default_rng(RESAMPLE_SEED)
    thresholds = {method: {} for method in POLICIES}
    tails = []
    for radius in RADII:
        pool = parent_results['calibration_pools'][str(radius)][precision.FULL]
        lower, upper = pool['radius_range']
        eligible = [record for record in baseline_rows
                    if record['trial']['name'] in calibration_names and
                    lower <= record['trial']['nominal_radius'] <= upper and
                    all(record['models'][parent_method]['valid']
                        for parent_method in POLICIES.values())]
        inner.full.radial.require(len(eligible) >= 20,
                                  'fewer than 20 common-support calibration trials')
        if len(eligible) == 20:
            indices = np.tile(np.arange(20), (RESAMPLES, 1))
        else:
            random_keys = rng.random((RESAMPLES, len(eligible)))
            indices = np.argpartition(random_keys, 19, axis=1)[:, :20]
        selected_names = pool['trials']
        for key, parent_method in POLICIES.items():
            values = np.array([record['models'][parent_method]['search_score']
                               for record in eligible])
            selected = np.array([rows[name]['models'][parent_method]['search_score']
                                 for name in selected_names])
            ordered = np.sort(selected)
            thresholds[key][radius] = np.max(values[indices], axis=1)
            tails.append({'radius': radius, 'policy': key,
                'eligible_common_trials': len(eligible),
                'eligible_distribution': distribution(values),
                'selected_distribution': distribution(selected),
                'selected_second_largest': float(ordered[-2]),
                'selected_maximum_gap': float(ordered[-1]-ordered[-2]),
                'completed_threshold': float(ordered[-1]),
                'resampled_threshold': distribution(thresholds[key][radius])})
    positive_by_level = {level: [record for record in measurements
        if record['task']['job']['brightness_multiplier'] == level]
        for level in protocol['brightnesses']}
    sensitivity = []
    for key in POLICIES:
        for mode, method in (('refit', REFIT[key]), ('frozen', FROZEN[key])):
            null_counts = np.zeros(RESAMPLES, dtype=int)
            for record in nulls:
                radius = record['trial']['nominal_radius']
                null_counts += record['models'][method]['search_score'] > thresholds[key][radius]
            levels = []
            for level, records_at_level in positive_by_level.items():
                detections = np.zeros(RESAMPLES, dtype=int)
                for record in records_at_level:
                    radius = record['task']['job']['nominal_radius']
                    detections += record['models'][method]['search_score'] > thresholds[key][radius]
                levels.append({'brightness_multiplier': level,
                               'detections': distribution(detections.astype(float))})
            sensitivity.append({'policy': key, 'mode': mode,
                                'heldout_null_exceedances': distribution(null_counts.astype(float)),
                                'levels': levels})
    return {'draws': RESAMPLES, 'draw_size': 20, 'seed': RESAMPLE_SEED,
        'paired_common_support': True, 'tails': tails, 'sensitivity': sensitivity,
        'limitation': 'Radius 6 and radii 16–24 have exactly 20 common-support candidates, so their subset thresholds cannot vary without widening the frozen radial band.'}


def planet_analysis(root: Path) -> dict:
    """Measure the known planet with the selected precision and reference filters."""
    directory = root/'planet'
    complete = directory/'complete.json'
    if complete.exists():
        inner.full.verify(inner.full.read(complete)['products'])
        return inner.full.read(directory/'results.json')
    protocol = inner.full.read(root/'protocol.json')
    study = Path(protocol['parent_study'])
    study_protocol = inner.full.read(study/'protocol.json')
    science, header = fits.getdata(study/'payload/baseline.fits', header=True)
    science = science.squeeze().astype(float)
    settings = protocol['planet']
    yy, xx = np.indices(science.shape)
    radius = np.hypot(xx-127.5, yy-127.5).astype('f4')
    distance = np.hypot(xx-settings['source_x'], yy-settings['source_y']).astype('f4')
    aperture = distance <= settings['aperture_radius']+.5
    bins = set()
    for y, x in np.argwhere(aperture):
        lower = math.floor(float(radius[y, x])-.5)
        bins.update((lower, lower+1))
    selected = np.zeros(science.shape, bool)
    for lower in bins:
        selected |= (radius > lower) & (radius <= lower+1)
    forbidden = distance <= study_protocol['known_source_circle'][2]+.5
    for y, x in np.argwhere(aperture):
        inner.full.mark_footprint(forbidden, int(x), int(y))
    directory.mkdir()
    fits.writeto(directory/'training_exclusion.fits', forbidden.astype('u1'), header)
    templates = inner.full.load_templates(study/'payload/response')
    maps = np.full((len(PLANET_METHODS), *science.shape), np.nan, dtype='f4')
    diagnostics = {}
    fitted = 0
    for (x, y), template in sorted(templates.items()):
        if not selected[y, x]:
            continue
        data = science[y-5:y+6, x-5:x+6].ravel()
        if data.size != 121 or not np.isfinite(data).all():
            continue
        weights, fitted_mean, detail = model_weights(science, (x, y), template, forbidden)
        for key in POLICIES:
            maps[PLANET_METHODS.index(key), y, x] = float(weights[key] @ (data-fitted_mean))
        filtered = smoothing.smooth_template(template, 1.8)
        maps[PLANET_METHODS.index('identity_response_lpf1p8'), y, x] = float(
            filtered @ data/(filtered @ filtered))
        if aperture[y, x]:
            diagnostics[f'{x},{y}'] = detail
        fitted += 1
    gaussian = fits.getdata(study/'references/baseline/gaussian.fits').squeeze().astype(float)
    maps[PLANET_METHODS.index('gaussian')] = np.where(selected, gaussian, np.nan)
    inner.full.radial.require(all(np.isfinite(maps[index][aperture]).all()
                                  for index in range(len(PLANET_METHODS))),
                              'planet aperture lacks complete method support')
    header['HCI FILTER LABELS'] = ','.join(PLANET_METHODS)
    header['HCI RADIAL BINS'] = ','.join(map(str, sorted(bins)))
    fits.writeto(directory/'amplitudes.fits', maps, header)
    command = [str(study/'software/hciAnalyze'), '--config', str(root/'analyze.conf'),
        '--file='+str(directory/'amplitudes.fits'), '--filter.psfResponse=',
        '--filter.lpfGaussFW=0', '--filter.hpfGaussFW=0', '--noise.model=identity',
        '--noise.only=false', '--noise.outputDiagnostics=false']
    environment = inner.full.binary_environment(study, reduction=False)
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (directory/'analysis.log').open('w') as log:
        subprocess.run(command, cwd=directory, env=environment, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    snr, snr_header = fits.getdata(directory/'amplitudes_snr.fits', header=True)
    snr = snr.reshape(maps.shape)
    inner.full.radial.require(snr_header['SNRAPER'] == 3 and snr_header['SNRMINR'] == 6 and
                              snr_header['SNRMAXR'] == 60 and
                              snr_header['SNRMEAN'] == snr_header['SNRSMALL'] == 1,
                              'planet production SNR settings changed')
    summaries, oracle_errors = {}, {}
    nearest_y, nearest_x = min(np.argwhere(aperture),
        key=lambda point: float(distance[point[0], point[1]]))
    oracle_settings = {'source_x': settings['source_x'], 'source_y': settings['source_y'],
        'source_radius': settings['source_radius'], 'lambda_d': settings['lambda_d'],
        'min_radius': settings['min_radius'], 'max_radius': settings['max_radius']}
    for index, name in enumerate(PLANET_METHODS):
        expected, profile = annular_oracle(maps[index], oracle_settings)
        valid = np.isfinite(maps[index]) & np.isfinite(expected)
        oracle_errors[name] = float(np.max(np.abs(snr[index][valid]-expected[valid])))
        inner.full.radial.require(np.allclose(snr[index][valid], expected[valid],
                                  rtol=2e-6, atol=2e-6),
                                  'planet annular oracle mismatch for '+name)
        values = np.where(aperture, snr[index], np.nan)
        peak = int(np.nanargmax(values))
        py, px = np.unravel_index(peak, values.shape)
        summaries[name] = {'label': PLANET_LABELS[name],
            'annular_snr_nearest_pixel': float(snr[index, nearest_y, nearest_x]),
            'nearest_pixel_x_y': [int(nearest_x), int(nearest_y)],
            'annular_snr_aperture_maximum': float(snr[index, py, px]),
            'peak_x_y': [int(px), int(py)],
            'valid_aperture_pixels': int(np.isfinite(values).sum()),
            'minimum_annular_pixels': min(row['pixels'] for row in profile
                if int(row['radius']-.5) in bins)}
    result = {'settings': settings, 'methods': summaries, 'required_bins': sorted(bins),
        'fitted_pixels': fitted, 'aperture_pixels': int(aperture.sum()),
        'search_diagnostics': diagnostics, 'annular_oracle_max_errors': oracle_errors,
        'interpretation': 'descriptive single known planet; not a method-selection sample'}
    write_json(directory/'results.json', result)
    write_json(directory/'command.json', command)
    products = [fingerprint(directory/name) for name in ('training_exclusion.fits',
        'amplitudes.fits', 'amplitudes_snr.fits', 'results.json', 'command.json')]
    write_json(complete, {'products': products})
    return result


def summarize(root: Path, measurements: list[dict], planet: dict) -> None:
    """Summarize frozen/refit recovery, threshold sensitivity, and planet SNR."""
    protocol = inner.full.read(root/'protocol.json')
    parent_root = Path(protocol['parent_comparison'])
    parent_results = inner.full.read(parent_root/'results.json')
    parent_thresholds = inner.full.read(parent_root/'thresholds.json')
    thresholds = output_thresholds(parent_thresholds)
    nulls = parent_baseline_models(parent_results, thresholds)
    baseline_by_site = {record['trial']['name']: record for record in nulls}
    for record in measurements:
        job = record['task']['job']
        for method, value in record['models'].items():
            value['detected'] = value['search_score'] > thresholds[str(job['nominal_radius'])][method]
            if method != 'gaussian':
                baseline = baseline_by_site[job['site']]['models'][method]['center_amplitude']
                value['paired_amplitude_increment'] = value['center_amplitude']-baseline
                value['paired_throughput'] = value['paired_amplitude_increment']/job['contrast']
    levels = protocol['brightnesses']
    groups, aggregate = [], []
    for radius in RADII:
        for method in METHODS:
            null_rows = [record for record in nulls if record['trial']['nominal_radius'] == radius]
            row = {'radius': radius, 'method': method,
                'null_exceedances': sum(record['models'][method]['detected'] for record in null_rows),
                'levels': []}
            for level in levels:
                selected = [record['models'][method] for record in measurements
                    if record['task']['job']['nominal_radius'] == radius and
                    record['task']['job']['brightness_multiplier'] == level]
                throughputs = [value['paired_throughput'] for value in selected
                               if 'paired_throughput' in value]
                row['levels'].append({'brightness_multiplier': level,
                    'trials': len(selected), 'detections': sum(value['detected'] for value in selected),
                    'mean_search_snr': float(np.mean([value['search_score'] for value in selected])),
                    'mean_center_snr': float(np.mean([value['snr_pixels'][0] for value in selected])),
                    'mean_paired_throughput': float(np.mean(throughputs)) if throughputs else None})
            groups.append(row)
    for method in METHODS:
        method_groups = [group for group in groups if group['method'] == method]
        aggregate_levels = []
        for level in levels:
            selected = [record['models'][method] for record in measurements
                        if record['task']['job']['brightness_multiplier'] == level]
            throughputs = [value['paired_throughput'] for value in selected
                           if 'paired_throughput' in value]
            aggregate_levels.append({'brightness_multiplier': level, 'trials': len(selected),
                'detections': sum(value['detected'] for value in selected),
                'mean_search_snr': float(np.mean([value['search_score'] for value in selected])),
                'mean_center_snr': float(np.mean([value['snr_pixels'][0] for value in selected])),
                'mean_paired_throughput': float(np.mean(throughputs)) if throughputs else None})
        aggregate.append({'method': method,
            'null_exceedances': sum(group['null_exceedances'] for group in method_groups),
            'levels': aggregate_levels})
    adaptation = []
    for key in POLICIES:
        for radius in RADII:
            for level in levels:
                selected = [record for record in measurements
                    if record['task']['job']['nominal_radius'] == radius and
                    record['task']['job']['brightness_multiplier'] == level]
                refit_values = [record['models'][REFIT[key]] for record in selected]
                frozen_values = [record['models'][FROZEN[key]] for record in selected]
                adaptation.append({'policy': key, 'radius': radius,
                    'brightness_multiplier': level,
                    'mean_search_snr_change_frozen_minus_refit': float(np.mean([
                        frozen['search_score']-refit['search_score']
                        for frozen, refit in zip(frozen_values, refit_values)])),
                    'mean_center_snr_change_frozen_minus_refit': float(np.mean([
                        frozen['snr_pixels'][0]-refit['snr_pixels'][0]
                        for frozen, refit in zip(frozen_values, refit_values)])),
                    'mean_throughput_change_frozen_minus_refit': float(np.mean([
                        frozen['paired_throughput']-refit['paired_throughput']
                        for frozen, refit in zip(frozen_values, refit_values)])),
                    'frozen_only': [record['task']['job']['name'] for record in selected
                        if record['models'][FROZEN[key]]['detected'] and
                        not record['models'][REFIT[key]]['detected']],
                    'refit_only': [record['task']['job']['name'] for record in selected
                        if record['models'][REFIT[key]]['detected'] and
                        not record['models'][FROZEN[key]]['detected']]})
    audit = threshold_audit(root, measurements, nulls)
    verification = {'new_reductions': 0, 'copied_parent_snr_reproduced': True,
        'maximum_copied_parent_snr_difference': max(value
            for record in measurements for value in record['copied_parent_snr_max_errors'].values()
            if value is not None),
        'maximum_annular_oracle_difference': max(value
            for record in measurements for value in record['annular_oracle_max_errors'].values()),
        'maximum_frozen_baseline_parent_amplitude_difference': max(value
            for site in protocol['sites'] for value in inner.full.read(
                root/'sites'/site/'baseline_control.json')[
                    'maximum_parent_amplitude_difference'].values())}
    write_json(root/'results.json', {'aggregate': aggregate, 'groups': groups,
        'adaptation': adaptation, 'threshold_audit': audit, 'planet': planet,
        'nulls': nulls, 'measurements': measurements, 'thresholds': thresholds,
        'verification': verification,
        'caveats': [protocol['dependence'], audit['limitation'],
            'The known planet is a descriptive single-source measurement.']})
    lines = ['# P4 precision closure: frozen covariance, null tails, and known planet', '',
        'All positive images and thresholds are reused. Frozen methods fit the covariance, mean, and '
        'unit-response weight only on the uninjected baseline for the matching trial exclusion.', '',
        '## Aggregate refit versus frozen recovery', '',
        '| Method | '+' | '.join(protocol['level_labels'])+' detections | Nulls | '
        + ' | '.join(protocol['level_labels'])+' mean max SNR |',
        '| --- | '+' | '.join('---:' for _ in levels)+' | ---: | '
        + ' | '.join('---:' for _ in levels)+' |']
    for row in aggregate:
        lines.append(f'| {LABELS[row["method"]]} | '
            + ' | '.join(str(level['detections']) for level in row['levels'])
            + f' | {row["null_exceedances"]}/36 | '
            + ' | '.join(f'{level["mean_search_snr"]:.4f}' for level in row['levels'])+' |')
    lines += ['', '## Mean frozen-minus-refit SNR by radius', '',
        '| Policy | Radius | '+' | '.join(protocol['level_labels'])+' maximum-SNR change |',
        '| --- | ---: | '+' | '.join('---:' for _ in levels)+' |']
    for key in POLICIES:
        for radius in RADII:
            rows = [row for row in adaptation if row['policy'] == key and row['radius'] == radius]
            lines.append(f'| {key} | {radius} | '
                + ' | '.join(f'{row["mean_search_snr_change_frozen_minus_refit"]:+.4f}'
                             for row in rows)+' |')
    lines += ['', '## Known planet', '',
        '| Method | Nearest-pixel SNR | Aperture maximum SNR | Peak x,y |',
        '| --- | ---: | ---: | --- |']
    for name in PLANET_METHODS:
        row = planet['methods'][name]
        lines.append(f'| {row["label"]} | {row["annular_snr_nearest_pixel"]:.4f} | '
            f'{row["annular_snr_aperture_maximum"]:.4f} | {row["peak_x_y"]} |')
    lines += ['', 'Threshold-tail and resampling diagnostics are recorded in `results.json`. '
        'No P4 reduction was run.', '']
    (root/'results.md').write_text('\n'.join(lines))
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), layout='constrained')
    for axis, key in zip(axes.flat, POLICIES):
        for mode, method, marker in (('Refit', REFIT[key], 'o-'),
                                     ('Baseline-frozen', FROZEN[key], 's--')):
            row = next(item for item in aggregate if item['method'] == method)
            axis.plot(protocol['target_source_snrs'],
                      [level['mean_search_snr'] for level in row['levels']], marker, label=mode)
        gaussian = next(item for item in aggregate if item['method'] == 'gaussian')
        axis.plot(protocol['target_source_snrs'],
                  [level['mean_search_snr'] for level in gaussian['levels']], 'k:', label='Gaussian')
        axis.set(title=PLANET_LABELS[key], xlabel='Nominal identity source SNR',
                 ylabel='Mean five-pixel maximum SNR')
        axis.legend()
    fig.suptitle('Per-image covariance refit versus baseline-frozen precision')
    fig.savefig(root/'comparison.png', dpi=170)
    plt.close(fig)


def self_check() -> None:
    """Verify fixed-weight algebra and deterministic subset-threshold sampling."""
    rng = np.random.default_rng(RESAMPLE_SEED)
    matrix = rng.normal(size=(121, 121))
    covariance = matrix @ matrix.T + .3*np.trace(matrix @ matrix.T)/121*np.eye(121)
    eigenvalues, eigenvectors = eigh(covariance, check_finite=False)
    template = rng.normal(size=121)
    mean = rng.normal(size=121)
    baseline = rng.normal(size=121)
    coefficients = eigenvectors.T @ template
    numerator = eigenvectors @ (coefficients/eigenvalues)
    weight = numerator/(template @ numerator)
    direct = np.linalg.solve(covariance, template)
    direct /= template @ direct
    inner.full.radial.require(np.allclose(weight, direct, rtol=1e-10, atol=1e-12) and
                              np.isclose(weight @ template, 1, rtol=1e-12),
                              'frozen full weight algebra failed')
    first = float(weight @ (baseline-mean))
    second = float(weight @ (baseline+template-mean))
    inner.full.radial.require(np.isclose(second-first, 1, rtol=1e-12, atol=1e-12),
                              'frozen weight lost unit response')
    random_keys = rng.random((100, 24))
    indices = np.argpartition(random_keys, 19, axis=1)[:, :20]
    inner.full.radial.require(indices.shape == (100, 20) and
                              all(len(set(row)) == 20 for row in indices),
                              'threshold subset sampler is not without replacement')
    print('precision-closure checks passed', flush=True)


def run(root: Path) -> None:
    """Run saved-positive replay, planet measurement, and final closure summary."""
    manifest = inner.full.read(root/'manifest.json')
    inner.full.verify([*manifest['inputs'], manifest['protocol'], manifest['config']])
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (root/'complete.json').exists():
            inner.full.verify(inner.full.read(root/'complete.json')['products'])
            print('Precision closure already complete.', flush=True)
            return
        try:
            write_json(root/'state.json', {'status': 'positive_replay', 'pid': os.getpid(),
                'new_reductions': 0})
            measurements = run_sites(root)
            write_json(root/'state.json', {'status': 'planet', 'pid': os.getpid(),
                'new_reductions': 0})
            planet = planet_analysis(root)
            write_json(root/'state.json', {'status': 'summarize', 'pid': os.getpid(),
                'new_reductions': 0})
            summarize(root, measurements, planet)
            inner.full.verify([*manifest['inputs'], manifest['protocol'], manifest['config']])
            products = [fingerprint(root/name) for name in ('results.json', 'results.md',
                'comparison.png')]
            products += [fingerprint(root/'planet'/name) for name in ('results.json',
                'amplitudes.fits', 'amplitudes_snr.fits', 'complete.json')]
            write_json(root/'complete.json', {'sites': len(inner.full.read(root/'protocol.json')['sites']),
                'positive_analyses': len(measurements), 'new_reductions': 0,
                'frozen_inputs_unchanged': True, 'products': products})
            write_json(root/'state.json', {'status': 'complete', 'new_reductions': 0,
                'results': str(root/'results.json')})
        except Exception as error:
            write_json(root/'state.json', {'status': 'failed', 'pid': os.getpid(),
                'new_reductions': 0, 'error': str(error)})
            raise


def main() -> None:
    """Expose algebra validation, immutable preparation, and resumable execution."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('check', 'prepare', 'run'))
    parser.add_argument('--root', type=Path)
    parser.add_argument('--parent', type=Path,
                        default=Path('working/roc/p4_psd_precision_20260920'))
    parser.add_argument('--config', type=Path, default=Path('working/analyze.conf'))
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
