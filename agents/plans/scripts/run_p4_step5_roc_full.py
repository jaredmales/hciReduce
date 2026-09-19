#!/usr/bin/env python3
"""Prepare and run the full, unattended Step-5 PSD injection study on ROC.

Design selects new native sites by geometry only. Each site has its own held-out
footprint and calibration; the same field therefore supplies correlated trials.
Prepare freezes native binaries and dependencies. Run produces one baseline,
freezes all thresholds/contrasts, then reduces and analyzes all 90 injections.
"""
from __future__ import annotations

import argparse
import fcntl
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.ndimage import binary_dilation
from astropy.io import fits
import astropy

import compare_p4_step5_radial_pooling as radial
import compare_p4_step5_shrinkage as shrinkage
import compare_p4_step5_welch_psd as psd
import run_p4_step5_gaussian as gaussian
from run_p4_step5_full_injections import fingerprint, replacement, write_json

COVARIANCE_MODELS = ('psd_hann_b0_m0.1', 'psd_rectangular_b5_m0.3', 'pca_b5_f1', 'isotropic_b5')
REFERENCES = ('identity', 'gaussian_snr', 'gaussian_raw', 'identity_snr')
MODELS = (*COVARIANCE_MODELS, *REFERENCES)
LEVELS = (.5, .75, 1.)
SCRIPTS = ('run_p4_step5_roc_full.py', 'compare_p4_step5_radial_pooling.py',
           'compare_p4_step5_shrinkage.py', 'compare_p4_step5_welch_psd.py',
           'compare_p4_step5_variance_floor.py', 'diagnose_p4_step5_projected_noise.py',
           'analyze_p4_step5_development.py', 'run_p4_step5_full_injections.py', 'run_p4_step5_gaussian.py')


def read(path: Path):
    """Read one frozen protocol or atomic progress record."""
    return radial.read_json(path)


def verify(records: list) -> None:
    """Stop if a declared source, dependency, result, or protocol has changed."""
    for record in records:
        radial.require(fingerprint(Path(record['path'])) == record, 'changed frozen file: '+record['path'])


def mark_footprint(mask: np.ndarray, x: int, y: int, half_width: int = 7) -> None:
    """Mark the exact native square-kernel union for the five-pixel search aperture."""
    for ox, oy in radial.OFFSETS:
        mask[y+oy-half_width:y+oy+half_width+1, x+ox-half_width:x+ox+half_width+1] = True


def load_templates(directory: Path) -> dict:
    """Read unchanged, finite unit-source stamps in native candidate coordinates."""
    coordinates = fits.getdata(directory/'p4PSF_coordinates.fits').T
    models = fits.getdata(directory/'p4PSF_model_0000.fits').astype(float)
    valid = fits.getdata(directory/'p4PSF_validity_0000.fits').ravel()
    return {(int(c[0]), int(c[1])): m.ravel() for c, m, v in zip(coordinates, models, valid)
            if v == 1 and np.isfinite(m).all()}


def holdout_mask(shape: tuple, protocol: dict, site: dict) -> np.ndarray:
    """Exclude the known source, all calibration footprints, and the current tested site."""
    yy, xx = np.indices(shape)
    x, y, radius = protocol['known_source_circle']
    mask = np.hypot(xx-x, yy-y) <= radius
    for trial in [*protocol['calibration_trials'], site]:
        mark_footprint(mask, trial['row'], trial['column'])
    return mask


def support_counts(science: np.ndarray, position: tuple, mask: np.ndarray) -> list:
    """Audit same-radius and five-pixel-band sample counts without evaluating scores."""
    x, y = position
    result = []
    for ox, oy in radial.OFFSETS:
        rings = radial.geometry(science, (x+ox, y+oy), mask)
        result.append({'position': [x+ox, y+oy], 'same_radius': len(rings[0]['halves']),
                       'band5': sum(len(r['halves']) for o, r in rings.items() if abs(o) <= 5)})
    return result


def design(args: argparse.Namespace) -> None:
    """Freeze 30 geometry-selected sites and stage portable source data and analysis scripts."""
    root = args.root.resolve()
    old = read(args.previous_study/'protocol.json')
    old_nulls = read(args.previous_study/'null_results.json')
    sampling = read(args.previous_sampling/'manifest.json')
    batch = read(args.source_batch/'manifest.json')
    verify([*batch['input_records'], *batch['frozen_records'], old['science'], *old['source_response_inputs']])
    science = fits.getdata(old['science']['path']).squeeze().astype(float)
    radial.require(science.shape == (256, 256), 'unexpected baseline shape')
    field = Path(old['source_response_inputs'][0]['path']).parent
    templates = load_templates(field)
    calibration = [{k: t[k] for k in ('name', 'row', 'column', 'nominal_radius', 'angular_block')}
                   for t in old_nulls['trials'] if t['role'] == 'calibration' and t['common_eligible']]
    radial.require(len(calibration) == 28, 'changed original calibration set')
    protocol = {'schema': 1, 'purpose': '30 new injection sites on the existing field; sitewise held-out calibration',
        'source_cpp_commit': args.source_commit, 'pilot': 'separate ROC numerical/throughput pilot waived by user',
        'models': list(MODELS), 'brightness_multipliers': LEVELS, 'nominal_radii': [26, 32, 38, 44, 50],
        'sites_per_radius': 6, 'full_positive_reductions': 90, 'baseline_reductions': 1,
        'calibration_trials': calibration, 'known_source_circle': sampling['exclusions'][0],
        'holdout': 'per site: known source plus exact 15x15 Gaussian/five-pixel search footprints of all 28 calibration trials and this site; other new sites and historical development pixels may train',
        'threshold_rule': 'per site and model, ceil((28+1)*0.95) calibration order statistic; maximum; strict exceedance',
        'reference_contrast': 'site-specific isotropic_b5 calibration threshold times its candidate sigma; target stamp values are not used',
        'selection': 'radii within 0.65 pixel, azimuth 90..270 degrees, no overlap with prior Gaussian/development measurement footprints; six evenly spaced angular ranks with minimum four-pixel same-ring center spacing; finite response and >=8 training patches',
        'dependence': 'sites, their calibration sets, and training pixels overlap; these are not independent trials or a new residual field',
        'covariance': 'raw pixels, 11x11 untapered candidates, radial/angular step 5, minimum 8; same-radius Hann/mix .1 and band5 rectangular/mix .3; matching band5 PCA floor1 and mean-variance isotropic controls',
        'gaussian': 'production FWHM3.6, plus ordinary application SNR; legacy annular normalization includes trial neighborhoods',
        'source': 'same 621 frames, P4-M32D64 mode fraction .15, full image, mean combine, existing cropped 12x12 injection PSF without renormalization',
        'reduction_cpu_ids': list(range(24)), 'reduction_openmp_threads': 24, 'blas_threads': 1,
        'analysis_cpu_ids': [12, 13], 'analysis_openmp_threads': 2}
    yy, xx = np.indices(science.shape)
    old_measurements = np.zeros(science.shape, bool)
    for trial in old['trials']:
        mark_footprint(old_measurements, trial['row'], trial['column'])
    for x, y in ((120, 137), (109, 143), (98, 98)):
        mark_footprint(old_measurements, x, y, 5)
    x, y, radius = protocol['known_source_circle']
    old_measurements |= np.hypot(xx-x, yy-y) <= radius
    aperture = np.zeros((17, 17), bool)
    mark_footprint(aperture, 8, 8)
    excluded_centers = binary_dilation(old_measurements, structure=aperture)
    separation = np.hypot(xx-127.5, yy-127.5)
    azimuth = np.degrees(np.arctan2(yy-127.5, xx-127.5)) % 360
    sites, audits, pool_counts = [], [], {}
    for nominal in protocol['nominal_radii']:
        permitted = (np.abs(separation-nominal) < .65) & ~excluded_centers & (azimuth >= 90) & (azimuth <= 270)
        coordinates = sorted([(int(x), int(y)) for y, x in np.argwhere(permitted)], key=lambda p: azimuth[p[1], p[0]])
        viable = []
        for x, y in coordinates:
            if not all((x+ox, y+oy) in templates for ox, oy in radial.OFFSETS):
                continue
            site = {'row': x, 'column': y}
            counts = support_counts(science, (x, y), holdout_mask(science.shape, protocol, site))
            if min(c['same_radius'] for c in counts) >= 8:
                viable.append((x, y))
        pool_counts[str(nominal)] = len(viable)
        radial.require(len(viable) >= 6, 'insufficient score-free candidate pool')
        chosen = []
        for target in np.linspace(0, len(viable)-1, 6):
            ranked = sorted(range(len(viable)), key=lambda i: abs(i-target))
            choice = next((viable[i] for i in ranked if all(math.dist(viable[i], p) >= 4 for p in chosen)), None)
            radial.require(choice is not None, 'insufficient separated sites')
            chosen.append(choice)
        for index, (x, y) in enumerate(chosen):
            site = {'name': f'fresh_r{nominal}_b{index}', 'row': x, 'column': y, 'nominal_radius': nominal,
                    'azimuth_rank': index, 'azimuth_degrees': float(azimuth[y, x])}
            mask = holdout_mask(science.shape, protocol, site)
            candidate = support_counts(science, (x, y), mask)
            calibration_counts = {t['name']: support_counts(science, (t['row'], t['column']), mask) for t in calibration}
            radial.require(all(c['same_radius'] >= 8 for records in calibration_counts.values() for c in records),
                           'site cannot retain the complete 28-search calibration set')
            fresh_mask = np.zeros(science.shape, bool)
            mark_footprint(fresh_mask, x, y)
            radial.require(not np.any(fresh_mask & old_measurements), 'new measurement footprint overlaps previous measurement')
            sites.append(site)
            audits.append({'site': site['name'], 'candidate': candidate, 'calibration': calibration_counts})
    protocol['sites'] = sites
    radial.require(len(sites) == 30 and len({(s['row'], s['column']) for s in sites}) == 30, 'changed study size')
    root.mkdir(parents=True, exist_ok=False)
    (root/'payload/response').mkdir(parents=True)
    (root/'software').mkdir()
    for source in field.glob('p4PSF_*.fits'):
        shutil.copy2(source, root/'payload/response'/source.name)
    for name in ('reduction.conf', 'injection_psf.fits', 'common_command.json'):
        shutil.copy2(args.source_batch/name, root/'payload'/name)
    shutil.copy2(old['science']['path'], root/'payload/previous_baseline.fits')
    for name in SCRIPTS:
        shutil.copy2(Path(__file__).with_name(name), root/'software'/name)
    write_json(root/'protocol.json', protocol)
    write_json(root/'geometry.json', {'selection_used_scores': False, 'minimum_same_ring_center_spacing': 4,
        'new_and_previous_measurement_footprints_disjoint': True, 'viable_pool_counts': pool_counts, 'sites': audits})
    write_json(root/'input_records.json', batch['input_records'])
    staged = [fingerprint(p) for directory in ('payload', 'software') for p in sorted((root/directory).rglob('*')) if p.is_file()]
    staged += [fingerprint(root/n) for n in ('protocol.json', 'geometry.json', 'input_records.json')]
    write_json(root/'staging.json', {'root': str(root), 'records': staged,
        'prior_science': old['science'], 'prior_response': old['source_response_inputs'],
        'prior_protocol': fingerprint(args.previous_study/'protocol.json'),
        'prior_sampling': fingerprint(args.previous_sampling/'manifest.json'),
        'source_batch': fingerprint(args.source_batch/'manifest.json')})
    print(f'Prepared {len(sites)} sites and {len(sites)*len(LEVELS)} positive reductions in {root}', flush=True)


def prepare(args: argparse.Namespace) -> None:
    """Freeze ROC-native executables, shared libraries, and the complete immutable run contract."""
    root = args.root.resolve()
    radial.require(not (root/'manifest.json').exists(), 'prepared run already exists')
    staging = read(root/'staging.json')
    verify(staging['records'])
    inputs = read(root/'input_records.json')
    verify(inputs)
    executables = {'p4ReductionPrecisionBenchmark': args.build/'benchmarks/p4ReductionPrecisionBenchmark',
                   'hciGaussianReference': args.build/'benchmarks/hciGaussianReference', 'hciAnalyze': args.build/'src/hciAnalyze'}
    dependency_sources = []
    for name, source in executables.items():
        shutil.copy2(source, root/'software'/name)
        linked = subprocess.check_output(['ldd', str(source)], text=True)
        radial.require('not found' not in linked, 'missing native binary dependency')
        (root/'payload'/f'{name}.ldd.txt').write_text(linked)
        for soname, filename in re.findall(r'^\s*(\S+)\s+=>\s+(/\S+)', linked, re.M):
            destination = root/'software'/soname
            source_record = fingerprint(Path(filename))
            if destination.exists():
                radial.require(fingerprint(destination)['sha256'] == source_record['sha256'], 'conflicting dependency '+soname)
            else:
                shutil.copy2(filename, destination)
            dependency_sources.append(source_record)
    input_list = root/'payload/inputs.txt'
    input_list.write_text(''.join(record['path']+'\n' for record in inputs))
    command = read(root/'payload/common_command.json')
    command[0] = str(root/'software/p4ReductionPrecisionBenchmark')
    command[command.index('--config')+1] = str(root/'payload/reduction.conf')
    command = replacement(command, {'input.fileList': input_list, 'fake.fileName': root/'payload/injection_psf.fits'})
    write_json(root/'reduction_command.json', command)
    protocol = read(root/'protocol.json')
    radial.require(set(protocol['reduction_cpu_ids']).issubset(os.sched_getaffinity(0)), 'unavailable CPU affinity')
    environment = binary_environment(root, reduction=False)
    for name in ('p4ReductionPrecisionBenchmark', 'hciAnalyze'):
        with (root/'payload'/f'{name}.help.txt').open('w') as stream:
            subprocess.run([str(root/'software'/name), '--help'], env=environment, stdout=stream,
                           stderr=subprocess.STDOUT, check=True)
    records = [fingerprint(p) for directory in ('payload', 'software') for p in sorted((root/directory).rglob('*'))
               if p.is_file() and '__pycache__' not in p.parts]
    records += [fingerprint(root/n) for n in ('protocol.json', 'geometry.json', 'input_records.json', 'staging.json', 'reduction_command.json')]
    write_json(root/'manifest.json', {'schema': 1, 'source_cpp_commit': protocol['source_cpp_commit'],
        'frozen_records': records, 'input_records': inputs, 'native_dependency_sources': dependency_sources,
        'build_cache': fingerprint(args.build/'CMakeCache.txt'),
        'source_archive': fingerprint(root/'source.tar.gz'), 'host': os.uname().nodename,
        'python': {'executable': fingerprint(Path(sys.executable)), 'numpy': np.__version__, 'scipy': scipy.__version__,
                   'astropy': astropy.__version__, 'matplotlib': matplotlib.__version__},
        'separate_roc_pilot': False, 'positive_reductions': 90, 'baseline_reductions': 1})
    write_json(root/'state.json', {'status': 'prepared', 'positive_reductions': 90, 'finished': []})
    print('Native dependencies, all 621 input hashes, and 30-site geometry are ready.', flush=True)


def binary_environment(root: Path, reduction: bool) -> dict:
    """Use frozen binary dependencies and bounded, non-nested CPU threading."""
    protocol = read(root/'protocol.json')
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(protocol['reduction_openmp_threads'] if reduction else protocol['analysis_openmp_threads']),
                       OMP_PROC_BIND='true', OMP_PLACES='cores', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       LD_LIBRARY_PATH=str(root/'software'))
    environment.pop('P4REDUCE_GLOBAL_CONFIG', None)
    environment.pop('HCIANALYZE_GLOBAL_CONFIG', None)
    return environment


def reduce_image(root: Path, trial: dict) -> Path:
    """Run or reuse a complete full-image reduction without overwriting partial work."""
    directory = root/'reductions'/trial['name']
    complete = directory/'complete.json'
    if complete.exists():
        verify(read(complete)['products'])
        return directory/'finim.fits'
    directory.mkdir(parents=True)
    protocol = read(root/'protocol.json')
    x, y = trial['row']-127.5, trial['column']-127.5
    command = replacement(read(root/'reduction_command.json'), {'fake.sep': math.hypot(x, y),
        'fake.PA': math.degrees(-math.atan2(x, y)) % 360, 'fake.contrast': trial['contrast'], 'output.directory': directory})
    command = ['taskset', '-c', ','.join(map(str, protocol['reduction_cpu_ids'])), *command]
    write_json(directory/'command.json', command)
    start = time.monotonic()
    with (directory/'run.log').open('w') as log:
        subprocess.run(['/usr/bin/time', '-f', 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M',
                        '-o', str(directory/'resource_usage.txt'), *command], env=binary_environment(root, True), cwd=directory,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    data, header = fits.getdata(directory/'finim.fits', header=True)
    radial.require(data.shape in ((256, 256), (1, 256, 256)) and header['P4 LOCAL STAMP SIZE'] == 0 and
                   header['COMBINATION METHOD'].strip() == 'mean', 'changed full-image reduction contract')
    write_json(complete, {'trial': trial, 'elapsed_seconds': time.monotonic()-start,
                          'products': [fingerprint(p) for p in sorted(directory.glob('*.fits'))]})
    return directory/'finim.fits'


def reference_maps(root: Path, source: Path, name: str) -> dict:
    """Generate production Gaussian and ordinary annular-SNR references, with exact replay checks."""
    directory = root/'references'/name
    complete = directory/'complete.json'
    if complete.exists():
        verify(read(complete)['products'])
    else:
        directory.parent.mkdir(exist_ok=True)
        gaussian.analyze_image(directory, source, root/'software', root/'payload/response/p4PSF_manifest.fits',
                               binary_environment(root, False))
    return gaussian.maps(directory)


def candidate_models(science: np.ndarray, position: tuple, mask: np.ndarray) -> dict:
    """Estimate the four fixed covariance models using only accepted training pixels."""
    rings = radial.geometry(science, position, mask)
    matrices = {o: radial.extract(science, r) for o, r in rings.items() if abs(o) <= 5}
    for o, ring in rings.items():
        if abs(o) <= 5:
            used = ring['indices'][ring['weights'] != 0]
            radial.require(not np.any(mask.ravel()[used]), 'held-out native pixel enters covariance training')
    band = np.vstack(list(matrices.values()))
    models = {'psd_hann_b0_m0.1': psd.fit_psd(matrices[0], 'hann', .1),
              'psd_rectangular_b5_m0.3': psd.fit_psd(band, 'rectangular', .3),
              'pca_b5_f1': radial.fit(band, 1), 'isotropic_b5': shrinkage.fit_shrinkage(band, 1)}
    for name, model in models.items():
        if model is not None:
            model['samples'] = len(matrices[0]) if name == 'psd_hann_b0_m0.1' else len(band)
    return models


def measure_trial(science: np.ndarray, trial: dict, templates: dict, mask: np.ndarray, maps: dict) -> dict:
    """Keep raw exact-center photometry separate from the five-pixel search statistic."""
    x, y = trial['row'], trial['column']
    pixels = {name: [] for name in (*COVARIANCE_MODELS, 'identity')}
    for ox, oy in radial.OFFSETS:
        position = (x+ox, y+oy)
        template = templates[position]
        data = science[y+oy-5:y+oy+6, x+ox-5:x+ox+6].ravel()
        radial.require(np.isfinite(data).all(), 'invalid candidate science support')
        models = candidate_models(science, position, mask)
        for name, model in models.items():
            if model is None:
                pixels[name].append({'valid': False})
                continue
            result = radial.filter_stamp(data, template, model, np.ones(121))
            w = result['physical_weight']
            radial.require(np.isclose(w @ template, 1, rtol=1e-12, atol=0) and
                           np.isclose(w @ model['covariance'] @ w, result['sigma']**2, rtol=1e-10, atol=0), 'invalid covariance normalization')
            pixels[name].append({'valid': True, 'training_samples': model['samples'],
                                 **{k: result[k] for k in ('amplitude', 'sigma', 'score')}})
        energy = float(template @ template)
        amplitude = float(template @ data/energy)
        sigma = 1/math.sqrt(energy)
        pixels['identity'].append({'valid': True, 'amplitude': amplitude, 'sigma': sigma, 'score': amplitude/sigma})
    result = {}
    for name, values in pixels.items():
        valid = all(p['valid'] for p in values)
        result[name] = {'valid': valid, 'center': values[0], 'pixels': values,
                        'search_score': max(p['score'] for p in values) if valid else None}
    for name in ('gaussian_snr', 'gaussian_raw', 'identity_snr'):
        result[name] = {'valid': True, **gaussian.score(maps[name], trial)}
    return result


def calibrate(root: Path, science: np.ndarray, templates: dict, maps: dict) -> None:
    """Freeze every threshold and source contrast before inspecting new-site baseline scores or positives."""
    complete = root/'calibration_complete.json'
    if complete.exists():
        verify(read(complete)['products'])
        return
    radial.require(not (root/'thresholds.json').exists(), 'partial calibration exists; preserve it before any restart')
    protocol = read(root/'protocol.json')
    thresholds, jobs, calibration = {}, [], []
    for site in protocol['sites']:
        write_json(root/'state.json', {'status': 'calibrating', 'site': site['name'], 'pid': os.getpid(), 'finished': []})
        mask = holdout_mask(science.shape, protocol, site)
        rows = [{'trial': t, 'models': measure_trial(science, t, templates, mask, maps)} for t in protocol['calibration_trials']]
        one = {}
        for name in MODELS:
            radial.require(all(r['models'][name]['valid'] for r in rows), 'changed 28-search calibration support')
            scores = sorted(r['models'][name]['search_score'] for r in rows)
            rank = math.ceil((len(scores)+1)*.95)
            one[name] = scores[rank-1]
        thresholds[site['name']] = one
        position = (site['row'], site['column'])
        model = candidate_models(science, position, mask)['isotropic_b5']
        radial.require(model is not None, 'invalid contrast-reference covariance')
        # Supply zero data here: the source scale uses only covariance and template, never target pixels.
        sigma = radial.filter_stamp(np.zeros(121), templates[position], model, np.ones(121))['sigma']
        reference = one['isotropic_b5']*sigma
        radial.require(np.isfinite(reference) and reference > 0, 'nonpositive calibration-based source scale')
        calibration.append({'site': site['name'], 'calibration_trials': rows, 'reference_contrast_scale': reference,
                            'reference_sigma': sigma, 'target_values_used_for_contrast': False})
        for index, level in enumerate(LEVELS):
            jobs.append({**site, 'site': site['name'], 'name': site['name']+f'_l{index}',
                         'brightness_multiplier': level, 'reference_contrast_scale': reference, 'contrast': level*reference})
        print('calibrated '+site['name'], flush=True)
    write_json(root/'thresholds.json', thresholds)
    write_json(root/'jobs.json', jobs)
    write_json(root/'calibration.json', calibration)
    write_json(complete, {'products': [fingerprint(root/name) for name in ('thresholds.json', 'jobs.json', 'calibration.json')],
                          'all_thresholds_and_contrasts_frozen_before_new_site_measurements': True})


def summarize(root: Path, nulls: list, rows: list, thresholds: dict) -> None:
    """Persist recovery, false positives, contrast error, and descriptive plots for every fixed method."""
    groups = []
    model_stats = {}
    for name in MODELS:
        model_stats[name] = {'evaluation_trials': len(nulls), 'evaluation_exceedances': sum(n['models'][name]['detected'] for n in nulls),
                             'invalid_null_searches': sum(not n['models'][name]['valid'] for n in nulls)}
        for level in LEVELS:
            selected = [r for r in rows if r['trial']['brightness_multiplier'] == level]
            measured = [r['models'][name] for r in selected]
            errors = [v['raw_contrast_error'] for v in measured if 'raw_contrast_error' in v]
            groups.append({'model': name, 'brightness_multiplier': level, 'trials': len(selected),
                'detections': sum(v['detected'] for v in measured), 'invalid_searches': sum(not v['valid'] for v in measured),
                'median_raw_contrast_error': float(np.median(errors)) if errors else None,
                'median_absolute_raw_contrast_error': float(np.median(np.abs(errors))) if errors else None,
                'rms_raw_contrast_error': float(np.sqrt(np.mean(np.array(errors)**2))) if errors else None,
                'within_one_conditional_sigma': sum(v.get('inside_one_conditional_sigma', False) for v in measured) if errors else None})
    write_json(root/'results.json', {'models': model_stats, 'groups': groups, 'measurements': rows, 'nulls': nulls,
        'thresholds': thresholds, 'caveats': ['30 new source positions in a reused residual field; overlapping sites and shared calibration/training are correlated.',
        'Each source has its own full-footprint holdout and recalibration; this changes the former global exclusion union.',
        'Gaussian/identity application-SNR references retain legacy whole-annulus normalization including trial neighborhoods.',
        'Identity conditional sigma uses unit pixel covariance and is not a fitted physical uncertainty.']})
    lines = ['# Full ROC PSD study', '', 'One baseline and 90 full-image injections; no separate ROC pilot.', '',
             '| Method | 0.5× | 0.75× | 1× | Evaluation null exceedances |', '| --- | --- | --- | --- | --- |']
    for name in MODELS:
        selected = [g for g in groups if g['model'] == name]
        lines.append('| '+name+' | '+' | '.join(f'{g["detections"]}/{g["trials"]}' for g in selected)+
                     f' | {model_stats[name]["evaluation_exceedances"]}/{len(nulls)} |')
    lines += ['', 'All thresholds and injection contrasts were fixed from calibration before new-site scores or positives. '
              'Invalid searches remain nondetections. The 30 sites are correlated; counts do not establish independent-trial confidence intervals.', '',
              'See `results.json` for individual raw photometry and conditional errors, and `protocol.json` for the fixed holdouts and controls.', '']
    (root/'results.md').write_text('\n'.join(lines))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout='constrained')
    for name in MODELS:
        selected = [g for g in groups if g['model'] == name]
        axes[0].plot(LEVELS, [g['detections']/g['trials'] for g in selected], 'o-', label=name)
    axes[0].set(xlabel='Brightness / calibration-based reference contrast', ylabel='Recovery fraction', ylim=(0, 1.05))
    axes[0].legend(fontsize=7)
    axes[1].barh(MODELS, [model_stats[name]['evaluation_exceedances'] for name in MODELS])
    axes[1].set(xlabel='Null exceedances out of 30', title='Site-specific calibration thresholds')
    fig.suptitle('Full ROC study: 30 new sites × three brightnesses; one correlated field')
    fig.savefig(root/'comparison.png', dpi=170)
    plt.close(fig)


def run(root: Path) -> None:
    """Complete the immutable baseline/calibration/90-injection queue under a single supervisor lock."""
    manifest = read(root/'manifest.json')
    protocol = read(root/'protocol.json')
    verify([*manifest['frozen_records'], *manifest['input_records']])
    templates = load_templates(root/'payload/response')
    finished = []
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (root/'complete.json').exists():
            verify(read(root/'complete.json')['products'])
            print('Full study already complete.', flush=True)
            return
        try:
            write_json(root/'state.json', {'status': 'reducing_baseline', 'pid': os.getpid(), 'finished': finished})
            source = reduce_image(root, {'name': 'baseline', 'row': 127, 'column': 127, 'contrast': 0.})
            science = fits.getdata(source).squeeze().astype(float)
            previous = fits.getdata(root/'payload/previous_baseline.fits').squeeze()
            radial.require(np.array_equal(np.isfinite(science), np.isfinite(previous)), 'full-study baseline changed finite support')
            maps = reference_maps(root, source, 'baseline')
            calibrate(root, science, templates, maps)
            calibration_record = read(root/'calibration_complete.json')
            thresholds = read(root/'thresholds.json')
            jobs = read(root/'jobs.json')
            if (root/'baseline_nulls_complete.json').exists():
                verify(read(root/'baseline_nulls_complete.json')['products'])
                nulls = read(root/'baseline_nulls.json')
            else:
                radial.require(not (root/'baseline_nulls.json').exists(), 'partial baseline-null record exists')
                nulls = []
                for site in protocol['sites']:
                    measured = measure_trial(science, site, templates, holdout_mask(science.shape, protocol, site), maps)
                    for name, value in measured.items():
                        value['detected'] = value['valid'] and value['search_score'] > thresholds[site['name']][name]
                    nulls.append({'site': site, 'models': measured})
                write_json(root/'baseline_nulls.json', nulls)
                write_json(root/'baseline_nulls_complete.json', {'products': [fingerprint(root/'baseline_nulls.json')]})
            null_record = fingerprint(root/'baseline_nulls.json')
            measurements = []
            for job in jobs:
                record_path = root/'measurements'/(job['name']+'.json')
                if record_path.exists():
                    record = read(record_path)
                    verify(record['products'])
                else:
                    write_json(root/'state.json', {'status': 'reducing', 'current_trial': job['name'], 'pid': os.getpid(), 'finished': finished})
                    source = reduce_image(root, job)
                    positive = fits.getdata(source).squeeze().astype(float)
                    radial.require(np.array_equal(np.isfinite(science), np.isfinite(positive)), 'injection changed finite support')
                    write_json(root/'state.json', {'status': 'analyzing', 'current_trial': job['name'], 'pid': os.getpid(), 'finished': finished})
                    maps = reference_maps(root, source, job['name'])
                    site = next(s for s in protocol['sites'] if s['name'] == job['site'])
                    measured = measure_trial(positive, job, templates, holdout_mask(science.shape, protocol, site), maps)
                    baseline = next(n for n in nulls if n['site']['name'] == job['site'])
                    for name, value in measured.items():
                        value['detected'] = value['valid'] and value['search_score'] > thresholds[job['site']][name]
                        if name in (*COVARIANCE_MODELS, 'identity') and value['center']['valid']:
                            center = value['center']
                            value['raw_contrast_error'] = center['amplitude']/job['contrast']-1
                            value['inside_one_conditional_sigma'] = abs(center['amplitude']-job['contrast']) <= center['sigma']
                            before = baseline['models'][name]['center']
                            if before['valid']:
                                value['paired_increment_error'] = (center['amplitude']-before['amplitude'])/job['contrast']-1
                    record = {'trial': job, 'models': measured, 'products': [fingerprint(source),
                        *read(root/'references'/job['name']/'complete.json')['products']]}
                    record_path.parent.mkdir(exist_ok=True)
                    write_json(record_path, record)
                measurements.append(record)
                finished.append(job['name'])
                write_json(root/'state.json', {'status': 'running', 'pid': os.getpid(), 'finished': finished, 'total': len(jobs)})
                print(f'completed {len(finished)}/{len(jobs)}: {job["name"]}', flush=True)
            summarize(root, nulls, measurements, thresholds)
            verify([*manifest['frozen_records'], *manifest['input_records'], *calibration_record['products'], null_record])
            write_json(root/'complete.json', {'baseline_reductions': 1, 'positive_reductions': len(finished),
                'finished': finished, 'frozen_inputs_unchanged': True,
                'products': [fingerprint(root/name) for name in ('results.json', 'results.md', 'comparison.png', 'baseline_nulls.json')]})
            write_json(root/'state.json', {'status': 'complete', 'finished': finished, 'results': str(root/'results.json')})
        except Exception as error:
            write_json(root/'state.json', {'status': 'failed', 'pid': os.getpid(), 'finished': finished, 'error': str(error)})
            raise


def main() -> None:
    """Expose separate score-free design, native preparation, and unattended execution phases."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('design', 'prepare', 'run'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--previous-study', type=Path)
    parser.add_argument('--previous-sampling', type=Path)
    parser.add_argument('--source-batch', type=Path)
    parser.add_argument('--source-commit', default='882cb15')
    parser.add_argument('--build', type=Path)
    args = parser.parse_args()
    if args.action == 'design':
        if any(getattr(args, key) is None for key in ('previous_study', 'previous_sampling', 'source_batch')):
            parser.error('design requires previous-study, previous-sampling, and source-batch')
        design(args)
    elif args.action == 'prepare':
        if args.build is None:
            parser.error('prepare requires build')
        prepare(args)
    else:
        run(args.root.resolve())


if __name__ == '__main__':
    main()
