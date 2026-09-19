#!/usr/bin/env python3
"""Prepare and run the inner-radius pooled-rectangular PSD injection study.

The score-free design compares radial training half-widths 5, 10, and 20
pixels at six radii.  Every filter input and noise sample is kept clear of the
known planet.  Setup freezes the design, software, and reduction inputs; run
freezes baseline thresholds and injection contrasts before reducing positives.
"""
from __future__ import annotations

import argparse
import fcntl
import math
import os
from pathlib import Path
import shutil
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.ndimage import binary_dilation

import run_p4_step5_roc_full as full
from measure_p4_step5_planet import annular_oracle
from run_p4_step5_full_injections import fingerprint, replacement, write_json

RADII = (6, 8, 12, 16, 20, 24)
WIDTHS = (5, 10, 20)
LEVELS = (.5, .75, 1.)
RECTANGULAR = {width: f'psd_rectangular_b{width}_m0.3' for width in WIDTHS}
METHODS = (*RECTANGULAR.values(), 'identity', 'gaussian')
METHOD_WIDTH = {RECTANGULAR[width]: width for width in WIDTHS}
CALIBRATION_SEARCHES = 20
SEARCH_OFFSETS = tuple(full.radial.OFFSETS)
SCRIPT_NAMES = (*full.SCRIPTS, 'measure_p4_step5_planet.py', 'run_p4_step5_inner_rectangular.py')


def planet_masks(shape: tuple, circle: list) -> tuple[np.ndarray, np.ndarray]:
    """Return the production source mask and centers with a clean 15-square kernel."""
    yy, xx = np.indices(shape)
    x, y, radius = circle
    source = np.hypot(xx-x, yy-y) <= radius+.5
    contaminated = binary_dilation(source, structure=np.ones((15, 15), bool))
    return source, ~contaminated


def search_is_planet_clean(clean_centers: np.ndarray, x: int, y: int) -> bool:
    """Require all five detection kernels to avoid the buffered planet exclusion."""
    return all(clean_centers[y+dy, x+dx] for dx, dy in SEARCH_OFFSETS)


def active_methods(trial: dict) -> list[str]:
    """Return references and covariance models with a complete five-pixel search."""
    return [name for name in METHODS if name not in METHOD_WIDTH or
            trial['support'][str(METHOD_WIDTH[name])]['valid']]


def support(science: np.ndarray, position: tuple, source_mask: np.ndarray) -> dict:
    """Count accepted patches for each width across the complete five-pixel search."""
    x, y = position
    forbidden = source_mask.copy()
    full.mark_footprint(forbidden, x, y)
    counts = {width: [] for width in WIDTHS}
    for dx, dy in SEARCH_OFFSETS:
        rings = full.radial.geometry(science, (x+dx, y+dy), forbidden)
        for width in WIDTHS:
            counts[width].append(sum(len(ring['halves']) for offset, ring in rings.items()
                                     if abs(offset) <= width))
    return {str(width): {'pixels': values, 'minimum': min(values), 'valid': min(values) >= 8}
            for width, values in counts.items()}


def maximin(rows: list, count: int) -> list:
    """Choose a deterministic geometry-only set with large pairwise separations."""
    full.radial.require(len(rows) >= count, 'insufficient candidates for maximin selection')
    ordered = sorted(rows, key=lambda row: (row['azimuth_degrees'], row['row'], row['column']))
    best, best_key = None, None
    for seed in ordered:
        chosen = [seed]
        while len(chosen) < count:
            remaining = [row for row in ordered if row not in chosen]
            choice = max(remaining, key=lambda row: (min(math.dist((row['row'], row['column']),
                (prior['row'], prior['column'])) for prior in chosen), -row['azimuth_degrees']))
            chosen.append(choice)
        distances = sorted(math.dist((a['row'], a['column']), (b['row'], b['column']))
                           for index, a in enumerate(chosen) for b in chosen[index+1:])
        key = tuple(distances)
        if best_key is None or key > best_key:
            best, best_key = chosen, key
    return sorted(best, key=lambda row: row['azimuth_degrees'])


def candidate_rows(science: np.ndarray, templates: dict, circle: list) -> dict:
    """Enumerate score-free, planet-clean candidates and exact training support."""
    source_mask, clean_centers = planet_masks(science.shape, circle)
    yy, xx = np.indices(science.shape)
    radius = np.hypot(xx-127.5, yy-127.5)
    azimuth = np.degrees(np.arctan2(yy-127.5, xx-127.5)) % 360
    result = {}
    for nominal in range(5, 31):
        rows = []
        for y, x in np.argwhere(abs(radius-nominal) < .5):
            x, y = int(x), int(y)
            if not search_is_planet_clean(clean_centers, x, y):
                continue
            if not all((x+dx, y+dy) in templates for dx, dy in SEARCH_OFFSETS):
                continue
            if any(not np.isfinite(science[y+dy-5:y+dy+6, x+dx-5:x+dx+6]).all()
                   for dx, dy in SEARCH_OFFSETS):
                continue
            rows.append({'row': x, 'column': y, 'nominal_radius': nominal,
                'actual_radius': float(radius[y, x]), 'azimuth_degrees': float(azimuth[y, x]),
                'support': support(science, (x, y), source_mask)})
        result[nominal] = rows
        print(f'design radius {nominal}: {len(rows)} planet-clean candidates', flush=True)
    return result


def calibration_pool(candidates: dict, sites: list, nominal: int, method: str) -> tuple[list, int]:
    """Choose 20 null searches from the narrowest score-free radial band with support."""
    site_coordinates = {(site['row'], site['column']) for site in sites}
    width = METHOD_WIDTH.get(method, 20)
    for half_width in range(0, 26):
        eligible = [row for radius, rows in candidates.items() if abs(radius-nominal) <= half_width
                    for row in rows if (row['row'], row['column']) not in site_coordinates and
                    row['support'][str(width)]['valid']]
        if len(eligible) >= CALIBRATION_SEARCHES:
            return maximin(eligible, CALIBRATION_SEARCHES), half_width
    raise RuntimeError(f'cannot calibrate {method} at radius {nominal}')


def design(parent: Path, cpu_ids: list[int]) -> tuple[dict, dict]:
    """Build the complete geometry-only protocol without reading any detection score."""
    full.radial.require(bool(cpu_ids) and len(cpu_ids) == len(set(cpu_ids)) and min(cpu_ids) >= 0,
                        'reduction CPU IDs must be nonnegative and distinct')
    parent_protocol = full.read(parent/'protocol.json')
    science_path = parent/'reductions/baseline/finim.fits'
    science = fits.getdata(science_path).squeeze().astype(float)
    templates = full.load_templates(parent/'payload/response')
    candidates = candidate_rows(science, templates, parent_protocol['known_source_circle'])
    sites = []
    for nominal in RADII:
        viable = [row for row in candidates[nominal] if row['support']['20']['valid']]
        chosen = maximin(viable, 6)
        for index, row in enumerate(chosen):
            sites.append({**row, 'name': f'inner_r{nominal}_b{index}', 'azimuth_rank': index})
    pools, unique_nulls = {}, {}
    for nominal in RADII:
        active = {'identity', 'gaussian', RECTANGULAR[20]}
        for width in (5, 10):
            if any(site['nominal_radius'] == nominal and site['support'][str(width)]['valid'] for site in sites):
                active.add(RECTANGULAR[width])
        pools[str(nominal)] = {}
        for method in METHODS:
            if method not in active:
                pools[str(nominal)][method] = None
                continue
            selected, half_width = calibration_pool(candidates, sites, nominal, method)
            names = []
            for row in selected:
                key = (row['row'], row['column'])
                name = f'null_x{row["row"]}_y{row["column"]}'
                unique_nulls.setdefault(key, {**row, 'name': name})
                names.append(name)
            pools[str(nominal)][method] = {'trials': names, 'radial_half_width': half_width,
                'radius_range': [nominal-half_width, nominal+half_width]}
    by_name = {row['name']: row for row in unique_nulls.values()}
    full.radial.require(all(name in by_name for radius in pools.values() for pool in radius.values()
                            if pool is not None for name in pool['trials']), 'missing calibration trial')
    spacing = {str(radius): min(math.dist((a['row'], a['column']), (b['row'], b['column']))
        for index, a in enumerate([s for s in sites if s['nominal_radius'] == radius])
        for b in [s for s in sites if s['nominal_radius'] == radius][index+1:]) for radius in RADII}
    protocol = {'schema': 1, 'purpose': 'inner-radius recovery comparison for pooled rectangular PSD widths',
        'parent_study': str(parent.resolve()), 'models': list(METHODS), 'nominal_radii': list(RADII),
        'sites_per_radius': 6, 'sites': sites, 'brightness_multipliers': list(LEVELS),
        'full_positive_reductions': len(sites)*len(LEVELS), 'baseline_reductions': 0,
        'known_source_circle': parent_protocol['known_source_circle'], 'lambda_d_pixels': 3.6,
        'planet_guard': 'production planet circle plus 0.5-pixel boundary; the 15x15 kernel at every one of the five search pixels must be disjoint',
        'annular_noise': 'complete required one-pixel annuli, but only filter centers whose 15x15 kernel is disjoint from the buffered known-planet circle; production mean/stddev interpolation and small-sample correction; injection neighborhoods remain in normalization',
        'training': 'raw aligned 11x11 patches; five-pixel angular/radial center step; the trial full five-search 15x15 kernel union and buffered known-planet circle excluded; rectangular PSD; mixture0.3; minimum8',
        'holdout': 'per search trial: buffered known-planet circle plus the full five-search 15x15 kernel union; no union of unrelated calibration footprints',
        'search': 'native center plus four axial one-pixel neighbors; all five required; invalid searches are nondetections',
        'calibration_searches_per_active_radius_method': CALIBRATION_SEARCHES,
        'calibration_pools': pools, 'calibration_trials': sorted(unique_nulls.values(), key=lambda r: r['name']),
        'threshold_rule': 'maximum of 20 geometry-selected planet-clean baseline search scores; strict exceedance',
        'reference_contrast': 'identity matched-filter annular-SNR threshold times interpolated identity-amplitude radial sigma divided by the small-sample correction; target pixels unused',
        'selection': 'six deterministic maximin sites per requested radius from planet-clean candidates with complete band20 support; scores and positive images unused; narrower bands may be invalid',
        'dependence': 'same residual field, overlapping training patches, radial null pools, and close inner sites are correlated; counts are descriptive',
        'reduction_cpu_ids': cpu_ids, 'reduction_openmp_threads': len(cpu_ids), 'blas_threads': 1,
        'analysis_cpu_ids': [12, 13], 'analysis_openmp_threads': 2,
        'source_cpp_commit': parent_protocol['source_cpp_commit']}
    geometry = {'scores_or_recovery_used': False, 'candidate_radii': [5, 30],
        'planet_clean_candidates': {str(radius): len(rows) for radius, rows in candidates.items()},
        'minimum_site_spacing_by_radius': spacing, 'sites': sites, 'calibration_pools': pools,
        'unique_calibration_trials': len(unique_nulls)}
    return protocol, geometry


def audit(args: argparse.Namespace) -> None:
    """Write only the score-free design for local review."""
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    protocol, geometry = design(args.parent.resolve(), args.cpus)
    write_json(root/'protocol.json', protocol)
    write_json(root/'geometry.json', geometry)
    print(f'audited {len(protocol["sites"])} sites and {len(protocol["calibration_trials"])} unique null centers', flush=True)


def setup(args: argparse.Namespace) -> None:
    """Freeze the score-free design, native software, response, and reduction inputs."""
    root, parent = args.root.resolve(), args.parent.resolve()
    full.radial.require(full.read(parent/'state.json')['status'] == 'complete', 'parent study is incomplete')
    parent_manifest = full.read(parent/'manifest.json')
    baseline_complete = full.read(parent/'reductions/baseline/complete.json')
    full.radial.require(set(args.cpus).issubset(os.sched_getaffinity(0)), 'requested reduction CPUs unavailable')
    full.verify([*parent_manifest['frozen_records'], *parent_manifest['input_records'],
                 *full.read(parent/'complete.json')['products'], *baseline_complete['products']])
    protocol, geometry = design(parent, args.cpus)
    root.mkdir(parents=True, exist_ok=False)
    (root/'payload/response').mkdir(parents=True)
    (root/'software').mkdir()
    for source in sorted((parent/'payload/response').glob('*.fits')):
        shutil.copy2(source, root/'payload/response'/source.name)
    for name in ('reduction.conf', 'injection_psf.fits', 'inputs.txt'):
        shutil.copy2(parent/'payload'/name, root/'payload'/name)
    shutil.copy2(parent/'reductions/baseline/finim.fits', root/'payload/baseline.fits')
    for source in sorted((parent/'software').iterdir()):
        if source.is_file() and source.suffix != '.py':
            shutil.copy2(source, root/'software'/source.name)
    for name in SCRIPT_NAMES:
        shutil.copy2(Path(__file__).with_name(name), root/'software'/name)
    command = full.read(parent/'reduction_command.json')
    command[0] = str(root/'software'/Path(command[0]).name)
    command[command.index('--config')+1] = str(root/'payload/reduction.conf')
    command = replacement(command, {'input.fileList': root/'payload/inputs.txt',
        'fake.fileName': root/'payload/injection_psf.fits'})
    write_json(root/'protocol.json', protocol)
    write_json(root/'geometry.json', geometry)
    write_json(root/'reduction_command.json', command)
    write_json(root/'input_records.json', parent_manifest['input_records'])
    records = [fingerprint(path) for directory in ('payload', 'software')
               for path in sorted((root/directory).rglob('*')) if path.is_file()]
    records += [fingerprint(root/name) for name in ('protocol.json', 'geometry.json',
        'reduction_command.json', 'input_records.json')]
    write_json(root/'manifest.json', {'schema': 1, 'parent_manifest': fingerprint(parent/'manifest.json'),
        'parent_complete': fingerprint(parent/'complete.json'), 'frozen_records': records,
        'parent_baseline_complete': fingerprint(parent/'reductions/baseline/complete.json'),
        'input_records': parent_manifest['input_records'], 'host': os.uname().nodename})
    write_json(root/'state.json', {'status': 'prepared', 'positive_reductions': len(protocol['sites'])*len(LEVELS),
        'finished': []})
    print(f'prepared {len(protocol["sites"])} sites and {len(protocol["sites"])*len(LEVELS)} injections in {root}', flush=True)


def repair(args: argparse.Namespace) -> None:
    """Update a pre-calibration frozen runner after the unsupported-layer abort."""
    root = args.root.resolve()
    manifest = full.read(root/'manifest.json')
    state = full.read(root/'state.json')
    full.radial.require(state['status'] == 'failed' and 'hciAnalyze' in state.get('error', '') and
                        'SIGABRT' in state.get('error', '') and not (root/'calibration_complete.json').exists() and
                        not (root/'thresholds.json').exists() and not (root/'jobs.json').exists(),
                        'repair is allowed only after the pre-calibration unsupported-layer failure')
    reductions = root/'reductions'
    full.radial.require(not reductions.exists() or not any(reductions.iterdir()),
                        'positive reductions already exist; do not repair this study in place')
    target = (root/'software'/Path(__file__).name).resolve()
    old = next((record for record in manifest['frozen_records'] if Path(record['path']).resolve() == target), None)
    full.radial.require(old is not None and fingerprint(target) == old, 'frozen runner changed before repair')
    full.radial.require(fingerprint(Path(__file__).resolve())['sha256'] != old['sha256'],
                        'repair source is not a newer runner')
    full.verify([record for record in manifest['frozen_records'] if record != old])
    full.verify([*manifest['input_records'], *manifest.get('repair_records', [])])
    archived = []
    analysis = root/'analysis'
    if analysis.exists():
        archive = root/'pre_repair_failures'
        for directory in sorted(path for path in analysis.iterdir() if path.is_dir() and
                                not (path/'complete.json').exists()):
            archive.mkdir(exist_ok=True)
            destination = archive/directory.name
            full.radial.require(not destination.exists(), 'pre-repair archive already exists: '+str(destination))
            shutil.move(directory, destination)
            archived.append(str(destination))
    previous_state = fingerprint(root/'state.json')
    shutil.copy2(Path(__file__).resolve(), target)
    updated = fingerprint(target)
    manifest['frozen_records'] = [updated if record == old else record for record in manifest['frozen_records']]
    record = {'schema': 1, 'reason': 'omit geometrically unsupported covariance layers from production SNR interpolation',
        'pre_calibration': True, 'positive_reductions_before_repair': 0, 'previous_runner': old,
        'updated_runner': updated, 'previous_state': previous_state, 'archived_partial_analysis': archived}
    write_json(root/'repair.json', record)
    manifest.setdefault('repair_records', []).append(fingerprint(root/'repair.json'))
    write_json(root/'manifest.json', manifest)
    write_json(root/'state.json', {'status': 'repaired', 'positive_reductions': 0, 'finished': [],
        'repair': str(root/'repair.json')})
    print(f'repaired frozen runner and archived {len(archived)} partial analysis directories', flush=True)


def required_bins(trials: list, radius: np.ndarray) -> tuple[list, set, np.ndarray]:
    """Select complete one-pixel annuli bracketing every search-pixel radius."""
    bins, positions = set(), set()
    for trial in trials:
        for dx, dy in SEARCH_OFFSETS:
            x, y = trial['row']+dx, trial['column']+dy
            positions.add((x, y))
            lower = math.floor(float(radius[y, x])-.5)
            bins.update((lower, lower+1))
    selected = np.zeros(radius.shape, bool)
    for lower in bins:
        selected |= (radius > lower) & (radius <= lower+1)
    return sorted(bins), positions, selected


def filter_pixel(science: np.ndarray, position: tuple, template: np.ndarray,
                 forbidden: np.ndarray) -> tuple[dict, dict]:
    """Measure identity and all three fixed rectangular-PSD amplitudes at one pixel."""
    x, y = position
    data = science[y-5:y+6, x-5:x+6].ravel()
    energy = float(template @ template)
    amplitudes = {'identity': float(template @ data/energy)}
    diagnostics = {'identity': {'valid': True}}
    rings = full.radial.geometry(science, position, forbidden)
    matrices = {offset: full.radial.extract(science, ring) for offset, ring in rings.items()}
    for width in WIDTHS:
        samples = np.vstack([matrix for offset, matrix in matrices.items() if abs(offset) <= width])
        model = full.psd.fit_psd(samples, 'rectangular', .3)
        name = RECTANGULAR[width]
        diagnostics[name] = {'valid': model is not None, 'samples': len(samples)}
        if model is None:
            continue
        result = full.radial.filter_stamp(data, template, model, np.ones(121))
        amplitudes[name] = result['amplitude']
        diagnostics[name].update({key: result[key] for key in ('amplitude', 'sigma', 'score')})
    return amplitudes, diagnostics


def score_trial(trial: dict, maps: np.ndarray, amplitudes: np.ndarray,
                profiles: dict, radius: np.ndarray) -> dict:
    """Record the fixed five-pixel maximum and center amplitude for each method."""
    result = {}
    for index, name in enumerate(METHODS):
        values = [float(maps[index, trial['column']+dy, trial['row']+dx]) for dx, dy in SEARCH_OFFSETS]
        valid = all(np.isfinite(values))
        one = {'valid': valid, 'search_score': max(values) if valid else None,
               'pixels': values, 'center_amplitude': float(amplitudes[index, trial['column'], trial['row']])
                    if np.isfinite(amplitudes[index, trial['column'], trial['row']]) else None}
        if name == 'identity' and valid:
            r = float(radius[trial['column'], trial['row']])
            profile = profiles[name]
            one['center_noise_sigma'] = float(np.interp(r, [p['radius'] for p in profile],
                [p['stddev'] if p['stddev'] is not None else np.nan for p in profile]))
            count = 2*math.pi*r/3.6-1
            one['small_sample_correction'] = 1/math.sqrt(1+1/count)
        result[name] = one
    return result


def analyze_image(root: Path, source: Path, analysis_name: str, trial: dict,
                  reference_name: str) -> dict:
    """Build clean amplitude maps, apply production annular SNR, and score trials."""
    directory = root/'analysis'/analysis_name
    complete = directory/'complete.json'
    if complete.exists():
        full.verify(full.read(complete)['products'])
        return full.read(directory/'measurements.json')
    directory.mkdir(parents=True, exist_ok=False)
    science, header = fits.getdata(source, header=True)
    science = science.squeeze().astype(float)
    protocol = full.read(root/'protocol.json')
    templates = full.load_templates(root/'payload/response')
    source_mask, clean_centers = planet_masks(science.shape, protocol['known_source_circle'])
    forbidden = source_mask.copy()
    full.mark_footprint(forbidden, trial['row'], trial['column'])
    yy, xx = np.indices(science.shape)
    radius = np.hypot(xx-127.5, yy-127.5).astype('f4')
    bins, positions, selected = required_bins([trial], radius)
    maps = np.full((len(METHODS), *science.shape), np.nan, dtype='f4')
    enabled = active_methods(trial)
    enabled_indices = [METHODS.index(name) for name in enabled]
    references = full.reference_maps(root, source, reference_name)
    maps[METHODS.index('gaussian')][selected & clean_centers] = references['gaussian_raw'][selected & clean_centers]
    details = {}
    for index, ((x, y), template) in enumerate(sorted(templates.items())):
        if not selected[y, x] or not clean_centers[y, x]:
            continue
        data = science[y-5:y+6, x-5:x+6]
        if data.size != 121 or not np.isfinite(data).all():
            continue
        amplitudes, diagnostic = filter_pixel(science, (x, y), template, forbidden)
        for name, amplitude in amplitudes.items():
            maps[METHODS.index(name), y, x] = amplitude
        if (x, y) in positions:
            details[f'{x},{y}'] = diagnostic
        if index % 1000 == 0:
            print(f'{analysis_name}: filtered {index}/{len(templates)} response positions', flush=True)
    header['HCI FILTER LABELS'] = ','.join(METHODS)
    header['HCI RADIAL BINS'] = ','.join(map(str, bins))
    fits.writeto(directory/'amplitudes.fits', maps, header)
    active_header = header.copy()
    active_header['HCI FILTER LABELS'] = ','.join(enabled)
    fits.writeto(directory/'active_amplitudes.fits', maps[enabled_indices], active_header)
    command = [str(root/'software/hciAnalyze'), '--file='+str(directory/'active_amplitudes.fits'), '--lambdaD=3.6',
        '--planet.sep=11.782', '--planet.PA=262.051', '--planet.R=7.3', '--snr.apertureR=60',
        '--snr.minRad=0', '--snr.maxRad=60', '--filter.psfResponse=', '--filter.lpfGaussFW=0',
        '--filter.hpfGaussFW=0', '--noise.model=identity', '--noise.only=false', '--noise.outputDiagnostics=false']
    write_json(directory/'command.json', command)
    environment = full.binary_environment(root, reduction=False)
    environment.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with (directory/'analysis.log').open('w') as log:
        subprocess.run(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
    active_snr, snr_header = fits.getdata(directory/'active_amplitudes_snr.fits', header=True)
    active_snr = active_snr.reshape((len(enabled), *science.shape))
    full.radial.require(active_snr.shape == maps[enabled_indices].shape and
                        snr_header['SNRMEAN'] == 1 and snr_header['SNRSMALL'] == 1,
                        'changed production annular SNR contract')
    snr = np.full(maps.shape, np.nan, dtype='f4')
    snr[enabled_indices] = active_snr
    fits.writeto(directory/'amplitudes_snr.fits', snr, snr_header)
    settings = {'source_x': protocol['known_source_circle'][0], 'source_y': protocol['known_source_circle'][1],
        'source_radius': protocol['known_source_circle'][2], 'lambda_d': 3.6, 'min_radius': 0, 'max_radius': 60}
    profiles, errors = {}, {}
    for index, name in enumerate(METHODS):
        if name not in enabled:
            profiles[name], errors[name] = [], None
            continue
        expected, profiles[name] = annular_oracle(maps[index], settings)
        checked = np.array([[y, x] for x, y in positions if np.isfinite(maps[index, y, x])])
        if len(checked):
            difference = np.abs(snr[index, checked[:, 0], checked[:, 1]]-expected[checked[:, 0], checked[:, 1]])
            errors[name] = float(np.nanmax(difference))
            full.radial.require(np.allclose(snr[index, checked[:, 0], checked[:, 1]],
                                expected[checked[:, 0], checked[:, 1]], rtol=2e-6, atol=2e-6, equal_nan=True),
                                'annular oracle mismatch for '+name)
        else:
            errors[name] = None
    rows = [{'trial': trial, 'models': score_trial(trial, snr, maps, profiles, radius)}]
    write_json(directory/'measurements.json', {'image': reference_name, 'analysis': analysis_name,
        'source': fingerprint(source), 'rows': rows,
        'required_bins': bins, 'active_methods': enabled, 'search_details': details,
        'annular_oracle_max_errors': errors})
    products = [fingerprint(directory/name) for name in ('amplitudes.fits', 'active_amplitudes.fits',
        'active_amplitudes_snr.fits', 'amplitudes_snr.fits', 'measurements.json', 'command.json')]
    write_json(complete, {'products': products})
    return full.read(directory/'measurements.json')


def calibrate(root: Path) -> tuple[dict, list, dict]:
    """Freeze null thresholds and identity-based injection contrasts before positives."""
    record = root/'calibration_complete.json'
    if record.exists():
        full.verify(full.read(record)['products'])
        return full.read(root/'thresholds.json'), full.read(root/'jobs.json'), full.read(root/'baseline.json')
    protocol = full.read(root/'protocol.json')
    trials = [*protocol['calibration_trials'], *protocol['sites']]
    baseline_rows = []
    for index, trial in enumerate(trials):
        analyzed = analyze_image(root, root/'payload/baseline.fits', 'baseline__'+trial['name'], trial, 'baseline')
        baseline_rows.append(analyzed['rows'][0])
        print(f'calibrated baseline search {index+1}/{len(trials)}: {trial["name"]}', flush=True)
    baseline = {'rows': baseline_rows, 'source': fingerprint(root/'payload/baseline.fits'),
        'per_search_holdout': True}
    rows = {row['trial']['name']: row for row in baseline_rows}
    thresholds = {}
    for nominal in RADII:
        thresholds[str(nominal)] = {}
        for method in METHODS:
            pool = protocol['calibration_pools'][str(nominal)][method]
            if pool is None:
                thresholds[str(nominal)][method] = None
                continue
            values = [rows[name]['models'][method] for name in pool['trials']]
            full.radial.require(len(values) == CALIBRATION_SEARCHES and all(value['valid'] for value in values),
                                f'invalid calibration pool for radius {nominal} {method}')
            thresholds[str(nominal)][method] = max(value['search_score'] for value in values)
    jobs = []
    for site in protocol['sites']:
        identity = rows[site['name']]['models']['identity']
        threshold = thresholds[str(site['nominal_radius'])]['identity']
        reference = threshold*identity['center_noise_sigma']/identity['small_sample_correction']
        full.radial.require(np.isfinite(reference) and reference > 0, 'invalid identity-based contrast scale')
        for index, level in enumerate(LEVELS):
            jobs.append({**{key: site[key] for key in ('name', 'row', 'column', 'nominal_radius')},
                'site': site['name'], 'name': site['name']+f'_l{index}', 'brightness_multiplier': level,
                'reference_contrast_scale': reference, 'contrast': level*reference,
                'target_values_used_for_contrast': False})
    write_json(root/'baseline.json', baseline)
    write_json(root/'thresholds.json', thresholds)
    write_json(root/'jobs.json', jobs)
    products = [fingerprint(root/name) for name in ('baseline.json', 'thresholds.json', 'jobs.json')]
    write_json(record, {'all_thresholds_and_contrasts_frozen_before_positive_analysis': True,
                        'products': products})
    return thresholds, jobs, baseline


def summarize(root: Path, thresholds: dict, baseline: dict, measurements: list) -> None:
    """Write per-radius recovery, null exceedances, invalid counts, and a compact plot."""
    protocol = full.read(root/'protocol.json')
    baseline_rows = {row['trial']['name']: row for row in baseline['rows']}
    nulls, groups = [], []
    for site in protocol['sites']:
        models = baseline_rows[site['name']]['models']
        for method, value in models.items():
            threshold = thresholds[str(site['nominal_radius'])][method]
            value['detected'] = value['valid'] and threshold is not None and value['search_score'] > threshold
        nulls.append({'site': site, 'models': models})
    for record in measurements:
        trial = record['trial']
        for method, value in record['models'].items():
            threshold = thresholds[str(trial['nominal_radius'])][method]
            value['detected'] = value['valid'] and threshold is not None and value['search_score'] > threshold
            if value['center_amplitude'] is not None and method != 'gaussian':
                value['raw_contrast_error'] = value['center_amplitude']/trial['contrast']-1
    for nominal in RADII:
        for method in METHODS:
            one = {'radius': nominal, 'method': method,
                'null_exceedances': sum(row['models'][method]['detected'] for row in nulls
                                        if row['site']['nominal_radius'] == nominal),
                'invalid_nulls': sum(not row['models'][method]['valid'] for row in nulls
                                     if row['site']['nominal_radius'] == nominal), 'levels': []}
            for level in LEVELS:
                selected = [row for row in measurements if row['trial']['nominal_radius'] == nominal and
                            row['trial']['brightness_multiplier'] == level]
                values = [row['models'][method] for row in selected]
                errors = [value['raw_contrast_error'] for value in values if 'raw_contrast_error' in value]
                one['levels'].append({'brightness_multiplier': level, 'trials': len(values),
                    'detections': sum(value['detected'] for value in values),
                    'invalid_searches': sum(not value['valid'] for value in values),
                    'median_raw_contrast_error': float(np.median(errors)) if errors else None})
            groups.append(one)
    write_json(root/'results.json', {'groups': groups, 'nulls': nulls, 'measurements': measurements,
        'thresholds': thresholds, 'calibration_pools': protocol['calibration_pools'],
        'caveats': [protocol['dependence'], 'Radius-6 planet-clean sites are necessarily close together.',
                    'Method-specific calibration bands are recorded and can span adjacent radii.']})
    lines = ['# Inner-radius rectangular-PSD recovery', '',
        '| Radius | Method | 0.5× | 0.75× | 1× | Nulls | Invalid nulls |',
        '| ---: | --- | ---: | ---: | ---: | ---: | ---: |']
    for group in groups:
        lines.append(f'| {group["radius"]} | {group["method"]} | '+
            ' | '.join(f'{row["detections"]}/{row["trials"]}' for row in group['levels'])+
            f' | {group["null_exceedances"]}/6 | {group["invalid_nulls"]}/6 |')
    lines += ['', 'Thresholds and injection contrasts were frozen from the baseline before positive reductions. '
              'Invalid searches are nondetections. Counts are correlated within this one residual field.', '']
    (root/'results.md').write_text('\n'.join(lines))
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True, sharey=True, layout='constrained')
    for axis, nominal in zip(axes.flat, RADII):
        for method in METHODS:
            group = next(row for row in groups if row['radius'] == nominal and row['method'] == method)
            axis.plot(LEVELS, [row['detections']/6 for row in group['levels']], 'o-', label=method)
        axis.set(title=f'r = {nominal} px', xlabel='Brightness / identity threshold scale', ylabel='Recovery fraction',
                 ylim=(-.03, 1.03))
    axes[0, 0].legend(fontsize=6)
    fig.suptitle('Inner-radius full P4 injections: pooled rectangular PSD widths')
    fig.savefig(root/'comparison.png', dpi=170)
    plt.close(fig)


def run(root: Path) -> None:
    """Calibrate once, then resume the unattended queue of 108 full reductions."""
    manifest = full.read(root/'manifest.json')
    full.verify([*manifest['frozen_records'], *manifest['input_records'], *manifest.get('repair_records', [])])
    with (root/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if (root/'complete.json').exists():
            full.verify(full.read(root/'complete.json')['products'])
            print('Inner-radius study already complete.', flush=True)
            return
        finished = []
        try:
            write_json(root/'state.json', {'status': 'calibrating', 'pid': os.getpid(), 'finished': finished})
            thresholds, jobs, baseline = calibrate(root)
            calibration_record = full.read(root/'calibration_complete.json')
            measurements = []
            for job in jobs:
                record_path = root/'measurements'/(job['name']+'.json')
                if record_path.exists():
                    record = full.read(record_path)
                    full.verify(record['products'])
                else:
                    write_json(root/'state.json', {'status': 'reducing', 'current_trial': job['name'],
                        'pid': os.getpid(), 'finished': finished, 'total': len(jobs)})
                    source = full.reduce_image(root, job)
                    write_json(root/'state.json', {'status': 'analyzing', 'current_trial': job['name'],
                        'pid': os.getpid(), 'finished': finished, 'total': len(jobs)})
                    site = next(site for site in full.read(root/'protocol.json')['sites'] if site['name'] == job['site'])
                    analyzed = analyze_image(root, source, job['name'], site, job['name'])
                    record = {'trial': job, 'models': analyzed['rows'][0]['models'],
                        'products': [fingerprint(source), fingerprint(root/'analysis'/job['name']/'measurements.json')]}
                    record_path.parent.mkdir(exist_ok=True)
                    write_json(record_path, record)
                measurements.append(record)
                finished.append(job['name'])
                write_json(root/'state.json', {'status': 'running', 'pid': os.getpid(), 'finished': finished,
                    'total': len(jobs)})
                print(f'completed {len(finished)}/{len(jobs)}: {job["name"]}', flush=True)
            summarize(root, thresholds, baseline, measurements)
            full.verify([*manifest['frozen_records'], *manifest['input_records'],
                         *manifest.get('repair_records', []), *calibration_record['products']])
            products = [fingerprint(root/name) for name in ('results.json', 'results.md', 'comparison.png',
                'thresholds.json', 'jobs.json', 'baseline.json')]
            write_json(root/'complete.json', {'positive_reductions': len(finished), 'finished': finished,
                'frozen_inputs_unchanged': True, 'products': products})
            write_json(root/'state.json', {'status': 'complete', 'finished': finished,
                'results': str(root/'results.json')})
        except Exception as error:
            write_json(root/'state.json', {'status': 'failed', 'pid': os.getpid(),
                'finished': finished, 'error': str(error)})
            raise


def main() -> None:
    """Expose score-free audit, immutable setup, and resumable execution actions."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('audit', 'setup', 'repair', 'run'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--parent', type=Path, default=Path('working/roc/p4_psd_full_20260918'))
    parser.add_argument('--cpus', type=int, nargs='+', default=list(range(24)))
    args = parser.parse_args()
    if args.action == 'audit':
        audit(args)
    elif args.action == 'setup':
        setup(args)
    elif args.action == 'repair':
        repair(args)
    else:
        run(args.root.resolve())


if __name__ == '__main__':
    main()
