#!/usr/bin/env python3
"""Separate background subtraction, PSD weighting, and SNR normalization on saved injections.

Retain the full ROC study's sites, masks, contrasts, and five-pixel searches.
Compute complete radial bins needed by each search, then use production annular
SNR on eight amplitude maps. Freeze all baseline thresholds before positives.
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
from scipy.linalg import cho_solve

import run_p4_step5_roc_full as full
from measure_p4_step5_planet import annular_oracle
from run_p4_step5_full_injections import fingerprint, write_json

METHODS = ('gaussian', 'identity', 'hann_mean_only', 'hann_psd_no_mean', 'hann_psd_mean',
           'rect_mean_only', 'rect_psd_no_mean', 'rect_psd_mean')
FAMILIES = (('hann', 'hann', .1, 0, 'psd_hann_b0_m0.1'),
            ('rect', 'rectangular', .3, 5, 'psd_rectangular_b5_m0.3'))
LABELS = ('Gaussian 3.6', 'Identity', 'Same-ring mean only', 'Hann PSD, no mean', 'Hann PSD + mean',
          'Pooled mean only', 'Rect. PSD, no mean', 'Rect. PSD + mean')
_CONTEXT = None


def narrow_samples(science: np.ndarray, position: tuple, forbidden: np.ndarray) -> dict:
    """Extract the original native interpolation stencils for only the three needed rings."""
    xq, yq = position
    cy, cx = (np.array(science.shape)-1)/2
    radius, angle = math.hypot(xq-cx, yq-cy), math.atan2(yq-cy, xq-cx)
    dy, dx = np.mgrid[-5:6, -5:6]
    matrices = {}
    for offset in (-5, 0, 5):
        radial = radius+offset
        if radial <= 0:
            continue
        count = math.ceil(2*math.pi*radial/5)
        theta = np.array([2*math.pi*i/count for i in range(count)])
        centers = np.array([[cx+radial*math.cos(t), cy+radial*math.sin(t)] for t in theta])
        cosine = np.array([math.cos(t-angle) for t in theta])[:, None]
        sine = np.array([math.sin(t-angle) for t in theta])[:, None]
        xx = centers[:, 0, None]+cosine*dx.ravel()-sine*dy.ravel()
        yy = centers[:, 1, None]+sine*dx.ravel()+cosine*dy.ravel()
        x0, y0 = np.floor(xx).astype(int), np.floor(yy).astype(int)
        fx, fy = xx-x0, yy-y0
        indices, weights = [], []
        rejected = np.zeros(count, bool)
        for oy in (0, 1):
            for ox in (0, 1):
                w = (fx if ox else 1-fx)*(fy if oy else 1-fy)
                used = w != 0
                px, py = x0+ox, y0+oy
                outside = (px < 0) | (px >= science.shape[1]) | (py < 0) | (py >= science.shape[0])
                sx, sy = np.clip(px, 0, science.shape[1]-1), np.clip(py, 0, science.shape[0]-1)
                bad = forbidden[sy, sx] & ~outside
                bad |= (np.abs(px-xq) <= 5) & (np.abs(py-yq) <= 5)
                rejected |= np.any(used & (bad | outside | ~np.isfinite(science[sy, sx])), axis=1)
                indices.append(sy*science.shape[1]+sx)
                weights.append(w)
        indices = np.stack(indices, axis=-1)[~rejected]
        weights = np.stack(weights, axis=-1)[~rejected]
        full.radial.require(not np.any(forbidden.ravel()[indices[weights != 0]]), 'excluded pixel enters training')
        values = np.where(weights != 0, science.ravel()[indices], 0)
        matrices[offset] = np.sum(weights*values, axis=2)
    return matrices


def filter_pixel(science: np.ndarray, position: tuple, template: np.ndarray, forbidden: np.ndarray) -> tuple:
    """Hold covariance fixed within each mean-on/off pair and use its exact mean in the identity control."""
    x, y = position
    data = science[y-5:y+6, x-5:x+6].ravel()
    matrices = narrow_samples(science, position, forbidden)
    energy = float(template @ template)
    amplitudes = {'identity': float(template @ data/energy)}
    diagnostics = {}
    for prefix, window, mixing, width, old_name in FAMILIES:
        samples = matrices[0] if width == 0 else np.vstack(list(matrices.values()))
        model = full.psd.fit_psd(samples, window, mixing)
        diagnostics[prefix] = {'valid': model is not None, 'samples': len(samples)}
        if model is None:
            continue
        inverse_template = cho_solve(model['factorization'], template, check_finite=False)
        energy_psd = float(template @ inverse_template)
        weight = inverse_template/energy_psd
        sigma = 1/math.sqrt(energy_psd)
        mean_only = float(template @ (data-model['mean'])/energy)
        off, on = float(weight @ data), float(weight @ (data-model['mean']))
        full.radial.require(np.isclose(weight @ template, 1, rtol=1e-12) and
                            np.isclose(off-on, weight @ model['mean'], rtol=1e-9, atol=1e-13), 'mean-toggle identity failed')
        amplitudes.update({prefix+'_mean_only': mean_only, prefix+'_psd_no_mean': off, prefix+'_psd_mean': on})
        diagnostics[prefix].update(amplitude=on, sigma=sigma, score=on/sigma, mean_projection=float(weight @ model['mean']))
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
    """Freeze the eight-method comparison and source records before any new measurements."""
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
    scripts = [fingerprint(p) for p in sorted(Path(__file__).parent.glob('*.py'))]
    records = list({r['path']: r for r in [*records, *native, *scripts]}.values())
    full.verify(records)
    protocol = {'schema': 1, 'purpose': 'development ablation on inspected images: mean subtraction, PSD weights, common annular SNR',
        'study': str(study), 'methods': list(METHODS), 'families': FAMILIES, 'new_reductions': 0,
        'sites': parent['sites'], 'positive_images': 90, 'brightnesses': list(full.LEVELS),
        'search': 'unchanged native center plus four axial one-pixel neighbors; all five required; invalid search is nondetection',
        'training': 'unchanged sitewise parent holdout and complete finite 11x11 raw patches; offsets -5,0,+5; step5; minimum8',
        'mean_toggle': 'same centered training covariance and weights for both PSD variants; toggle only subtraction of the fitted mean from candidate data',
        'mean_only': 'identity weights with the exact same mean patch and support as the corresponding PSD family',
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
    _CONTEXT = (path, study, parent, radius, full.load_templates(study/'payload/response'), calibration, nulls, positives)


def measure(task: dict) -> str:
    """Produce one site's baseline or positive amplitude/SNR cubes and audited search records."""
    root, study, parent, radius, templates, old_calibration, old_nulls, old_positives = _CONTEXT
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
    previous = fits.getdata(study/'reductions/baseline/finim.fits').squeeze()
    full.radial.require(np.array_equal(np.isfinite(science), np.isfinite(previous)), 'changed science support')
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
            comparisons.append((pixel['amplitudes']['rect_mean_only'], old['isotropic_b5']['pixels'][pixel_index]['amplitude']))
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
        for reference in (prefix+'_mean_only', prefix+'_psd_no_mean', 'identity', 'gaussian'):
            tested = prefix+'_psd_mean'
            pairs[tested+'/'+reference] = [{'level': level,
                'test_only': [r['trial']['name'] for r in rows if r['job']['brightness_multiplier'] == level and
                              r['models'][tested]['detected'] and not r['models'][reference]['detected']],
                'reference_only': [r['trial']['name'] for r in rows if r['job']['brightness_multiplier'] == level and
                                   r['models'][reference]['detected'] and not r['models'][tested]['detected']]}
                for level in full.LEVELS]
    verification = {'archived_scalar_checks': sum(r['archived_scalar_checks'] for r in [*baseline, *positives]),
        'fitted_psd_models': sum(r['fitted_psd_models'] for r in [*baseline, *positives]),
        'maximum_archived_scalar_difference': max(r['maximum_archived_scalar_difference'] for r in [*baseline, *positives]),
        'maximum_annular_oracle_difference': max(v for r in [*baseline, *positives] for v in r['oracle_errors'].values() if v is not None)}
    write_json(root/'results.json', {'models': summary, 'paired_decisions': pairs, 'nulls': nulls, 'injections': rows,
        'thresholds': thresholds, 'verification': verification, 'baseline_records': baseline, 'new_reductions': 0,
        'original_study_summary': {'models': parent['models'], 'groups': parent['groups']}, 'reference_decisions_unchanged': True})
    lines = ['# Mean subtraction and common annular SNR', '', '| Method | 0.5× /30 | 0.75× /30 | 1× /30 | Nulls /30 |',
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
    fig.suptitle('Saved 90 injections: mean subtraction, PSD weighting, common annular SNR')
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
