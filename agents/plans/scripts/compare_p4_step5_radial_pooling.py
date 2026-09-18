#!/usr/bin/env python3
"""Compare radial pooling and variance normalization on saved Step-5 development data.

The existing calibration/evaluation labels retain the original spatial split, but
all these inspected data now serve development. The same-radius raw control is
checked against production before any positive-image comparison. No P4 reduction
is rerun, and no production filtering policy is changed by this NumPy prototype.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.linalg import cho_factor, cho_solve, eigh

from analyze_p4_step5_development import samples
from run_p4_step5_full_injections import fingerprint, write_json

WIDTHS = (0, 5, 10, 20)
OFFSETS = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))
POLICIES = [(f'{"normalized" if normalized else "raw"}_b{width}', width, normalized)
            for normalized in (False, True) for width in WIDTHS]


def read_json(path: Path):
    """Read a saved experiment record."""
    return json.loads(path.read_text())


def require(condition: bool, message: str) -> None:
    """Stop on an inconsistent source, geometry, or numerical control."""
    if not condition:
        raise RuntimeError(message)


def variance_profile(science: np.ndarray, forbidden: np.ndarray) -> tuple[np.ndarray, list]:
    """Estimate the previously specified native radial scale outside all held-out pixels."""
    yy, xx = np.indices(science.shape)
    radius = np.hypot(xx - (science.shape[1] - 1) / 2, yy - (science.shape[0] - 1) / 2)
    allowed = np.isfinite(science) & ~forbidden & (radius < 60)
    profile = []
    for lower in np.arange(0, 60, 3.6):
        upper = min(float(lower + 3.6), 60)
        pixels = science[allowed & (radius >= lower) & (radius < upper)]
        require(len(pixels) >= 20, 'unsupported radial variance bin')
        variance = float(np.var(pixels, ddof=1))
        require(np.isfinite(variance) and variance > 0, 'invalid radial variance')
        profile.append({'lower': float(lower), 'upper': upper, 'pixels': len(pixels), 'variance': variance})
    scale = np.sqrt(np.exp(np.interp(radius, [(p['lower'] + p['upper']) / 2 for p in profile],
                                    np.log([p['variance'] for p in profile]))))
    scale[radius >= 60] = np.nan
    return scale, profile


def geometry(science: np.ndarray, position: tuple[int, int], forbidden: np.ndarray) -> dict:
    """Cache exact nonzero interpolation stencils across the fixed radial grid."""
    xq, yq = position
    cy, cx = (np.array(science.shape) - 1) / 2
    radius, angle = math.hypot(xq - cx, yq - cy), math.atan2(yq - cy, xq - cx)
    dy, dx = np.mgrid[-5:6, -5:6]
    rings = {}
    for offset in range(-20, 21, 5):
        radial = radius + offset
        if radial <= 0:
            continue
        count = math.ceil(2 * math.pi * radial / 5)
        theta = np.array([2 * math.pi * i / count for i in range(count)])
        centers = np.array([[cx + radial * math.cos(t), cy + radial * math.sin(t)] for t in theta])
        cosine = np.array([math.cos(t - angle) for t in theta])[:, None]
        sine = np.array([math.sin(t - angle) for t in theta])[:, None]
        xx = centers[:, 0, None] + cosine * dx.ravel() - sine * dy.ravel()
        yy = centers[:, 1, None] + sine * dx.ravel() + cosine * dy.ravel()
        x0, y0 = np.floor(xx).astype(int), np.floor(yy).astype(int)
        fx, fy = xx - x0, yy - y0
        indices, weights, input_y = [], [], []
        excluded = np.zeros(count, bool)
        incomplete = np.zeros(count, bool)
        for oy in (0, 1):
            for ox in (0, 1):
                w = (fx if ox else 1 - fx) * (fy if oy else 1 - fy)
                used = w != 0
                px, py = x0 + ox, y0 + oy
                outside = (px < 0) | (px >= science.shape[1]) | (py < 0) | (py >= science.shape[0])
                sx, sy = np.clip(px, 0, science.shape[1] - 1), np.clip(py, 0, science.shape[0] - 1)
                bad = forbidden[sy, sx] & ~outside
                bad |= (np.abs(px - xq) <= 5) & (np.abs(py - yq) <= 5)
                excluded |= np.any(used & bad, axis=1)
                incomplete |= np.any(used & (outside | ~np.isfinite(science[sy, sx])), axis=1)
                indices.append(sy * science.shape[1] + sx)
                weights.append(w)
                input_y.append(np.where(used, py, np.nan))
        accepted = ~(excluded | incomplete)
        yy_used = np.stack(input_y, axis=-1)[accepted]
        # Reject straddlers for split-covariance diagnostics: the two sets read disjoint y half-planes.
        halves = np.where(np.nanmax(yy_used, axis=(1, 2)) < cy, 0,
                          np.where(np.nanmin(yy_used, axis=(1, 2)) > cy, 1, -1))
        rings[offset] = {'indices': np.stack(indices, axis=-1)[accepted],
            'weights': np.stack(weights, axis=-1)[accepted], 'centers': centers[accepted], 'halves': halves,
            'counts': {'attempted': count, 'excluded': int(excluded.sum()),
                       'incomplete': int((incomplete & ~excluded).sum()), 'accepted': int(accepted.sum())}}
    return rings


def extract(science: np.ndarray, ring: dict) -> np.ndarray:
    """Read a cached stencil without allowing unused NaN neighbors to contaminate a sample."""
    values = np.where(ring['weights'] != 0, science.ravel()[ring['indices']], 0)
    return np.sum(ring['weights'] * values, axis=2)


def fit(samples: np.ndarray, floor_fraction: float = 0.1) -> dict | None:
    """Fit at most three modes with a configurable median-variance floor; default to the original policy."""
    require(np.isfinite(floor_fraction) and floor_fraction > 0, 'floor fraction must be positive')
    if len(samples) < 8:
        return None
    mean = samples.mean(axis=0)
    centered = samples - mean
    covariance = centered.T @ centered / (len(samples) - 1)
    floor = floor_fraction * float(np.median(np.diag(covariance)))
    if not np.isfinite(floor) or floor <= 0:
        return None
    eigenvalues, eigenvectors = eigh(covariance, subset_by_index=(118, 120), check_finite=False)
    retained = eigenvalues > floor
    modes = eigenvectors[:, retained]
    model = floor * np.eye(121) + (modes * (eigenvalues[retained] - floor)) @ modes.T
    return {'mean': mean, 'covariance': model, 'factorization': cho_factor(model, lower=True, check_finite=False),
            'floor': floor, 'modes': int(retained.sum()), 'samples': len(samples),
            'sample_trace': float(np.trace(covariance)), 'top_three_fraction': float(eigenvalues.sum() / np.trace(covariance))}


def filter_stamp(data: np.ndarray, template: np.ndarray, model: dict, scale: np.ndarray) -> dict:
    """Scale data and response consistently, retaining amplitude in physical contrast units."""
    t = template / scale
    z = data / scale - model['mean']
    weight = cho_solve(model['factorization'], t, check_finite=False)
    energy = float(t @ weight)
    require(energy > 0 and np.isfinite(energy), 'invalid filter normalization')
    amplitude = float(weight @ z / energy)
    sigma = 1 / math.sqrt(energy)
    return {'amplitude': amplitude, 'sigma': sigma, 'score': amplitude / sigma,
            'physical_weight': weight / (energy * scale)}


def analyze_position(science: np.ndarray, scale_map: np.ndarray, position: tuple[int, int],
                     template: np.ndarray, rings: dict, split: bool = False,
                     floor_fraction: float = 0.1) -> tuple[dict, dict]:
    """Evaluate the fixed sampling/normalization grid at one candidate and variance-floor fraction."""
    x, y = position
    data = science[y - 5:y + 6, x - 5:x + 6].ravel()
    scale = scale_map[y - 5:y + 6, x - 5:x + 6].ravel()
    require(np.isfinite(data).all() and np.isfinite(template).all() and np.isfinite(scale).all(), 'invalid candidate support')
    matrices = {}
    for normalized in (False, True):
        image = science / scale_map if normalized else science
        matrices[normalized] = {offset: extract(image, ring) for offset, ring in rings.items()}
    rows, details = {}, {}
    for name, width, normalized in POLICIES:
        selected = [offset for offset in rings if abs(offset) <= width]
        samples = np.vstack([matrices[normalized][offset] for offset in selected])
        model = fit(samples, floor_fraction)
        rows[name] = {'valid': model is not None, 'samples': len(samples)}
        if model is None:
            continue
        factors = scale if normalized else np.ones(121)
        measured = filter_stamp(data, template, model, factors)
        physical_covariance = factors[:, None] * model['covariance'] * factors[None, :]
        physical_mean = factors * model['mean']
        rows[name].update({k: measured[k] for k in ('amplitude', 'sigma', 'score')})
        rows[name].update({'floor_in_fit_units': model['floor'], 'modes': model['modes'],
                          'top_three_sample_variance_fraction': model['top_three_fraction'],
                          'modeled_sample_trace_fraction': float(np.trace(model['covariance']) / model['sample_trace'])})
        details[name] = {'model': model, 'physical_covariance': physical_covariance,
                         'physical_weight': measured['physical_weight'], 'samples': samples}
        if normalized:
            # An independent solve in original units must give the same contrast and sigma.
            factorization = cho_factor(physical_covariance, lower=True, check_finite=False)
            weight = cho_solve(factorization, template, check_finite=False)
            energy = float(template @ weight)
            physical_amplitude = float(weight @ (data - physical_mean) / energy)
            require(np.isclose(physical_amplitude, measured['amplitude'], rtol=1e-9, atol=1e-12) and
                    np.isclose(1 / math.sqrt(energy), measured['sigma'], rtol=1e-9), 'normalization changes contrast units')
        if split:
            halves = np.concatenate([rings[offset]['halves'] for offset in selected])
            parts = [fit(samples[halves == half], floor_fraction) for half in (0, 1)]
            diagnostic = {'samples_by_half': [int((halves == half).sum()) for half in (0, 1)],
                          'straddling_patches_removed': int((halves < 0).sum()), 'valid': all(p is not None for p in parts)}
            if diagnostic['valid']:
                weights = [filter_stamp(data, template, p, factors)['physical_weight'] for p in parts]
                covariance = [factors[:, None] * p['covariance'] * factors[None, :] for p in parts]
                diagnostic.update({'weight_cosine': float(weights[0] @ weights[1] / (np.linalg.norm(weights[0]) * np.linalg.norm(weights[1]))),
                    'relative_covariance_difference': float(2 * np.linalg.norm(covariance[0] - covariance[1]) /
                        (np.linalg.norm(covariance[0]) + np.linalg.norm(covariance[1])))})
            rows[name]['split_stability'] = diagnostic
    return rows, details


def trial_measurement(trial: dict, pixels: dict, policy: str) -> dict:
    """Keep search selection separate from exact-center photometry and uncertainty."""
    x, y = trial['row'], trial['column']
    selected = [pixels[(x + dx, y + dy)][policy] for dx, dy in OFFSETS]
    valid = all(p['valid'] for p in selected)
    center = selected[0]
    return {'valid': valid, 'search_score': max(p['score'] for p in selected) if valid else None,
            'minimum_search_samples': min(p['samples'] for p in selected), 'center': center}


def summarize(nulls: list, injections: list, thresholds: dict, original: dict, gaussian: dict) -> dict:
    """Describe the reused development sample without changing its eligibility or thresholds."""
    groups = {}
    common_stability = [n for n in nulls if n['original_common_eligible'] and
        all(n['policies'][policy]['center'].get('split_stability', {}).get('valid', False)
            for policy, _, _ in POLICIES)]
    for policy, _, _ in POLICIES:
        evaluated = [n for n in nulls if n['role'] == 'evaluation' and n['original_common_eligible']]
        centers = [n['policies'][policy]['center'] for n in evaluated]
        positives = [p for p in injections if p['role'] == 'original_evaluation']
        level_groups = []
        for level in (0.5, 1, 2):
            chosen = [p for p in positives if p['trial']['brightness_multiplier'] == level]
            measured = [p['policies'][policy] for p in chosen]
            errors = [m['center']['amplitude'] / p['trial']['contrast'] - 1 for p, m in zip(chosen, measured) if m['center']['valid']]
            level_groups.append({'brightness_multiplier': level, 'trials': len(chosen),
                'detections': sum(m['detected'] for m in measured), 'invalid': sum(not m['valid'] for m in measured),
                'median_raw_contrast_error': float(np.median(errors)) if errors else None,
                'inside_one_conditional_sigma': sum(m['inside_one_conditional_sigma'] for m in measured)})
        stability = [n['policies'][policy]['center']['split_stability'] for n in nulls
                     if n['original_common_eligible'] and n['policies'][policy]['center']['valid']]
        usable = [s for s in stability if s['valid']]
        common = [n['policies'][policy]['center']['split_stability'] for n in common_stability]
        groups[policy] = {'threshold': thresholds[policy], 'evaluation_nulls': len(evaluated),
            'evaluation_exceedances': sum(n['policies'][policy]['valid'] and n['policies'][policy]['search_score'] > thresholds[policy] for n in evaluated),
            'evaluation_invalid': sum(not n['policies'][policy]['valid'] for n in evaluated),
            'null_center_one_sigma': sum(c['valid'] and abs(c['amplitude']) <= c['sigma'] for c in centers),
            'positive_one_sigma': sum(p['policies'][policy]['inside_one_conditional_sigma'] for p in positives),
            'new_radius20_valid_searches': sum(n['nominal_radius'] == 20 and n['policies'][policy]['valid'] for n in nulls),
            'groups': level_groups, 'split_stability_eligible_centers': len(usable),
            'median_split_weight_cosine': float(np.median([s['weight_cosine'] for s in usable])) if usable else None,
            'median_split_relative_covariance_difference': float(np.median([s['relative_covariance_difference'] for s in usable])) if usable else None,
            'common_split_stability_centers': len(common),
            'common_median_split_weight_cosine': float(np.median([s['weight_cosine'] for s in common])) if common else None,
            'common_median_split_relative_covariance_difference': float(np.median([s['relative_covariance_difference'] for s in common])) if common else None,
            'median_modeled_sample_trace_fraction': float(np.median([c['modeled_sample_trace_fraction'] for c in centers if c['valid']]))}
    return {'policies': groups,
            'common_split_stability_trials': [n['name'] for n in common_stability],
            'references': {name: gaussian['models'][name] for name in ('identity', 'gaussian_raw', 'gaussian_snr')},
            'reference_original_groups': original['groups']}


def plots(summary: dict, output: Path) -> None:
    """Show sample-reuse recovery, uncertainty coverage, and descriptive split stability."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), layout='constrained')
    colors = ('#b66d16', '#436eae')
    for normalized, color in zip((False, True), colors):
        prefix = 'normalized' if normalized else 'raw'
        rows = [summary['policies'][f'{prefix}_b{width}'] for width in WIDTHS]
        label = 'radial normalization' if normalized else 'raw pixels'
        axes[0, 0].plot(WIDTHS, [r['groups'][1]['detections'] for r in rows], 'o-', color=color, label=label)
        axes[0, 1].plot(WIDTHS, [r['evaluation_exceedances'] for r in rows], 'o-', color=color, label=label)
        axes[1, 0].plot(WIDTHS, [r['positive_one_sigma'] for r in rows], 'o-', color=color, label=label)
        axes[1, 1].plot(WIDTHS, [r['common_median_split_weight_cosine'] for r in rows], 'o-', color=color, label=label)
    axes[0, 0].axhline(6, color='0.5', ls='--', label='identity / Gaussian-SNR reference')
    axes[0, 1].axhline(2, color='0.5', ls='--', label='identity / Gaussian-SNR reference')
    axes[0, 0].set(title='Middle-brightness recovery', ylabel='Recovered injections out of 6', ylim=(0, 6.4))
    axes[0, 1].set(title='Original evaluation nulls', ylabel='Threshold exceedances out of 28')
    axes[1, 0].set(title='Raw photometric conditional coverage', ylabel='Within ±1 conditional sigma out of 18', ylim=(0, 18))
    axes[1, 1].set(title=f'Split-half weight stability: {len(summary["common_split_stability_trials"])} common sites',
                   ylabel='Median cosine between two fitted weight vectors')
    for axis in axes.flat:
        axis.set(xlabel='Radial band half-width [pixels]', xticks=WIDTHS)
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)
    fig.suptitle('Step 5 development: radial pooling × variance normalization\nRank 3 and floor fraction 0.1 fixed; previously inspected images reused', fontsize=12)
    fig.savefig(output / 'comparison.png', dpi=170)
    plt.close(fig)


def main() -> None:
    """Verify the original production control, freeze thresholds, then reuse saved positive images."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluation', type=Path, required=True)
    parser.add_argument('--gaussian', type=Path, required=True)
    parser.add_argument('--development', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    evaluation = args.evaluation.resolve()
    protocol = read_json(evaluation / 'protocol.json')
    original_nulls = read_json(evaluation / 'null_results.json')
    original = read_json(evaluation / 'results.json')
    gaussian = read_json(args.gaussian / 'results.json')
    baseline_path = Path(protocol['science']['path'])
    science = fits.getdata(baseline_path).squeeze().astype(float)
    # Match the actual float-resolved source exclusion in the original production maps.
    holdout = Path(protocol['software'][0]['path']).parent.parent
    production_header = fits.getheader(holdout / 'analysis/pca/science_noise_samples.fits')
    exclusions = [list(map(float, circle.split(','))) for circle in production_header['HCIA NOISE EXCLUSIONS'].split(';')]
    yy, xx = np.indices(science.shape)
    forbidden = np.zeros(science.shape, bool)
    for x, y, radius in exclusions:
        forbidden |= np.hypot(xx - x, yy - y) <= radius
    baseline_scale, baseline_profile = variance_profile(science, forbidden)
    field = Path(protocol['source_response_inputs'][0]['path']).parent
    coordinates = fits.getdata(field / 'p4PSF_coordinates.fits').T
    responses = fits.getdata(field / 'p4PSF_model_0000.fits').astype(float)
    validity = fits.getdata(field / 'p4PSF_validity_0000.fits').ravel()
    templates = {(int(c[0]), int(c[1])): t.ravel() for c, t, v in zip(coordinates, responses, validity)
                 if v == 1 and np.isfinite(t).all()}
    jobs = [{'trial': j, 'role': 'original_evaluation', 'image': evaluation / 'reductions' / j['name'] / 'finim.fits'}
            for j in read_json(evaluation / 'jobs.json')]
    development_manifest = read_json(args.development / 'manifest.json')
    jobs += [{'trial': j, 'role': 'original_development', 'image': args.development.resolve() / j['name'] / 'finim.fits'}
             for j in development_manifest['trials'] if j['contrast'] > 0]
    centers = {(t['row'], t['column']) for t in protocol['trials']} | {(j['trial']['row'], j['trial']['column']) for j in jobs}
    positions = sorted({(x + dx, y + dy) for x, y in centers for dx, dy in OFFSETS})
    positive_positions = {(j['trial']['row'] + dx, j['trial']['column'] + dy) for j in jobs for dx, dy in OFFSETS}
    inputs = [protocol['science'], *protocol['source_response_inputs'], fingerprint(evaluation / 'protocol.json'),
              fingerprint(evaluation / 'null_results.json'), fingerprint(evaluation / 'results.json'),
              fingerprint(args.gaussian / 'results.json'), fingerprint(args.development / 'manifest.json')]
    for job in jobs:
        completed = read_json(job['image'].parent / 'complete.json')
        inputs += completed['products']
    for record in inputs:
        require(fingerprint(Path(record['path'])) == record, 'changed original input')
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = {'purpose': 'development-only radial pooling/normalization comparison on already inspected images',
        'policies': [{'name': n, 'half_width': w, 'normalized': u} for n, w, u in POLICIES],
        'maximum_modes': 3, 'floor_fraction': 0.1, 'minimum_patches': 8, 'radial_and_angular_step': 5,
        'profile': '3.6-pixel variance bins after union exclusions; log interpolation; native normalization before patch interpolation; re-estimated for every image',
        'comparison_support': 'original 28+28 common-eligible null searches; newly available radius-20 searches reported separately',
        'split_stability': 'native stencil entirely above/below image-center horizontal line; common full-image variance profile; descriptive, not independent cross-validation',
        'exclusions': exclusions, 'inputs': inputs,
        'scripts': [fingerprint(Path(__file__)), fingerprint(Path(__file__).with_name('analyze_p4_step5_development.py'))]}
    write_json(args.output / 'manifest.json', manifest)
    baseline_pixels, cached_geometry, baseline_details = {}, {}, {}
    production_maps = {role: fits.getdata(holdout / f'analysis/pca/science_{role}.fits').squeeze()
                       for role in ('psf_amplitude', 'psf_sigma', 'psf_score', 'noise_samples', 'noise_status')}
    original_positions = {(t['row'] + dx, t['column'] + dy) for t in protocol['trials'] for dx, dy in OFFSETS}
    comparisons, control_errors = 0, []
    for index, position in enumerate(positions):
        rings = geometry(science, position, forbidden)
        rows, details = analyze_position(science, baseline_scale, position, templates[position], rings, split=position in centers)
        baseline_pixels[position] = rows
        if position in positive_positions:
            cached_geometry[position] = rings
        if position in centers:
            baseline_details[position] = details
        if position in original_positions:
            x, y = position
            control = rows['raw_b0']
            require(control['samples'] == int(production_maps['noise_samples'][y, x]), 'cached geometry differs from production')
            require(control['valid'] == (production_maps['noise_status'][y, x] == 0), 'validity differs from production')
            if control['valid']:
                for key, role in (('amplitude', 'psf_amplitude'), ('sigma', 'psf_sigma'), ('score', 'psf_score')):
                    observed = float(production_maps[role][y, x])
                    require(np.isclose(control[key], observed, rtol=2e-5, atol=1e-10), 'same-radius production numerical mismatch')
                    control_errors.append(abs(control[key] - observed) / max(abs(observed), 1e-12))
            comparisons += 1
        if index % 50 == 0:
            print(f'baseline candidate {index + 1}/{len(positions)}', flush=True)
    # Cross-check the new cached stencils against the prior independent sampler at all development rings.
    stencil_checks = 0
    for job in jobs:
        if job['role'] != 'original_development':
            continue
        q = job['trial']; position = (q['row'], q['column'])
        radius = math.hypot(q['row'] - 127.5, q['column'] - 127.5)
        for offset, ring in cached_geometry[position].items():
            expected, expected_centers, expected_counts = samples(science, position, exclusions, training_radius=radius + offset)
            require(expected_counts == ring['counts'] and np.allclose(np.array(expected_centers).reshape((-1, 2)), ring['centers'], rtol=0, atol=1e-12), 'independent stencil geometry mismatch')
            require(np.allclose(expected, extract(science, ring), rtol=1e-12, atol=1e-12), 'independent stencil sample mismatch')
            stencil_checks += 1
    nulls = [{**{k: t[k] for k in ('name', 'role', 'nominal_radius', 'angular_block', 'row', 'column')},
              'original_common_eligible': t['common_eligible'],
              'policies': {name: trial_measurement(t, baseline_pixels, name) for name, _, _ in POLICIES}}
             for t in original_nulls['trials']]
    thresholds = {}
    for name, _, _ in POLICIES:
        calibration = [t['policies'][name] for t in nulls if t['role'] == 'calibration' and t['original_common_eligible']]
        require(len(calibration) == 28 and all(c['valid'] for c in calibration), 'policy lacks common calibration support')
        thresholds[name] = max(c['search_score'] for c in calibration)
    write_json(args.output / 'thresholds.json', thresholds)
    threshold_record = fingerprint(args.output / 'thresholds.json')
    write_json(args.output / 'nulls.json', nulls)
    positive_rows, profiles, positive_checks = [], {'baseline': baseline_profile}, 0
    for index, job in enumerate(jobs):
        trial = job['trial']; name = trial['name']
        positive = fits.getdata(job['image']).squeeze().astype(float)
        require(np.array_equal(np.isfinite(science), np.isfinite(positive)), 'injection changes finite support')
        scale, profiles[name] = variance_profile(positive, forbidden)
        pixels, central_details = {}, {}
        position = (trial['row'], trial['column'])
        for dx, dy in OFFSETS:
            p = (position[0] + dx, position[1] + dy)
            pixels[p], details = analyze_position(positive, scale, p, templates[p], cached_geometry[p])
            if dx == dy == 0:
                central_details = details
        measured = {}
        for policy, _, _ in POLICIES:
            value = trial_measurement(trial, pixels, policy)
            value['detected'] = value['valid'] and value['search_score'] > thresholds[policy]
            center = value['center']
            value['inside_one_conditional_sigma'] = center['valid'] and abs(center['amplitude'] - trial['contrast']) <= center['sigma']
            base = baseline_pixels[position][policy]
            if center['valid']:
                value['raw_contrast_error'] = center['amplitude'] / trial['contrast'] - 1
            if center['valid'] and base['valid']:
                value['paired_increment_error'] = (center['amplitude'] - base['amplitude']) / trial['contrast'] - 1
                before, after = baseline_details[position][policy], central_details[policy]
                centered = before['samples'] - before['samples'].mean(axis=0)
                value['training_change_over_centered_baseline_norm'] = float(np.linalg.norm(after['samples'] - before['samples']) / np.linalg.norm(centered))
                value['physical_covariance_relative_change'] = float(np.linalg.norm(after['physical_covariance'] - before['physical_covariance']) / np.linalg.norm(before['physical_covariance']))
            measured[policy] = value
        if job['role'] == 'original_evaluation':
            saved = next(r for r in original['measurements'] if r['model'] == 'pca' and r['trial']['name'] == name)
            control = measured['raw_b0']
            for actual, expected in ((control['search_score'], saved['search_score']),
                (control['center']['amplitude'], saved['center_amplitude']),
                (control['center']['sigma'], saved['center_conditional_sigma'])):
                require(np.isclose(actual, expected, rtol=2e-5, atol=1e-10), 'positive production control mismatch')
            require(control['detected'] == saved['detected'] and control['inside_one_conditional_sigma'] == saved['inside_one_conditional_sigma'], 'production decision mismatch')
            positive_checks += 1
        positive_rows.append({'trial': trial, 'role': job['role'], 'policies': measured})
        print(f'positive image {index + 1}/{len(jobs)}: {name}', flush=True)
    summary = summarize(nulls, positive_rows, thresholds, original, gaussian)
    write_json(args.output / 'summary.json', summary)
    write_json(args.output / 'injections.json', positive_rows)
    write_json(args.output / 'variance_profiles.json', profiles)
    verification = {'production_baseline_pixels_checked': comparisons, 'production_positive_trials_checked': positive_checks,
        'maximum_baseline_control_relative_difference': max(control_errors), 'independent_ring_stencil_checks': stencil_checks,
        'normalized_and_original_unit_solves_agree': True, 'original_production_decisions_preserved': True,
        'variance_profiles_and_geometry': 'same native forbidden mask; finite support unchanged in all positive images'}
    write_json(args.output / 'verification.json', verification)
    for record in [*inputs, *manifest['scripts'], threshold_record]:
        require(fingerprint(Path(record['path'])) == record, 'frozen input or threshold changed')
    plots(summary, args.output)
    write_json(args.output / 'complete.json', {'baseline_images': 1, 'positive_images': len(jobs), 'policies': len(POLICIES),
        'new_reductions': 0, 'purpose': 'development reuse, not blind validation', 'verification': verification})
    print(json.dumps(summary['policies'], indent=2))


if __name__ == '__main__':
    main()
