#!/usr/bin/env python3
"""Prepare a spatially held-out Step-5 pilot and freeze null thresholds.

Trial locations, shared exclusions, the five-pixel search aperture, and noise
settings are recorded before scores are calculated. This small single-field
experiment reports achieved rates and angular-block sensitivity; overlapping
trials do not provide an independent-sample or five-sigma guarantee.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from astropy.io import fits

from analyze_p4_step5_development import samples, source_circle
from run_p4_step5_full_injections import fingerprint, write_json

MODELS = ('identity', 'diagonal', 'pca0', 'pca')
OFFSETS = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))


def prepare(args: argparse.Namespace) -> None:
    """Freeze disjoint trial neighborhoods and a lossless subset of the response field."""
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    science = fits.getdata(args.science).squeeze()
    if science.shape != (256, 256) or fits.getheader(args.science)['COMBINATION METHOD'].strip() != 'mean':
        raise ValueError('this experiment requires the full 256x256 mean-combined baseline')
    trials = []
    for role, angles in (('calibration', (-60, -50, -40, -30)), ('evaluation', (30, 40, 50, 60))):
        for radius in (20, 26, 30, 34, 38, 42, 46, 50):
            for block, angle in enumerate(angles):
                row = round(127.5 + radius * math.cos(math.radians(angle)))
                column = round(127.5 + radius * math.sin(math.radians(angle)))
                trials.append({'name': f'{role}_r{radius}_b{block}', 'role': role,
                    'nominal_radius': radius, 'angle_degrees_from_positive_first_axis': angle,
                    'angular_block': block, 'row': row, 'column': column})
    # sqrt(6^2+5^2) encloses every 11x11 stamp in the radius-one, five-pixel search aperture.
    holdout_radius = math.sqrt(61)
    circles = [[trial['row'], trial['column'], holdout_radius] for trial in trials]
    circles += [[x, y, holdout_radius] for x, y in ((120, 137), (109, 143), (98, 98))]
    all_exclusions = [source_circle(), *circles]
    geometry = []
    for trial in trials:
        support = []
        for dx, dy in OFFSETS:
            position = (trial['row']+dx, trial['column']+dy)
            _, _, counts = samples(science, position, all_exclusions)
            support.append({'row': position[0], 'column': position[1], **counts})
        geometry.append({'trial': trial['name'], 'candidates': support})
    # Verify the development/calibration/evaluation support unions are pairwise disjoint.
    yy, xx = np.indices(science.shape)
    masks = {}
    for role in ('calibration', 'evaluation'):
        mask = np.zeros(science.shape, bool)
        for trial in trials:
            if trial['role'] == role:
                mask |= np.hypot(xx-trial['row'], yy-trial['column']) <= holdout_radius
        masks[role] = mask
    masks['development'] = np.logical_or.reduce([
        np.hypot(xx-x, yy-y) <= holdout_radius for x, y in ((120, 137), (109, 143), (98, 98))])
    for a, b in (('development', 'calibration'), ('development', 'evaluation'), ('calibration', 'evaluation')):
        if (masks[a] & masks[b]).any():
            raise RuntimeError(f'{a}/{b} holdout pixel neighborhoods overlap; revise geometry before scoring')
    field = args.products.resolve() / 'finim_outputs'
    coordinates = fits.getdata(field / 'p4PSF_coordinates.fits')
    requested = {(trial['row']+dx, trial['column']+dy) for trial in trials for dx, dy in OFFSETS}
    selection = np.array([(int(x), int(y)) in requested for x, y in coordinates[:2].T])
    subset = root / 'response'
    subset.mkdir()
    for name in ('manifest', 'coordinates', 'model_0000', 'validity_0000'):
        data, header = fits.getdata(field / f'p4PSF_{name}.fits', header=True)
        if name == 'coordinates' or name.startswith('validity'):
            data = data[:, selection]
        elif name.startswith('model'):
            data = data[selection]
        header['HIERARCH P4 PSF SOURCE COUNT'] = str(int(selection.sum()))
        header['HIERARCH P4 STEP5 SUBSET'] = 'unchanged native response values at preregistered search pixels'
        fits.writeto(subset / f'p4PSF_{name}.fits', data, header)
    (root / 'software').mkdir()
    for source in (args.binary, args.library):
        shutil.copy2(source, root / 'software' / source.name)
    policy = {'schema': 1, 'purpose': 'single-field held-out covariance pilot, not independent high-significance calibration',
        'trials': trials, 'search_offsets_row_column': OFFSETS, 'known_source_circle': source_circle(),
        'extra_exclusions': circles, 'holdout_radius_pixels': holdout_radius,
        'models': list(MODELS), 'pca_modes': 3, 'pca0_modes': 0, 'floor_fraction': 0.1,
        'minimum_training_patches': 8, 'arc_step_pixels': 5, 'guard_pixels': 0,
        'exact_exclusion': True, 'lambdaD_pixels': 3.6, 'target_false_positive_fraction': 0.05,
        'selection': 'fixed initial policy; no rank/floor tuning on any recovered score',
        'threshold_rule': 'calibration order statistic ceil((n+1)*0.95), strict exceedance; unresolvable if rank>n',
        'uncertainty': 'descriptive resampling of four angular blocks; blocks may remain correlated',
        'science': fingerprint(args.science), 'cpu_ids': args.cpus,
        'source_response_inputs': [fingerprint(p) for p in sorted(field.glob('p4PSF_*.fits'))],
        'frozen_response_inputs': [fingerprint(p) for p in sorted(subset.glob('*.fits'))],
        'software': [fingerprint(p) for p in sorted((root / 'software').iterdir())],
        'script': fingerprint(Path(__file__))}
    write_json(root / 'protocol.json', policy)
    write_json(root / 'geometry.json', {'pairwise_holdout_masks_disjoint': True, 'trials': geometry})
    print(root)


def analyze(args: argparse.Namespace) -> None:
    """Run the production filters, audit the holdouts, and freeze separate model thresholds."""
    root = args.output.resolve()
    protocol = json.loads((root / 'protocol.json').read_text())
    for expected in [protocol['science'], *protocol['source_response_inputs'],
                     *protocol['frozen_response_inputs'], *protocol['software']]:
        if fingerprint(Path(expected['path'])) != expected:
            raise RuntimeError('a preregistered input changed: ' + expected['path'])
    if (root / 'null_results.json').exists():
        raise RuntimeError('null results already exist; do not overwrite frozen thresholds')
    pilot = Path(__file__).with_name('run_p4_step5_pilot.py')
    command = [sys.executable, str(pilot), '--science', protocol['science']['path'],
        '--manifest', str(root / 'response/p4PSF_manifest.fits'),
        '--binary', str(root / 'software/hciAnalyze'), '--library', str(root / 'software/libhcireduce.so'),
        '--output', str(root / 'analysis'), '--lambda-d', '3.6', '--source-radius', '7.3',
        '--exact-exclusion', '--zero-mode-control', '--conditional-only', '--cpus', *map(str, protocol['cpu_ids'])]
    for circle in protocol['extra_exclusions']:
        command += ['--exclusion', *map(str, circle)]
    write_json(root / 'analysis_command.json', command)
    with (root / 'analysis.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    maps = {model: {role: fits.getdata(root / f'analysis/{model}/science_{role}.fits').squeeze()
                   for role in ('psf_score', 'psf_amplitude', 'psf_sigma', 'noise_status', 'noise_samples',
                                'noise_attempted', 'noise_excluded', 'noise_incomplete')}
            for model in MODELS}
    geometry = json.loads((root / 'geometry.json').read_text())
    records = []
    for trial, expected in zip(protocol['trials'], geometry['trials']):
        assert trial['name'] == expected['trial']
        locations = [(trial['row']+dx, trial['column']+dy) for dx, dy in OFFSETS]
        # Cross-check every candidate's actual production training against the independent heldout-stencil audit.
        for candidate, (x, y) in zip(expected['candidates'], locations):
            status = maps['pca']['noise_status'][y, x]
            if status in (0, 1, 2):
                for key, role in (('accepted', 'samples'), ('attempted', 'attempted'),
                                  ('excluded', 'excluded'), ('incomplete', 'incomplete')):
                    if candidate[key] != maps['pca']['noise_' + role][y, x]:
                        raise RuntimeError('production/independent heldout geometry mismatch: ' + trial['name'])
        eligible = {model: all(maps[model]['noise_status'][y, x] == 0 for x, y in locations) for model in MODELS}
        sx, sy, sr = protocol['known_source_circle']
        source_overlap = math.hypot(trial['row']-sx, trial['column']-sy) <= sr + protocol['holdout_radius_pixels']
        record = {**trial, 'eligible_by_model': eligible, 'source_overlap': source_overlap,
                  'common_eligible': all(eligible.values()) and not source_overlap, 'models': {}}
        for model in MODELS:
            if not eligible[model]:
                record['models'][model] = None
                continue
            scores = [float(maps[model]['psf_score'][y, x]) for x, y in locations]
            peak = int(np.argmax(scores))
            x, y = trial['row'], trial['column']
            record['models'][model] = {'search_score': scores[peak], 'peak_row_column': locations[peak],
                'center_amplitude': float(maps[model]['psf_amplitude'][y, x]),
                'center_conditional_sigma': float(maps[model]['psf_sigma'][y, x]),
                'center_training_samples': int(maps[model]['noise_samples'][y, x])}
        records.append(record)
    summaries = {}
    for model in MODELS:
        calibration = [r['models'][model]['search_score'] for r in records
                       if r['role'] == 'calibration' and r['common_eligible']]
        rank = math.ceil((len(calibration)+1) * (1-protocol['target_false_positive_fraction']))
        threshold = sorted(calibration)[rank-1] if rank <= len(calibration) else None
        evaluation = [r for r in records if r['role'] == 'evaluation' and r['common_eligible']]
        blocks = []
        for block in range(4):
            selected = [r for r in evaluation if r['angular_block'] == block]
            blocks.append({'angular_block': block, 'trials': len(selected),
                'exceedances': sum(r['models'][model]['search_score'] > threshold for r in selected)
                if threshold is not None else None})
        if threshold is not None:
            rng = np.random.default_rng(51018)
            block_ids = rng.integers(0, 4, size=(10000, 4))
            denominators = np.array([b['trials'] for b in blocks])[block_ids].sum(axis=1)
            fractions = np.array([b['exceedances'] for b in blocks])[block_ids].sum(axis=1) / denominators
            bootstrap = np.percentile(fractions, [2.5, 97.5]).tolist()
        else:
            bootstrap = None
        summaries[model] = {'calibration_trials': len(calibration), 'threshold_rank': rank,
            'threshold': threshold, 'calibration_exceedances': sum(s > threshold for s in calibration)
            if threshold is not None else None,
            'evaluation_trials': len(evaluation), 'evaluation_exceedances': sum(b['exceedances'] for b in blocks)
            if threshold is not None else None, 'angular_blocks': blocks,
            'block_resampling_2p5_97p5_percentiles': bootstrap,
            'uncertainty_caveat': 'four spatial blocks are few and may be correlated; zero observed events cannot bound the tail'}
    write_json(root / 'null_results.json', {'protocol': fingerprint(root / 'protocol.json'),
        'production_geometry_audit_passed': True, 'models': summaries, 'trials': records,
        'response_subset_values': 'copied without modification from the full analytic field',
        'independent_trials_assumed': False, 'detection_gain_established': False})
    print(json.dumps(summaries, indent=2))


def main() -> None:
    """Separate score-free protocol preparation from threshold/evaluation analysis."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'analyze'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--science', type=Path)
    parser.add_argument('--products', type=Path)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--library', type=Path)
    parser.add_argument('--cpus', type=int, nargs='+', default=[12, 13])
    args = parser.parse_args()
    if args.action == 'prepare' and any(getattr(args, name) is None for name in ('science', 'products', 'binary', 'library')):
        parser.error('prepare requires --science, --products, --binary, and --library')
    (prepare if args.action == 'prepare' else analyze)(args)


if __name__ == '__main__':
    main()
