#!/usr/bin/env python3
"""Audit inner-radius support for wider pooled rectangular PSD training bands."""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import run_p4_step5_roc_full as full
from run_p4_step5_full_injections import fingerprint, write_json

WIDTHS = (0, 5, 10, 20)
RADII = tuple(range(6, 27))
INJECTION_RADII = (8, 12, 16, 20, 24)


def support(science: np.ndarray, position: tuple, mask: np.ndarray) -> dict:
    """Count accepted patches for every fixed half-width at all five search pixels."""
    x, y = position
    counts = {width: [] for width in WIDTHS}
    for dx, dy in full.radial.OFFSETS:
        rings = full.radial.geometry(science, (x+dx, y+dy), mask)
        for width in WIDTHS:
            counts[width].append(sum(len(ring['halves']) for offset, ring in rings.items()
                                     if abs(offset) <= width))
    return {str(width): {'pixels': values, 'minimum': min(values), 'valid': min(values) >= 8}
            for width, values in counts.items()}


def choose_sites(rows: list, radius: int, prior_centers: list) -> list:
    """Choose six score-free azimuth ranks with complete ±10- and ±20-pixel support."""
    viable = [row for row in rows if 45 <= row['azimuth_degrees'] <= 315 and
              row['support']['10']['valid'] and row['support']['20']['valid'] and
              all(math.dist((row['row'], row['column']), center) >= 4 for center in prior_centers)]
    if len(viable) < 6:
        raise RuntimeError(f'insufficient inner candidates at radius {radius}')
    best = []
    for seed in viable:
        chosen = [seed]
        while True:
            available = [row for row in viable if row not in chosen and
                         all(math.dist((row['row'], row['column']), (prior['row'], prior['column'])) >= 4
                             for prior in chosen)]
            if not available:
                break
            chosen.append(max(available, key=lambda row: min(
                math.dist((row['row'], row['column']), (prior['row'], prior['column'])) for prior in chosen)))
        if len(chosen) > len(best):
            best = chosen
    if len(best) < 6:
        raise RuntimeError(f'cannot separate six candidates at radius {radius}')
    chosen = sorted(best[:6], key=lambda row: row['azimuth_degrees'])
    return [{**row, 'name': f'inner_r{radius}_b{index}', 'azimuth_rank': index}
            for index, row in enumerate(chosen)]


def main() -> None:
    """Map coverage, select a geometry-only inner grid, and verify its PSD fits."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    study, output = args.study.resolve(), args.output.resolve()
    if output.exists():
        raise RuntimeError('output already exists')
    output.mkdir(parents=True)
    protocol = full.read(study/'protocol.json')
    staging = full.read(study/'staging.json')
    previous_protocol = full.read(Path(staging['prior_protocol']['path']))
    prior_centers = [(trial['row'], trial['column']) for trial in previous_protocol['trials']]
    prior_centers += [(120, 137), (109, 143), (98, 98)]
    science_path = study/'reductions/baseline/finim.fits'
    science = fits.getdata(science_path).squeeze().astype(float)
    templates = full.load_templates(study/'payload/response')
    yy, xx = np.indices(science.shape)
    radius_map = np.hypot(xx-127.5, yy-127.5)
    azimuth_map = np.degrees(np.arctan2(yy-127.5, xx-127.5)) % 360
    records, candidates_by_radius = [], {}
    for nominal in RADII:
        candidates = []
        for y, x in np.argwhere(abs(radius_map-nominal) < .5):
            x, y = int(x), int(y)
            if not all((x+dx, y+dy) in templates for dx, dy in full.radial.OFFSETS):
                continue
            if any(not np.isfinite(science[y+dy-5:y+dy+6, x+dx-5:x+dx+6]).all()
                   for dx, dy in full.radial.OFFSETS):
                continue
            site = {'row': x, 'column': y}
            one = {'row': x, 'column': y, 'nominal_radius': nominal,
                   'actual_radius': float(radius_map[y, x]), 'azimuth_degrees': float(azimuth_map[y, x]),
                   'support': support(science, (x, y), full.holdout_mask(science.shape, protocol, site))}
            candidates.append(one)
        candidates.sort(key=lambda row: row['azimuth_degrees'])
        candidates_by_radius[nominal] = candidates
        summary = {'radius': nominal, 'candidate_centers': len(candidates), 'widths': {}}
        for width in WIDTHS:
            valid = [row for row in candidates if row['support'][str(width)]['valid']]
            minima = [row['support'][str(width)]['minimum'] for row in valid]
            summary['widths'][str(width)] = {'valid_searches': len(valid),
                'valid_fraction': len(valid)/len(candidates) if candidates else None,
                'minimum_samples': min(minima) if minima else None,
                'median_minimum_samples': float(np.median(minima)) if minima else None}
        records.append(summary)
        print(f'radius {nominal}: '+', '.join(f'±{width}={summary["widths"][str(width)]["valid_searches"]}'
                                              for width in WIDTHS), flush=True)
    selected = [site for nominal in INJECTION_RADII
                for site in choose_sites(candidates_by_radius[nominal], nominal, prior_centers)]
    fit_checks = []
    for site in selected:
        x, y = site['row'], site['column']
        mask = full.holdout_mask(science.shape, protocol, site)
        rings = full.radial.geometry(science, (x, y), mask)
        matrices = {offset: full.radial.extract(science, ring) for offset, ring in rings.items()}
        one = {'site': site['name'], 'fits': {}}
        for width in (5, 10, 20):
            samples = np.vstack([matrix for offset, matrix in matrices.items() if abs(offset) <= width])
            model = full.psd.fit_psd(samples, 'rectangular', .3)
            one['fits'][str(width)] = {'samples': len(samples), 'valid': model is not None,
                'condition_number': model['condition_number'] if model is not None else None}
            if width >= 10 and model is None:
                raise RuntimeError('selected inner site has invalid fixed PSD fit')
        fit_checks.append(one)
    inputs = [fingerprint(study/name) for name in ('protocol.json', 'geometry.json', 'complete.json')]
    inputs += [fingerprint(science_path), *[fingerprint(path) for path in sorted((study/'payload/response').glob('*.fits'))]]
    write_json(output/'coverage.json', {'schema': 1, 'purpose': 'score-free inner-radius rectangular PSD coverage',
        'widths': list(WIDTHS), 'radii': list(RADII), 'minimum_patches': 8,
        'candidate': 'all five response templates and 11x11 science stamps finite; exact current known-source, 28-calibration, and candidate-footprint exclusions',
        'selection': 'radii 8,12,16,20,24; azimuth 45..315 excludes the known-source sector; ±10 and ±20 complete; six maximin-spaced ranks; minimum four-pixel distance between selected and previously inspected centers',
        'scores_or_recovery_used_for_selection': False, 'records': records, 'selected_sites': selected,
        'selected_center_fit_checks': fit_checks, 'inputs': inputs,
        'script': fingerprint(Path(__file__).resolve()), 'prior_protocol': staging['prior_protocol'], 'new_reductions': 0})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), layout='constrained')
    for width in WIDTHS:
        axes[0].plot(RADII, [row['widths'][str(width)]['valid_fraction'] for row in records], 'o-', label=f'±{width} px')
        axes[1].plot(RADII, [row['widths'][str(width)]['median_minimum_samples'] or np.nan for row in records],
                     'o-', label=f'±{width} px')
    axes[0].set(xlabel='Candidate radius (pixels)', ylabel='Fraction of complete five-pixel searches', ylim=(-.03, 1.03))
    axes[1].set(xlabel='Candidate radius (pixels)', ylabel='Median minimum training patches')
    axes[1].axhline(8, color='0.5', linestyle=':', linewidth=.8)
    for axis in axes:
        axis.legend()
    fig.suptitle('Inner-radius support: rectangular PSD radial pooling')
    fig.savefig(output/'coverage.png', dpi=170)
    plt.close(fig)
    write_json(output/'complete.json', {'all_selected_b10_b20_fits_valid': True,
        'products': [fingerprint(output/name) for name in ('coverage.json', 'coverage.png')], 'new_reductions': 0})


if __name__ == '__main__':
    main()
