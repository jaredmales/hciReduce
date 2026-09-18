#!/usr/bin/env python3
"""Audit radial pooling geometry and variance normalization before filter tuning.

This uses only the saved baseline. It reports sample counts and covariance
diagnostics, not source recovery or independent effective sample counts.
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

from analyze_p4_step5_development import samples
from run_p4_step5_full_injections import fingerprint, write_json


def covariance_summary(matrix: np.ndarray) -> dict:
    """Describe the centered sample covariance without interpreting patch count as independence."""
    if len(matrix) < 2:
        return {'samples': len(matrix), 'mean_pixel_variance': None, 'top_three_variance_fraction': None}
    centered = matrix - matrix.mean(axis=0)
    eigenvalues = np.linalg.eigvalsh(centered.T @ centered / (len(matrix) - 1))
    trace = float(eigenvalues.sum())
    return {'samples': len(matrix), 'mean_pixel_variance': trace / matrix.shape[1],
            'top_three_variance_fraction': float(eigenvalues[-3:].sum() / trace) if trace > 0 else None}


def main() -> None:
    """Compare radius bands with and without an independently withheld native variance profile."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--holdout', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    protocol_path = args.holdout / 'protocol.json'
    protocol = json.loads(protocol_path.read_text())
    science_path = Path(protocol['science']['path'])
    if fingerprint(science_path) != protocol['science']:
        raise RuntimeError('baseline changed since the held-out protocol')
    science = fits.getdata(science_path).squeeze().astype(float)
    yy, xx = np.indices(science.shape)
    cy, cx = (np.array(science.shape) - 1) / 2
    radii = np.hypot(xx - cx, yy - cy)
    exclusions = [protocol['known_source_circle'], *protocol['extra_exclusions']]
    profile_mask = np.isfinite(science) & (radii < 60)
    for x, y, radius in exclusions:
        profile_mask &= np.hypot(xx - x, yy - y) > radius
    profile = []
    for lower in np.arange(0, 60, 3.6):
        upper = min(float(lower + 3.6), 60)
        pixels = science[profile_mask & (radii >= lower) & (radii < upper)]
        if len(pixels) < 20:
            raise RuntimeError('variance-profile bin has too few training pixels; revise policy before filtering')
        variance = float(np.var(pixels, ddof=1))
        if not np.isfinite(variance) or variance <= 0:
            raise RuntimeError('variance profile is not positive and finite')
        profile.append({'lower': float(lower), 'upper': upper, 'center': float((lower + upper) / 2),
                        'pixels': len(pixels), 'variance': variance})
    # Interpolate log variance. Keep endpoint values within the first/last bins;
    # do not extrapolate past the declared 0--60 pixel reduction range.
    scale = np.sqrt(np.exp(np.interp(radii, [p['center'] for p in profile],
                                    np.log([p['variance'] for p in profile]))))
    scale[radii >= 60] = np.nan
    normalized = science / scale
    positions = ((120, 137), (109, 143), (98, 98), (145, 138))
    widths = (0, 5, 10, 20)
    records = []
    for position in positions:
        radius = math.hypot(position[0] - cx, position[1] - cy)
        rings = {}
        for offset in range(-20, 21, 5):
            if radius + offset <= 0:
                continue
            raw, centers, counts = samples(science, position, exclusions, training_radius=radius + offset)
            unit, unit_centers, unit_counts = samples(normalized, position, exclusions, training_radius=radius + offset)
            if centers != unit_centers or counts != unit_counts:
                raise RuntimeError('normalization changes support; compare on common support before proceeding')
            if offset == 0:
                previous, previous_centers, previous_counts = samples(science, position, exclusions)
                if not np.array_equal(previous, raw) or previous_centers != centers or previous_counts != counts:
                    raise RuntimeError('same-radius sampler behavior changed')
            rings[offset] = (raw, unit, counts)
        bands = []
        for width in widths:
            selected = [offset for offset in rings if abs(offset) <= width]
            raw = np.vstack([rings[offset][0] for offset in selected])
            unit = np.vstack([rings[offset][1] for offset in selected])
            bands.append({'half_width_pixels': width, 'raw': covariance_summary(raw),
                          'normalized': covariance_summary(unit),
                          'rings': [{'center_radius': radius + offset, **rings[offset][2]} for offset in selected]})
        records.append({'position_row_column': position, 'radius': radius, 'bands': bands})
    args.output.mkdir(parents=True, exist_ok=False)
    result = {'purpose': 'development geometry and variance-scale audit; no detection-gain test',
        'inputs': [fingerprint(science_path), fingerprint(protocol_path)],
        'scripts': [fingerprint(Path(__file__)), fingerprint(Path(__file__).with_name('analyze_p4_step5_development.py'))],
        'normalization': 'native pixels divided by sigma(radius) before patch interpolation; annular variance estimated after excluding all known-source/development/calibration/evaluation circles',
        'profile_bin_width_pixels': 3.6, 'profile_minimum_pixels': 20,
        'profile_interpolation': 'log variance between bin centers; constant within endpoint half-bins; undefined outside [0,60)',
        'radial_step_pixels': 5, 'angular_arc_step_pixels': 5, 'patch_size_pixels': 11,
        'patch_weighting': 'equal weight per accepted patch; outer rings offer more centers',
        'same_radius_backward_compatibility': True, 'raw_normalized_geometry_identical': True,
        'variance_profile': profile, 'candidates': records,
        'caveats': ['Counts are overlapping patches, not independent samples.',
                    'A positive variance profile does not guarantee transferable covariance shape.',
                    'Profile scale uncertainty and interpolation effects require subsequent calibration.',
                    'The outer variance rise must be investigated before choosing a pooling band.']}
    write_json(args.output / 'audit.json', result)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), layout='constrained')
    centers = [p['center'] for p in profile]
    axes[0, 0].semilogy(centers, np.sqrt([p['variance'] for p in profile]), 'o-')
    axes[0, 0].set(title='Native radial noise scale', xlabel='Stellar radius [pixels]', ylabel='Annular standard deviation [image units]')
    axes[0, 1].plot(centers, [p['pixels'] for p in profile], 'o-')
    axes[0, 1].set(title='Profile training pixels after exclusions', xlabel='Stellar radius [pixels]', ylabel='Finite unmasked pixels per bin')
    for record in records:
        label = f'r={record["radius"]:.1f} px'
        axes[1, 0].plot(widths, [b['raw']['samples'] for b in record['bands']], 'o-', label=label)
        axes[1, 1].plot(widths, [b['normalized']['mean_pixel_variance'] for b in record['bands']], 'o-', label=label)
    axes[1, 0].axhline(8, color='0.5', ls='--', lw=1)
    axes[1, 0].set(title='Training support from additional radii', xlabel='Radial band half-width [pixels]', ylabel='Accepted overlapping patches')
    axes[1, 1].set(title='Sample variance after radial normalization', xlabel='Radial band half-width [pixels]', ylabel='Mean diagonal of centered patch covariance')
    for axis in axes[1]:
        axis.legend(fontsize=8)
    for axis in axes.flat:
        axis.grid(alpha=0.2)
    fig.suptitle('Step 5 development: radial pooling and variance-profile normalization\nBaseline only; no recovery result or independent-sample claim', fontsize=12)
    fig.savefig(args.output / 'diagnostics.png', dpi=170)
    plt.close(fig)
    for expected in result['inputs']:
        if fingerprint(Path(expected['path'])) != expected:
            raise RuntimeError('input changed during audit')
    print(json.dumps([{'radius': r['radius'], 'patches': [b['raw']['samples'] for b in r['bands']]} for r in records], indent=2))


if __name__ == '__main__':
    main()
