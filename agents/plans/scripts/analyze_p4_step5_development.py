#!/usr/bin/env python3
"""Audit response support and full-image source leakage for Step-5 development.

This independent NumPy sampler checks the exact native-pixel exclusion contract.
Difference images are used only to diagnose source support and training leakage;
they are not measurements of completeness, false-positive rate, or contrast bias.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

from run_p4_step5_full_injections import fingerprint, write_json


def samples(science: np.ndarray, position: tuple[int, int], exclusions: list,
            half: int = 5, arc_step: float = 5) -> tuple[np.ndarray, list, dict]:
    """Independently gather complete stencils and audit every nonzero input pixel."""
    row, column = position
    cy, cx = (np.array(science.shape) - 1) / 2
    radius = math.hypot(row - cx, column - cy)
    angle = math.atan2(column - cy, row - cx)
    count = math.ceil(2 * math.pi * radius / arc_step)
    dy, dx = np.mgrid[-half:half+1, -half:half+1]
    accepted, centers = [], []
    counts = {'attempted': count, 'excluded': 0, 'incomplete': 0}
    for index in range(count):
        theta = 2 * math.pi * index / count
        xx, yy = cx + radius * math.cos(theta), cy + radius * math.sin(theta)
        cosine, sine = math.cos(theta - angle), math.sin(theta - angle)
        x, y = xx + cosine * dx - sine * dy, yy + sine * dx + cosine * dy
        x0, y0 = np.floor(x).astype(int), np.floor(y).astype(int)
        fx, fy = x - x0, y - y0
        stencil = []
        excluded, incomplete = False, False
        for oy in (0, 1):
            for ox in (0, 1):
                weight = (fx if ox else 1-fx) * (fy if oy else 1-fy)
                used = weight != 0
                px, py = x0 + ox, y0 + oy
                forbidden = (abs(px-row) <= half) & (abs(py-column) <= half)
                for er, ec, size in exclusions:
                    forbidden |= np.hypot(px-er, py-ec) <= size
                excluded |= bool((used & forbidden).any())
                outside = (px < 0) | (px >= science.shape[1]) | (py < 0) | (py >= science.shape[0])
                incomplete |= bool((used & outside).any())
                stencil.append((px, py, weight, used))
        if excluded:
            counts['excluded'] += 1
            continue
        if incomplete:
            counts['incomplete'] += 1
            continue
        values = np.zeros(dx.shape)
        for px, py, weight, used in stencil:
            pixels = science[py[used], px[used]]
            incomplete |= bool((~np.isfinite(pixels)).any())
            values[used] += weight[used] * pixels
        if incomplete:
            counts['incomplete'] += 1
            continue
        accepted.append(values.ravel())  # FITS/NumPy order agrees with Eigen column-major native coordinates.
        centers.append([xx, yy])
    counts['accepted'] = len(accepted)
    return np.array(accepted).reshape((-1, (2*half+1)**2)), centers, counts


def source_circle() -> list[float]:
    """Return the known-source coordinates, with a full-response enclosing mask."""
    angle = math.radians(262.051)
    return [127.5 - 11.782 * math.sin(angle), 127.5 + 11.782 * math.cos(angle), 7.3]


def response_support(products: Path) -> dict:
    """Measure contained energy relative to the stored stamp, without assuming absent wings are zero."""
    field = products / 'finim_outputs'
    coordinates = fits.getdata(field / 'p4PSF_coordinates.fits').T
    response = fits.getdata(field / 'p4PSF_model_0000.fits').astype(float)
    available = fits.getdata(field / 'p4PSF_validity_0000.fits').ravel() == 1
    dy, dx = np.mgrid[-5:6, -5:6]
    rho = np.hypot(dx, dy)
    radius = np.hypot(coordinates[:, 0]-127.5, coordinates[:, 1]-127.5)
    result = {}
    for lower, upper in ((10, 20), (20, 30), (30, 40), (40, 50)):
        selection = available & (radius >= lower) & (radius < upper) & np.isfinite(response).all(axis=(1, 2))
        stamps = response[selection]
        energy = np.sum(stamps**2, axis=(1, 2))
        negative = np.minimum(0, stamps)**2
        negative_energy = negative.sum(axis=(1, 2))
        fractions = []
        for cutoff in (3.6, 4, 5, 6, math.sqrt(50)):
            retained = rho <= cutoff
            fractions.append({'radius_pixels': cutoff,
                'energy_fraction_p05_median_p95': np.percentile(
                    np.sum(stamps[:, retained]**2, axis=1) / energy, [5, 50, 95]).tolist(),
                'negative_energy_fraction_p05_median_p95': np.percentile(
                    negative[:, retained].sum(axis=1)[negative_energy > 0] / negative_energy[negative_energy > 0],
                    [5, 50, 95]).tolist()})
        result[f'{lower}-{upper}'] = {'templates': int(selection.sum()), 'fractions': fractions}
    return result


def audit_pilot(pilot: Path) -> dict:
    """Cross-check training geometry against production maps at deterministic representative candidates."""
    science = fits.getdata(pilot / 'pca/science.fits').squeeze()
    status = fits.getdata(pilot / 'pca/science_noise_status.fits').squeeze()
    header = fits.getheader(pilot / 'pca/science_noise_samples.fits')
    if header['HCIA NOISE EXCLUSION'].strip() != 'NONZERO_STENCIL_PIXEL_CENTERS':
        raise ValueError('independent audit requires exact stencil exclusion')
    exclusions = [list(map(float, circle.split(','))) for circle in header['HCIA NOISE EXCLUSIONS'].split(';')]
    maps = {key: fits.getdata(pilot / f'pca/science_noise_{key}.fits').squeeze()
            for key in ('samples', 'attempted', 'excluded', 'incomplete')}
    candidates = [(int(x), int(y)) for y, x in zip(*np.where((status == 0) | (status == 1)))][::137]
    candidates.append((139, 126))
    records = []
    for candidate in candidates:
        _, _, counts = samples(science, candidate, exclusions)
        x, y = candidate
        for key, role in (('accepted', 'samples'), ('attempted', 'attempted'), ('excluded', 'excluded'),
                          ('incomplete', 'incomplete')):
            if counts[key] != int(maps[role][y, x]):
                raise RuntimeError(f'production geometry disagrees at {candidate}: {key}')
        records.append({'row': x, 'column': y, **counts})
    return {'production_checked_candidates': len(records), 'matches': True, 'records': records}


def leakage(root: Path) -> list:
    """Quantify finite-source changes to accepted training data in full-image reductions."""
    if not (root / 'complete.json').exists():
        raise RuntimeError('full development batch must finish before leakage analysis')
    manifest = json.loads((root / 'manifest.json').read_text())
    baseline = fits.getdata(root / 'baseline/finim.fits').squeeze().astype(float)
    y, x = np.indices(baseline.shape)
    rows = []
    for trial in manifest['trials']:
        if trial['contrast'] == 0:
            continue
        positive = fits.getdata(root / trial['name'] / 'finim.fits').squeeze().astype(float)
        row, column = trial['row'], trial['column']
        distance = np.hypot(x-row, y-column)
        difference = positive - baseline
        energy = np.nansum(difference**2)
        record = {'trial': trial, 'source_difference_energy': float(energy), 'masks': []}
        for cutoff in (3.6, 5, 6, 7.3, 10):
            exclusions = [source_circle(), [row, column, cutoff]]
            before, centers, counts = samples(baseline, (row, column), exclusions)
            after, after_centers, after_counts = samples(positive, (row, column), exclusions)
            if centers != after_centers or counts != after_counts:
                raise RuntimeError('injection changes complete training support; requires separate support analysis')
            b = before - before.mean(axis=0)
            a = after - after.mean(axis=0)
            covariance_before = b.T @ b
            covariance_after = a.T @ a
            record['masks'].append({'radius_pixels': cutoff, **counts,
                'difference_energy_fraction_inside_circle': float(np.nansum(difference[distance <= cutoff]**2)/energy),
                'training_change_over_centered_noise_norm': float(np.linalg.norm(after-before)/np.linalg.norm(b)),
                'covariance_relative_frobenius_change': float(
                    np.linalg.norm(covariance_after-covariance_before)/np.linalg.norm(covariance_before))})
        rows.append(record)
    return rows


def main() -> None:
    """Write reproducible geometry/support audits, optionally including finished full-image trials."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--products', type=Path, required=True)
    parser.add_argument('--pilot', type=Path, required=True)
    parser.add_argument('--development', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('--output must be new')
    records = [fingerprint(args.pilot / 'manifest.json'), fingerprint(args.pilot / 'summary.json')]
    result = {'purpose': 'development geometry/source support; no detection-gain claim',
              'provenance': records, 'lambdaD_pixels': 3.6,
              'support_denominator': 'energy inside stored 11x11 stamp, not total processed source energy',
              'response_support': response_support(args.products), 'production_audit': audit_pilot(args.pilot)}
    if args.development:
        result['full_image_leakage'] = leakage(args.development)
        result['provenance'].append(fingerprint(args.development / 'manifest.json'))
    write_json(args.output, result)
    print(args.output)


if __name__ == '__main__':
    main()
