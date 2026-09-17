#!/usr/bin/env python3
"""Measure finite-amplitude P4 injection recovery with fixed sparse response fields.

Use exact local reductions of the full frame list, at six integer sky positions
and three brightnesses. Baseline subtraction isolates response calibration from
realization-specific residual noise; this is not a completeness/false-alarm test.
Mean and sigmaMean runs expose clipping changes relative to the same mean-combined
response fields. Finite amplitude and sparse spatial interpolation both contribute.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess

import numpy as np
from astropy.io import fits

from fit_p4_matched_response import (bounded_quadratic_maximum, quadratic_coefficients,
                                     quadratic_value, read_coordinates)
from run_p4_response_convergence import fingerprint, read_stamp

POSITIONS = ((120, 137), (137, 120), (109, 143), (146, 112), (98, 98), (157, 157))
AMPLITUDES = (0.0011909814823390977, 0.004763925929356391, 0.019055703717425564)


def write_json(path: Path, value: object) -> None:
    """Persist a human-readable experiment record."""
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def replacement(command: list[str], overrides: dict) -> list[str]:
    """Override recorded equals-style flags without relying on duplicate-option order."""
    return [v for v in command if v.split('=', 1)[0].removeprefix('--') not in overrides] + [
        f'--{k}={v}' for k, v in overrides.items()]


def run(args: argparse.Namespace) -> None:
    """Execute the signed local reductions with fixed software and complete provenance."""
    root = args.output.resolve()
    root.mkdir(exist_ok=False)
    stage = args.products.resolve()
    product_manifest = json.loads((stage / 'manifest.json').read_text())
    binary = stage / 'software' / Path(product_manifest['binary']['path']).name
    template = Path(product_manifest['template']['path'])
    if (fingerprint(binary)['sha256'] != product_manifest['binary']['sha256'] or
        fingerprint(template) != product_manifest['template']):
        raise RuntimeError('product software or template changed')
    library = stage / 'software/libhcireduce.so'
    if fingerprint(library)['sha256'] != product_manifest['libraries']['libhcireduce.so']['sha256']:
        raise RuntimeError('product hciReduce library changed')
    for record in product_manifest['inputs']:
        if fingerprint(Path(record['path'])) != record:
            raise RuntimeError(f"input changed: {record['path']}")
    # The analytic 11-pixel response uses a parity-preserving 12-pixel source crop.
    # A larger local output stamp must not silently broaden the injected source support.
    original_template, template_header = fits.getdata(template, header=True)
    injection_template = np.zeros_like(original_template)
    first_column, first_row = (original_template.shape[0]-12)//2, (original_template.shape[1]-12)//2
    injection_template[first_column:first_column+12, first_row:first_row+12] = \
        original_template[first_column:first_column+12, first_row:first_row+12]
    prepared_template = root / 'injection_psf.fits'
    template_header['HIERARCH P4 STEP3 TEMPLATE CROP'] = 12
    fits.writeto(prepared_template, injection_template, template_header)
    base = json.loads((args.experiment / 'common_command.json').read_text())
    base[0] = str(binary)
    base = replacement(base, {'input.fileList': args.experiment.resolve() / 'inputs.txt', 'p4.localStampSize': 15, 'p4.writeDiagnostics': 'false', 'p4.memoryFraction': 0,
                             'psfResponse.file': '', 'psfResponse.outputModels': 'false', 'psfResponse.filter': 'false',
                             'fake.method': 'single', 'fake.fileName': prepared_template, 'fake.scaleFileName': '',
                             'fake.subtractPlanet': 'false', 'combine.minGoodFract': 0,
                             'combine.weightFile': '', 'output.outputPSFSub': 'false'})
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       LD_LIBRARY_PATH=str(binary.parent) + ':' + environment.get('LD_LIBRARY_PATH', ''))
    environment.pop('P4REDUCE_GLOBAL_CONFIG', None)
    write_json(root / 'manifest.json', {'schema': 1, 'product_provenance': product_manifest,
        'script': fingerprint(Path(__file__).resolve()), 'positions_row_column': POSITIONS,
        'positive_amplitudes': AMPLITUDES, 'combination_methods': ['mean', 'sigmaMean'],
        'threads': args.threads, 'injection_template': fingerprint(prepared_template),
        'source_crop_size': 12, 'local_output_stamp_size': 15, 'response_stamp_size': 11, 'interpretation': 'baseline-subtracted conditional response calibration'})
    for combination in ('mean', 'sigmaMean'):
        for index, (row, column) in enumerate(POSITIONS):
            radius = math.hypot(row - 127.5, column - 127.5)
            angle = math.degrees(-math.atan2(row - 127.5, column - 127.5)) % 360
            for level, amplitude in enumerate((0., *AMPLITUDES, *(-a for a in AMPLITUDES))):
                directory = root / f'{combination}_p{index}_a{level}'
                directory.mkdir()
                command = replacement(base, {'combine.method': combination, 'output.directory': directory,
                    'fake.sep': radius, 'fake.PA': angle, 'fake.contrast': amplitude})
                write_json(directory / 'command.json', command)
                print(directory.name, flush=True)
                with (directory / 'run.log').open('w') as log:
                    subprocess.run(['/usr/bin/time', '-f', 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M',
                        '-o', str(directory / 'resource_usage.txt'), *command], cwd=directory, env=environment,
                        stdout=log, stderr=subprocess.STDOUT, check=True)
                read_stamp(directory, [0.15], 15, (radius, angle, amplitude), 'P4-M32D64')
    write_json(root / 'runs_complete.json', {'runs': 84})


def matched_surface(data: np.ndarray, support: np.ndarray, origin: tuple, position: tuple,
                    templates: dict) -> tuple[np.ndarray, np.ndarray]:
    """Calculate amplitudes and profiled identity-noise likelihood on a 5x5 candidate grid."""
    amplitudes = np.full((5, 5), np.nan)
    likelihood = np.full((5, 5), np.nan)
    row, column = position
    for j, dc in enumerate(range(-2, 3)):
        for i, dr in enumerate(range(-2, 3)):
            template = templates[(row + dr, column + dc)]
            rr, cc = row + dr - origin[0], column + dc - origin[1]
            patch = data[cc-5:cc+6, rr-5:rr+6]
            valid = support[cc-5:cc+6, rr-5:rr+6] & np.isfinite(template) & np.isfinite(patch)
            if valid.sum() != 121:
                raise RuntimeError('injection fit requires common full 11x11 support')
            t, d = template[valid], patch[valid]
            energy = float(t @ t)
            if energy <= 0:
                raise RuntimeError('injection fit has a zero-energy template')
            amplitudes[j, i] = float(t @ d) / energy
            likelihood[j, i] = 0.5 * max(0., float(t @ d)) ** 2 / energy
    return amplitudes, likelihood


def recovered(data: np.ndarray, support: np.ndarray, origin: tuple, position: tuple, templates: dict) -> dict:
    """Locate a bounded quadratic peak within one pixel of the known injected position."""
    amplitudes, likelihood = matched_surface(data, support, origin, position, templates)
    cc, rr = np.unravel_index(int(np.argmax(likelihood[1:4, 1:4])), (3, 3))
    cc, rr = cc + 1, rr + 1
    coefficients = quadratic_coefficients(likelihood, cc, rr)
    try:
        dr, dc, _, boundary = bounded_quadratic_maximum(coefficients,
            (max(-1., 1. - rr), min(1., 3. - rr)), (max(-1., 1. - cc), min(1., 3. - cc)))
        contrast = quadratic_value(quadratic_coefficients(amplitudes, cc, rr), dr, dc)
        error = math.hypot(rr + dr - 2, cc + dc - 2)
        return {'status': 'bounded_peak' if boundary else 'converged', 'contrast': contrast,
                'position_error_pixels': error, 'row_error_pixels': rr+dr-2, 'column_error_pixels': cc+dc-2,
                'fixed_position_contrast': float(amplitudes[2, 2])}
    except RuntimeError as exception:
        return {'status': str(exception), 'contrast': None, 'position_error_pixels': None,
                'row_error_pixels': None, 'column_error_pixels': None,
                'fixed_position_contrast': float(amplitudes[2, 2])}


def analyze(args: argparse.Namespace) -> None:
    """Compare published analytic and paired-refit templates against all signed injection trials."""
    root = args.output.resolve()
    if not args.completed_mean_only and not (root / 'runs_complete.json').is_file():
        raise RuntimeError('injection reductions are incomplete')
    fields = {}
    for name, path in (('analytic', args.products / 'finim_outputs'), ('refit', args.refit_products)):
        coords = read_coordinates(path / 'p4PSF_coordinates.fits')
        model = fits.getdata(path / 'p4PSF_model_0000.fits').astype(float)
        fields[name] = {(int(c[0]), int(c[1])): m for c, m in zip(coords, model)}
    rows = []
    for combination in ('mean',) if args.completed_mean_only else ('mean', 'sigmaMean'):
        for index, position in enumerate(POSITIONS):
            for level in range(7):
                header = fits.getheader(root / f'{combination}_p{index}_a{level}/finim.fits')
                if str(header['COMBINATION METHOD']).strip() != combination:
                    raise RuntimeError('injection combination differs from its experiment label')
                if combination == 'sigmaMean' and header['SIGMA THRESHOLD'] != 5:
                    raise RuntimeError('injection clipping threshold differs from the protocol')
            row, column = position
            radius = math.hypot(row - 127.5, column - 127.5)
            angle = math.degrees(-math.atan2(row - 127.5, column - 127.5)) % 360
            baseline, base_valid, origin = read_stamp(root / f'{combination}_p{index}_a0', [0.15], 15,
                                                     (radius, angle, 0), 'P4-M32D64')
            if origin != (row-7, column-7):
                raise RuntimeError('local stamp is not centered on the requested injection pixel')
            for level, amplitude in enumerate(AMPLITUDES, 1):
                positive, pos_valid, pos_origin = read_stamp(root / f'{combination}_p{index}_a{level}', [0.15], 15,
                                                            (radius, angle, amplitude), 'P4-M32D64')
                negative, neg_valid, neg_origin = read_stamp(root / f'{combination}_p{index}_a{level+3}', [0.15], 15,
                                                            (radius, angle, -amplitude), 'P4-M32D64')
                if pos_origin != origin or neg_origin != origin:
                    raise RuntimeError('signed local lattices differ')
                difference = positive[0] - baseline[0]
                valid = base_valid[0] & pos_valid[0] & neg_valid[0]
                unit = difference[2:-2, 2:-2] / amplitude
                paired = (positive[0] - negative[0])[2:-2, 2:-2] / (2 * amplitude)
                for name, templates in fields.items():
                    template = templates[position]
                    retained = valid[2:-2, 2:-2] & np.isfinite(template)
                    if retained.sum() != 121:
                        raise RuntimeError('injection diagnostics require full common support')
                    t, d, p = template[retained], unit[retained], paired[retained]
                    eta = float(t @ d / (np.linalg.norm(t) * np.linalg.norm(d)))
                    scale = float(t @ d / (t @ t))
                    fit = recovered(difference, valid, origin, position, templates)
                    raw = recovered(positive[0], valid, origin, position, templates)
                    rows.append({'combination': combination, 'position': index, 'row': row, 'column': column,
                        'radius': radius, 'pa_degrees': angle, 'amplitude': amplitude, 'template': name,
                        'cosine': eta, 'fixed_position_bias_fraction': scale - 1,
                        'relative_template_error': float(np.linalg.norm(t-d) / np.linalg.norm(d)),
                        'best_scaled_template_error': float(np.linalg.norm(d-scale*t) / np.linalg.norm(d)),
                        'paired_relative_template_error': float(np.linalg.norm(t-p) / np.linalg.norm(p)),
                        'one_sided_vs_paired_relative_error': float(np.linalg.norm(d-p) / np.linalg.norm(p)),
                        'fit_status': fit['status'], 'contrast_bias_fraction': None if fit['contrast'] is None else fit['contrast']/amplitude-1,
                        'position_error_pixels': fit['position_error_pixels'],
                        'row_error_pixels': fit['row_error_pixels'], 'column_error_pixels': fit['column_error_pixels'],
                        'raw_fit_status': raw['status'],
                        'raw_contrast_bias_fraction': None if raw['contrast'] is None else raw['contrast']/amplitude-1,
                        'raw_position_error_pixels': raw['position_error_pixels'],
                        'raw_row_error_pixels': raw['row_error_pixels'], 'raw_column_error_pixels': raw['column_error_pixels']})
    stem = 'mean_summary' if args.completed_mean_only else 'summary'
    write_json(root / f'{stem}.json', rows)
    with (root / f'{stem}.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({'rows': len(rows), 'summary': str(root / f'{stem}.json')}))


def main() -> None:
    """Dispatch independently reproducible reductions and response-field analysis."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('run', 'analyze'))
    parser.add_argument('--experiment', type=Path, required=True)
    parser.add_argument('--products', type=Path, required=True)
    parser.add_argument('--refit-products', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=20)
    parser.add_argument('--completed-mean-only', action='store_true',
                        help='analyze all 42 completed mean trials into separate mean_summary files')
    args = parser.parse_args()
    if args.action == 'analyze' and args.refit_products is None:
        parser.error('--refit-products is required for analysis')
    if args.action == 'run' and args.completed_mean_only:
        parser.error('--completed-mean-only is only available for analysis')
    (run if args.action == 'run' else analyze)(args)


if __name__ == '__main__':
    main()
