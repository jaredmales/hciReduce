#!/usr/bin/env python3
"""Summarize the full AF Lep analytic/reference comparison and known-source fit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from astropy.io import fits

from compare_p4_refit_difference import (optimizer_summary_path, position_difference, read_exact_fit,
                                         read_response_fit, resource_usage)
from fit_p4_matched_response import read_coordinates, selected_mode_index
from run_p4_response_convergence import fingerprint


def metrics(actual: np.ndarray, reference: np.ndarray) -> dict:
    """Measure scale, shape and total error on finite common support."""
    valid = np.isfinite(actual) & np.isfinite(reference)
    a, b = actual[valid].astype(float), reference[valid].astype(float)
    scale = float(a @ b / (a @ a))
    return {'pixels': int(valid.sum()), 'cosine': float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b))),
            'reference_projection_on_template': scale,
            'relative_error': float(np.linalg.norm(a-b) / np.linalg.norm(b)),
            'best_scaled_error': float(np.linalg.norm(scale*a-b) / np.linalg.norm(b))}


def main() -> None:
    """Reproduce the known-source fit and compare current products with frozen accepted references."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', type=Path, required=True)
    parser.add_argument('--stage', default='analytic_sparse_batch32')
    parser.add_argument('--refit-reference', type=Path, required=True)
    parser.add_argument('--exact-reference', type=Path, required=True)
    args = parser.parse_args()
    root, refit, exact_root = args.experiment.resolve(), args.refit_reference.resolve(), args.exact_reference.resolve()
    stage = root / args.stage
    completion = json.loads((stage / 'complete.json').read_text())
    if not completion['science_bitwise_equal']:
        raise RuntimeError('analytic run failed science invariance')
    control_stage = root / 'baseline_same_build'
    control = json.loads((control_stage / 'complete.json').read_text())
    control_manifest = json.loads((control_stage / 'manifest.json').read_text())
    product_manifest = json.loads((stage / 'manifest.json').read_text())
    if not control['science_bitwise_equal']:
        raise RuntimeError('same-build response-free control changed science')
    if (control_manifest['binary']['sha256'] != product_manifest['binary']['sha256'] or
        {k: v['sha256'] for k, v in control_manifest['libraries'].items()} !=
        {k: v['sha256'] for k, v in product_manifest['libraries'].items()}):
        raise RuntimeError('science-invariance control used different software')
    products = stage / 'finim_outputs'
    previous = refit / 'refit_response' / 'finim_outputs'
    output = root / 'scientific_comparison'
    output.mkdir(exist_ok=False)
    fit_dir = output / 'analytic_fit'
    subprocess.run([sys.executable, str(Path(__file__).with_name('fit_p4_matched_response.py')),
        str(stage / 'finim.fits'), str(products / 'p4PSF_manifest.fits'), str(fit_dir),
        '--initial-separation', '11.782', '--initial-pa', '262.051'], check=True)
    coordinates = read_coordinates(products / 'p4PSF_coordinates.fits')
    if not np.array_equal(coordinates, read_coordinates(previous / 'p4PSF_coordinates.fits')):
        raise RuntimeError('accepted and analytic response coordinates differ')
    model = fits.getdata(products / 'p4PSF_model_0000.fits')
    reference = fits.getdata(previous / 'p4PSF_model_0000.fits')
    if not np.array_equal(np.isfinite(model), np.isfinite(reference)):
        raise RuntimeError('accepted and analytic finite support differ')
    header = fits.getheader(products / 'p4PSF_manifest.fits')
    diagnostics, diagnostic_header = fits.getdata(products / 'p4PSF_measurement_diagnostics.fits', header=True)
    diagnostics = diagnostics.T
    if int(header['P4 PSF MEASUREMENT COUNT']) != 232 or int(header['P4 PSF SAMPLE EXCLUDED COUNT']) != 221:
        raise RuntimeError('full accepted sampling/avoidance geometry was not reproduced')
    exact = read_exact_fit(optimizer_summary_path(exact_root))
    summaries = {}
    for name, path in [('analytic', fit_dir / 'summary.json'), ('accepted_refit', refit / 'refit_fit/summary.json'),
                       ('replayed_refit', root / 'replayed_refit_fit/summary.json')]:
        fit = read_response_fit(path)
        summaries[name] = {**fit, 'contrast_bias_fraction_vs_exact': fit['contrast']/exact['contrast']-1,
                           'position_error_pixels_vs_exact': position_difference(fit, exact)}
    # Preserve the accepted stamp, source-support and subpixel differences in the comparison's interpretation.
    accepted_fit = json.loads((refit / 'refit_fit/summary.json').read_text())['fit']
    row, column = int(accepted_fit['integer_peak_row']), int(accepted_fit['integer_peak_column'])
    source = np.flatnonzero((coordinates[:, 0] == row) & (coordinates[:, 1] == column))
    if source.size != 1:
        raise RuntimeError('known-source template is not unique')
    original_path = exact_root / 'sparse_response/finim.fits'
    original, original_header = fits.getdata(original_path, header=True)
    original = original[selected_mode_index(original_header, .15, original_path)]
    removed = fits.getdata(exact_root / 'signal_free_oracle/finim.fits')
    empirical = ((original - removed)/exact['contrast'])[column-5:column+6, row-5:row+6]
    norms = np.sqrt(np.nansum(reference.astype(float)**2, axis=(1, 2)))
    errors = np.sqrt(np.nansum((model.astype(float)-reference)**2, axis=(1, 2)))
    positive = norms > 0
    ratios = errors[positive] / norms[positive]
    radial_groups = []
    radii = np.hypot(coordinates[:, 0]-127.5, coordinates[:, 1]-127.5)
    for lower, upper in ((0, 3.5), (3.5, 6), (6, 12), (12, 24), (24, 42), (42, 60)):
        selected = positive & (radii >= lower) & (radii < upper)
        relative = errors[selected]/norms[selected]
        radial_groups.append({'minimum_radius': lower, 'maximum_radius': upper, 'count': int(selected.sum()),
                              'median': float(np.median(relative)), 'p95': float(np.quantile(relative, .95)),
                              'maximum': float(relative.max())})
    positive_indices = np.flatnonzero(positive)
    largest_errors = []
    for index in positive_indices[np.argsort(ratios)[-8:][::-1]]:
        largest_errors.append({'row': int(coordinates[index, 0]), 'column': int(coordinates[index, 1]),
            'radius': float(radii[index]), 'finite_pixels': int(np.isfinite(model[index]).sum()),
            'relative_error': float(errors[index]/norms[index]), 'reference_norm': float(norms[index])})
    psf = fits.getdata(product_manifest['template']['path']).astype(float)
    source_support = []
    for size in (12, 14, 16):
        first_column, first_row = (psf.shape[0]-size)//2, (psf.shape[1]-size)//2
        crop = psf[first_column:first_column+size, first_row:first_row+size]
        source_support.append({'crop_pixels': size,
            'fraction_of_full_squared_energy': float(np.sum(crop**2)/np.sum(psf**2)),
            'fraction_of_full_signed_sum': float(crop.sum()/psf.sum())})
    summary = {'schema': 1, 'science_bitwise_equal': True, 'same_build_control': control, 'frames': 621, 'mode_fraction': .15,
        'combination': {'science': header['P4 SCIENCE COMBINATION'], 'response': header['P4 PSF COMBINATION'],
                        'science_sigma_threshold': header['P4 SCIENCE SIGMA THRESHOLD']},
        'sampling': {k: header[k] for k in ['P4 PSF MEASUREMENT COUNT', 'P4 PSF SAMPLE EXCLUDED COUNT',
                    'P4 PSF ANALYTIC FACTOR COUNT', 'P4 PSF ANALYTIC BATCH SIZE']},
        'outcome_counts': dict(zip(str(diagnostic_header['P4 PSF DIAGNOSTIC COLUMNS']).split(',')[4:],
                                  diagnostics[:, 4:].sum(axis=0).astype(int).tolist())),
        'all_templates': metrics(model, reference),
        'per_template_relative_error': {'count': int(positive.sum()), 'median': float(np.median(ratios)),
                                       'p95': float(np.quantile(ratios, .95)), 'maximum': float(ratios.max())},
        'template_error_radial_groups': radial_groups, 'largest_template_errors': largest_errors,
        'exact_optimizer': exact, 'fits': summaries, 'source_support': source_support,
        'known_source_stamp': {'row': row, 'column': column, 'reference': 'finite negative-planet removal with sigmaMean',
            'interpretation': ('Total mismatch includes finite amplitude, subpixel registration, clipping, and source '
                               'support: the archived full-image removal shifts the full 256-pixel PSF, while '
                               'the 11-pixel sparse response crops it to 12 pixels. The accepted optimizer used '
                               'a 13-pixel local window with a 14-pixel source crop. These metrics do not isolate '
                               'analytic-derivative error; the new injection study matches source support.'),
            'analytic': metrics(model[source[0]], empirical), 'accepted_refit': metrics(reference[source[0]], empirical)},
        'resources': {'baseline': resource_usage(control_stage / 'resource_usage.txt'),
                      'analytic': resource_usage(stage / 'resource_usage.txt'),
                      'archived_refit': resource_usage(refit / 'refit_response/resource_usage.txt')},
        'provenance': [fingerprint(p) for p in [products / 'p4PSF_manifest.fits', previous / 'p4PSF_manifest.fits',
                        optimizer_summary_path(exact_root), original_path, exact_root / 'signal_free_oracle/finim.fits']]}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
