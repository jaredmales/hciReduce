#!/usr/bin/env python3
"""Run a fixed-policy covariance-filter integration pilot, not a completeness study.

All models use the same existing final science image and response manifest.
Outputs go to a new directory, with frozen software and input hashes. No parameter
is selected from the known planet's score. Spatial training is candidate-excluded;
held-out threshold calibration and fresh full-image injections remain separate work.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from astropy.io import fits


def fingerprint(path: Path) -> dict:
    """Identify one exact file without depending on its modification time."""
    path = path.resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def write_json(path: Path, value: object) -> None:
    """Write strict JSON for reproducible experiment records."""
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def run(args: argparse.Namespace) -> None:
    """Freeze software, run three unmodified policies, and summarize training support."""
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    software = root / 'software'
    software.mkdir()
    binary = software / 'hciAnalyze'
    shutil.copy2(args.binary, binary)
    shutil.copy2(args.library, software / 'libhcireduce.so')
    manifest = args.manifest.resolve()
    prefix = manifest.name.removesuffix('manifest.fits')
    inputs = [fingerprint(args.science)] + [fingerprint(p) for p in sorted(manifest.parent.glob(prefix + '*.fits'))]
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       LD_LIBRARY_PATH=str(software) + ':' + environment.get('LD_LIBRARY_PATH', ''))
    environment.pop('HCIANALYZE_GLOBAL_CONFIG', None)
    policy = {'arc_step_pixels': 0, 'extra_guard_pixels': 0, 'minimum_samples': 8,
              'maximum_modes': 3, 'floor_fraction': 0.1, 'known_source_exclusion_pixels': args.source_radius,
              'exact_exclusion': args.exact_exclusion,
              'conditional_maps_only': args.conditional_only,
              'extra_exclusions_row_column_radius': args.exclusion,
              'lambdaD_pixels': args.lambda_d, 'cpu_ids': args.cpus,
              'science_combination': fits.getheader(args.science)['COMBINATION METHOD'].strip(),
              'response_combination': fits.getheader(manifest)['P4 PSF COMBINATION'].strip()}
    repository = Path(__file__).resolve().parents[3]
    sources = software / 'sources'
    sources.mkdir()
    source_records = []
    for relative in ('src/apps/hciAnalyze.cpp', 'src/common/PSFNoiseTraining.hpp',
                     'src/common/PSFNoiseTraining.cpp', 'src/common/PSFNoiseModel.hpp',
                     'src/common/PSFNoiseModel.cpp', 'src/common/P4PSFFilter.hpp', 'src/common/P4PSFFilter.cpp'):
        source = repository / relative
        shutil.copy2(source, sources / source.name)
        source_records.append(fingerprint(source))
    write_json(software / 'source_records.json', source_records)
    (software / 'git_head.txt').write_text(subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=repository, text=True))
    write_json(root / 'manifest.json', {'schema': 1, 'purpose': 'integration and training-coverage pilot; no detection-gain claim',
        'inputs': inputs, 'policy': policy, 'binary': fingerprint(binary),
        'library': fingerprint(software / 'libhcireduce.so'), 'script': fingerprint(Path(__file__))})
    results = {}
    for model in (('identity', 'diagonal', 'pca0', 'pca') if args.zero_mode_control else
                  ('identity', 'diagonal', 'pca')):
        directory = root / model
        directory.mkdir()
        science = directory / 'science.fits'
        shutil.copy2(args.science, science)
        command = ['taskset', '-c', ','.join(map(str, args.cpus)), str(binary),
            f'--file={science}', f'--filter.psfResponse={manifest}',
            f'--noise.model={"pca" if model == "pca0" else model}',
            '--noise.outputDiagnostics=true', '--noise.minimumSamples=8',
            f'--noise.maximumModes={0 if model == "pca0" else 3}',
            '--noise.floorFraction=0.1', '--noise.arcStep=0', '--noise.guardRadius=0',
            f'--noise.exactExclusion={str(args.exact_exclusion).lower()}',
            f'--noise.only={str(args.conditional_only).lower()}',
            f'--lambdaD={args.lambda_d}', '--planet.sep=11.782', '--planet.PA=262.051',
            f'--planet.R={args.source_radius}',
            '--snr.apertureR=2', '--filter.hpfGaussFW=0', '--filter.lpfGaussFW=0']
        if args.exclusion:
            command += [f'--noise.exclude{key}=' + ','.join(str(circle[i]) for circle in args.exclusion)
                        for i, key in enumerate(('Rows', 'Columns', 'Radii'))]
        write_json(directory / 'command.json', command)
        start = time.monotonic()
        with (directory / 'run.log').open('w') as log:
            completed = subprocess.run(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT)
        elapsed = time.monotonic() - start
        record = {'returncode': completed.returncode, 'elapsed_seconds': elapsed}
        write_json(directory / 'process.json', record)
        if completed.returncode:
            raise RuntimeError(f'{model} analysis failed; inspect {directory / "run.log"}')
        status = np.asarray(fits.getdata(directory / 'science_noise_status.fits'))
        counts = np.asarray(fits.getdata(directory / 'science_noise_samples.fits'))
        row, column = np.meshgrid(np.arange(status.shape[-1]), np.arange(status.shape[-2]))
        radius = np.hypot(row - (status.shape[-1] - 1) / 2, column - (status.shape[-2] - 1) / 2)
        radial = []
        for lower, upper in ((0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60)):
            selected = np.isfinite(status) & (radius >= lower) & (radius < upper)
            valid = selected & (status == 0)
            radial.append({'radius_pixels': [lower, upper], 'attempted_candidates': int(selected.sum()),
                'valid_candidates': int(valid.sum()), 'too_few_samples': int((selected & (status == 1)).sum()),
                'zero_variance': int((selected & (status == 2)).sum()),
                'invalid_filter': int((selected & (status == 3)).sum()),
                'median_training_patches': float(np.median(counts[valid])) if valid.any() else None})
        record['radial_training_coverage'] = radial
        record['status_counts'] = {str(code): int((status == code).sum()) for code in range(4)}
        record['no_response'] = int(np.isnan(status).sum())
        record['products'] = [fingerprint(p) for p in sorted(directory.glob('science_*.fits'))]
        results[model] = record
        write_json(root / 'summary.json', {'policy': policy, 'models': results})
    # Confirm the original image and all response inputs remained unchanged.
    for before in inputs:
        if fingerprint(Path(before['path'])) != before:
            raise RuntimeError('pilot input changed: ' + before['path'])
    write_json(root / 'complete.json', {'models': list(results), 'inputs_unchanged': True,
        'held_out_detection_calibrated': False})
    print(json.dumps({'output': str(root), 'models': {k: v['status_counts'] for k, v in results.items()}}, indent=2))


def main() -> None:
    """Parse explicit input, software, and new-output paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--science', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpus', type=int, nargs='+', default=[0, 2])
    parser.add_argument('--lambda-d', type=float, required=True, help='explicit dataset scale in pixels per lambda/D')
    parser.add_argument('--source-radius', type=float, required=True, help='known-source mask radius in pixels')
    parser.add_argument('--exact-exclusion', action='store_true')
    parser.add_argument('--conditional-only', action='store_true', help='skip the separate empirical annular SNR stage')
    parser.add_argument('--zero-mode-control', action='store_true')
    parser.add_argument('--exclusion', type=float, nargs=3, action='append', default=[],
                        metavar=('ROW', 'COLUMN', 'RADIUS'), help='extra held-out circle; may be repeated')
    args = parser.parse_args()
    if not args.manifest.name.endswith('manifest.fits'):
        parser.error('--manifest must name a manifest.fits product')
    if not np.isfinite(args.lambda_d) or args.lambda_d <= 0 or not np.isfinite(args.source_radius) or args.source_radius < 0:
        parser.error('require a positive finite lambda/D scale and a nonnegative finite source radius')
    run(args)


if __name__ == '__main__':
    main()
