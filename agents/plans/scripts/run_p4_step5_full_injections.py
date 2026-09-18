#!/usr/bin/env python3
"""Run fresh full-image Step-5 reductions with fixed source support and CPU placement.

The default batch measures development-source leakage, not detection completeness.
A later evaluation batch can supply a frozen JSON trial list. Every trial reruns the
whole reduction; no local residual is pasted into an existing science image.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from astropy.io import fits


def fingerprint(path: Path) -> dict:
    """Record the bytes of a file as well as its resolved location."""
    path = path.resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': digest.hexdigest()}


def write_json(path: Path, value: object) -> None:
    """Replace an experiment record atomically with strict JSON."""
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def replacement(command: list[str], overrides: dict) -> list[str]:
    """Replace equals-style options without relying on duplicate-option precedence."""
    return [v for v in command if v.split('=', 1)[0].removeprefix('--') not in overrides] + [
        f'--{key}={value}' for key, value in overrides.items()]


def prepare(args: argparse.Namespace) -> None:
    """Freeze a new full-image batch before examining any of its recovered scores."""
    root = args.output.resolve()
    products = args.products.resolve()
    provenance = json.loads((products / 'manifest.json').read_text())
    trials = json.loads(args.trials.read_text()) if args.trials else [
        {'name': 'baseline', 'row': 120, 'column': 137, 'contrast': 0.0},
        {'name': 'development_r12', 'row': 120, 'column': 137, 'contrast': 0.0011909814823390977},
        {'name': 'development_r24', 'row': 109, 'column': 143, 'contrast': 0.0011909814823390977},
        {'name': 'development_r42', 'row': 98, 'column': 98, 'contrast': 0.0011909814823390977},
    ]
    names = set()
    for trial in trials:
        name = trial['name']
        if (not name or Path(name).name != name or name in ('.', '..') or name in names or
                not all(math.isfinite(trial[key]) for key in ('row', 'column', 'contrast')) or
                not all(0 <= trial[key] < 256 for key in ('row', 'column')) or trial['contrast'] < 0):
            raise ValueError('trials require unique safe names, finite in-image coordinates, and nonnegative contrasts')
        names.add(name)
    if not trials or any(cpu < 0 for cpu in args.cpus) or len(set(args.cpus)) != len(args.cpus):
        raise ValueError('provide trials and distinct nonnegative CPU IDs')
    root.mkdir(parents=True, exist_ok=False)
    software = root / 'software'
    software.mkdir()
    records = []
    for name, expected in ((Path(provenance['binary']['path']).name, provenance['binary']),
                           ('libhcireduce.so', provenance['libraries']['libhcireduce.so'])):
        source = products / 'software' / name
        if fingerprint(source)['sha256'] != expected['sha256']:
            raise RuntimeError('source product software changed: ' + str(source))
        shutil.copy2(source, software / name)
        records.append(fingerprint(software / name))
    # Also freeze mxlib and BLAS so a library rebuild cannot change a later trial.
    for name, expected in provenance['libraries'].items():
        if name == 'libhcireduce.so':
            continue
        source = Path(expected['path'])
        if fingerprint(source)['sha256'] != expected['sha256']:
            raise RuntimeError('source product dependency changed: ' + str(source))
        shutil.copy2(source, software / name)
        records.append(fingerprint(software / name))
    shutil.copy2(Path(__file__), software / Path(__file__).name)
    records.append(fingerprint(software / Path(__file__).name))
    inputs = []
    for expected in provenance['inputs']:
        observed = fingerprint(Path(expected['path']))
        if observed['sha256'] != expected['sha256']:
            raise RuntimeError('input changed: ' + expected['path'])
        inputs.append(observed)
    (root / 'inputs.txt').write_text(''.join(record['path'] + '\n' for record in inputs))
    template = Path(provenance['template']['path'])
    if fingerprint(template)['sha256'] != provenance['template']['sha256']:
        raise RuntimeError('source PSF changed')
    original, header = fits.getdata(template, header=True)
    cropped = np.zeros_like(original)
    first, second = ((size - 12) // 2 for size in original.shape)
    cropped[first:first+12, second:second+12] = original[first:first+12, second:second+12]
    header['HIERARCH P4 STEP5 SOURCE CROP'] = 12
    fits.writeto(root / 'injection_psf.fits', cropped, header)
    command = json.loads((args.experiment / 'common_command.json').read_text())
    config_index = command.index('--config') + 1
    config = Path(command[config_index])
    shutil.copy2(config, root / 'reduction.conf')
    command[config_index] = str(root / 'reduction.conf')
    command[0] = str(software / Path(provenance['binary']['path']).name)
    command = replacement(command, {
        'input.fileList': root / 'inputs.txt', 'p4.localStampSize': 0,
        'p4.writeDiagnostics': 'false', 'p4.memoryFraction': 0,
        'psfResponse.file': '', 'psfResponse.outputModels': 'false', 'psfResponse.filter': 'false',
        'fake.method': 'single', 'fake.fileName': root / 'injection_psf.fits',
        'fake.scaleFileName': '', 'fake.subtractPlanet': 'false', 'combine.method': 'mean',
        'combine.minGoodFract': 0, 'combine.weightFile': '', 'output.outputPSFSub': 'false'})
    write_json(root / 'common_command.json', command)
    write_json(root / 'trials.json', trials)
    records.extend(fingerprint(root / name) for name in (
        'reduction.conf', 'injection_psf.fits', 'inputs.txt', 'common_command.json', 'trials.json'))
    write_json(root / 'manifest.json', {
        'schema': 1, 'purpose': args.purpose, 'input_records': inputs, 'frozen_records': records,
        'source_product_manifest': fingerprint(products / 'manifest.json'),
        'response_manifest': fingerprint(products / 'finim_outputs/p4PSF_manifest.fits'),
        'source_crop_size': 12, 'response_stamp_size': 11, 'source_normalization': 'stored; no renormalization',
        'combination': 'mean', 'cpu_ids': args.cpus, 'openmp_threads': len(args.cpus),
        'blas_threads': 1, 'trials': trials, 'lambdaD_pixels': 3.6,
        'lambdaD_provenance': 'working/analyze.conf; descriptive scale, not a reduction parameter'})
    print(root, flush=True)


def run(args: argparse.Namespace) -> None:
    """Run or resume the frozen batch, stopping rather than hiding any failed trial."""
    root = args.output.resolve()
    manifest = json.loads((root / 'manifest.json').read_text())
    for expected in manifest['input_records'] + manifest['frozen_records']:
        if fingerprint(Path(expected['path'])) != expected:
            raise RuntimeError('frozen file changed: ' + expected['path'])
    command = json.loads((root / 'common_command.json').read_text())
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(manifest['openmp_threads']), OMP_PROC_BIND='true', OMP_PLACES='cores',
                       OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       LD_LIBRARY_PATH=str(root / 'software') + ':' + environment.get('LD_LIBRARY_PATH', ''))
    environment.pop('P4REDUCE_GLOBAL_CONFIG', None)
    finished = []
    write_json(root / 'state.json', {'status': 'running', 'pid': os.getpid(), 'finished': finished})
    try:
        for trial in manifest['trials']:
            directory = root / trial['name']
            complete = directory / 'complete.json'
            if complete.exists():
                previous = json.loads(complete.read_text())
                for expected in previous['products']:
                    if fingerprint(Path(expected['path'])) != expected:
                        raise RuntimeError('completed product changed: ' + expected['path'])
                finished.append(trial['name'])
                continue
            # A failed/incomplete trial is deliberately not overwritten on resume.
            directory.mkdir()
            row, column = trial['row'] - 127.5, trial['column'] - 127.5
            options = replacement(command, {'fake.sep': math.hypot(row, column),
                'fake.PA': math.degrees(-math.atan2(row, column)) % 360,
                'fake.contrast': trial['contrast'], 'output.directory': directory})
            options = ['taskset', '-c', ','.join(map(str, manifest['cpu_ids'])), *options]
            write_json(directory / 'command.json', options)
            write_json(root / 'state.json', {'status': 'running', 'pid': os.getpid(),
                       'current_trial': trial['name'], 'finished': finished})
            print(trial['name'], flush=True)
            start = time.monotonic()
            with (directory / 'run.log').open('w') as log:
                subprocess.run(['/usr/bin/time', '-f',
                    'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M',
                    '-o', str(directory / 'resource_usage.txt'), *options],
                    cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
            data, header = fits.getdata(directory / 'finim.fits', header=True)
            if (data.shape not in ((256, 256), (1, 256, 256)) or header['P4 LOCAL STAMP SIZE'] != 0 or
                    header['COMBINATION METHOD'].strip() != manifest['combination']):
                raise RuntimeError('trial does not have the declared full-image combination contract')
            write_json(complete, {'trial': trial, 'elapsed_seconds': time.monotonic() - start,
                'products': [fingerprint(p) for p in sorted(directory.glob('*.fits'))]})
            finished.append(trial['name'])
        for expected in manifest['input_records'] + manifest['frozen_records']:
            if fingerprint(Path(expected['path'])) != expected:
                raise RuntimeError('frozen file changed during batch: ' + expected['path'])
        write_json(root / 'complete.json', {'finished': finished, 'inputs_unchanged': True,
                                          'detection_calibrated': False})
        write_json(root / 'state.json', {'status': 'complete', 'finished': finished})
    except Exception as error:
        write_json(root / 'state.json', {'status': 'failed', 'finished': finished, 'error': str(error)})
        raise


def main() -> None:
    """Separate batch preparation from execution so prepared work can run unattended."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--products', type=Path)
    parser.add_argument('--experiment', type=Path)
    parser.add_argument('--trials', type=Path)
    parser.add_argument('--cpus', type=int, nargs='+', default=[0, 2, 4, 6, 8, 10])
    parser.add_argument('--purpose', default='development-only full-image source leakage and covariance-mask validation')
    args = parser.parse_args()
    if args.action == 'prepare' and (args.products is None or args.experiment is None):
        parser.error('prepare requires --products and --experiment')
    (prepare if args.action == 'prepare' else run)(args)


if __name__ == '__main__':
    main()
