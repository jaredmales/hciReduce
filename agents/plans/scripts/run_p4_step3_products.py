#!/usr/bin/env python3
"""Run a timed Step-3 AF Lep product comparison from a recorded baseline command.

Snapshot the executable and hciReduce library so subsequent builds cannot change a
running experiment. Every run retains commands, input and software fingerprints,
resource usage, and the exact same-build science-invariance comparison.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
from astropy.io import fits

from run_p4_response_convergence import fingerprint


def main() -> None:
    """Run one response method and compare its science with the recorded baseline."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--experiment', type=Path, required=True)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--template', type=Path, required=True)
    parser.add_argument('--method', choices=('baseline', 'analytic', 'refitDifference'), default='analytic')
    parser.add_argument('--batch-size', type=int, default=32)
    parser.add_argument('--contrast', type=float, default=0)
    parser.add_argument('--threads', type=int, default=20)
    parser.add_argument('--name', required=True)
    parser.add_argument('--science-reference', type=Path, help='science image for the exact invariance check')
    args = parser.parse_args()
    root = args.experiment.resolve()
    stage = root / args.name
    stage.mkdir(exist_ok=False)
    software = stage / 'software'
    software.mkdir()
    executable = software / args.binary.name
    shutil.copy2(args.binary.resolve(), executable)
    runner = software / Path(__file__).name
    shutil.copy2(Path(__file__).resolve(), runner)
    environment = os.environ.copy()
    adjacent_library = args.binary.resolve().parent / 'libhcireduce.so'
    if adjacent_library.is_file():
        environment['LD_LIBRARY_PATH'] = str(adjacent_library.parent) + ':' + environment.get('LD_LIBRARY_PATH', '')
    libraries = {}
    for line in subprocess.check_output(['ldd', str(args.binary.resolve())], text=True, env=environment).splitlines():
        if '=>' not in line:
            continue
        name, value = line.split('=>', 1)
        path = Path(value.split(' (', 1)[0].strip()).resolve()
        if name.strip().startswith(('libhcireduce.', 'libmxlib.', 'libopenblas.', 'liblapack.')):
            libraries[name.strip()] = fingerprint(path)
            if name.strip().startswith('libhcireduce.'):
                shutil.copy2(path, software / name.strip())
    command = json.loads((root / 'common_command.json').read_text())
    command[0] = str(executable)
    command = [value for value in command if not value.startswith('--input.fileList=')]
    command.append(f'--input.fileList={root / "inputs.txt"}')
    if args.method == 'baseline':
        command += [f'--output.directory={stage}', '--psfResponse.file=', '--psfResponse.outputModels=false',
                    '--psfResponse.filter=false']
    else:
        command += [f'--output.directory={stage}', f'--psfResponse.file={args.template.resolve()}',
                '--psfResponse.stampSize=11', '--psfResponse.outputModels=true', '--psfResponse.filter=false',
                '--psfResponse.outputPrefix=p4PSF_', f'--psfResponse.method={args.method}',
                f'--psfResponse.refitContrast={args.contrast}', f'--psfResponse.analyticBatchSize={args.batch_size}',
                '--psfResponse.radiiPerRegion=2', '--psfResponse.samplesPerRadius=4',
                '--psfResponse.sampleAvoidRadius=5']
    environment.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       LD_LIBRARY_PATH=str(software) + ':' + environment.get('LD_LIBRARY_PATH', ''))
    environment.pop('P4REDUCE_GLOBAL_CONFIG', None)
    provenance = {'arguments': {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                  'binary': fingerprint(executable), 'libraries': libraries, 'runner': fingerprint(runner),
                  'config': fingerprint(Path(command[command.index('--config')+1])),
                  'template': fingerprint(args.template.resolve()),
                  'inputs': [fingerprint(Path(p)) for p in (root / 'inputs.txt').read_text().splitlines()],
                  'command': command, 'environment': {k: environment[k] for k in
                                                    ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'LD_LIBRARY_PATH')}}
    (stage / 'manifest.json').write_text(json.dumps(provenance, indent=2) + '\n')
    (stage / 'command.json').write_text(json.dumps(command, indent=2) + '\n')
    with (stage / 'run.log').open('w') as stream:
        subprocess.run(['/usr/bin/time', '-f', 'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M',
                        '-o', str(stage / 'resource_usage.txt'), *command], env=environment, cwd=stage,
                       stdout=stream, stderr=subprocess.STDOUT, check=True)
    reference_path = args.science_reference.resolve() if args.science_reference else root / 'baseline' / 'finim.fits'
    baseline = fits.getdata(reference_path)
    science = fits.getdata(stage / 'finim.fits')
    unchanged = (science.shape == baseline.shape and science.dtype == baseline.dtype and
                 science.tobytes() == baseline.tobytes())
    header = fits.getheader(stage / 'finim_outputs' / 'p4PSF_manifest.fits') if args.method != 'baseline' else {}
    result = {'science_bitwise_equal': unchanged, 'science_reference': fingerprint(reference_path),
              'science_max_absolute_difference': float(np.nanmax(np.abs(science - baseline))),
              'response_header': {k: v for k, v in header.items() if k.startswith('P4 PSF')},
              'resources': (stage / 'resource_usage.txt').read_text()}
    (stage / 'complete.json').write_text(json.dumps(result, indent=2) + '\n')
    if not unchanged:
        raise RuntimeError('response generation changed science')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
