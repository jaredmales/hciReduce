#!/usr/bin/env python3
"""Freeze and run the positive-injection evaluation at preassigned Step-5 sites.

The run phase is self-contained in its software directory. It completes all full
reductions, applies four frozen production filters, and writes raw recovery and
conditional-coverage summaries without subtracting the baseline science image.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from astropy.io import fits

from run_p4_step5_full_injections import fingerprint, write_json

MODELS = ('identity', 'diagonal', 'pca0', 'pca')
OFFSETS = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))


def prepare(args: argparse.Namespace) -> None:
    """Freeze software, thresholds, exclusions, and 18 full-image trials before running them."""
    root = args.output.resolve()
    holdout = args.holdout.resolve()
    nulls = json.loads((holdout / 'null_results.json').read_text())
    protocol = json.loads((holdout / 'protocol.json').read_text())
    design = json.loads((holdout / 'injection_design.json').read_text())
    if nulls['protocol'] != fingerprint(holdout / 'protocol.json'):
        raise RuntimeError('protocol changed after calibration')
    jobs = []
    for radius, block in design['sites']:
        site = next(trial for trial in nulls['trials'] if trial['role'] == 'evaluation' and
                    trial['nominal_radius'] == radius and trial['angular_block'] == block)
        if not site['common_eligible']:
            raise RuntimeError('a preassigned injection site is ineligible; do not replace it based on scores')
        scale = nulls['models']['pca0']['threshold'] * site['models']['pca0']['center_conditional_sigma']
        if not np.isfinite(scale) or scale <= 0:
            raise RuntimeError('calibration cannot supply a positive finite injection scale')
        for level, multiplier in enumerate(design['brightness_multipliers']):
            jobs.append({'name': f'evaluation_r{radius}_b{block}_l{level}', 'row': site['row'],
                'column': site['column'], 'contrast': multiplier * scale,
                'nominal_radius': radius, 'angular_block': block, 'brightness_multiplier': multiplier,
                'reference_contrast_scale': scale, 'null_trial': site['name']})
    root.mkdir(parents=True, exist_ok=False)
    software = root / 'software'
    software.mkdir()
    for name in ('hciAnalyze', 'libhcireduce.so'):
        expected = next(r for r in protocol['software'] if Path(r['path']).name == name)
        if fingerprint(Path(expected['path'])) != expected:
            raise RuntimeError('calibration software changed')
        shutil.copy2(expected['path'], software / name)
    for name in (Path(__file__).name, 'run_p4_step5_full_injections.py'):
        shutil.copy2(Path(__file__).with_name(name), software / name)
    shutil.copytree(holdout / 'response', root / 'response')
    write_json(root / 'jobs.json', jobs)
    write_json(root / 'null_results.json', nulls)
    write_json(root / 'protocol.json', protocol)
    write_json(root / 'injection_design.json', design)
    command = [sys.executable, str(software / 'run_p4_step5_full_injections.py'), 'prepare',
        '--output', str(root / 'reductions'), '--products', str(args.products.resolve()),
        '--experiment', str(args.experiment.resolve()), '--trials', str(root / 'jobs.json'),
        '--purpose', 'held-out positive-injection evaluation with frozen null thresholds']
    subprocess.run(command, check=True)
    # Match the dependent libraries used by the frozen reductions and original response software.
    reduction_software = root / 'reductions/software'
    for name in ('libmxlib.so', 'liblapack.so.3', 'libopenblas.so.0'):
        shutil.copy2(reduction_software / name, software / name)
    records = [fingerprint(p) for p in sorted(software.iterdir())]
    records += [fingerprint(p) for p in sorted((root / 'response').glob('*.fits'))]
    records += [fingerprint(root / name) for name in ('jobs.json', 'null_results.json', 'protocol.json', 'injection_design.json')]
    write_json(root / 'manifest.json', {'schema': 1, 'frozen_records': records,
        'holdout_null_results': fingerprint(holdout / 'null_results.json'),
        'source_frame_count': 621, 'full_reductions': len(jobs), 'models': list(MODELS),
        'contrast_and_coverage': 'raw positive science, fixed injected position; no baseline subtraction',
        'detection': 'maximum signed score in the preassigned five-pixel search aperture, strict frozen-threshold exceedance',
        'analysis_cpu_ids': [12, 13], 'analysis_openmp_threads': 2, 'analysis_blas_threads': 1})
    print(root, flush=True)


def analyze(root: Path, manifest: dict) -> None:
    """Apply frozen production filters and summarize every preassigned trial, including failures."""
    jobs = json.loads((root / 'jobs.json').read_text())
    protocol = json.loads((root / 'protocol.json').read_text())
    nulls = json.loads((root / 'null_results.json').read_text())
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       LD_LIBRARY_PATH=str(root / 'software') + ':' + environment.get('LD_LIBRARY_PATH', ''))
    environment.pop('HCIANALYZE_GLOBAL_CONFIG', None)
    rows = []
    for job in jobs:
        for model in MODELS:
            directory = root / 'analysis' / job['name'] / model
            record_path = directory / 'measurement.json'
            if record_path.exists():
                record = json.loads(record_path.read_text())
                for expected in record['products']:
                    if fingerprint(Path(expected['path'])) != expected:
                        raise RuntimeError('a completed analysis product changed')
                rows.append(record)
                continue
            directory.mkdir(parents=True)
            science = directory / 'science.fits'
            shutil.copy2(root / 'reductions' / job['name'] / 'finim.fits', science)
            command = ['taskset', '-c', ','.join(map(str, manifest['analysis_cpu_ids'])),
                str(root / 'software/hciAnalyze'), f'--file={science}',
                f'--filter.psfResponse={root / "response/p4PSF_manifest.fits"}',
                f'--noise.model={"pca" if model == "pca0" else model}', '--noise.only=true',
                '--noise.outputDiagnostics=true', '--noise.minimumSamples=8',
                f'--noise.maximumModes={0 if model == "pca0" else 3}', '--noise.floorFraction=0.1',
                '--noise.arcStep=5', '--noise.guardRadius=0', '--noise.exactExclusion=true',
                '--lambdaD=3.6', '--planet.sep=11.782', '--planet.PA=262.051', '--planet.R=7.3',
                '--filter.hpfGaussFW=0', '--filter.lpfGaussFW=0']
            command += [f'--noise.exclude{key}=' + ','.join(str(circle[index]) for circle in protocol['extra_exclusions'])
                        for index, key in enumerate(('Rows', 'Columns', 'Radii'))]
            write_json(directory / 'command.json', command)
            write_json(root / 'state.json', {'status': 'analyzing', 'trial': job['name'], 'model': model,
                                           'measurements_finished': len(rows)})
            with (directory / 'run.log').open('w') as log:
                subprocess.run(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
            maps = {role: fits.getdata(directory / f'science_{role}.fits').squeeze()
                    for role in ('psf_score', 'psf_amplitude', 'psf_sigma', 'noise_status', 'noise_samples')}
            x, y = job['row'], job['column']
            statuses = [float(maps['noise_status'][y+dy, x+dx]) for dx, dy in OFFSETS]
            valid = all(s == 0 for s in statuses)
            amplitude = float(maps['psf_amplitude'][y, x]) if maps['noise_status'][y, x] == 0 else None
            sigma = float(maps['psf_sigma'][y, x]) if amplitude is not None else None
            score = max(float(maps['psf_score'][y+dy, x+dx]) for dx, dy in OFFSETS) if valid else None
            threshold = nulls['models'][model]['threshold']
            record = {'trial': job, 'model': model, 'search_valid': valid,
                'statuses': [s if np.isfinite(s) else None for s in statuses],
                'search_score': score, 'threshold': threshold,
                'detected': valid and score > threshold,
                'center_amplitude': amplitude, 'center_conditional_sigma': sigma,
                'center_training_samples': int(maps['noise_samples'][y, x])
                if np.isfinite(maps['noise_samples'][y, x]) else None,
                'contrast_bias_fraction': amplitude / job['contrast'] - 1 if amplitude is not None else None,
                'inside_one_conditional_sigma': abs(amplitude-job['contrast']) <= sigma if amplitude is not None else False,
                'products': [fingerprint(p) for p in sorted(directory.glob('*.fits'))]}
            write_json(record_path, record)
            rows.append(record)
    groups = []
    for multiplier in sorted({job['brightness_multiplier'] for job in jobs}):
        for model in MODELS:
            selected = [r for r in rows if r['model'] == model and r['trial']['brightness_multiplier'] == multiplier]
            biases = [r['contrast_bias_fraction'] for r in selected if r['contrast_bias_fraction'] is not None]
            groups.append({'model': model, 'brightness_multiplier': multiplier, 'trials': len(selected),
                'detections': sum(r['detected'] for r in selected),
                'invalid_searches': sum(not r['search_valid'] for r in selected),
                'median_contrast_bias_fraction': float(np.median(biases)) if biases else None,
                'bias_available_trials': len(biases),
                'inside_one_conditional_sigma': sum(r['inside_one_conditional_sigma'] for r in selected)})
    write_json(root / 'results.json', {'groups': groups, 'measurements': rows, 'null_models': nulls['models'],
        'caveats': ['Six fixed spatial sites; overlapping field noise is not an independent ensemble.',
                    'Identity conditional sigma assumes unit pixel covariance and is not a fitted noise uncertainty.',
                    'No baseline subtraction; reported contrast errors include residual noise and finite-source response bias.',
                    'Null thresholds target a pooled 5% operating point; achieved evaluation rates and block sensitivity are separate.']})
    lines = ['# Step 5 positive-injection results', '',
        'Each brightness level has six preassigned positions. Results use raw positive-injection images, without '
        'baseline subtraction. Invalid searches count as nondetections. Contrast errors and conditional intervals '
        'use the exact injected position, without selecting a fitted peak.', '',
        '| Model | Brightness / reference threshold | Detections | Invalid searches | Median contrast error | Within one conditional sigma |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for group in groups:
        bias = group['median_contrast_bias_fraction']
        bias_text = f'{100*bias:+.1f}%' if bias is not None else 'unavailable'
        lines.append(f"| {group['model']} | {group['brightness_multiplier']:g} | "
            f"{group['detections']}/{group['trials']} | {group['invalid_searches']} | {bias_text} | "
            f"{group['inside_one_conditional_sigma']}/{group['trials']} |")
    lines += ['', 'Identity sigma assumes unit pixel covariance; it is not a fitted noise uncertainty. '
        'The spatial trials overlap and provide a small correlated sample. Consult `results.json` for the separate '
        'null exceedance rates, block sensitivity, individual measurements, and file hashes. These results require '
        'scientific review before selecting a preferred noise model.', '']
    (root / 'results.md').write_text('\n'.join(lines))


def run_locked(args: argparse.Namespace) -> None:
    """Complete the frozen reduction/analysis queue and persist final status or a visible failure."""
    root = args.output.resolve()
    manifest = json.loads((root / 'manifest.json').read_text())
    for expected in manifest['frozen_records']:
        if fingerprint(Path(expected['path'])) != expected:
            raise RuntimeError('frozen evaluation input changed: ' + expected['path'])
    try:
        write_json(root / 'state.json', {'status': 'reducing', 'pid': os.getpid()})
        command = [sys.executable, str(root / 'reductions/software/run_p4_step5_full_injections.py'),
                   'run', '--output', str(root / 'reductions')]
        subprocess.run(command, check=True)
        analyze(root, manifest)
        for expected in manifest['frozen_records']:
            if fingerprint(Path(expected['path'])) != expected:
                raise RuntimeError('frozen evaluation input changed during the run: ' + expected['path'])
        write_json(root / 'complete.json', {'reductions': manifest['full_reductions'],
            'model_measurements': manifest['full_reductions'] * len(MODELS),
            'frozen_inputs_unchanged': True, 'results': fingerprint(root / 'results.json')})
        write_json(root / 'state.json', {'status': 'complete', 'results': str(root / 'results.json')})
    except Exception as error:
        write_json(root / 'state.json', {'status': 'failed', 'error': str(error)})
        raise


def run(args: argparse.Namespace) -> None:
    """Prevent duplicate supervisors from running the same expensive frozen queue."""
    with (args.output.resolve() / 'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        run_locked(args)


def main() -> None:
    """Prepare or execute an unattended, self-contained scientific evaluation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--holdout', type=Path)
    parser.add_argument('--products', type=Path)
    parser.add_argument('--experiment', type=Path)
    args = parser.parse_args()
    if args.action == 'prepare' and any(getattr(args, name) is None for name in ('holdout', 'products', 'experiment')):
        parser.error('prepare requires --holdout, --products, and --experiment')
    (prepare if args.action == 'prepare' else run)(args)


if __name__ == '__main__':
    main()
