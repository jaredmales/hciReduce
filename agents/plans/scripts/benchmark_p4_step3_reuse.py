#!/usr/bin/env python3
"""Compare P4 factor reuse, source caching, and sparse temporal products.

Use separate frozen software directories on fixed 24/96-frame experiments.
Each serial trial retains /usr/bin/time results and checks every template against
the FP64 paired-refit reference. Report medians of three runs, not selected minima.
Pin a homogeneous CPU set on hybrid processors to keep native arithmetic and
performance comparable between separately launched processes.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys

import numpy as np
from astropy.io import fits

from compare_p4_refit_difference import resource_usage
from run_p4_response_convergence import fingerprint


def timed_run(command: list[str], stage: Path, environment: dict) -> None:
    """Run one preserved command with resource accounting and a separate log."""
    stage.mkdir()
    (stage / 'command.json').write_text(json.dumps(command, indent=2)+'\n')
    print(stage.name, flush=True)
    with (stage / 'run.log').open('w') as log:
        subprocess.run(['/usr/bin/time', '-f',
            'wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M',
            '-o', str(stage / 'resource_usage.txt'), *command], cwd=stage, env=environment,
            stdout=log, stderr=subprocess.STDOUT, check=True)


def replaced(command: list[str], overrides: dict) -> list[str]:
    """Replace equals-style options without relying on duplicate-option order."""
    return [arg for arg in command if arg.split('=', 1)[0].removeprefix('--') not in overrides] + [
        f'--{key}={value}' for key, value in overrides.items()]


def refit_controls(args: argparse.Namespace, output: Path) -> list[dict]:
    """Time production-precision and FP64 paired refits on the same subsets."""
    records = []
    software = args.cached_software.resolve()
    executable = str(software / 'p4ReductionPrecisionBenchmark')
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       LD_LIBRARY_PATH=str(software)+':'+environment.get('LD_LIBRARY_PATH', ''))
    environment.pop('P4REDUCE_GLOBAL_CONFIG', None)
    policies = [('native_refit', 'P4-M32D64', 0.0047639259293563909), ('fp64_refit', 'P4-D64', 1e-5)]
    for frames in (24, 96):
        reference = args.reference_root.resolve() / f'refit{frames}_D64'
        previous = json.loads((reference / 'manifest.json').read_text())['reference_manifest']
        last_amplitude = len(previous['amplitudes_descending'])-1
        common = json.loads((reference / f'amplitude_{last_amplitude:03d}/command.json').read_text())
        common[0] = executable
        analytic = output / f'frames{frames}_cached_sparse_batch32_repeat0'
        baseline = {'P4-M32D64': fits.getdata(analytic / 'baseline/finim.fits')}
        baseline_stage = output / f'frames{frames}_fp64_baseline'
        baseline_command = json.loads((reference / 'baseline/command.json').read_text())
        baseline_command[0] = executable
        baseline_command = replaced(baseline_command, {'precision': 'P4-D64',
            'output.directory': baseline_stage, 'input.fileList': analytic / 'inputs.txt'})
        timed_run(baseline_command, baseline_stage, environment)
        baseline['P4-D64'] = fits.getdata(baseline_stage / 'finim.fits')
        templates = [fits.getdata(path) for path in sorted(
            (analytic / 'analytic/finim_outputs').glob('response_model_*.fits'))]
        for label, precision, amplitude in policies:
            for repeat in range(args.repeats):
                stage = output / f'frames{frames}_{label}_repeat{repeat}'
                command = replaced(common, {'precision': precision, 'output.directory': stage,
                    'input.fileList': analytic / 'inputs.txt', 'psfResponse.refitContrast': amplitude})
                timed_run(command, stage, environment)
                science = fits.getdata(stage / 'finim.fits')
                if (science.shape != baseline[precision].shape or science.dtype != baseline[precision].dtype or
                    science.tobytes() != baseline[precision].tobytes()):
                    raise RuntimeError('paired-refit control changed same-precision science')
                paths = sorted((stage / 'finim_outputs').glob('response_model_*.fits'))
                if len(paths) != len(templates):
                    raise RuntimeError('paired-refit control has different mode counts')
                errors = []
                for path, template in zip(paths, templates):
                    model = fits.getdata(path)
                    if not np.array_equal(np.isfinite(model), np.isfinite(template)):
                        raise RuntimeError('paired-refit control changed response support')
                    norms = np.sqrt(np.nansum(template.astype(float)**2, axis=(1, 2)))
                    difference = np.sqrt(np.nansum((model.astype(float)-template)**2, axis=(1, 2)))
                    errors.extend((difference[norms > 0]/norms[norms > 0]).tolist())
                records.append({'frames': frames, 'repeat': repeat, 'precision': precision,
                    'cpu_affinity': sorted(os.sched_getaffinity(0)),
                    'half_contrast': amplitude, 'science_bitwise_equal': True,
                    'maximum_stamp_relative_error_vs_analytic': max(errors),
                    'resources': resource_usage(stage / 'resource_usage.txt')})
                (output / 'refit_trials.json').write_text(json.dumps(records, indent=2)+'\n')
    summary = []
    for frames in (24, 96):
        for _, precision, _ in policies:
            runs = [row for row in records if (row['frames'], row['precision']) == (frames, precision)]
            summary.append({'frames': frames, 'precision': precision, 'trials': len(runs),
                'cpu_affinity': sorted(os.sched_getaffinity(0)),
                'half_contrast': runs[0]['half_contrast'], 'science_bitwise_equal': True,
                'maximum_stamp_relative_error_vs_analytic': max(
                    row['maximum_stamp_relative_error_vs_analytic'] for row in runs),
                'median_wall_seconds': statistics.median(row['resources']['wall_seconds'] for row in runs),
                'median_user_seconds': statistics.median(row['resources']['user_seconds'] for row in runs),
                'median_peak_rss_kib': statistics.median(row['resources']['maximum_rss_kib'] for row in runs)})
    (output / 'refit_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


def main() -> None:
    """Execute controlled variants serially and require identical float products."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-root', type=Path, required=True)
    parser.add_argument('--uncached-software', type=Path, required=True)
    parser.add_argument('--cached-software', type=Path, required=True)
    parser.add_argument('--dense-cached-software', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--cpus', type=int, nargs='+', help='CPU IDs inherited by every trial, e.g. --cpus 0 2')
    args = parser.parse_args()
    if args.cpus is not None:
        os.sched_setaffinity(0, args.cpus)
    output = args.output.resolve()
    output.mkdir(exist_ok=False)
    (output / 'benchmark_manifest.json').write_text(json.dumps({
        'cpu_affinity': sorted(os.sched_getaffinity(0)), 'threads': args.threads,
        'blas_threads': 1, 'repeats': args.repeats, 'script': fingerprint(Path(__file__).resolve()),
        'software_directories': {name: str(path.resolve()) for name, path in
                                [('uncached_dense', args.uncached_software),
                                 ('cached_dense', args.dense_cached_software),
                                 ('cached_sparse', args.cached_software)]}}, indent=2)+'\n')
    records = []
    expected = {}
    variants = [('uncached_dense', args.uncached_software), ('cached_dense', args.dense_cached_software),
                ('cached_sparse', args.cached_software)]
    for frames in (24, 96):
        for variant, directory in variants:
            software = directory.resolve()
            executable = software / 'p4ReductionPrecisionBenchmark'
            environment = os.environ.copy()
            environment.update(OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                               LD_LIBRARY_PATH=str(software)+':'+environment.get('LD_LIBRARY_PATH', ''))
            for batch in (1, 32):
                for repeat in range(args.repeats):
                    name = f'frames{frames}_{variant}_batch{batch}_repeat{repeat}'
                    stage = output / name
                    command = [sys.executable, str(Path(__file__).with_name('run_p4_analytic_products.py')),
                        '--reference', str(args.reference_root.resolve() / f'refit{frames}_D64'),
                        '--binary', str(executable), '--precision', 'P4-M32D64', '--batch-size', str(batch),
                        '--output', str(stage)]
                    print(name, flush=True)
                    subprocess.run(command, env=environment, check=True)
                    complete = json.loads((stage / 'complete.json').read_text())
                    model_paths = sorted((stage / 'analytic/finim_outputs').glob('response_model_*.fits'))
                    models = [fits.getdata(path) for path in model_paths]
                    key = (frames, variant)
                    identical = True
                    if key not in expected:
                        expected[key] = models
                    elif any(not np.array_equal(a, b, equal_nan=True) for a, b in zip(models, expected[key])):
                        raise RuntimeError('batching changed response templates within an arithmetic variant')
                    reference_models = expected[(frames, 'uncached_dense')]
                    relative_errors = []
                    for actual, reference in zip(models, reference_models):
                        if not np.array_equal(np.isfinite(actual), np.isfinite(reference)):
                            raise RuntimeError('arithmetic variants changed response support')
                        identical &= np.array_equal(actual, reference, equal_nan=True)
                        numerator = np.sqrt(np.nansum((actual.astype(float)-reference)**2, axis=(1,2)))
                        denominator = np.sqrt(np.nansum(reference.astype(float)**2, axis=(1,2)))
                        relative_errors.extend((numerator[denominator>0]/denominator[denominator>0]).tolist())
                    maximum_error = max(relative_errors)
                    if maximum_error > 8*np.finfo(np.float32).eps:
                        raise RuntimeError('sparse arithmetic differs beyond the FP32 product tolerance')
                    records.append({'frames': frames, 'variant': variant, 'requested_batch': batch,
                        'cpu_affinity': sorted(os.sched_getaffinity(0)),
                        'bitwise_identical_vs_dense': bool(identical), 'maximum_stamp_relative_error_vs_dense': maximum_error,
                        'repeat': repeat, 'binary': fingerprint(executable), **complete,
                        'resources': resource_usage(stage / 'analytic/resource_usage.txt')})
                    (output / 'trials.json').write_text(json.dumps(records, indent=2)+'\n')
    summary = []
    for frames in (24, 96):
        for variant, _ in variants:
            for batch in (1, 32):
                runs = [r for r in records if (r['frames'], r['variant'], r['requested_batch']) ==
                        (frames, variant, batch)]
                summary.append({'frames': frames, 'variant': variant, 'requested_batch': batch,
                    'cpu_affinity': sorted(os.sched_getaffinity(0)),
                    'bitwise_identical_vs_dense': all(r['bitwise_identical_vs_dense'] for r in runs),
                    'maximum_stamp_relative_error_vs_dense': max(r['maximum_stamp_relative_error_vs_dense'] for r in runs),
                    'trials': len(runs), 'baseline_factor_count': runs[0]['baseline_factor_count'],
                    'realized_batch': runs[0]['realized_batch_size'], 'bitwise_identical_templates': True,
                    'median_wall_seconds': statistics.median(r['resources']['wall_seconds'] for r in runs),
                    'median_user_seconds': statistics.median(r['resources']['user_seconds'] for r in runs),
                    'median_peak_rss_kib': statistics.median(r['resources']['maximum_rss_kib'] for r in runs)})
    refit_controls(args, output)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
