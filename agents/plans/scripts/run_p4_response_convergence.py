#!/usr/bin/env python3
"""Sweep paired source amplitudes through production P4 local, mean-combined stamps.

Inputs must already be preprocessed P4 images. This measures final-stamp amplitude
dependence, not a known exact derivative or detector-level eigengaps. Use the
P4PCAResponse Catch2 tests for the independent FP64 numerical oracle.
"""

from __future__ import annotations

import argparse
import csv
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
    """Identify a preserved experiment input by absolute path and content hash."""
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "sha256": digest}


def read_stamp(directory: Path, modes: list[float], size: int, trial: tuple, precision: str | None = None) -> tuple:
    """Validate the requested local product and its explicit validity map."""
    data, header = fits.getdata(directory / "finim.fits", header=True)
    validity, validity_header = fits.getdata(directory / "finim_outputs/finim_local_validity.fits", header=True)
    if data.ndim == 2:
        data = data[None, :, :]
    if validity.ndim == 2:
        validity = validity[None, :, :]
    expected_shape = (len(modes), size, size)
    if data.shape != expected_shape or validity.shape != expected_shape:
        raise ValueError(f"unexpected local product shape in {directory}")
    if str(header["P4 PRODUCT ROLE"]).strip() != "LOCAL_RESIDUAL":
        raise ValueError(f"not a local residual: {directory}")
    if precision is not None and str(header["P4POLCY"]).strip() != precision.removeprefix("P4-"):
        raise ValueError(f"PCA precision policy does not match: {directory}")
    if str(validity_header["P4 PRODUCT ROLE"]).strip() != "LOCAL_VALIDITY":
        raise ValueError(f"not a local validity map: {directory}")
    actual_modes = [float(value) for value in str(header["P4 MODE FRACTIONS"]).split(",")]
    if len(actual_modes) != len(modes) or not np.allclose(actual_modes, modes, rtol=2e-6, atol=0):
        raise ValueError(f"mode fractions do not match: {directory}")
    for key, expected in zip(("FAKESEP", "FAKEPA", "FAKECONT"), trial):
        if not np.isclose(float(header[key]), expected, rtol=2e-6, atol=0):
            raise ValueError(f"{key} does not match: {directory}")
    if not np.all(np.isfinite(validity)) or not np.all((validity == 0) | (validity == 1)):
        raise ValueError(f"invalid validity map: {directory}")
    support = validity != 0
    if not np.all(np.isfinite(data[support])):
        raise ValueError(f"nonfinite residual on declared valid support: {directory}")
    origin_keys = ("P4 LOCAL ORIGIN ROW", "P4 LOCAL ORIGIN COLUMN")
    origin = tuple(int(header[key]) for key in origin_keys)
    if origin != tuple(int(validity_header[key]) for key in origin_keys):
        raise ValueError(f"residual/validity origins differ: {directory}")
    return np.asarray(data, dtype=np.float64), support, origin


def ratio(numerator: float, denominator: float) -> float | None:
    """Represent a zero-denominator diagnostic as missing rather than spurious agreement."""
    return float(numerator / denominator) if denominator > 0 else None


def summarize(position: tuple, modes: list, amplitudes: list, baseline: tuple, pairs: list) -> list[dict]:
    """Compare adjacent amplitudes on common support without declaring either one exact."""
    rows = []
    for mode, fraction in enumerate(modes):
        for index in range(len(amplitudes) - 1):
            large, small = pairs[index:index + 2]
            support = baseline[1][mode] & large[1][mode] & small[1][mode]
            coarse, fine = large[0][mode][support], small[0][mode][support]
            coarse_norm, fine_norm = np.linalg.norm(coarse), np.linalg.norm(fine)
            rows.append({
                "radius_pixels": position[0],
                "position_angle_degrees": position[1],
                "mode_fraction": fraction,
                "amplitude_large": amplitudes[index],
                "amplitude_small": amplitudes[index + 1],
                "common_pixels": int(support.sum()),
                "total_pixels": int(support.size),
                "support_changed_pixels": int(np.count_nonzero(large[1][mode] != small[1][mode])),
                "sign_support_changed_large": int(large[3][mode].sum()),
                "sign_support_changed_small": int(small[3][mode].sum()),
                "relative_change": ratio(np.linalg.norm(coarse - fine), fine_norm),
                "cosine": ratio(np.dot(coarse, fine), coarse_norm * fine_norm),
                "projection_scale_large_to_small": ratio(np.dot(coarse, fine), fine_norm**2),
                "small_response_norm": float(fine_norm),
                "relative_even_large": ratio(np.linalg.norm(large[2][mode][support]), coarse_norm),
                "relative_even_small": ratio(np.linalg.norm(small[2][mode][support]), fine_norm),
            })
    return rows


def main() -> None:
    """Run a fresh, fully logged experiment using fixed direct-P4 and mean-combination settings."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path,
                        help="base P4 geometry and angle configuration; external file paths must be absolute")
    parser.add_argument("--inputs", required=True, type=Path, help="directory of preprocessed image FITS files")
    parser.add_argument("--pattern", default="*.fits", help="input filename glob, sorted lexically")
    parser.add_argument("--psf", required=True, type=Path, help="stored unit-contrast PSF in P4 input units")
    parser.add_argument("--binary", default="p4Reduce")
    parser.add_argument("--precision", choices=["P4-D64", "P4-M32D64"],
                        help="requires p4ReductionPrecisionBenchmark; image storage remains FP32")
    parser.add_argument("--capture-detectors", action="store_true",
                        help="capture detector inputs at zero/unit amplitude and residuals at every amplitude")
    parser.add_argument("--output", required=True, type=Path,
                        help="new experiment directory; existing paths are refused")
    parser.add_argument("--radii", nargs="+", type=float, required=True)
    parser.add_argument("--angles", nargs="+", type=float, required=True)
    parser.add_argument("--modes", nargs="+", type=float, default=[0.05, 0.1, 0.2])
    parser.add_argument("--amplitudes", nargs="+", type=float, required=True, help="positive half-amplitudes")
    parser.add_argument("--stamp-size", type=int, default=11)
    args = parser.parse_args()
    if args.capture_detectors and args.precision is None:
        parser.error("--capture-detectors requires --precision and the experimental benchmark binary")
    if args.stamp_size <= 0 or args.stamp_size % 2 != 1:
        parser.error("--stamp-size must be positive and odd")
    numeric_groups = (args.radii, args.angles, args.modes, args.amplitudes)
    if not all(np.isfinite(value) for values in numeric_groups for value in values):
        parser.error("numeric arguments must be finite")
    if any(value < 0 for value in args.radii) or any(value <= 0 for value in args.amplitudes):
        parser.error("radii must be nonnegative and amplitudes positive")
    if any(not 0 < value <= 1 for value in args.modes) or args.modes != sorted(set(args.modes)):
        parser.error("mode fractions must be distinct, increasing, and in (0,1]")
    amplitudes = sorted(set(args.amplitudes), reverse=True)
    if len(amplitudes) < 2:
        parser.error("at least two distinct amplitudes are required")
    if any(value > np.finfo(np.float32).max or np.float32(value) == 0 for value in amplitudes):
        parser.error("amplitudes must remain finite and nonzero in production FP32 input")
    resolved_binary = shutil.which(args.binary)
    if resolved_binary is None:
        parser.error(f"executable not found: {args.binary}")
    binary = Path(resolved_binary).resolve()
    config, psf, output = args.config.resolve(), args.psf.resolve(), args.output.resolve()
    inputs = sorted(path.resolve() for path in args.inputs.glob(args.pattern))
    if len(inputs) < 2:
        parser.error("at least two input FITS images are required")
    for path in inputs:
        if any(character.isspace() for character in str(path)):
            parser.error("the production input file-list reader requires paths without whitespace")
        if fits.getheader(path).get("HCIREDUCE PREPROCESSED") != 1:
            parser.error(f"input is not a marked preprocessed product: {path}")
    manifest = {
        "schema": 1,
        "scope": "direct detector P4, post-preprocessing injection, fixed mean combination",
        "reference": "adjacent amplitudes only; no exact final-stamp derivative assumed",
        "implicit_configs": "disabled: unset P4REDUCE_GLOBAL_CONFIG and run in each new trial directory",
        "arguments": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
        "amplitudes_descending": amplitudes,
        "binary": fingerprint(binary), "config": fingerprint(config), "psf": fingerprint(psf),
        "runner": fingerprint(Path(__file__).resolve()),
        "inputs": [fingerprint(path) for path in inputs],
        "thread_environment": {
            key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        },
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.copyfile(config, output / "base.conf")
    file_list = output / "inputs.txt"
    file_list.write_text("".join(f"{path}\n" for path in inputs))
    fixed = {
        "mode": "normal", "input.directory": "", "input.fileList": str(file_list),
        "input.deleteFront": 0, "input.deleteBack": 0, "input.qualityFile": "",
        "input.qualityThreshold": 0, "input.imSize": 0,
        "coadd.method": "none", "preProcess.skip": "true", "preProcess.only": "false",
        "preProcess.outputPrefix": "", "adi.postMedSub": "false", "adi.excludeMethod": "none",
        "solver.exclusionSolver": "explicitRefit", "solver.deletionBackend": "leadingCovariance",
        "p4.regressionFrame": "detector", "p4.temporalPredictor": "none", "p4.numberImages": 0,
        "p4.modeFractions": ",".join(map(str, args.modes)), "p4.localStampSize": args.stamp_size,
        "p4.writeDiagnostics": "false", "p4Optimize.enabled": "false",
        "psfResponse.file": "", "psfResponse.outputModels": "false", "psfResponse.filter": "false",
        "fake.method": "single", "fake.fileName": str(psf), "fake.scaleFileName": "",
        "fake.subtractPlanet": "false", "combine.method": "mean", "combine.weightFile": "",
        "combine.minGoodFract": 0, "combine.noDerotate": "false",
        "output.fileName": "finim.fits", "output.exactFName": "true", "output.outputPSFSub": "false",
        "showTiming": "true",
    }
    common = [str(binary), "--config", str(config)]
    for key, value in fixed.items():
        common.append(f"--{key}={value}")
    if args.precision is not None:
        common.extend([f"--precision={args.precision}", "--p4.memoryFraction=0"])
    environment = os.environ.copy()
    environment.pop("P4REDUCE_GLOBAL_CONFIG", None)

    def run(directory: Path, radius: float, angle: float, amplitude: float) -> tuple:
        """Run and validate one baseline or signed trial, preserving its exact command and log."""
        directory.mkdir()
        command = common + ["--output.directory", str(directory), "--fake.sep", str(radius),
                            "--fake.PA", str(angle), "--fake.contrast", str(amplitude)]
        if args.capture_detectors:
            command.append(f"--capture-detectors={directory / 'detectors'}")
            command.append(f"--capture-inputs={'true' if amplitude in (0, 1) else 'false'}")
        (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        print(f"{directory.name}: radius={radius:g} PA={angle:g} contrast={amplitude:g}", flush=True)
        start = time.monotonic()
        with (directory / "run.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, cwd=directory, env=environment, check=True)
        (directory / "timing.json").write_text(json.dumps({"wall_seconds": time.monotonic() - start}) + "\n")
        return read_stamp(directory, args.modes, args.stamp_size, (radius, angle, amplitude), args.precision)

    rows = []
    position_index = 0
    for radius in args.radii:
        for angle in args.angles:
            position_dir = output / f"position_{position_index:03d}"
            position_index += 1
            position_dir.mkdir()
            baseline = run(position_dir / "baseline", radius, angle, 0)
            if args.capture_detectors:
                run(position_dir / "unit_source", radius, angle, 1)
            pairs = []
            for index, amplitude in enumerate(amplitudes):
                positive = run(position_dir / f"amplitude_{index:02d}_positive", radius, angle, amplitude)
                negative = run(position_dir / f"amplitude_{index:02d}_negative", radius, angle, -amplitude)
                if positive[2] != baseline[2] or negative[2] != baseline[2]:
                    raise ValueError("stamp origin changed across signed trials")
                support = positive[1] & negative[1]
                response = (positive[0] - negative[0]) / (2 * amplitude)
                even = (positive[0] + negative[0] - 2 * baseline[0]) / (2 * amplitude)
                response[~support] = np.nan
                pairs.append((response, support, even, positive[1] != negative[1]))
            rows.extend(summarize((radius, angle), args.modes, amplitudes, baseline, pairs))
            (output / "summary.json").write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")
            with (output / "summary.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
    (output / "complete.json").write_text(json.dumps({"comparisons": len(rows)}) + "\n")
    print(f"Wrote {len(rows)} adjacent-amplitude comparisons to {output / 'summary.csv'}")


if __name__ == "__main__":
    main()
