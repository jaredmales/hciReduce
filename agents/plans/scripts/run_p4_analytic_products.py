#!/usr/bin/env python3
"""Compare full analytic P4 products with a completed FP64 refit-product sweep.

Run an independent science baseline and an analytic response reduction. Preserve
commands, input/library fingerprints, per-mode diagnostics, and per-stamp errors.
The reference is the smallest finite amplitude, not an exact mathematical oracle.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from astropy.io import fits

from run_p4_refit_convergence import read_products
from run_p4_response_convergence import fingerprint, ratio


def main() -> None:
    """Run both reductions and require unchanged science and identical response support."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="completed P4-D64 refit-product sweep")
    parser.add_argument("--binary", required=True, type=Path, help="p4ReductionPrecisionBenchmark executable")
    parser.add_argument("--precision", choices=("P4-D64", "P4-M32D64"), required=True)
    parser.add_argument("--output", required=True, type=Path, help="new experiment directory")
    parser.add_argument("--batch-size", type=int, default=32, help="maximum analytic source measurements per batch")
    args = parser.parse_args()
    root, binary, output = args.reference.resolve(), args.binary.resolve(), args.output.resolve()
    previous = json.loads((root / "manifest.json").read_text())
    original = previous["reference_manifest"]
    settings = original["arguments"]
    if not (root / "complete.json").is_file() or settings["precision"] != "P4-D64":
        parser.error("a completed FP64 refit sweep is required")
    for record in [original["config"], original["psf"]] + original["inputs"]:
        if fingerprint(Path(record["path"])) != record:
            parser.error(f"reference input changed: {record['path']}")
    amplitudes = original["amplitudes_descending"]
    reference = root / f"amplitude_{len(amplitudes) - 1:03d}"
    expected, support, coordinates = read_products(
        reference, settings["modes"], settings["stamp_size"], "P4-D64", amplitudes[-1])
    libraries = {}
    for line in subprocess.check_output(["ldd", str(binary)], text=True).splitlines():
        if "=>" in line:
            name, resolved = line.split("=>", 1)
            if name.strip().startswith(("libhcireduce.", "libmxlib.", "libopenblas.", "liblapack.")):
                libraries[name.strip()] = fingerprint(Path(resolved.split(" (", 1)[0].strip()).resolve())
    output.mkdir(parents=True, exist_ok=False)
    (output / "manifest.json").write_text(json.dumps({
        "schema": 1, "reference": fingerprint(root / "manifest.json"), "reference_manifest": previous,
        "reference_products": [fingerprint(path) for path in sorted((reference / "finim_outputs").glob("*.fits"))],
        "binary": fingerprint(binary), "libraries": libraries, "precision": args.precision,
        "scripts": [fingerprint(Path(__file__).resolve().with_name(name)) for name in
                    (Path(__file__).name, "run_p4_refit_convergence.py", "run_p4_response_convergence.py")],
        "thread_environment": {key: os.environ.get(key) for key in
                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    }, indent=2) + "\n")
    environment = os.environ.copy()
    environment.pop("P4REDUCE_GLOBAL_CONFIG", None)
    inputs = output / "inputs.txt"
    inputs.write_text("".join(f"{record['path']}\n" for record in original["inputs"]))
    baseline, timings = None, {}
    for stage in ("baseline", "analytic"):
        source = root / "baseline" if stage == "baseline" else reference
        command = json.loads((source / "command.json").read_text())
        command[0] = str(binary)
        directory = output / stage
        directory.mkdir()
        overrides = {"precision": args.precision, "output.directory": directory, "input.fileList": inputs}
        if stage == "analytic":
            overrides.update({"psfResponse.method": "analytic", "psfResponse.refitContrast": 0,
                              "psfResponse.analyticGapTolerance": 0, "psfResponse.analyticBatchSize": args.batch_size})
        command = [arg for arg in command if arg.split("=", 1)[0].removeprefix("--") not in overrides]
        command += [f"--{key}={value}" for key, value in overrides.items()]
        (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        print(f"{args.precision} {stage}: {output}", flush=True)
        start = time.monotonic()
        with (directory / "run.log").open("w") as stream:
            subprocess.run(["/usr/bin/time", "-f", "wall_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmaximum_rss_kib=%M",
                            "-o", str(directory / "resource_usage.txt"), *command],
                           cwd=directory, env=environment, stdout=stream, stderr=subprocess.STDOUT, check=True)
        timings[stage] = time.monotonic() - start
        science = fits.getdata(directory / "finim.fits")
        if stage == "baseline":
            baseline = science.copy()
        elif not np.array_equal(science, baseline, equal_nan=True):
            raise ValueError("analytic response calculation changed science")
    products = output / "analytic" / "finim_outputs"
    actual_coordinates, coordinate_header = fits.getdata(products / "response_coordinates.fits", header=True)
    expected_coordinates = coordinates.copy()
    expected_coordinates[0] -= coordinate_header["P4 PSF COORDINATE ORIGIN ROW"]
    expected_coordinates[1] -= coordinate_header["P4 PSF COORDINATE ORIGIN COLUMN"]
    if not np.array_equal(actual_coordinates, expected_coordinates):
        raise ValueError("source coordinates differ from reference")
    diagnostics, diagnostic_header = fits.getdata(products / "response_measurement_diagnostics.fits", header=True)
    diagnostics = diagnostics.T.astype(np.int64)
    columns = str(diagnostic_header["P4 PSF DIAGNOSTIC COLUMNS"]).strip().split(",")
    if diagnostics.shape[1] != 11 or len(columns) != 11:
        raise ValueError("unexpected diagnostic columns")
    rows = []
    for mode, fraction in enumerate(settings["modes"]):
        actual, header = fits.getdata(products / f"response_model_{mode:04d}.fits", header=True)
        actual = actual.astype(np.float64)
        validity = fits.getdata(products / f"response_validity_{mode:04d}.fits").reshape(-1)
        if (header["P4 PSF PRODUCT SCHEMA"] != 8 or str(header["P4 PSF RESPONSE"]).strip() != "ANALYTIC_PROJECTOR"
                or str(header["P4 PSF RESPONSE PRECISION"]).strip() != "D64"
                or str(header["P4POLCY"]).strip() != args.precision.removeprefix("P4-")):
            raise ValueError("unexpected analytic provenance")
        actual_support = np.isfinite(actual) & (validity[:, None, None] != 0)
        if actual.shape != expected[mode].shape or not np.array_equal(actual_support, support[mode]):
            raise ValueError("analytic/reference dimensions or support differ")
        counts = diagnostics[diagnostics[:, 3] == mode, 4:].sum(axis=0)
        if counts[0] <= 0 or counts[1:].any():
            raise ValueError("reference dataset has unavailable or fallback responses")
        count_keys = ("ANALYTIC", "RANK INSUFFICIENT", "RANK BOUNDARY", "CUTOFF UNRESOLVED",
                      "FALLBACK ATTEMPTED", "FALLBACK ACCEPTED", "UNAVAILABLE")
        for key, count in zip(count_keys, counts):
            if int(header[f"P4 PSF {key} FITS"]) != count:
                raise ValueError("model header and measurement diagnostics differ")
        mask = actual_support
        a, b = actual[mask], expected[mode][mask]
        errors = [ratio(np.linalg.norm((x - y)[valid]), np.linalg.norm(y[valid]))
                  for x, y, valid in zip(actual, expected[mode], mask)]
        finite_errors = [error for error in errors if error is not None]
        if not finite_errors:
            raise ValueError("no nonzero response stamps to compare")
        rows.append({
            "mode_fraction": fraction, "reference_amplitude": amplitudes[-1],
            "common_pixels": int(mask.sum()), "compared_sources": len(finite_errors),
            "relative_error": ratio(np.linalg.norm(a - b), np.linalg.norm(b)),
            "max_source_relative_error": max(finite_errors), "source_relative_errors": errors,
            "cosine": ratio(np.dot(a, b), np.linalg.norm(a) * np.linalg.norm(b)),
            "projection_scale_analytic_to_refit": ratio(np.dot(a, b), np.dot(b, b)),
            "detector_counts": dict(zip(columns[4:], counts.tolist())),
        })
    (output / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    (output / "complete.json").write_text(json.dumps({
        "science_exactly_unchanged": True, "identical_support": True, "timing_seconds": timings,
        "requested_batch_size": args.batch_size,
        "baseline_factor_count": int(diagnostic_header["P4 PSF ANALYTIC FACTOR COUNT"]),
        "realized_batch_size": int(diagnostic_header["P4 PSF ANALYTIC BATCH SIZE"]),
        "maximum_source_relative_error": max(row["max_source_relative_error"] for row in rows),
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
