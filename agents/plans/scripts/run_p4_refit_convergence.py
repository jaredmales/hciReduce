#!/usr/bin/env python3
"""Sweep pre-storage refit differences through reconstruction and radial PSF products.

Reuse the input/geometry provenance of a completed run_p4_response_convergence.py
experiment. Comparisons cover the published radial model, including its angular
averaging; they do not designate the smallest amplitude as an exact derivative.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from astropy.io import fits

from run_p4_response_convergence import fingerprint, ratio


def read_products(directory: Path, modes: list, size: int, precision: str, amplitude: float) -> tuple:
    """Validate response provenance, dimensions and explicit per-source validity."""
    products = directory / "finim_outputs"
    coordinates = fits.getdata(products / "response_coordinates.fits")
    responses, validities = [], []
    for mode in range(len(modes)):
        data, header = fits.getdata(products / f"response_model_{mode:04d}.fits", header=True)
        validity = fits.getdata(products / f"response_validity_{mode:04d}.fits").reshape(-1)
        if data.shape != (coordinates.shape[1], size, size) or validity.size != data.shape[0]:
            raise ValueError(f"unexpected response dimensions in {directory}")
        if not str(header["P4 PSF RESPONSE"]).startswith("REFIT_CENTRAL_DIFFERENCE"):
            raise ValueError(f"not a refit-difference response in {directory}")
        if str(header["P4POLCY"]).strip() != precision.removeprefix("P4-"):
            raise ValueError(f"incorrect precision policy in {directory}")
        if not np.isclose(float(header["P4 PSF REFIT CONTRAST"]), amplitude, rtol=2e-6, atol=0):
            raise ValueError(f"incorrect refit amplitude in {directory}")
        actual_modes = [float(value) for value in str(header["P4 MODE FRACTIONS"]).split(",")]
        if len(actual_modes) != len(modes) or not np.allclose(actual_modes, modes, rtol=2e-6, atol=0):
            raise ValueError(f"incorrect mode fractions in {directory}")
        if not np.all((validity == 0) | (validity == 1)):
            raise ValueError(f"invalid response validity in {directory}")
        support = np.isfinite(data) & (validity[:, None, None] != 0)
        responses.append(data.astype(np.float64))
        validities.append(support)
    return np.asarray(responses), np.asarray(validities), coordinates


def main() -> None:
    """Run a fresh amplitude sweep, retaining commands, inputs, products and comparison metrics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", required=True, type=Path, help="completed local response sweep")
    parser.add_argument("--output", required=True, type=Path, help="new output directory")
    parser.add_argument("--binary", type=Path, help="explicitly use and fingerprint a rebuilt benchmark")
    parser.add_argument("--samples-per-radius", type=int, default=4)
    args = parser.parse_args()
    if args.samples_per_radius < 1:
        parser.error("samples per radius must be positive")
    experiment, output = args.experiment.resolve(), args.output.resolve()
    if not (experiment / "complete.json").is_file():
        parser.error("reference experiment is incomplete")
    manifest_path = experiment / "manifest.json"
    previous = json.loads(manifest_path.read_text())
    settings = previous["arguments"]
    precision = settings["precision"]
    if precision not in ("P4-D64", "P4-M32D64"):
        parser.error("reference experiment must use an experimental precision policy")
    binary = fingerprint(args.binary.resolve()) if args.binary else previous["binary"]
    for record in [binary, previous["config"], previous["psf"]] + previous["inputs"]:
        if fingerprint(Path(record["path"])) != record:
            parser.error(f"reference input changed: {record['path']}")
    amplitudes = previous["amplitudes_descending"]
    modes, size = settings["modes"], settings["stamp_size"]
    output.mkdir(parents=True, exist_ok=False)
    (output / "manifest.json").write_text(json.dumps({
        "schema": 1, "reference_experiment": fingerprint(manifest_path),
        "reference_manifest": previous, "binary": binary, "runner": fingerprint(Path(__file__).resolve()),
        "samples_per_radius": args.samples_per_radius,
        "comparison": "adjacent amplitudes of published radial models, including angular averaging",
        "thread_environment": {key: os.environ.get(key) for key in
                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    }, indent=2) + "\n")
    file_list = output / "inputs.txt"
    file_list.write_text("".join(f"{record['path']}\n" for record in previous["inputs"]))
    fixed = {
        "mode": "normal", "precision": precision, "input.directory": "", "input.fileList": file_list,
        "input.deleteFront": 0, "input.deleteBack": 0, "input.qualityFile": "", "input.qualityThreshold": 0,
        "input.imSize": 0, "coadd.method": "none", "preProcess.skip": "true", "preProcess.only": "false",
        "preProcess.outputPrefix": "", "adi.postMedSub": "false", "adi.excludeMethod": "none",
        "solver.exclusionSolver": "explicitRefit", "solver.deletionBackend": "leadingCovariance",
        "p4.regressionFrame": "detector", "p4.temporalPredictor": "none", "p4.numberImages": 0,
        "p4.modeFractions": ",".join(map(str, modes)), "p4.localStampSize": 0, "p4.memoryFraction": 0,
        "p4.writeDiagnostics": "false", "p4Optimize.enabled": "false",
        "psfResponse.file": previous["psf"]["path"], "psfResponse.stampSize": size,
        "psfResponse.method": "refitDifference", "psfResponse.filter": "false",
        "psfResponse.outputPrefix": "response_", "psfResponse.sampleRadii": ",".join(map(str, settings["radii"])),
        "psfResponse.radiiPerRegion": 0, "psfResponse.samplesPerRadius": args.samples_per_radius,
        "psfResponse.sampleArcStep": 0, "psfResponse.sampleAvoidRadius": 0,
        "fake.fileName": "", "combine.method": "mean", "combine.weightFile": "",
        "combine.minGoodFract": 0, "combine.noDerotate": "false",
        "output.fileName": "finim.fits", "output.exactFName": "true", "output.outputPSFSub": "false",
    }
    common = [binary["path"], "--config", previous["config"]["path"]]
    common.extend(f"--{key}={value}" for key, value in fixed.items())
    environment = os.environ.copy()
    environment.pop("P4REDUCE_GLOBAL_CONFIG", None)
    baseline, last, rows = None, None, []
    for index, amplitude in enumerate([None] + amplitudes):
        directory = output / ("baseline" if amplitude is None else f"amplitude_{index - 1:03d}")
        directory.mkdir()
        if amplitude is None:
            overrides = ["--psfResponse.file=", "--psfResponse.sampleRadii=", "--psfResponse.samplesPerRadius=0",
                         "--psfResponse.method=skyExact", "--psfResponse.refitContrast=0",
                         "--psfResponse.outputModels=false"]
        else:
            overrides = [f"--psfResponse.refitContrast={amplitude}", "--psfResponse.outputModels=true"]
        command = common + overrides + [f"--output.directory={directory}"]
        (directory / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        print(f"{precision} {directory.name}: half-amplitude={amplitude}", flush=True)
        start = time.monotonic()
        with (directory / "run.log").open("w") as log:
            subprocess.run(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
        (directory / "timing.json").write_text(json.dumps({"wall_seconds": time.monotonic() - start}) + "\n")
        science = fits.getdata(directory / "finim.fits")
        if baseline is None:
            baseline = science.copy()
            continue
        if not np.array_equal(science, baseline, equal_nan=True):
            raise ValueError(f"enabling response calculation changed science in {directory}")
        current = read_products(directory, modes, size, precision, amplitude)
        if last is not None:
            if not np.array_equal(last[2], current[2]):
                raise ValueError("source coordinate order changed between amplitudes")
            for mode, fraction in enumerate(modes):
                common_support = last[1][mode] & current[1][mode]
                a, b = last[0][mode][common_support], current[0][mode][common_support]
                an, bn = np.linalg.norm(a), np.linalg.norm(b)
                changes = []
                for source in range(current[0].shape[1]):
                    mask = common_support[source]
                    x, y = last[0][mode, source][mask], current[0][mode, source][mask]
                    change = ratio(np.linalg.norm(x - y), np.linalg.norm(y))
                    if change is not None:
                        changes.append(change)
                rows.append({
                    "mode_fraction": fraction, "amplitude_large": amplitudes[index - 2],
                    "amplitude_small": amplitude, "common_pixels": int(common_support.sum()),
                    "compared_sources": len(changes),
                    "support_changed_pixels": int(np.count_nonzero(last[1][mode] != current[1][mode])),
                    "relative_change": ratio(np.linalg.norm(a - b), bn),
                    "max_source_relative_change": max(changes, default=None),
                    "cosine": ratio(np.dot(a, b), an * bn),
                    "projection_scale_large_to_small": ratio(np.dot(a, b), bn**2),
                })
        last = current
    with (output / "summary.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "summary.json").write_text(json.dumps(rows, indent=2) + "\n")
    (output / "complete.json").write_text(json.dumps({
        "reductions": len(amplitudes) + 1, "science_exactly_unchanged": True,
        "comparisons": len(rows),
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
