#!/usr/bin/env python3
"""Compare captured P4 detector derivatives and final stamps across PCA precision policies.

The independent reference is the FP64 SVD projector derivative of baseline inputs
in the source direction measured by a unit-amplitude capture. Source interpolation
and image storage remain FP32 in both production-layout reductions.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits

sys.dont_write_bytecode = True
from run_p4_response_convergence import fingerprint, ratio, read_stamp


def read_capture(path: Path) -> dict:
    """Read checked column-major FP64 inputs and pre-storage residuals from a capture pair."""
    metadata = json.loads(path.read_text())
    if metadata["schema"] != 1 or metadata["order"] != "F" or metadata["dtype"] not in ("<f8", ">f8"):
        raise ValueError(f"unsupported detector capture: {path}")
    rows, columns, modes = metadata["rows"], metadata["columns"], metadata["modes"]
    if min(rows, columns) <= 0 or not modes or modes != sorted(set(modes)):
        raise ValueError(f"invalid dimensions or modes: {path}")
    if modes[0] <= 0 or modes[-1] > min(rows, columns):
        raise ValueError(f"retained count exceeds structural rank: {path}")
    if len(metadata["supported"]) != len(modes):
        raise ValueError(f"invalid support flags: {path}")
    values = np.fromfile(path.with_suffix(".bin"), dtype=metadata["dtype"])
    input_count = rows * (columns + 1) if metadata["inputs"] else 0
    if values.size != input_count + rows * len(modes):
        raise ValueError(f"truncated or oversized capture: {path}")
    result = {"metadata": metadata, "residual": values[input_count:].reshape((rows, len(modes)), order="F")}
    if metadata["inputs"]:
        result["predictors"] = values[:rows * columns].reshape((rows, columns), order="F")
        result["target"] = values[rows * columns:input_count]
        if not np.all(np.isfinite(values[:input_count])):
            raise ValueError(f"nonfinite captured inputs: {path}")
    if not np.all(np.isfinite(result["residual"][:, metadata["supported"]])):
        raise ValueError(f"nonfinite supported residual: {path}")
    return result


def same_geometry(left: dict, right: dict) -> None:
    """Reject mismatched coordinates, temporal dimensions, retained counts, or rank thresholds."""
    for key in ("row", "column", "rows", "columns", "modes", "rank_tolerance"):
        if left["metadata"][key] != right["metadata"][key]:
            raise ValueError(f"detector capture mismatch in {key}")


def reference_derivative(baseline: dict, unit: dict) -> tuple:
    """Differentiate the complete temporal projector using an independent SVD and sampled source direction."""
    same_geometry(baseline, unit)
    predictors, target = baseline["predictors"], baseline["target"]
    source = unit["predictors"] - predictors
    target_source = unit["target"] - target
    rows, columns = predictors.shape
    basis, singular, _ = np.linalg.svd(predictors, full_matrices=rows > columns)
    eigenvalues = np.zeros(rows)
    eigenvalues[:singular.size] = singular**2
    gram_derivative = source @ predictors.T + predictors @ source.T
    transformed = basis.T @ gram_derivative @ basis
    derivatives, diagnostics = [], []
    for index, modes in enumerate(baseline["metadata"]["modes"]):
        retained, discarded = basis[:, :modes], basis[:, modes:]
        gap = eigenvalues[modes - 1] - eigenvalues[modes] if modes < rows else None
        resolved = baseline["metadata"]["supported"][index] and (
            gap is None or gap > 64 * np.finfo(float).eps * max(rows, columns) * eigenvalues[0]
        )
        svd_residual = target - retained @ (retained.T @ target)
        diagnostics.append({
            "retained_modes": modes,
            "cutoff_gap_over_leading": None if gap is None else ratio(gap, eigenvalues[0]),
            "cutoff_gap_over_retained": None if gap is None else ratio(gap, eigenvalues[modes - 1]),
            "baseline_double_to_svd_relative_error": ratio(
                np.linalg.norm(baseline["residual"][:, index] - svd_residual), np.linalg.norm(target)
            ) if baseline["metadata"]["supported"][index] else None,
            "reference_resolved": bool(resolved),
            "baseline_rank": baseline["metadata"]["rank"],
        })
        if not resolved:
            derivatives.append(np.full(rows, np.nan))
            continue
        coupling = transformed[modes:, :modes] / (eigenvalues[:modes][None, :] - eigenvalues[modes:, None])
        derivative = target_source - retained @ (retained.T @ target_source)
        derivative -= discarded @ (coupling @ (retained.T @ target))
        derivative -= retained @ (coupling.T @ (discarded.T @ target))
        derivatives.append(derivative)
    return np.column_stack(derivatives), diagnostics


def write_csv(path: Path, rows: list[dict]) -> None:
    """Persist a nonempty comparison table with explicit column names."""
    if not rows:
        raise ValueError(f"empty comparison table: {path}")
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def capture_paths(trial: Path) -> dict:
    """Require the number of unique captured fits advertised by the local product."""
    paths = {path.stem: path for path in (trial / "detectors").glob("fit_*.json")}
    header = fits.getheader(trial / "finim.fits")
    attempted = sum(int(value) for value in str(header["P4 VALID FIT COUNT"]).split(","))
    if not paths or len(paths) != attempted:
        raise ValueError(f"incomplete detector capture: {trial}")
    return paths


def main() -> None:
    """Analyze two complete sweeps with matching data, geometry, trial positions, and amplitudes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--double", required=True, type=Path)
    parser.add_argument("--mixed", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="new analysis directory")
    args = parser.parse_args()
    roots = [args.double.resolve(), args.mixed.resolve()]
    manifests = []
    for root, precision in zip(roots, ("P4-D64", "P4-M32D64")):
        if not (root / "complete.json").exists():
            parser.error(f"incomplete sweep: {root}")
        manifest = json.loads((root / "manifest.json").read_text())
        if manifest["arguments"].get("precision") != precision or not manifest["arguments"].get("capture_detectors"):
            parser.error(f"incorrect precision or missing captures: {root}")
        manifests.append(manifest)
    for key in ("config", "psf"):
        if manifests[0][key]["sha256"] != manifests[1][key]["sha256"]:
            parser.error(f"mismatched {key}")
    if [item["sha256"] for item in manifests[0]["inputs"]] != [item["sha256"] for item in manifests[1]["inputs"]]:
        parser.error("different input sequences")
    settings = manifests[0]["arguments"]
    for key in ("radii", "angles", "modes", "amplitudes", "stamp_size"):
        if settings[key] != manifests[1]["arguments"][key]:
            parser.error(f"mismatched {key}")
    args.output.mkdir(parents=True, exist_ok=False)
    provenance = {"analysis": fingerprint(Path(__file__).resolve()),
                  "sweeps": [fingerprint(root / "manifest.json") for root in roots]}
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    amplitudes = manifests[0]["amplitudes_descending"]
    detector_rows, gap_rows, summary_rows = [], [], []
    positions = [(radius, angle) for radius in settings["radii"] for angle in settings["angles"]]
    for position_index, (radius, angle) in enumerate(positions):
        directories = [root / f"position_{position_index:03d}" for root in roots]
        baseline_paths = [capture_paths(directory / "baseline") for directory in directories]
        unit_paths = [capture_paths(directory / "unit_source") for directory in directories]
        names = sorted(baseline_paths[0])
        if any(set(paths) != set(names) for paths in baseline_paths + unit_paths):
            raise ValueError("baseline/unit capture coordinates differ")
        references = {}
        for name in names:
            baseline, mixed_baseline = [read_capture(paths[name]) for paths in baseline_paths]
            unit, mixed_unit = [read_capture(paths[name]) for paths in unit_paths]
            for left, right in ((baseline, mixed_baseline), (unit, mixed_unit)):
                same_geometry(left, right)
                for key in ("predictors", "target"):
                    if not np.array_equal(left[key], right[key]):
                        raise ValueError(f"precision runs did not receive identical {key}")
            reference, gaps = reference_derivative(baseline, unit)
            references[name] = (reference, gaps, baseline)
            for gap in gaps:
                gap_rows.append({"position": position_index, "fit": name, **gap})
        for amplitude_index, amplitude in enumerate(amplitudes):
            trial_dirs = [directory / f"amplitude_{amplitude_index:02d}_{sign}"
                          for directory in directories for sign in ("positive", "negative")]
            paths = [capture_paths(directory) for directory in trial_dirs]
            if any(set(item) != set(names) for item in paths):
                raise ValueError("signed capture coordinates differ")
            errors = np.zeros((len(settings["modes"]), 4))
            counts = np.zeros(len(settings["modes"]), dtype=int)
            for name in names:
                captures = [read_capture(item[name]) for item in paths]
                reference, gaps, baseline = references[name]
                for capture in captures:
                    same_geometry(baseline, capture)
                double = (captures[0]["residual"] - captures[1]["residual"]) / (2 * amplitude)
                mixed = (captures[2]["residual"] - captures[3]["residual"]) / (2 * amplitude)
                for index, fraction in enumerate(settings["modes"]):
                    supported = gaps[index]["reference_resolved"] and all(
                        capture["metadata"]["supported"][index] for capture in captures
                    )
                    row = {"position": position_index, "fit": name, "amplitude": amplitude,
                           "mode_fraction": fraction, "common_supported": bool(supported),
                           "double_relative_error": None, "mixed_relative_error": None}
                    if supported:
                        norm = np.linalg.norm(reference[:, index])
                        double_error = np.linalg.norm(double[:, index] - reference[:, index])
                        mixed_error = np.linalg.norm(mixed[:, index] - reference[:, index])
                        difference = np.linalg.norm(mixed[:, index] - double[:, index])
                        errors[index] += np.array([norm, double_error, mixed_error, difference])**2
                        counts[index] += 1
                        row.update(double_relative_error=ratio(double_error, norm),
                                   mixed_relative_error=ratio(mixed_error, norm))
                    detector_rows.append(row)
            stamps = [read_stamp(directory, settings["modes"], settings["stamp_size"],
                                 (radius, angle, amplitude if index % 2 == 0 else -amplitude),
                                 "P4-D64" if index < 2 else "P4-M32D64")
                      for index, directory in enumerate(trial_dirs)]
            if any(stamp[2] != stamps[0][2] for stamp in stamps):
                raise ValueError("final stamp origins differ")
            for index, fraction in enumerate(settings["modes"]):
                support = np.logical_and.reduce([stamp[1][index] for stamp in stamps])
                double = (stamps[0][0][index][support] - stamps[1][0][index][support]) / (2 * amplitude)
                mixed = (stamps[2][0][index][support] - stamps[3][0][index][support]) / (2 * amplitude)
                norms = np.sqrt(errors[index])
                summary_rows.append({
                    "position": position_index, "radius_pixels": radius, "angle_degrees": angle,
                    "mode_fraction": fraction, "amplitude": amplitude,
                    "common_detector_fits": int(counts[index]), "total_detector_fits": len(names),
                    "double_detector_relative_error": ratio(norms[1], norms[0]),
                    "mixed_detector_relative_error": ratio(norms[2], norms[0]),
                    "detector_precision_difference_over_reference": ratio(norms[3], norms[0]),
                    "common_final_pixels": int(support.sum()),
                    "final_support_disagreement_pixels": int(np.count_nonzero(
                        np.logical_or.reduce([stamp[1][index] for stamp in stamps]) != support)),
                    "final_mixed_to_double_relative_difference": ratio(np.linalg.norm(mixed - double),
                                                                      np.linalg.norm(double)),
                })
        print(f"Analyzed position {position_index}: {len(names)} detector fits", flush=True)
    write_csv(args.output / "eigengaps.csv", gap_rows)
    write_csv(args.output / "detectors.csv", detector_rows)
    write_csv(args.output / "summary.csv", summary_rows)
    (args.output / "summary.json").write_text(json.dumps(summary_rows, indent=2, allow_nan=False) + "\n")
    (args.output / "complete.json").write_text(json.dumps({"comparisons": len(summary_rows)}) + "\n")


if __name__ == "__main__":
    main()
