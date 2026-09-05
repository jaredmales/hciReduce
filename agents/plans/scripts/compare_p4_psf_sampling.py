#!/usr/bin/env python3
"""Compare maintained sparse P4 radial PSF experiments with a dense case."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np
from astropy.io import fits


MODEL_PATTERN = re.compile(r"p4PSF_model_(\d+)\.fits$")


def completed_product_directory(case_directory: Path) -> Path:
    """Return a completed case product directory or raise a useful error."""
    product_directory = case_directory / "finim_outputs"
    manifest = product_directory / "p4PSF_manifest.fits"
    if not manifest.is_file():
        raise RuntimeError(f"missing completion manifest: {manifest}")
    header = fits.getheader(manifest)
    if int(header.get("P4 PSF COMPLETE", 0)) != 1:
        raise RuntimeError(f"incomplete PSF product set: {manifest}")
    return product_directory


def mode_files(product_directory: Path) -> dict[int, Path]:
    """Map mode indices to response-model paths."""
    paths: dict[int, Path] = {}
    for path in product_directory.glob("p4PSF_model_*.fits"):
        match = MODEL_PATTERN.match(path.name)
        if match:
            paths[int(match.group(1))] = path
    if not paths:
        raise RuntimeError(f"no response models found in {product_directory}")
    return paths


def wall_seconds(case_directory: Path) -> float:
    """Read the externally measured wall time when available."""
    path = case_directory / "resource_usage.txt"
    if not path.is_file():
        return math.nan
    for line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition("=")
        if separator and key == "wall_seconds":
            return float(value)
    return math.nan


def logged_seconds(case_directory: Path, label: str) -> float:
    """Read a named timing value from the p4Reduce timing report."""
    path = case_directory / "run.log"
    if not path.is_file():
        return math.nan
    pattern = re.compile(rf"{re.escape(label)}\s+([0-9.eE+-]+)")
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            return float(match.group(1))
    return math.nan


def safe_ratio(numerator: float, denominator: float) -> float:
    """Return a finite ratio or NaN when its denominator is unusable."""
    if denominator <= 0 or not math.isfinite(denominator):
        return math.nan
    return numerator / denominator


def compare_mode(
    dense_model_path: Path,
    dense_validity_path: Path,
    sparse_model_path: Path,
    sparse_validity_path: Path,
) -> dict[str, float | int]:
    """Calculate response-field fidelity metrics for one P4 output mode."""
    dense = np.asarray(fits.getdata(dense_model_path), dtype=np.float64)
    sparse = np.asarray(fits.getdata(sparse_model_path), dtype=np.float64)
    dense_validity = np.asarray(fits.getdata(dense_validity_path)).reshape(-1) > 0.5
    sparse_validity = np.asarray(fits.getdata(sparse_validity_path)).reshape(-1) > 0.5
    if dense.shape != sparse.shape:
        raise RuntimeError(
            f"model shape mismatch: {dense_model_path} {dense.shape} != {sparse_model_path} {sparse.shape}"
        )
    if dense.ndim != 3 or dense.shape[0] != dense_validity.size or dense.shape[0] != sparse_validity.size:
        raise RuntimeError(f"unexpected model/validity shape for {dense_model_path}")

    both_centers = dense_validity & sparse_validity
    dense_finite = np.isfinite(dense) & dense_validity[:, None, None]
    overlap = np.isfinite(dense) & np.isfinite(sparse) & both_centers[:, None, None]
    dense_values = dense[overlap]
    sparse_values = sparse[overlap]
    differences = sparse_values - dense_values

    dense_energy = float(np.dot(dense_values, dense_values))
    sparse_energy = float(np.dot(sparse_values, sparse_values))
    error_energy = float(np.dot(differences, differences))
    cross_energy = float(np.dot(dense_values, sparse_values))

    per_source_dense_energy = np.where(overlap, dense * dense, 0.0).sum(axis=(1, 2))
    per_source_error_energy = np.where(overlap, (sparse - dense) ** 2, 0.0).sum(axis=(1, 2))
    source_eligible = both_centers & (per_source_dense_energy > 0) & overlap.any(axis=(1, 2))
    source_relative_l2 = np.sqrt(
        per_source_error_energy[source_eligible] / per_source_dense_energy[source_eligible]
    )

    center_row = dense.shape[1] // 2
    center_column = dense.shape[2] // 2
    dense_center = dense[:, center_row, center_column]
    sparse_center = sparse[:, center_row, center_column]
    center_overlap = both_centers & np.isfinite(dense_center) & np.isfinite(sparse_center)
    center_error = sparse_center[center_overlap] - dense_center[center_overlap]
    center_denominator = float(np.dot(dense_center[center_overlap], dense_center[center_overlap]))

    return {
        "source_count": dense.shape[0],
        "center_validity_mismatch_fraction": float(np.mean(dense_validity != sparse_validity)),
        "finite_overlap_fraction": safe_ratio(float(overlap.sum()), float(dense_finite.sum())),
        "relative_l2": math.sqrt(safe_ratio(error_energy, dense_energy)),
        "median_source_relative_l2": float(np.median(source_relative_l2)),
        "p95_source_relative_l2": float(np.quantile(source_relative_l2, 0.95)),
        "cosine_similarity": safe_ratio(cross_energy, math.sqrt(dense_energy * sparse_energy)),
        "best_fit_scale": safe_ratio(cross_energy, dense_energy),
        "central_relative_l2": math.sqrt(
            safe_ratio(float(np.dot(center_error, center_error)), center_denominator)
        ),
    }


def compare_filtered_image(dense_case: Path, sparse_case: Path) -> dict[str, float]:
    """Compare the optional spatially filtered science images."""
    dense_path = dense_case / "finim_filtered.fits"
    sparse_path = sparse_case / "finim_filtered.fits"
    if not dense_path.is_file() or not sparse_path.is_file():
        return {
            "filtered_overlap_fraction": math.nan,
            "filtered_relative_l2": math.nan,
            "filtered_cosine_similarity": math.nan,
        }
    dense = np.asarray(fits.getdata(dense_path), dtype=np.float64)
    sparse = np.asarray(fits.getdata(sparse_path), dtype=np.float64)
    if dense.shape != sparse.shape:
        raise RuntimeError(f"filtered image shape mismatch: {dense_path} != {sparse_path}")
    dense_finite = np.isfinite(dense)
    overlap = dense_finite & np.isfinite(sparse)
    dense_values = dense[overlap]
    sparse_values = sparse[overlap]
    differences = sparse_values - dense_values
    dense_energy = float(np.dot(dense_values, dense_values))
    sparse_energy = float(np.dot(sparse_values, sparse_values))
    return {
        "filtered_overlap_fraction": safe_ratio(float(overlap.sum()), float(dense_finite.sum())),
        "filtered_relative_l2": math.sqrt(
            safe_ratio(float(np.dot(differences, differences)), dense_energy)
        ),
        "filtered_cosine_similarity": safe_ratio(
            float(np.dot(dense_values, sparse_values)), math.sqrt(dense_energy * sparse_energy)
        ),
    }


def format_number(value: float | int) -> str:
    """Format a numeric table value compactly while preserving NaN."""
    if isinstance(value, int):
        return str(value)
    if not math.isfinite(value):
        return "nan"
    return f"{value:.8g}"


def write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    """Write dictionaries with a stable shared field order."""
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    """Compare every completed sparse case under an experiment directory."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment_directory", type=Path)
    arguments = parser.parse_args()
    experiment_directory = arguments.experiment_directory.resolve()
    dense_case = experiment_directory / "dense"
    dense_products = completed_product_directory(dense_case)
    dense_modes = mode_files(dense_products)
    dense_coordinates = np.asarray(fits.getdata(dense_products / "p4PSF_coordinates.fits"))
    dense_wall = wall_seconds(dense_case)
    dense_psf_worker = logged_seconds(dense_case, "Local PSF calculation")

    comparison_rows: list[dict[str, float | int | str]] = []
    summary_rows: list[dict[str, float | int | str]] = []
    for sparse_case in sorted(path for path in experiment_directory.iterdir() if path.is_dir() and path != dense_case):
        try:
            sparse_products = completed_product_directory(sparse_case)
        except RuntimeError as error:
            print(f"Skipping {sparse_case.name}: {error}")
            continue
        sparse_modes = mode_files(sparse_products)
        if set(sparse_modes) != set(dense_modes):
            raise RuntimeError(f"mode set mismatch for {sparse_case}")
        sparse_coordinates = np.asarray(fits.getdata(sparse_products / "p4PSF_coordinates.fits"))
        if not np.array_equal(dense_coordinates, sparse_coordinates):
            raise RuntimeError(f"coordinate table mismatch for {sparse_case}")

        sparse_wall = wall_seconds(sparse_case)
        sparse_psf_worker = logged_seconds(sparse_case, "Local PSF calculation")
        sparse_reconstruction = logged_seconds(sparse_case, "PSF field reconstruction/filtering")
        case_mode_rows: list[dict[str, float | int | str]] = []
        for mode_index, dense_model_path in sorted(dense_modes.items()):
            sparse_model_path = sparse_modes[mode_index]
            header = fits.getheader(sparse_model_path)
            metrics = compare_mode(
                dense_model_path,
                dense_products / f"p4PSF_validity_{mode_index:04d}.fits",
                sparse_model_path,
                sparse_products / f"p4PSF_validity_{mode_index:04d}.fits",
            )
            measurements = int(header.get("P4 PSF MEASUREMENT COUNT", metrics["source_count"]))
            row: dict[str, float | int | str] = {
                "case": sparse_case.name,
                "mode_index": mode_index,
                "measurement_count": measurements,
                "measurement_reduction": safe_ratio(int(metrics["source_count"]), measurements),
                "wall_seconds": sparse_wall,
                "speedup_vs_dense": safe_ratio(dense_wall, sparse_wall),
                "psf_worker_seconds": sparse_psf_worker,
                "psf_worker_speedup_vs_dense": safe_ratio(dense_psf_worker, sparse_psf_worker),
                "reconstruction_seconds": sparse_reconstruction,
                **metrics,
            }
            comparison_rows.append(row)
            case_mode_rows.append(row)

        filtered_metrics = compare_filtered_image(dense_case, sparse_case)
        first = case_mode_rows[0]
        relative_l2_values = np.array([float(row["relative_l2"]) for row in case_mode_rows])
        cosine_values = np.array([float(row["cosine_similarity"]) for row in case_mode_rows])
        summary_rows.append(
            {
                "case": sparse_case.name,
                "measurement_count": int(first["measurement_count"]),
                "measurement_reduction": float(first["measurement_reduction"]),
                "wall_seconds": sparse_wall,
                "speedup_vs_dense": float(first["speedup_vs_dense"]),
                "psf_worker_seconds": sparse_psf_worker,
                "psf_worker_speedup_vs_dense": safe_ratio(dense_psf_worker, sparse_psf_worker),
                "reconstruction_seconds": sparse_reconstruction,
                "median_mode_relative_l2": float(np.median(relative_l2_values)),
                "worst_mode_relative_l2": float(np.max(relative_l2_values)),
                "mean_mode_cosine_similarity": float(np.mean(cosine_values)),
                **filtered_metrics,
            }
        )

    if not comparison_rows:
        raise RuntimeError(f"no completed sparse cases found in {experiment_directory}")

    comparison_path = experiment_directory / "response_comparison.csv"
    summary_path = experiment_directory / "comparison_summary.csv"
    report_path = experiment_directory / "comparison_summary.md"
    write_csv(comparison_path, comparison_rows)
    write_csv(summary_path, summary_rows)

    with report_path.open("w", encoding="utf-8") as stream:
        stream.write("# P4 sparse radial PSF sampling comparison\n\n")
        stream.write(f"Dense reference wall time: {format_number(dense_wall)} seconds.\n\n")
        stream.write(f"Dense local-PSF worker time: {format_number(dense_psf_worker)} seconds.\n\n")
        stream.write(
            "| Case | Measurements | Reduction | Wall (s) | Wall speedup | PSF worker (s) | PSF speedup | "
            "Median mode rel. L2 | Worst mode rel. L2 | Mean cosine | Filtered rel. L2 |\n"
        )
        stream.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in summary_rows:
            stream.write(
                f"| {row['case']} | {row['measurement_count']} | "
                f"{format_number(float(row['measurement_reduction']))}x | "
                f"{format_number(float(row['wall_seconds']))} | "
                f"{format_number(float(row['speedup_vs_dense']))}x | "
                f"{format_number(float(row['psf_worker_seconds']))} | "
                f"{format_number(float(row['psf_worker_speedup_vs_dense']))}x | "
                f"{format_number(float(row['median_mode_relative_l2']))} | "
                f"{format_number(float(row['worst_mode_relative_l2']))} | "
                f"{format_number(float(row['mean_mode_cosine_similarity']))} | "
                f"{format_number(float(row['filtered_relative_l2']))} |\n"
            )
        stream.write("\nPer-mode metrics are in `response_comparison.csv`.\n")

    print(f"Wrote {comparison_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
