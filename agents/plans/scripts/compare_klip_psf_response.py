#!/usr/bin/env python3
"""Compare maintained sparse KLIP radial-response experiments with a fine radial reference."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import numpy as np
from astropy.io import fits


RESPONSE_PATTERN = re.compile(r"klipPSF_mode(\d+)_radial_response\.fits$")


def safe_ratio(numerator: float, denominator: float) -> float:
    """Return a finite ratio or NaN when its denominator is unusable."""
    if denominator <= 0 or not math.isfinite(denominator):
        return math.nan
    return numerator / denominator


def case_complete(case_directory: Path) -> bool:
    """Return whether the runner marked a case complete and its final image remains present."""
    return (case_directory / "complete").is_file() and (case_directory / "finim.fits").is_file()


def response_mode_files(case_directory: Path) -> dict[int, Path]:
    """Map KLIP output-mode indices to response-product paths."""
    product_directory = case_directory / "finim_outputs"
    paths: dict[int, Path] = {}
    for path in product_directory.glob("klipPSF_mode*_radial_response.fits"):
        match = RESPONSE_PATTERN.match(path.name)
        if match:
            paths[int(match.group(1))] = path
    if not paths:
        raise RuntimeError(f"no KLIP response products found in {product_directory}")
    return paths


def resource_value(case_directory: Path, key: str) -> float:
    """Read a numeric value from the runner's external resource report."""
    path = case_directory / "resource_usage.txt"
    if not path.is_file():
        return math.nan
    for line in path.read_text(encoding="utf-8").splitlines():
        name, separator, value = line.partition("=")
        if separator and name == key:
            return float(value)
    return math.nan


def logged_seconds(case_directory: Path, label: str) -> float:
    """Read one timing value from the klipReduce timing report."""
    path = case_directory / "run.log"
    if not path.is_file():
        return math.nan
    pattern = re.compile(rf"{re.escape(label)}\s+([0-9.eE+-]+)")
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            return float(match.group(1))
    return math.nan


def header_radii(header: fits.Header, path: Path) -> np.ndarray:
    """Read the ordered KLIP radial sample coordinates from one response header."""
    text = str(header.get("KLIP PSF SAMPLE RADII", ""))
    try:
        radii = np.asarray([float(token) for token in text.split(",") if token.strip()], dtype=np.float64)
    except ValueError as error:
        raise RuntimeError(f"invalid KLIP response radii in {path}") from error
    if radii.size == 0 or not np.isfinite(radii).all() or np.any(np.diff(radii) <= 0):
        raise RuntimeError(f"missing or unordered KLIP response radii in {path}")
    return radii


def load_response(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, fits.Header]:
    """Load one canonical radial response cube and its paired validity cube."""
    header = fits.getheader(path)
    response = np.asarray(fits.getdata(path), dtype=np.float64)
    validity_path = path.with_name(path.name.replace("_radial_response.fits", "_radial_validity.fits"))
    if not validity_path.is_file():
        raise RuntimeError(f"missing KLIP response validity product: {validity_path}")
    validity = np.asarray(fits.getdata(validity_path)) > 0.5
    radii = header_radii(header, path)
    if response.ndim != 3 or response.shape != validity.shape or response.shape[0] != radii.size:
        raise RuntimeError(f"unexpected response/validity geometry in {path}")
    if int(header.get("KLIP PSF PRODUCT SCHEMA", 0)) != 1:
        raise RuntimeError(f"unsupported KLIP response schema in {path}")
    return response, validity, radii, header


def compare_mode(reference_path: Path, candidate_path: Path) -> dict[str, float | int | str]:
    """Compare one sparse radial response with nearest-radius evaluation of the fine reference."""
    reference, reference_validity, reference_radii, reference_header = load_response(reference_path)
    candidate, candidate_validity, candidate_radii, candidate_header = load_response(candidate_path)
    if reference.shape[1:] != candidate.shape[1:]:
        raise RuntimeError(f"response stamp mismatch: {reference_path} != {candidate_path}")
    if int(reference_header.get("KLIP PSF MODE COUNT", -1)) != int(
        candidate_header.get("KLIP PSF MODE COUNT", -2)
    ):
        raise RuntimeError(f"KLIP mode-count mismatch: {reference_path} != {candidate_path}")

    nearest_indices = np.abs(candidate_radii[:, None] - reference_radii[None, :]).argmin(axis=0)
    evaluated = candidate[nearest_indices]
    evaluated_validity = candidate_validity[nearest_indices]
    overlap = reference_validity & evaluated_validity & np.isfinite(reference) & np.isfinite(evaluated)
    reference_finite = reference_validity & np.isfinite(reference)
    reference_values = reference[overlap]
    evaluated_values = evaluated[overlap]
    differences = evaluated_values - reference_values

    reference_energy = float(np.dot(reference_values, reference_values))
    candidate_energy = float(np.dot(evaluated_values, evaluated_values))
    error_energy = float(np.dot(differences, differences))
    cross_energy = float(np.dot(reference_values, evaluated_values))

    per_radius_reference_energy = np.where(overlap, reference * reference, 0.0).sum(axis=(1, 2))
    per_radius_error_energy = np.where(overlap, (evaluated - reference) ** 2, 0.0).sum(axis=(1, 2))
    eligible = (per_radius_reference_energy > 0) & overlap.any(axis=(1, 2))
    per_radius_relative_l2 = np.sqrt(per_radius_error_energy[eligible] / per_radius_reference_energy[eligible])
    radius_errors = np.abs(candidate_radii[nearest_indices] - reference_radii)

    return {
        "mode_count": int(candidate_header.get("KLIP PSF MODE COUNT", -1)),
        "science_combination": str(candidate_header.get("KLIP SCIENCE COMBINATION", "unknown")).strip(),
        "response_combination": str(candidate_header.get("KLIP PSF COMBINATION", "unknown")).strip(),
        "accumulation": str(candidate_header.get("KLIP PSF ACCUMULATION", "unknown")).strip(),
        "measurement_count": int(candidate_header.get("KLIP PSF MEASUREMENT COUNT", 0)),
        "retained_bytes": int(candidate_header.get("KLIP PSF RETAINED BYTES", 0)),
        "radius_count": int(candidate_radii.size),
        "max_nearest_radius_error": float(radius_errors.max()),
        "mean_nearest_radius_error": float(radius_errors.mean()),
        "validity_mismatch_fraction": float(np.mean(reference_validity != evaluated_validity)),
        "finite_overlap_fraction": safe_ratio(float(overlap.sum()), float(reference_finite.sum())),
        "relative_l2": math.sqrt(safe_ratio(error_energy, reference_energy)),
        "median_radius_relative_l2": (
            float(np.median(per_radius_relative_l2)) if per_radius_relative_l2.size else math.nan
        ),
        "p95_radius_relative_l2": (
            float(np.quantile(per_radius_relative_l2, 0.95)) if per_radius_relative_l2.size else math.nan
        ),
        "cosine_similarity": safe_ratio(cross_energy, math.sqrt(reference_energy * candidate_energy)),
        "unit_signal_amplitude_ratio": safe_ratio(cross_energy, candidate_energy),
    }


def compare_science(reference_path: Path, candidate_path: Path) -> dict[str, float]:
    """Compare final science cubes to verify response measurement does not alter the reduction."""
    if not reference_path.is_file() or not candidate_path.is_file():
        return {"science_relative_l2": math.nan, "science_max_abs_error": math.nan}
    reference = np.asarray(fits.getdata(reference_path), dtype=np.float64)
    candidate = np.asarray(fits.getdata(candidate_path), dtype=np.float64)
    if reference.shape != candidate.shape:
        raise RuntimeError(f"science-product shape mismatch: {reference_path} != {candidate_path}")
    overlap = np.isfinite(reference) & np.isfinite(candidate)
    differences = candidate[overlap] - reference[overlap]
    reference_values = reference[overlap]
    return {
        "science_relative_l2": math.sqrt(
            safe_ratio(float(np.dot(differences, differences)), float(np.dot(reference_values, reference_values)))
        ),
        "science_max_abs_error": float(np.max(np.abs(differences))) if differences.size else math.nan,
    }


def product_bytes(case_directory: Path) -> int:
    """Return the stored byte size of radial response and validity products."""
    product_directory = case_directory / "finim_outputs"
    return sum(path.stat().st_size for path in product_directory.glob("klipPSF_mode*_radial_*.fits"))


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
    reference_case = experiment_directory / "reference_dr1_a16"
    if not case_complete(reference_case):
        raise RuntimeError(f"fine response reference is incomplete: {reference_case}")
    reference_modes = response_mode_files(reference_case)
    reference_first_path = reference_modes[min(reference_modes)]
    reference_measurement_count = int(load_response(reference_first_path)[3].get("KLIP PSF MEASUREMENT COUNT", 0))
    reference_wall = resource_value(reference_case, "wall_seconds")
    reference_projection = logged_seconds(reference_case, "PSF calc/sub")

    science_case = experiment_directory / "science_only"
    science_wall = resource_value(science_case, "wall_seconds") if case_complete(science_case) else math.nan
    science_projection = logged_seconds(science_case, "PSF calc/sub") if case_complete(science_case) else math.nan
    science_rss = resource_value(science_case, "maximum_rss_kib") if case_complete(science_case) else math.nan

    comparison_rows: list[dict[str, float | int | str]] = []
    summary_rows: list[dict[str, float | int | str]] = []
    excluded = {reference_case, science_case}
    candidate_cases = sorted(
        path for path in experiment_directory.iterdir() if path.is_dir() and path not in excluded
    )
    for candidate_case in candidate_cases:
        if not case_complete(candidate_case):
            print(f"Skipping incomplete case: {candidate_case.name}")
            continue
        candidate_modes = response_mode_files(candidate_case)
        if set(candidate_modes) != set(reference_modes):
            raise RuntimeError(f"mode set mismatch for {candidate_case}")

        candidate_wall = resource_value(candidate_case, "wall_seconds")
        candidate_projection = logged_seconds(candidate_case, "PSF calc/sub")
        candidate_rss = resource_value(candidate_case, "maximum_rss_kib")
        science_metrics = compare_science(science_case / "finim.fits", candidate_case / "finim.fits")
        mode_rows: list[dict[str, float | int | str]] = []
        for mode_index, reference_path in sorted(reference_modes.items()):
            metrics = compare_mode(reference_path, candidate_modes[mode_index])
            row: dict[str, float | int | str] = {
                "case": candidate_case.name,
                "mode_index": mode_index,
                "wall_seconds": candidate_wall,
                "wall_overhead_vs_science": candidate_wall - science_wall,
                "response_overhead_speedup_vs_reference": safe_ratio(
                    reference_wall - science_wall, candidate_wall - science_wall
                ),
                "projection_worker_seconds": candidate_projection,
                "projection_overhead_vs_science": candidate_projection - science_projection,
                "maximum_rss_kib": candidate_rss,
                "rss_delta_vs_science_kib": candidate_rss - science_rss,
                "product_bytes": product_bytes(candidate_case),
                **science_metrics,
                **metrics,
            }
            comparison_rows.append(row)
            mode_rows.append(row)

        relative_l2 = np.asarray([float(row["relative_l2"]) for row in mode_rows])
        cosine = np.asarray([float(row["cosine_similarity"]) for row in mode_rows])
        amplitude = np.asarray([float(row["unit_signal_amplitude_ratio"]) for row in mode_rows])
        first = mode_rows[0]
        summary_rows.append(
            {
                "case": candidate_case.name,
                "measurement_count": int(first["measurement_count"]),
                "measurement_reduction": safe_ratio(reference_measurement_count, int(first["measurement_count"])),
                "radius_count": int(first["radius_count"]),
                "wall_seconds": candidate_wall,
                "wall_overhead_vs_science": float(first["wall_overhead_vs_science"]),
                "response_overhead_speedup_vs_reference": float(first["response_overhead_speedup_vs_reference"]),
                "projection_overhead_vs_science": float(first["projection_overhead_vs_science"]),
                "maximum_rss_kib": candidate_rss,
                "rss_delta_vs_science_kib": float(first["rss_delta_vs_science_kib"]),
                "retained_bytes": int(first["retained_bytes"]),
                "product_bytes": int(first["product_bytes"]),
                "science_relative_l2": float(first["science_relative_l2"]),
                "science_max_abs_error": float(first["science_max_abs_error"]),
                "median_mode_relative_l2": float(np.median(relative_l2)),
                "worst_mode_relative_l2": float(np.max(relative_l2)),
                "mean_mode_cosine_similarity": float(np.mean(cosine)),
                "worst_unit_signal_amplitude_error": float(np.max(np.abs(amplitude - 1.0))),
            }
        )

    if not comparison_rows:
        raise RuntimeError(f"no completed sparse cases found in {experiment_directory}")

    comparison_path = experiment_directory / "klip_response_comparison.csv"
    summary_path = experiment_directory / "klip_response_summary.csv"
    report_path = experiment_directory / "klip_response_summary.md"
    write_csv(comparison_path, comparison_rows)
    write_csv(summary_path, summary_rows)

    with report_path.open("w", encoding="utf-8") as stream:
        stream.write("# KLIP sparse radial PSF-response comparison\n\n")
        stream.write(f"Fine-reference wall time: {format_number(reference_wall)} seconds.\n\n")
        stream.write(f"Science-only wall time: {format_number(science_wall)} seconds.\n\n")
        stream.write(
            "The cosine and unit-signal amplitude metrics are template-level matched-filter proxies. They do not "
            "replace the later negative-fit/zero-signal-injection inference comparison.\n\n"
        )
        stream.write(
            "| Case | Measurements | Reduction | Wall overhead (s) | Response speedup | RSS delta (KiB) | "
            "Retained (bytes) | Science rel. L2 | Median response rel. L2 | Worst response rel. L2 | "
            "Mean cosine | Worst amplitude error |\n"
        )
        stream.write("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in summary_rows:
            stream.write(
                f"| {row['case']} | {row['measurement_count']} | "
                f"{format_number(float(row['measurement_reduction']))}x | "
                f"{format_number(float(row['wall_overhead_vs_science']))} | "
                f"{format_number(float(row['response_overhead_speedup_vs_reference']))}x | "
                f"{format_number(float(row['rss_delta_vs_science_kib']))} | {row['retained_bytes']} | "
                f"{format_number(float(row['science_relative_l2']))} | "
                f"{format_number(float(row['median_mode_relative_l2']))} | "
                f"{format_number(float(row['worst_mode_relative_l2']))} | "
                f"{format_number(float(row['mean_mode_cosine_similarity']))} | "
                f"{format_number(float(row['worst_unit_signal_amplitude_error']))} |\n"
            )
        stream.write("\nPer-mode metrics are in `klip_response_comparison.csv`.\n")

    print(f"Wrote {comparison_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
