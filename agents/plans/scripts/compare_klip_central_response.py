#!/usr/bin/env python3
"""Analyze paired central-difference KLIP response measurements."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

from compare_klip_finite_response import (
    evaluate_response,
    extract_stamp,
    filter_amplitude,
    mode_counts,
    read_cube,
    response_paths,
    rotate_response,
)


def fraction_tag(fraction: float) -> str:
    """Return the directory-safe amplitude tag written by the runner."""
    return format(float(fraction), ".12g").replace(".", "p")


def finite_json(value: object) -> object:
    """Replace non-finite numeric values before strict JSON serialization."""
    if isinstance(value, dict):
        return {str(key): finite_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite_json(item) for item in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def format_number(value: object) -> str:
    """Format a finite number compactly for the Markdown report."""
    if value is None:
        return "nan"
    number = float(value)
    return f"{number:.7g}" if math.isfinite(number) else "nan"


def shape_metrics(
    measured: np.ndarray,
    reference: np.ndarray,
    validity: np.ndarray,
) -> dict[str, float]:
    """Measure scale, cosine, and residual of one response against another."""
    retained = validity & np.isfinite(measured) & np.isfinite(reference)
    measured_values = measured[retained]
    reference_values = reference[retained]
    measured_energy = float(np.dot(measured_values, measured_values))
    reference_energy = float(np.dot(reference_values, reference_values))
    if measured_energy <= 0 or reference_energy <= 0:
        raise RuntimeError("response comparison has no positive-energy finite support")
    cross_energy = float(np.dot(measured_values, reference_values))
    projection = cross_energy / reference_energy
    cosine = cross_energy / math.sqrt(measured_energy * reference_energy)
    residual = measured_values - projection * reference_values
    return {
        "valid_samples": int(np.count_nonzero(retained)),
        "projection": projection,
        "cosine": cosine,
        "best_scaled_relative_residual": math.sqrt(float(np.dot(residual, residual)) / measured_energy),
        "relative_l2_difference": math.sqrt(
            float(np.dot(measured_values - reference_values, measured_values - reference_values))
            / reference_energy
        ),
    }


def response_coordinates(
    separation: float,
    position_angle: float,
    detector_rows: int,
    detector_columns: int,
) -> dict[str, float | int]:
    """Resolve one configured sky sample to the response lattice coordinate."""
    center_row = 0.5 * (detector_rows - 1)
    center_column = 0.5 * (detector_columns - 1)
    pa_radians = math.radians(position_angle)
    source_row = center_row - separation * math.sin(pa_radians)
    source_column = center_column + separation * math.cos(pa_radians)
    # Match the production response lattice's floor(source + 0.5) anchor;
    # Python's round-to-even differs at exact half-pixel coordinates.
    target_row = int(math.floor(source_row + 0.5))
    target_column = int(math.floor(source_column + 0.5))
    return {
        "source_row": source_row,
        "source_column": source_column,
        "target_row": target_row,
        "target_column": target_column,
        "target_radius": math.hypot(target_row - center_row, target_column - center_column),
        "target_angle": math.atan2(target_row - center_row, target_column - center_column),
    }


def read_resource_summary(experiment_directory: Path) -> dict[str, float | int]:
    """Summarize timing records from all completed perturbation reductions."""
    records: list[dict[str, float]] = []
    for path in sorted((experiment_directory / "runs").glob("**/resource_usage.txt")):
        values: dict[str, float] = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            key, separator, text = line.partition("=")
            if separator:
                values[key] = float(text)
        if "wall_seconds" in values and "maximum_rss_kib" in values:
            records.append(values)
    return {
        "completed_reductions": len(records),
        "total_wall_seconds": sum(record["wall_seconds"] for record in records),
        "mean_wall_seconds": (
            sum(record["wall_seconds"] for record in records) / len(records) if records else math.nan
        ),
        "maximum_rss_kib": max((record["maximum_rss_kib"] for record in records), default=math.nan),
    }


def write_response_cube(
    path: Path,
    data: list[np.ndarray],
    modes: list[int],
    role: str,
) -> None:
    """Write a row-major diagnostic response stack with JSON sidecar mapping."""
    header = fits.Header()
    header["HIERARCH KLIP CENTRAL RESPONSE PRODUCT"] = role
    header["HIERARCH KLIP CENTRAL RESPONSE SCHEMA"] = 1
    header["HIERARCH KLIP CENTRAL RESPONSE ROW COUNT"] = len(data)
    header["NMODES"] = ",".join(str(mode) for mode in modes)
    fits.writeto(path, np.asarray(data, dtype=np.float32), header=header, overwrite=False)


def average(values: list[float]) -> float:
    """Return the arithmetic mean of a nonempty numeric list."""
    return sum(values) / len(values)


def main() -> int:
    """Analyze all paired KLIP perturbations and write comparison products."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_experiment", type=Path)
    parser.add_argument("experiment_directory", type=Path)
    parser.add_argument("--reference-case", default="radial_ld_fixed16_filter")
    parser.add_argument("--finite-response-experiment", type=Path, required=True)
    arguments = parser.parse_args()

    reference_experiment = arguments.reference_experiment.resolve()
    experiment_directory = arguments.experiment_directory.resolve()
    finite_response_experiment = arguments.finite_response_experiment.resolve()
    manifest_path = experiment_directory / "samples.json"
    if not manifest_path.is_file():
        raise RuntimeError(f"sample manifest is missing: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if int(manifest.get("schema", 0)) != 1:
        raise RuntimeError("unsupported central-response sample-manifest schema")
    samples = manifest["samples"]
    fractions = [float(value) for value in manifest["amplitude_fractions"]]
    if not samples or not fractions or fractions != sorted(set(fractions)):
        raise RuntimeError("sample manifest has no ordered samples or amplitude fractions")
    fiducial_contrast = float(manifest["planet"]["contrast"])

    reference_case = reference_experiment / arguments.reference_case
    original_path = reference_case / "finim.fits"
    filtered_path = reference_case / "finim_filtered.fits"
    original, original_header = read_cube(original_path)
    filtered, filtered_header = read_cube(filtered_path)
    original_modes = mode_counts(original_header, original_path)
    if mode_counts(filtered_header, filtered_path) != original_modes or filtered.shape != original.shape:
        raise RuntimeError("reference science and frozen-filter cubes are inconsistent")
    paths = response_paths(reference_case)
    if set(paths) != set(original_modes):
        raise RuntimeError("reference radial responses do not match the final-image KL modes")

    finite_path = finite_response_experiment / "empirical_response.fits"
    finite_cube, finite_header = read_cube(finite_path)
    finite_modes = mode_counts(finite_header, finite_path)
    if finite_modes != original_modes or finite_cube.shape[0] != len(original_modes):
        raise RuntimeError("one-sided finite-response modes do not match the central-response reference")

    detector_columns, detector_rows = original.shape[1:]
    row_products: list[np.ndarray] = []
    frozen_products: list[np.ndarray] = []
    residual_products: list[np.ndarray] = []
    rows: list[dict[str, object]] = []
    response_cache: dict[tuple[str, float, int], tuple[np.ndarray, np.ndarray, dict[str, float | int]]] = {}
    smallest_cache: dict[tuple[str, int], tuple[np.ndarray, np.ndarray]] = {}
    analysis_modes: list[int] | None = None

    for sample in samples:
        label = str(sample["label"])
        separation = float(sample["separation"])
        position_angle = float(sample["position_angle"])
        coordinates = response_coordinates(separation, position_angle, detector_rows, detector_columns)
        target_row = int(coordinates["target_row"])
        target_column = int(coordinates["target_column"])
        for fraction in fractions:
            epsilon = fiducial_contrast * fraction
            run_root = experiment_directory / "runs" / label / f"fraction_{fraction_tag(fraction)}"
            plus_path = run_root / "plus" / "finim.fits"
            minus_path = run_root / "minus" / "finim.fits"
            plus, plus_header = read_cube(plus_path)
            minus, minus_header = read_cube(minus_path)
            plus_modes = mode_counts(plus_header, plus_path)
            minus_modes = mode_counts(minus_header, minus_path)
            if plus_modes != minus_modes:
                raise RuntimeError(f"paired perturbation modes differ for sample {label}")
            if analysis_modes is None:
                analysis_modes = plus_modes
                missing_modes = sorted(set(analysis_modes) - set(original_modes))
                if missing_modes:
                    raise RuntimeError(f"perturbation modes are absent from the reference: {missing_modes}")
            elif plus_modes != analysis_modes:
                raise RuntimeError(f"perturbation mode ordering differs for sample {label}")
            if plus.shape[1:] != original.shape[1:] or minus.shape[1:] != original.shape[1:]:
                raise RuntimeError(f"perturbation image dimensions differ for sample {label}")

            for mode_index, retained_modes in enumerate(original_modes):
                if retained_modes not in plus_modes:
                    continue
                perturbation_index = plus_modes.index(retained_modes)
                finite_index = finite_modes.index(retained_modes)
                frozen, frozen_validity = evaluate_response(
                    paths[retained_modes],
                    float(coordinates["target_radius"]),
                    float(coordinates["target_angle"]),
                )
                stamp_size = frozen.shape[0]
                original_stamp = extract_stamp(
                    original[mode_index], target_row, target_column, stamp_size
                )
                plus_stamp = extract_stamp(
                    plus[perturbation_index], target_row, target_column, stamp_size
                )
                minus_stamp = extract_stamp(
                    minus[perturbation_index], target_row, target_column, stamp_size
                )
                central = (plus_stamp - minus_stamp) / (2 * epsilon)
                positive_secant = (plus_stamp - original_stamp) / epsilon
                negative_secant = (original_stamp - minus_stamp) / epsilon
                validity = (
                    frozen_validity
                    & np.isfinite(original_stamp)
                    & np.isfinite(plus_stamp)
                    & np.isfinite(minus_stamp)
                    & np.isfinite(central)
                )
                central_values = central[validity]
                central_energy = float(np.dot(central_values, central_values))
                if central_energy <= 0:
                    raise RuntimeError(
                        f"central response has no positive energy for {label}, fraction {fraction}, "
                        f"mode {retained_modes}"
                    )
                asymmetry = 0.5 * (positive_secant - negative_secant)
                asymmetry_values = asymmetry[validity]
                secant_asymmetry = math.sqrt(
                    float(np.dot(asymmetry_values, asymmetry_values)) / central_energy
                )
                frozen_metrics = shape_metrics(central, frozen, validity)

                cache_key = (label, retained_modes)
                if fraction == fractions[0]:
                    smallest_cache[cache_key] = (central.copy(), validity.copy())
                smallest, smallest_validity = smallest_cache[cache_key]
                smallest_metrics = shape_metrics(
                    central, smallest, validity & smallest_validity
                )

                finite_metrics = {
                    "projection": math.nan,
                    "cosine": math.nan,
                    "best_scaled_relative_residual": math.nan,
                    "relative_l2_difference": math.nan,
                }
                if str(sample["role"]) == "candidate":
                    finite = finite_cube[finite_index].T
                    if finite.shape != central.shape:
                        raise RuntimeError("finite and central response stamps have different dimensions")
                    finite_metrics = shape_metrics(central, finite, validity & np.isfinite(finite))

                central_amplitude, central_normalization, central_support = filter_amplitude(
                    original_stamp, central, validity
                )
                frozen_amplitude, _, frozen_support = filter_amplitude(
                    original_stamp, frozen, frozen_validity
                )
                production_frozen_amplitude = float(
                    filtered[mode_index, target_column, target_row]
                )
                production_frozen_valid = math.isfinite(production_frozen_amplitude)
                if production_frozen_valid and not math.isclose(
                    frozen_amplitude,
                    production_frozen_amplitude,
                    rel_tol=2e-6,
                    abs_tol=2e-10,
                ):
                    raise RuntimeError(
                        f"frozen response reconstruction failed for {label}, mode {retained_modes}"
                    )

                row_index = len(rows)
                rows.append(
                    {
                        "row_index": row_index,
                        "sample_index": int(sample["index"]),
                        "sample_label": label,
                        "sample_role": str(sample["role"]),
                        "configured_separation": separation,
                        "configured_position_angle": position_angle,
                        "candidate_distance": float(sample["candidate_distance"]),
                        "target_row": target_row,
                        "target_column": target_column,
                        "target_radius": float(coordinates["target_radius"]),
                        "target_angle_radians": float(coordinates["target_angle"]),
                        "amplitude_fraction": fraction,
                        "epsilon": epsilon,
                        "mode_count": retained_modes,
                        "valid_stamp_samples": int(frozen_metrics["valid_samples"]),
                        "central_projection_on_frozen": frozen_metrics["projection"],
                        "central_frozen_cosine": frozen_metrics["cosine"],
                        "central_frozen_best_scaled_residual": frozen_metrics[
                            "best_scaled_relative_residual"
                        ],
                        "central_frozen_relative_l2_difference": frozen_metrics[
                            "relative_l2_difference"
                        ],
                        "central_projection_on_smallest_epsilon": smallest_metrics["projection"],
                        "central_smallest_epsilon_cosine": smallest_metrics["cosine"],
                        "central_smallest_epsilon_best_scaled_residual": smallest_metrics[
                            "best_scaled_relative_residual"
                        ],
                        "central_smallest_epsilon_relative_l2_difference": smallest_metrics[
                            "relative_l2_difference"
                        ],
                        "paired_secant_asymmetry": secant_asymmetry,
                        "central_projection_on_finite_secant": finite_metrics["projection"],
                        "central_finite_secant_cosine": finite_metrics["cosine"],
                        "central_finite_secant_best_scaled_residual": finite_metrics[
                            "best_scaled_relative_residual"
                        ],
                        "central_finite_secant_relative_l2_difference": finite_metrics[
                            "relative_l2_difference"
                        ],
                        "central_filter_amplitude": central_amplitude,
                        "central_filter_normalization": central_normalization,
                        "central_filter_support": central_support,
                        "frozen_filter_amplitude": frozen_amplitude,
                        "production_frozen_filter_amplitude": production_frozen_amplitude,
                        "production_frozen_filter_valid": production_frozen_valid,
                        "frozen_filter_support": frozen_support,
                    }
                )
                response_cache[(label, fraction, retained_modes)] = (
                    central.copy(),
                    validity.copy(),
                    coordinates,
                )
                row_products.append(np.where(validity, central, np.nan).T)
                frozen_products.append(np.where(frozen_validity, frozen, np.nan).T)
                residual_products.append(
                    np.where(
                        validity,
                        central - float(frozen_metrics["projection"]) * frozen,
                        np.nan,
                    ).T
                )

    radial_groups: dict[float, list[dict[str, object]]] = defaultdict(list)
    for sample in samples:
        if str(sample["role"]) == "clear":
            radial_groups[float(sample["separation"])].append(sample)

    if analysis_modes is None:
        raise RuntimeError("no paired perturbation modes were analyzed")

    radial_rows: list[dict[str, object]] = []
    for separation, radial_samples in sorted(radial_groups.items()):
        if len(radial_samples) < 2:
            continue
        for fraction in fractions:
            for retained_modes in analysis_modes:
                canonical_responses: list[np.ndarray] = []
                canonical_validities: list[np.ndarray] = []
                target_radii: list[float] = []
                for sample in radial_samples:
                    response, validity, coordinates = response_cache[
                        (str(sample["label"]), fraction, retained_modes)
                    ]
                    canonical, canonical_validity = rotate_response(
                        response, validity, float(coordinates["target_angle"])
                    )
                    canonical_responses.append(canonical)
                    canonical_validities.append(canonical_validity)
                    target_radii.append(float(coordinates["target_radius"]))
                common_validity = np.logical_and.reduce(canonical_validities)
                stack = np.asarray(canonical_responses)
                mean_response = np.mean(stack, axis=0)
                mean_values = mean_response[common_validity]
                mean_energy = float(np.dot(mean_values, mean_values))
                if mean_energy <= 0:
                    raise RuntimeError(
                        f"radial mean has no energy at radius {separation}, fraction {fraction}, "
                        f"mode {retained_modes}"
                    )
                deviations = [
                    math.sqrt(
                        float(
                            np.dot(
                                (response - mean_response)[common_validity],
                                (response - mean_response)[common_validity],
                            )
                        )
                        / mean_energy
                    )
                    for response in canonical_responses
                ]
                representative_radius = average(target_radii)
                frozen_canonical, frozen_validity = evaluate_response(
                    paths[retained_modes], representative_radius, 0.0
                )
                mean_metrics = shape_metrics(
                    mean_response, frozen_canonical, common_validity & frozen_validity
                )
                radial_rows.append(
                    {
                        "configured_separation": separation,
                        "sample_count": len(radial_samples),
                        "sample_labels": [str(sample["label"]) for sample in radial_samples],
                        "amplitude_fraction": fraction,
                        "epsilon": fiducial_contrast * fraction,
                        "mode_count": retained_modes,
                        "representative_lattice_radius": representative_radius,
                        "lattice_radius_range": [min(target_radii), max(target_radii)],
                        "common_valid_samples": int(np.count_nonzero(common_validity)),
                        "angular_scatter_rms": math.sqrt(
                            average([deviation * deviation for deviation in deviations])
                        ),
                        "maximum_angular_relative_deviation": max(deviations),
                        "radial_mean_projection_on_frozen": mean_metrics["projection"],
                        "radial_mean_frozen_cosine": mean_metrics["cosine"],
                        "radial_mean_frozen_best_scaled_residual": mean_metrics[
                            "best_scaled_relative_residual"
                        ],
                    }
                )

    managed = (
        "klip_central_response.json",
        "klip_central_response.csv",
        "klip_central_response.md",
        "central_response.fits",
        "frozen_response.fits",
        "central_minus_scaled_frozen.fits",
    )
    existing = [name for name in managed if (experiment_directory / name).exists()]
    if existing:
        raise RuntimeError("central-response output files already exist: " + ", ".join(existing))

    resource_summary = read_resource_summary(experiment_directory)
    result = {
        "schema": 1,
        "reference_experiment": str(reference_experiment),
        "reference_case": arguments.reference_case,
        "finite_response_experiment": str(finite_response_experiment),
        "sample_manifest": manifest,
        "mode_counts": analysis_modes,
        "smallest_amplitude_fraction": fractions[0],
        "resource_summary": resource_summary,
        "rows": rows,
        "radial_averages": radial_rows,
    }
    json_path = experiment_directory / managed[0]
    csv_path = experiment_directory / managed[1]
    markdown_path = experiment_directory / managed[2]
    json_path.write_text(
        json.dumps(finite_json(result), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(finite_json(rows))
    write_response_cube(
        experiment_directory / managed[3], row_products, analysis_modes, "PAIRED_CENTRAL_DIFFERENCE"
    )
    write_response_cube(
        experiment_directory / managed[4], frozen_products, analysis_modes, "FROZEN_BASIS_COMPARISON"
    )
    write_response_cube(
        experiment_directory / managed[5],
        residual_products,
        analysis_modes,
        "CENTRAL_MINUS_BEST_SCALED_FROZEN",
    )

    clear_rows = [row for row in rows if row["sample_role"] == "clear"]
    candidate_rows = [row for row in rows if row["sample_role"] == "candidate"]
    with markdown_path.open("w", encoding="utf-8") as stream:
        stream.write("# KLIP paired central-response comparison\n\n")
        stream.write(
            f"Measured `{len(samples)}` positions at `{len(fractions)}` perturbation amplitudes and "
            f"`{len(analysis_modes)}` KL mode counts using `{resource_summary['completed_reductions']}` "
            "complete KLIP reductions. "
        )
        if int(resource_summary["completed_reductions"]) > 0:
            stream.write(
                f"Their summed wall time was `{float(resource_summary['total_wall_seconds']):.3f}` seconds, "
                f"with mean `{float(resource_summary['mean_wall_seconds']):.3f}` seconds per reduction.\n\n"
            )
        else:
            stream.write("\n\n")

        stream.write("## Clear-sample response by perturbation amplitude\n\n")
        stream.write(
            "| epsilon/fiducial | Central projection on frozen | Central/frozen cosine | "
            "Best-scaled residual | Change from smallest epsilon | Paired-secant asymmetry |\n"
        )
        stream.write("|---:|---:|---:|---:|---:|---:|\n")
        for fraction in fractions:
            selected = [row for row in clear_rows if row["amplitude_fraction"] == fraction]
            if not selected:
                continue
            stream.write(
                f"| {format_number(fraction)} | "
                f"{format_number(average([float(row['central_projection_on_frozen']) for row in selected]))} | "
                f"{format_number(average([float(row['central_frozen_cosine']) for row in selected]))} | "
                f"{format_number(average([float(row['central_frozen_best_scaled_residual']) for row in selected]))} | "
                f"{format_number(average([float(row['central_smallest_epsilon_relative_l2_difference']) for row in selected]))} | "
                f"{format_number(average([float(row['paired_secant_asymmetry']) for row in selected]))} |\n"
            )

        stream.write("\n## Candidate response\n\n")
        stream.write(
            "| KL modes | epsilon/fiducial | Central projection on frozen | Central/frozen cosine | "
            "Change from smallest epsilon | Central projection on finite secant | Central/finite cosine |\n"
        )
        stream.write("|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in candidate_rows:
            stream.write(
                f"| {row['mode_count']} | {format_number(row['amplitude_fraction'])} | "
                f"{format_number(row['central_projection_on_frozen'])} | "
                f"{format_number(row['central_frozen_cosine'])} | "
                f"{format_number(row['central_smallest_epsilon_relative_l2_difference'])} | "
                f"{format_number(row['central_projection_on_finite_secant'])} | "
                f"{format_number(row['central_finite_secant_cosine'])} |\n"
            )

        stream.write("\n## Azimuthally averaged clear samples\n\n")
        stream.write(
            "| Radius | Samples | epsilon/fiducial | Mean angular scatter | Worst angular deviation | "
            "Radial-mean projection on frozen | Radial-mean/frozen cosine |\n"
        )
        stream.write("|---:|---:|---:|---:|---:|---:|---:|\n")
        for separation in sorted(radial_groups):
            for fraction in fractions:
                selected = [
                    row
                    for row in radial_rows
                    if row["configured_separation"] == separation
                    and row["amplitude_fraction"] == fraction
                ]
                if not selected:
                    continue
                stream.write(
                    f"| {format_number(separation)} | {selected[0]['sample_count']} | "
                    f"{format_number(fraction)} | "
                    f"{format_number(average([float(row['angular_scatter_rms']) for row in selected]))} | "
                    f"{format_number(max(float(row['maximum_angular_relative_deviation']) for row in selected))} | "
                    f"{format_number(average([float(row['radial_mean_projection_on_frozen']) for row in selected]))} | "
                    f"{format_number(average([float(row['radial_mean_frozen_cosine']) for row in selected]))} |\n"
                )

        stream.write(
            "\nThe smallest-epsilon central response is the local numerical-derivative oracle. Stability "
            "with increasing epsilon establishes the usable linear range. Paired-secant asymmetry measures "
            "the even nonlinear component, while the azimuthal rows test whether clear responses can still be "
            "rotated and averaged into a radial model after basis adaptation.\n"
        )

    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
