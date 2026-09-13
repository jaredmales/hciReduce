#!/usr/bin/env python3
"""Compare a finite negative-injection KLIP response with the frozen-basis response."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np
from astropy.io import fits


RESPONSE_PATTERN = re.compile(r"klipPSF_mode(\d+)_radial_response\.fits$")


def read_cube(path: Path) -> tuple[np.ndarray, fits.Header]:
    """Read a two- or three-dimensional FITS image as a mode-major cube."""
    header = fits.getheader(path)
    cube = np.asarray(fits.getdata(path), dtype=np.float64)
    if cube.ndim == 2:
        cube = cube[np.newaxis, :, :]
    if cube.ndim != 3:
        raise RuntimeError(f"expected a two- or three-dimensional FITS image: {path}")
    return cube, header


def mode_counts(header: fits.Header, path: Path) -> list[int]:
    """Read the exact ordered KL mode-count vector from a FITS header."""
    try:
        values = [int(token) for token in str(header["NMODES"]).split(",")]
    except (KeyError, ValueError) as error:
        raise RuntimeError(f"FITS product has no valid NMODES vector: {path}") from error
    if not values or len(set(values)) != len(values) or any(value <= 0 for value in values):
        raise RuntimeError(f"FITS product has invalid or duplicate KL mode counts: {path}")
    return values


def response_paths(case_directory: Path) -> dict[int, Path]:
    """Map retained KL mode counts to canonical radial-response products."""
    paths: dict[int, Path] = {}
    for path in (case_directory / "finim_outputs").glob("klipPSF_mode*_radial_response.fits"):
        if RESPONSE_PATTERN.fullmatch(path.name) is None:
            continue
        header = fits.getheader(path)
        if int(header.get("KLIP PSF PRODUCT SCHEMA", 0)) != 1:
            raise RuntimeError(f"unsupported KLIP response schema: {path}")
        if str(header.get("KLIP PSF PRODUCT", "")).strip() != "RADIAL_RESPONSE":
            raise RuntimeError(f"KLIP response product has the wrong declared role: {path}")
        retained_modes = int(header.get("KLIP PSF MODE COUNT", 0))
        if retained_modes <= 0 or retained_modes in paths:
            raise RuntimeError(f"KLIP response mode count is missing or duplicate: {path}")
        paths[retained_modes] = path
    if not paths:
        raise RuntimeError(f"no canonical KLIP radial-response products found in {case_directory}")
    return paths


def response_radii(header: fits.Header, path: Path) -> np.ndarray:
    """Read and validate the ordered radial nodes from a response product."""
    try:
        radii = np.asarray(
            [float(token) for token in str(header["KLIP PSF SAMPLE RADII"]).split(",")],
            dtype=np.float64,
        )
    except (KeyError, ValueError) as error:
        raise RuntimeError(f"KLIP response has no valid radial-node vector: {path}") from error
    if radii.size == 0 or not np.isfinite(radii).all() or np.any(np.diff(radii) <= 0):
        raise RuntimeError(f"KLIP response radial nodes are empty or unordered: {path}")
    return radii


def cubic_weight(distance: float) -> float:
    """Evaluate mxlib's default negative-half cubic-convolution kernel."""
    distance = abs(distance)
    cubic = -0.5
    if distance <= 1:
        return (cubic + 2) * distance**3 - (cubic + 3) * distance**2 + 1
    if distance < 2:
        return cubic * distance**3 - 5 * cubic * distance**2 + 8 * cubic * distance - 4 * cubic
    return 0.0


def rotate_response(
    response: np.ndarray,
    validity: np.ndarray,
    angle: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce RadialPSFModel's cubic response rotation and validity rule."""
    if response.ndim != 2 or response.shape != validity.shape or not math.isfinite(angle):
        raise RuntimeError("response rotation requires matching two-dimensional arrays and a finite angle")
    rows, columns = response.shape
    center_row = 0.5 * (rows - 1)
    center_column = 0.5 * (columns - 1)
    cosine = math.cos(angle)
    sine = math.sin(angle)
    output = np.zeros_like(response, dtype=np.float64)
    output_validity = np.zeros_like(validity, dtype=bool)
    for output_column in range(columns):
        delta_column = output_column - center_column
        for output_row in range(rows):
            delta_row = output_row - center_row
            input_row = center_row + delta_row * cosine + delta_column * sine
            input_column = center_column - delta_row * sine + delta_column * cosine
            floor_row = math.floor(input_row)
            floor_column = math.floor(input_column)
            footprint_row = floor_row - 1
            footprint_column = floor_column - 1
            row_fraction = input_row - floor_row
            column_fraction = input_column - floor_column
            row_weights = (
                cubic_weight(1 + row_fraction),
                cubic_weight(row_fraction),
                cubic_weight(1 - row_fraction),
                cubic_weight(2 - row_fraction),
            )
            column_weights = (
                cubic_weight(1 + column_fraction),
                cubic_weight(column_fraction),
                cubic_weight(1 - column_fraction),
                cubic_weight(2 - column_fraction),
            )
            value = 0.0
            valid = True
            for column_offset, column_weight in enumerate(column_weights):
                for row_offset, row_weight in enumerate(row_weights):
                    weight = row_weight * column_weight
                    if weight == 0:
                        continue
                    input_sample_row = footprint_row + row_offset
                    input_sample_column = footprint_column + column_offset
                    if (
                        input_sample_row < 0
                        or input_sample_row >= rows
                        or input_sample_column < 0
                        or input_sample_column >= columns
                    ):
                        continue
                    if not validity[input_sample_row, input_sample_column]:
                        valid = False
                        break
                    value += response[input_sample_row, input_sample_column] * weight
                if not valid:
                    break
            if valid and math.isfinite(value):
                output[output_row, output_column] = value
                output_validity[output_row, output_column] = True
    return output, output_validity


def evaluate_response(path: Path, radius: float, angle: float) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate one canonical radial model using production interpolation conventions."""
    header = fits.getheader(path)
    response_cube = np.asarray(fits.getdata(path), dtype=np.float64)
    validity_path = path.with_name(path.name.replace("_radial_response.fits", "_radial_validity.fits"))
    if not validity_path.is_file():
        raise RuntimeError(f"missing paired KLIP response validity product: {validity_path}")
    validity_cube = np.asarray(fits.getdata(validity_path), dtype=np.float64) > 0.5
    radii = response_radii(header, path)
    if response_cube.ndim != 3 or response_cube.shape != validity_cube.shape or response_cube.shape[0] != radii.size:
        raise RuntimeError(f"KLIP radial response and validity dimensions are inconsistent: {path}")

    upper = int(np.searchsorted(radii, radius, side="left"))
    if upper == 0:
        lower = upper = 0
        upper_fraction = 0.0
    elif upper == radii.size:
        lower = upper = radii.size - 1
        upper_fraction = 0.0
    elif radii[upper] == radius:
        lower = upper
        upper_fraction = 0.0
    else:
        lower = upper - 1
        upper_fraction = (radius - radii[lower]) / (radii[upper] - radii[lower])

    # hciReduce FITS images appear in NumPy as (column, row), so transpose into
    # the Eigen (row, column) convention used by RadialPSFModel.
    if lower == upper:
        canonical = response_cube[lower].T
        canonical_validity = validity_cube[lower].T
    else:
        canonical = ((1 - upper_fraction) * response_cube[lower] + upper_fraction * response_cube[upper]).T
        canonical_validity = (validity_cube[lower] & validity_cube[upper]).T
    return rotate_response(canonical, canonical_validity, -angle)


def extract_stamp(image: np.ndarray, source_row: int, source_column: int, size: int) -> np.ndarray:
    """Extract an Eigen-oriented square stamp from a FITS-oriented image."""
    half_width = size // 2
    stamp = image[
        source_column - half_width : source_column + half_width + 1,
        source_row - half_width : source_row + half_width + 1,
    ].T
    if stamp.shape != (size, size):
        raise RuntimeError("comparison response stamp extends outside the final image")
    return stamp


def filter_amplitude(image: np.ndarray, response: np.ndarray, validity: np.ndarray) -> tuple[float, float, float]:
    """Apply the production signed normalized filter to one extracted science stamp."""
    retained = validity & np.isfinite(response) & np.isfinite(image)
    support = float(np.count_nonzero(retained)) / response.size
    response_values = response[retained]
    science_values = image[retained]
    normalization = float(np.dot(response_values, response_values))
    if support <= 0 or normalization <= 0:
        raise RuntimeError("comparison response has no finite positive-energy support")
    amplitude = float(np.dot(response_values, science_values) / normalization)
    return amplitude, normalization, support


def optional_fit(reference_experiment: Path, retained_modes: int) -> dict[str, object] | None:
    """Read the prior response-backed fit for one mode when it is available."""
    path = reference_experiment / f"matched_fit_mode{retained_modes}" / "summary.json"
    if not path.is_file():
        return None
    summary = json.loads(path.read_text(encoding="utf-8"))
    if int(summary.get("mode_count", 0)) != retained_modes or not isinstance(summary.get("fit"), dict):
        raise RuntimeError(f"prior KLIP matched-response summary is inconsistent: {path}")
    return summary["fit"]


def finite_json(value: object) -> object:
    """Replace non-finite NumPy and Python values before strict JSON serialization."""
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
    """Format a finite numeric value or an unavailable diagnostic."""
    if value is None:
        return "nan"
    number = float(value)
    return f"{number:.8g}" if math.isfinite(number) else "nan"


def write_stamp_cube(
    path: Path,
    data: np.ndarray,
    modes: list[int],
    role: str,
    source_row: float,
    source_column: float,
) -> None:
    """Write one mode-major diagnostic stamp cube with compact provenance."""
    header = fits.Header()
    header["HIERARCH KLIP FINITE RESPONSE PRODUCT"] = role
    header["NMODES"] = ",".join(str(mode) for mode in modes)
    header["HIERARCH KLIP FINITE RESPONSE ROW"] = source_row
    header["HIERARCH KLIP FINITE RESPONSE COLUMN"] = source_column
    fits.writeto(path, np.asarray(data, dtype=np.float32), header=header, overwrite=False)


def main() -> int:
    """Compare the frozen response with an end-to-end finite negative injection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_experiment", type=Path)
    parser.add_argument("negative_final_image", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--separation", type=float, required=True)
    parser.add_argument("--pa", type=float, required=True)
    parser.add_argument("--contrast", type=float, required=True)
    parser.add_argument("--reference-case", default="radial_ld_fixed16_filter")
    arguments = parser.parse_args()
    if (
        not math.isfinite(arguments.separation)
        or arguments.separation < 0
        or not math.isfinite(arguments.pa)
        or not math.isfinite(arguments.contrast)
        or arguments.contrast <= 0
    ):
        raise RuntimeError("comparison separation, PA, and positive contrast must be finite and valid")

    reference_experiment = arguments.reference_experiment.resolve()
    reference_case = reference_experiment / arguments.reference_case
    original_path = reference_case / "finim.fits"
    filtered_path = reference_case / "finim_filtered.fits"
    negative_path = arguments.negative_final_image.resolve()
    for path in (original_path, filtered_path, negative_path):
        if not path.is_file():
            raise RuntimeError(f"required KLIP product is missing: {path}")

    original, original_header = read_cube(original_path)
    filtered, filtered_header = read_cube(filtered_path)
    negative, negative_header = read_cube(negative_path)
    original_modes = mode_counts(original_header, original_path)
    filtered_modes = mode_counts(filtered_header, filtered_path)
    negative_modes = mode_counts(negative_header, negative_path)
    if original.shape != filtered.shape or original_modes != filtered_modes:
        raise RuntimeError("reference KLIP science and filtered cubes have inconsistent modes or dimensions")
    if negative.shape[1:] != original.shape[1:]:
        raise RuntimeError("negative-injection and reference KLIP images have different spatial dimensions")
    missing_modes = sorted(set(negative_modes) - set(original_modes))
    if missing_modes:
        raise RuntimeError(f"negative-injection modes are absent from the response reference: {missing_modes}")

    paths = response_paths(reference_case)
    missing_responses = sorted(set(negative_modes) - set(paths))
    if missing_responses:
        raise RuntimeError(f"negative-injection modes lack frozen response products: {missing_responses}")

    detector_columns, detector_rows = original.shape[1:]
    center_row = 0.5 * (detector_rows - 1)
    center_column = 0.5 * (detector_columns - 1)
    pa_radians = math.radians(arguments.pa)
    source_row = center_row - arguments.separation * math.sin(pa_radians)
    source_column = center_column + arguments.separation * math.cos(pa_radians)
    target_row = int(round(source_row))
    target_column = int(round(source_column))
    target_radius = math.hypot(target_row - center_row, target_column - center_column)
    target_angle = math.atan2(target_row - center_row, target_column - center_column)

    rows: list[dict[str, object]] = []
    empirical_cubes: list[np.ndarray] = []
    frozen_cubes: list[np.ndarray] = []
    residual_cubes: list[np.ndarray] = []
    for negative_index, retained_modes in enumerate(negative_modes):
        original_index = original_modes.index(retained_modes)
        frozen, frozen_validity = evaluate_response(paths[retained_modes], target_radius, target_angle)
        stamp_size = frozen.shape[0]
        original_stamp = extract_stamp(original[original_index], target_row, target_column, stamp_size)
        negative_stamp = extract_stamp(negative[negative_index], target_row, target_column, stamp_size)
        empirical = (original_stamp - negative_stamp) / arguments.contrast
        comparison_validity = frozen_validity & np.isfinite(empirical) & np.isfinite(frozen)
        empirical_values = empirical[comparison_validity]
        frozen_values = frozen[comparison_validity]
        empirical_energy = float(np.dot(empirical_values, empirical_values))
        frozen_energy = float(np.dot(frozen_values, frozen_values))
        cross_energy = float(np.dot(empirical_values, frozen_values))
        if empirical_energy <= 0 or frozen_energy <= 0:
            raise RuntimeError(f"mode {retained_modes} has no positive-energy finite response comparison")
        projection_ratio = cross_energy / frozen_energy
        cosine = cross_energy / math.sqrt(empirical_energy * frozen_energy)
        best_scaled_residual = math.sqrt(
            float(np.dot(empirical_values - projection_ratio * frozen_values,
                         empirical_values - projection_ratio * frozen_values))
            / empirical_energy
        )

        original_amplitude, normalization, support = filter_amplitude(
            original_stamp, frozen, frozen_validity
        )
        negative_amplitude, negative_normalization, _ = filter_amplitude(
            negative_stamp, frozen, frozen_validity
        )
        if not math.isclose(normalization, negative_normalization, rel_tol=0, abs_tol=1e-9 * normalization):
            raise RuntimeError("the fixed frozen response changed between science-stamp filter applications")
        production_amplitude = float(filtered[original_index, target_column, target_row])
        reconstruction_error = original_amplitude - production_amplitude
        if not math.isfinite(production_amplitude) or abs(reconstruction_error) > 2e-6 * abs(production_amplitude):
            raise RuntimeError(
                f"mode {retained_modes} Python response evaluation does not reproduce production filtering"
            )

        empirical_original, empirical_normalization, empirical_support = filter_amplitude(
            original_stamp, empirical, comparison_validity
        )
        empirical_negative, _, _ = filter_amplitude(negative_stamp, empirical, comparison_validity)
        prior_fit = optional_fit(reference_experiment, retained_modes)
        prior_contrast = float(prior_fit["contrast"]) if prior_fit is not None else math.nan
        rows.append(
            {
                "mode_count": retained_modes,
                "valid_stamp_samples": int(np.count_nonzero(comparison_validity)),
                "frozen_filter_support": support,
                "empirical_filter_support": empirical_support,
                "cosine_similarity": cosine,
                "finite_response_projection_on_frozen": projection_ratio,
                "best_scaled_relative_residual": best_scaled_residual,
                "frozen_normalization": normalization,
                "empirical_normalization": empirical_normalization,
                "production_frozen_amplitude": production_amplitude,
                "reconstructed_frozen_amplitude": original_amplitude,
                "reconstruction_error": reconstruction_error,
                "signal_cancelled_frozen_amplitude": negative_amplitude,
                "removed_frozen_amplitude": original_amplitude - negative_amplitude,
                "removed_frozen_amplitude_fraction": (original_amplitude - negative_amplitude)
                / arguments.contrast,
                "prior_response_fit_contrast": prior_contrast,
                "prior_response_fit_fraction": prior_contrast / arguments.contrast,
                "empirical_filter_original_amplitude": empirical_original,
                "empirical_filter_signal_cancelled_amplitude": empirical_negative,
                "empirical_filter_removed_amplitude": empirical_original - empirical_negative,
            }
        )
        empirical_cubes.append(empirical.T)
        frozen_cubes.append(frozen.T)
        best_scaled = projection_ratio * frozen
        residual_cubes.append((empirical - best_scaled).T)

    output_directory = arguments.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    managed = (
        "klip_finite_response_comparison.csv",
        "klip_finite_response_comparison.json",
        "klip_finite_response_comparison.md",
        "empirical_response.fits",
        "frozen_response.fits",
        "empirical_minus_scaled_frozen.fits",
    )
    existing = [name for name in managed if (output_directory / name).exists()]
    if existing:
        raise RuntimeError("finite-response output files already exist: " + ", ".join(existing))

    result = {
        "schema": 1,
        "reference_experiment": str(reference_experiment),
        "reference_case": arguments.reference_case,
        "original_final_image": str(original_path),
        "negative_final_image": str(negative_path),
        "negative_injection": {
            "separation": arguments.separation,
            "position_angle": arguments.pa % 360,
            "contrast": arguments.contrast,
            "source_row": source_row,
            "source_column": source_column,
        },
        "comparison_coordinate": {
            "row": target_row,
            "column": target_column,
            "radius": target_radius,
            "angle_radians": target_angle,
        },
        "modes": rows,
    }
    json_path = output_directory / managed[1]
    csv_path = output_directory / managed[0]
    markdown_path = output_directory / managed[2]
    json_path.write_text(
        json.dumps(finite_json(result), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    write_stamp_cube(
        output_directory / managed[3],
        np.asarray(empirical_cubes),
        negative_modes,
        "END_TO_END_REMOVED_SIGNAL_PER_CONTRAST",
        source_row,
        source_column,
    )
    write_stamp_cube(
        output_directory / managed[4],
        np.asarray(frozen_cubes),
        negative_modes,
        "FROZEN_BASIS_RESPONSE",
        target_row,
        target_column,
    )
    write_stamp_cube(
        output_directory / managed[5],
        np.asarray(residual_cubes),
        negative_modes,
        "EMPIRICAL_MINUS_BEST_SCALED_FROZEN",
        target_row,
        target_column,
    )

    with markdown_path.open("w", encoding="utf-8") as stream:
        stream.write("# KLIP finite-response comparison\n\n")
        stream.write(
            f"A companion with separation `{arguments.separation:.10g}` pixels, PA "
            f"`{arguments.pa % 360:.10g}` degrees, and contrast `{arguments.contrast:.10g}` was subtracted "
            "before a complete KLIP refit. The empirical response is the original-minus-subtracted final image "
            "divided by that contrast.\n\n"
        )
        stream.write(
            "| KL modes | Cosine | Empirical projection on frozen | Best-scaled residual | "
            "Frozen amplitude before | Frozen amplitude after | Prior fitted contrast/fiducial |\n"
        )
        stream.write("|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in rows:
            stream.write(
                f"| {row['mode_count']} | {format_number(row['cosine_similarity'])} | "
                f"{format_number(row['finite_response_projection_on_frozen'])} | "
                f"{format_number(row['best_scaled_relative_residual'])} | "
                f"{format_number(row['production_frozen_amplitude'])} | "
                f"{format_number(row['signal_cancelled_frozen_amplitude'])} | "
                f"{format_number(row['prior_response_fit_fraction'])} |\n"
            )
        stream.write(
            "\nIf the projection factors reproduce the prior fitted-contrast fractions near 0.245, the "
            "end-to-end response is missing from the frozen-basis model and radial sparsity is not the cause.\n"
        )

    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
