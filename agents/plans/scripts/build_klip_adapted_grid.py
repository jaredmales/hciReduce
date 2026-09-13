#!/usr/bin/env python3
"""Build and apply a radial KLIP response model from paired-refit grid measurements."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

from compare_klip_central_response import shape_metrics
from compare_klip_finite_response import extract_stamp, filter_amplitude, mode_counts, rotate_response


def read_cube(path: Path) -> tuple[np.ndarray, fits.Header]:
    """Read a two- or three-dimensional FITS image as a mode-major cube."""
    header = fits.getheader(path)
    cube = np.asarray(fits.getdata(path), dtype=np.float64)
    if cube.ndim == 2:
        cube = cube[np.newaxis, :, :]
    if cube.ndim != 3:
        raise RuntimeError(f"expected a two- or three-dimensional FITS image: {path}")
    return cube, header


def format_number(value: object) -> str:
    """Format one finite number compactly for Markdown output."""
    if value is None:
        return "nan"
    number = float(value)
    return f"{number:.8g}" if math.isfinite(number) else "nan"


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


def cubic_weight(distance: np.ndarray) -> np.ndarray:
    """Evaluate the production negative-half cubic-convolution kernel."""
    absolute = np.abs(distance)
    result = np.zeros_like(absolute, dtype=np.float64)
    inner = absolute <= 1
    outer = (absolute > 1) & (absolute < 2)
    result[inner] = 1.5 * absolute[inner] ** 3 - 2.5 * absolute[inner] ** 2 + 1
    result[outer] = (
        -0.5 * absolute[outer] ** 3
        + 2.5 * absolute[outer] ** 2
        - 4 * absolute[outer]
        + 2
    )
    return result


def rotate_response_batch(
    responses: np.ndarray,
    validities: np.ndarray,
    angles: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Rotate a batch of response stamps using production cubic conventions."""
    if (
        responses.ndim != 3
        or responses.shape != validities.shape
        or angles.shape != (responses.shape[0],)
        or responses.shape[1] <= 0
        or responses.shape[2] <= 0
    ):
        raise RuntimeError("batched response rotation has inconsistent dimensions")
    batch, rows, columns = responses.shape
    center_row = 0.5 * (rows - 1)
    center_column = 0.5 * (columns - 1)
    output_rows, output_columns = np.meshgrid(
        np.arange(rows, dtype=np.float64),
        np.arange(columns, dtype=np.float64),
        indexing="ij",
    )
    delta_rows = (output_rows - center_row).reshape(1, -1)
    delta_columns = (output_columns - center_column).reshape(1, -1)
    cosine = np.cos(angles).reshape(-1, 1)
    sine = np.sin(angles).reshape(-1, 1)
    input_rows = center_row + delta_rows * cosine + delta_columns * sine
    input_columns = center_column - delta_rows * sine + delta_columns * cosine
    floor_rows = np.floor(input_rows).astype(np.int64)
    floor_columns = np.floor(input_columns).astype(np.int64)
    footprint_rows = floor_rows - 1
    footprint_columns = floor_columns - 1
    row_fraction = input_rows - floor_rows
    column_fraction = input_columns - floor_columns
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

    values = np.zeros((batch, rows * columns), dtype=np.float64)
    output_validity = np.ones((batch, rows * columns), dtype=bool)
    batch_indices = np.arange(batch).reshape(-1, 1)
    for column_offset, column_weight in enumerate(column_weights):
        for row_offset, row_weight in enumerate(row_weights):
            weight = row_weight * column_weight
            input_row = footprint_rows + row_offset
            input_column = footprint_columns + column_offset
            in_bounds = (
                (input_row >= 0)
                & (input_row < rows)
                & (input_column >= 0)
                & (input_column < columns)
            )
            active = in_bounds & (weight != 0)
            safe_row = np.clip(input_row, 0, rows - 1)
            safe_column = np.clip(input_column, 0, columns - 1)
            sampled_values = responses[batch_indices, safe_row, safe_column]
            sampled_validity = validities[batch_indices, safe_row, safe_column]
            values += np.where(active, sampled_values * weight, 0)
            output_validity &= ~(active & ~sampled_validity)
    output_validity &= np.isfinite(values)
    return values.reshape(batch, rows, columns), output_validity.reshape(batch, rows, columns)


def evaluate_radial_batch(
    responses: np.ndarray,
    validities: np.ndarray,
    radii: np.ndarray,
    requested_radii: np.ndarray,
    requested_angles: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate and rotate a canonical radial model at a batch of source coordinates."""
    if responses.ndim != 3 or responses.shape != validities.shape or responses.shape[0] != radii.size:
        raise RuntimeError("radial response model dimensions are inconsistent")
    upper = np.searchsorted(radii, requested_radii, side="left")
    lower = np.empty_like(upper)
    upper_fraction = np.zeros(requested_radii.shape, dtype=np.float64)
    below = upper == 0
    above = upper == radii.size
    exact_or_below = below | (~above & (radii[np.minimum(upper, radii.size - 1)] == requested_radii))
    lower[below] = 0
    upper[below] = 0
    lower[above] = radii.size - 1
    upper[above] = radii.size - 1
    exact = exact_or_below & ~below
    lower[exact] = upper[exact]
    interpolate = ~(below | above | exact)
    lower[interpolate] = upper[interpolate] - 1
    upper_fraction[interpolate] = (
        requested_radii[interpolate] - radii[lower[interpolate]]
    ) / (radii[upper[interpolate]] - radii[lower[interpolate]])
    lower_fraction = 1 - upper_fraction
    canonical = (
        lower_fraction[:, np.newaxis, np.newaxis] * responses[lower]
        + upper_fraction[:, np.newaxis, np.newaxis] * responses[upper]
    )
    canonical_validity = validities[lower] & validities[upper]
    return rotate_response_batch(canonical, canonical_validity, -requested_angles)


def build_radial_models(
    summary: dict[str, object],
    central_cube: np.ndarray,
    amplitude_fraction: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int], list[int]]:
    """Rotate and average clear central responses at every measured radius."""
    modes = [int(value) for value in summary["mode_counts"]]
    rows = summary["rows"]
    selected_rows = [
        row
        for row in rows
        if row["sample_role"] == "clear"
        and math.isclose(float(row["amplitude_fraction"]), amplitude_fraction, rel_tol=0, abs_tol=1e-15)
    ]
    radii = np.asarray(sorted({float(row["configured_separation"]) for row in selected_rows}))
    if radii.size == 0:
        raise RuntimeError("central-response result contains no clear samples at the requested amplitude")
    stamp_size = central_cube.shape[-1]
    responses = np.zeros((len(modes), radii.size, stamp_size, stamp_size), dtype=np.float64)
    validities = np.zeros_like(responses, dtype=bool)
    sample_counts: list[int] = []
    for radius_index, radius in enumerate(radii):
        radius_rows = [row for row in selected_rows if float(row["configured_separation"]) == radius]
        labels = {str(row["sample_label"]) for row in radius_rows}
        sample_counts.append(len(labels))
        if not labels:
            raise RuntimeError(f"radial node {radius} has no clear response measurements")
        for mode_index, retained_modes in enumerate(modes):
            measured_rows = [row for row in radius_rows if int(row["mode_count"]) == retained_modes]
            if len(measured_rows) != len(labels):
                raise RuntimeError(f"radial node {radius}, mode {retained_modes} has incomplete samples")
            sums = np.zeros((stamp_size, stamp_size), dtype=np.float64)
            counts = np.zeros((stamp_size, stamp_size), dtype=np.int64)
            for row in measured_rows:
                measured = central_cube[int(row["row_index"])].T
                measured_validity = np.isfinite(measured)
                canonical, canonical_validity = rotate_response(
                    np.nan_to_num(measured),
                    measured_validity,
                    float(row["target_angle_radians"]),
                )
                sums[canonical_validity] += canonical[canonical_validity]
                counts[canonical_validity] += 1
            validities[mode_index, radius_index] = counts > 0
            responses[mode_index, radius_index] = np.divide(
                sums,
                counts,
                out=np.zeros_like(sums),
                where=counts > 0,
            )
    return responses, validities, radii, modes, sample_counts


def filter_science_cube(
    science: np.ndarray,
    responses: np.ndarray,
    validities: np.ndarray,
    radii: np.ndarray,
    minimum_radius: float,
    maximum_radius: float,
    minimum_support_fraction: float,
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Apply a radial response model to every configured final-image location."""
    if science.ndim != 3 or responses.shape[0] != science.shape[0]:
        raise RuntimeError("science and response mode dimensions differ")
    mode_count, detector_columns, detector_rows = science.shape
    center_row = 0.5 * (detector_rows - 1)
    center_column = 0.5 * (detector_columns - 1)
    row_grid, column_grid = np.meshgrid(
        np.arange(detector_rows),
        np.arange(detector_columns),
        indexing="ij",
    )
    grid_radius = np.hypot(row_grid - center_row, column_grid - center_column)
    inside = (grid_radius >= minimum_radius) & (grid_radius < maximum_radius)
    source_rows = row_grid[inside].astype(np.int64)
    source_columns = column_grid[inside].astype(np.int64)
    source_radii = grid_radius[inside]
    source_angles = np.arctan2(source_rows - center_row, source_columns - center_column)

    filtered = np.full_like(science, np.nan, dtype=np.float64)
    normalization = np.full_like(science, np.nan, dtype=np.float64)
    support = np.full_like(science, np.nan, dtype=np.float64)
    filter_validity = np.zeros_like(science, dtype=np.float64)
    stamp_size = responses.shape[-1]
    stamp_center = stamp_size // 2
    row_offsets, column_offsets = np.meshgrid(
        np.arange(stamp_size) - stamp_center,
        np.arange(stamp_size) - stamp_center,
        indexing="ij",
    )
    for mode_index in range(mode_count):
        science_image = science[mode_index]
        for begin in range(0, source_rows.size, chunk_size):
            end = min(begin + chunk_size, source_rows.size)
            rows = source_rows[begin:end]
            columns = source_columns[begin:end]
            radial_response, response_validity = evaluate_radial_batch(
                responses[mode_index],
                validities[mode_index],
                radii,
                source_radii[begin:end],
                source_angles[begin:end],
            )
            image_rows = rows[:, np.newaxis, np.newaxis] + row_offsets
            image_columns = columns[:, np.newaxis, np.newaxis] + column_offsets
            in_bounds = (
                (image_rows >= 0)
                & (image_rows < detector_rows)
                & (image_columns >= 0)
                & (image_columns < detector_columns)
            )
            safe_rows = np.clip(image_rows, 0, detector_rows - 1)
            safe_columns = np.clip(image_columns, 0, detector_columns - 1)
            science_values = science_image[safe_columns, safe_rows]
            retained = response_validity & in_bounds & np.isfinite(science_values)
            response_values = np.where(retained, radial_response, 0)
            retained_science = np.where(retained, science_values, 0)
            correlations = np.sum(response_values * retained_science, axis=(1, 2), dtype=np.float64)
            normalizations = np.sum(response_values * response_values, axis=(1, 2), dtype=np.float64)
            supports = np.count_nonzero(retained, axis=(1, 2)) / float(stamp_size * stamp_size)
            center_valid = response_validity[:, stamp_center, stamp_center]
            valid = (
                center_valid
                & (supports >= minimum_support_fraction)
                & (normalizations > 0)
                & np.isfinite(correlations)
                & np.isfinite(normalizations)
            )
            amplitudes = np.divide(
                correlations,
                normalizations,
                out=np.full_like(correlations, np.nan),
                where=valid,
            )
            filtered[mode_index, columns, rows] = amplitudes
            normalization[mode_index, columns, rows] = np.where(
                center_valid & np.isfinite(normalizations) & (normalizations >= 0),
                normalizations,
                np.nan,
            )
            support[mode_index, columns, rows] = np.where(center_valid, supports, np.nan)
            filter_validity[mode_index, columns, rows] = valid.astype(np.float64)
    return filtered, normalization, support, filter_validity


def candidate_comparison(
    candidate_experiment: Path,
    amplitude_fraction: float,
    modes: list[int],
    science: np.ndarray,
    responses: np.ndarray,
    validities: np.ndarray,
    radii: np.ndarray,
) -> list[dict[str, object]]:
    """Compare the clear radial model with an independent candidate response."""
    summary_path = candidate_experiment / "klip_central_response.json"
    cube_path = candidate_experiment / "central_response.fits"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    cube = np.asarray(fits.getdata(cube_path), dtype=np.float64)
    candidate_rows = [
        row
        for row in summary["rows"]
        if row["sample_role"] == "candidate"
        and math.isclose(float(row["amplitude_fraction"]), amplitude_fraction, rel_tol=0, abs_tol=1e-15)
    ]
    if not candidate_rows:
        raise RuntimeError("candidate reference has no response at the requested amplitude fraction")
    results: list[dict[str, object]] = []
    for mode_index, retained_modes in enumerate(modes):
        matches = [row for row in candidate_rows if int(row["mode_count"]) == retained_modes]
        if len(matches) != 1:
            raise RuntimeError(f"candidate reference does not uniquely contain mode {retained_modes}")
        row = matches[0]
        candidate = cube[int(row["row_index"])].T
        candidate_validity = np.isfinite(candidate)
        model, model_validity = evaluate_radial_batch(
            responses[mode_index],
            validities[mode_index],
            radii,
            np.asarray([float(row["target_radius"])]),
            np.asarray([float(row["target_angle_radians"])]),
        )
        model = model[0]
        model_validity = model_validity[0]
        metrics = shape_metrics(candidate, model, candidate_validity & model_validity)
        stamp = extract_stamp(
            science[mode_index],
            int(row["target_row"]),
            int(row["target_column"]),
            model.shape[0],
        )
        adapted_amplitude, normalization, support = filter_amplitude(
            stamp, model, model_validity
        )
        results.append(
            {
                "mode_count": retained_modes,
                "candidate_projection_on_clear_grid": metrics["projection"],
                "candidate_clear_grid_cosine": metrics["cosine"],
                "candidate_clear_grid_best_scaled_residual": metrics[
                    "best_scaled_relative_residual"
                ],
                "adapted_filter_amplitude": adapted_amplitude,
                "adapted_filter_normalization": normalization,
                "adapted_filter_support": support,
                "candidate_local_filter_amplitude": float(row["central_filter_amplitude"]),
            }
        )
    return results


def product_header(
    base: fits.Header,
    role: str,
    radii: np.ndarray,
    sample_counts: list[int],
    amplitude_fraction: float,
    epsilon: float,
    minimum_support_fraction: float,
) -> fits.Header:
    """Create adapted-response FITS provenance from the reference final header."""
    header = base.copy()
    header["HIERARCH KLIP PSF PRODUCT SCHEMA"] = (1, "sparse radial response-product schema")
    header["HIERARCH KLIP PSF PRODUCT"] = (role, "product role")
    header["HIERARCH KLIP PSF COMPLETE"] = (0, "complete product set available")
    header["HIERARCH KLIP PSF FILTER"] = (1, "normalized filtering enabled")
    header["HIERARCH KLIP PSF FILTER MIN GOOD FRACTION"] = minimum_support_fraction
    header["HIERARCH KLIP PSF RESPONSE METHOD"] = "PAIRED_CENTRAL_REFIT"
    header["HIERARCH KLIP PSF SAMPLE RADII"] = ",".join(format(value, ".12g") for value in radii)
    header["HIERARCH KLIP PSF ACTUAL SAMPLES PER RADIUS"] = ",".join(
        str(value) for value in sample_counts
    )
    header["HIERARCH KLIP PSF CENTRAL AMPLITUDE FRACTION"] = amplitude_fraction
    header["HIERARCH KLIP PSF CENTRAL EPSILON"] = epsilon
    return header


def main() -> int:
    """Build the adapted radial model, filter the science cube, and write products."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_case", type=Path)
    parser.add_argument("measurement_experiment", type=Path)
    parser.add_argument("output_case", type=Path)
    parser.add_argument("--candidate-experiment", type=Path, required=True)
    parser.add_argument("--amplitude-fraction", type=float, default=1.0)
    parser.add_argument("--candidate-amplitude-fraction", type=float, default=1.0)
    parser.add_argument("--minimum-radius", type=float, default=6.0)
    parser.add_argument("--maximum-radius", type=float, default=60.0)
    parser.add_argument("--minimum-support-fraction", type=float, default=1.0)
    parser.add_argument("--filter-chunk-size", type=int, default=256)
    arguments = parser.parse_args()
    if (
        not math.isfinite(arguments.amplitude_fraction)
        or arguments.amplitude_fraction <= 0
        or not math.isfinite(arguments.candidate_amplitude_fraction)
        or arguments.candidate_amplitude_fraction <= 0
        or not math.isfinite(arguments.minimum_radius)
        or not math.isfinite(arguments.maximum_radius)
        or arguments.minimum_radius < 0
        or arguments.minimum_radius >= arguments.maximum_radius
        or not math.isfinite(arguments.minimum_support_fraction)
        or arguments.minimum_support_fraction < 0
        or arguments.minimum_support_fraction > 1
        or arguments.filter_chunk_size <= 0
    ):
        raise RuntimeError("adapted-grid amplitude, radius, support, or chunk setting is invalid")

    reference_case = arguments.reference_case.resolve()
    measurement_experiment = arguments.measurement_experiment.resolve()
    candidate_experiment = arguments.candidate_experiment.resolve()
    output_case = arguments.output_case.resolve()
    original_path = reference_case / "finim.fits"
    science, science_header = read_cube(original_path)
    science_modes = mode_counts(science_header, original_path)
    summary_path = measurement_experiment / "klip_central_response.json"
    central_path = measurement_experiment / "central_response.fits"
    for path in (summary_path, central_path):
        if not path.is_file():
            raise RuntimeError(f"paired central-response product is missing: {path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    central_cube = np.asarray(fits.getdata(central_path), dtype=np.float64)
    responses, validities, radii, modes, sample_counts = build_radial_models(
        summary, central_cube, arguments.amplitude_fraction
    )
    missing_modes = [mode for mode in modes if mode not in science_modes]
    if missing_modes:
        raise RuntimeError(
            f"measurement modes {missing_modes} are absent from reference science modes {science_modes}"
        )
    science_indices = [science_modes.index(mode) for mode in modes]
    science = science[np.asarray(science_indices)]
    science_header["NMODES"] = ",".join(str(mode) for mode in modes)

    reference_responses = sorted(
        (reference_case / "finim_outputs").glob("klipPSF_mode*_radial_response.fits")
    )
    if len(reference_responses) != len(science_modes):
        raise RuntimeError("reference case does not contain one radial response per science mode")
    expected_radii = np.asarray(
        [
            float(token)
            for token in str(fits.getheader(reference_responses[0])["KLIP PSF SAMPLE RADII"]).split(",")
        ]
    )
    if not np.array_equal(radii, expected_radii):
        raise RuntimeError(
            f"adapted grid radii {radii.tolist()} do not match the reference grid {expected_radii.tolist()}"
        )

    filtered, normalization, support, filter_validity = filter_science_cube(
        science,
        responses,
        validities,
        radii,
        arguments.minimum_radius,
        arguments.maximum_radius,
        arguments.minimum_support_fraction,
        arguments.filter_chunk_size,
    )
    candidate_rows = candidate_comparison(
        candidate_experiment,
        arguments.candidate_amplitude_fraction,
        modes,
        science,
        responses,
        validities,
        radii,
    )

    if output_case.exists() and any(output_case.iterdir()):
        raise RuntimeError(f"adapted-grid output directory is not empty: {output_case}")
    output_directory = output_case / "finim_outputs"
    output_directory.mkdir(parents=True, exist_ok=True)
    epsilon = float(summary["sample_manifest"]["planet"]["contrast"]) * arguments.amplitude_fraction
    base_header = product_header(
        science_header,
        "SCIENCE",
        radii,
        sample_counts,
        arguments.amplitude_fraction,
        epsilon,
        arguments.minimum_support_fraction,
    )
    base_header["HIERARCH KLIP PSF GRID CANDIDATE AVOID RADIUS"] = float(
        summary["sample_manifest"]["candidate_avoid_radius"]
    )
    fits.writeto(output_case / "finim.fits", science.astype(np.float32), header=base_header, overwrite=False)

    products = (
        (output_case / "finim_filtered.fits", filtered, "FILTERED"),
        (output_directory / "finim_filter_normalization.fits", normalization, "FILTER_NORMALIZATION"),
        (output_directory / "finim_filter_support.fits", support, "FILTER_SUPPORT"),
        (output_directory / "finim_filter_validity.fits", filter_validity, "FILTER_VALIDITY"),
    )
    for path, data, role in products:
        header = base_header.copy()
        header["HIERARCH KLIP PSF PRODUCT"] = role
        fits.writeto(path, np.asarray(data, dtype=np.float32), header=header, overwrite=False)

    for mode_index, retained_modes in enumerate(modes):
        response_header = base_header.copy()
        response_header["HIERARCH KLIP PSF PRODUCT"] = "RADIAL_RESPONSE"
        response_header["HIERARCH KLIP PSF MODE COUNT"] = retained_modes
        validity_header = response_header.copy()
        validity_header["HIERARCH KLIP PSF PRODUCT"] = "RADIAL_VALIDITY"
        fits.writeto(
            output_directory / f"klipPSF_mode{mode_index:03d}_radial_response.fits",
            responses[mode_index].transpose(0, 2, 1).astype(np.float32),
            header=response_header,
            overwrite=False,
        )
        fits.writeto(
            output_directory / f"klipPSF_mode{mode_index:03d}_radial_validity.fits",
            validities[mode_index].transpose(0, 2, 1).astype(np.float32),
            header=validity_header,
            overwrite=False,
        )

    manifest_header = base_header.copy()
    manifest_header["HIERARCH KLIP PSF PRODUCT"] = "MANIFEST"
    manifest_header["HIERARCH KLIP PSF COMPLETE"] = 1
    fits.writeto(
        output_directory / "klipPSF_manifest.fits",
        np.ones((1, 1), dtype=np.float32),
        header=manifest_header,
        overwrite=False,
    )

    summary_result = {
        "schema": 1,
        "reference_case": str(reference_case),
        "measurement_experiment": str(measurement_experiment),
        "candidate_experiment": str(candidate_experiment),
        "output_case": str(output_case),
        "mode_counts": modes,
        "radii": radii.tolist(),
        "actual_samples_per_radius": sample_counts,
        "amplitude_fraction": arguments.amplitude_fraction,
        "epsilon": epsilon,
        "minimum_support_fraction": arguments.minimum_support_fraction,
        "valid_filter_pixels_per_mode": [
            int(np.count_nonzero(filter_validity[index] > 0.5)) for index in range(len(modes))
        ],
        "candidate_comparison": candidate_rows,
    }
    summary_json = output_case / "adapted_grid_summary.json"
    summary_csv = output_case / "adapted_grid_candidate_comparison.csv"
    summary_markdown = output_case / "adapted_grid_summary.md"
    summary_json.write_text(
        json.dumps(finite_json(summary_result), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    with summary_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(candidate_rows[0]))
        writer.writeheader()
        writer.writerows(candidate_rows)
    with summary_markdown.open("w", encoding="utf-8") as stream:
        stream.write("# KLIP candidate-avoiding adapted radial grid\n\n")
        stream.write(
            f"Built `{len(radii)}` radial nodes from `{sum(sample_counts)}` clear paired-refit samples at "
            f"epsilon/fiducial `{arguments.amplitude_fraction:.8g}`. Actual angular sample counts were "
            f"`{','.join(str(value) for value in sample_counts)}`.\n\n"
        )
        stream.write(
            "| KL modes | Candidate projection on clear grid | Candidate/clear-grid cosine | "
            "Best-scaled residual | Adapted filter amplitude | Candidate-local amplitude |\n"
        )
        stream.write("|---:|---:|---:|---:|---:|---:|\n")
        for row in candidate_rows:
            stream.write(
                f"| {row['mode_count']} | {format_number(row['candidate_projection_on_clear_grid'])} | "
                f"{format_number(row['candidate_clear_grid_cosine'])} | "
                f"{format_number(row['candidate_clear_grid_best_scaled_residual'])} | "
                f"{format_number(row['adapted_filter_amplitude'])} | "
                f"{format_number(row['candidate_local_filter_amplitude'])} |\n"
            )
    print(f"Wrote {summary_json}")
    print(f"Wrote {summary_markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
