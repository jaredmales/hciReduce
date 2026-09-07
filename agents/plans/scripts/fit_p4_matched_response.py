#!/usr/bin/env python3
"""Fit one source with a persisted P4 response field and matched filter."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits


MODE_KEYWORDS = ("P4 MODE FRACTIONS", "P4MODFR", "FRACT NMODES", "NMODES")


def header_vector(header: fits.Header, keywords: tuple[str, ...]) -> np.ndarray:
    """Read the first available comma-separated numeric FITS vector."""
    for keyword in keywords:
        if keyword not in header:
            continue
        text = str(header[keyword]).strip()
        if not text:
            return np.empty(0, dtype=np.float64)
        return np.asarray([float(value) for value in text.split(",")], dtype=np.float64)
    raise RuntimeError(f"none of the required FITS keywords are present: {', '.join(keywords)}")


def selected_mode_index(header: fits.Header, requested_mode: float, source: Path) -> int:
    """Resolve one requested mode against a product header without nearest-mode substitution."""
    modes = header_vector(header, MODE_KEYWORDS)
    matches = np.flatnonzero(
        np.abs(modes - requested_mode)
        <= 8
        * np.finfo(np.float32).eps
        * np.maximum.reduce((np.ones_like(modes), np.abs(modes), np.full_like(modes, abs(requested_mode))))
    )
    if matches.size != 1:
        raise RuntimeError(f"mode {requested_mode:.17g} is not an exact unique mode in {source}")
    return int(matches[0])


def product_prefix(manifest: Path) -> tuple[Path, str]:
    """Return the product directory and filename prefix represented by a manifest path."""
    suffix = "manifest.fits"
    if not manifest.name.endswith(suffix):
        raise RuntimeError(f"P4 response manifest name must end with {suffix}: {manifest}")
    return manifest.parent, manifest.name[: -len(suffix)]


def read_science(science_path: Path, detector_rows: int, detector_columns: int) -> tuple[np.ndarray, fits.Header]:
    """Read a science cube and restore an automatic centered crop when necessary."""
    science_header = fits.getheader(science_path)
    science = np.asarray(fits.getdata(science_path), dtype=np.float64)
    if science.ndim == 2:
        science = science[np.newaxis, :, :]
    if science.ndim != 3:
        raise RuntimeError(f"science input must be a two- or three-dimensional FITS image: {science_path}")

    # FITS/NumPy axes are (plane, detector column, detector row) for the Eigen images written here.
    target_shape = (science.shape[0], detector_columns, detector_rows)
    if science.shape == target_shape:
        return science, science_header
    if science.shape[1] > detector_columns or science.shape[2] > detector_rows:
        raise RuntimeError(
            f"science dimensions {science.shape[1:]} exceed response detector dimensions "
            f"{(detector_columns, detector_rows)}"
        )
    column_difference = detector_columns - science.shape[1]
    row_difference = detector_rows - science.shape[2]
    if column_difference % 2 != 0 or row_difference % 2 != 0:
        raise RuntimeError("a cropped science product cannot be centered on the response detector dimensions")
    restored = np.full(target_shape, np.nan, dtype=np.float64)
    first_column = column_difference // 2
    first_row = row_difference // 2
    restored[
        :,
        first_column : first_column + science.shape[1],
        first_row : first_row + science.shape[2],
    ] = science
    return restored, science_header


def read_coordinates(path: Path) -> np.ndarray:
    """Read the coordinate product as rows of row, column, region, and regional index."""
    coordinates = np.asarray(fits.getdata(path), dtype=np.float64)
    if coordinates.ndim != 2:
        raise RuntimeError(f"P4 coordinate product must be two-dimensional: {path}")
    if coordinates.shape[0] == 4:
        coordinates = coordinates.T
    if coordinates.shape[1] != 4:
        raise RuntimeError(f"P4 coordinate product must contain four columns: {path}")
    if not np.all(np.isfinite(coordinates[:, :2])) or not np.all(coordinates[:, :2] == np.trunc(coordinates[:, :2])):
        raise RuntimeError(f"P4 coordinate product contains non-integer detector coordinates: {path}")
    return coordinates


def apply_response(
    science: np.ndarray,
    science_mode: int,
    model: np.ndarray,
    source_validity: np.ndarray,
    coordinates: np.ndarray,
    minimum_support: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply the normalized P4 response at every persisted source coordinate."""
    if model.ndim != 3 or model.shape[0] != coordinates.shape[0]:
        raise RuntimeError("P4 model and coordinate source counts differ")
    if model.shape[1] != model.shape[2] or model.shape[1] % 2 == 0:
        raise RuntimeError("P4 response stamps must be odd and square")
    if source_validity.size != coordinates.shape[0]:
        raise RuntimeError("P4 source validity and coordinate counts differ")
    if not math.isfinite(minimum_support) or minimum_support < 0 or minimum_support > 1:
        raise RuntimeError("P4 response minimum support must be finite and between zero and one")

    detector_columns = science.shape[1]
    detector_rows = science.shape[2]
    amplitude = np.full((detector_columns, detector_rows), np.nan, dtype=np.float64)
    normalization = np.full_like(amplitude, np.nan)
    support = np.full_like(amplitude, np.nan)
    half_width = model.shape[1] // 2
    stamp_pixels = model.shape[1] * model.shape[2]

    for source, coordinate in enumerate(coordinates):
        if not math.isfinite(float(source_validity[source])) or source_validity[source] <= 0:
            continue
        row = int(coordinate[0])
        column = int(coordinate[1])
        if (
            row - half_width < 0
            or row + half_width >= detector_rows
            or column - half_width < 0
            or column + half_width >= detector_columns
        ):
            continue
        response = model[source]
        patch = science[
            science_mode,
            column - half_width : column + half_width + 1,
            row - half_width : row + half_width + 1,
        ]
        retained = np.isfinite(response) & np.isfinite(patch)
        support_fraction = float(np.count_nonzero(retained)) / stamp_pixels
        support[column, row] = support_fraction
        if support_fraction < minimum_support:
            continue
        response_values = response[retained]
        science_values = patch[retained]
        denominator = float(np.dot(response_values, response_values))
        if not math.isfinite(denominator) or denominator <= 0:
            continue
        numerator = float(np.dot(response_values, science_values))
        if not math.isfinite(numerator):
            continue
        normalization[column, row] = denominator
        amplitude[column, row] = numerator / denominator

    return amplitude, normalization, support


def radial_noise_maps(
    amplitude: np.ndarray,
    source_row: float,
    source_column: float,
    exclusion_radius: float,
    minimum_radius: float,
    maximum_radius: float,
    lambda_d: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reproduce hciAnalyze's one-pixel radial mean, deviation, and small-sample SNR profiles."""
    detector_columns, detector_rows = amplitude.shape
    column_indices, row_indices = np.indices(amplitude.shape, dtype=np.float64)
    center_row = 0.5 * (detector_rows - 1)
    center_column = 0.5 * (detector_columns - 1)
    radii = np.hypot(row_indices - center_row, column_indices - center_column)
    signal_distance = np.hypot(row_indices - source_row, column_indices - source_column)
    noise_mask = np.isfinite(amplitude) & (signal_distance > exclusion_radius + 0.5)

    profile_radii: list[float] = []
    profile_means: list[float] = []
    profile_deviations: list[float] = []
    upper_radius = int(math.ceil(float(np.max(radii))))
    for lower in range(upper_radius):
        selected = noise_mask & (radii >= lower) & (radii <= lower + 1)
        values = amplitude[selected]
        if values.size < 2:
            continue
        profile_radii.append(lower + 0.5)
        profile_means.append(float(np.mean(values)))
        profile_deviations.append(float(np.std(values, ddof=1)))
    if len(profile_radii) < 2:
        raise RuntimeError("matched-filter amplitude map has too few radial noise samples")

    mean_map = np.interp(radii, profile_radii, profile_means)
    deviation_map = np.interp(radii, profile_radii, profile_deviations)
    comparison_samples = 2 * np.pi * radii / lambda_d - 1
    correction = np.zeros_like(radii)
    positive = comparison_samples > 0
    correction[positive] = 1 / np.sqrt(1 + 1 / comparison_samples[positive])
    uncertainty = np.full_like(amplitude, np.nan)
    usable_deviation = np.isfinite(deviation_map) & (deviation_map > 0) & (correction > 0)
    uncertainty[usable_deviation] = deviation_map[usable_deviation] / correction[usable_deviation]
    signal = amplitude - mean_map
    snr = np.full_like(amplitude, np.nan)
    usable = (
        np.isfinite(amplitude)
        & usable_deviation
        & (radii >= minimum_radius)
        & (radii <= maximum_radius)
    )
    snr[usable] = signal[usable] / uncertainty[usable]
    return signal, uncertainty, snr, radii


def quadratic_coefficients(values: np.ndarray, center_column: int, center_row: int) -> np.ndarray:
    """Fit a complete quadratic to the three-by-three neighborhood of one detector pixel."""
    design: list[list[float]] = []
    samples: list[float] = []
    for column in range(center_column - 1, center_column + 2):
        for row in range(center_row - 1, center_row + 2):
            value = float(values[column, row])
            if not math.isfinite(value):
                raise RuntimeError("matched-filter peak does not have a complete finite three-by-three neighborhood")
            delta_row = float(row - center_row)
            delta_column = float(column - center_column)
            design.append(
                [
                    1,
                    delta_row,
                    delta_column,
                    delta_row * delta_row,
                    delta_row * delta_column,
                    delta_column * delta_column,
                ]
            )
            samples.append(value)
    coefficients, _, rank, _ = np.linalg.lstsq(np.asarray(design), np.asarray(samples), rcond=None)
    if rank != 6:
        raise RuntimeError("matched-filter quadratic surface is rank deficient")
    return coefficients


def quadratic_value(coefficients: np.ndarray, delta_row: float, delta_column: float) -> float:
    """Evaluate a fitted two-dimensional quadratic surface."""
    return float(
        coefficients
        @ np.asarray(
            [
                1,
                delta_row,
                delta_column,
                delta_row * delta_row,
                delta_row * delta_column,
                delta_column * delta_column,
            ]
        )
    )


def bounded_quadratic_maximum(
    coefficients: np.ndarray,
    row_bounds: tuple[float, float],
    column_bounds: tuple[float, float],
) -> tuple[float, float, np.ndarray, bool]:
    """Maximize a concave quadratic over a rectangular subpixel domain."""
    gradient = coefficients[1:3]
    hessian = np.asarray(
        [[2 * coefficients[3], coefficients[4]], [coefficients[4], 2 * coefficients[5]]], dtype=np.float64
    )
    eigenvalues = np.linalg.eigvalsh(hessian)
    if not np.all(eigenvalues < 0):
        raise RuntimeError("matched-filter likelihood peak is not locally concave")

    candidates: list[tuple[float, float]] = []
    stationary = -np.linalg.solve(hessian, gradient)
    if row_bounds[0] <= stationary[0] <= row_bounds[1] and column_bounds[0] <= stationary[1] <= column_bounds[1]:
        candidates.append((float(stationary[0]), float(stationary[1])))
    for row in row_bounds:
        column = -(gradient[1] + hessian[1, 0] * row) / hessian[1, 1]
        candidates.append((row, float(np.clip(column, *column_bounds))))
    for column in column_bounds:
        row = -(gradient[0] + hessian[0, 1] * column) / hessian[0, 0]
        candidates.append((float(np.clip(row, *row_bounds)), column))
    for row in row_bounds:
        for column in column_bounds:
            candidates.append((row, column))
    best = max(candidates, key=lambda point: quadratic_value(coefficients, *point))
    tolerance = 32 * np.finfo(np.float64).eps
    at_boundary = (
        abs(best[0] - row_bounds[0]) <= tolerance
        or abs(best[0] - row_bounds[1]) <= tolerance
        or abs(best[1] - column_bounds[0]) <= tolerance
        or abs(best[1] - column_bounds[1]) <= tolerance
    )
    return best[0], best[1], hessian, at_boundary


def polar_coordinate(row: float, column: float, center_row: float, center_column: float) -> tuple[float, float]:
    """Convert P4 detector row/column offsets to separation and east-of-north PA."""
    row_delta = row - center_row
    column_delta = column - center_column
    separation = math.hypot(row_delta, column_delta)
    position_angle = (-math.degrees(math.atan2(row_delta, column_delta))) % 360
    return separation, position_angle


def covariance_summary(
    covariance: np.ndarray,
    row: float,
    column: float,
    center_row: float,
    center_column: float,
) -> dict[str, float | list[list[float]]]:
    """Propagate local likelihood curvature into Cartesian, separation, and PA diagnostics."""
    row_delta = row - center_row
    column_delta = column - center_column
    radius_squared = row_delta * row_delta + column_delta * column_delta
    if radius_squared <= 0:
        raise RuntimeError("cannot propagate position covariance at zero separation")
    separation = math.sqrt(radius_squared)
    separation_jacobian = np.asarray([row_delta / separation, column_delta / separation])
    pa_jacobian = (180 / np.pi) * np.asarray([-column_delta / radius_squared, row_delta / radius_squared])
    separation_variance = float(separation_jacobian @ covariance @ separation_jacobian)
    pa_variance = float(pa_jacobian @ covariance @ pa_jacobian)
    return {
        "row_standard_error": math.sqrt(max(0.0, float(covariance[0, 0]))),
        "column_standard_error": math.sqrt(max(0.0, float(covariance[1, 1]))),
        "separation_standard_error": math.sqrt(max(0.0, separation_variance)),
        "position_angle_standard_error": math.sqrt(max(0.0, pa_variance)),
        "row_column_covariance": covariance.tolist(),
    }


def fit_response(
    amplitude: np.ndarray,
    normalization: np.ndarray,
    support: np.ndarray,
    initial_separation: float,
    initial_pa: float,
    position_bound: float,
    exclusion_radius: float,
    minimum_radius: float,
    maximum_radius: float,
    lambda_d: float,
) -> tuple[dict[str, object], list[dict[str, float | int]]]:
    """Fit a bounded subpixel matched-filter peak and return its diagnostic surface."""
    detector_columns, detector_rows = amplitude.shape
    center_row = 0.5 * (detector_rows - 1)
    center_column = 0.5 * (detector_columns - 1)
    initial_radians = math.radians(initial_pa)
    initial_row = center_row - initial_separation * math.sin(initial_radians)
    initial_column = center_column + initial_separation * math.cos(initial_radians)
    signal, uncertainty, snr, radii = radial_noise_maps(
        amplitude,
        initial_row,
        initial_column,
        exclusion_radius,
        minimum_radius,
        maximum_radius,
        lambda_d,
    )
    likelihood = 0.5 * snr * snr

    column_indices, row_indices = np.indices(amplitude.shape)
    bounded = (
        np.isfinite(likelihood)
        & np.isfinite(signal)
        & (signal > 0)
        & (row_indices >= initial_row - position_bound)
        & (row_indices <= initial_row + position_bound)
        & (column_indices >= initial_column - position_bound)
        & (column_indices <= initial_column + position_bound)
    )
    if not np.any(bounded):
        raise RuntimeError("no positive valid matched-filter sample lies inside the requested position bounds")
    bounded_likelihood = np.where(bounded, likelihood, -np.inf)
    peak_column, peak_row = np.unravel_index(int(np.argmax(bounded_likelihood)), bounded_likelihood.shape)
    likelihood_coefficients = quadratic_coefficients(likelihood, peak_column, peak_row)

    row_bounds = (
        max(-1.0, initial_row - position_bound - peak_row),
        min(1.0, initial_row + position_bound - peak_row),
    )
    column_bounds = (
        max(-1.0, initial_column - position_bound - peak_column),
        min(1.0, initial_column + position_bound - peak_column),
    )
    if row_bounds[0] > row_bounds[1] or column_bounds[0] > column_bounds[1]:
        raise RuntimeError("the matched-filter interpolation neighborhood does not intersect the position bounds")
    delta_row, delta_column, hessian, at_boundary = bounded_quadratic_maximum(
        likelihood_coefficients, row_bounds, column_bounds
    )
    fitted_row = peak_row + delta_row
    fitted_column = peak_column + delta_column
    fitted_separation, fitted_pa = polar_coordinate(fitted_row, fitted_column, center_row, center_column)

    signal_coefficients = quadratic_coefficients(signal, peak_column, peak_row)
    amplitude_coefficients = quadratic_coefficients(amplitude, peak_column, peak_row)
    uncertainty_coefficients = quadratic_coefficients(uncertainty, peak_column, peak_row)
    snr_coefficients = quadratic_coefficients(snr, peak_column, peak_row)
    normalization_coefficients = quadratic_coefficients(normalization, peak_column, peak_row)
    support_coefficients = quadratic_coefficients(support, peak_column, peak_row)
    fitted_signal = quadratic_value(signal_coefficients, delta_row, delta_column)
    fitted_uncertainty = quadratic_value(uncertainty_coefficients, delta_row, delta_column)
    fitted_likelihood = quadratic_value(likelihood_coefficients, delta_row, delta_column)
    covariance = np.linalg.inv(-hessian)

    fit: dict[str, object] = {
        "status": "bounded_peak" if at_boundary else "converged",
        "at_bound": bool(at_boundary),
        "detector_center_row": center_row,
        "detector_center_column": center_column,
        "initial_row": initial_row,
        "initial_column": initial_column,
        "row": fitted_row,
        "column": fitted_column,
        "row_delta": fitted_row - initial_row,
        "column_delta": fitted_column - initial_column,
        "separation": fitted_separation,
        "position_angle": fitted_pa,
        "contrast": fitted_signal,
        "raw_amplitude": quadratic_value(amplitude_coefficients, delta_row, delta_column),
        "contrast_standard_error": fitted_uncertainty,
        "snr": quadratic_value(snr_coefficients, delta_row, delta_column),
        "likelihood_snr": math.sqrt(max(0.0, 2 * fitted_likelihood)),
        "response_normalization": quadratic_value(normalization_coefficients, delta_row, delta_column),
        "support_fraction": quadratic_value(support_coefficients, delta_row, delta_column),
        "integer_peak_row": int(peak_row),
        "integer_peak_column": int(peak_column),
        "curvature_diagnostic": covariance_summary(
            covariance, fitted_row, fitted_column, center_row, center_column
        ),
    }

    surface: list[dict[str, float | int]] = []
    margin = int(math.ceil(position_bound)) + 2
    for column in range(max(0, math.floor(initial_column) - margin), min(detector_columns, math.ceil(initial_column) + margin + 1)):
        for row in range(max(0, math.floor(initial_row) - margin), min(detector_rows, math.ceil(initial_row) + margin + 1)):
            separation, position_angle = polar_coordinate(row, column, center_row, center_column)
            surface.append(
                {
                    "row": row,
                    "column": column,
                    "separation": separation,
                    "position_angle": position_angle,
                    "amplitude": float(amplitude[column, row]),
                    "signal": float(signal[column, row]),
                    "uncertainty": float(uncertainty[column, row]),
                    "snr": float(snr[column, row]),
                    "likelihood": float(likelihood[column, row]),
                    "normalization": float(normalization[column, row]),
                    "support": float(support[column, row]),
                    "inside_position_bound": int(
                        abs(row - initial_row) <= position_bound and abs(column - initial_column) <= position_bound
                    ),
                }
            )
    fit["finite_amplitude_count"] = int(np.count_nonzero(np.isfinite(amplitude)))
    fit["finite_snr_count"] = int(np.count_nonzero(np.isfinite(snr)))
    fit["minimum_noise_radius"] = minimum_radius
    fit["maximum_noise_radius"] = maximum_radius
    fit["noise_exclusion_radius"] = exclusion_radius
    fit["lambda_d"] = lambda_d
    fit["source_radius_at_fit"] = float(
        np.hypot(fitted_row - center_row, fitted_column - center_column)
    )
    return fit, surface


def write_image(path: Path, data: np.ndarray, role: str, mode_fraction: float) -> None:
    """Write one fit diagnostic image with compact, unambiguous provenance."""
    header = fits.Header()
    header["HIERARCH P4 RESPONSE FIT PRODUCT"] = role
    header["HIERARCH P4 RESPONSE FIT MODE FRACTION"] = mode_fraction
    fits.writeto(path, np.asarray(data, dtype=np.float32), header=header, overwrite=False)


def main() -> int:
    """Apply a P4 response field, fit its local likelihood peak, and persist diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("science", type=Path, help="ordinary P4 final science image")
    parser.add_argument("manifest", type=Path, help="P4 response manifest used as the matched-filter template")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--mode-fraction", type=float, default=0.15)
    parser.add_argument("--initial-separation", type=float, required=True)
    parser.add_argument("--initial-pa", type=float, required=True)
    parser.add_argument("--position-bound", type=float, default=1.0)
    parser.add_argument("--noise-exclusion-radius", type=float, default=5.0)
    parser.add_argument("--noise-min-radius", type=float, default=6.0)
    parser.add_argument("--noise-max-radius", type=float, default=60.0)
    parser.add_argument("--lambda-d", type=float, default=3.6)
    arguments = parser.parse_args()

    for name in (
        "mode_fraction",
        "initial_separation",
        "initial_pa",
        "position_bound",
        "noise_exclusion_radius",
        "noise_min_radius",
        "noise_max_radius",
        "lambda_d",
    ):
        value = float(getattr(arguments, name))
        if not math.isfinite(value):
            raise RuntimeError(f"--{name.replace('_', '-')} must be finite")
    if arguments.initial_separation < 0 or arguments.position_bound <= 0 or arguments.noise_exclusion_radius < 0:
        raise RuntimeError("separation/exclusion must be nonnegative and position bound must be positive")
    if arguments.noise_min_radius < 0 or arguments.noise_min_radius >= arguments.noise_max_radius:
        raise RuntimeError("noise radii must be nonnegative and strictly ordered")
    if arguments.lambda_d <= 0:
        raise RuntimeError("--lambda-d must be positive")

    science_path = arguments.science.resolve()
    manifest_path = arguments.manifest.resolve()
    output_directory = arguments.output_directory.resolve()
    if not science_path.is_file() or not manifest_path.is_file():
        raise RuntimeError("science and manifest inputs must be readable files")
    output_directory.mkdir(parents=True, exist_ok=True)
    managed_outputs = {
        "amplitude.fits",
        "signal.fits",
        "uncertainty.fits",
        "snr.fits",
        "normalization.fits",
        "support.fits",
        "summary.json",
        "surface.csv",
    }
    existing_outputs = sorted(name for name in managed_outputs if (output_directory / name).exists())
    if existing_outputs:
        raise RuntimeError(
            f"output directory already contains matched-response products: {', '.join(existing_outputs)}"
        )

    manifest_header = fits.getheader(manifest_path)
    if int(manifest_header.get("P4 PSF COMPLETE", 0)) != 1:
        raise RuntimeError(f"P4 response manifest is incomplete: {manifest_path}")
    detector_rows = int(manifest_header.get("P4 PSF TEMPLATE ROWS", 0))
    detector_columns = int(manifest_header.get("P4 PSF TEMPLATE COLUMNS", 0))
    mode_count = int(manifest_header.get("P4 PSF MODE COUNT", 0))
    if detector_rows <= 0 or detector_columns <= 0 or mode_count <= 0:
        raise RuntimeError(f"P4 response manifest lacks detector or mode dimensions: {manifest_path}")
    manifest_mode = selected_mode_index(manifest_header, arguments.mode_fraction, manifest_path)
    if manifest_mode >= mode_count:
        raise RuntimeError("selected manifest mode exceeds the declared mode count")

    science, science_header = read_science(science_path, detector_rows, detector_columns)
    science_mode = selected_mode_index(science_header, arguments.mode_fraction, science_path)
    if science_mode >= science.shape[0]:
        raise RuntimeError("selected science mode exceeds the science cube plane count")
    directory, prefix = product_prefix(manifest_path)
    coordinate_path = directory / f"{prefix}coordinates.fits"
    model_path = directory / f"{prefix}model_{manifest_mode:04d}.fits"
    validity_path = directory / f"{prefix}validity_{manifest_mode:04d}.fits"
    for path in (coordinate_path, model_path, validity_path):
        if not path.is_file():
            raise RuntimeError(f"missing P4 response product: {path}")

    coordinates = read_coordinates(coordinate_path)
    model = np.asarray(fits.getdata(model_path), dtype=np.float64)
    source_validity = np.asarray(fits.getdata(validity_path), dtype=np.float64).reshape(-1)
    minimum_support = float(manifest_header.get("P4 PSF FILTER MIN GOOD FRACTION", 1.0))
    amplitude, normalization, support = apply_response(
        science, science_mode, model, source_validity, coordinates, minimum_support
    )
    fit, surface = fit_response(
        amplitude,
        normalization,
        support,
        arguments.initial_separation,
        arguments.initial_pa,
        arguments.position_bound,
        arguments.noise_exclusion_radius,
        arguments.noise_min_radius,
        arguments.noise_max_radius,
        arguments.lambda_d,
    )

    signal, uncertainty, snr, _ = radial_noise_maps(
        amplitude,
        float(fit["initial_row"]),
        float(fit["initial_column"]),
        arguments.noise_exclusion_radius,
        arguments.noise_min_radius,
        arguments.noise_max_radius,
        arguments.lambda_d,
    )
    write_image(output_directory / "amplitude.fits", amplitude, "MATCHED_FILTER_AMPLITUDE", arguments.mode_fraction)
    write_image(output_directory / "signal.fits", signal, "RADIAL_MEAN_SUBTRACTED_AMPLITUDE", arguments.mode_fraction)
    write_image(output_directory / "uncertainty.fits", uncertainty, "RADIAL_NOISE_UNCERTAINTY", arguments.mode_fraction)
    write_image(output_directory / "snr.fits", snr, "SMALL_SAMPLE_CORRECTED_SNR", arguments.mode_fraction)
    write_image(
        output_directory / "normalization.fits",
        normalization,
        "MATCHED_FILTER_RESPONSE_ENERGY",
        arguments.mode_fraction,
    )
    write_image(output_directory / "support.fits", support, "MATCHED_FILTER_SUPPORT", arguments.mode_fraction)

    summary = {
        "schema": 1,
        "science": str(science_path),
        "manifest": str(manifest_path),
        "manifest_schema": int(manifest_header.get("P4 PSF PRODUCT SCHEMA", 0)),
        "spatial_model": str(manifest_header.get("P4 PSF SPATIAL MODEL", "PER_PIXEL")).strip(),
        "composition": str(manifest_header.get("P4 PSF COMPOSITION", "SOURCE_PIXEL")).strip(),
        "mode_fraction": arguments.mode_fraction,
        "science_mode_index": science_mode,
        "manifest_mode_index": manifest_mode,
        "position_bound": arguments.position_bound,
        "initial_separation": arguments.initial_separation,
        "initial_position_angle": arguments.initial_pa,
        "fit": fit,
    }
    summary_path = output_directory / "summary.json"
    with summary_path.open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    surface_path = output_directory / "surface.csv"
    with surface_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(surface[0]))
        writer.writeheader()
        writer.writerows(surface)

    print(f"Wrote {summary_path}")
    print(
        "Matched-response fit: "
        f"sep={fit['separation']:.12g}, PA={fit['position_angle']:.12g}, "
        f"contrast={fit['contrast']:.12g}, SNR={fit['snr']:.8g}, status={fit['status']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
