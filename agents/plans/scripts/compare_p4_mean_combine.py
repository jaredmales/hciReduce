#!/usr/bin/env python3
"""Compare sigma-clipped and mean-combined finite P4 responses."""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


def load_response_fit_module():
    """Load the maintained matched-response implementation without writing bytecode."""
    sys.dont_write_bytecode = True
    script = Path(__file__).resolve().with_name("fit_p4_matched_response.py")
    specification = importlib.util.spec_from_file_location("p4_matched_response_fit", script)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load matched-response implementation: {script}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def optimizer_summary_path(reference: Path) -> Path:
    """Return the unique P4 optimizer summary in a matched-response experiment."""
    candidates = sorted((reference / "exact_optimizer" / "finim_outputs").glob("p4Negative*_summary.yaml"))
    if len(candidates) != 1:
        raise RuntimeError(f"expected one exact-optimizer summary under {reference}")
    return candidates[0]


def read_exact_planet(summary: Path) -> dict[str, float]:
    """Read and validate the converged negative-planet point estimate."""
    text = summary.read_text(encoding="utf-8")
    status = re.search(r'^  status: "([^"]+)"$', text, flags=re.MULTILINE)
    converged = re.search(r"^  converged: (true|false)$", text, flags=re.MULTILINE)
    dense = re.search(r"^  denseAgreement: (true|false)$", text, flags=re.MULTILINE)
    fitted = re.search(r"^  fitted:\n(?P<body>(?:    [^\n]*\n)+)", text, flags=re.MULTILINE)
    if (
        status is None
        or status.group(1) != "converged"
        or converged is None
        or converged.group(1) != "true"
        or dense is None
        or dense.group(1) != "true"
        or fitted is None
    ):
        raise RuntimeError(f"exact optimizer point fit did not converge: {summary}")

    values: dict[str, float] = {}
    for key in ("separation", "positionAngle", "contrast"):
        match = re.search(rf"^    {key}: ([+\-0-9.eE]+)$", fitted.group("body"), flags=re.MULTILINE)
        if match is None:
            raise RuntimeError(f"exact optimizer summary lacks fitted {key}: {summary}")
        values[key] = float(match.group(1))
    if (
        not all(math.isfinite(value) for value in values.values())
        or values["separation"] < 0
        or values["contrast"] >= 0
    ):
        raise RuntimeError(f"exact optimizer point estimate is invalid: {summary}")
    return {
        "separation": values["separation"],
        "position_angle": values["positionAngle"] % 360,
        "contrast": -values["contrast"],
    }


def validate_combination(path: Path, expected: str) -> None:
    """Require one final image to advertise the expected combination method."""
    actual = str(fits.getheader(path).get("COMBINATION METHOD", "")).strip()
    if actual != expected:
        raise RuntimeError(f"expected {expected} combination in {path}, found {actual or 'unset'}")


def quadratic_at(response_fit, values: np.ndarray, row: float, column: float) -> float:
    """Interpolate a map at one subpixel coordinate with the maintained quadratic convention."""
    center_row = int(round(row))
    center_column = int(round(column))
    coefficients = response_fit.quadratic_coefficients(values, center_column, center_row)
    return response_fit.quadratic_value(
        coefficients,
        row - center_row,
        column - center_column,
    )


def pair_metrics(
    response_fit,
    label: str,
    combination: str,
    original_path: Path,
    signal_free_path: Path,
    detector_rows: int,
    detector_columns: int,
    mode_fraction: float,
    model: np.ndarray,
    source_validity: np.ndarray,
    coordinates: np.ndarray,
    minimum_support: float,
    source_row: float,
    source_column: float,
    exact_contrast: float,
) -> dict[str, object]:
    """Measure the finite original-minus-signal-free response for one final-image combiner."""
    validate_combination(original_path, combination)
    validate_combination(signal_free_path, combination)
    original, original_header = response_fit.read_science(original_path, detector_rows, detector_columns)
    signal_free, signal_free_header = response_fit.read_science(
        signal_free_path, detector_rows, detector_columns
    )
    original_mode = response_fit.selected_mode_index(original_header, mode_fraction, original_path)
    signal_free_mode = response_fit.selected_mode_index(signal_free_header, mode_fraction, signal_free_path)
    if original.shape[1:] != signal_free.shape[1:]:
        raise RuntimeError(f"original and signal-free detector dimensions differ for {label}")

    original_plane = original[original_mode]
    signal_free_plane = signal_free[signal_free_mode]
    difference = original_plane - signal_free_plane
    amplitude_maps: dict[str, np.ndarray] = {}
    for name, plane in (
        ("original", original_plane),
        ("signal_free", signal_free_plane),
        ("difference", difference),
    ):
        amplitude, _, _ = response_fit.apply_response(
            plane[np.newaxis, :, :],
            0,
            model,
            source_validity,
            coordinates,
            minimum_support,
        )
        amplitude_maps[name] = amplitude

    nearest_row = int(round(source_row))
    nearest_column = int(round(source_column))
    source_matches = np.flatnonzero(
        (coordinates[:, 0] == nearest_row) & (coordinates[:, 1] == nearest_column)
    )
    if source_matches.size != 1:
        raise RuntimeError(f"response field does not uniquely contain ({nearest_row}, {nearest_column})")
    source = int(source_matches[0])
    if not math.isfinite(float(source_validity[source])) or source_validity[source] <= 0:
        raise RuntimeError(f"response field is invalid at ({nearest_row}, {nearest_column})")

    response = model[source]
    half_width = response.shape[0] // 2
    difference_patch = difference[
        nearest_column - half_width : nearest_column + half_width + 1,
        nearest_row - half_width : nearest_row + half_width + 1,
    ]
    retained = np.isfinite(response) & np.isfinite(difference_patch)
    response_values = response[retained]
    difference_values = difference_patch[retained]
    response_energy = float(np.dot(response_values, response_values))
    difference_energy = float(np.dot(difference_values, difference_values))
    correlation = float(np.dot(response_values, difference_values))
    if response_energy <= 0 or difference_energy <= 0:
        raise RuntimeError(f"finite response comparison is degenerate for {label}")
    integer_amplitude = correlation / response_energy
    cosine = correlation / math.sqrt(response_energy * difference_energy)
    best_residual = difference_values - integer_amplitude * response_values
    expected_residual = difference_values - exact_contrast * response_values

    projected = {
        name: quadratic_at(response_fit, amplitude, source_row, source_column)
        for name, amplitude in amplitude_maps.items()
    }
    return {
        "label": label,
        "combination": combination,
        "original": str(original_path),
        "signal_free": str(signal_free_path),
        "source_row": source_row,
        "source_column": source_column,
        "nearest_source_row": nearest_row,
        "nearest_source_column": nearest_column,
        "original_projected_amplitude": projected["original"],
        "signal_free_projected_amplitude": projected["signal_free"],
        "difference_projected_amplitude": projected["difference"],
        "difference_ratio_to_exact": projected["difference"] / exact_contrast,
        "integer_difference_amplitude": integer_amplitude,
        "response_difference_cosine": cosine,
        "best_scale_residual_fraction": float(np.linalg.norm(best_residual)) / math.sqrt(difference_energy),
        "exact_scale_residual_fraction": float(np.linalg.norm(expected_residual)) / math.sqrt(difference_energy),
    }


def format_number(value: object) -> str:
    """Format a finite diagnostic value compactly."""
    number = float(value)
    if not math.isfinite(number):
        raise RuntimeError("mean-combination diagnostic produced a non-finite value")
    return f"{number:.8g}"


def main() -> int:
    """Write the sigmaMean-versus-mean finite-response comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_experiment", type=Path)
    parser.add_argument("mean_original", type=Path)
    parser.add_argument("mean_signal_free", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--mode-fraction", type=float, default=0.15)
    parser.add_argument("--candidate-combination", default="mean", choices=("mean", "sigmaMean"))
    arguments = parser.parse_args()

    reference = arguments.reference_experiment.resolve()
    mean_original = arguments.mean_original.resolve()
    mean_signal_free = arguments.mean_signal_free.resolve()
    output_directory = arguments.output_directory.resolve()
    if not all(path.is_file() for path in (mean_original, mean_signal_free)):
        raise RuntimeError("candidate original and signal-free final images must be readable files")

    response_fit = load_response_fit_module()
    optimizer_summary = optimizer_summary_path(reference)
    exact = read_exact_planet(optimizer_summary)
    manifest = reference / "signal_free_oracle" / "finim_outputs" / "p4PSF_manifest.fits"
    if not manifest.is_file():
        raise RuntimeError(f"dense signal-free response manifest is missing: {manifest}")
    manifest_header = fits.getheader(manifest)
    if int(manifest_header.get("P4 PSF COMPLETE", 0)) != 1:
        raise RuntimeError(f"dense signal-free response manifest is incomplete: {manifest}")
    if str(manifest_header.get("P4 PSF SPATIAL MODEL", "")).strip() != "PER_PIXEL":
        raise RuntimeError(f"mean-combination diagnostic requires a dense response: {manifest}")
    if str(manifest_header.get("P4 PSF COMBINATION", "")).strip() != "mean":
        raise RuntimeError(f"reference response was not mean-combined: {manifest}")
    detector_rows = int(manifest_header.get("P4 PSF TEMPLATE ROWS", 0))
    detector_columns = int(manifest_header.get("P4 PSF TEMPLATE COLUMNS", 0))
    if detector_rows <= 0 or detector_columns <= 0:
        raise RuntimeError(f"response manifest lacks detector dimensions: {manifest}")
    manifest_mode = response_fit.selected_mode_index(manifest_header, arguments.mode_fraction, manifest)
    directory, prefix = response_fit.product_prefix(manifest)
    coordinates = response_fit.read_coordinates(directory / f"{prefix}coordinates.fits")
    model = np.asarray(
        fits.getdata(directory / f"{prefix}model_{manifest_mode:04d}.fits"), dtype=np.float64
    )
    source_validity = np.asarray(
        fits.getdata(directory / f"{prefix}validity_{manifest_mode:04d}.fits"), dtype=np.float64
    ).reshape(-1)
    minimum_support = float(manifest_header.get("P4 PSF FILTER MIN GOOD FRACTION", 1.0))

    center_row = 0.5 * (detector_rows - 1)
    center_column = 0.5 * (detector_columns - 1)
    radians = math.radians(exact["position_angle"])
    source_row = center_row - exact["separation"] * math.sin(radians)
    source_column = center_column + exact["separation"] * math.cos(radians)
    common = (
        response_fit,
        detector_rows,
        detector_columns,
        arguments.mode_fraction,
        model,
        source_validity,
        coordinates,
        minimum_support,
        source_row,
        source_column,
        exact["contrast"],
    )
    baseline = pair_metrics(
        common[0],
        "sigma_clipped_baseline",
        "sigmaMean",
        reference / "sparse_response" / "finim.fits",
        reference / "signal_free_oracle" / "finim.fits",
        *common[1:],
    )
    candidate = pair_metrics(
        common[0],
        "mean_candidate",
        arguments.candidate_combination,
        mean_original,
        mean_signal_free,
        *common[1:],
    )

    result = {
        "schema": 1,
        "mode_fraction": arguments.mode_fraction,
        "response_manifest": str(manifest),
        "response_combination": "mean",
        "optimizer_summary": str(optimizer_summary),
        "exact_planet": exact,
        "pairs": [baseline, candidate],
    }
    output_directory.mkdir(parents=True, exist_ok=True)
    json_path = output_directory / "mean_combine_comparison.json"
    with json_path.open("w", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")

    markdown_path = output_directory / "mean_combine_comparison.md"
    with markdown_path.open("w", encoding="utf-8") as stream:
        stream.write("# P4 mean-combination finite-response diagnostic\n\n")
        stream.write(
            f"Exact negative-fit contrast: `{exact['contrast']:.12g}` at separation "
            f"`{exact['separation']:.12g}` pixels and PA `{exact['position_angle']:.12g}` degrees.\n\n"
        )
        stream.write(
            "Both rows project the original-minus-signal-free final-image difference onto the same dense, "
            "signal-free, mean-combined analytic response. The projection is evaluated at the exact fitted "
            "subpixel position.\n\n"
        )
        stream.write(
            "| Final-image combination | Original projection | Signal-free projection | Difference projection | "
            "Difference / exact | Shape cosine | Best-scale residual |\n"
        )
        stream.write("|---|---:|---:|---:|---:|---:|---:|\n")
        for pair in (baseline, candidate):
            stream.write(
                f"| {pair['combination']} | {format_number(pair['original_projected_amplitude'])} | "
                f"{format_number(pair['signal_free_projected_amplitude'])} | "
                f"{format_number(pair['difference_projected_amplitude'])} | "
                f"{format_number(pair['difference_ratio_to_exact'])} | "
                f"{format_number(pair['response_difference_cosine'])} | "
                f"{format_number(pair['best_scale_residual_fraction'])} |\n"
            )
        stream.write(
            "\nIf the mean-combined difference approaches the exact contrast, science sigma-clipping caused the "
            "calibration loss. If it remains near the sigmaMean ratio, refitting the P4 regression is the dominant "
            "missing response term. An intermediate result indicates that both effects matter.\n"
        )

    print(f"Wrote {json_path}")
    print(f"Wrote {markdown_path}")
    print(
        "Mean-combination finite response: "
        f"contrast={candidate['difference_projected_amplitude']:.12g}, "
        f"ratio_to_exact={candidate['difference_ratio_to_exact']:.8g}, "
        f"sigmaMean_ratio={baseline['difference_ratio_to_exact']:.8g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
