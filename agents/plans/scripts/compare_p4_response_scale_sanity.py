#!/usr/bin/env python3
"""Compare negative-injection and scaled-template P4 response sanity checks."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


def load_response_fit_module():
    """Load the maintained P4 matched-response implementation without writing bytecode."""
    sys.dont_write_bytecode = True
    script = Path(__file__).resolve().with_name("fit_p4_matched_response.py")
    specification = importlib.util.spec_from_file_location("p4_matched_response_fit", script)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load matched-response implementation: {script}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def optimizer_summary_path(reference: Path) -> Path:
    """Return the unique exact-optimizer summary in a matched-response experiment."""
    candidates = sorted((reference / "exact_optimizer" / "finim_outputs").glob("p4Negative*_summary.yaml"))
    if len(candidates) != 1:
        raise RuntimeError(f"expected one exact-optimizer summary under {reference}")
    return candidates[0]


def read_exact_planet(path: Path) -> dict[str, float]:
    """Read the converged exact negative-planet point estimate."""
    text = path.read_text(encoding="utf-8")
    fitted = re.search(r"^  fitted:\n(?P<body>(?:    [^\n]*\n)+)", text, flags=re.MULTILINE)
    if fitted is None or not re.search(r"^  converged: true$", text, flags=re.MULTILINE):
        raise RuntimeError(f"exact optimizer did not converge: {path}")
    values: dict[str, float] = {}
    for key in ("separation", "positionAngle", "contrast"):
        match = re.search(rf"^    {key}: ([+\-0-9.eE]+)$", fitted.group("body"), flags=re.MULTILINE)
        if match is None:
            raise RuntimeError(f"exact optimizer summary lacks fitted {key}: {path}")
        values[key] = float(match.group(1))
    if values["separation"] < 0 or values["contrast"] >= 0 or not all(map(math.isfinite, values.values())):
        raise RuntimeError(f"invalid exact optimizer point estimate: {path}")
    return {
        "separation": values["separation"],
        "position_angle": values["positionAngle"] % 360,
        "contrast": -values["contrast"],
    }


def read_response_fit(path: Path) -> dict[str, float]:
    """Read a converged matched-response position and contrast."""
    summary = json.loads(path.read_text(encoding="utf-8"))
    fit = summary.get("fit", {})
    result = {
        "separation": float(fit["separation"]),
        "position_angle": float(fit["position_angle"]) % 360,
        "contrast": float(fit["contrast"]),
    }
    if fit.get("status") != "converged" or result["separation"] < 0 or result["contrast"] <= 0:
        raise RuntimeError(f"matched-response fit did not converge to a positive source: {path}")
    if not all(map(math.isfinite, result.values())):
        raise RuntimeError(f"matched-response fit is non-finite: {path}")
    return result


def quadratic_at(response_fit, values: np.ndarray, row: float, column: float) -> float:
    """Evaluate one map with the maintained local quadratic convention."""
    center_row = int(round(row))
    center_column = int(round(column))
    coefficients = response_fit.quadratic_coefficients(values, center_column, center_row)
    return response_fit.quadratic_value(coefficients, row - center_row, column - center_column)


def response_field(response_fit, manifest: Path, mode_fraction: float) -> dict[str, object]:
    """Read one complete response field selected by exact mode fraction."""
    header = fits.getheader(manifest)
    if int(header.get("P4 PSF COMPLETE", 0)) != 1:
        raise RuntimeError(f"response manifest is incomplete: {manifest}")
    detector_rows = int(header.get("P4 PSF TEMPLATE ROWS", 0))
    detector_columns = int(header.get("P4 PSF TEMPLATE COLUMNS", 0))
    if detector_rows <= 0 or detector_columns <= 0:
        raise RuntimeError(f"response manifest lacks detector dimensions: {manifest}")
    mode_index = response_fit.selected_mode_index(header, mode_fraction, manifest)
    directory, prefix = response_fit.product_prefix(manifest)
    coordinates = response_fit.read_coordinates(directory / f"{prefix}coordinates.fits")
    model = np.asarray(fits.getdata(directory / f"{prefix}model_{mode_index:04d}.fits"), dtype=np.float64)
    validity = np.asarray(
        fits.getdata(directory / f"{prefix}validity_{mode_index:04d}.fits"), dtype=np.float64
    ).reshape(-1)
    if model.ndim != 3 or model.shape[0] != coordinates.shape[0] or validity.size != coordinates.shape[0]:
        raise RuntimeError(f"response field products have inconsistent dimensions: {manifest}")
    return {
        "header": header,
        "detector_rows": detector_rows,
        "detector_columns": detector_columns,
        "coordinates": coordinates,
        "model": model,
        "validity": validity,
        "minimum_support": float(header.get("P4 PSF FILTER MIN GOOD FRACTION", 1.0)),
    }


def amplitude_map(
    response_fit,
    science_path: Path,
    mode_fraction: float,
    field: dict[str, object],
) -> np.ndarray:
    """Apply one response field to one selected final-image plane."""
    science, header = response_fit.read_science(
        science_path,
        int(field["detector_rows"]),
        int(field["detector_columns"]),
    )
    mode_index = response_fit.selected_mode_index(header, mode_fraction, science_path)
    amplitude, _, _ = response_fit.apply_response(
        science,
        mode_index,
        np.asarray(field["model"]),
        np.asarray(field["validity"]),
        np.asarray(field["coordinates"]),
        float(field["minimum_support"]),
    )
    return amplitude


def proportionality_metrics(
    reference: dict[str, object],
    scaled: dict[str, object],
    scale: float,
) -> dict[str, float]:
    """Measure whether a scaled-template response is a scalar multiple of the reference response."""
    reference_coordinates = np.asarray(reference["coordinates"])
    scaled_coordinates = np.asarray(scaled["coordinates"])
    if not np.array_equal(reference_coordinates, scaled_coordinates):
        raise RuntimeError("reference and scaled response coordinate products differ")
    reference_model = np.asarray(reference["model"])
    scaled_model = np.asarray(scaled["model"])
    if reference_model.shape != scaled_model.shape:
        raise RuntimeError("reference and scaled response models have different dimensions")
    valid = (
        (np.asarray(reference["validity"]) > 0)
        & (np.asarray(scaled["validity"]) > 0)
    )
    retained = valid[:, None, None] & np.isfinite(reference_model) & np.isfinite(scaled_model)
    reference_values = reference_model[retained]
    scaled_values = scaled_model[retained]
    expected_values = scale * reference_values
    reference_energy = float(np.dot(reference_values, reference_values))
    scaled_energy = float(np.dot(scaled_values, scaled_values))
    expected_energy = float(np.dot(expected_values, expected_values))
    cross_energy = float(np.dot(reference_values, scaled_values))
    if reference_energy <= 0 or scaled_energy <= 0 or expected_energy <= 0:
        raise RuntimeError("response proportionality comparison has zero energy")
    realized_scale = cross_energy / reference_energy
    return {
        "finite_overlap_count": int(reference_values.size),
        "cosine_similarity": cross_energy / math.sqrt(reference_energy * scaled_energy),
        "realized_scale": realized_scale,
        "realized_scale_ratio_to_requested": realized_scale / scale,
        "relative_l2_from_requested_scale": float(np.linalg.norm(scaled_values - expected_values))
        / math.sqrt(expected_energy),
    }


def format_number(value: object) -> str:
    """Format one finite diagnostic value compactly."""
    number = float(value)
    if not math.isfinite(number):
        raise RuntimeError("sanity-check diagnostic produced a non-finite value")
    return f"{number:.8g}"


def main() -> int:
    """Write negative-injection and response-scale comparisons."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_experiment", type=Path)
    parser.add_argument("response_fit_subtracted", type=Path)
    parser.add_argument("over_subtracted", type=Path)
    parser.add_argument("scaled_manifest", type=Path)
    parser.add_argument("scaled_fit_summary", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--mode-fraction", type=float, default=0.15)
    parser.add_argument("--template-scale", type=float, required=True)
    parser.add_argument("--over-subtraction-factor", type=float, default=4.0)
    arguments = parser.parse_args()
    if arguments.template_scale <= 0 or arguments.over_subtraction_factor <= 1:
        raise RuntimeError("template scale must be positive and over-subtraction factor must exceed one")

    reference = arguments.reference_experiment.resolve()
    response_fit_subtracted = arguments.response_fit_subtracted.resolve()
    over_subtracted = arguments.over_subtracted.resolve()
    scaled_manifest = arguments.scaled_manifest.resolve()
    scaled_fit_summary = arguments.scaled_fit_summary.resolve()
    output_directory = arguments.output_directory.resolve()
    required = (
        response_fit_subtracted,
        over_subtracted,
        scaled_manifest,
        scaled_fit_summary,
        reference / "sparse_response" / "finim.fits",
        reference / "signal_free_oracle" / "finim.fits",
        reference / "sparse_response" / "finim_outputs" / "p4PSF_manifest.fits",
        reference / "sparse_fit" / "summary.json",
        reference / "sparse_fit" / "uncertainty.fits",
    )
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise RuntimeError("missing required sanity-check products: " + ", ".join(map(str, missing)))

    response_fit_module = load_response_fit_module()
    exact = read_exact_planet(optimizer_summary_path(reference))
    original_response_fit = read_response_fit(reference / "sparse_fit" / "summary.json")
    scaled_response_fit = read_response_fit(scaled_fit_summary)
    reference_manifest = reference / "sparse_response" / "finim_outputs" / "p4PSF_manifest.fits"
    reference_field = response_field(response_fit_module, reference_manifest, arguments.mode_fraction)
    scaled_field = response_field(response_fit_module, scaled_manifest, arguments.mode_fraction)
    proportionality = proportionality_metrics(reference_field, scaled_field, arguments.template_scale)

    detector_rows = int(reference_field["detector_rows"])
    detector_columns = int(reference_field["detector_columns"])
    center_row = 0.5 * (detector_rows - 1)
    center_column = 0.5 * (detector_columns - 1)
    radians = math.radians(exact["position_angle"])
    source_row = center_row - exact["separation"] * math.sin(radians)
    source_column = center_column + exact["separation"] * math.cos(radians)

    original_path = reference / "sparse_response" / "finim.fits"
    exact_subtracted_path = reference / "signal_free_oracle" / "finim.fits"
    case_specs = (
        ("original", original_path, 0.0),
        ("response_fit_subtracted", response_fit_subtracted, original_response_fit["contrast"]),
        ("exact_subtracted", exact_subtracted_path, exact["contrast"]),
        (
            "four_times_exact_subtracted",
            over_subtracted,
            arguments.over_subtraction_factor * exact["contrast"],
        ),
    )
    case_amplitudes = {
        name: quadratic_at(
            response_fit_module,
            amplitude_map(response_fit_module, path, arguments.mode_fraction, reference_field),
            source_row,
            source_column,
        )
        for name, path, _ in case_specs
    }
    exact_subtracted_amplitude = case_amplitudes["exact_subtracted"]
    uncertainty = np.asarray(fits.getdata(reference / "sparse_fit" / "uncertainty.fits"), dtype=np.float64)
    source_uncertainty = quadratic_at(response_fit_module, uncertainty, source_row, source_column)
    if source_uncertainty <= 0:
        raise RuntimeError("reference response-fit uncertainty is not positive at the exact source")

    cases: list[dict[str, object]] = []
    for name, path, subtracted_contrast in case_specs:
        delta = case_amplitudes[name] - exact_subtracted_amplitude
        cases.append(
            {
                "case": name,
                "final_image": str(path),
                "negative_injection_contrast": subtracted_contrast,
                "expected_residual_input_contrast": exact["contrast"] - subtracted_contrast,
                "fixed_position_response_amplitude": case_amplitudes[name],
                "amplitude_relative_to_exact_subtracted": delta,
                "relative_signal_to_noise": delta / source_uncertainty,
            }
        )

    scaled_physical_contrast = scaled_response_fit["contrast"] * arguments.template_scale
    result = {
        "schema": 1,
        "reference_experiment": str(reference),
        "mode_fraction": arguments.mode_fraction,
        "exact_planet": {**exact, "row": source_row, "column": source_column},
        "original_response_fit": original_response_fit,
        "negative_injection_cases": cases,
        "reference_source_uncertainty": source_uncertainty,
        "scaled_template": {
            "requested_scale": arguments.template_scale,
            "fit_in_scaled_template_units": scaled_response_fit,
            "fitted_physical_contrast": scaled_physical_contrast,
            "physical_contrast_ratio_to_original_response_fit": scaled_physical_contrast
            / original_response_fit["contrast"],
            "field_proportionality": proportionality,
        },
    }

    output_directory.mkdir(parents=True, exist_ok=True)
    json_path = output_directory / "response_scale_sanity.json"
    csv_path = output_directory / "negative_injection_cases.csv"
    markdown_path = output_directory / "response_scale_sanity.md"
    json_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cases[0]))
        writer.writeheader()
        writer.writerows(cases)

    markdown: list[str] = [
        "# P4 response scale sanity checks",
        "",
        f"- Exact negative-fit contrast: `{format_number(exact['contrast'])}`",
        f"- Original matched-response contrast: `{format_number(original_response_fit['contrast'])}`",
        f"- Fixed-position reference uncertainty: `{format_number(source_uncertainty)}`",
        "",
        "## Negative-injection bracket",
        "",
        "| Case | Negative contrast | Expected residual input contrast | "
        "Response amplitude vs exact-zero | Relative S/N |",
        "|---|---:|---:|---:|---:|",
    ]
    for case in cases:
        markdown.append(
            f"| {case['case']} | {format_number(case['negative_injection_contrast'])} | "
            f"{format_number(case['expected_residual_input_contrast'])} | "
            f"{format_number(case['amplitude_relative_to_exact_subtracted'])} | "
            f"{format_number(case['relative_signal_to_noise'])} |"
        )
    markdown.extend(
        [
            "",
            (
                "The response-fit subtraction is expected to remain positive because its injected contrast "
                "is smaller than the exact negative-fit contrast. The four-times-exact control is expected "
                "to be strongly negative."
            ),
            "",
            "## Pre-scaled response template",
            "",
            f"- Input PSF scale: `{format_number(arguments.template_scale)}`",
            f"- Fit in scaled-template units: `{format_number(scaled_response_fit['contrast'])}`",
            f"- Converted physical contrast: `{format_number(scaled_physical_contrast)}`",
            "- Converted/original response-fit contrast: "
            f"`{format_number(scaled_physical_contrast / original_response_fit['contrast'])}`",
            f"- Realized/requested response-field scale: "
            f"`{format_number(proportionality['realized_scale_ratio_to_requested'])}`",
            f"- Response-field relative L2 error from exact scaling: "
            f"`{format_number(proportionality['relative_l2_from_requested_scale'])}`",
            f"- Response-field cosine similarity: `{format_number(proportionality['cosine_similarity'])}`",
            "",
        ]
    )
    markdown_path.write_text("\n".join(markdown), encoding="utf-8")
    print(f"wrote {json_path}")
    print(f"wrote {csv_path}")
    print(f"wrote {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
