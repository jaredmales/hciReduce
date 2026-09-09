#!/usr/bin/env python3
"""Compare sparse and signal-free P4 matched-response fits with the exact negative fit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np


def read_response_fit(path: Path) -> dict[str, object]:
    """Read and minimally validate one matched-response fit summary."""
    with path.open(encoding="utf-8") as stream:
        summary = json.load(stream)
    if int(summary.get("schema", 0)) != 1 or not isinstance(summary.get("fit"), dict):
        raise RuntimeError(f"unsupported matched-response fit summary: {path}")
    return summary


def read_optimizer_fit(products: Path) -> dict[str, object]:
    """Read a converged exact point estimate and its optional complete jackknife."""
    candidates = sorted(products.glob("p4Negative*_summary.yaml"))
    if len(candidates) != 1:
        raise RuntimeError(f"expected one exact optimizer summary: {products}")
    summary_path = candidates[0]
    prefix = summary_path.name[: -len("summary.yaml")]
    jackknife_path = products / f"{prefix}jackknife.csv"
    summary_text = summary_path.read_text(encoding="utf-8")
    status_match = re.search(r'^  status: "([^"]+)"$', summary_text, flags=re.MULTILINE)
    converged_match = re.search(r"^  converged: (true|false)$", summary_text, flags=re.MULTILINE)
    dense_match = re.search(r"^  denseAgreement: (true|false)$", summary_text, flags=re.MULTILINE)
    fitted_match = re.search(
        r"^  fitted:\n(?P<body>(?:    [^\n]*\n)+)", summary_text, flags=re.MULTILINE
    )
    if (
        status_match is None
        or status_match.group(1) != "converged"
        or converged_match is None
        or converged_match.group(1) != "true"
        or dense_match is None
        or dense_match.group(1) != "true"
        or fitted_match is None
    ):
        raise RuntimeError(f"exact optimizer point fit did not converge: {summary_path}")

    fitted_values: dict[str, float] = {}
    for key in ("separation", "positionAngle", "contrast"):
        match = re.search(
            rf"^    {key}: ([+\-0-9.eE]+)$", fitted_match.group("body"), flags=re.MULTILINE
        )
        if match is None:
            raise RuntimeError(f"exact optimizer summary lacks fitted {key}: {summary_path}")
        fitted_values[key] = float(match.group(1))
    separation = fitted_values["separation"]
    position_angle = fitted_values["positionAngle"] % 360
    signed_contrast = fitted_values["contrast"]
    if (
        not all(math.isfinite(value) for value in (separation, position_angle, signed_contrast))
        or separation < 0
        or signed_contrast >= 0
    ):
        raise RuntimeError(f"exact optimizer point estimate is invalid: {summary_path}")

    evaluation_match = re.search(r"^  evaluationCount: ([0-9]+)$", summary_text, flags=re.MULTILINE)
    elapsed_match = re.search(
        r"^  evaluationElapsedSeconds: ([+\-0-9.eE]+)$", summary_text, flags=re.MULTILINE
    )
    result: dict[str, object] = {
        "status": status_match.group(1) if status_match else "unknown",
        "separation": separation,
        "position_angle": position_angle,
        "contrast": -signed_contrast,
        "signed_contrast": signed_contrast,
        "evaluation_count": int(evaluation_match.group(1)) if evaluation_match else 0,
        "evaluation_elapsed_seconds": float(elapsed_match.group(1)) if elapsed_match else math.nan,
    }

    jackknife_match = re.search(
        r"^  jackknife:\n(?P<body>(?:    [^\n]*(?:\n|$))+)", summary_text, flags=re.MULTILINE
    )
    if jackknife_match is not None:
        jackknife_body = jackknife_match.group("body")
        requested_match = re.search(r"^    requestedBlocks: ([0-9]+)$", jackknife_body, flags=re.MULTILINE)
        complete_match = re.search(r"^    complete: (true|false)$", jackknife_body, flags=re.MULTILINE)
        jackknife_status_match = re.search(r'^    status: "([^"]+)"$', jackknife_body, flags=re.MULTILINE)
        requested = int(requested_match.group(1)) if requested_match else 0
        complete = complete_match is not None and complete_match.group(1) == "true"
        jackknife: dict[str, object] = {
            "requested_blocks": requested,
            "complete": complete,
            "status": jackknife_status_match.group(1) if jackknife_status_match else "unknown",
            "sample_count": 0,
            "converged_sample_count": 0,
        }
        samples: list[dict[str, str]] = []
        converged_samples: list[dict[str, str]] = []
        if jackknife_path.is_file():
            with jackknife_path.open(encoding="utf-8", newline="") as stream:
                samples = list(csv.DictReader(stream))
            converged_samples = [sample for sample in samples if sample.get("converged") == "1"]
        jackknife["sample_count"] = len(samples)
        jackknife["converged_sample_count"] = len(converged_samples)

        if complete:
            if requested < 2 or len(samples) != requested or len(converged_samples) != requested:
                raise RuntimeError(f"completed optimizer jackknife table is inconsistent: {jackknife_path}")
            values = {
                name: np.asarray([float(sample[name]) for sample in converged_samples], dtype=np.float64)
                for name in ("row_delta", "column_delta", "separation", "position_angle", "contrast")
            }
            angle_reference = position_angle
            values["position_angle"] = angle_reference + (
                values["position_angle"] - angle_reference + 180
            ) % 360 - 180
            count = len(converged_samples)
            standard_error = {
                name: float(math.sqrt((count - 1) / count * np.sum((data - np.mean(data)) ** 2)))
                for name, data in values.items()
            }
            jackknife.update(
                {
                    "row_standard_error": standard_error["row_delta"],
                    "column_standard_error": standard_error["column_delta"],
                    "separation_standard_error": standard_error["separation"],
                    "position_angle_standard_error": standard_error["position_angle"],
                    "contrast_standard_error": standard_error["contrast"],
                }
            )
        result["jackknife"] = jackknife
    return result


def response_method(name: str, summary: dict[str, object]) -> dict[str, object]:
    """Flatten the common comparison fields from one matched-response fit."""
    fit = summary["fit"]
    if not isinstance(fit, dict):
        raise RuntimeError(f"response fit for {name} is not an object")
    curvature = fit["curvature_diagnostic"]
    if not isinstance(curvature, dict):
        raise RuntimeError(f"curvature diagnostic for {name} is not an object")
    return {
        "method": name,
        "status": str(fit["status"]),
        "separation": float(fit["separation"]),
        "position_angle": float(fit["position_angle"]),
        "contrast": float(fit["contrast"]),
        "snr": float(fit["snr"]),
        "row_standard_error": float(curvature["row_standard_error"]),
        "column_standard_error": float(curvature["column_standard_error"]),
        "separation_standard_error": float(curvature["separation_standard_error"]),
        "position_angle_standard_error": float(curvature["position_angle_standard_error"]),
        "contrast_standard_error": float(fit["contrast_standard_error"]),
        "uncertainty_method": "local matched-likelihood curvature and radial filtered-image noise",
        "spatial_model": str(summary["spatial_model"]),
        "composition": str(summary["composition"]),
    }


def optimizer_method(summary: dict[str, object]) -> dict[str, object]:
    """Flatten exact negative-optimizer values into the comparison schema."""
    jackknife = summary.get("jackknife", {})
    if not isinstance(jackknife, dict):
        raise RuntimeError("exact optimizer jackknife summary is not an object")
    requested_blocks = int(jackknife.get("requested_blocks", 0))
    jackknife_complete = bool(jackknife.get("complete", False))
    if jackknife_complete:
        uncertainty_method = "delete-one-time-block jackknife"
    elif requested_blocks:
        uncertainty_method = (
            f"incomplete delete-one-time-block jackknife "
            f"({int(jackknife.get('converged_sample_count', 0))}/{requested_blocks} converged)"
        )
    else:
        uncertainty_method = "disabled"
    return {
        "method": "exact_negative_optimizer",
        "status": str(summary["status"]),
        "separation": float(summary["separation"]),
        "position_angle": float(summary["position_angle"]),
        "contrast": float(summary["contrast"]),
        "snr": math.nan,
        "row_standard_error": float(jackknife.get("row_standard_error", math.nan)),
        "column_standard_error": float(jackknife.get("column_standard_error", math.nan)),
        "separation_standard_error": float(jackknife.get("separation_standard_error", math.nan)),
        "position_angle_standard_error": float(jackknife.get("position_angle_standard_error", math.nan)),
        "contrast_standard_error": float(jackknife.get("contrast_standard_error", math.nan)),
        "uncertainty_method": uncertainty_method,
        "spatial_model": "finite-amplitude local P4",
        "composition": "refitted each evaluation",
    }


def coordinate(method: dict[str, object], center: tuple[float, float]) -> tuple[float, float]:
    """Convert a method's polar location to P4 detector row and column."""
    separation = float(method["separation"])
    radians = math.radians(float(method["position_angle"]))
    return center[0] - separation * math.sin(radians), center[1] + separation * math.cos(radians)


def wrapped_angle_difference(left: float, right: float) -> float:
    """Return the signed shortest angular difference in degrees."""
    return (left - right + 180) % 360 - 180


def format_number(value: object) -> str:
    """Format finite comparison values while retaining missing diagnostics."""
    number = float(value)
    return f"{number:.8g}" if math.isfinite(number) else "nan"


def json_finite(value: object) -> object:
    """Replace non-finite floating-point diagnostics with JSON null values."""
    if isinstance(value, dict):
        return {str(key): json_finite(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_finite(item) for item in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    return value


def main() -> int:
    """Write machine-readable and Markdown comparisons for the complete experiment."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sparse_summary", type=Path)
    parser.add_argument("oracle_summary", type=Path)
    parser.add_argument("optimizer_products", type=Path)
    parser.add_argument("output_directory", type=Path)
    arguments = parser.parse_args()

    sparse_summary = read_response_fit(arguments.sparse_summary.resolve())
    oracle_summary = read_response_fit(arguments.oracle_summary.resolve())
    optimizer_summary = read_optimizer_fit(arguments.optimizer_products.resolve())
    if float(sparse_summary["mode_fraction"]) != float(oracle_summary["mode_fraction"]):
        raise RuntimeError("sparse and signal-free oracle fits use different mode fractions")
    if sparse_summary["science"] != oracle_summary["science"]:
        raise RuntimeError("sparse and signal-free responses were not applied to the same science product")

    methods = [
        response_method("sparse_original_response", sparse_summary),
        response_method("dense_signal_free_response", oracle_summary),
        optimizer_method(optimizer_summary),
    ]
    exact = methods[-1]
    sparse_fit = sparse_summary["fit"]
    if not isinstance(sparse_fit, dict):
        raise RuntimeError("sparse response fit is not an object")
    center = (
        float(sparse_fit["detector_center_row"]),
        float(sparse_fit["detector_center_column"]),
    )
    exact_coordinate = coordinate(exact, center)
    for method in methods:
        location = coordinate(method, center)
        method["position_difference_from_exact"] = math.hypot(
            location[0] - exact_coordinate[0], location[1] - exact_coordinate[1]
        )
        method["separation_difference_from_exact"] = float(method["separation"]) - float(exact["separation"])
        method["pa_difference_from_exact"] = wrapped_angle_difference(
            float(method["position_angle"]), float(exact["position_angle"])
        )
        method["contrast_difference_from_exact"] = float(method["contrast"]) - float(exact["contrast"])
        method["contrast_ratio_to_exact"] = float(method["contrast"]) / float(exact["contrast"])

    output_directory = arguments.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    csv_path = output_directory / "fit_comparison.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(methods[0]))
        writer.writeheader()
        writer.writerows(methods)

    json_path = output_directory / "fit_comparison.json"
    with json_path.open("w", encoding="utf-8") as stream:
        json.dump(
            json_finite(
                {
                    "schema": 1,
                    "mode_fraction": sparse_summary["mode_fraction"],
                    "science": sparse_summary["science"],
                    "methods": methods,
                    "optimizer_evaluation_count": optimizer_summary["evaluation_count"],
                    "optimizer_evaluation_elapsed_seconds": optimizer_summary["evaluation_elapsed_seconds"],
                    "optimizer_jackknife": optimizer_summary.get("jackknife"),
                }
            ),
            stream,
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        stream.write("\n")

    markdown_path = output_directory / "fit_comparison.md"
    with markdown_path.open("w", encoding="utf-8") as stream:
        stream.write("# P4 matched-response fit validation\n\n")
        stream.write(f"Selected P4 mode fraction: {float(sparse_summary['mode_fraction']):.8g}.\n\n")
        stream.write(
            "The sparse and dense signal-free response fields were applied to the same original science image. "
            "The exact negative optimizer refits finite-amplitude P4 for every trial and is the point-estimate "
            "oracle. Response-fit position errors are local likelihood-curvature diagnostics; the exact optimizer "
            "uncertainties are reported only when every requested delete-one-time-block jackknife refit "
            "converged.\n\n"
        )
        jackknife = optimizer_summary.get("jackknife", {})
        if isinstance(jackknife, dict) and int(jackknife.get("requested_blocks", 0)):
            requested = int(jackknife["requested_blocks"])
            converged = int(jackknife.get("converged_sample_count", 0))
            if bool(jackknife.get("complete", False)):
                stream.write(f"Exact-optimizer jackknife: complete ({converged}/{requested} blocks).\n\n")
            else:
                stream.write(
                    f"Exact-optimizer jackknife: incomplete ({converged}/{requested} blocks converged); "
                    "no jackknife uncertainty is reported.\n\n"
                )
        else:
            stream.write("Exact-optimizer jackknife: disabled for this point-estimate validation.\n\n")
        stream.write(
            "| Method | Status | Separation | PA | Contrast | SNR | Position delta (pix) | Contrast / exact | "
            "sigma contrast | sigma row | sigma column |\n"
        )
        stream.write("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for method in methods:
            stream.write(
                f"| {method['method']} | {method['status']} | {format_number(method['separation'])} | "
                f"{format_number(method['position_angle'])} | {format_number(method['contrast'])} | "
                f"{format_number(method['snr'])} | {format_number(method['position_difference_from_exact'])} | "
                f"{format_number(method['contrast_ratio_to_exact'])} | "
                f"{format_number(method['contrast_standard_error'])} | "
                f"{format_number(method['row_standard_error'])} | "
                f"{format_number(method['column_standard_error'])} |\n"
            )
        stream.write("\nDetailed fit surfaces and FITS diagnostics are in the two response-fit directories.\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {markdown_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
