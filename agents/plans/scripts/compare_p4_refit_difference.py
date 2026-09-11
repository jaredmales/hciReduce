#!/usr/bin/env python3
"""Compare a P4 paired-refit response fit with frozen and exact fits."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

from astropy.io import fits


def read_response_fit(path: Path) -> dict[str, float]:
    """Read a converged matched-response fit summary."""
    summary = json.loads(path.read_text(encoding="utf-8"))
    fit = summary.get("fit", {})
    result = {
        "separation": float(fit["separation"]),
        "position_angle": float(fit["position_angle"]) % 360,
        "contrast": float(fit["contrast"]),
        "contrast_standard_error": float(fit["contrast_standard_error"]),
        "snr": float(fit["snr"]),
    }
    if fit.get("status") != "converged" or result["separation"] < 0 or result["contrast"] <= 0:
        raise RuntimeError(f"matched-response fit did not converge: {path}")
    if not all(math.isfinite(value) for value in result.values()):
        raise RuntimeError(f"matched-response fit contains non-finite values: {path}")
    return result


def optimizer_summary_path(reference: Path) -> Path:
    """Return the unique exact-optimizer summary in a matched-response experiment."""
    candidates = sorted((reference / "exact_optimizer" / "finim_outputs").glob("p4Negative*_summary.yaml"))
    if len(candidates) != 1:
        raise RuntimeError(f"expected one exact optimizer summary under {reference}")
    return candidates[0]


def read_exact_fit(path: Path) -> dict[str, float]:
    """Read the converged exact negative-planet point estimate."""
    text = path.read_text(encoding="utf-8")
    fitted = re.search(r"^  fitted:\n(?P<body>(?:    [^\n]*\n)+)", text, flags=re.MULTILINE)
    if fitted is None or not re.search(r"^  converged: true$", text, flags=re.MULTILINE):
        raise RuntimeError(f"exact optimizer did not converge: {path}")
    values: dict[str, float] = {}
    for key in ("separation", "positionAngle", "contrast"):
        match = re.search(rf"^    {key}: ([+\-0-9.eE]+)$", fitted.group("body"), flags=re.MULTILINE)
        if match is None:
            raise RuntimeError(f"exact optimizer summary lacks {key}: {path}")
        values[key] = float(match.group(1))
    result = {
        "separation": values["separation"],
        "position_angle": values["positionAngle"] % 360,
        "contrast": -values["contrast"],
    }
    if result["separation"] < 0 or result["contrast"] <= 0 or not all(map(math.isfinite, result.values())):
        raise RuntimeError(f"exact optimizer point is invalid: {path}")
    return result


def resource_usage(path: Path) -> dict[str, float]:
    """Read optional key-value resource measurements."""
    result: dict[str, float] = {}
    if not path.is_file():
        return result
    for line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            result[key] = float(value)
    return result


def position_difference(left: dict[str, float], right: dict[str, float]) -> float:
    """Return detector-coordinate separation between two polar fits."""
    left_angle = math.radians(left["position_angle"])
    right_angle = math.radians(right["position_angle"])
    left_row = -left["separation"] * math.sin(left_angle)
    left_column = left["separation"] * math.cos(left_angle)
    right_row = -right["separation"] * math.sin(right_angle)
    right_column = right["separation"] * math.cos(right_angle)
    return math.hypot(left_row - right_row, left_column - right_column)


def format_number(value: float | int) -> str:
    """Format one finite diagnostic compactly."""
    return f"{value:.8g}" if isinstance(value, float) else str(value)


def main() -> int:
    """Write JSON, CSV, and Markdown comparisons."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference_experiment", type=Path)
    parser.add_argument("refit_manifest", type=Path)
    parser.add_argument("refit_fit_summary", type=Path)
    parser.add_argument("output_directory", type=Path)
    arguments = parser.parse_args()

    reference = arguments.reference_experiment.resolve()
    manifest = arguments.refit_manifest.resolve()
    fit_summary = arguments.refit_fit_summary.resolve()
    output = arguments.output_directory.resolve()
    frozen_summary = reference / "sparse_fit" / "summary.json"
    required = (manifest, fit_summary, frozen_summary, optimizer_summary_path(reference))
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise RuntimeError("missing comparison products: " + ", ".join(map(str, missing)))

    exact = read_exact_fit(required[-1])
    frozen = read_response_fit(frozen_summary)
    refit = read_response_fit(fit_summary)
    header = fits.getheader(manifest)
    if int(header.get("P4 PSF PRODUCT SCHEMA", 0)) != 7:
        raise RuntimeError(f"refit response manifest has the wrong schema: {manifest}")
    if str(header.get("P4 PSF RESPONSE", "")).strip() != "REFIT_CENTRAL_DIFFERENCE":
        raise RuntimeError(f"manifest does not describe a refit-difference response: {manifest}")

    methods = []
    for name, fit in (("frozen_sparse_response", frozen), ("refit_difference_response", refit)):
        methods.append(
            {
                "method": name,
                **fit,
                "contrast_ratio_to_exact": fit["contrast"] / exact["contrast"],
                "contrast_fractional_difference_from_exact": (fit["contrast"] - exact["contrast"])
                / exact["contrast"],
                "position_difference_from_exact_pixels": position_difference(fit, exact),
            }
        )
    methods.append(
        {
            "method": "exact_negative_optimizer",
            **exact,
            "contrast_standard_error": math.nan,
            "snr": math.nan,
            "contrast_ratio_to_exact": 1.0,
            "contrast_fractional_difference_from_exact": 0.0,
            "position_difference_from_exact_pixels": 0.0,
        }
    )
    response_usage = resource_usage(manifest.parents[1] / "resource_usage.txt")
    response = {
        "schema": 1,
        "reference_experiment": str(reference),
        "manifest": str(manifest),
        "refit_half_amplitude": float(header["P4 PSF REFIT CONTRAST"]),
        "paired_detector_fit_count": int(header["P4 PSF REFIT FIT COUNT"]),
        "measurement_count": int(header["P4 PSF MEASUREMENT COUNT"]),
        "excluded_candidate_count": int(header["P4 PSF SAMPLE EXCLUDED COUNT"]),
        "response_resource_usage": response_usage,
        "methods": methods,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "refit_difference_comparison.json").write_text(
        json.dumps(response, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    columns = tuple(methods[0])
    with (output / "refit_difference_comparison.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(methods)

    lines = [
        "# P4 paired-refit response comparison",
        "",
        f"- Central-difference half-amplitude: {format_number(response['refit_half_amplitude'])}",
        f"- Sparse response measurements: {response['measurement_count']}",
        f"- Positive-plus-negative detector fits: {response['paired_detector_fit_count']}",
        f"- Avoided detector candidates: {response['excluded_candidate_count']}",
    ]
    if "wall_seconds" in response_usage:
        lines.append(f"- Response-stage wall time: {format_number(response_usage['wall_seconds'])} s")
    lines.extend(
        [
            "",
            "| method | separation | PA (deg) | contrast | ratio to exact | fractional contrast difference | position difference (px) | SNR |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for method in methods:
        lines.append(
            "| {method} | {separation} | {position_angle} | {contrast} | {contrast_ratio_to_exact} | "
            "{contrast_fractional_difference_from_exact} | {position_difference_from_exact_pixels} | {snr} |".format(
                **{key: format_number(value) if isinstance(value, (float, int)) else value for key, value in method.items()}
            )
        )
    (output / "refit_difference_comparison.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
