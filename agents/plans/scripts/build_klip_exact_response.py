#!/usr/bin/env python3
"""Package one exact-pixel paired KLIP response for external hciAnalyze filtering."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

from compare_klip_finite_response import extract_stamp, mode_counts, read_cube


def exact_coordinates(
    separation: float,
    position_angle: float,
    rows: int,
    columns: int,
) -> tuple[int, int, float, float]:
    """Resolve polar coordinates that must land on one exact integer detector pixel."""
    center_row = 0.5 * (rows - 1)
    center_column = 0.5 * (columns - 1)
    angle = math.radians(position_angle)
    source_row = center_row - separation * math.sin(angle)
    source_column = center_column + separation * math.cos(angle)
    target_row = int(math.floor(source_row + 0.5))
    target_column = int(math.floor(source_column + 0.5))
    if not math.isclose(source_row, target_row, rel_tol=0, abs_tol=2e-6) or not math.isclose(
        source_column, target_column, rel_tol=0, abs_tol=2e-6
    ):
        raise RuntimeError(
            "configured exact-response coordinates do not land on an integer detector pixel: "
            f"({source_row:.12g}, {source_column:.12g})"
        )
    delta_row = target_row - center_row
    delta_column = target_column - center_column
    return target_row, target_column, math.hypot(delta_row, delta_column), math.atan2(delta_row, delta_column)


def product_header(
    reference: fits.Header,
    role: str,
    mode_count: int | None,
    stamp_size: int,
    radius: float,
    angle: float,
    row: int,
    column: int,
    perturbation: float,
    minimum_support: float,
) -> fits.Header:
    """Return exact-response provenance compatible with the KLIP manifest consumer."""
    header = reference.copy()
    header["HIERARCH KLIP PSF STAMP SIZE"] = (stamp_size, "square response-stamp size")
    header["HIERARCH KLIP PSF SAMPLE RADII"] = (f"{radius:.12g}", "exact response sample radius")
    header["HIERARCH KLIP PSF SAMPLES PER RADIUS"] = (1, "one exact detector-pixel sample")
    header["HIERARCH KLIP PSF REQUESTED SAMPLES PER RADIUS"] = ("1", "count")
    header["HIERARCH KLIP PSF SPATIAL MODEL"] = (
        "EXACT_AZIMUTHAL",
        "model",
    )
    header["HIERARCH KLIP PSF SAMPLE ANGLE"] = (angle, "radians")
    header["HIERARCH KLIP PSF EXACT ROW"] = (row, "integer final-image anchor row")
    header["HIERARCH KLIP PSF EXACT COLUMN"] = (column, "integer final-image anchor column")
    header["HIERARCH KLIP PSF RESPONSE METHOD"] = ("refitDifference", "method")
    header["HIERARCH KLIP PSF SAMPLE AVOID RADIUS"] = (0.0, "pixels")
    header["HIERARCH KLIP PSF REFIT CONTRAST"] = (perturbation, "contrast")
    header["HIERARCH KLIP PSF REFIT TRIAL COUNT"] = (2, "positive and negative KLIP reductions")
    header["HIERARCH KLIP PSF MEASUREMENT COUNT"] = (1, "exact detector-pixel measurements")
    header["HIERARCH KLIP PSF ACCUMULATION"] = (
        "PAIRED_FINAL_DIFFERENCE",
        "response accumulation storage",
    )
    header["HIERARCH KLIP PSF RETAINED BYTES"] = (0, "external pair")
    header["HIERARCH KLIP PSF FILTER"] = (0, "external filtering only")
    header["HIERARCH KLIP PSF FILTER MIN GOOD FRACTION"] = (
        minimum_support,
        "fraction",
    )
    header["HIERARCH KLIP PSF PRODUCT SCHEMA"] = (1, "sparse response-product schema")
    header["HIERARCH KLIP PSF PRODUCT"] = (role, "product role")
    if mode_count is None:
        header["HIERARCH KLIP PSF COMPLETE"] = (1, "complete product set available")
        if "KLIP PSF MODE COUNT" in header:
            del header["KLIP PSF MODE COUNT"]
    else:
        header["HIERARCH KLIP PSF MODE COUNT"] = (mode_count, "requested retained KL modes")
        if "KLIP PSF COMPLETE" in header:
            del header["KLIP PSF COMPLETE"]
    return header


def main() -> int:
    """Build exact response, validity, manifest, and machine-readable summary products."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("science_file", type=Path)
    parser.add_argument("plus_file", type=Path)
    parser.add_argument("minus_file", type=Path)
    parser.add_argument("reference_manifest", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--separation", type=float, required=True)
    parser.add_argument("--pa", type=float, required=True)
    parser.add_argument("--perturbation", type=float, required=True)
    parser.add_argument("--stamp-size", type=int, default=11)
    parser.add_argument("--minimum-support", type=float, default=1.0)
    parser.add_argument("--output-prefix", default="klipExact_")
    arguments = parser.parse_args()

    paths = (
        arguments.science_file,
        arguments.plus_file,
        arguments.minus_file,
        arguments.reference_manifest,
    )
    for path in paths:
        if not path.is_file():
            raise RuntimeError(f"required exact-response input is missing: {path}")
    if (
        not math.isfinite(arguments.separation)
        or arguments.separation < 0
        or not math.isfinite(arguments.pa)
        or not math.isfinite(arguments.perturbation)
        or arguments.perturbation <= 0
        or arguments.stamp_size <= 0
        or arguments.stamp_size % 2 == 0
        or not math.isfinite(arguments.minimum_support)
        or not 0 <= arguments.minimum_support <= 1
        or not arguments.output_prefix
    ):
        raise RuntimeError("exact-response geometry, perturbation, support, or output prefix is invalid")

    science, science_header = read_cube(arguments.science_file)
    plus, plus_header = read_cube(arguments.plus_file)
    minus, minus_header = read_cube(arguments.minus_file)
    science_modes = mode_counts(science_header, arguments.science_file)
    plus_modes = mode_counts(plus_header, arguments.plus_file)
    minus_modes = mode_counts(minus_header, arguments.minus_file)
    if (
        science.shape != plus.shape
        or science.shape != minus.shape
        or science_modes != plus_modes
        or science_modes != minus_modes
    ):
        raise RuntimeError("science and paired exact-response cubes have inconsistent dimensions or modes")
    if science.ndim != 3 or len(science_modes) != science.shape[0]:
        raise RuntimeError("exact-response input must be a mode-labelled three-dimensional cube")
    detector_columns, detector_rows = science.shape[1:]
    target_row, target_column, target_radius, target_angle = exact_coordinates(
        arguments.separation,
        arguments.pa,
        detector_rows,
        detector_columns,
    )

    reference_header = fits.getheader(arguments.reference_manifest)
    if int(reference_header.get("KLIP PSF COMPLETE", 0)) != 1 or str(
        reference_header.get("KLIP PSF PRODUCT", "")
    ).strip() != "MANIFEST":
        raise RuntimeError("reference response manifest is incomplete")
    reference_modes = mode_counts(reference_header, arguments.reference_manifest)
    if reference_modes != science_modes:
        raise RuntimeError("reference manifest and science cube have different KL mode labels")
    for keyword in ("REGMINR", "REGMAXR"):
        if keyword not in science_header or keyword not in reference_header:
            raise RuntimeError(f"science and reference manifest must provide {keyword}")
        if str(science_header[keyword]).strip() != str(reference_header[keyword]).strip():
            raise RuntimeError(f"science and reference manifest have different {keyword}")

    arguments.output_directory.mkdir(parents=True, exist_ok=False)
    manifest_header = product_header(
        reference_header,
        "MANIFEST",
        None,
        arguments.stamp_size,
        target_radius,
        target_angle,
        target_row,
        target_column,
        arguments.perturbation,
        arguments.minimum_support,
    )
    manifest_path = arguments.output_directory / f"{arguments.output_prefix}manifest.fits"
    fits.writeto(manifest_path, np.ones((1, 1), dtype=np.float32), header=manifest_header, overwrite=False)

    mode_products = []
    for mode_index, retained_modes in enumerate(science_modes):
        plus_stamp = extract_stamp(plus[mode_index], target_row, target_column, arguments.stamp_size)
        minus_stamp = extract_stamp(minus[mode_index], target_row, target_column, arguments.stamp_size)
        validity = np.isfinite(plus_stamp) & np.isfinite(minus_stamp)
        response = np.zeros_like(plus_stamp, dtype=np.float64)
        response[validity] = (plus_stamp[validity] - minus_stamp[validity]) / (2 * arguments.perturbation)
        if not validity[arguments.stamp_size // 2, arguments.stamp_size // 2]:
            raise RuntimeError(f"mode {retained_modes} has an invalid exact-response anchor")
        response_values = response[validity]
        if response_values.size == 0 or float(np.dot(response_values, response_values)) <= 0:
            raise RuntimeError(f"mode {retained_modes} exact response has no positive energy")

        response_header = product_header(
            reference_header,
            "RADIAL_RESPONSE",
            retained_modes,
            arguments.stamp_size,
            target_radius,
            target_angle,
            target_row,
            target_column,
            arguments.perturbation,
            arguments.minimum_support,
        )
        validity_header = product_header(
            reference_header,
            "RADIAL_VALIDITY",
            retained_modes,
            arguments.stamp_size,
            target_radius,
            target_angle,
            target_row,
            target_column,
            arguments.perturbation,
            arguments.minimum_support,
        )
        mode_tag = f"mode{mode_index:03d}"
        response_path = arguments.output_directory / f"{arguments.output_prefix}{mode_tag}_radial_response.fits"
        validity_path = arguments.output_directory / f"{arguments.output_prefix}{mode_tag}_radial_validity.fits"
        fits.writeto(
            response_path,
            np.asarray([response.T], dtype=np.float32),
            header=response_header,
            overwrite=False,
        )
        fits.writeto(
            validity_path,
            np.asarray([validity.T], dtype=np.float32),
            header=validity_header,
            overwrite=False,
        )
        mode_products.append(
            {
                "plane": mode_index,
                "mode_count": retained_modes,
                "valid_samples": int(np.count_nonzero(validity)),
                "response_energy": float(np.dot(response_values, response_values)),
                "response_file": response_path.name,
                "validity_file": validity_path.name,
            }
        )

    summary = {
        "schema": 1,
        "science_file": str(arguments.science_file.resolve()),
        "plus_file": str(arguments.plus_file.resolve()),
        "minus_file": str(arguments.minus_file.resolve()),
        "reference_manifest": str(arguments.reference_manifest.resolve()),
        "manifest": manifest_path.name,
        "spatial_model": "EXACT_AZIMUTHAL",
        "target_row": target_row,
        "target_column": target_column,
        "target_radius": target_radius,
        "target_angle_radians": target_angle,
        "separation": arguments.separation,
        "position_angle": arguments.pa % 360,
        "perturbation": arguments.perturbation,
        "stamp_size": arguments.stamp_size,
        "minimum_support": arguments.minimum_support,
        "modes": mode_products,
    }
    summary_path = arguments.output_directory / "exact_response_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {manifest_path}")
    print(
        f"Exact anchor: row={target_row} column={target_column} "
        f"radius={target_radius:.12g} angle={target_angle:.12g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
