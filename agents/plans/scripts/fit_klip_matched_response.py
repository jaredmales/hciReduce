#!/usr/bin/env python3
"""Fit one source from production KLIP normalized-response filter products."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

from fit_p4_matched_response import fit_response, radial_noise_maps


def product_with_role(case_directory: Path, role: str) -> Path:
    """Return the unique KLIP PSF product declaring the requested role."""
    candidates: list[Path] = []
    for directory in (case_directory, *case_directory.glob("*_outputs")):
        if not directory.is_dir():
            continue
        for path in directory.glob("*.fits"):
            try:
                declared_role = str(fits.getheader(path).get("KLIP PSF PRODUCT", "")).strip()
            except OSError:
                continue
            if declared_role == role:
                candidates.append(path)
    if len(candidates) != 1:
        raise RuntimeError(
            f"expected one KLIP PSF {role} product in {case_directory}, found {len(candidates)}"
        )
    return candidates[0]


def read_cube(path: Path) -> tuple[np.ndarray, fits.Header]:
    """Read a two- or three-dimensional FITS product as a mode-major cube."""
    header = fits.getheader(path)
    cube = np.asarray(fits.getdata(path), dtype=np.float64)
    if cube.ndim == 2:
        cube = cube[np.newaxis, :, :]
    if cube.ndim != 3:
        raise RuntimeError(f"KLIP filter product must be two- or three-dimensional: {path}")
    return cube, header


def mode_index(header: fits.Header, requested_mode: int, path: Path) -> int:
    """Resolve an exact retained-mode count from the ordered NMODES header vector."""
    try:
        modes = [int(token) for token in str(header["NMODES"]).split(",")]
    except (KeyError, ValueError) as error:
        raise RuntimeError(f"KLIP filter product has no valid NMODES vector: {path}") from error
    matches = [index for index, mode in enumerate(modes) if mode == requested_mode]
    if len(matches) != 1:
        raise RuntimeError(f"KL mode count {requested_mode} is not exact and unique in {path}")
    return matches[0]


def write_image(path: Path, data: np.ndarray, role: str, mode_count: int) -> None:
    """Write one KLIP response-fit diagnostic image with compact provenance."""
    header = fits.Header()
    header["HIERARCH KLIP RESPONSE FIT PRODUCT"] = role
    header["HIERARCH KLIP RESPONSE FIT MODE COUNT"] = mode_count
    fits.writeto(path, np.asarray(data, dtype=np.float32), header=header, overwrite=False)


def main() -> int:
    """Fit a bounded local likelihood peak from persisted KLIP filter products."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path, help="completed klipReduce output directory")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--mode-count", type=int, default=200)
    parser.add_argument("--initial-separation", type=float, required=True)
    parser.add_argument("--initial-pa", type=float, required=True)
    parser.add_argument("--position-bound", type=float, default=1.0)
    parser.add_argument("--noise-exclusion-radius", type=float, default=5.0)
    parser.add_argument("--noise-min-radius", type=float, default=6.0)
    parser.add_argument("--noise-max-radius", type=float, default=60.0)
    parser.add_argument("--lambda-d", type=float, default=3.6)
    arguments = parser.parse_args()

    for name in (
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
    if arguments.mode_count <= 0:
        raise RuntimeError("--mode-count must be positive")
    if arguments.initial_separation < 0 or arguments.position_bound <= 0 or arguments.noise_exclusion_radius < 0:
        raise RuntimeError("separation/exclusion must be nonnegative and position bound must be positive")
    if arguments.noise_min_radius < 0 or arguments.noise_min_radius >= arguments.noise_max_radius:
        raise RuntimeError("noise radii must be nonnegative and strictly ordered")
    if arguments.lambda_d <= 0:
        raise RuntimeError("--lambda-d must be positive")

    case_directory = arguments.case_directory.resolve()
    output_directory = arguments.output_directory.resolve()
    if not case_directory.is_dir():
        raise RuntimeError(f"KLIP case directory is not readable: {case_directory}")
    output_directory.mkdir(parents=True, exist_ok=True)
    managed_outputs = {"signal.fits", "uncertainty.fits", "snr.fits", "summary.json", "surface.csv"}
    existing_outputs = sorted(name for name in managed_outputs if (output_directory / name).exists())
    if existing_outputs:
        raise RuntimeError(
            f"output directory already contains KLIP fit products: {', '.join(existing_outputs)}"
        )

    paths = {
        role: product_with_role(case_directory, role)
        for role in ("FILTERED", "FILTER_NORMALIZATION", "FILTER_SUPPORT", "FILTER_VALIDITY", "MANIFEST")
    }
    manifest_header = fits.getheader(paths["MANIFEST"])
    if int(manifest_header.get("KLIP PSF COMPLETE", 0)) != 1:
        raise RuntimeError(f"KLIP response manifest is incomplete: {paths['MANIFEST']}")
    if int(manifest_header.get("KLIP PSF PRODUCT SCHEMA", 0)) != 1:
        raise RuntimeError(f"unsupported KLIP response schema: {paths['MANIFEST']}")
    if int(manifest_header.get("KLIP PSF FILTER", 0)) != 1:
        raise RuntimeError(f"KLIP response manifest does not declare normalized filtering: {paths['MANIFEST']}")

    amplitude_cube, amplitude_header = read_cube(paths["FILTERED"])
    normalization_cube, _ = read_cube(paths["FILTER_NORMALIZATION"])
    support_cube, _ = read_cube(paths["FILTER_SUPPORT"])
    validity_cube, _ = read_cube(paths["FILTER_VALIDITY"])
    if not (
        amplitude_cube.shape
        == normalization_cube.shape
        == support_cube.shape
        == validity_cube.shape
    ):
        raise RuntimeError("KLIP filter amplitude, normalization, support, and validity dimensions differ")

    selected_mode = mode_index(amplitude_header, arguments.mode_count, paths["FILTERED"])
    valid = np.isfinite(validity_cube[selected_mode]) & (validity_cube[selected_mode] > 0.5)
    amplitude = np.where(valid, amplitude_cube[selected_mode], np.nan)
    normalization = np.where(valid, normalization_cube[selected_mode], np.nan)
    support = np.where(valid, support_cube[selected_mode], np.nan)
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
    write_image(output_directory / "signal.fits", signal, "RADIAL_MEAN_SUBTRACTED_AMPLITUDE", arguments.mode_count)
    write_image(output_directory / "uncertainty.fits", uncertainty, "RADIAL_NOISE_UNCERTAINTY", arguments.mode_count)
    write_image(output_directory / "snr.fits", snr, "SMALL_SAMPLE_CORRECTED_SNR", arguments.mode_count)

    summary = {
        "schema": 1,
        "case_directory": str(case_directory),
        "mode_count": arguments.mode_count,
        "mode_index": selected_mode,
        "manifest": str(paths["MANIFEST"]),
        "filtered": str(paths["FILTERED"]),
        "normalization": str(paths["FILTER_NORMALIZATION"]),
        "support": str(paths["FILTER_SUPPORT"]),
        "validity": str(paths["FILTER_VALIDITY"]),
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
        "KLIP matched-response fit: "
        f"Nmodes={arguments.mode_count}, sep={fit['separation']:.12g}, "
        f"PA={fit['position_angle']:.12g}, contrast={fit['contrast']:.12g}, "
        f"SNR={fit['snr']:.8g}, status={fit['status']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
