#!/usr/bin/env python3
"""Fit a KLIP negative companion with repeated end-to-end reductions."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shlex
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.optimize import minimize


def polar_from_offset(row_offset: float, column_offset: float) -> tuple[float, float]:
    """Convert a detector offset to separation and PA east of north."""
    separation = math.hypot(row_offset, column_offset)
    position_angle = math.degrees(math.atan2(-row_offset, column_offset)) % 360
    return separation, position_angle


def header_vector(header: fits.Header, keyword: str) -> list[float]:
    """Read one comma-separated numeric FITS-header vector."""
    try:
        return [float(token) for token in str(header[keyword]).split(",")]
    except (KeyError, ValueError) as error:
        raise RuntimeError(f"KLIP optimizer output has no valid {keyword} vector") from error


def read_final_image(path: Path, mode_count: int) -> tuple[np.ndarray, fits.Header]:
    """Read and validate one single-mode optimizer result."""
    header = fits.getheader(path)
    image = np.asarray(fits.getdata(path), dtype=np.float64)
    if image.ndim == 3 and image.shape[0] == 1:
        image = image[0]
    if image.ndim != 2 or image.size == 0:
        raise RuntimeError(f"KLIP optimizer result is not one nonempty image: {path}")
    modes = header_vector(header, "NMODES")
    if modes != [float(mode_count)]:
        raise RuntimeError(f"KLIP optimizer result has unexpected mode vector: {modes}")
    return image, header


def aperture_merit(
    image: np.ndarray,
    aperture_row: float,
    aperture_column: float,
    aperture_radius: float,
) -> tuple[float, int]:
    """Calculate uniform mean-square residual in one fixed sky aperture."""
    detector_columns, detector_rows = image.shape
    column_indices, row_indices = np.indices(image.shape, dtype=np.float64)
    selected = (
        (row_indices - aperture_row) ** 2 + (column_indices - aperture_column) ** 2
        <= aperture_radius**2
    )
    values = image[selected]
    if values.size == 0:
        raise RuntimeError("KLIP negative-planet aperture contains no pixel centers")
    if not np.isfinite(values).all():
        raise RuntimeError("KLIP negative-planet aperture contains invalid output pixels")
    if (
        aperture_row - aperture_radius < -0.5
        or aperture_row + aperture_radius > detector_rows - 0.5
        or aperture_column - aperture_radius < -0.5
        or aperture_column + aperture_radius > detector_columns - 0.5
    ):
        raise RuntimeError("KLIP negative-planet aperture crosses the final-image boundary")
    merit = float(np.mean(values * values))
    if not math.isfinite(merit):
        raise RuntimeError("KLIP negative-planet merit is non-finite")
    return merit, int(values.size)


def main() -> int:
    """Run a bounded joint position/contrast negative-companion fit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base_config", type=Path)
    parser.add_argument("psf_file", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--klipreduce-bin", default="klipReduce")
    parser.add_argument("--mode-count", type=int, default=200)
    parser.add_argument("--initial-separation", type=float, required=True)
    parser.add_argument("--initial-pa", type=float, required=True)
    parser.add_argument("--initial-contrast", type=float, required=True)
    parser.add_argument("--position-bound", type=float, default=1.0)
    parser.add_argument("--contrast-lower", type=float, default=0.0)
    parser.add_argument("--contrast-upper", type=float, default=0.01)
    parser.add_argument("--aperture-radius", type=float, default=5.0)
    parser.add_argument("--maximum-evaluations", type=int, default=192)
    parser.add_argument("--parameter-tolerance", type=float, default=5e-4)
    parser.add_argument("--merit-tolerance", type=float, default=1e-5)
    arguments = parser.parse_args()

    base_config = arguments.base_config.resolve()
    psf_file = arguments.psf_file.resolve()
    output_directory = arguments.output_directory.resolve()
    klipreduce_path = shutil.which(arguments.klipreduce_bin)
    if not base_config.is_file() or not psf_file.is_file():
        raise RuntimeError("KLIP optimizer base configuration or PSF file is missing")
    if klipreduce_path is None:
        raise RuntimeError(f"klipReduce executable was not found: {arguments.klipreduce_bin}")
    numeric_values = (
        arguments.initial_separation,
        arguments.initial_pa,
        arguments.initial_contrast,
        arguments.position_bound,
        arguments.contrast_lower,
        arguments.contrast_upper,
        arguments.aperture_radius,
        arguments.parameter_tolerance,
        arguments.merit_tolerance,
    )
    if not all(math.isfinite(value) for value in numeric_values):
        raise RuntimeError("KLIP optimizer numeric controls must be finite")
    if (
        arguments.mode_count <= 0
        or arguments.initial_separation < 0
        or arguments.initial_contrast <= 0
        or arguments.position_bound <= 0
        or arguments.contrast_lower < 0
        or arguments.contrast_upper <= arguments.contrast_lower
        or not arguments.contrast_lower < arguments.initial_contrast < arguments.contrast_upper
        or arguments.aperture_radius <= 0
        or arguments.maximum_evaluations < 16
        or arguments.parameter_tolerance <= 0
        or arguments.merit_tolerance <= 0
    ):
        raise RuntimeError("KLIP optimizer controls do not define a valid bounded fit")

    settings = {
        "schema": 1,
        "base_config": str(base_config),
        "psf_file": str(psf_file),
        "klipreduce": str(Path(klipreduce_path).resolve()),
        "mode_count": arguments.mode_count,
        "initial_separation": arguments.initial_separation,
        "initial_position_angle": arguments.initial_pa % 360,
        "initial_contrast": arguments.initial_contrast,
        "position_bound": arguments.position_bound,
        "contrast_lower": arguments.contrast_lower,
        "contrast_upper": arguments.contrast_upper,
        "aperture_radius": arguments.aperture_radius,
        "maximum_evaluations": arguments.maximum_evaluations,
        "parameter_tolerance": arguments.parameter_tolerance,
        "merit_tolerance": arguments.merit_tolerance,
    }
    settings_path = output_directory / "settings.json"
    summary_path = output_directory / "summary.json"
    if settings_path.is_file():
        existing_settings = json.loads(settings_path.read_text(encoding="utf-8"))
        if existing_settings != settings:
            raise RuntimeError(f"existing KLIP optimizer settings differ: {settings_path}")
        if summary_path.is_file():
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            if (
                int(summary.get("schema", 0)) == 1
                and summary.get("complete") is True
                and summary.get("status") == "converged"
            ):
                print(f"Completed KLIP negative-planet fit exists; skipping: {summary_path}")
                return 0
            raise RuntimeError(f"existing KLIP negative-planet fit is not converged: {summary_path}")
    else:
        output_directory.mkdir(parents=True, exist_ok=True)
        settings_path.write_text(json.dumps(settings, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    evaluations_directory = output_directory / "evaluations"
    evaluations_directory.mkdir(exist_ok=True)
    evaluations: list[dict[str, object]] = []
    cache: dict[tuple[str, str, str], dict[str, object]] = {}
    for evaluation_path in sorted(evaluations_directory.glob("eval_*/evaluation.json")):
        evaluation = json.loads(evaluation_path.read_text(encoding="utf-8"))
        if evaluation.get("complete") is not True:
            continue
        point = tuple(f"{float(value):.17g}" for value in evaluation["normalized_point"])
        if len(point) != 3 or point in cache:
            raise RuntimeError(f"invalid or duplicate cached optimizer point: {evaluation_path}")
        cache[point] = evaluation
        evaluations.append(evaluation)

    pa_radians = math.radians(arguments.initial_pa)
    initial_row_offset = -arguments.initial_separation * math.sin(pa_radians)
    initial_column_offset = arguments.initial_separation * math.cos(pa_radians)
    contrast_range = arguments.contrast_upper - arguments.contrast_lower
    aperture_coordinates: tuple[float, float] | None = None
    expected_shape: tuple[int, int] | None = None
    merit_scale: float | None = None
    if evaluations:
        cached_image, _ = read_final_image(Path(str(evaluations[0]["final_image"])), arguments.mode_count)
        expected_shape = cached_image.shape
        detector_columns, detector_rows = cached_image.shape
        aperture_coordinates = (
            0.5 * (detector_rows - 1) + initial_row_offset,
            0.5 * (detector_columns - 1) + initial_column_offset,
        )

    def physical_point(normalized_point: np.ndarray) -> tuple[float, float, float, float, float]:
        """Map normalized bounded coordinates to sky geometry and positive contrast."""
        row_delta = float(normalized_point[0]) * arguments.position_bound
        column_delta = float(normalized_point[1]) * arguments.position_bound
        contrast = arguments.contrast_lower + 0.5 * (float(normalized_point[2]) + 1) * contrast_range
        separation, position_angle = polar_from_offset(
            initial_row_offset + row_delta,
            initial_column_offset + column_delta,
        )
        return row_delta, column_delta, separation, position_angle, contrast

    def evaluate(normalized_point: np.ndarray) -> float:
        """Run or reuse one KLIP reduction and return its scaled aperture merit."""
        nonlocal aperture_coordinates, expected_shape, merit_scale
        bounded = np.clip(np.asarray(normalized_point, dtype=np.float64), -1, 1)
        key = tuple(f"{float(value):.17g}" for value in bounded)
        row_delta, column_delta, separation, position_angle, contrast = physical_point(bounded)
        if key in cache:
            evaluation = cache[key]
            merit = float(evaluation["merit"])
            if merit_scale is None:
                merit_scale = merit
            return merit / merit_scale

        evaluation_index = len(list(evaluations_directory.glob("eval_*")))
        evaluation_directory = evaluations_directory / f"eval_{evaluation_index:04d}"
        evaluation_directory.mkdir(exist_ok=False)
        final_path = evaluation_directory / "finim.fits"
        command = [
            str(Path(klipreduce_path).resolve()),
            "--config",
            str(base_config),
            "--klip.Nmodes",
            str(arguments.mode_count),
            "--psfResponse.file",
            "",
            "--psfResponse.outputModels=false",
            "--psfResponse.filter=false",
            "--planet.sep",
            f"{separation:.17g}",
            "--planet.PA",
            f"{position_angle:.17g}",
            "--planet.contrast",
            f"{contrast:.17g}",
            "--fake.method",
            "single",
            "--fake.fileName",
            str(psf_file),
            "--fake.sep",
            f"{separation:.17g}",
            "--fake.PA",
            f"{position_angle:.17g}",
            "--fake.contrast",
            f"{-contrast:.17g}",
            "--fake.subtractPlanet=false",
            "--output.directory",
            str(evaluation_directory),
            "--output.fileName",
            "finim.fits",
            "--output.exactFName=true",
            "--showTiming=true",
        ]
        (evaluation_directory / "command.txt").write_text(shlex.join(command) + "\n", encoding="utf-8")
        begin = time.monotonic()
        with (evaluation_directory / "run.log").open("w", encoding="utf-8") as log_stream:
            completed = subprocess.run(command, stdout=log_stream, stderr=subprocess.STDOUT, check=False)
        elapsed = time.monotonic() - begin
        if completed.returncode != 0 or not final_path.is_file():
            raise RuntimeError(
                f"KLIP optimizer evaluation {evaluation_index} failed with status {completed.returncode}; "
                f"see {evaluation_directory / 'run.log'}"
            )

        image, header = read_final_image(final_path, arguments.mode_count)
        if expected_shape is None:
            expected_shape = image.shape
            detector_columns, detector_rows = image.shape
            center_row = 0.5 * (detector_rows - 1)
            center_column = 0.5 * (detector_columns - 1)
            aperture_coordinates = (
                center_row + initial_row_offset,
                center_column + initial_column_offset,
            )
        if image.shape != expected_shape or aperture_coordinates is None:
            raise RuntimeError("KLIP optimizer output dimensions changed between evaluations")
        merit, aperture_samples = aperture_merit(
            image,
            aperture_coordinates[0],
            aperture_coordinates[1],
            arguments.aperture_radius,
        )
        fake_separation = header_vector(header, "FAKESEP")
        fake_pa = header_vector(header, "FAKEPA")
        fake_contrast = header_vector(header, "FAKECONT")
        if (
            len(fake_separation) != 1
            or len(fake_pa) != 1
            or len(fake_contrast) != 1
            or not math.isclose(fake_separation[0], separation, rel_tol=0, abs_tol=5e-4)
            or abs((fake_pa[0] - position_angle + 180) % 360 - 180) > 5e-4
            or not math.isclose(fake_contrast[0], -contrast, rel_tol=0, abs_tol=5e-9)
        ):
            raise RuntimeError("KLIP optimizer output does not record the requested negative fake")

        evaluation = {
            "complete": True,
            "evaluation": evaluation_index,
            "normalized_point": [float(value) for value in bounded],
            "row_delta": row_delta,
            "column_delta": column_delta,
            "separation": separation,
            "position_angle": position_angle,
            "contrast": contrast,
            "injected_contrast": -contrast,
            "merit": merit,
            "aperture_samples": aperture_samples,
            "elapsed_seconds": elapsed,
            "final_image": str(final_path),
        }
        (evaluation_directory / "evaluation.json").write_text(
            json.dumps(evaluation, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        cache[key] = evaluation
        evaluations.append(evaluation)
        if merit_scale is None:
            merit_scale = merit
        print(
            f"evaluation {evaluation_index}: rowDelta={row_delta:.8g}, columnDelta={column_delta:.8g}, "
            f"sep={separation:.10g}, PA={position_angle:.10g}, contrast={contrast:.10g}, merit={merit:.10g}",
            flush=True,
        )
        return merit / merit_scale

    normalized_contrast = 2 * (arguments.initial_contrast - arguments.contrast_lower) / contrast_range - 1
    initial_point = np.asarray([0.0, 0.0, normalized_contrast], dtype=np.float64)
    evaluate(initial_point)
    result = minimize(
        evaluate,
        initial_point,
        method="Powell",
        bounds=[(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)],
        options={
            "maxfev": arguments.maximum_evaluations,
            "xtol": arguments.parameter_tolerance,
            "ftol": arguments.merit_tolerance,
            "disp": True,
        },
    )

    if not evaluations:
        raise RuntimeError("KLIP optimizer completed without any evaluated reductions")
    best = min(evaluations, key=lambda evaluation: float(evaluation["merit"]))
    best_point = np.asarray(best["normalized_point"], dtype=np.float64)
    bound_coordinates = [index for index, value in enumerate(best_point) if abs(abs(value) - 1) <= 1e-3]
    status = "converged"
    if not bool(result.success):
        status = "optimizer-not-converged"
    elif bound_coordinates:
        status = "best-fit-at-bound"

    evaluations.sort(key=lambda evaluation: int(evaluation["evaluation"]))
    table_path = output_directory / "evaluations.csv"
    with table_path.open("w", encoding="utf-8", newline="") as stream:
        fieldnames = [
            "evaluation",
            "row_delta",
            "column_delta",
            "separation",
            "position_angle",
            "contrast",
            "injected_contrast",
            "merit",
            "aperture_samples",
            "elapsed_seconds",
            "final_image",
        ]
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(evaluations)

    config_path = output_directory / "best_signal_free.conf"
    config_path.write_text(
        "[planet]\n"
        f"sep={float(best['separation']):.17g}\n"
        f"PA={float(best['position_angle']):.17g}\n"
        f"contrast={float(best['contrast']):.17g}\n"
        "\n[fake]\n"
        "method=single\n"
        f"fileName={psf_file}\n"
        "subtractPlanet=true\n",
        encoding="utf-8",
    )
    summary = {
        "schema": 1,
        "complete": True,
        "status": status,
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "optimizer_evaluations": int(result.nfev),
        "stored_evaluations": len(evaluations),
        "aperture_frame": "fixed-initial-sky",
        "aperture_row": aperture_coordinates[0] if aperture_coordinates is not None else None,
        "aperture_column": aperture_coordinates[1] if aperture_coordinates is not None else None,
        "best": best,
        "best_signal_free_config": str(config_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        "KLIP negative-planet fit: "
        f"status={status}, separation={float(best['separation']):.12g}, "
        f"PA={float(best['position_angle']):.12g}, contrast={float(best['contrast']):.12g}, "
        f"merit={float(best['merit']):.12g}, evaluations={len(evaluations)}"
    )
    if status != "converged":
        raise RuntimeError(f"KLIP negative-planet fit did not establish an interior converged optimum: {status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
