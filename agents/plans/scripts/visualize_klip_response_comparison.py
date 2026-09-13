#!/usr/bin/env python3
"""Compare a native sparse KLIP response estimate with the measured planet response."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.colors import Normalize

from compare_klip_central_response import shape_metrics
from compare_klip_finite_response import evaluate_response, response_paths


COLORMAP = "RdBu_r"


def parse_modes(text: str) -> list[int]:
    """Parse a comma-separated list of unique positive KL mode counts."""
    try:
        modes = [int(token) for token in text.split(",")]
    except ValueError as error:
        raise argparse.ArgumentTypeError("modes must be comma-separated integers") from error
    if not modes or any(mode <= 0 for mode in modes) or len(set(modes)) != len(modes):
        raise argparse.ArgumentTypeError("modes must be unique positive integers")
    return modes


def finite_limit(arrays: list[np.ndarray]) -> float:
    """Return a positive symmetric display limit from finite array values."""
    finite = np.concatenate([np.abs(array[np.isfinite(array)]) for array in arrays])
    if finite.size == 0:
        raise RuntimeError("response comparison contains no finite values")
    limit = float(np.max(finite))
    if not math.isfinite(limit) or limit <= 0:
        raise RuntimeError("response comparison has no positive finite display scale")
    return limit


def draw_response(
    axis: plt.Axes,
    response: np.ndarray,
    title: str,
    normalization: Normalize,
) -> matplotlib.image.AxesImage:
    """Draw one Eigen-oriented response stamp in Cartesian detector coordinates."""
    half_rows = 0.5 * response.shape[0]
    half_columns = 0.5 * response.shape[1]
    artist = axis.imshow(
        response,
        origin="lower",
        interpolation="nearest",
        extent=(-half_columns, half_columns, -half_rows, half_rows),
        cmap=COLORMAP,
        norm=normalization,
    )
    axis.axhline(0, color="0.2", linewidth=0.35, alpha=0.35)
    axis.axvline(0, color="0.2", linewidth=0.35, alpha=0.35)
    axis.set_title(title, fontsize=10)
    axis.set_xlabel(r"$\Delta$ column [pix]")
    axis.set_ylabel(r"$\Delta$ row [pix]")
    axis.set_aspect("equal")
    return artist


def main() -> int:
    """Write the requested KLIP planet-response comparison PNG."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native_experiment", type=Path, help="completed native KLIP response experiment")
    parser.add_argument("planet_experiment", type=Path, help="paired-refit experiment containing a candidate sample")
    parser.add_argument("output", type=Path, help="output PNG path")
    parser.add_argument("--case", default="radial_ld_refit4_filter", help="native response case directory")
    parser.add_argument("--modes", type=parse_modes, default=parse_modes("125,200,350"))
    parser.add_argument("--amplitude-fraction", type=float, default=1.0)
    parser.add_argument("--overwrite", action="store_true")
    arguments = parser.parse_args()

    if not math.isfinite(arguments.amplitude_fraction) or arguments.amplitude_fraction <= 0:
        raise RuntimeError("--amplitude-fraction must be positive and finite")
    output = arguments.output.resolve()
    if output.exists() and not arguments.overwrite:
        raise RuntimeError(f"output already exists: {output}")

    case_directory = arguments.native_experiment.resolve() / arguments.case
    model_paths = response_paths(case_directory)
    missing_modes = sorted(set(arguments.modes) - set(model_paths))
    if missing_modes:
        raise RuntimeError(f"native response products do not contain KL modes {missing_modes}")

    planet_directory = arguments.planet_experiment.resolve()
    summary_path = planet_directory / "klip_central_response.json"
    cube_path = planet_directory / "central_response.fits"
    if not summary_path.is_file() or not cube_path.is_file():
        raise RuntimeError("planet experiment lacks its central-response summary or cube")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    planet_cube = np.asarray(fits.getdata(cube_path), dtype=np.float64)
    if planet_cube.ndim != 3:
        raise RuntimeError("planet central-response product is not a response cube")

    selected_rows: dict[int, dict[str, object]] = {}
    for row in summary["rows"]:
        mode = int(row["mode_count"])
        if (
            row["sample_role"] == "candidate"
            and mode in arguments.modes
            and math.isclose(
                float(row["amplitude_fraction"]),
                arguments.amplitude_fraction,
                rel_tol=0,
                abs_tol=1e-15,
            )
        ):
            if mode in selected_rows:
                raise RuntimeError(f"planet response is not unique for KL mode {mode}")
            selected_rows[mode] = row
    if set(selected_rows) != set(arguments.modes):
        raise RuntimeError("planet experiment does not contain every requested candidate response")

    comparisons: list[dict[str, object]] = []
    for mode in arguments.modes:
        row = selected_rows[mode]
        row_index = int(row["row_index"])
        if row_index < 0 or row_index >= planet_cube.shape[0]:
            raise RuntimeError(f"planet response row is outside the response cube for KL mode {mode}")
        planet = planet_cube[row_index].T
        estimate, estimate_validity = evaluate_response(
            model_paths[mode],
            float(row["target_radius"]),
            float(row["target_angle_radians"]),
        )
        retained = np.isfinite(planet) & estimate_validity & np.isfinite(estimate)
        metrics = shape_metrics(planet, estimate, retained)
        scale = float(metrics["projection"])
        scaled_estimate = scale * estimate
        comparisons.append(
            {
                "mode": mode,
                "planet": np.where(retained, planet, np.nan),
                "estimate": np.where(retained, estimate, np.nan),
                "residual": np.where(retained, planet - scaled_estimate, np.nan),
                "metrics": metrics,
                "row": row,
            }
        )

    figure, axes = plt.subplots(
        len(comparisons),
        3,
        figsize=(10.9, 3.15 * len(comparisons)),
        constrained_layout=True,
        squeeze=False,
    )
    for row_index, comparison in enumerate(comparisons):
        planet = comparison["planet"]
        estimate = comparison["estimate"]
        residual = comparison["residual"]
        metrics = comparison["metrics"]
        limit = finite_limit([planet, estimate])
        normalization = Normalize(vmin=-limit, vmax=limit)
        titles = (
            f"KL {comparison['mode']}: paired-refit planet response",
            "Sparse radial estimate",
            f"Planet − {metrics['projection']:.3f} × estimate",
        )
        for column, (image, title) in enumerate(zip((planet, estimate, residual), titles)):
            artist = draw_response(axes[row_index, column], image, title, normalization)
        colorbar = figure.colorbar(artist, ax=axes[row_index], shrink=0.88, pad=0.01)
        colorbar.set_label("final-image response per unit input contrast")
        axes[row_index, 2].text(
            0.03,
            0.04,
            f"cosine = {metrics['cosine']:.4f}\n"
            f"best-scaled residual = {metrics['best_scaled_relative_residual']:.3f}",
            transform=axes[row_index, 2].transAxes,
            fontsize=8.5,
            color="black",
            bbox={"facecolor": "white", "edgecolor": "0.5", "alpha": 0.82, "pad": 3},
        )

    first_row = comparisons[0]["row"]
    figure.suptitle(
        "Native sparse KLIP response estimate versus independent planet response\n"
        f"candidate sep={float(first_row['configured_separation']):.3f} pix, "
        f"PA={float(first_row['configured_position_angle']):.3f}°, "
        f"paired perturbation={arguments.amplitude_fraction:g} × fiducial",
        fontsize=13,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=190)
    plt.close(figure)
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
