#!/usr/bin/env python3
"""Diagnose registration and stamp support in a KLIP planet-response estimate."""

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
from matplotlib.patches import Rectangle
from scipy.ndimage import shift
from scipy.optimize import differential_evolution

from compare_klip_central_response import fraction_tag, shape_metrics
from compare_klip_finite_response import evaluate_response, mode_counts, read_cube, response_paths


COLORMAP = "RdBu_r"


def crop_response(image: np.ndarray, row: int, column: int, size: int) -> np.ndarray:
    """Extract an odd square stamp from an Eigen-oriented response image."""
    if size <= 0 or size % 2 == 0:
        raise RuntimeError("response stamp size must be positive and odd")
    half_width = size // 2
    result = image[
        row - half_width : row + half_width + 1,
        column - half_width : column + half_width + 1,
    ]
    if result.shape != (size, size):
        raise RuntimeError("response stamp extends outside the paired final image")
    return result


def registered_metrics(
    planet: np.ndarray,
    estimate: np.ndarray,
    displacement: np.ndarray,
) -> tuple[np.ndarray, dict[str, float]]:
    """Shift an estimate and return fixed-support response metrics."""
    registered = shift(
        estimate,
        displacement,
        order=3,
        mode="nearest",
        prefilter=True,
    )
    validity = np.isfinite(planet) & np.isfinite(registered)
    return registered, shape_metrics(planet, registered, validity)


def optimize_registration(
    planet: np.ndarray,
    estimate: np.ndarray,
    maximum_shift: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Find the bounded subpixel displacement that maximizes response cosine."""
    def objective(displacement: np.ndarray) -> float:
        """Return negative cosine for the registration optimizer."""
        _, metrics = registered_metrics(planet, estimate, displacement)
        return -float(metrics["cosine"])

    result = differential_evolution(
        objective,
        [(-maximum_shift, maximum_shift), (-maximum_shift, maximum_shift)],
        seed=7,
        tol=1e-10,
        polish=True,
    )
    if not result.success or not np.isfinite(result.x).all():
        raise RuntimeError(f"response registration did not converge: {result.message}")
    registered, metrics = registered_metrics(planet, estimate, result.x)
    return np.asarray(result.x), registered, metrics


def scaled_residual(
    planet: np.ndarray,
    estimate: np.ndarray,
    metrics: dict[str, float],
) -> np.ndarray:
    """Subtract the best amplitude-scaled estimate from the planet response."""
    return planet - float(metrics["projection"]) * estimate


def draw_response(
    axis: plt.Axes,
    response: np.ndarray,
    title: str,
    normalization: Normalize,
) -> matplotlib.image.AxesImage:
    """Draw one Eigen-oriented response stamp in detector coordinates."""
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
    """Write a registration- and support-aware planet-response comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native_experiment", type=Path, help="completed native KLIP response experiment")
    parser.add_argument("planet_experiment", type=Path, help="paired-refit experiment containing a candidate sample")
    parser.add_argument("output", type=Path, help="output PNG path")
    parser.add_argument("--case", default="radial_ld_refit4_filter", help="native response case directory")
    parser.add_argument("--mode-count", type=int, default=200)
    parser.add_argument("--amplitude-fraction", type=float, default=1.0)
    parser.add_argument("--context-size", type=int, default=31)
    parser.add_argument("--maximum-shift", type=float, default=0.75)
    parser.add_argument("--overwrite", action="store_true")
    arguments = parser.parse_args()

    if arguments.mode_count <= 0:
        raise RuntimeError("--mode-count must be positive")
    if not math.isfinite(arguments.amplitude_fraction) or arguments.amplitude_fraction <= 0:
        raise RuntimeError("--amplitude-fraction must be positive and finite")
    if arguments.context_size <= 0 or arguments.context_size % 2 == 0:
        raise RuntimeError("--context-size must be positive and odd")
    if not math.isfinite(arguments.maximum_shift) or arguments.maximum_shift <= 0:
        raise RuntimeError("--maximum-shift must be positive and finite")
    output = arguments.output.resolve()
    if output.exists() and not arguments.overwrite:
        raise RuntimeError(f"output already exists: {output}")

    case_directory = arguments.native_experiment.resolve() / arguments.case
    paths = response_paths(case_directory)
    if arguments.mode_count not in paths:
        raise RuntimeError(f"native response does not contain KL mode {arguments.mode_count}")

    planet_directory = arguments.planet_experiment.resolve()
    summary_path = planet_directory / "klip_central_response.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    matching_rows = [
        row
        for row in summary["rows"]
        if row["sample_role"] == "candidate"
        and int(row["mode_count"]) == arguments.mode_count
        and math.isclose(
            float(row["amplitude_fraction"]),
            arguments.amplitude_fraction,
            rel_tol=0,
            abs_tol=1e-15,
        )
    ]
    if len(matching_rows) != 1:
        raise RuntimeError("paired experiment does not contain a unique requested planet response")
    row = matching_rows[0]

    run_directory = (
        planet_directory
        / "runs"
        / str(row["sample_label"])
        / f"fraction_{fraction_tag(arguments.amplitude_fraction)}"
    )
    plus_path = run_directory / "plus" / "finim.fits"
    minus_path = run_directory / "minus" / "finim.fits"
    plus, plus_header = read_cube(plus_path)
    minus, minus_header = read_cube(minus_path)
    plus_modes = mode_counts(plus_header, plus_path)
    minus_modes = mode_counts(minus_header, minus_path)
    if plus.shape != minus.shape or plus_modes != minus_modes or arguments.mode_count not in plus_modes:
        raise RuntimeError("paired positive and negative planet reductions are inconsistent")
    mode_index = plus_modes.index(arguments.mode_count)
    full_response = ((plus[mode_index] - minus[mode_index]) / (2 * float(row["epsilon"]))).T

    target_row = int(row["target_row"])
    target_column = int(row["target_column"])
    estimate, estimate_validity = evaluate_response(
        paths[arguments.mode_count],
        float(row["target_radius"]),
        float(row["target_angle_radians"]),
    )
    stamp_size = estimate.shape[0]
    if estimate.shape[0] != estimate.shape[1] or stamp_size % 2 == 0:
        raise RuntimeError("native response estimate must be square with odd dimensions")
    planet = crop_response(full_response, target_row, target_column, stamp_size)
    context = crop_response(full_response, target_row, target_column, arguments.context_size)
    retained = np.isfinite(planet) & estimate_validity & np.isfinite(estimate)
    raw_metrics = shape_metrics(planet, estimate, retained)
    displacement, registered, registered_values = optimize_registration(
        planet,
        np.where(estimate_validity, estimate, 0),
        arguments.maximum_shift,
    )
    raw_residual = scaled_residual(planet, estimate, raw_metrics)
    registered_residual = scaled_residual(planet, registered, registered_values)

    stamp_energy = float(np.sum(planet[retained] ** 2))
    context_energy = float(np.sum(context[np.isfinite(context)] ** 2))
    if stamp_energy <= 0 or context_energy <= 0:
        raise RuntimeError("planet response has no positive finite energy")
    retained_energy_fraction = stamp_energy / context_energy
    border = np.zeros(planet.shape, dtype=bool)
    border[[0, -1], :] = True
    border[:, [0, -1]] = True
    border_energy_fraction = float(np.sum(planet[border & retained] ** 2)) / stamp_energy

    display_limit = float(
        np.max(
            np.abs(
                np.concatenate(
                    (
                        planet[retained],
                        estimate[retained],
                        registered[np.isfinite(registered)],
                    )
                )
            )
        )
    )
    normalization = Normalize(vmin=-display_limit, vmax=display_limit)
    figure, axes = plt.subplots(2, 3, figsize=(11.4, 7.7), constrained_layout=True)
    images = (
        context,
        planet,
        estimate,
        raw_residual,
        registered,
        registered_residual,
    )
    titles = (
        f"{arguments.context_size}×{arguments.context_size} paired planet response",
        f"Planet response: {stamp_size}×{stamp_size} crop",
        f"Raw sparse estimate\ncos={raw_metrics['cosine']:.4f}",
        f"Raw best-scaled residual\nrelative norm={raw_metrics['best_scaled_relative_residual']:.3f}",
        f"Registered estimate\nshift=({displacement[0]:+.3f}, {displacement[1]:+.3f}) pix",
        f"Registered best-scaled residual\ncos={registered_values['cosine']:.4f}, "
        f"relative norm={registered_values['best_scaled_relative_residual']:.3f}",
    )
    for axis, image, title in zip(axes.flat, images, titles):
        artist = draw_response(axis, image, title, normalization)
    half_stamp = 0.5 * stamp_size
    axes[0, 0].add_patch(
        Rectangle(
            (-half_stamp, -half_stamp),
            stamp_size,
            stamp_size,
            fill=False,
            edgecolor="#f6c85f",
            linewidth=1.6,
        )
    )
    colorbar = figure.colorbar(artist, ax=axes, shrink=0.88, pad=0.015)
    colorbar.set_label("final-image response per unit input contrast")
    figure.suptitle(
        "KLIP response registration and stamp-support diagnostic\n"
        f"KL modes={arguments.mode_count}; candidate sep={float(row['configured_separation']):.3f} pix, "
        f"PA={float(row['configured_position_angle']):.3f}°; "
        f"{stamp_size}×{stamp_size} contains {retained_energy_fraction:.1%} of "
        f"{arguments.context_size}×{arguments.context_size} energy "
        f"({border_energy_fraction:.1%} on outer ring)",
        fontsize=12.5,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=190)
    plt.close(figure)
    print(f"Wrote {output}")
    print(
        f"raw cosine={raw_metrics['cosine']:.8g}, registered cosine={registered_values['cosine']:.8g}, "
        f"shift=({displacement[0]:.8g},{displacement[1]:.8g}), "
        f"stamp/context energy={retained_energy_fraction:.8g}, border/stamp energy={border_energy_fraction:.8g}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
