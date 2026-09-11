#!/usr/bin/env python3
"""Visualize the final-image cube and PSF model in a P4 response experiment."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.colors import Normalize


MODE_KEYWORDS = ("P4 MODE FRACTIONS", "P4MODFR", "FRACT NMODES", "NMODES")
COLORMAP = "RdBu_r"


def header_vector(header: fits.Header, keywords: tuple[str, ...]) -> np.ndarray:
    """Read the first available comma-separated numeric FITS vector."""
    for keyword in keywords:
        if keyword in header:
            return np.asarray([float(value) for value in str(header[keyword]).split(",")])
    raise RuntimeError(f"none of the required FITS keywords are present: {', '.join(keywords)}")


def selected_mode_index(header: fits.Header, requested_mode: float, source: Path) -> int:
    """Resolve one requested mode without nearest-mode substitution."""
    modes = header_vector(header, MODE_KEYWORDS)
    tolerance = 8 * np.finfo(np.float32).eps * np.maximum.reduce(
        (np.ones_like(modes), np.abs(modes), np.full_like(modes, abs(requested_mode)))
    )
    matches = np.flatnonzero(np.abs(modes - requested_mode) <= tolerance)
    if matches.size != 1:
        raise RuntimeError(f"mode {requested_mode:.17g} is not an exact unique mode in {source}")
    return int(matches[0])


def read_exact_planet(path: Path) -> dict[str, float]:
    """Read the converged exact negative-planet point estimate."""
    text = path.read_text(encoding="utf-8")
    fitted = re.search(r"^  fitted:\n(?P<body>(?:    [^\n]*\n)+)", text, flags=re.MULTILINE)
    if fitted is None or not re.search(r"^  converged: true$", text, flags=re.MULTILINE):
        raise RuntimeError(f"exact optimizer did not converge: {path}")

    result: dict[str, float] = {}
    for key in ("separation", "positionAngle", "contrast"):
        match = re.search(rf"^    {key}: ([+\-0-9.eE]+)$", fitted.group("body"), flags=re.MULTILINE)
        if match is None:
            raise RuntimeError(f"exact optimizer summary lacks fitted {key}: {path}")
        result[key] = float(match.group(1))
    if result["contrast"] >= 0:
        raise RuntimeError(f"exact optimizer contrast is not a negative injection: {path}")
    return {
        "separation": result["separation"],
        "position_angle": result["positionAngle"] % 360,
        "contrast": -result["contrast"],
    }


def read_cube(path: Path) -> tuple[np.ndarray, fits.Header]:
    """Read a two- or three-dimensional final image as a three-dimensional cube."""
    header = fits.getheader(path)
    data = np.asarray(fits.getdata(path), dtype=np.float64)
    if data.ndim == 2:
        data = data[np.newaxis, :, :]
    if data.ndim != 3:
        raise RuntimeError(f"expected a two- or three-dimensional FITS image: {path}")
    return data, header


def read_coordinates(path: Path) -> np.ndarray:
    """Read response coordinates as row, column, region, and regional-index records."""
    coordinates = np.asarray(fits.getdata(path), dtype=np.float64)
    if coordinates.ndim != 2:
        raise RuntimeError(f"response coordinates are not two-dimensional: {path}")
    if coordinates.shape[0] == 4:
        coordinates = coordinates.T
    if coordinates.shape[1] != 4:
        raise RuntimeError(f"response coordinates do not have four fields: {path}")
    return coordinates


def source_index(coordinates: np.ndarray, row: int, column: int, path: Path) -> int:
    """Return the unique response-model index at one detector coordinate."""
    matches = np.flatnonzero((coordinates[:, 0] == row) & (coordinates[:, 1] == column))
    if matches.size != 1:
        raise RuntimeError(f"{path} does not uniquely contain response coordinate ({row}, {column})")
    return int(matches[0])


def crop_at(
    image: np.ndarray,
    center_row: int,
    center_column: int,
    size: int,
) -> np.ndarray:
    """Extract a square detector stamp from an image stored as column by row."""
    if size <= 0 or size % 2 == 0:
        raise RuntimeError("display stamp size must be positive and odd")
    half_width = size // 2
    stamp = image[
        center_column - half_width : center_column + half_width + 1,
        center_row - half_width : center_row + half_width + 1,
    ]
    if stamp.shape != (size, size):
        raise RuntimeError("display stamp extends beyond the final image")
    return stamp


def local_crop(image: np.ndarray, header: fits.Header, shape: tuple[int, int]) -> np.ndarray:
    """Extract the detector region described by a local P4 product header."""
    first_row = int(header["P4 LOCAL ORIGIN ROW"])
    first_column = int(header["P4 LOCAL ORIGIN COLUMN"])
    columns, rows = shape
    stamp = image[first_column : first_column + columns, first_row : first_row + rows]
    if stamp.shape != shape:
        raise RuntimeError("local final-image product extends beyond the complete final image")
    return stamp


def response_metrics(empirical: np.ndarray, model: np.ndarray) -> dict[str, float]:
    """Measure shape agreement and the projection of a finite response onto a model."""
    retained = np.isfinite(empirical) & np.isfinite(model)
    empirical_values = empirical[retained]
    model_values = model[retained]
    empirical_energy = float(np.dot(empirical_values, empirical_values))
    model_energy = float(np.dot(model_values, model_values))
    cross_energy = float(np.dot(empirical_values, model_values))
    if empirical_energy <= 0 or model_energy <= 0:
        raise RuntimeError("finite and analytic responses must have positive energy")
    return {
        "cosine_similarity": cross_energy / math.sqrt(empirical_energy * model_energy),
        "projection_ratio_to_exact_contrast": cross_energy / model_energy,
        "relative_l2": math.sqrt(
            float(np.dot(empirical_values - model_values, empirical_values - model_values))
            / empirical_energy
        ),
    }


def finite_limit(arrays: list[np.ndarray], percentile: float = 100.0) -> float:
    """Return a positive symmetric color limit from finite values."""
    values = np.concatenate([np.abs(array[np.isfinite(array)]) for array in arrays])
    if values.size == 0:
        raise RuntimeError("cannot scale a collection without finite image values")
    limit = float(np.percentile(values, percentile))
    if not math.isfinite(limit) or limit <= 0:
        raise RuntimeError("image collection does not have a positive finite color limit")
    return limit


def image_extent(shape: tuple[int, int]) -> tuple[float, float, float, float]:
    """Return a center-relative image extent for a detector stamp."""
    columns, rows = shape
    return (-0.5 * columns, 0.5 * columns, -0.5 * rows, 0.5 * rows)


def draw_stamp(
    axis: plt.Axes,
    stamp: np.ndarray,
    title: str,
    normalization: Normalize,
    source_offset: tuple[float, float] | None = None,
) -> matplotlib.image.AxesImage:
    """Draw one column-by-row detector stamp in conventional Cartesian orientation."""
    image = axis.imshow(
        stamp.T,
        origin="lower",
        interpolation="nearest",
        extent=image_extent(stamp.shape),
        cmap=COLORMAP,
        norm=normalization,
    )
    if source_offset is not None:
        axis.plot(
            source_offset[0],
            source_offset[1],
            marker="+",
            markersize=8,
            markeredgewidth=1.2,
            color="#f6c85f",
        )
    axis.set_title(title, fontsize=10)
    axis.set_xlabel(r"$\Delta$ column [pix]")
    axis.set_ylabel(r"$\Delta$ row [pix]")
    axis.set_aspect("equal")
    return image


def write_overview(
    path: Path,
    original: np.ndarray,
    exact_best: np.ndarray,
    full_best: np.ndarray,
    empirical: np.ndarray,
    sparse_model: np.ndarray,
    dense_model: np.ndarray,
    source_offset: tuple[float, float],
    mode_fraction: float,
    exact_contrast: float,
    sparse_metrics: dict[str, float],
    dense_metrics: dict[str, float],
    full_local_fraction: float,
) -> None:
    """Write the selected-mode final-image and response diagnostic."""
    science_limit = finite_limit([original, exact_best, full_best], percentile=99.5)
    response_limit = finite_limit([empirical, sparse_model, dense_model])
    science_norm = Normalize(vmin=-science_limit, vmax=science_limit)
    response_norm = Normalize(vmin=-response_limit, vmax=response_limit)

    figure, axes = plt.subplots(2, 3, figsize=(12.4, 8.0), constrained_layout=True)
    science_titles = (
        "Original final image",
        "Exact optimizer: best negative planet",
        "Full rerun: same negative planet",
    )
    science_images = (original, exact_best, full_best)
    for axis, image, title in zip(axes[0], science_images, science_titles):
        artist = draw_stamp(axis, image, title, science_norm, source_offset)
    science_bar = figure.colorbar(artist, ax=axes[0], shrink=0.84, pad=0.01)
    science_bar.set_label("final-image value")

    response_titles = (
        "End-to-end removed signal / exact contrast",
        "Sparse analytic response used by fit",
        "Dense signal-free analytic response",
    )
    response_images = (empirical, sparse_model, dense_model)
    for axis, image, title in zip(axes[1], response_images, response_titles):
        artist = draw_stamp(axis, image, title, response_norm, source_offset)
    response_bar = figure.colorbar(artist, ax=axes[1], shrink=0.84, pad=0.01)
    response_bar.set_label("response per unit input contrast")

    figure.suptitle(
        "P4 final images and fitted PSF response\n"
        f"mode fraction {mode_fraction:.3g}; exact contrast {exact_contrast:.7g}; "
        f"best-local/full difference {full_local_fraction:.3%} of removed signal\n"
        "finite vs sparse: "
        f"cos={sparse_metrics['cosine_similarity']:.3f}, "
        f"projection/exact={sparse_metrics['projection_ratio_to_exact_contrast']:.3f}; "
        "finite vs dense: "
        f"cos={dense_metrics['cosine_similarity']:.3f}, "
        f"projection/exact={dense_metrics['projection_ratio_to_exact_contrast']:.3f}",
        fontsize=12,
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)


def write_cube_animation(
    path: Path,
    modes: np.ndarray,
    original: np.ndarray,
    exact_best: np.ndarray,
    models: list[np.ndarray],
    exact_contrast: float,
    source_offset: tuple[float, float],
) -> None:
    """Write an animation that walks the final-image cube and its PSF estimate by mode."""
    figure, axes = plt.subplots(1, 4, figsize=(12.2, 3.55), constrained_layout=True)

    def draw_frame(index: int) -> None:
        """Render one output-mode frame into the animation figure."""
        for axis in axes:
            axis.clear()
        empirical = (original[index] - exact_best[index]) / exact_contrast
        science_limit = finite_limit([original[index], exact_best[index]], percentile=99.5)
        response_limit = finite_limit([empirical, models[index]])
        science_norm = Normalize(vmin=-science_limit, vmax=science_limit)
        response_norm = Normalize(vmin=-response_limit, vmax=response_limit)
        draw_stamp(axes[0], original[index], "Original final", science_norm, source_offset)
        draw_stamp(axes[1], exact_best[index], "Best negative subtracted", science_norm, source_offset)
        draw_stamp(axes[2], empirical, "Removed signal / contrast", response_norm, source_offset)
        draw_stamp(axes[3], models[index], "Sparse analytic PSF", response_norm, source_offset)
        figure.suptitle(
            f"P4 mode fraction {modes[index]:.3g}   "
            f"final-image scale ±{science_limit:.3g}   response scale ±{response_limit:.3g}",
            fontsize=12,
        )

    movie = animation.FuncAnimation(figure, draw_frame, frames=len(modes), interval=900, repeat=True)
    movie.save(path, writer=animation.PillowWriter(fps=1.15), dpi=105)
    plt.close(figure)


def main() -> int:
    """Write a selected-mode diagnostic, cube animation, and machine-readable summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", type=Path, help="completed p4_matched_response experiment")
    parser.add_argument(
        "output_directory",
        type=Path,
        nargs="?",
        help="output directory (default: EXPERIMENT/visualization)",
    )
    parser.add_argument("--mode-fraction", type=float, default=0.15)
    parser.add_argument("--skip-animation", action="store_true")
    arguments = parser.parse_args()

    experiment = arguments.experiment.resolve()
    output_directory = (
        arguments.output_directory.resolve()
        if arguments.output_directory is not None
        else experiment / "visualization"
    )
    original_path = experiment / "sparse_response" / "finim.fits"
    exact_best_path = experiment / "exact_optimizer" / "finim_outputs" / "p4Negative_best.fits"
    exact_summary_path = experiment / "exact_optimizer" / "finim_outputs" / "p4Negative_summary.yaml"
    full_best_path = experiment / "signal_free_oracle" / "finim.fits"
    sparse_directory = experiment / "sparse_response" / "finim_outputs"
    dense_directory = experiment / "signal_free_oracle" / "finim_outputs"
    sparse_fit_path = experiment / "sparse_fit" / "summary.json"
    required = (
        original_path,
        exact_best_path,
        exact_summary_path,
        full_best_path,
        sparse_directory / "p4PSF_coordinates.fits",
        dense_directory / "p4PSF_coordinates.fits",
        sparse_fit_path,
    )
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise RuntimeError("missing required matched-response products: " + ", ".join(map(str, missing)))

    original_cube, original_header = read_cube(original_path)
    exact_best_cube, exact_best_header = read_cube(exact_best_path)
    full_best_cube, full_best_header = read_cube(full_best_path)
    modes = header_vector(original_header, MODE_KEYWORDS)
    if modes.size != original_cube.shape[0] or exact_best_cube.shape[0] != original_cube.shape[0]:
        raise RuntimeError("original and exact-best final-image cubes do not have the same mode axis")
    mode_index = selected_mode_index(original_header, arguments.mode_fraction, original_path)
    full_mode_index = selected_mode_index(full_best_header, arguments.mode_fraction, full_best_path)

    exact = read_exact_planet(exact_summary_path)
    detector_columns, detector_rows = original_cube.shape[1:]
    center_row = 0.5 * (detector_rows - 1)
    center_column = 0.5 * (detector_columns - 1)
    angle = math.radians(exact["position_angle"])
    exact_row = center_row - exact["separation"] * math.sin(angle)
    exact_column = center_column + exact["separation"] * math.cos(angle)

    sparse_fit = json.loads(sparse_fit_path.read_text(encoding="utf-8"))
    target_row = int(sparse_fit["fit"]["integer_peak_row"])
    target_column = int(sparse_fit["fit"]["integer_peak_column"])
    source_offset = (exact_column - target_column, exact_row - target_row)

    sparse_coordinates_path = sparse_directory / "p4PSF_coordinates.fits"
    sparse_coordinates = read_coordinates(sparse_coordinates_path)
    sparse_source = source_index(sparse_coordinates, target_row, target_column, sparse_coordinates_path)
    sparse_models: list[np.ndarray] = []
    for index in range(len(modes)):
        model_path = sparse_directory / f"p4PSF_model_{index:04d}.fits"
        validity_path = sparse_directory / f"p4PSF_validity_{index:04d}.fits"
        if not model_path.is_file() or not validity_path.is_file():
            raise RuntimeError(f"missing sparse response mode {index} in {sparse_directory}")
        validity = np.asarray(fits.getdata(validity_path)).reshape(-1)
        if not math.isfinite(float(validity[sparse_source])) or validity[sparse_source] <= 0:
            raise RuntimeError(f"sparse response is invalid at ({target_row}, {target_column}), mode {index}")
        sparse_models.append(np.asarray(fits.getdata(model_path)[sparse_source], dtype=np.float64))

    dense_coordinates_path = dense_directory / "p4PSF_coordinates.fits"
    dense_coordinates = read_coordinates(dense_coordinates_path)
    dense_source = source_index(dense_coordinates, target_row, target_column, dense_coordinates_path)
    dense_model_path = dense_directory / f"p4PSF_model_{full_mode_index:04d}.fits"
    dense_validity_path = dense_directory / f"p4PSF_validity_{full_mode_index:04d}.fits"
    dense_validity = np.asarray(fits.getdata(dense_validity_path)).reshape(-1)
    if not math.isfinite(float(dense_validity[dense_source])) or dense_validity[dense_source] <= 0:
        raise RuntimeError(f"dense response is invalid at ({target_row}, {target_column})")
    dense_model = np.asarray(fits.getdata(dense_model_path)[dense_source], dtype=np.float64)

    local_columns, local_rows = exact_best_cube.shape[1:]
    original_local_cube = np.stack(
        [local_crop(image, exact_best_header, (local_columns, local_rows)) for image in original_cube]
    )
    full_best_local = local_crop(
        full_best_cube[full_mode_index], exact_best_header, (local_columns, local_rows)
    )
    response_size = sparse_models[mode_index].shape[0]
    if sparse_models[mode_index].shape != dense_model.shape or response_size % 2 == 0:
        raise RuntimeError("sparse and dense selected-mode responses must be matching odd square stamps")
    original_response_stamp = crop_at(
        original_cube[mode_index], target_row, target_column, response_size
    )
    full_best_response_stamp = crop_at(
        full_best_cube[full_mode_index], target_row, target_column, response_size
    )
    empirical_response = (original_response_stamp - full_best_response_stamp) / exact["contrast"]
    sparse_metrics = response_metrics(empirical_response, sparse_models[mode_index])
    dense_metrics = response_metrics(empirical_response, dense_model)
    removed_signal_norm = float(
        np.linalg.norm(original_local_cube[mode_index] - full_best_local)
    )
    full_local_fraction = (
        float(np.linalg.norm(exact_best_cube[mode_index] - full_best_local)) / removed_signal_norm
    )

    output_directory.mkdir(parents=True, exist_ok=True)
    overview_path = output_directory / "p4_response_overview.png"
    animation_path = output_directory / "p4_final_cube_and_psf.gif"
    summary_path = output_directory / "p4_response_visualization.json"
    write_overview(
        overview_path,
        original_local_cube[mode_index],
        exact_best_cube[mode_index],
        full_best_local,
        empirical_response,
        sparse_models[mode_index],
        dense_model,
        source_offset,
        float(modes[mode_index]),
        exact["contrast"],
        sparse_metrics,
        dense_metrics,
        full_local_fraction,
    )
    if not arguments.skip_animation:
        write_cube_animation(
            animation_path,
            modes,
            original_local_cube,
            exact_best_cube,
            sparse_models,
            exact["contrast"],
            source_offset,
        )

    summary = {
        "schema": 1,
        "experiment": str(experiment),
        "mode_fraction": float(modes[mode_index]),
        "mode_index": mode_index,
        "exact_planet": {
            **exact,
            "row": exact_row,
            "column": exact_column,
        },
        "response_coordinate": {"row": target_row, "column": target_column},
        "sparse_response": sparse_metrics,
        "dense_signal_free_response": dense_metrics,
        "exact_best_vs_full_best_fraction_of_removed_signal": full_local_fraction,
        "outputs": {
            "overview": str(overview_path),
            "cube_animation": None if arguments.skip_animation else str(animation_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {overview_path}")
    if not arguments.skip_animation:
        print(f"wrote {animation_path}")
    print(f"wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
