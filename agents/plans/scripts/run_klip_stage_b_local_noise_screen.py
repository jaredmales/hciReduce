#!/usr/bin/env python3
"""Run the raw-pixel mode-200 KLIP Stage-B noise-only screen locally."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_klip_stage_b_psd_extension as extension  # noqa: E402
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402


TRAINING_SUPPORT = 11
SUPPORTS = (11, 31, 47)
RADII = (7.5, 10.0, 12.0, 16.0, 20.0, 24.0)
PRIMARY_RADII = (7.5, 10.0, 12.0)
VARIANTS = (("identity", "rectangular", 1.0),
            ("rectangular_m0.1", "rectangular", 0.1),
            ("rectangular_m0.3", "rectangular", 0.3),
            ("hann_m0.1", "hann", 0.1),
            ("hann_m0.3", "hann", 0.3))


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def verify_bundle(path: Path) -> tuple[dict[str, object], dict[str, np.ndarray], dict[str, object]]:
    """Verify the transferred bundle receipt and return copied arrays."""
    receipt_path = path.with_suffix(path.suffix + ".receipt.json")
    receipt = read(receipt_path)
    observed = stage.fingerprint(path)
    expected = receipt["bundle"]
    stage.require(receipt["status"] == "complete" and observed["bytes"] == expected["bytes"] and
                  observed["sha256"] == expected["sha256"], "local bundle does not match its ROC receipt")
    with np.load(path, allow_pickle=False) as bundle:
        metadata = json.loads(bytes(bundle["metadata_json"]).decode("utf-8"))
        arrays = {name: np.array(bundle[name], copy=True)
                  for name in ("baseline", "positions", "source_indices", "response", "validity")}
    stage.require(metadata["mode"] == 200 and metadata["supports"] == list(SUPPORTS) and
                  metadata["radii"] == list(RADII) and arrays["baseline"].shape == (128, 128) and
                  arrays["response"].shape == arrays["validity"].shape ==
                  (metadata["position_count"], 47, 47), "unexpected local bundle schema")
    return metadata, arrays, receipt


def planet_position(shape: tuple[int, int], metadata: dict[str, object]) -> tuple[float, float]:
    """Convert the frozen planet polar coordinates to row and column coordinates."""
    center_row = 0.5 * (shape[1] - 1)
    center_column = 0.5 * (shape[0] - 1)
    planet = metadata["known_planet"]
    angle = math.radians(float(planet["position_angle"]))
    separation = float(planet["separation"])
    return center_row - separation * math.sin(angle), center_column + separation * math.cos(angle)


def training_stencils(finite: np.ndarray, query: tuple[int, int],
                      searches: list[tuple[int, int]], planet: tuple[float, float],
                      exclusion_support: int) -> dict[int, dict[str, object]]:
    """Cache the exact accepted 11-pixel interpolation stencils used by the preflight."""
    row_query, column_query = query
    center_row = 0.5 * (finite.shape[1] - 1)
    center_column = 0.5 * (finite.shape[0] - 1)
    radius = math.hypot(row_query - center_row, column_query - center_column)
    angle = math.atan2(column_query - center_column, row_query - center_row)
    half = TRAINING_SUPPORT // 2
    exclusion_half = exclusion_support // 2
    delta_column, delta_row = np.mgrid[-half:half + 1, -half:half + 1]
    minimum_k = math.ceil((6.0 - radius) / half)
    maximum_k = math.floor((np.nextafter(60.0, 0.0) - radius) / half)
    rings = {}
    for multiplier in range(minimum_k, maximum_k + 1):
        offset = multiplier * half
        ring_radius = radius + offset
        count = math.ceil(2 * math.pi * ring_radius / half)
        theta = np.arange(count, dtype=np.float64) * (2 * math.pi / count)
        centers_row = center_row + ring_radius * np.cos(theta)
        centers_column = center_column + ring_radius * np.sin(theta)
        cosine = np.cos(theta - angle)[:, None]
        sine = np.sin(theta - angle)[:, None]
        sample_row = centers_row[:, None] + cosine * delta_row.ravel() - sine * delta_column.ravel()
        sample_column = centers_column[:, None] + sine * delta_row.ravel() + cosine * delta_column.ravel()
        row0 = np.floor(sample_row).astype(int)
        column0 = np.floor(sample_column).astype(int)
        row_fraction = sample_row - row0
        column_fraction = sample_column - column0
        excluded = np.zeros(count, dtype=bool)
        incomplete = np.zeros(count, dtype=bool)
        indices, weights, used_columns = [], [], []
        for column_delta in (0, 1):
            for row_delta in (0, 1):
                weight = ((row_fraction if row_delta else 1 - row_fraction) *
                          (column_fraction if column_delta else 1 - column_fraction))
                used = weight != 0
                native_row = row0 + row_delta
                native_column = column0 + column_delta
                outside = ((native_row < 0) | (native_row >= finite.shape[1]) |
                           (native_column < 0) | (native_column >= finite.shape[0]))
                safe_row = np.clip(native_row, 0, finite.shape[1] - 1)
                safe_column = np.clip(native_column, 0, finite.shape[0] - 1)
                bad = outside | ~finite[safe_column, safe_row]
                bad |= (np.hypot(native_row - planet[0], native_column - planet[1]) <=
                        stage.PLANET_EXCLUSION_RADIUS)
                for search_row, search_column in searches:
                    bad |= ((np.abs(native_row - search_row) <= exclusion_half) &
                            (np.abs(native_column - search_column) <= exclusion_half))
                excluded |= np.any(used & bad & ~outside, axis=1)
                incomplete |= np.any(used & (outside | ~finite[safe_column, safe_row]), axis=1)
                indices.append(safe_column * finite.shape[1] + safe_row)
                weights.append(weight)
                used_columns.append(np.where(used, native_column, np.nan))
        accepted = ~(excluded | incomplete)
        columns = np.stack(used_columns, axis=-1)[accepted]
        halves = np.where(np.nanmax(columns, axis=(1, 2)) < center_column, 0,
                          np.where(np.nanmin(columns, axis=(1, 2)) > center_column, 1, -1))
        rings[offset] = {
            "indices": np.stack(indices, axis=-1)[accepted],
            "weights": np.stack(weights, axis=-1)[accepted],
            "halves": halves,
            "counts": {"attempted": count, "accepted": int(np.count_nonzero(accepted)),
                       "top": int(np.count_nonzero(halves == 0)),
                       "bottom": int(np.count_nonzero(halves == 1)),
                       "straddling": int(np.count_nonzero(halves < 0))},
        }
    return rings


def extract(image: np.ndarray, ring: dict[str, object]) -> np.ndarray:
    """Extract all accepted patches from a cached bilinear stencil."""
    weights = np.asarray(ring["weights"])
    indices = np.asarray(ring["indices"])
    values = np.where(weights != 0, image.ravel()[indices], 0)
    return np.sum(weights * values, axis=2)


def counts_for_rings(rings: dict[int, dict[str, object]]) -> dict[str, dict[str, object]]:
    """Return preflight-compatible band counts for cached stencils."""
    return footprint.band_counts({offset: dict(ring["counts"]) for offset, ring in rings.items()},
                                 TRAINING_SUPPORT)


def band_samples(image: np.ndarray, rings: dict[int, dict[str, object]],
                 width: int, half: int) -> tuple[np.ndarray, float]:
    """Extract one detector half of one radial training band and its interpolation gain."""
    samples, gains = [], []
    for offset, ring in rings.items():
        if abs(offset) > width:
            continue
        selected = np.asarray(ring["halves"]) == half
        if not np.any(selected):
            continue
        samples.append(extract(image, ring)[selected])
        gains.append(np.sum(np.square(np.asarray(ring["weights"])[selected]), axis=2).ravel())
    stage.require(samples, "training half contains no accepted patches")
    return np.vstack(samples), float(np.mean(np.concatenate(gains)))


def fit_periodogram_base(samples: np.ndarray, response_size: int, window_name: str) -> dict[str, object]:
    """Fit the reusable unmixed part of one extended Welch spectrum."""
    samples = np.asarray(samples, dtype=np.float64)
    stage.require(samples.ndim == 2 and samples.shape[1] == 121 and len(samples) >= 8 and
                  np.all(np.isfinite(samples)), "invalid local PSD training samples")
    mean = np.mean(samples, axis=0)
    centered = samples - mean
    target = float(np.sum(np.square(centered)) / ((len(samples) - 1) * 121))
    stage.require(np.isfinite(target) and target > 0, "invalid local training variance")
    window = extension.estimation_window(window_name)
    window_energy = float(np.sum(np.square(window)))
    fft_size = 2 * response_size - 1
    transformed = np.fft.fft2(centered.reshape((-1, 11, 11)) * window,
                              s=(fft_size, fft_size), axes=(-2, -1))
    raw_power = np.sum(np.square(np.abs(transformed)), axis=0) / ((len(samples) - 1) * window_energy)
    raw_zero_lag = float(np.mean(raw_power))
    stage.require(np.isfinite(raw_zero_lag) and raw_zero_lag > 0, "invalid local periodogram")
    return {"mean": mean, "target_variance": target, "window": window_name,
            "window_energy": window_energy, "mixing": None, "fft_size": fft_size,
            "raw_windowed_zero_lag": raw_zero_lag, "psd_rescaling": target / raw_zero_lag,
            "rescaled_power": raw_power * (target / raw_zero_lag)}


def regularize(base: dict[str, object], mixing: float) -> dict[str, object]:
    """Apply an isotropic spectral mixture to a cached periodogram."""
    target = float(base["target_variance"])
    power = (1 - mixing) * np.asarray(base["rescaled_power"]) + mixing * target
    lag = np.fft.ifft2(power)
    stage.require(np.max(np.abs(lag.imag)) <= 1e-12 * target and
                  np.min(power) >= mixing * target * (1 - 1e-12),
                  "local PSD regularization lost its positive real contract")
    return {key: value for key, value in base.items() if key != "rescaled_power"} | {
        "mixing": mixing, "power": power, "lag": lag.real}


def score_stats(values: np.ndarray) -> dict[str, object]:
    """Summarize one finite standardized-score vector."""
    values = np.asarray(values, dtype=np.float64)
    stage.require(values.ndim == 1 and len(values) >= 2 and np.all(np.isfinite(values)),
                  "invalid score vector")
    return {"count": len(values), "mean": float(np.mean(values)),
            "variance": float(np.var(values, ddof=1)), "mse": float(np.mean(np.square(values))),
            "within_one_fraction": float(np.mean(np.abs(values) <= 1)),
            "within_two_fraction": float(np.mean(np.abs(values) <= 2)),
            "quantiles": np.percentile(values, [0, 10, 50, 90, 100]).tolist(),
            "maximum_absolute": float(np.max(np.abs(values)))}


def analyze_model(model: dict[str, object], template: np.ndarray, candidate: np.ndarray,
                  heldout: np.ndarray | None) -> tuple[dict[str, object], dict[str, object]]:
    """Solve and score one candidate, retaining vectors only for split comparison."""
    solved = extension.solve_template(model, template)
    weight = np.asarray(solved["weight"])
    vector = candidate.ravel()
    sigma = float(solved["sigma"])
    stage.require(np.isclose(weight @ template.ravel(), 1, rtol=2e-10, atol=2e-12),
                  "local matched-filter weight lost unit response")
    amplitude = float(weight @ vector)
    record = {"amplitude_no_mean": amplitude, "score_no_mean": amplitude / sigma,
              "sigma": sigma, "iterations": solved["iterations"],
              "relative_residual": solved["relative_residual"],
              "spectral_condition_number": float(np.max(model["power"]) / np.min(model["power"])),
              "target_variance": model["target_variance"],
              "psd_rescaling": model["psd_rescaling"]}
    if template.shape == (11, 11):
        mean = np.asarray(model["mean"])
        mean_amplitude = float(weight @ (vector - mean))
        record.update({"amplitude_fitted_mean": mean_amplitude,
                       "score_fitted_mean": mean_amplitude / sigma})
        if heldout is not None:
            record["heldout_no_mean"] = score_stats(heldout @ weight / sigma)
            record["heldout_fitted_mean"] = score_stats((heldout - mean) @ weight / sigma)
    detail = {"weight": weight, "score_no_mean": record["score_no_mean"],
              "amplitude_no_mean": amplitude, "sigma": sigma}
    if "score_fitted_mean" in record:
        detail["score_fitted_mean"] = record["score_fitted_mean"]
    return record, detail


def crop(array: np.ndarray, support: int) -> np.ndarray:
    """Return the central square support from a 47-pixel stamp."""
    center = array.shape[0] // 2
    half = support // 2
    return array[center - half:center + half + 1, center - half:center + half + 1]


def summarize_records(records: list[dict[str, object]]) -> dict[str, object]:
    """Aggregate candidate calibration and split stability by fixed policy."""
    groups = {}
    for support in SUPPORTS:
        for radius in RADII:
            for band_role in ("narrow", "full"):
                for method, _, _ in VARIANTS:
                    chosen = [row for row in records if row["support"] == support and
                              row["radius"] == radius and row["band_role"] == band_role and
                              row["method"] == method]
                    stage.require(len(chosen) == 60, "fixed local summary lost query records")
                    key = f"s{support}_r{radius:g}_{band_role}_{method}"
                    directional = np.asarray([value for row in chosen for value in row["score_no_mean"]])
                    centers = np.asarray([value for row in chosen if row["search_index"] == 0
                                          for value in row["score_no_mean"]])
                    maximums = []
                    for site in range(12):
                        site_rows = sorted((row for row in chosen if row["site_index"] == site),
                                           key=lambda row: row["search_index"])
                        stage.require(len(site_rows) == 5, "five-pixel search is incomplete")
                        for half in (0, 1):
                            maximums.append(max(row["score_no_mean"][half] for row in site_rows))
                    first = np.asarray([row["score_no_mean"][0] for row in chosen])
                    second = np.asarray([row["score_no_mean"][1] for row in chosen])
                    entry = {"support": support, "radius": radius, "band_role": band_role,
                             "band_half_width": chosen[0]["band_half_width"], "method": method,
                             "directional_candidate_scores_no_mean": score_stats(directional),
                             "directional_center_scores_no_mean": score_stats(centers),
                             "five_pixel_maximum_no_mean": score_stats(np.asarray(maximums)),
                             "median_split_weight_cosine": float(np.median([row["weight_cosine"] for row in chosen])),
                             "split_weight_cosine_10_90": np.percentile(
                                 [row["weight_cosine"] for row in chosen], [10, 90]).tolist(),
                             "split_score_correlation": float(np.corrcoef(first, second)[0, 1]),
                             "median_absolute_split_score_difference": float(np.median(np.abs(first - second))),
                             "median_sigma_ratio": float(np.median([
                                 max(row["sigma"]) / min(row["sigma"]) for row in chosen])),
                             "training_patch_count_range": [min(value for row in chosen for value in row["samples"]),
                                                            max(value for row in chosen for value in row["samples"])],
                             "median_interpolation_variance_gain": float(np.median([
                                 value for row in chosen for value in row["interpolation_variance_gain"]])),
                             "maximum_solver_relative_residual": float(max(
                                 value for row in chosen for value in row["relative_residual"])),
                             "maximum_solver_iterations": int(max(
                                 value for row in chosen for value in row["iterations"]))}
                    if support == 11:
                        directional_mean = np.asarray([value for row in chosen
                                                       for value in row["score_fitted_mean"]])
                        entry["directional_candidate_scores_fitted_mean"] = score_stats(directional_mean)
                        for label in ("heldout_no_mean", "heldout_fitted_mean"):
                            projections = [item for row in chosen for item in row[label]]
                            entry[label] = {
                                "directional_fits": len(projections),
                                "median_variance_over_prediction": float(np.median([
                                    item["variance"] for item in projections])),
                                "variance_10_90": np.percentile([
                                    item["variance"] for item in projections], [10, 90]).tolist(),
                                "median_mse_over_prediction": float(np.median([
                                    item["mse"] for item in projections])),
                                "median_within_one_fraction": float(np.median([
                                    item["within_one_fraction"] for item in projections])),
                                "median_within_two_fraction": float(np.median([
                                    item["within_two_fraction"] for item in projections])),
                            }
                    groups[key] = entry
    return {"groups": groups, "records": len(records)}


def pooled_primary(summary: dict[str, object]) -> list[dict[str, object]]:
    """Pool the controlling radii for a compact primary comparison."""
    rows = []
    groups = summary["groups"]
    for support in SUPPORTS:
        for band_role in ("narrow", "full"):
            for method, _, _ in VARIANTS:
                selected = [groups[f"s{support}_r{radius:g}_{band_role}_{method}"]
                            for radius in PRIMARY_RADII]
                rows.append({"support": support, "band_role": band_role, "method": method,
                             "score_variance_median": float(np.median([
                                 row["directional_candidate_scores_no_mean"]["variance"] for row in selected])),
                             "absolute_score_mean_median": float(np.median([
                                 abs(row["directional_candidate_scores_no_mean"]["mean"]) for row in selected])),
                             "within_one_fraction_median": float(np.median([
                                 row["directional_candidate_scores_no_mean"]["within_one_fraction"] for row in selected])),
                             "split_weight_cosine_median": float(np.median([
                                 row["median_split_weight_cosine"] for row in selected])),
                             "split_score_correlation_median": float(np.median([
                                 row["split_score_correlation"] for row in selected])),
                             "five_pixel_maximum_median": float(np.median([
                                 row["five_pixel_maximum_no_mean"]["quantiles"][2] for row in selected])),
                             "heldout_patch_variance_median": (float(np.median([
                                 row["heldout_no_mean"]["median_variance_over_prediction"]
                                 for row in selected])) if support == 11 else None)})
    return rows


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the compact human-readable comparison and one diagnostic figure."""
    lines = ["# KLIP Stage-B local mode-200 noise screen", "",
             "This baseline-only screen uses exact response templates, raw 11-by-11 Welch training patches, "
             "candidate-specific five-footprint exclusions, and independent detector-half fits. Scores do not "
             "subtract a fitted candidate mean in the common 11/31/47 comparison.", "",
             "## Controlling radii: 7.5, 10, and 12 pixels", "",
             "| Response | Band | Method | Score variance | |Score mean| | Within ±1σ | Split weight cosine | Split score corr. | Median 5-pixel max | 11px opposite-half variance |",
             "| ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in result["primary"]:
        heldout = "—" if row["heldout_patch_variance_median"] is None else f"{row['heldout_patch_variance_median']:.3f}"
        lines.append(f"| {row['support']} | {row['band_role']} | {row['method']} | "
                     f"{row['score_variance_median']:.3f} | {row['absolute_score_mean_median']:.3f} | "
                     f"{row['within_one_fraction_median']:.3f} | {row['split_weight_cosine_median']:.3f} | "
                     f"{row['split_score_correlation_median']:.3f} | {row['five_pixel_maximum_median']:.3f} | {heldout} |")
    lines.extend(["", "The score-variance columns use 120 directional scores per radius (60 fixed queries × two "
                  "fits) and remain spatially correlated. The five-pixel maximum uses 24 directional searches per "
                  "radius. The 11-pixel opposite-half column is the median, across query-specific fits, of the "
                  "variance of genuinely disjoint held-out patches divided by the conditional prediction.", "",
                  "For 31- and 47-pixel responses, complete response-sized held-out training stamps do not exist "
                  "under the frozen exclusion geometry. Their validation here is therefore limited to fixed "
                  "candidate-null calibration, split-weight stability, and split-score agreement.", "",
                  "This first screen covers raw pixels and no candidate-mean subtraction. Radial-variance "
                  "standardization, post-mean patch-RMS normalization, the direct 11-pixel PCA/diagonal grid, and "
                  "the support-independent radial-mean control remain separate Stage-B arms."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    methods = [item[0] for item in VARIANTS]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), layout="constrained")
    for column, support in enumerate(SUPPORTS):
        chosen = [row for row in result["primary"] if row["support"] == support]
        variance = np.asarray([[next(row["score_variance_median"] for row in chosen
                                     if row["band_role"] == band and row["method"] == method)
                                for method in methods] for band in ("narrow", "full")])
        cosine = np.asarray([[next(row["split_weight_cosine_median"] for row in chosen
                                   if row["band_role"] == band and row["method"] == method)
                              for method in methods] for band in ("narrow", "full")])
        for axis, values, title, limits, cmap in ((axes[0, column], variance, "score variance", (0, 3), "RdBu_r"),
                                                  (axes[1, column], cosine, "split weight cosine", (0, 1), "Blues")):
            mesh = axis.imshow(values, vmin=limits[0], vmax=limits[1], aspect="auto", cmap=cmap)
            axis.set(title=f"{support}px response: {title}", yticks=(0, 1), yticklabels=("narrow", "full"),
                     xticks=np.arange(len(methods)), xticklabels=methods)
            axis.tick_params(axis="x", rotation=35, labelsize=8)
            for y in range(2):
                for x in range(len(methods)):
                    axis.text(x, y, f"{values[y, x]:.2f}", ha="center", va="center", fontsize=8)
            fig.colorbar(mesh, ax=axis, shrink=0.72)
    fig.suptitle("KLIP mode 200 baseline: controlling radii pooled by median")
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, root: Path) -> None:
    """Run the frozen local screen and write replayable diagnostics."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = verify_bundle(bundle_path)
    root.mkdir(parents=True)
    scripts = [stage.fingerprint(Path(__file__)),
               stage.fingerprint(Path(__file__).with_name("check_klip_stage_b_psd_extension.py")),
               stage.fingerprint(Path(__file__).with_name("run_klip_stage_b_footprint_preflight.py"))]
    protocol = {"schema": 1, "purpose": "local raw-pixel mode-200 KLIP Stage-B noise-only screen",
                "mode": 200, "supports": list(SUPPORTS), "radii": list(RADII),
                "primary_radii": list(PRIMARY_RADII), "training_support": TRAINING_SUPPORT,
                "training_bands": "narrowest split-supported width per support/radius plus full width 60",
                "variants": [{"name": name, "window": window, "isotropic_mixing": mixing}
                             for name, window, mixing in VARIANTS],
                "candidate_mean": "off for common support; fitted ensemble mean also reported at support 11",
                "normalization": "raw pixels; no radial or per-patch RMS normalization",
                "split": "disjoint native detector-column halves; each fit scores the same held-out candidate",
                "bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
                "source_response_complete_sha256": metadata["response_complete_sha256"],
                "source_preflight_complete_sha256": metadata["preflight_complete_sha256"],
                "scripts": scripts}
    stage.write_json(root / "protocol.json", protocol)
    stage.write_json(root / "state.json", {"status": "running", "started_unix": time.time()})

    baseline = np.asarray(arrays["baseline"], dtype=np.float64)
    finite = np.isfinite(baseline)
    planet = planet_position(baseline.shape, metadata)
    lookup = {(int(row), int(column)): index for index, (row, column) in enumerate(arrays["positions"])}
    records = []
    geometry_checks = 0
    started = time.monotonic()
    for support in SUPPORTS:
        for radius in RADII:
            geometry = metadata["geometry"][str(support)][str(radius)]
            narrow = int(geometry["training"]["narrowest_split_supported_half_width"])
            count_records = []
            for site_index, site in enumerate(geometry["selected_sites"]):
                searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                            for delta_row, delta_column in footprint.SEARCH_OFFSETS]
                for search_index, query in enumerate(searches):
                    rings = training_stencils(finite, query, searches, planet, support)
                    observed_bands = counts_for_rings(rings)
                    oracle = footprint.band_counts(footprint.training_rings(
                        finite, query, searches, TRAINING_SUPPORT, planet,
                        exclusion_support=support), TRAINING_SUPPORT)
                    stage.require(observed_bands == oracle, "cached training stencil differs from preflight oracle")
                    count_records.append({"bands": observed_bands})
                    geometry_checks += 1

                    source = lookup[query]
                    canonical_response = np.asarray(arrays["response"][source], dtype=np.float64).T
                    canonical_validity = np.asarray(arrays["validity"][source], dtype=bool).T
                    template = crop(canonical_response, support)
                    stage.require(np.all(crop(canonical_validity, support)) and np.all(np.isfinite(template)) and
                                  np.sum(np.square(template)) > 0, "invalid exact response crop")
                    half_support = support // 2
                    row, column = query
                    candidate = baseline[column - half_support:column + half_support + 1,
                                         row - half_support:row + half_support + 1].T
                    stage.require(candidate.shape == template.shape and np.all(np.isfinite(candidate)),
                                  "invalid candidate data crop")
                    for band_role, width in (("narrow", narrow), ("full", 60)):
                        halves = []
                        for detector_half in (0, 1):
                            samples, interpolation_gain = band_samples(
                                baseline, rings, width, detector_half)
                            rectangular = fit_periodogram_base(samples, support, "rectangular")
                            hann = fit_periodogram_base(samples, support, "hann")
                            models = {"identity": regularize(rectangular, 1.0),
                                      "rectangular_m0.1": regularize(rectangular, 0.1),
                                      "rectangular_m0.3": regularize(rectangular, 0.3),
                                      "hann_m0.1": regularize(hann, 0.1),
                                      "hann_m0.3": regularize(hann, 0.3)}
                            opposite, _ = band_samples(baseline, rings, width, 1 - detector_half)
                            model_results = {}
                            model_details = {}
                            for method, _, _ in VARIANTS:
                                result, detail = analyze_model(models[method], template, candidate,
                                                               opposite if support == 11 else None)
                                model_results[method] = result
                                model_details[method] = detail
                            halves.append({"samples": len(samples),
                                           "interpolation_variance_gain": interpolation_gain,
                                           "results": model_results, "details": model_details})
                        for method, _, _ in VARIANTS:
                            first = halves[0]["details"][method]
                            second = halves[1]["details"][method]
                            weight_cosine = float(first["weight"] @ second["weight"] /
                                                  (np.linalg.norm(first["weight"]) *
                                                   np.linalg.norm(second["weight"])))
                            row_record = {"support": support, "radius": radius,
                                          "site_index": site_index, "search_index": search_index,
                                          "row": query[0], "column": query[1],
                                          "band_role": band_role, "band_half_width": width,
                                          "method": method,
                                          "samples": [halves[index]["samples"] for index in (0, 1)],
                                          "interpolation_variance_gain": [halves[index]["interpolation_variance_gain"]
                                                                           for index in (0, 1)],
                                          "weight_cosine": weight_cosine,
                                          "score_no_mean": [halves[index]["results"][method]["score_no_mean"]
                                                            for index in (0, 1)],
                                          "amplitude_no_mean": [halves[index]["results"][method]["amplitude_no_mean"]
                                                                for index in (0, 1)],
                                          "sigma": [halves[index]["results"][method]["sigma"]
                                                    for index in (0, 1)],
                                          "iterations": [halves[index]["results"][method]["iterations"]
                                                         for index in (0, 1)],
                                          "relative_residual": [halves[index]["results"][method]["relative_residual"]
                                                                for index in (0, 1)],
                                          "spectral_condition_number": [
                                              halves[index]["results"][method]["spectral_condition_number"]
                                              for index in (0, 1)]}
                            if support == 11:
                                row_record["score_fitted_mean"] = [
                                    halves[index]["results"][method]["score_fitted_mean"] for index in (0, 1)]
                                row_record["heldout_no_mean"] = [
                                    halves[index]["results"][method]["heldout_no_mean"] for index in (0, 1)]
                                row_record["heldout_fitted_mean"] = [
                                    halves[index]["results"][method]["heldout_fitted_mean"] for index in (0, 1)]
                            records.append(row_record)
            observed_summary = footprint.geometry_summary(count_records, TRAINING_SUPPORT)
            stage.require(observed_summary == geometry["training"],
                          f"saved preflight geometry did not replay for support {support}, radius {radius}")
            print(f"completed support {support}, radius {radius:g}: {len(records)} records", flush=True)

    summary = summarize_records(records)
    primary = pooled_primary(summary)
    result = {"purpose": protocol["purpose"], "mode": 200, "records": len(records),
              "geometry_query_checks": geometry_checks, "elapsed_seconds": time.monotonic() - started,
              "primary": primary, "groups": summary["groups"],
              "limitations": ["single signal-free residual field", "spatially correlated fixed null queries",
                              "raw-pixel arm only", "31/47-pixel held-out evidence uses candidate nulls rather than complete off-candidate stamps"]}
    stage.write_json(root / "records.json", records)
    stage.write_json(root / "results.json", result)
    write_report(root, result)
    completion = {"status": "complete", "protocol": stage.fingerprint(root / "protocol.json"),
                  "records": stage.fingerprint(root / "records.json"),
                  "results": stage.fingerprint(root / "results.json"),
                  "report": stage.fingerprint(root / "results.md"),
                  "figure": stage.fingerprint(root / "comparison.png")}
    stage.write_json(root / "complete.json", completion)
    stage.write_json(root / "state.json", {"status": "complete", "elapsed_seconds": result["elapsed_seconds"]})
    print(root / "results.md", flush=True)


def check() -> None:
    """Check cached PSD construction against the numerical extension contract."""
    generator = np.random.default_rng(92031)
    samples = generator.normal(size=(32, 121))
    template = generator.normal(size=(31, 31))
    for window in ("rectangular", "hann"):
        base = fit_periodogram_base(samples, 31, window)
        for mixing in (0.1, 0.3, 1.0):
            observed = regularize(base, mixing)
            expected = extension.fit_extended_psd(samples, 31, window, mixing)
            stage.require(np.array_equal(observed["mean"], expected["mean"]) and
                          np.allclose(observed["power"], expected["power"], rtol=2e-15, atol=1e-15) and
                          np.allclose(observed["lag"], expected["lag"], rtol=2e-15, atol=1e-15),
                          "cached local PSD differs from the verified extension")
            first = extension.solve_template(observed, template)
            second = extension.solve_template(expected, template)
            stage.require(np.allclose(first["weight"], second["weight"], rtol=2e-12, atol=2e-14),
                          "cached local PSD changes matched-filter weights")
    print("KLIP Stage-B local noise-screen checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the local noise-screen command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the requested local noise-screen action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        run(args.bundle.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
