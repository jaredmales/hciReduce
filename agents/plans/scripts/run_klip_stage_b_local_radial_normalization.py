#!/usr/bin/env python3
"""Test strict candidate-excluded radial normalization in the local KLIP screen."""
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
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402
import run_klip_stage_b_local_noise_screen as raw  # noqa: E402


SUPPORT = 11
RADIAL_BIN_WIDTH = 3.6
MINIMUM_PROFILE_PIXELS = 20


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def radius_map(shape: tuple[int, int]) -> np.ndarray:
    """Return native-pixel radius in the baseline's column-by-row layout."""
    column, row = np.indices(shape)
    return np.hypot(row - 0.5 * (shape[1] - 1), column - 0.5 * (shape[0] - 1))


def exclusion_mask(shape: tuple[int, int], searches: list[tuple[int, int]],
                   support: int, planet: tuple[float, float]) -> np.ndarray:
    """Mask the known planet and the union of five complete candidate footprints."""
    column, row = np.indices(shape)
    excluded = np.hypot(row - planet[0], column - planet[1]) <= stage.PLANET_EXCLUSION_RADIUS
    half = support // 2
    for search_row, search_column in searches:
        excluded |= ((np.abs(row - search_row) <= half) &
                     (np.abs(column - search_column) <= half))
    return excluded


def profile_counts(image: np.ndarray, excluded: np.ndarray) -> np.ndarray:
    """Count available pixels in the fixed radial-variance bins."""
    radius = radius_map(image.shape)
    allowed = np.isfinite(image) & ~excluded & (radius < 60)
    return np.asarray([np.count_nonzero(allowed & (radius >= lower) &
                                                  (radius < min(float(lower + RADIAL_BIN_WIDTH), 60)))
                       for lower in np.arange(0, 60, RADIAL_BIN_WIDTH)], dtype=int)


def variance_profile(image: np.ndarray, excluded: np.ndarray) -> tuple[np.ndarray, list[dict[str, object]]]:
    """Fit the frozen 3.6-pixel log-variance interpolation outside one held-out site."""
    radius = radius_map(image.shape)
    allowed = np.isfinite(image) & ~excluded & (radius < 60)
    profile = []
    for lower in np.arange(0, 60, RADIAL_BIN_WIDTH):
        upper = min(float(lower + RADIAL_BIN_WIDTH), 60)
        pixels = image[allowed & (radius >= lower) & (radius < upper)]
        stage.require(len(pixels) >= MINIMUM_PROFILE_PIXELS, "unsupported radial variance bin")
        variance = float(np.var(pixels, ddof=1))
        stage.require(np.isfinite(variance) and variance >= 0, "invalid radial variance bin")
        structural_zero = variance == 0
        if structural_zero:
            stage.require(lower == 0 and np.all(pixels == 0),
                          "zero variance is not the exact KLIP structural core")
        profile.append({"lower": float(lower), "upper": upper, "pixels": len(pixels),
                        "variance": variance, "standard_deviation": math.sqrt(variance),
                        "structural_zero": structural_zero})
    positive = [item for item in profile if not item["structural_zero"]]
    stage.require(positive and all(item["variance"] > 0 for item in positive),
                  "radial profile has no positive supported bins")
    centers = [(item["lower"] + item["upper"]) / 2 for item in positive]
    scale = np.sqrt(np.exp(np.interp(radius, centers,
                                     np.log([item["variance"] for item in positive]))))
    scale[radius >= 60] = np.nan
    return scale, profile


def score_stats(values: np.ndarray) -> dict[str, object]:
    """Reuse the raw screen's finite standardized-score summary."""
    return raw.score_stats(np.asarray(values, dtype=np.float64))


def summarize(records: list[dict[str, object]], raw_result: dict[str, object]) -> dict[str, object]:
    """Compare normalized narrow-band results with the exact raw checkpoint."""
    groups = {}
    raw_groups = raw_result["groups"]
    for radius in raw.RADII:
        for method, _, _ in raw.VARIANTS:
            chosen = [row for row in records if row["radius"] == radius and row["method"] == method]
            stage.require(len(chosen) == 60, "normalized summary lost fixed query records")
            directional = np.asarray([value for row in chosen for value in row["score_no_mean"]])
            centers = np.asarray([value for row in chosen if row["search_index"] == 0
                                  for value in row["score_no_mean"]])
            maximums = []
            for site in range(12):
                site_rows = sorted((row for row in chosen if row["site_index"] == site),
                                   key=lambda row: row["search_index"])
                stage.require(len(site_rows) == 5, "normalized five-pixel search is incomplete")
                for detector_half in (0, 1):
                    maximums.append(max(row["score_no_mean"][detector_half] for row in site_rows))
            first = np.asarray([row["score_no_mean"][0] for row in chosen])
            second = np.asarray([row["score_no_mean"][1] for row in chosen])
            heldout = [item for row in chosen for item in row["heldout_no_mean"]]
            raw_group = raw_groups[f"s11_r{radius:g}_narrow_{method}"]
            normalized_stats = score_stats(directional)
            entry = {"support": SUPPORT, "radius": radius, "method": method,
                     "band_half_width": chosen[0]["band_half_width"],
                     "directional_candidate_scores_no_mean": normalized_stats,
                     "directional_center_scores_no_mean": score_stats(centers),
                     "five_pixel_maximum_no_mean": score_stats(np.asarray(maximums)),
                     "median_split_physical_weight_cosine": float(np.median([
                         row["physical_weight_cosine"] for row in chosen])),
                     "split_physical_weight_cosine_10_90": np.percentile([
                         row["physical_weight_cosine"] for row in chosen], [10, 90]).tolist(),
                     "split_score_correlation": float(np.corrcoef(first, second)[0, 1]),
                     "median_absolute_split_score_difference": float(np.median(np.abs(first - second))),
                     "heldout_no_mean": {
                         "directional_fits": len(heldout),
                         "median_variance_over_prediction": float(np.median([
                             item["variance"] for item in heldout])),
                         "variance_10_90": np.percentile([
                             item["variance"] for item in heldout], [10, 90]).tolist(),
                         "median_mse_over_prediction": float(np.median([item["mse"] for item in heldout])),
                         "median_within_one_fraction": float(np.median([
                             item["within_one_fraction"] for item in heldout]))},
                     "raw_comparison": {
                         "candidate_variance": raw_group["directional_candidate_scores_no_mean"]["variance"],
                         "candidate_variance_ratio": (normalized_stats["variance"] /
                                                       raw_group["directional_candidate_scores_no_mean"]["variance"]),
                         "center_variance": raw_group["directional_center_scores_no_mean"]["variance"],
                         "heldout_variance_median": raw_group["heldout_no_mean"]["median_variance_over_prediction"],
                         "split_weight_cosine_median": raw_group["median_split_weight_cosine"]}}
            groups[f"r{radius:g}_{method}"] = entry

    primary = []
    for method, _, _ in raw.VARIANTS:
        selected = [groups[f"r{radius:g}_{method}"] for radius in raw.PRIMARY_RADII]
        primary.append({"method": method,
                        "normalized_score_variance_median": float(np.median([
                            row["directional_candidate_scores_no_mean"]["variance"] for row in selected])),
                        "raw_score_variance_median": float(np.median([
                            row["raw_comparison"]["candidate_variance"] for row in selected])),
                        "normalized_over_raw_variance_median": float(np.median([
                            row["raw_comparison"]["candidate_variance_ratio"] for row in selected])),
                        "normalized_center_variance_median": float(np.median([
                            row["directional_center_scores_no_mean"]["variance"] for row in selected])),
                        "normalized_heldout_variance_median": float(np.median([
                            row["heldout_no_mean"]["median_variance_over_prediction"] for row in selected])),
                        "raw_heldout_variance_median": float(np.median([
                            row["raw_comparison"]["heldout_variance_median"] for row in selected])),
                        "normalized_split_weight_cosine_median": float(np.median([
                            row["median_split_physical_weight_cosine"] for row in selected])),
                        "raw_split_weight_cosine_median": float(np.median([
                            row["raw_comparison"]["split_weight_cosine_median"] for row in selected])),
                        "normalized_within_one_fraction_median": float(np.median([
                            row["directional_candidate_scores_no_mean"]["within_one_fraction"] for row in selected]))})
    return {"groups": groups, "primary": primary}


def profile_geometry(metadata: dict[str, object], baseline: np.ndarray,
                     planet: tuple[float, float]) -> dict[str, object]:
    """Audit strict profile support for all three response footprints."""
    result = {}
    for support in raw.SUPPORTS:
        records = []
        for radius in raw.RADII:
            for site_index, site in enumerate(metadata["geometry"][str(support)][str(radius)]["selected_sites"]):
                searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                            for delta_row, delta_column in footprint.SEARCH_OFFSETS]
                counts = profile_counts(baseline, exclusion_mask(baseline.shape, searches, support, planet))
                records.append({"radius": radius, "site_index": site_index,
                                "minimum_bin_pixels": int(np.min(counts)),
                                "supported": bool(np.all(counts >= MINIMUM_PROFILE_PIXELS)),
                                "counts": counts.tolist()})
        minimum = np.min(np.asarray([record["counts"] for record in records]), axis=0)
        result[str(support)] = {"sites": len(records),
                                "supported_sites": sum(record["supported"] for record in records),
                                "minimum_bin_counts": minimum.tolist(), "records": records}
    return result


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the radial-normalization result and comparison figure."""
    lines = ["# KLIP Stage-B strict radial-normalization screen", "",
             "This mode-200 test changes only the input coordinates of the 11-pixel narrow-band PSD arm. "
             "For each fixed site, a 3.6-pixel-bin radial variance profile excludes the known planet and the "
             "union of all five complete candidate footprints. The same leave-site-out profile scales native "
             "training pixels, candidate data, and exact responses by standard deviation before interpolation "
             "and filtering.", "", "## Controlling radii", "",
             "| Method | Raw score variance | Normalized score variance | Normalized/raw | Center variance | Opposite-half variance | Split physical-weight cosine | Within ±1σ |",
             "| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in result["primary"]:
        lines.append(f"| {row['method']} | {row['raw_score_variance_median']:.3f} | "
                     f"{row['normalized_score_variance_median']:.3f} | "
                     f"{row['normalized_over_raw_variance_median']:.3f} | "
                     f"{row['normalized_center_variance_median']:.3f} | "
                     f"{row['normalized_heldout_variance_median']:.3f} | "
                     f"{row['normalized_split_weight_cosine_median']:.3f} | "
                     f"{row['normalized_within_one_fraction_median']:.3f} |")
    geometry = result["profile_geometry"]
    lines.extend(["", "## Normalization geometry", "",
                  "| Response support | Supported fixed sites | Minimum pixels in any radial bin |",
                  "| ---: | ---: | ---: |"])
    for support in raw.SUPPORTS:
        row = geometry[str(support)]
        lines.append(f"| {support} | {row['supported_sites']}/{row['sites']} | "
                     f"{min(row['minimum_bin_counts'])} |")
    lines.extend(["", "Strict candidate-excluded profiles are supported for every 11-pixel site, with at "
                  "least 21 pixels in every bin. They fail at 50 of 72 31-pixel sites and all 47-pixel sites "
                  "because the larger held-out response footprints cover complete inner annuli. Those response "
                  "supports were not normalized using a profile that reads the candidate or extrapolates across "
                  "an unsupported radial interval.", "",
                  "The innermost 0--3.6-pixel KLIP bin is a structural zero: all 44 baseline pixels and all "
                  "selected exact-response occurrences are bitwise zero. It receives the first positive-bin "
                  "scale only after this gate; scaling a zero data/response coordinate cannot change its filter "
                  "contribution."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    methods = [item[0] for item in raw.VARIANTS]
    raw_values = [row["raw_score_variance_median"] for row in result["primary"]]
    normalized = [row["normalized_score_variance_median"] for row in result["primary"]]
    heldout_raw = [row["raw_heldout_variance_median"] for row in result["primary"]]
    heldout_normalized = [row["normalized_heldout_variance_median"] for row in result["primary"]]
    x = np.arange(len(methods))
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    width = 0.36
    axes[0].bar(x - width / 2, raw_values, width, label="raw")
    axes[0].bar(x + width / 2, normalized, width, label="radial standardized")
    axes[0].axhline(1, color="black", linestyle="--", linewidth=1)
    axes[0].set(title="Fixed candidate score variance", ylabel="median over radii 7.5, 10, 12",
                xticks=x, xticklabels=methods)
    axes[1].bar(x - width / 2, heldout_raw, width, label="raw")
    axes[1].bar(x + width / 2, heldout_normalized, width, label="radial standardized")
    axes[1].axhline(1, color="black", linestyle="--", linewidth=1)
    axes[1].set(title="Opposite-half patch variance / prediction", xticks=x, xticklabels=methods)
    for axis in axes:
        axis.tick_params(axis="x", rotation=30, labelsize=8)
        axis.grid(axis="y", alpha=0.2)
        axis.legend()
    fig.suptitle("KLIP mode 200, 11-pixel response, narrow candidate-excluded training bands")
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, raw_root: Path, root: Path) -> None:
    """Run the strict 11-pixel radial-normalization comparison."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = raw.verify_bundle(bundle_path)
    raw_completion = read(raw_root / "complete.json")
    stage.verify([raw_completion[key] for key in ("protocol", "records", "results", "report", "figure")])
    raw_result = read(raw_root / "results.json")
    root.mkdir(parents=True)
    scripts = [stage.fingerprint(Path(__file__)),
               stage.fingerprint(Path(__file__).with_name("run_klip_stage_b_local_noise_screen.py"))]
    protocol = {"schema": 1, "purpose": "strict candidate-excluded radial normalization of the local KLIP screen",
                "mode": 200, "response_support": SUPPORT, "radii": list(raw.RADII),
                "training_bands": "narrowest split-supported width only",
                "variance_profile": {"native_pixels": True, "bin_width": RADIAL_BIN_WIDTH,
                                     "radial_range": [0, 60], "minimum_pixels_per_bin": MINIMUM_PROFILE_PIXELS,
                                     "interpolation": "linear in log variance between bin centers",
                                     "endpoint": "constant inside endpoint half-bins; undefined at radius >= 60",
                                     "structural_core": "the exact-zero 0--3.6 pixel bin receives the first positive-bin scale after a zero-data/zero-response gate",
                                     "exclusion": "known planet plus union of five complete candidate footprints",
                                     "scope": "one leave-site-out profile shared by the site's five queries"},
                "coordinate_transform": "divide native image, candidate, and response by interpolated standard deviation before bilinear extraction",
                "bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
                "raw_parent": stage.fingerprint(raw_root / "complete.json"), "scripts": scripts}
    stage.write_json(root / "protocol.json", protocol)
    stage.write_json(root / "state.json", {"status": "running", "started_unix": time.time()})

    baseline = np.asarray(arrays["baseline"], dtype=np.float64)
    finite = np.isfinite(baseline)
    planet = raw.planet_position(baseline.shape, metadata)
    lookup = {(int(row), int(column)): index for index, (row, column) in enumerate(arrays["positions"])}
    geometry = profile_geometry(metadata, baseline, planet)
    stage.require(geometry["11"]["supported_sites"] == geometry["11"]["sites"] and
                  geometry["31"]["supported_sites"] < geometry["31"]["sites"] and
                  geometry["47"]["supported_sites"] == 0,
                  "strict profile-support boundary changed")
    absolute_radius = radius_map(baseline.shape)
    stage.require(np.array_equal(baseline[absolute_radius < RADIAL_BIN_WIDTH],
                                 np.zeros(np.count_nonzero(absolute_radius < RADIAL_BIN_WIDTH))),
                  "KLIP structural core is not exactly zero in the baseline")
    structural_response_pixels = 0
    for radius in raw.RADII:
        for site in metadata["geometry"][str(SUPPORT)][str(radius)]["selected_sites"]:
            for delta_row, delta_column in footprint.SEARCH_OFFSETS:
                query = (int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                source = lookup[query]
                template = raw.crop(np.asarray(arrays["response"][source], dtype=np.float64).T, SUPPORT)
                stamp_row, stamp_column = np.mgrid[query[0] - SUPPORT // 2:query[0] + SUPPORT // 2 + 1,
                                                   query[1] - SUPPORT // 2:query[1] + SUPPORT // 2 + 1]
                core = np.hypot(stamp_row - 63.5, stamp_column - 63.5) < RADIAL_BIN_WIDTH
                stage.require(np.array_equal(template[core], np.zeros(np.count_nonzero(core))),
                              "KLIP structural core is not exactly zero in an exact response")
                structural_response_pixels += int(np.count_nonzero(core))

    records, profiles = [], []
    started = time.monotonic()
    for radius in raw.RADII:
        fixed = metadata["geometry"][str(SUPPORT)][str(radius)]
        width = int(fixed["training"]["narrowest_split_supported_half_width"])
        for site_index, site in enumerate(fixed["selected_sites"]):
            searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                        for delta_row, delta_column in footprint.SEARCH_OFFSETS]
            excluded = exclusion_mask(baseline.shape, searches, SUPPORT, planet)
            scale, profile = variance_profile(baseline, excluded)
            standardized = baseline / scale
            profiles.append({"radius": radius, "site_index": site_index,
                             "row": int(site["row"]), "column": int(site["column"]), "bins": profile})
            for search_index, query in enumerate(searches):
                rings = raw.training_stencils(finite, query, searches, planet, SUPPORT)
                source = lookup[query]
                template = raw.crop(np.asarray(arrays["response"][source], dtype=np.float64).T, SUPPORT)
                validity = raw.crop(np.asarray(arrays["validity"][source], dtype=bool).T, SUPPORT)
                row, column = query
                half_support = SUPPORT // 2
                candidate = baseline[column - half_support:column + half_support + 1,
                                     row - half_support:row + half_support + 1].T
                candidate_scale = scale[column - half_support:column + half_support + 1,
                                        row - half_support:row + half_support + 1].T
                stage.require(np.all(validity) and np.all(np.isfinite(candidate_scale)) and
                              np.all(candidate_scale > 0), "invalid standardized candidate support")
                standardized_template = template / candidate_scale
                standardized_candidate = candidate / candidate_scale
                halves = []
                for detector_half in (0, 1):
                    samples, interpolation_gain = raw.band_samples(
                        standardized, rings, width, detector_half)
                    opposite, _ = raw.band_samples(standardized, rings, width, 1 - detector_half)
                    rectangular = raw.fit_periodogram_base(samples, SUPPORT, "rectangular")
                    hann = raw.fit_periodogram_base(samples, SUPPORT, "hann")
                    models = {"identity": raw.regularize(rectangular, 1.0),
                              "rectangular_m0.1": raw.regularize(rectangular, 0.1),
                              "rectangular_m0.3": raw.regularize(rectangular, 0.3),
                              "hann_m0.1": raw.regularize(hann, 0.1),
                              "hann_m0.3": raw.regularize(hann, 0.3)}
                    results, details = {}, {}
                    for method, _, _ in raw.VARIANTS:
                        result, detail = raw.analyze_model(models[method], standardized_template,
                                                          standardized_candidate, opposite)
                        results[method] = result
                        details[method] = detail | {"physical_weight": detail["weight"] /
                                                                     candidate_scale.ravel()}
                    halves.append({"samples": len(samples), "interpolation_variance_gain": interpolation_gain,
                                   "results": results, "details": details})
                for method, _, _ in raw.VARIANTS:
                    first = halves[0]["details"][method]
                    second = halves[1]["details"][method]
                    physical_cosine = float(first["physical_weight"] @ second["physical_weight"] /
                                            (np.linalg.norm(first["physical_weight"]) *
                                             np.linalg.norm(second["physical_weight"])))
                    records.append({"radius": radius, "site_index": site_index,
                                    "search_index": search_index, "row": row, "column": column,
                                    "band_half_width": width, "method": method,
                                    "samples": [halves[index]["samples"] for index in (0, 1)],
                                    "physical_weight_cosine": physical_cosine,
                                    "score_no_mean": [halves[index]["results"][method]["score_no_mean"]
                                                      for index in (0, 1)],
                                    "score_fitted_mean": [halves[index]["results"][method]["score_fitted_mean"]
                                                          for index in (0, 1)],
                                    "sigma": [halves[index]["results"][method]["sigma"] for index in (0, 1)],
                                    "heldout_no_mean": [halves[index]["results"][method]["heldout_no_mean"]
                                                        for index in (0, 1)],
                                    "heldout_fitted_mean": [
                                        halves[index]["results"][method]["heldout_fitted_mean"]
                                        for index in (0, 1)]})
        print(f"completed radius {radius:g}: {len(records)} records", flush=True)

    summarized = summarize(records, raw_result)
    profile_std = np.asarray([[item["standard_deviation"] for item in profile["bins"]]
                              for profile in profiles])
    result = {"purpose": protocol["purpose"], "mode": 200, "support": SUPPORT,
              "records": len(records), "elapsed_seconds": time.monotonic() - started,
              "profile_geometry": geometry,
              "structural_zero_core": {"radius_limit": RADIAL_BIN_WIDTH,
                                       "baseline_pixels": int(np.count_nonzero(absolute_radius < RADIAL_BIN_WIDTH)),
                                       "selected_response_pixel_occurrences": structural_response_pixels,
                                       "all_values_bitwise_zero": True,
                                       "assigned_scale": "first positive radial bin"},
              "profile_site_standard_deviation": {
                  "profiles": len(profiles), "minimum": np.min(profile_std, axis=0).tolist(),
                  "median": np.median(profile_std, axis=0).tolist(),
                  "maximum": np.max(profile_std, axis=0).tolist()},
              "primary": summarized["primary"], "groups": summarized["groups"]}
    stage.write_json(root / "profiles.json", profiles)
    stage.write_json(root / "records.json", records)
    stage.write_json(root / "results.json", result)
    write_report(root, result)
    completion = {"status": "complete", "protocol": stage.fingerprint(root / "protocol.json"),
                  "profiles": stage.fingerprint(root / "profiles.json"),
                  "records": stage.fingerprint(root / "records.json"),
                  "results": stage.fingerprint(root / "results.json"),
                  "report": stage.fingerprint(root / "results.md"),
                  "figure": stage.fingerprint(root / "comparison.png")}
    stage.write_json(root / "complete.json", completion)
    stage.write_json(root / "state.json", {"status": "complete", "elapsed_seconds": result["elapsed_seconds"]})
    print(root / "results.md", flush=True)


def check() -> None:
    """Check profile masking and consistent data/template standardization."""
    generator = np.random.default_rng(19723)
    image = generator.normal(size=(128, 128))
    searches = [(70, 60), (69, 60), (71, 60), (70, 59), (70, 61)]
    planet = (76.0, 62.0)
    excluded = exclusion_mask(image.shape, searches, SUPPORT, planet)
    scale, profile = variance_profile(image, excluded)
    stage.require(len(profile) == 17 and np.all(profile_counts(image, excluded) >= 20) and
                  np.all(np.isfinite(scale[radius_map(image.shape) < 60])),
                  "radial profile check failed")
    samples = generator.normal(size=(32, 121))
    template = generator.normal(size=(11, 11))
    factors = np.exp(generator.normal(scale=0.2, size=(11, 11)))
    model = raw.regularize(raw.fit_periodogram_base(samples, 11, "hann"), 0.1)
    amplitude = 0.037
    result, _ = raw.analyze_model(model, template / factors, amplitude * template / factors, samples)
    stage.require(np.isclose(result["amplitude_no_mean"], amplitude, rtol=2e-12, atol=2e-14),
                  "standardizing data and template changed physical amplitude")
    print("KLIP Stage-B local radial-normalization checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the local radial-normalization command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("raw", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the requested radial-normalization action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        run(args.bundle.resolve(), args.raw.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
