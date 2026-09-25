#!/usr/bin/env python3
"""Test patch-RMS PSD weighting after corrected KLIP radial standardization."""
from __future__ import annotations

import argparse
import json
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
import run_klip_stage_b_local_patch_rms as patch  # noqa: E402
import run_klip_stage_b_local_planet_masked_screen as masked  # noqa: E402
import run_klip_stage_b_local_radial_normalization as radial  # noqa: E402


SUPPORT = 11
METHODS = patch.METHODS


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def summarize(records: list[dict[str, object]], parent: dict[str, object],
              parent_records: list[dict[str, object]]) -> dict[str, object]:
    """Aggregate radially standardized patch-RMS scores and paired parent scores."""
    parent_lookup = {(row["radius"], row["site_index"], row["search_index"], row["method"]): row
                     for row in parent_records if row["method"] in {
                         method[0] for method in METHODS}}
    stage.require(len(parent_lookup) == len(records),
                  "radial parent does not cover every patch-RMS record")
    groups = {}
    for radius in raw.RADII:
        for method, _, _ in METHODS:
            chosen = [row for row in records if row["radius"] == radius and
                      row["method"] == method]
            stage.require(len(chosen) == 60, "radial patch-RMS summary lost query records")
            no_mean = np.asarray([value for row in chosen for value in row["score_no_mean"]])
            fitted_mean = np.asarray([value for row in chosen
                                      for value in row["score_fitted_mean"]])
            parent_no_mean = np.asarray([
                value for row in chosen for value in parent_lookup[(
                    row["radius"], row["site_index"], row["search_index"],
                    row["method"])]["score_no_mean"]])
            parent_fitted_mean = np.asarray([
                value for row in chosen for value in parent_lookup[(
                    row["radius"], row["site_index"], row["search_index"],
                    row["method"])]["score_fitted_mean"]])
            centers = np.asarray([value for row in chosen if row["search_index"] == 0
                                  for value in row["score_no_mean"]])
            maximums = []
            for site_index in range(12):
                site = sorted((row for row in chosen if row["site_index"] == site_index),
                              key=lambda row: row["search_index"])
                stage.require(len(site) == 5, "radial patch-RMS search is incomplete")
                for detector_half in (0, 1):
                    maximums.append(max(row["score_no_mean"][detector_half] for row in site))
            first = np.asarray([row["score_no_mean"][0] for row in chosen])
            second = np.asarray([row["score_no_mean"][1] for row in chosen])
            delta = no_mean - parent_no_mean
            entry = {"support": SUPPORT, "radius": radius, "method": method,
                     "band_half_width": chosen[0]["band_half_width"],
                     "directional_candidate_scores_no_mean": raw.score_stats(no_mean),
                     "directional_candidate_scores_fitted_mean": raw.score_stats(fitted_mean),
                     "directional_center_scores_no_mean": raw.score_stats(centers),
                     "five_pixel_maximum_no_mean": raw.score_stats(np.asarray(maximums)),
                     "parent_candidate_scores_no_mean": raw.score_stats(parent_no_mean),
                     "parent_candidate_scores_fitted_mean": raw.score_stats(parent_fitted_mean),
                     "median_split_physical_weight_cosine": float(np.median([
                         row["physical_weight_cosine"] for row in chosen])),
                     "split_score_correlation": float(np.corrcoef(first, second)[0, 1]),
                     "median_patch_rms_max_min_ratio": float(np.median([
                         maximum / minimum for row in chosen
                         for minimum, maximum in zip(row["patch_rms_min"],
                                                    row["patch_rms_max"])])),
                     "median_normalized_mean_rms": float(np.median([
                         value for row in chosen for value in row["normalized_mean_rms"]])),
                     "paired_parent_scores": {
                         "correlation": float(np.corrcoef(no_mean, parent_no_mean)[0, 1]),
                         "median_absolute_difference": float(np.median(np.abs(delta))),
                         "rms_difference": float(np.sqrt(np.mean(np.square(delta)))),
                         "maximum_absolute_difference": float(np.max(np.abs(delta)))},
                     "maximum_solver_relative_residual": float(max(
                         value for row in chosen for value in row["relative_residual"]))}
            for label in ("heldout_no_mean", "heldout_fitted_mean"):
                values = [item for row in chosen for item in row[label]]
                entry[label] = {
                    "directional_fits": len(values),
                    "median_variance_over_prediction": float(np.median([
                        item["variance"] for item in values])),
                    "variance_10_90": np.percentile([
                        item["variance"] for item in values], [10, 90]).tolist(),
                    "median_mse_over_prediction": float(np.median([
                        item["mse"] for item in values])),
                    "median_within_one_fraction": float(np.median([
                        item["within_one_fraction"] for item in values]))}
            groups[f"r{radius:g}_{method}"] = entry

    primary = []
    for method, _, _ in METHODS:
        selected = [groups[f"r{radius:g}_{method}"] for radius in raw.PRIMARY_RADII]
        primary_rows = [row for row in records if row["radius"] in raw.PRIMARY_RADII and
                        row["method"] == method]
        scores = np.asarray([value for row in primary_rows for value in row["score_no_mean"]])
        parent_scores = np.asarray([
            value for row in primary_rows for value in parent_lookup[(
                row["radius"], row["site_index"], row["search_index"],
                row["method"])]["score_no_mean"]])
        delta = scores - parent_scores
        primary.append({
            "method": method,
            "patch_rms_score_variance_median": float(np.median([
                row["directional_candidate_scores_no_mean"]["variance"] for row in selected])),
            "parent_score_variance_median": float(np.median([
                row["parent_candidate_scores_no_mean"]["variance"] for row in selected])),
            "patch_rms_over_parent_variance_median": float(np.median([
                row["directional_candidate_scores_no_mean"]["variance"] /
                row["parent_candidate_scores_no_mean"]["variance"] for row in selected])),
            "patch_rms_fitted_mean_variance_median": float(np.median([
                row["directional_candidate_scores_fitted_mean"]["variance"]
                for row in selected])),
            "parent_fitted_mean_variance_median": float(np.median([
                row["parent_candidate_scores_fitted_mean"]["variance"]
                for row in selected])),
            "patch_rms_heldout_variance_median": float(np.median([
                row["heldout_no_mean"]["median_variance_over_prediction"]
                for row in selected])),
            "patch_rms_split_physical_weight_cosine_median": float(np.median([
                row["median_split_physical_weight_cosine"] for row in selected])),
            "paired_parent_score_correlation": float(np.corrcoef(scores, parent_scores)[0, 1]),
            "paired_parent_score_median_absolute_difference": float(np.median(np.abs(delta))),
            "paired_parent_score_rms_difference": float(np.sqrt(np.mean(np.square(delta)))),
            "paired_parent_score_maximum_absolute_difference": float(np.max(np.abs(delta))),
            "median_patch_rms_max_min_ratio": float(np.median([
                row["median_patch_rms_max_min_ratio"] for row in selected]))})
    return {"groups": groups, "primary": primary}


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the radial-plus-patch-RMS comparison."""
    lines = ["# KLIP Stage-B radial plus post-mean patch-RMS control", "",
             "This test applies the post-ensemble-mean per-patch RMS control after the strict leave-site-out "
             "radial standardization. Candidate data and exact responses use the same radial scale and fixed "
             "planet mask as the corrected radial parent. The PSD trace remains on the standardized raw "
             "post-mean variance scale.", "", "## Controlling radii", "",
             "| Method | Radial parent variance | + Patch-RMS variance | Ratio | + Fitted mean | Opposite-half variance | Split physical-weight cosine | Paired score corr. | Median abs. Δscore |",
             "| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in result["primary"]:
        lines.append(f"| {row['method']} | {row['parent_score_variance_median']:.3f} | "
                     f"{row['patch_rms_score_variance_median']:.3f} | "
                     f"{row['patch_rms_over_parent_variance_median']:.3f} | "
                     f"{row['patch_rms_fitted_mean_variance_median']:.3f} | "
                     f"{row['patch_rms_heldout_variance_median']:.3f} | "
                     f"{row['patch_rms_split_physical_weight_cosine_median']:.3f} | "
                     f"{row['paired_parent_score_correlation']:.6f} | "
                     f"{row['paired_parent_score_median_absolute_difference']:.3f} |")
    lines.extend(["", "Every comparison is paired on the same frozen query, candidate mask, radial profile, "
                  "training band, and detector-half split."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for axis, (method, _, _) in zip(axes, METHODS):
        parent_values = []
        patch_values = []
        for radius in raw.RADII:
            row = result["groups"][f"r{radius:g}_{method}"]
            parent_values.append(row["parent_candidate_scores_no_mean"]["variance"])
            patch_values.append(row["directional_candidate_scores_no_mean"]["variance"])
        x = np.arange(len(raw.RADII)); width = 0.36
        axis.bar(x - width / 2, parent_values, width, label="radial parent")
        axis.bar(x + width / 2, patch_values, width, label="radial + patch RMS")
        axis.axhline(1, color="black", linestyle="--", linewidth=1)
        axis.set(title=method, ylabel="candidate score variance", xticks=x,
                 xticklabels=[f"{radius:g}" for radius in raw.RADII], xlabel="radius (pixels)")
        axis.grid(axis="y", alpha=0.2)
        axis.legend()
    fig.suptitle("KLIP mode 200, corrected 11-pixel radial coordinates")
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, parent_root: Path, root: Path) -> None:
    """Run patch-RMS weighting in the verified radial-standardized coordinates."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = raw.verify_bundle(bundle_path)
    completion = read(parent_root / "complete.json")
    stage.verify([completion[key] for key in ("protocol", "records", "results", "report", "figure")])
    parent_result = read(parent_root / "results.json")
    parent_records = read(parent_root / "records.json")
    root.mkdir(parents=True)
    protocol = {"schema": 1,
                "purpose": "patch-RMS control after corrected strict radial standardization",
                "mode": 200, "support": SUPPORT, "radii": list(raw.RADII),
                "primary_radii": list(raw.PRIMARY_RADII),
                "methods": [{"name": name, "window": window, "mixing": mixing}
                            for name, window, mixing in METHODS],
                "radial_profile": "parent leave-site-out 3.6-pixel profile",
                "normalization": "in radial-standardized coordinates, subtract ensemble pixelwise mean; divide every residual patch by its own RMS; do not recenter; preserve standardized post-mean variance",
                "candidate_mask": "fixed seven-pixel planet disk removed from standardized candidate and exact response through the covariance principal-submatrix solve",
                "bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
                "parent": stage.fingerprint(parent_root / "complete.json"),
                "scripts": [stage.fingerprint(Path(__file__)),
                            stage.fingerprint(Path(__file__).with_name(
                                "run_klip_stage_b_local_patch_rms.py")),
                            stage.fingerprint(Path(__file__).with_name(
                                "run_klip_stage_b_local_planet_masked_normalization.py"))]}
    stage.write_json(root / "protocol.json", protocol)
    stage.write_json(root / "state.json", {"status": "running", "started_unix": time.time()})

    baseline = np.asarray(arrays["baseline"], dtype=np.float64)
    finite = np.isfinite(baseline)
    planet = raw.planet_position(baseline.shape, metadata)
    lookup = {(int(row), int(column)): index
              for index, (row, column) in enumerate(arrays["positions"])}
    records = []
    started = time.monotonic()
    for radius in raw.RADII:
        fixed = metadata["geometry"]["11"][str(radius)]
        width = int(fixed["training"]["narrowest_split_supported_half_width"])
        for site_index, site in enumerate(fixed["selected_sites"]):
            searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                        for delta_row, delta_column in footprint.SEARCH_OFFSETS]
            excluded = radial.exclusion_mask(baseline.shape, searches, SUPPORT, planet)
            scale, _ = radial.variance_profile(baseline, excluded)
            standardized = baseline / scale
            for search_index, query in enumerate(searches):
                rings = raw.training_stencils(finite, query, searches, planet, SUPPORT)
                source = lookup[query]
                template = raw.crop(np.asarray(arrays["response"][source], dtype=np.float64).T,
                                    SUPPORT)
                validity = raw.crop(np.asarray(arrays["validity"][source], dtype=bool).T, SUPPORT)
                planet_mask = masked.candidate_mask(query, SUPPORT, planet) & validity
                row, column = query
                candidate = baseline[column - 5:column + 6, row - 5:row + 6].T
                candidate_scale = scale[column - 5:column + 6, row - 5:row + 6].T
                standardized_template = template / candidate_scale
                standardized_candidate = candidate / candidate_scale
                halves = []
                for detector_half in (0, 1):
                    samples, interpolation_gain = raw.band_samples(
                        standardized, rings, width, detector_half)
                    opposite, _ = raw.band_samples(standardized, rings, width, 1 - detector_half)
                    results = {}
                    details = {}
                    for method, window, mixing in METHODS:
                        base = patch.fit_patch_rms_base(samples, window)
                        reference = raw.fit_periodogram_base(samples, SUPPORT, window)
                        stage.require(np.array_equal(base["mean"], reference["mean"]) and
                                      np.isclose(base["target_variance"],
                                                 reference["target_variance"],
                                                 rtol=1e-14, atol=0),
                                      "radial patch RMS changed mean or variance")
                        model = raw.regularize(base, mixing)
                        result, detail = masked.analyze_model(
                            model, standardized_template, standardized_candidate,
                            planet_mask, opposite)
                        result.update({"patch_rms_min": float(np.min(base["patch_rms"])),
                                       "patch_rms_median": float(np.median(base["patch_rms"])),
                                       "patch_rms_max": float(np.max(base["patch_rms"])),
                                       "normalized_mean_rms": base["normalized_mean_rms"]})
                        results[method] = result
                        details[method] = detail | {
                            "physical_weight": detail["weight"] / candidate_scale.ravel()}
                    halves.append({"samples": len(samples),
                                   "interpolation_variance_gain": interpolation_gain,
                                   "results": results, "details": details})
                for method, _, _ in METHODS:
                    first = halves[0]["details"][method]
                    second = halves[1]["details"][method]
                    cosine = float(first["physical_weight"] @ second["physical_weight"] /
                                   (np.linalg.norm(first["physical_weight"]) *
                                    np.linalg.norm(second["physical_weight"])))
                    records.append({
                        "support": SUPPORT, "radius": radius, "site_index": site_index,
                        "search_index": search_index, "row": row, "column": column,
                        "band_half_width": width, "method": method,
                        "samples": [halves[index]["samples"] for index in (0, 1)],
                        "physical_weight_cosine": cosine,
                        "valid_pixels": int(np.count_nonzero(planet_mask)),
                        "score_no_mean": [halves[index]["results"][method]["score_no_mean"]
                                          for index in (0, 1)],
                        "score_fitted_mean": [halves[index]["results"][method][
                            "score_fitted_mean"] for index in (0, 1)],
                        "heldout_no_mean": [halves[index]["results"][method][
                            "heldout_no_mean"] for index in (0, 1)],
                        "heldout_fitted_mean": [halves[index]["results"][method][
                            "heldout_fitted_mean"] for index in (0, 1)],
                        "patch_rms_min": [halves[index]["results"][method]["patch_rms_min"]
                                          for index in (0, 1)],
                        "patch_rms_median": [halves[index]["results"][method][
                            "patch_rms_median"] for index in (0, 1)],
                        "patch_rms_max": [halves[index]["results"][method]["patch_rms_max"]
                                          for index in (0, 1)],
                        "normalized_mean_rms": [halves[index]["results"][method][
                            "normalized_mean_rms"] for index in (0, 1)],
                        "iterations": [halves[index]["results"][method]["iterations"]
                                       for index in (0, 1)],
                        "relative_residual": [halves[index]["results"][method][
                            "relative_residual"] for index in (0, 1)]})
        print(f"completed radius {radius:g}: {len(records)} records", flush=True)

    summarized = summarize(records, parent_result, parent_records)
    result = {"purpose": protocol["purpose"], "mode": 200, "support": SUPPORT,
              "records": len(records), "elapsed_seconds": time.monotonic() - started,
              "primary": summarized["primary"], "groups": summarized["groups"]}
    stage.write_json(root / "records.json", records)
    stage.write_json(root / "results.json", result)
    write_report(root, result)
    completion = {"status": "complete", "protocol": stage.fingerprint(root / "protocol.json"),
                  "records": stage.fingerprint(root / "records.json"),
                  "results": stage.fingerprint(root / "results.json"),
                  "report": stage.fingerprint(root / "results.md"),
                  "figure": stage.fingerprint(root / "comparison.png")}
    stage.write_json(root / "complete.json", completion)
    stage.write_json(root / "state.json", {"status": "complete",
                                             "elapsed_seconds": result["elapsed_seconds"]})
    print(root / "results.md", flush=True)


def check() -> None:
    """Check simultaneous radial scaling, patch RMS, and masking preserve amplitude."""
    generator = np.random.default_rng(250926)
    samples = generator.normal(size=(41, 121))
    samples[:6] *= 3
    scale = np.exp(generator.normal(scale=0.25, size=(11, 11)))
    template = generator.normal(size=(11, 11))
    mask = generator.random((11, 11)) > 0.2
    amplitude = 0.029
    for _, window, mixing in METHODS:
        model = raw.regularize(patch.fit_patch_rms_base(samples, window), mixing)
        result, _ = masked.analyze_model(model, template / scale,
                                         amplitude * template / scale, mask, samples)
        stage.require(np.isclose(result["amplitude_no_mean"], amplitude,
                                 rtol=2e-12, atol=2e-14),
                      "radial patch-RMS solve changed physical amplitude")
    print("KLIP Stage-B radial plus patch-RMS checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the radial patch-RMS command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("parent", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the radial patch-RMS action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        run(args.bundle.resolve(), args.parent.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
