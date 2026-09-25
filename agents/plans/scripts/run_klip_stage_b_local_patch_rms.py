#!/usr/bin/env python3
"""Test post-ensemble-mean patch-RMS PSD weighting on corrected KLIP geometry."""
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
import run_klip_stage_b_local_noise_screen as raw  # noqa: E402
import run_klip_stage_b_local_planet_masked_screen as masked  # noqa: E402


SUPPORT = 11
METHODS = (("rectangular_m0.3", "rectangular", 0.3),
           ("hann_m0.1", "hann", 0.1))
BANDS = ("narrow", "full")


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def fit_patch_rms_base(samples: np.ndarray, window_name: str) -> dict[str, object]:
    """Fit one unmixed PSD after equalizing post-mean training-patch RMS."""
    samples = np.asarray(samples, dtype=np.float64)
    stage.require(samples.ndim == 2 and samples.shape[1] == 121 and len(samples) >= 8 and
                  np.all(np.isfinite(samples)), "invalid patch-RMS training samples")
    mean = np.mean(samples, axis=0)
    centered = samples - mean
    rms = np.sqrt(np.mean(np.square(centered), axis=1))
    roundoff = 16 * np.finfo(np.float64).eps * float(np.max(np.abs(samples)))
    stage.require(np.all(np.isfinite(rms)) and np.all(rms > roundoff),
                  "patch-RMS training contains a degenerate residual")
    normalized = centered / rms[:, None]
    stage.require(np.allclose(np.mean(np.square(normalized), axis=1), 1,
                              rtol=1e-12, atol=1e-14),
                  "post-mean patch-RMS normalization failed")

    target = float(np.sum(np.square(centered)) / ((len(samples) - 1) * 121))
    stage.require(np.isfinite(target) and target > 0, "invalid physical training variance")
    window = extension.estimation_window(window_name)
    window_energy = float(np.sum(np.square(window)))
    transformed = np.fft.fft2(normalized.reshape((-1, 11, 11)) * window,
                              s=(21, 21), axes=(-2, -1))
    raw_power = np.sum(np.square(np.abs(transformed)), axis=0) / (
        (len(samples) - 1) * window_energy)
    raw_zero_lag = float(np.mean(raw_power))
    stage.require(np.isfinite(raw_zero_lag) and raw_zero_lag > 0,
                  "invalid patch-RMS periodogram")
    return {"mean": mean, "target_variance": target, "window": window_name,
            "window_energy": window_energy, "mixing": None, "fft_size": 21,
            "raw_windowed_zero_lag": raw_zero_lag,
            "psd_rescaling": target / raw_zero_lag,
            "rescaled_power": raw_power * (target / raw_zero_lag),
            "patch_rms": rms,
            "normalized_mean_rms": float(np.sqrt(np.mean(np.square(np.mean(normalized, axis=0)))))}


def score_stats(values: np.ndarray) -> dict[str, object]:
    """Delegate finite standardized-score summaries to the raw-screen contract."""
    return raw.score_stats(np.asarray(values, dtype=np.float64))


def summarize(records: list[dict[str, object]], parent: dict[str, object],
              parent_records: list[dict[str, object]]) -> dict[str, object]:
    """Aggregate patch-RMS candidate calibration and paired raw-parent comparisons."""
    parent_groups = parent["groups"]
    parent_lookup = {
        (row["radius"], row["site_index"], row["search_index"], row["band_role"], row["method"]): row
        for row in parent_records if row["support"] == SUPPORT and
        row["method"] in {method[0] for method in METHODS}}
    stage.require(len(parent_lookup) == len(records),
                  "paired raw parent does not cover every patch-RMS record")
    groups = {}
    for radius in raw.RADII:
        for band_role in BANDS:
            for method, _, _ in METHODS:
                chosen = [row for row in records if row["radius"] == radius and
                          row["band_role"] == band_role and row["method"] == method]
                stage.require(len(chosen) == 60, "patch-RMS summary lost fixed query records")
                no_mean = np.asarray([value for row in chosen for value in row["score_no_mean"]])
                fitted_mean = np.asarray([value for row in chosen
                                          for value in row["score_fitted_mean"]])
                centers = np.asarray([value for row in chosen if row["search_index"] == 0
                                      for value in row["score_no_mean"]])
                maximums = []
                for site_index in range(12):
                    site = sorted((row for row in chosen if row["site_index"] == site_index),
                                  key=lambda row: row["search_index"])
                    stage.require(len(site) == 5, "patch-RMS five-pixel search is incomplete")
                    for detector_half in (0, 1):
                        maximums.append(max(row["score_no_mean"][detector_half] for row in site))
                first = np.asarray([row["score_no_mean"][0] for row in chosen])
                second = np.asarray([row["score_no_mean"][1] for row in chosen])
                paired_raw = np.asarray([
                    value for row in chosen for value in parent_lookup[(
                        row["radius"], row["site_index"], row["search_index"],
                        row["band_role"], row["method"])]["score_no_mean"]])
                delta = no_mean - paired_raw
                raw_group = parent_groups[f"s11_r{radius:g}_{band_role}_{method}"]
                entry = {"support": SUPPORT, "radius": radius, "band_role": band_role,
                         "band_half_width": chosen[0]["band_half_width"], "method": method,
                         "directional_candidate_scores_no_mean": score_stats(no_mean),
                         "directional_candidate_scores_fitted_mean": score_stats(fitted_mean),
                         "directional_center_scores_no_mean": score_stats(centers),
                         "five_pixel_maximum_no_mean": score_stats(np.asarray(maximums)),
                         "median_split_weight_cosine": float(np.median([
                             row["weight_cosine"] for row in chosen])),
                         "split_weight_cosine_10_90": np.percentile([
                             row["weight_cosine"] for row in chosen], [10, 90]).tolist(),
                         "split_score_correlation": float(np.corrcoef(first, second)[0, 1]),
                         "median_absolute_split_score_difference": float(np.median(
                             np.abs(first - second))),
                         "training_patch_count_range": [
                             min(value for row in chosen for value in row["samples"]),
                             max(value for row in chosen for value in row["samples"])],
                         "median_patch_rms_max_min_ratio": float(np.median([
                             maximum / minimum for row in chosen
                             for minimum, maximum in zip(row["patch_rms_min"],
                                                        row["patch_rms_max"])])),
                         "median_normalized_mean_rms": float(np.median([
                             value for row in chosen for value in row["normalized_mean_rms"]])),
                         "maximum_solver_relative_residual": float(max(
                             value for row in chosen for value in row["relative_residual"])),
                         "maximum_solver_iterations": int(max(
                             value for row in chosen for value in row["iterations"])),
                         "paired_raw_scores": {
                             "correlation": float(np.corrcoef(no_mean, paired_raw)[0, 1]),
                             "median_absolute_difference": float(np.median(np.abs(delta))),
                             "rms_difference": float(np.sqrt(np.mean(np.square(delta)))),
                             "maximum_absolute_difference": float(np.max(np.abs(delta)))},
                         "raw_comparison": {
                             "candidate_variance_no_mean": raw_group[
                                 "directional_candidate_scores_no_mean"]["variance"],
                             "candidate_variance_fitted_mean": raw_group[
                                 "directional_candidate_scores_fitted_mean"]["variance"],
                             "heldout_variance_no_mean": raw_group["heldout_no_mean"][
                                 "median_variance_over_prediction"],
                             "heldout_variance_fitted_mean": raw_group["heldout_fitted_mean"][
                                 "median_variance_over_prediction"],
                             "split_weight_cosine": raw_group["median_split_weight_cosine"]}}
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
                            item["within_two_fraction"] for item in projections]))}
                groups[f"r{radius:g}_{band_role}_{method}"] = entry

    primary = []
    for band_role in BANDS:
        for method, _, _ in METHODS:
            selected = [groups[f"r{radius:g}_{band_role}_{method}"]
                        for radius in raw.PRIMARY_RADII]
            primary_records = [row for row in records if row["radius"] in raw.PRIMARY_RADII and
                               row["band_role"] == band_role and row["method"] == method]
            patch_scores = np.asarray([value for row in primary_records
                                       for value in row["score_no_mean"]])
            raw_scores = np.asarray([
                value for row in primary_records for value in parent_lookup[(
                    row["radius"], row["site_index"], row["search_index"],
                    row["band_role"], row["method"])]["score_no_mean"]])
            score_delta = patch_scores - raw_scores
            primary.append({
                "band_role": band_role, "method": method,
                "patch_rms_score_variance_median": float(np.median([
                    row["directional_candidate_scores_no_mean"]["variance"]
                    for row in selected])),
                "raw_score_variance_median": float(np.median([
                    row["raw_comparison"]["candidate_variance_no_mean"] for row in selected])),
                "patch_rms_over_raw_variance_median": float(np.median([
                    row["directional_candidate_scores_no_mean"]["variance"] /
                    row["raw_comparison"]["candidate_variance_no_mean"] for row in selected])),
                "patch_rms_fitted_mean_variance_median": float(np.median([
                    row["directional_candidate_scores_fitted_mean"]["variance"]
                    for row in selected])),
                "raw_fitted_mean_variance_median": float(np.median([
                    row["raw_comparison"]["candidate_variance_fitted_mean"] for row in selected])),
                "patch_rms_heldout_variance_median": float(np.median([
                    row["heldout_no_mean"]["median_variance_over_prediction"]
                    for row in selected])),
                "raw_heldout_variance_median": float(np.median([
                    row["raw_comparison"]["heldout_variance_no_mean"] for row in selected])),
                "patch_rms_split_weight_cosine_median": float(np.median([
                    row["median_split_weight_cosine"] for row in selected])),
                "paired_raw_score_correlation": float(np.corrcoef(
                    patch_scores, raw_scores)[0, 1]),
                "paired_raw_score_median_absolute_difference": float(np.median(
                    np.abs(score_delta))),
                "paired_raw_score_rms_difference": float(np.sqrt(np.mean(
                    np.square(score_delta)))),
                "paired_raw_score_maximum_absolute_difference": float(np.max(
                    np.abs(score_delta))),
                "median_patch_rms_max_min_ratio": float(np.median([
                    row["median_patch_rms_max_min_ratio"] for row in selected]))})
    return {"groups": groups, "primary": primary}


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the patch-RMS comparison and diagnostic figure."""
    lines = ["# KLIP Stage-B post-mean patch-RMS control", "",
             "This mode-200 control uses the corrected 11-pixel planet-masked geometry. Each training patch "
             "first has the ensemble pixelwise mean removed and is then divided by its own RMS. The normalized "
             "patches enter the periodogram without a second ensemble centering. Candidate data and exact "
             "responses remain in physical contrast coordinates, and the spectrum is rescaled to the raw "
             "post-mean training variance.", "", "## Controlling radii", "",
             "| Band | Method | Raw variance | Patch-RMS variance | Patch/raw | Patch-RMS + fitted mean | Opposite-half variance | Split weight cosine | Paired score corr. | Median abs. Δscore |",
             "| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in result["primary"]:
        lines.append(f"| {row['band_role']} | {row['method']} | "
                     f"{row['raw_score_variance_median']:.3f} | "
                     f"{row['patch_rms_score_variance_median']:.3f} | "
                     f"{row['patch_rms_over_raw_variance_median']:.3f} | "
                     f"{row['patch_rms_fitted_mean_variance_median']:.3f} | "
                     f"{row['patch_rms_heldout_variance_median']:.3f} | "
                     f"{row['patch_rms_split_weight_cosine_median']:.3f} | "
                     f"{row['paired_raw_score_correlation']:.6f} | "
                     f"{row['paired_raw_score_median_absolute_difference']:.3f} |")
    spread = {row["band_role"]: row["median_patch_rms_max_min_ratio"]
              for row in result["primary"] if row["method"] == "hann_m0.1"}
    lines.extend(["", "The raw comparison is the verified planet-masked parent on identical frozen sites, "
                  "bands, candidate masks, and detector-half splits. Fitted mean changes only candidate "
                  "evaluation; patch RMS always follows training ensemble-mean subtraction. Median maximum-to-"
                  f"minimum patch-RMS ratios are {spread['narrow']:.2f} for the narrow bands and "
                  f"{spread['full']:.2f} for the full bands."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    for axis, (method, _, _) in zip(axes, METHODS):
        raw_values = []
        patch_values = []
        for radius in raw.RADII:
            row = result["groups"][f"r{radius:g}_narrow_{method}"]
            raw_values.append(row["raw_comparison"]["candidate_variance_no_mean"])
            patch_values.append(row["directional_candidate_scores_no_mean"]["variance"])
        x = np.arange(len(raw.RADII)); width = 0.36
        axis.bar(x - width / 2, raw_values, width, label="raw")
        axis.bar(x + width / 2, patch_values, width, label="post-mean patch RMS")
        axis.axhline(1, color="black", linestyle="--", linewidth=1)
        axis.set(title=method, ylabel="candidate score variance", xticks=x,
                 xticklabels=[f"{radius:g}" for radius in raw.RADII], xlabel="radius (pixels)")
        axis.grid(axis="y", alpha=0.2)
        axis.legend()
    fig.suptitle("KLIP mode 200, corrected 11-pixel planet-masked narrow bands")
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, parent_root: Path, root: Path) -> None:
    """Run the patch-RMS control from verified local products."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = raw.verify_bundle(bundle_path)
    parent_completion = read(parent_root / "complete.json")
    stage.verify([parent_completion[key] for key in ("protocol", "records", "results", "report", "figure")])
    parent_result = read(parent_root / "results.json")
    parent_records = read(parent_root / "records.json")
    root.mkdir(parents=True)
    protocol = {"schema": 1,
                "purpose": "post-ensemble-mean patch-RMS control on corrected KLIP geometry",
                "mode": 200, "support": SUPPORT, "radii": list(raw.RADII),
                "primary_radii": list(raw.PRIMARY_RADII),
                "methods": [{"name": name, "window": window, "mixing": mixing}
                            for name, window, mixing in METHODS],
                "bands": list(BANDS),
                "normalization": "subtract raw ensemble pixelwise mean; divide every residual patch by its own 121-pixel RMS; do not recenter normalized patches; rescale PSD trace to raw post-mean variance",
                "candidate_coordinates": "raw physical contrast; fixed seven-pixel planet disk removed from candidate and exact response through the covariance principal-submatrix solve",
                "bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
                "parent": stage.fingerprint(parent_root / "complete.json"),
                "scripts": [stage.fingerprint(Path(__file__)),
                            stage.fingerprint(Path(__file__).with_name(
                                "run_klip_stage_b_local_planet_masked_screen.py")),
                            stage.fingerprint(Path(__file__).with_name(
                                "run_klip_stage_b_local_noise_screen.py"))]}
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
        narrow = int(fixed["training"]["narrowest_split_supported_half_width"])
        for site_index, site in enumerate(fixed["selected_sites"]):
            searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                        for delta_row, delta_column in footprint.SEARCH_OFFSETS]
            for search_index, query in enumerate(searches):
                rings = raw.training_stencils(finite, query, searches, planet, SUPPORT)
                source = lookup[query]
                template = raw.crop(np.asarray(arrays["response"][source], dtype=np.float64).T,
                                    SUPPORT)
                validity = raw.crop(np.asarray(arrays["validity"][source], dtype=bool).T, SUPPORT)
                candidate_mask = masked.candidate_mask(query, SUPPORT, planet) & validity
                row, column = query
                candidate = baseline[column - 5:column + 6, row - 5:row + 6].T
                stage.require(np.all(np.isfinite(candidate[candidate_mask])) and
                              candidate_mask[5, 5], "patch-RMS candidate lost its source anchor")
                for band_role, width in (("narrow", narrow), ("full", 60)):
                    halves = []
                    for detector_half in (0, 1):
                        samples, interpolation_gain = raw.band_samples(
                            baseline, rings, width, detector_half)
                        opposite, _ = raw.band_samples(baseline, rings, width, 1 - detector_half)
                        results = {}
                        details = {}
                        for method, window, mixing in METHODS:
                            base = fit_patch_rms_base(samples, window)
                            raw_base = raw.fit_periodogram_base(samples, SUPPORT, window)
                            stage.require(np.array_equal(base["mean"], raw_base["mean"]) and
                                          np.isclose(base["target_variance"],
                                                     raw_base["target_variance"],
                                                     rtol=1e-14, atol=0),
                                          "patch RMS changed raw mean or physical variance")
                            model = raw.regularize(base, mixing)
                            result, detail = masked.analyze_model(
                                model, template, candidate, candidate_mask, opposite)
                            result.update({"patch_rms_min": float(np.min(base["patch_rms"])),
                                           "patch_rms_median": float(np.median(base["patch_rms"])),
                                           "patch_rms_max": float(np.max(base["patch_rms"])),
                                           "normalized_mean_rms": base["normalized_mean_rms"]})
                            results[method] = result
                            details[method] = detail
                        halves.append({"samples": len(samples),
                                       "interpolation_variance_gain": interpolation_gain,
                                       "results": results, "details": details})
                    for method, _, _ in METHODS:
                        first = halves[0]["details"][method]
                        second = halves[1]["details"][method]
                        cosine = float(first["weight"] @ second["weight"] /
                                       (np.linalg.norm(first["weight"]) *
                                        np.linalg.norm(second["weight"])))
                        records.append({
                            "support": SUPPORT, "radius": radius, "site_index": site_index,
                            "search_index": search_index, "row": row, "column": column,
                            "band_role": band_role, "band_half_width": width, "method": method,
                            "samples": [halves[index]["samples"] for index in (0, 1)],
                            "interpolation_variance_gain": [halves[index][
                                "interpolation_variance_gain"] for index in (0, 1)],
                            "weight_cosine": cosine,
                            "valid_pixels": int(np.count_nonzero(candidate_mask)),
                            "score_no_mean": [halves[index]["results"][method]["score_no_mean"]
                                              for index in (0, 1)],
                            "score_fitted_mean": [halves[index]["results"][method][
                                "score_fitted_mean"] for index in (0, 1)],
                            "sigma": [halves[index]["results"][method]["sigma"]
                                      for index in (0, 1)],
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
    """Check the patch-RMS coordinate definition and masked amplitude response."""
    generator = np.random.default_rng(250925)
    samples = generator.normal(size=(37, 121))
    samples[:5] *= 4
    for window, mixing in (("rectangular", 0.3), ("hann", 0.1)):
        base = fit_patch_rms_base(samples, window)
        reference = raw.fit_periodogram_base(samples, SUPPORT, window)
        stage.require(np.array_equal(base["mean"], reference["mean"]) and
                      np.isclose(base["target_variance"], reference["target_variance"],
                                 rtol=1e-14, atol=0) and
                      np.max(base["patch_rms"]) / np.min(base["patch_rms"]) > 3,
                      "patch-RMS check lost the raw physical-coordinate contract")
        model = raw.regularize(base, mixing)
        template = generator.normal(size=(SUPPORT, SUPPORT))
        mask = generator.random((SUPPORT, SUPPORT)) > 0.2
        amplitude = 0.037
        measured, _ = masked.analyze_model(model, template, amplitude * template,
                                           mask, samples)
        stage.require(np.isclose(measured["amplitude_no_mean"], amplitude,
                                 rtol=2e-12, atol=2e-14),
                      "patch-RMS masked solve changed physical amplitude")
    print("KLIP Stage-B post-mean patch-RMS checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the patch-RMS command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("parent", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the patch-RMS action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        run(args.bundle.resolve(), args.parent.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
