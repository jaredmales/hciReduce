#!/usr/bin/env python3
"""Test a support-independent radial mean with the corrected KLIP PSD filters."""
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
import run_klip_stage_b_local_planet_masked_screen as masked  # noqa: E402


METHODS = ("identity", "rectangular_m0.3", "hann_m0.1")
MAXIMUM_PROFILE_RADIUS = 60


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def radial_mean_map(image: np.ndarray, planet: tuple[float, float]) -> tuple[np.ndarray, list[dict[str, object]]]:
    """Fit and interpolate a one-pixel annular mean after excluding the known planet."""
    native_column, native_row = np.indices(image.shape, dtype=np.float64)
    center_row = 0.5 * (image.shape[1] - 1)
    center_column = 0.5 * (image.shape[0] - 1)
    radius = np.hypot(native_row - center_row, native_column - center_column)
    planet_distance = np.hypot(native_row - planet[0], native_column - planet[1])
    usable = np.isfinite(image) & (planet_distance > stage.PLANET_EXCLUSION_RADIUS + 0.5)
    centers, means, profile = [], [], []
    for lower in range(MAXIMUM_PROFILE_RADIUS):
        selected = usable & (radius >= lower) & (radius <= lower + 1)
        values = np.asarray(image[selected], dtype=np.float64)
        stage.require(len(values) >= 2 and np.all(np.isfinite(values)),
                      f"radial mean bin {lower} lacks finite support")
        mean = float(np.mean(values))
        centers.append(lower + 0.5)
        means.append(mean)
        profile.append({"lower_radius": lower, "center_radius": lower + 0.5,
                        "samples": len(values), "mean": mean,
                        "standard_deviation": float(np.std(values, ddof=1))})
    mean_map = np.interp(radius, np.asarray(centers), np.asarray(means))
    stage.require(np.all(np.isfinite(mean_map)), "radial mean interpolation is not finite")
    return mean_map, profile


def summarize_values(values: np.ndarray) -> dict[str, object]:
    """Summarize one finite score vector."""
    return raw.score_stats(np.asarray(values, dtype=np.float64))


def five_pixel_maximum(chosen: list[dict[str, object]], field: str) -> dict[str, object]:
    """Summarize the maximum score across each fixed five-pixel search."""
    maximums = []
    for site in range(12):
        rows = sorted((row for row in chosen if row["site_index"] == site),
                      key=lambda row: row["search_index"])
        stage.require(len(rows) == 5, "five-pixel search is incomplete")
        for detector_half in (0, 1):
            maximums.append(max(row[field][detector_half] for row in rows))
    return summarize_values(np.asarray(maximums))


def summarize(records: list[dict[str, object]]) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Aggregate raw and radial-mean scores by fixed Stage-B policy."""
    groups = {}
    for support in raw.SUPPORTS:
        for radius in raw.RADII:
            for band_role in ("narrow", "full"):
                for method in METHODS:
                    chosen = [row for row in records if row["support"] == support and
                              row["radius"] == radius and row["band_role"] == band_role and
                              row["method"] == method]
                    stage.require(len(chosen) == 60, "radial-mean summary lost fixed query records")
                    raw_scores = np.asarray([value for row in chosen for value in row["score_raw"]])
                    mean_scores = np.asarray([value for row in chosen for value in row["score_radial_mean"]])
                    projection = raw_scores - mean_scores
                    first = np.asarray([row["score_radial_mean"][0] for row in chosen])
                    second = np.asarray([row["score_radial_mean"][1] for row in chosen])
                    key = f"s{support}_r{radius:g}_{band_role}_{method}"
                    entry = {"support": support, "radius": radius, "band_role": band_role,
                             "band_half_width": chosen[0]["band_half_width"], "method": method,
                             "directional_candidate_scores_raw": summarize_values(raw_scores),
                             "directional_candidate_scores_radial_mean": summarize_values(mean_scores),
                             "directional_center_scores_radial_mean": summarize_values(np.asarray([
                                 value for row in chosen if row["search_index"] == 0
                                 for value in row["score_radial_mean"]])),
                             "five_pixel_maximum_raw": five_pixel_maximum(chosen, "score_raw"),
                             "five_pixel_maximum_radial_mean": five_pixel_maximum(
                                 chosen, "score_radial_mean"),
                             "radial_mean_projection_in_sigma": summarize_values(projection),
                             "raw_radial_score_correlation": float(np.corrcoef(raw_scores, mean_scores)[0, 1]),
                             "median_absolute_score_change": float(np.median(np.abs(projection))),
                             "median_split_weight_cosine": float(np.median([
                                 row["weight_cosine"] for row in chosen])),
                             "split_score_correlation_radial_mean": float(np.corrcoef(first, second)[0, 1]),
                             "median_absolute_split_score_difference_radial_mean": float(
                                 np.median(np.abs(first - second)))}
                    if support == 11:
                        heldout = [item for row in chosen for item in row["heldout_radial_mean"]]
                        entry["heldout_radial_mean"] = {
                            "directional_fits": len(heldout),
                            "median_variance_over_prediction": float(np.median([
                                item["variance"] for item in heldout])),
                            "variance_10_90": np.percentile([
                                item["variance"] for item in heldout], [10, 90]).tolist(),
                            "median_mse_over_prediction": float(np.median([
                                item["mse"] for item in heldout]))}
                    groups[key] = entry

    primary = []
    for support in raw.SUPPORTS:
        for band_role in ("narrow", "full"):
            for method in METHODS:
                selected = [groups[f"s{support}_r{radius:g}_{band_role}_{method}"]
                            for radius in raw.PRIMARY_RADII]
                primary.append({
                    "support": support, "band_role": band_role, "method": method,
                    "score_variance_raw_median": float(np.median([
                        row["directional_candidate_scores_raw"]["variance"] for row in selected])),
                    "score_variance_radial_mean_median": float(np.median([
                        row["directional_candidate_scores_radial_mean"]["variance"] for row in selected])),
                    "absolute_score_mean_raw_median": float(np.median([
                        abs(row["directional_candidate_scores_raw"]["mean"]) for row in selected])),
                    "absolute_score_mean_radial_mean_median": float(np.median([
                        abs(row["directional_candidate_scores_radial_mean"]["mean"]) for row in selected])),
                    "five_pixel_maximum_raw_median": float(np.median([
                        row["five_pixel_maximum_raw"]["quantiles"][2] for row in selected])),
                    "five_pixel_maximum_radial_mean_median": float(np.median([
                        row["five_pixel_maximum_radial_mean"]["quantiles"][2] for row in selected])),
                    "median_absolute_score_change": float(np.median([
                        row["median_absolute_score_change"] for row in selected])),
                    "raw_radial_score_correlation_median": float(np.median([
                        row["raw_radial_score_correlation"] for row in selected])),
                    "heldout_patch_variance_radial_mean_median": (float(np.median([
                        row["heldout_radial_mean"]["median_variance_over_prediction"]
                        for row in selected])) if support == 11 else None)})
    return groups, primary


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the support-independent radial-mean comparison and figure."""
    lines = ["# KLIP Stage-B support-independent radial-mean control", "",
             "A single one-pixel annular mean profile is fitted on the signal-free mode-200 image. Native pixels "
             "within 7.5 pixels of the known planet are excluded, matching the `R + 0.5` convention used by "
             "`hciAnalyze`. Linear interpolation evaluates that fixed profile at every retained native coordinate "
             "in the 11-, 31-, and 47-pixel candidate stamps.", "",
             "Only the candidate mean changes. The corrected fixed planet mask, exact response, raw 11-pixel "
             "Welch PSD estimate, covariance regularization, unit-response weights, and conditional uncertainty "
             "are identical to the corrected raw parent. Parent raw scores replay exactly.", "",
             "## Controlling radii: 7.5, 10, and 12 pixels", "",
             "| Response | Band | Method | Raw variance | Radial-mean variance | Ratio | Raw abs(mean) | Radial abs(mean) | Raw 5-pixel max | Radial 5-pixel max | Median abs(score change) | Raw/radial corr. | 11px opposite-half variance |",
             "| ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in result["primary"]:
        heldout = ("—" if row["heldout_patch_variance_radial_mean_median"] is None else
                   f"{row['heldout_patch_variance_radial_mean_median']:.3f}")
        ratio = row["score_variance_radial_mean_median"] / row["score_variance_raw_median"]
        lines.append(f"| {row['support']} | {row['band_role']} | {row['method']} | "
                     f"{row['score_variance_raw_median']:.3f} | "
                     f"{row['score_variance_radial_mean_median']:.3f} | {ratio:.3f} | "
                     f"{row['absolute_score_mean_raw_median']:.3f} | "
                     f"{row['absolute_score_mean_radial_mean_median']:.3f} | "
                     f"{row['five_pixel_maximum_raw_median']:.3f} | "
                     f"{row['five_pixel_maximum_radial_mean_median']:.3f} | "
                     f"{row['median_absolute_score_change']:.3f} | "
                     f"{row['raw_radial_score_correlation_median']:.3f} | {heldout} |")
    lines.extend(["", "The variance ratio is the radial-mean result divided by the unchanged raw result. "
                  "The 11-pixel opposite-half column subtracts the same radial profile from the disjoint generic "
                  "patches before projection. Complete 31- and 47-pixel held-out stamps do not exist under the "
                  "frozen exclusion geometry."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), layout="constrained")
    for axis, support in zip(axes, raw.SUPPORTS):
        rows = [row for row in result["primary"] if row["support"] == support and
                row["band_role"] == "narrow"]
        x = np.arange(len(METHODS)); width = 0.36
        raw_values = [next(row["score_variance_raw_median"] for row in rows
                           if row["method"] == method) for method in METHODS]
        mean_values = [next(row["score_variance_radial_mean_median"] for row in rows
                            if row["method"] == method) for method in METHODS]
        axis.bar(x - width / 2, raw_values, width, label="raw")
        axis.bar(x + width / 2, mean_values, width, label="radial mean")
        axis.axhline(1, color="black", linestyle="--", linewidth=1)
        axis.set(title=f"{support}-pixel response", ylabel="primary score variance",
                 xticks=x, xticklabels=METHODS)
        axis.tick_params(axis="x", rotation=28, labelsize=8)
        axis.grid(axis="y", alpha=0.2)
        axis.legend(fontsize=8)
    fig.suptitle("KLIP mode 200: support-independent radial-mean control")
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, parent_root: Path, root: Path) -> None:
    """Run the radial-mean control against the corrected raw parent."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = raw.verify_bundle(bundle_path)
    parent_completion = read(parent_root / "complete.json")
    stage.verify([parent_completion[key] for key in ("protocol", "records", "results", "report", "figure")])
    parent_records = read(parent_root / "records.json")
    parent_lookup = {(row["support"], row["radius"], row["site_index"], row["search_index"],
                      row["band_role"], row["method"]): row for row in parent_records}
    stage.require(len(parent_lookup) == len(parent_records), "corrected parent record keys are not unique")

    baseline = np.asarray(arrays["baseline"], dtype=np.float64)
    finite = np.isfinite(baseline)
    planet = raw.planet_position(baseline.shape, metadata)
    mean_map, profile = radial_mean_map(baseline, planet)
    positions = {(int(row), int(column)): index
                 for index, (row, column) in enumerate(arrays["positions"])}

    root.mkdir(parents=True)
    protocol = {
        "schema": 1,
        "purpose": "support-independent radial-mean control for corrected KLIP Stage-B PSD filters",
        "mode": 200, "supports": list(raw.SUPPORTS), "radii": list(raw.RADII),
        "primary_radii": list(raw.PRIMARY_RADII), "methods": list(METHODS),
        "radial_mean": {"source": "signal-free baseline", "bin_width_pixels": 1,
                        "profile_radius_range": [0, MAXIMUM_PROFILE_RADIUS],
                        "interpolation": "linear with endpoint extension as in numpy.interp",
                        "known_planet_exclusion_radius": stage.PLANET_EXCLUSION_RADIUS + 0.5,
                        "application": "subtract interpolated native-pixel mean map from candidate stamps only"},
        "fixed_parent_components": "candidate mask, exact response, PSD, regularization, weights, and sigma",
        "bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
        "parent": stage.fingerprint(parent_root / "complete.json"),
        "scripts": [stage.fingerprint(Path(__file__)),
                    stage.fingerprint(Path(__file__).with_name("run_klip_stage_b_local_planet_masked_screen.py")),
                    stage.fingerprint(Path(__file__).with_name("run_klip_stage_b_local_noise_screen.py"))]}
    stage.write_json(root / "protocol.json", protocol)
    stage.write_json(root / "radial_profile.json", profile)
    stage.write_json(root / "state.json", {"status": "running", "started_unix": time.time()})

    records = []
    started = time.monotonic()
    for support in raw.SUPPORTS:
        half_support = support // 2
        for radius in raw.RADII:
            geometry = metadata["geometry"][str(support)][str(radius)]
            narrow = int(geometry["training"]["narrowest_split_supported_half_width"])
            for site_index, site in enumerate(geometry["selected_sites"]):
                searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                            for delta_row, delta_column in footprint.SEARCH_OFFSETS]
                for search_index, query in enumerate(searches):
                    rings = raw.training_stencils(finite, query, searches, planet, support)
                    source = positions[query]
                    template = raw.crop(np.asarray(arrays["response"][source], dtype=np.float64).T, support)
                    validity = raw.crop(np.asarray(arrays["validity"][source], dtype=bool).T, support)
                    candidate_mask = masked.candidate_mask(query, support, planet) & validity
                    row, column = query
                    candidate = baseline[column - half_support:column + half_support + 1,
                                         row - half_support:row + half_support + 1].T
                    candidate_mean = mean_map[column - half_support:column + half_support + 1,
                                              row - half_support:row + half_support + 1].T
                    candidate_radius = math.hypot(row - 0.5 * (baseline.shape[1] - 1),
                                                  column - 0.5 * (baseline.shape[0] - 1))
                    stage.require(candidate.shape == template.shape and np.all(np.isfinite(candidate_mean)) and
                                  candidate_radius + math.sqrt(2) * half_support < MAXIMUM_PROFILE_RADIUS,
                                  "candidate stamp exceeds the fitted radial-mean domain")
                    for band_role, width in (("narrow", narrow), ("full", 60)):
                        halves = []
                        for detector_half in (0, 1):
                            samples, interpolation_gain = raw.band_samples(
                                baseline, rings, width, detector_half)
                            opposite, _ = raw.band_samples(baseline, rings, width, 1 - detector_half)
                            rectangular = raw.fit_periodogram_base(samples, support, "rectangular")
                            hann = raw.fit_periodogram_base(samples, support, "hann")
                            models = {"identity": raw.regularize(rectangular, 1.0),
                                      "rectangular_m0.3": raw.regularize(rectangular, 0.3),
                                      "hann_m0.1": raw.regularize(hann, 0.1)}
                            mean_opposite = None
                            if support == 11:
                                mean_opposite, _ = raw.band_samples(mean_map, rings, width, 1 - detector_half)
                            results = {}
                            for method in METHODS:
                                parent_result, detail = masked.analyze_model(
                                    models[method], template, candidate, candidate_mask,
                                    opposite if support == 11 else None)
                                weight = np.asarray(detail["weight"])
                                sigma = float(detail["sigma"])
                                amplitude = float(weight @ (candidate - candidate_mean).ravel())
                                result = {"score_raw": parent_result["score_no_mean"],
                                          "amplitude_raw": parent_result["amplitude_no_mean"],
                                          "score_radial_mean": amplitude / sigma,
                                          "amplitude_radial_mean": amplitude,
                                          "radial_mean_projection": float(weight @ candidate_mean.ravel()),
                                          "sigma": sigma, "weight": weight,
                                          "samples": len(samples),
                                          "interpolation_variance_gain": interpolation_gain,
                                          "iterations": parent_result["iterations"],
                                          "relative_residual": parent_result["relative_residual"]}
                                if support == 11:
                                    result["heldout_radial_mean"] = raw.score_stats(
                                        (opposite - mean_opposite) @ weight / sigma)
                                results[method] = result
                            halves.append(results)
                        for method in METHODS:
                            first, second = halves[0][method], halves[1][method]
                            key = (support, radius, site_index, search_index, band_role, method)
                            parent = parent_lookup[key]
                            observed_raw = [first["score_raw"], second["score_raw"]]
                            stage.require(np.allclose(observed_raw, parent["score_no_mean"],
                                                      rtol=2e-11, atol=2e-12),
                                          "corrected raw parent score did not replay")
                            records.append({
                                "support": support, "radius": radius, "site_index": site_index,
                                "search_index": search_index, "row": row, "column": column,
                                "band_role": band_role, "band_half_width": width, "method": method,
                                "samples": [first["samples"], second["samples"]],
                                "interpolation_variance_gain": [first["interpolation_variance_gain"],
                                                                second["interpolation_variance_gain"]],
                                "weight_cosine": float(first["weight"] @ second["weight"] /
                                                       (np.linalg.norm(first["weight"]) *
                                                        np.linalg.norm(second["weight"]))),
                                "score_raw": observed_raw,
                                "score_radial_mean": [first["score_radial_mean"],
                                                      second["score_radial_mean"]],
                                "amplitude_raw": [first["amplitude_raw"], second["amplitude_raw"]],
                                "amplitude_radial_mean": [first["amplitude_radial_mean"],
                                                          second["amplitude_radial_mean"]],
                                "radial_mean_projection": [first["radial_mean_projection"],
                                                           second["radial_mean_projection"]],
                                "sigma": [first["sigma"], second["sigma"]],
                                "iterations": [first["iterations"], second["iterations"]],
                                "relative_residual": [first["relative_residual"],
                                                      second["relative_residual"]],
                                **({"heldout_radial_mean": [first["heldout_radial_mean"],
                                                            second["heldout_radial_mean"]]}
                                   if support == 11 else {})})
            print(f"completed support {support}, radius {radius:g}: {len(records)} records", flush=True)

    groups, primary = summarize(records)
    result = {"purpose": protocol["purpose"], "mode": 200, "records": len(records),
              "elapsed_seconds": time.monotonic() - started, "primary": primary, "groups": groups,
              "radial_profile": {"bins": len(profile),
                                 "sample_count_range": [min(row["samples"] for row in profile),
                                                        max(row["samples"] for row in profile)],
                                 "mean_range": [min(row["mean"] for row in profile),
                                                max(row["mean"] for row in profile)]}}
    stage.write_json(root / "records.json", records)
    stage.write_json(root / "results.json", result)
    write_report(root, result)
    completion = {"status": "complete", "protocol": stage.fingerprint(root / "protocol.json"),
                  "radial_profile": stage.fingerprint(root / "radial_profile.json"),
                  "records": stage.fingerprint(root / "records.json"),
                  "results": stage.fingerprint(root / "results.json"),
                  "report": stage.fingerprint(root / "results.md"),
                  "figure": stage.fingerprint(root / "comparison.png")}
    stage.write_json(root / "complete.json", completion)
    stage.write_json(root / "state.json", {"status": "complete", "elapsed_seconds": result["elapsed_seconds"]})
    print(root / "results.md", flush=True)


def check() -> None:
    """Verify constant-profile fitting and support-independent stamp evaluation."""
    image = np.full((128, 128), 3.25)
    planet = (75.5, 61.5)
    mean_map, profile = radial_mean_map(image, planet)
    stage.require(len(profile) == MAXIMUM_PROFILE_RADIUS and
                  np.allclose(mean_map, 3.25, rtol=0, atol=1e-14),
                  "constant radial mean did not reproduce its input")
    for support in raw.SUPPORTS:
        half = support // 2
        stamp = mean_map[64 - half:64 + half + 1, 64 - half:64 + half + 1]
        stage.require(stamp.shape == (support, support) and np.all(stamp == 3.25),
                      "radial mean is not support independent")
    print("KLIP Stage-B radial-mean checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the radial-mean control command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("parent", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the radial-mean control action."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check()
    else:
        run(arguments.bundle.resolve(), arguments.parent.resolve(), arguments.output.resolve())


if __name__ == "__main__":
    main()
