#!/usr/bin/env python3
"""Run the direct 11-pixel diagonal/PCA KLIP Stage-B covariance screen."""
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
from scipy.linalg import cho_factor, cho_solve, qr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402
import run_klip_stage_b_local_noise_screen as raw  # noqa: E402
import run_klip_stage_b_local_planet_masked_screen as masked  # noqa: E402
import run_klip_stage_b_local_radial_normalization as radial  # noqa: E402


SUPPORT = 11
BAND_WIDTHS = (0, 5, 10, 20, 40, 60)
COORDINATES = ("raw", "radial")
BANDS = ("narrow", "wide")
FLOORS = (0.1, 0.3, 1.0)
RANKS: tuple[int | str, ...] = (0, 3, 8, "all")
METHODS = (("diagonal_f0.1", "diagonal", None, 0.1),) + tuple(
    (f"pca_r{rank}_f{floor:g}", "pca", rank, floor)
    for rank in RANKS for floor in FLOORS)


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def fit_direct_grid(samples: np.ndarray) -> dict[str, dict[str, object]]:
    """Fit the frozen diagonal and low-rank PCA grid from one sample matrix."""
    samples = np.asarray(samples, dtype=np.float64)
    stage.require(samples.ndim == 2 and samples.shape[1] == SUPPORT * SUPPORT and
                  len(samples) >= 8 and np.all(np.isfinite(samples)),
                  "invalid direct-covariance training samples")
    mean = np.mean(samples, axis=0)
    centered = samples - mean
    pixel_variances = np.sum(np.square(centered), axis=0) / (len(samples) - 1)
    median_variance = float(np.median(pixel_variances))
    sample_trace = float(np.sum(pixel_variances))
    stage.require(np.isfinite(median_variance) and median_variance > 0 and
                  np.isfinite(sample_trace) and sample_trace > 0,
                  "invalid direct-covariance sample variance")

    _, singular_values, right_vectors = np.linalg.svd(centered, full_matrices=False)
    maximum_estimable = min(len(samples) - 1, samples.shape[1])
    eigenvalues = np.square(singular_values[:maximum_estimable]) / (len(samples) - 1)
    eigenvectors = right_vectors[:maximum_estimable].T
    stage.require(np.all(np.isfinite(eigenvalues)) and
                  np.all(eigenvalues[:-1] >= eigenvalues[1:]) and
                  np.isclose(np.sum(eigenvalues), sample_trace, rtol=3e-12, atol=1e-14),
                  "direct PCA eigensystem lost the sample covariance trace")

    models = {}
    diagonal_floor = FLOORS[0] * median_variance
    diagonal = np.maximum(pixel_variances, diagonal_floor)
    models["diagonal_f0.1"] = {
        "kind": "diagonal", "mean": mean, "variances": diagonal,
        "factor": np.empty((samples.shape[1], 0)), "floor_fraction": FLOORS[0],
        "variance_floor": diagonal_floor, "requested_modes": None,
        "retained_modes": 0, "estimable_modes": maximum_estimable,
        "sample_trace": sample_trace, "model_trace": float(np.sum(diagonal)),
        "condition_number": float(np.max(diagonal) / np.min(diagonal))}

    for rank in RANKS:
        cap = maximum_estimable if rank == "all" else min(int(rank), maximum_estimable)
        for floor_fraction in FLOORS:
            floor = floor_fraction * median_variance
            retained = int(np.count_nonzero(eigenvalues[:cap] > floor))
            modes = eigenvectors[:, :retained]
            excess = np.sqrt(np.maximum(eigenvalues[:retained] - floor, 0))
            factor = modes * excess
            model_trace = samples.shape[1] * floor + float(np.sum(
                eigenvalues[:retained] - floor))
            largest = max(float(eigenvalues[0]) if retained else floor, floor)
            name = f"pca_r{rank}_f{floor_fraction:g}"
            models[name] = {
                "kind": "pca", "mean": mean,
                "variances": np.full(samples.shape[1], floor), "factor": factor,
                "floor_fraction": floor_fraction, "variance_floor": floor,
                "requested_modes": rank, "retained_modes": retained,
                "estimable_modes": maximum_estimable, "sample_trace": sample_trace,
                "model_trace": model_trace, "condition_number": largest / floor}
    stage.require(set(models) == {method[0] for method in METHODS},
                  "direct-covariance method grid changed")
    return models


def solve_model(model: dict[str, object], value: np.ndarray,
                support_mask: np.ndarray) -> np.ndarray:
    """Solve a covariance principal submatrix using production's scaled-QR algebra."""
    value = np.asarray(value, dtype=np.float64)
    mask = np.asarray(support_mask, dtype=bool).ravel()
    selected = np.flatnonzero(mask)
    stage.require(value.shape == (len(selected),) and np.all(np.isfinite(value)),
                  "direct-covariance solve has invalid support")
    variances = np.asarray(model["variances"], dtype=np.float64)[selected]
    factor = np.asarray(model["factor"], dtype=np.float64)[selected]
    scale = 1 / np.sqrt(variances)
    scaled = scale * value
    if factor.shape[1]:
        whitened = scale[:, None] * factor
        rotation, triangular = qr(whitened, mode="full", check_finite=False)
        count = min(whitened.shape)
        upper = np.triu(triangular[:count])
        reduced = np.eye(count) + upper @ upper.T
        decomposition = cho_factor(reduced, lower=True, check_finite=False)
        rotated = rotation.T @ scaled
        rotated[:count] = cho_solve(decomposition, rotated[:count], check_finite=False)
        scaled = rotation @ rotated
    result = scale * scaled
    stage.require(np.all(np.isfinite(result)), "direct-covariance solve is nonfinite")
    return result


def analyze_model(model: dict[str, object], template: np.ndarray, candidate: np.ndarray,
                  support_mask: np.ndarray, heldout: np.ndarray) -> tuple[dict[str, object], dict[str, object]]:
    """Score one candidate and the opposite detector-half patches."""
    template = np.asarray(template, dtype=np.float64)
    candidate = np.asarray(candidate, dtype=np.float64)
    mask = np.asarray(support_mask, dtype=bool)
    selected = np.flatnonzero(mask.ravel())
    inverse = solve_model(model, template.ravel()[selected], mask)
    energy = float(template.ravel()[selected] @ inverse)
    stage.require(np.isfinite(energy) and energy > 0, "direct covariance has invalid response energy")
    weight = np.zeros(template.size, dtype=np.float64)
    weight[selected] = inverse / energy
    sigma = 1 / math.sqrt(energy)
    stage.require(np.isclose(weight @ template.ravel(), 1, rtol=2e-11, atol=2e-13) and
                  np.all(weight[~mask.ravel()] == 0),
                  "direct covariance lost the masked unit-response contract")
    vector = candidate.ravel()
    mean = np.asarray(model["mean"], dtype=np.float64)
    amplitude_no_mean = float(weight @ vector)
    amplitude_fitted_mean = float(weight @ (vector - mean))
    record = {
        "amplitude_no_mean": amplitude_no_mean,
        "score_no_mean": amplitude_no_mean / sigma,
        "amplitude_fitted_mean": amplitude_fitted_mean,
        "score_fitted_mean": amplitude_fitted_mean / sigma,
        "sigma": sigma,
        "heldout_no_mean": raw.score_stats(np.asarray(heldout) @ weight / sigma),
        "heldout_fitted_mean": raw.score_stats((np.asarray(heldout) - mean) @ weight / sigma),
        "retained_modes": model["retained_modes"],
        "estimable_modes": model["estimable_modes"],
        "variance_floor": model["variance_floor"],
        "condition_number": model["condition_number"],
        "model_trace_fraction": model["model_trace"] / model["sample_trace"]}
    return record, {"weight": weight, "sigma": sigma}


def summarize(records: list[dict[str, object]]) -> dict[str, object]:
    """Aggregate candidate calibration, held-out prediction, and split stability."""
    groups = {}
    for coordinates in COORDINATES:
        for radius in raw.RADII:
            for band_role in BANDS:
                for method, _, _, _ in METHODS:
                    chosen = [row for row in records if row["coordinates"] == coordinates and
                              row["radius"] == radius and row["band_role"] == band_role and
                              row["method"] == method]
                    stage.require(len(chosen) == 60,
                                  "direct-covariance summary lost fixed query records")
                    no_mean = np.asarray([value for row in chosen for value in row["score_no_mean"]])
                    fitted_mean = np.asarray([value for row in chosen
                                              for value in row["score_fitted_mean"]])
                    centers = np.asarray([value for row in chosen if row["search_index"] == 0
                                          for value in row["score_no_mean"]])
                    maximums = []
                    for site_index in range(12):
                        site = sorted((row for row in chosen if row["site_index"] == site_index),
                                      key=lambda row: row["search_index"])
                        stage.require(len(site) == 5,
                                      "direct-covariance five-pixel search is incomplete")
                        for detector_half in (0, 1):
                            maximums.append(max(row["score_no_mean"][detector_half]
                                                for row in site))
                    first = np.asarray([row["score_no_mean"][0] for row in chosen])
                    second = np.asarray([row["score_no_mean"][1] for row in chosen])
                    entry = {
                        "coordinates": coordinates, "support": SUPPORT, "radius": radius,
                        "band_role": band_role, "band_half_width": chosen[0]["band_half_width"],
                        "method": method,
                        "directional_candidate_scores_no_mean": raw.score_stats(no_mean),
                        "directional_candidate_scores_fitted_mean": raw.score_stats(fitted_mean),
                        "directional_center_scores_no_mean": raw.score_stats(centers),
                        "five_pixel_maximum_no_mean": raw.score_stats(np.asarray(maximums)),
                        "median_split_physical_weight_cosine": float(np.median([
                            row["physical_weight_cosine"] for row in chosen])),
                        "split_physical_weight_cosine_10_90": np.percentile([
                            row["physical_weight_cosine"] for row in chosen], [10, 90]).tolist(),
                        "split_score_correlation": float(np.corrcoef(first, second)[0, 1]),
                        "median_absolute_split_score_difference": float(np.median(
                            np.abs(first - second))),
                        "training_patch_count_range": [
                            min(value for row in chosen for value in row["samples"]),
                            max(value for row in chosen for value in row["samples"])],
                        "retained_modes_median": float(np.median([
                            value for row in chosen for value in row["retained_modes"]])),
                        "retained_modes_range": [
                            min(value for row in chosen for value in row["retained_modes"]),
                            max(value for row in chosen for value in row["retained_modes"])],
                        "estimable_modes_range": [
                            min(value for row in chosen for value in row["estimable_modes"]),
                            max(value for row in chosen for value in row["estimable_modes"])],
                        "median_condition_number": float(np.median([
                            value for row in chosen for value in row["condition_number"]])),
                        "median_model_trace_fraction": float(np.median([
                            value for row in chosen for value in row["model_trace_fraction"]])),
                        "profile_undefined_training_patches_range": [
                            min(value for row in chosen for value in row[
                                "profile_undefined_training_patches"]),
                            max(value for row in chosen for value in row[
                                "profile_undefined_training_patches"])],
                        "profile_undefined_heldout_patches_range": [
                            min(value for row in chosen for value in row[
                                "profile_undefined_heldout_patches"]),
                            max(value for row in chosen for value in row[
                                "profile_undefined_heldout_patches"])]}
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
                                item["within_one_fraction"] for item in values])),
                            "median_within_two_fraction": float(np.median([
                                item["within_two_fraction"] for item in values]))}
                    groups[f"{coordinates}_r{radius:g}_{band_role}_{method}"] = entry

    primary = []
    for coordinates in COORDINATES:
        for band_role in BANDS:
            for method, kind, rank, floor in METHODS:
                selected = [groups[f"{coordinates}_r{radius:g}_{band_role}_{method}"]
                            for radius in raw.PRIMARY_RADII]
                primary.append({
                    "coordinates": coordinates, "band_role": band_role,
                    "method": method, "kind": kind, "requested_modes": rank,
                    "floor_fraction": floor,
                    "score_variance_no_mean_median": float(np.median([
                        row["directional_candidate_scores_no_mean"]["variance"]
                        for row in selected])),
                    "score_variance_fitted_mean_median": float(np.median([
                        row["directional_candidate_scores_fitted_mean"]["variance"]
                        for row in selected])),
                    "absolute_score_mean_no_mean_median": float(np.median([
                        abs(row["directional_candidate_scores_no_mean"]["mean"])
                        for row in selected])),
                    "heldout_variance_no_mean_median": float(np.median([
                        row["heldout_no_mean"]["median_variance_over_prediction"]
                        for row in selected])),
                    "heldout_variance_fitted_mean_median": float(np.median([
                        row["heldout_fitted_mean"]["median_variance_over_prediction"]
                        for row in selected])),
                    "split_physical_weight_cosine_median": float(np.median([
                        row["median_split_physical_weight_cosine"] for row in selected])),
                    "split_score_correlation_median": float(np.median([
                        row["split_score_correlation"] for row in selected])),
                    "retained_modes_median": float(np.median([
                        row["retained_modes_median"] for row in selected])),
                    "retained_modes_minimum": min(
                        row["retained_modes_range"][0] for row in selected),
                    "retained_modes_maximum": max(
                        row["retained_modes_range"][1] for row in selected),
                    "condition_number_median": float(np.median([
                        row["median_condition_number"] for row in selected])),
                    "model_trace_fraction_median": float(np.median([
                        row["median_model_trace_fraction"] for row in selected]))})
    return {"groups": groups, "primary": primary}


def primary_references(raw_parent: dict[str, object],
                       radial_parent: dict[str, object]) -> dict[str, object]:
    """Copy the corrected identity and PSD reference rows used for interpretation."""
    raw_rows = {(row["band_role"], row["method"]): row for row in raw_parent["primary"]
                if row["support"] == SUPPORT}
    radial_rows = {row["method"]: row for row in radial_parent["primary"]}
    return {
        "raw_narrow_identity": raw_rows[("narrow", "identity")],
        "raw_narrow_rectangular_m0.3": raw_rows[("narrow", "rectangular_m0.3")],
        "raw_narrow_hann_m0.1": raw_rows[("narrow", "hann_m0.1")],
        "radial_identity": radial_rows["identity"],
        "radial_rectangular_m0.3": radial_rows["rectangular_m0.3"],
        "radial_hann_m0.1": radial_rows["hann_m0.1"]}


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the complete direct-covariance primary grid and heatmap."""
    lines = ["# KLIP Stage-B direct diagonal/PCA covariance screen", "",
             "This mode-200 screen uses the corrected 11-pixel candidate planet mask and independent detector-half fits. Direct covariance is estimated in raw and strict radial-standardized coordinates. Each geometry uses its narrowest split-supported training band and the next wider fixed half-width.", "",
             "PCA covariance is an isotropic floor plus at most 0, 3, 8, or every estimable centered sample mode. The floor is 0.1, 0.3, or 1.0 times the median fitted pixel variance. Diagonal covariance uses independently fitted pixel variances clipped at the production 0.1 floor. Candidate mean subtraction is reported both off and on.", "",
             "## Controlling radii: narrow bands", "",
             "| Coordinates | Method | Floor | Modes retained | Variance, no mean | Variance, fitted mean | Opposite-half variance | Split weight cosine | Model/sample trace |",
             "| :--- | :--- | ---: | :--- | ---: | ---: | ---: | ---: | ---: |"]
    for coordinates in COORDINATES:
        rows = [row for row in result["primary"] if row["coordinates"] == coordinates and
                row["band_role"] == "narrow"]
        for row in rows:
            modes = ("—" if row["kind"] == "diagonal" else
                     f"{row['retained_modes_median']:.0f} [{row['retained_modes_minimum']}, {row['retained_modes_maximum']}]")
            lines.append(f"| {coordinates} | {row['method']} | {row['floor_fraction']:.1f} | {modes} | "
                         f"{row['score_variance_no_mean_median']:.3f} | "
                         f"{row['score_variance_fitted_mean_median']:.3f} | "
                         f"{row['heldout_variance_no_mean_median']:.3f} | "
                         f"{row['split_physical_weight_cosine_median']:.3f} | "
                         f"{row['model_trace_fraction_median']:.3f} |")
    refs = result["references"]
    lines.extend(["", "Corrected reference score variances are "
                  f"{refs['raw_narrow_identity']['score_variance_median']:.3f} for raw identity, "
                  f"{refs['raw_narrow_rectangular_m0.3']['score_variance_median']:.3f} for raw rectangular/mixing-0.3, "
                  f"{refs['raw_narrow_hann_m0.1']['score_variance_median']:.3f} for raw Hann/mixing-0.1, and "
                  f"{refs['radial_hann_m0.1']['normalized_score_variance_median']:.3f} for radial Hann/mixing-0.1."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), layout="constrained")
    rank_labels = ("0", "3", "8", "all")
    for row_index, coordinates in enumerate(COORDINATES):
        for column_index, mean_key in enumerate(("score_variance_no_mean_median",
                                                  "score_variance_fitted_mean_median")):
            values = np.full((len(RANKS), len(FLOORS)), np.nan)
            for y, rank in enumerate(RANKS):
                for x, floor in enumerate(FLOORS):
                    entry = next(row for row in result["primary"]
                                 if row["coordinates"] == coordinates and
                                 row["band_role"] == "narrow" and
                                 row["method"] == f"pca_r{rank}_f{floor:g}")
                    values[y, x] = entry[mean_key]
            axis = axes[row_index, column_index]
            image = axis.imshow(values, cmap="viridis", aspect="auto")
            axis.set(title=f"{coordinates}; {'no candidate mean' if column_index == 0 else 'fitted candidate mean'}",
                     xticks=range(len(FLOORS)), xticklabels=FLOORS,
                     yticks=range(len(RANKS)), yticklabels=rank_labels,
                     xlabel="floor fraction", ylabel="requested PCA modes")
            for y in range(len(RANKS)):
                for x in range(len(FLOORS)):
                    axis.text(x, y, f"{values[y, x]:.2f}", ha="center", va="center",
                              color="white" if values[y, x] < np.nanmedian(values) else "black")
            fig.colorbar(image, ax=axis, shrink=0.8, label="primary median score variance")
    fig.suptitle("KLIP mode 200 direct PCA covariance, corrected 11-pixel narrow bands")
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, raw_parent_root: Path, radial_parent_root: Path,
        root: Path) -> None:
    """Run the frozen direct-covariance screen from verified local products."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = raw.verify_bundle(bundle_path)
    parents = []
    for parent_root in (raw_parent_root, radial_parent_root):
        completion = read(parent_root / "complete.json")
        stage.verify([completion[key] for key in ("protocol", "records", "results", "report", "figure")])
        parents.append(read(parent_root / "results.json"))
    root.mkdir(parents=True)
    protocol = {
        "schema": 1, "purpose": "direct 11-pixel diagonal/PCA KLIP Stage-B covariance screen",
        "mode": 200, "support": SUPPORT, "radii": list(raw.RADII),
        "primary_radii": list(raw.PRIMARY_RADII), "coordinates": list(COORDINATES),
        "bands": "narrowest split-supported and next wider fixed half-width; radially standardized wide-band patches touching the undefined radius >=60 profile domain are rejected deterministically",
        "band_width_grid": list(BAND_WIDTHS),
        "diagonal": {"floor_fraction": 0.1, "model": "per-pixel sample variance clipped at floor"},
        "pca": {"requested_modes": list(RANKS), "floor_fractions": list(FLOORS),
                "model": "isotropic floor plus retained centered sample-eigenmode excess variance"},
        "fitted_mean": "candidate scores reported with and without the shared fitted training mean",
        "candidate_mask": "fixed seven-pixel planet disk removed from candidate and exact response through covariance principal-submatrix solve",
        "bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
        "raw_parent": stage.fingerprint(raw_parent_root / "complete.json"),
        "radial_parent": stage.fingerprint(radial_parent_root / "complete.json"),
        "scripts": [stage.fingerprint(Path(__file__)),
                    stage.fingerprint(Path(__file__).with_name(
                        "run_klip_stage_b_local_planet_masked_screen.py")),
                    stage.fingerprint(Path(__file__).with_name(
                        "run_klip_stage_b_local_planet_masked_normalization.py"))]}
    stage.write_json(root / "protocol.json", protocol)
    stage.write_json(root / "state.json", {"status": "running", "started_unix": time.time()})

    baseline = np.asarray(arrays["baseline"], dtype=np.float64)
    finite = np.isfinite(baseline)
    planet = raw.planet_position(baseline.shape, metadata)
    response_lookup = {(int(row), int(column)): index
                       for index, (row, column) in enumerate(arrays["positions"])}
    records = []
    started = time.monotonic()
    for radius in raw.RADII:
        fixed = metadata["geometry"]["11"][str(radius)]
        narrow = int(fixed["training"]["narrowest_split_supported_half_width"])
        wider = next(width for width in BAND_WIDTHS if width > narrow)
        for site_index, site in enumerate(fixed["selected_sites"]):
            searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                        for delta_row, delta_column in footprint.SEARCH_OFFSETS]
            excluded = radial.exclusion_mask(baseline.shape, searches, SUPPORT, planet)
            scale_map, _ = radial.variance_profile(baseline, excluded)
            coordinate_images = {"raw": baseline, "radial": baseline / scale_map}
            for search_index, query in enumerate(searches):
                rings = raw.training_stencils(finite, query, searches, planet, SUPPORT)
                source = response_lookup[query]
                template = raw.crop(np.asarray(arrays["response"][source], dtype=np.float64).T,
                                    SUPPORT)
                validity = raw.crop(np.asarray(arrays["validity"][source], dtype=bool).T, SUPPORT)
                planet_mask = masked.candidate_mask(query, SUPPORT, planet) & validity
                row, column = query
                candidate = baseline[column - 5:column + 6, row - 5:row + 6].T
                candidate_scale = scale_map[column - 5:column + 6, row - 5:row + 6].T
                stage.require(np.all(np.isfinite(candidate[planet_mask])) and planet_mask[5, 5],
                              "direct-covariance candidate lost its source anchor")
                for coordinates in COORDINATES:
                    scale = np.ones_like(candidate_scale) if coordinates == "raw" else candidate_scale
                    fit_template = template / scale
                    fit_candidate = candidate / scale
                    image = coordinate_images[coordinates]
                    for band_role, width in (("narrow", narrow), ("wide", wider)):
                        halves = []
                        for detector_half in (0, 1):
                            samples, interpolation_gain = raw.band_samples(
                                image, rings, width, detector_half)
                            opposite, _ = raw.band_samples(
                                image, rings, width, 1 - detector_half)
                            finite_training = np.all(np.isfinite(samples), axis=1)
                            finite_opposite = np.all(np.isfinite(opposite), axis=1)
                            rejected_training = int(np.count_nonzero(~finite_training))
                            rejected_opposite = int(np.count_nonzero(~finite_opposite))
                            if coordinates == "raw":
                                stage.require(rejected_training == 0 and rejected_opposite == 0,
                                              "raw direct covariance lost finite training support")
                            samples = samples[finite_training]
                            opposite = opposite[finite_opposite]
                            stage.require(len(samples) >= 8 and len(opposite) >= 2,
                                          "strict radial profile leaves insufficient wide-band support")
                            models = fit_direct_grid(samples)
                            results = {}
                            details = {}
                            for method, _, _, _ in METHODS:
                                try:
                                    result, detail = analyze_model(
                                        models[method], fit_template, fit_candidate,
                                        planet_mask, opposite)
                                except Exception as error:
                                    raise RuntimeError(
                                        f"direct model failed at radius={radius:g}, site={site_index}, "
                                        f"search={search_index}, coordinates={coordinates}, "
                                        f"band={band_role}, half={detector_half}, method={method}: {error}") from error
                                results[method] = result
                                details[method] = detail | {
                                    "physical_weight": detail["weight"] / scale.ravel()}
                            halves.append({"samples": len(samples),
                                           "interpolation_variance_gain": interpolation_gain,
                                           "profile_undefined_training_patches": rejected_training,
                                           "profile_undefined_heldout_patches": rejected_opposite,
                                           "results": results, "details": details})
                        for method, kind, rank, floor in METHODS:
                            first = halves[0]["details"][method]
                            second = halves[1]["details"][method]
                            cosine = float(first["physical_weight"] @ second["physical_weight"] /
                                           (np.linalg.norm(first["physical_weight"]) *
                                            np.linalg.norm(second["physical_weight"])))
                            records.append({
                                "coordinates": coordinates, "support": SUPPORT,
                                "radius": radius, "site_index": site_index,
                                "search_index": search_index, "row": row, "column": column,
                                "band_role": band_role, "band_half_width": width,
                                "method": method, "kind": kind, "requested_modes": rank,
                                "floor_fraction": floor,
                                "samples": [halves[index]["samples"] for index in (0, 1)],
                                "profile_undefined_training_patches": [halves[index][
                                    "profile_undefined_training_patches"] for index in (0, 1)],
                                "profile_undefined_heldout_patches": [halves[index][
                                    "profile_undefined_heldout_patches"] for index in (0, 1)],
                                "interpolation_variance_gain": [halves[index][
                                    "interpolation_variance_gain"] for index in (0, 1)],
                                "physical_weight_cosine": cosine,
                                "valid_pixels": int(np.count_nonzero(planet_mask)),
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
                                "retained_modes": [halves[index]["results"][method][
                                    "retained_modes"] for index in (0, 1)],
                                "estimable_modes": [halves[index]["results"][method][
                                    "estimable_modes"] for index in (0, 1)],
                                "condition_number": [halves[index]["results"][method][
                                    "condition_number"] for index in (0, 1)],
                                "model_trace_fraction": [halves[index]["results"][method][
                                    "model_trace_fraction"] for index in (0, 1)]})
        print(f"completed radius {radius:g}: {len(records)} records", flush=True)

    summarized = summarize(records)
    result = {"purpose": protocol["purpose"], "mode": 200, "support": SUPPORT,
              "records": len(records), "elapsed_seconds": time.monotonic() - started,
              "primary": summarized["primary"], "groups": summarized["groups"],
              "references": primary_references(parents[0], parents[1])}
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
    """Check model construction and restricted solves against dense covariance."""
    generator = np.random.default_rng(250927)
    samples = generator.normal(size=(43, SUPPORT * SUPPORT))
    samples[:, :4] += 2 * generator.normal(size=(43, 1))
    models = fit_direct_grid(samples)
    stage.require(models["pca_r0_f0.1"]["retained_modes"] == 0 and
                  models["pca_r3_f0.1"]["retained_modes"] <= 3 and
                  models["pca_r8_f0.1"]["retained_modes"] <= 8 and
                  models["pca_rall_f0.1"]["estimable_modes"] == len(samples) - 1,
                  "direct-covariance rank grid changed")
    template = generator.normal(size=(SUPPORT, SUPPORT))
    mask = generator.random((SUPPORT, SUPPORT)) > 0.2
    selected = np.flatnonzero(mask.ravel())
    for name, model in models.items():
        covariance = np.diag(model["variances"]) + model["factor"] @ model["factor"].T
        observed = solve_model(model, template.ravel()[selected], mask)
        expected = np.linalg.solve(covariance[np.ix_(selected, selected)],
                                   template.ravel()[selected])
        stage.require(np.allclose(observed, expected, rtol=2e-10, atol=2e-12),
                      f"restricted direct solve differs from dense covariance: {name}")
        amplitude = 0.031
        measured, _ = analyze_model(model, template, amplitude * template,
                                    mask, samples)
        stage.require(np.isclose(measured["amplitude_no_mean"], amplitude,
                                 rtol=2e-12, atol=2e-14),
                      f"direct covariance changed physical amplitude: {name}")
    print("KLIP Stage-B direct diagonal/PCA checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the direct-covariance command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("raw_parent", type=Path)
    run_parser.add_argument("radial_parent", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the direct-covariance action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        run(args.bundle.resolve(), args.raw_parent.resolve(),
            args.radial_parent.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
