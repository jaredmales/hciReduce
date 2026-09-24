#!/usr/bin/env python3
"""Rerun the local KLIP raw PSD screen with a fixed candidate planet mask."""
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
from scipy.sparse.linalg import LinearOperator, cg

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_klip_stage_b_psd_extension as extension  # noqa: E402
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402
import run_klip_stage_b_local_noise_screen as raw  # noqa: E402


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def candidate_mask(query: tuple[int, int], support: int,
                   planet: tuple[float, float]) -> np.ndarray:
    """Return candidate coordinates outside the fixed native planet disk."""
    half = support // 2
    row, column = query
    native_row, native_column = np.mgrid[row - half:row + half + 1,
                                         column - half:column + half + 1]
    return np.hypot(native_row - planet[0], native_column - planet[1]) > stage.PLANET_EXCLUSION_RADIUS


def solve_masked_template(model: dict[str, object], template: np.ndarray,
                          mask: np.ndarray, relative_tolerance: float = 1e-10) -> dict[str, object]:
    """Solve the positive covariance principal submatrix on unmasked coordinates."""
    template = np.asarray(template, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    stage.require(template.shape == mask.shape and template.ndim == 2 and np.any(mask) and
                  np.all(np.isfinite(template[mask])), "invalid masked response template")
    support = template.shape[0]
    full_operator = extension.covariance_operator(model, support)
    selected = np.flatnonzero(mask.ravel())

    def matvec(vector: np.ndarray) -> np.ndarray:
        full = np.zeros(template.size, dtype=np.float64)
        full[selected] = np.asarray(vector)
        return np.asarray(full_operator @ full)[selected]

    operator = LinearOperator((len(selected), len(selected)), matvec=matvec, rmatvec=matvec,
                              dtype=np.float64)
    target = float(model["target_variance"])
    preconditioner = LinearOperator(operator.shape, matvec=lambda value: np.asarray(value) / target,
                                    rmatvec=lambda value: np.asarray(value) / target, dtype=np.float64)
    vector = template.ravel()[selected]
    iterations = 0

    def count_iteration(_: np.ndarray) -> None:
        nonlocal iterations
        iterations += 1

    inverse, information = cg(operator, vector, M=preconditioner, rtol=relative_tolerance,
                              atol=0, maxiter=4 * len(selected), callback=count_iteration)
    residual = operator @ inverse - vector
    relative_residual = float(np.linalg.norm(residual) / np.linalg.norm(vector))
    stage.require(information == 0 and relative_residual <= 5 * relative_tolerance,
                  "masked finite PSD solve did not converge")
    energy = float(vector @ inverse)
    stage.require(np.isfinite(energy) and energy > 0, "masked solve has nonpositive template energy")
    weight = np.zeros(template.size, dtype=np.float64)
    weight[selected] = inverse / energy
    return {"weight": weight, "energy": energy, "sigma": 1 / math.sqrt(energy),
            "iterations": iterations, "relative_residual": relative_residual,
            "valid_pixels": len(selected)}


def analyze_model(model: dict[str, object], template: np.ndarray, candidate: np.ndarray,
                  mask: np.ndarray, heldout: np.ndarray | None) -> tuple[dict[str, object], dict[str, object]]:
    """Score one planet-masked candidate and optional opposite-half patches."""
    solved = solve_masked_template(model, template, mask)
    weight = np.asarray(solved["weight"])
    vector = candidate.ravel()
    sigma = float(solved["sigma"])
    stage.require(np.isclose(weight @ template.ravel(), 1, rtol=2e-10, atol=2e-12) and
                  np.all(weight[~mask.ravel()] == 0), "masked weight lost its support contract")
    amplitude = float(weight @ vector)
    record = {"amplitude_no_mean": amplitude, "score_no_mean": amplitude / sigma,
              "sigma": sigma, "iterations": solved["iterations"],
              "relative_residual": solved["relative_residual"],
              "spectral_condition_number": float(np.max(model["power"]) / np.min(model["power"])),
              "target_variance": model["target_variance"], "psd_rescaling": model["psd_rescaling"],
              "valid_pixels": solved["valid_pixels"]}
    if template.shape == (11, 11):
        mean = np.asarray(model["mean"])
        mean_amplitude = float(weight @ (vector - mean))
        record.update({"amplitude_fitted_mean": mean_amplitude,
                       "score_fitted_mean": mean_amplitude / sigma})
        if heldout is not None:
            record["heldout_no_mean"] = raw.score_stats(heldout @ weight / sigma)
            record["heldout_fitted_mean"] = raw.score_stats((heldout - mean) @ weight / sigma)
    detail = {"weight": weight, "score_no_mean": record["score_no_mean"],
              "amplitude_no_mean": amplitude, "sigma": sigma}
    if "score_fitted_mean" in record:
        detail["score_fitted_mean"] = record["score_fitted_mean"]
    return record, detail


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the corrected support comparison and compact figure."""
    lines = ["# KLIP Stage-B planet-masked raw PSD screen", "",
             "This rerun keeps every original score-blind site but removes native pixels inside the fixed "
             "seven-pixel fitted-planet disk from candidate data and exact responses. Each PSD covariance is "
             "restricted to the resulting principal submatrix before inversion. Training and radial geometry "
             "retain their existing planet exclusion.", "", "## Controlling radii", "",
             "| Response | Band | Method | Masked score variance | Previous variance | Ratio | Split weight cosine | Median retained pixels | Median response energy |",
             "| ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    parent = {(row["support"], row["band_role"], row["method"]): row
              for row in result["parent_primary"]}
    for row in result["primary"]:
        previous = parent[(row["support"], row["band_role"], row["method"])]
        support = result["mask_support"][f"s{row['support']}_{row['band_role']}"]
        lines.append(f"| {row['support']} | {row['band_role']} | {row['method']} | "
                     f"{row['score_variance_median']:.3f} | {previous['score_variance_median']:.3f} | "
                     f"{row['score_variance_median']/previous['score_variance_median']:.3f} | "
                     f"{row['split_weight_cosine_median']:.3f} | {support['median_valid_pixels']:.0f} | "
                     f"{support['median_template_energy_fraction']:.3f} |")
    lines.extend(["", "The mask is fixed by the known planet coordinates and never by a candidate score. "
                  "The full finite PSD covariance remains positive; selecting the same valid rows and columns "
                  "forms a positive principal submatrix. Unit response is imposed after every masked solve."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    methods = [item[0] for item in raw.VARIANTS]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), layout="constrained")
    for axis, support in zip(axes, raw.SUPPORTS):
        values = []
        previous_values = []
        for method in methods:
            values.append(next(row["score_variance_median"] for row in result["primary"]
                               if row["support"] == support and row["band_role"] == "narrow" and
                               row["method"] == method))
            previous_values.append(parent[(support, "narrow", method)]["score_variance_median"])
        x = np.arange(len(methods)); width = 0.36
        axis.bar(x - width / 2, previous_values, width, label="unmasked")
        axis.bar(x + width / 2, values, width, label="planet masked")
        axis.axhline(1, color="black", linestyle="--", linewidth=1)
        axis.set(title=f"{support}-pixel response", ylabel="primary median score variance",
                 xticks=x, xticklabels=methods)
        axis.tick_params(axis="x", rotation=32, labelsize=8)
        axis.grid(axis="y", alpha=0.2)
        axis.legend(fontsize=8)
    fig.suptitle("KLIP mode 200 raw PSD: fixed planet mask in candidate data and response")
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, parent_root: Path, root: Path) -> None:
    """Run the corrected planet-masked raw support comparison."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = raw.verify_bundle(bundle_path)
    parent_completion = read(parent_root / "complete.json")
    stage.verify([parent_completion[key] for key in ("protocol", "records", "results", "report", "figure")])
    parent_result = read(parent_root / "results.json")
    root.mkdir(parents=True)
    protocol = {"schema": 1, "purpose": "corrected local raw KLIP screen with fixed candidate planet mask",
                "mode": 200, "supports": list(raw.SUPPORTS), "radii": list(raw.RADII),
                "primary_radii": list(raw.PRIMARY_RADII),
                "planet_mask": {"center_row_column": list(raw.planet_position(arrays["baseline"].shape, metadata)),
                                "native_radius": stage.PLANET_EXCLUSION_RADIUS,
                                "application": "remove candidate data and exact-response coordinates before solving the covariance principal submatrix"},
                "training": "unchanged 11-pixel raw Welch patches and existing candidate/planet exclusions",
                "bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
                "parent": stage.fingerprint(parent_root / "complete.json"),
                "scripts": [stage.fingerprint(Path(__file__)),
                            stage.fingerprint(Path(__file__).with_name("run_klip_stage_b_local_noise_screen.py")),
                            stage.fingerprint(Path(__file__).with_name("check_klip_stage_b_psd_extension.py"))]}
    stage.write_json(root / "protocol.json", protocol)
    stage.write_json(root / "state.json", {"status": "running", "started_unix": time.time()})

    baseline = np.asarray(arrays["baseline"], dtype=np.float64)
    finite = np.isfinite(baseline)
    planet = tuple(protocol["planet_mask"]["center_row_column"])
    lookup = {(int(row), int(column)): index for index, (row, column) in enumerate(arrays["positions"])}
    records = []
    started = time.monotonic()
    for support in raw.SUPPORTS:
        for radius in raw.RADII:
            geometry = metadata["geometry"][str(support)][str(radius)]
            narrow = int(geometry["training"]["narrowest_split_supported_half_width"])
            for site_index, site in enumerate(geometry["selected_sites"]):
                searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                            for delta_row, delta_column in footprint.SEARCH_OFFSETS]
                for search_index, query in enumerate(searches):
                    rings = raw.training_stencils(finite, query, searches, planet, support)
                    source = lookup[query]
                    template = raw.crop(np.asarray(arrays["response"][source], dtype=np.float64).T, support)
                    validity = raw.crop(np.asarray(arrays["validity"][source], dtype=bool).T, support)
                    mask = candidate_mask(query, support, planet) & validity
                    half_support = support // 2
                    row, column = query
                    candidate = baseline[column - half_support:column + half_support + 1,
                                         row - half_support:row + half_support + 1].T
                    stage.require(np.all(np.isfinite(candidate[mask])) and mask[half_support, half_support],
                                  "masked candidate lost valid source anchor")
                    full_energy = float(np.sum(np.square(template[validity])))
                    retained_energy = float(np.sum(np.square(template[mask])))
                    stage.require(full_energy > 0 and retained_energy > 0, "planet mask removed response energy")
                    for band_role, width in (("narrow", narrow), ("full", 60)):
                        halves = []
                        for detector_half in (0, 1):
                            samples, interpolation_gain = raw.band_samples(baseline, rings, width, detector_half)
                            opposite, _ = raw.band_samples(baseline, rings, width, 1 - detector_half)
                            rectangular = raw.fit_periodogram_base(samples, support, "rectangular")
                            hann = raw.fit_periodogram_base(samples, support, "hann")
                            models = {"identity": raw.regularize(rectangular, 1.0),
                                      "rectangular_m0.1": raw.regularize(rectangular, 0.1),
                                      "rectangular_m0.3": raw.regularize(rectangular, 0.3),
                                      "hann_m0.1": raw.regularize(hann, 0.1),
                                      "hann_m0.3": raw.regularize(hann, 0.3)}
                            results, details = {}, {}
                            for method, _, _ in raw.VARIANTS:
                                result, detail = analyze_model(models[method], template, candidate, mask,
                                                               opposite if support == 11 else None)
                                results[method] = result; details[method] = detail
                            halves.append({"samples": len(samples),
                                           "interpolation_variance_gain": interpolation_gain,
                                           "results": results, "details": details})
                        for method, _, _ in raw.VARIANTS:
                            first = halves[0]["details"][method]; second = halves[1]["details"][method]
                            cosine = float(first["weight"] @ second["weight"] /
                                           (np.linalg.norm(first["weight"]) * np.linalg.norm(second["weight"])))
                            record = {"support": support, "radius": radius, "site_index": site_index,
                                      "search_index": search_index, "row": row, "column": column,
                                      "band_role": band_role, "band_half_width": width, "method": method,
                                      "samples": [halves[index]["samples"] for index in (0, 1)],
                                      "interpolation_variance_gain": [halves[index]["interpolation_variance_gain"]
                                                                       for index in (0, 1)],
                                      "weight_cosine": cosine,
                                      "valid_pixels": int(np.count_nonzero(mask)),
                                      "template_energy_fraction": retained_energy / full_energy,
                                      "score_no_mean": [halves[index]["results"][method]["score_no_mean"]
                                                        for index in (0, 1)],
                                      "amplitude_no_mean": [halves[index]["results"][method]["amplitude_no_mean"]
                                                            for index in (0, 1)],
                                      "sigma": [halves[index]["results"][method]["sigma"] for index in (0, 1)],
                                      "iterations": [halves[index]["results"][method]["iterations"]
                                                     for index in (0, 1)],
                                      "relative_residual": [halves[index]["results"][method]["relative_residual"]
                                                            for index in (0, 1)],
                                      "spectral_condition_number": [
                                          halves[index]["results"][method]["spectral_condition_number"]
                                          for index in (0, 1)]}
                            if support == 11:
                                record["score_fitted_mean"] = [
                                    halves[index]["results"][method]["score_fitted_mean"] for index in (0, 1)]
                                record["heldout_no_mean"] = [
                                    halves[index]["results"][method]["heldout_no_mean"] for index in (0, 1)]
                                record["heldout_fitted_mean"] = [
                                    halves[index]["results"][method]["heldout_fitted_mean"] for index in (0, 1)]
                            records.append(record)
            print(f"completed support {support}, radius {radius:g}: {len(records)} records", flush=True)

    summary = raw.summarize_records(records)
    primary = raw.pooled_primary(summary)
    support_summary = {}
    for support in raw.SUPPORTS:
        for band_role in ("narrow", "full"):
            chosen = [record for record in records if record["support"] == support and
                      record["band_role"] == band_role and record["radius"] in raw.PRIMARY_RADII and
                      record["method"] == "identity"]
            support_summary[f"s{support}_{band_role}"] = {
                "records": len(chosen), "minimum_valid_pixels": min(row["valid_pixels"] for row in chosen),
                "median_valid_pixels": float(np.median([row["valid_pixels"] for row in chosen])),
                "maximum_valid_pixels": max(row["valid_pixels"] for row in chosen),
                "minimum_template_energy_fraction": min(row["template_energy_fraction"] for row in chosen),
                "median_template_energy_fraction": float(np.median([
                    row["template_energy_fraction"] for row in chosen]))}
    result = {"purpose": protocol["purpose"], "mode": 200, "records": len(records),
              "elapsed_seconds": time.monotonic() - started, "primary": primary,
              "groups": summary["groups"], "mask_support": support_summary,
              "parent_primary": parent_result["primary"]}
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
    """Verify masked iterative solves against dense covariance principal submatrices."""
    generator = np.random.default_rng(240924)
    samples = generator.normal(size=(40, 121))
    for support in (11, 31):
        model = extension.fit_extended_psd(samples, support, "hann", 0.1)
        template = generator.normal(size=(support, support))
        mask = generator.random((support, support)) > 0.17
        observed = solve_masked_template(model, template, mask)
        covariance = extension.dense_covariance(model, support)
        selected = np.flatnonzero(mask.ravel())
        direct = np.linalg.solve(covariance[np.ix_(selected, selected)], template.ravel()[selected])
        direct /= float(template.ravel()[selected] @ direct)
        stage.require(np.allclose(observed["weight"][selected], direct, rtol=2e-8, atol=2e-10) and
                      np.all(observed["weight"][~mask.ravel()] == 0),
                      "masked iterative weights differ from dense principal-submatrix solve")
        isotropic = extension.fit_extended_psd(samples, support, "rectangular", 1.0)
        endpoint = solve_masked_template(isotropic, template, mask)
        expected = np.zeros(template.size); vector = template.ravel()[selected]
        expected[selected] = vector / float(vector @ vector)
        stage.require(np.allclose(endpoint["weight"], expected, rtol=2e-11, atol=2e-13),
                      "masked isotropic endpoint differs from identity weighting")
    print("KLIP Stage-B local planet-masked screen checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the planet-masked screen command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("parent", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the requested planet-masked screen action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        run(args.bundle.resolve(), args.parent.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
