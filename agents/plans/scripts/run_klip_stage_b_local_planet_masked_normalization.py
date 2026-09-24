#!/usr/bin/env python3
"""Apply strict radial normalization with the fixed KLIP candidate planet mask."""
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
import run_klip_stage_b_local_planet_masked_screen as masked  # noqa: E402
import run_klip_stage_b_local_radial_normalization as normalized  # noqa: E402


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the corrected normalized comparison."""
    lines = ["# KLIP Stage-B planet-masked radial-normalization screen", "",
             "This result combines the strict leave-site-out radial profile with the fixed native planet mask "
             "in candidate data and exact 11-pixel responses. Training already excludes both the candidate "
             "footprints and planet disk.", "", "## Controlling radii", "",
             "| Method | Masked raw variance | Masked normalized variance | Normalized/raw | Opposite-half variance | Split physical-weight cosine | Within ±1σ |",
             "| :--- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in result["primary"]:
        lines.append(f"| {row['method']} | {row['raw_score_variance_median']:.3f} | "
                     f"{row['normalized_score_variance_median']:.3f} | "
                     f"{row['normalized_over_raw_variance_median']:.3f} | "
                     f"{row['normalized_heldout_variance_median']:.3f} | "
                     f"{row['normalized_split_weight_cosine_median']:.3f} | "
                     f"{row['normalized_within_one_fraction_median']:.3f} |")
    lines.extend(["", "Every site retains its original score-blind selection. Candidate/template coordinates "
                  "inside the planet disk have zero weight, and every standardized solve is renormalized to unit "
                  "response on the remaining support."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    methods = [item[0] for item in raw.VARIANTS]
    raw_values = [row["raw_score_variance_median"] for row in result["primary"]]
    normalized_values = [row["normalized_score_variance_median"] for row in result["primary"]]
    x = np.arange(len(methods)); width = 0.36
    fig, axis = plt.subplots(figsize=(8, 5), layout="constrained")
    axis.bar(x - width / 2, raw_values, width, label="planet-masked raw")
    axis.bar(x + width / 2, normalized_values, width, label="planet-masked radial standardized")
    axis.axhline(1, color="black", linestyle="--", linewidth=1)
    axis.set(title="11-pixel response, controlling radii", ylabel="median candidate score variance",
             xticks=x, xticklabels=methods)
    axis.tick_params(axis="x", rotation=30, labelsize=8)
    axis.grid(axis="y", alpha=0.2); axis.legend()
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, masked_root: Path, unmasked_root: Path, root: Path) -> None:
    """Run the corrected masked and radially standardized 11-pixel comparison."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = raw.verify_bundle(bundle_path)
    for parent in (masked_root, unmasked_root):
        completion = read(parent / "complete.json")
        stage.verify([value for key, value in completion.items() if key != "status"])
    masked_result = read(masked_root / "results.json")
    root.mkdir(parents=True)
    protocol = {"schema": 1, "purpose": "corrected strict radial normalization with fixed candidate planet mask",
                "mode": 200, "support": 11, "radii": list(raw.RADII),
                "training_bands": "narrowest split-supported width only",
                "radial_profile": "leave-site-out 3.6-pixel bins from the verified unmasked-normalization contract",
                "planet_mask": "fixed seven-pixel native disk removed from candidate data and exact response before principal-submatrix solve",
                "bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
                "masked_raw_parent": stage.fingerprint(masked_root / "complete.json"),
                "unmasked_normalized_parent": stage.fingerprint(unmasked_root / "complete.json"),
                "scripts": [stage.fingerprint(Path(__file__)),
                            stage.fingerprint(Path(__file__).with_name("run_klip_stage_b_local_planet_masked_screen.py")),
                            stage.fingerprint(Path(__file__).with_name("run_klip_stage_b_local_radial_normalization.py"))]}
    stage.write_json(root / "protocol.json", protocol)
    stage.write_json(root / "state.json", {"status": "running", "started_unix": time.time()})

    baseline = np.asarray(arrays["baseline"], dtype=np.float64)
    finite = np.isfinite(baseline)
    planet = raw.planet_position(baseline.shape, metadata)
    lookup = {(int(row), int(column)): index for index, (row, column) in enumerate(arrays["positions"])}
    records = []
    started = time.monotonic()
    for radius in raw.RADII:
        fixed = metadata["geometry"]["11"][str(radius)]
        width = int(fixed["training"]["narrowest_split_supported_half_width"])
        for site_index, site in enumerate(fixed["selected_sites"]):
            searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                        for delta_row, delta_column in footprint.SEARCH_OFFSETS]
            excluded = normalized.exclusion_mask(baseline.shape, searches, 11, planet)
            scale, _ = normalized.variance_profile(baseline, excluded)
            standardized = baseline / scale
            for search_index, query in enumerate(searches):
                rings = raw.training_stencils(finite, query, searches, planet, 11)
                source = lookup[query]
                template = raw.crop(np.asarray(arrays["response"][source], dtype=np.float64).T, 11)
                validity = raw.crop(np.asarray(arrays["validity"][source], dtype=bool).T, 11)
                planet_mask = masked.candidate_mask(query, 11, planet) & validity
                row, column = query
                candidate = baseline[column - 5:column + 6, row - 5:row + 6].T
                candidate_scale = scale[column - 5:column + 6, row - 5:row + 6].T
                standardized_template = template / candidate_scale
                standardized_candidate = candidate / candidate_scale
                halves = []
                for detector_half in (0, 1):
                    samples, interpolation_gain = raw.band_samples(standardized, rings, width, detector_half)
                    opposite, _ = raw.band_samples(standardized, rings, width, 1 - detector_half)
                    rectangular = raw.fit_periodogram_base(samples, 11, "rectangular")
                    hann = raw.fit_periodogram_base(samples, 11, "hann")
                    models = {"identity": raw.regularize(rectangular, 1.0),
                              "rectangular_m0.1": raw.regularize(rectangular, 0.1),
                              "rectangular_m0.3": raw.regularize(rectangular, 0.3),
                              "hann_m0.1": raw.regularize(hann, 0.1),
                              "hann_m0.3": raw.regularize(hann, 0.3)}
                    results, details = {}, {}
                    for method, _, _ in raw.VARIANTS:
                        result, detail = masked.analyze_model(models[method], standardized_template,
                                                             standardized_candidate, planet_mask, opposite)
                        results[method] = result
                        details[method] = detail | {"physical_weight": detail["weight"] /
                                                                     candidate_scale.ravel()}
                    halves.append({"samples": len(samples), "results": results, "details": details})
                for method, _, _ in raw.VARIANTS:
                    first = halves[0]["details"][method]; second = halves[1]["details"][method]
                    cosine = float(first["physical_weight"] @ second["physical_weight"] /
                                   (np.linalg.norm(first["physical_weight"]) *
                                    np.linalg.norm(second["physical_weight"])))
                    records.append({"radius": radius, "site_index": site_index,
                                    "search_index": search_index, "row": row, "column": column,
                                    "band_half_width": width, "method": method,
                                    "samples": [halves[index]["samples"] for index in (0, 1)],
                                    "physical_weight_cosine": cosine,
                                    "valid_pixels": int(np.count_nonzero(planet_mask)),
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

    summarized = normalized.summarize(records, masked_result)
    result = {"purpose": protocol["purpose"], "mode": 200, "support": 11,
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
    stage.write_json(root / "state.json", {"status": "complete", "elapsed_seconds": result["elapsed_seconds"]})
    print(root / "results.md", flush=True)


def check() -> None:
    """Check simultaneous scaling and masking preserve physical amplitude."""
    generator = np.random.default_rng(924)
    samples = generator.normal(size=(32, 121))
    model = raw.regularize(raw.fit_periodogram_base(samples, 11, "hann"), 0.1)
    template = generator.normal(size=(11, 11))
    scale = np.exp(generator.normal(scale=0.2, size=(11, 11)))
    mask = generator.random((11, 11)) > 0.2
    amplitude = 0.041
    result, _ = masked.analyze_model(model, template / scale, amplitude * template / scale,
                                     mask, samples)
    stage.require(np.isclose(result["amplitude_no_mean"], amplitude, rtol=2e-12, atol=2e-14),
                  "combined scaling and masking changed physical amplitude")
    print("KLIP Stage-B masked radial-normalization checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the corrected normalization parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("masked", type=Path)
    run_parser.add_argument("unmasked", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the corrected normalization action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        run(args.bundle.resolve(), args.masked.resolve(), args.unmasked.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
