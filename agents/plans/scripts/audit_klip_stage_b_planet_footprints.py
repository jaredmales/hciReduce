#!/usr/bin/env python3
"""Audit known-planet overlap in the fixed KLIP Stage-B candidate footprints."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402
import run_klip_stage_b_local_noise_screen as raw  # noqa: E402


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def footprint_clear(query: tuple[int, int], support: int,
                    planet: tuple[float, float]) -> bool:
    """Return whether every native candidate pixel avoids the planet exclusion disk."""
    half = support // 2
    row, column = query
    native_row, native_column = np.mgrid[row - half:row + half + 1,
                                         column - half:column + half + 1]
    return not np.any(np.hypot(native_row - planet[0], native_column - planet[1]) <=
                      stage.PLANET_EXCLUSION_RADIUS)


def geometry(metadata: dict[str, object], planet: tuple[float, float]) -> dict[str, object]:
    """Count clear queries and complete five-pixel sites for each fixed geometry."""
    result = {}
    for support in raw.SUPPORTS:
        by_radius = {}
        for radius in raw.RADII:
            records = []
            for site_index, site in enumerate(metadata["geometry"][str(support)][str(radius)]["selected_sites"]):
                clear = []
                for search_index, (delta_row, delta_column) in enumerate(footprint.SEARCH_OFFSETS):
                    query = (int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                    value = footprint_clear(query, support, planet)
                    clear.append(value)
                    records.append({"site_index": site_index, "search_index": search_index,
                                    "row": query[0], "column": query[1], "clear": value})
            by_radius[str(radius)] = {"queries": len(records),
                                      "clear_queries": sum(record["clear"] for record in records),
                                      "sites": len(records) // len(footprint.SEARCH_OFFSETS),
                                      "clear_sites": sum(all(record["clear"] for record in records
                                                             if record["site_index"] == site_index)
                                                         for site_index in range(len(records) // 5)),
                                      "records": records}
        result[str(support)] = by_radius
    return result


def clear_site_lookup(audit: dict[str, object], support: int, radius: float) -> set[int]:
    """Return sites whose complete five-pixel search footprints avoid the planet."""
    row = audit[str(support)][str(radius)]
    return {site_index for site_index in range(row["sites"])
            if all(record["clear"] for record in row["records"]
                   if record["site_index"] == site_index)}


def score_stats(records: list[dict[str, object]]) -> dict[str, object]:
    """Summarize the two directional scores from one geometry-only subset."""
    values = np.asarray([value for record in records for value in record["score_no_mean"]])
    return raw.score_stats(values)


def score_audit(records: list[dict[str, object]], audit: dict[str, object],
                normalized: bool) -> dict[str, object]:
    """Recompute 11-pixel narrow-band score summaries on footprint-clear sites."""
    result = {}
    for radius in raw.RADII:
        clear_sites = clear_site_lookup(audit, 11, radius)
        for method, _, _ in raw.VARIANTS:
            selected = [record for record in records if record["radius"] == radius and
                        record["method"] == method and record["site_index"] in clear_sites and
                        (normalized or (record["support"] == 11 and record["band_role"] == "narrow"))]
            stage.require(len(selected) == 5 * len(clear_sites), "clear-site record count changed")
            result[f"r{radius:g}_{method}"] = {"radius": radius, "method": method,
                                               "clear_sites": len(clear_sites),
                                               "query_records": len(selected),
                                               "scores": score_stats(selected)}
    primary = []
    for method, _, _ in raw.VARIANTS:
        selected = [result[f"r{radius:g}_{method}"] for radius in raw.PRIMARY_RADII]
        primary.append({"method": method,
                        "median_clear_site_variance": float(np.median([
                            row["scores"]["variance"] for row in selected])),
                        "radius_variances": {str(row["radius"]): row["scores"]["variance"]
                                             for row in selected},
                        "clear_sites": {str(row["radius"]): row["clear_sites"] for row in selected}})
    return {"groups": result, "primary": primary}


def residual_diagnostics(image: np.ndarray, planet: tuple[float, float]) -> dict[str, object]:
    """Compare the subtracted-planet neighborhood with the clear radius-12 annulus."""
    native_column, native_row = np.indices(image.shape)
    radius = np.hypot(native_row - 63.5, native_column - 63.5)
    distance = np.hypot(native_row - planet[0], native_column - planet[1])
    annulus = image[(radius >= 11.5) & (radius < 13.5) &
                    (distance > stage.PLANET_EXCLUSION_RADIUS)]
    result = {"comparison_annulus": {"pixels": len(annulus),
                                      "rms": float(np.sqrt(np.mean(np.square(annulus)))),
                                      "maximum_absolute": float(np.max(np.abs(annulus)))}, "disks": {}}
    for disk_radius in (2.0, 3.6, 5.0, 7.0):
        pixels = image[distance <= disk_radius]
        result["disks"][str(disk_radius)] = {
            "pixels": len(pixels), "rms": float(np.sqrt(np.mean(np.square(pixels)))),
            "maximum_absolute": float(np.max(np.abs(pixels)))}
    return result


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the overlap audit and corrected interpretation."""
    lines = ["# KLIP Stage-B known-planet footprint audit", "",
             "The earlier preflight rejected candidate centers within seven pixels of the fitted planet and "
             "excluded that disk from covariance training and radial-profile estimation. It did not reject or "
             "mask candidate stamps whose outer pixels enter the disk.", "", "## Complete-footprint geometry", "",
             "Entries are clear five-pixel sites; parenthetical values are clear individual query footprints.", "",
             "| Radius | 11-pixel response | 31-pixel response | 47-pixel response |",
             "| ---: | ---: | ---: | ---: |"]
    for radius in raw.RADII:
        values = []
        for support in raw.SUPPORTS:
            row = result["geometry"][str(support)][str(radius)]
            values.append(f"{row['clear_sites']}/{row['sites']} ({row['clear_queries']}/{row['queries']})")
        lines.append(f"| {radius:g} | {values[0]} | {values[1]} | {values[2]} |")
    lines.extend(["", "## Footprint-clear 11-pixel scores", "",
                  "These geometry-only subsets retain all five queries at each selected site.", "",
                  "| Method | Raw all-site variance | Raw clear-site variance | Normalized all-site variance | Normalized clear-site variance |",
                  "| :--- | ---: | ---: | ---: | ---: |"])
    raw_all = {row["method"]: row for row in result["parent_primary"]["raw"]}
    normalized_all = {row["method"]: row for row in result["parent_primary"]["normalized"]}
    raw_clear = {row["method"]: row for row in result["clear_scores"]["raw"]["primary"]}
    normalized_clear = {row["method"]: row for row in result["clear_scores"]["normalized"]["primary"]}
    for method, _, _ in raw.VARIANTS:
        lines.append(f"| {method} | {raw_all[method]['score_variance_median']:.3f} | "
                     f"{raw_clear[method]['median_clear_site_variance']:.3f} | "
                     f"{normalized_all[method]['normalized_score_variance_median']:.3f} | "
                     f"{normalized_clear[method]['median_clear_site_variance']:.3f} |")
    hann_raw = raw_clear["hann_m0.1"]["radius_variances"]
    hann_normalized = normalized_clear["hann_m0.1"]["radius_variances"]
    lines.extend(["", "The clear-site restriction does not explain away the radius-12 mismatch. For "
                  f"Hann/mixing-0.1, radius-12 variance is {hann_raw['12.0']:.3f} raw and "
                  f"{hann_normalized['12.0']:.3f} after radial normalization, both larger than the all-site "
                  "values. The subtraction residual's seven-pixel RMS is comparable to the remaining radius-12 "
                  "annulus, but the conservative exclusion contract still requires candidate data to avoid it.", "",
                  "## Decision", "",
                  "Treat the earlier 31- and 47-pixel candidate-score comparisons as invalid under the intended "
                  "known-planet exclusion. The 11-pixel qualitative conclusions survive on the geometry-only "
                  "clear subset, with fewer and non-reselected angular sites. Before further Stage-B tuning, add "
                  "a fixed native planet mask to candidate data and exact responses and solve the corresponding "
                  "positive covariance submatrix. Then rerun the raw support comparison on the original frozen "
                  "sites. Do not use the actual planet or these baseline scores to choose replacement sites."])
    (root / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    methods = [item[0] for item in raw.VARIANTS]
    x = np.arange(len(methods))
    raw_all_values = [raw_all[name]["score_variance_median"] for name in methods]
    raw_clear_values = [raw_clear[name]["median_clear_site_variance"] for name in methods]
    norm_all_values = [normalized_all[name]["normalized_score_variance_median"] for name in methods]
    norm_clear_values = [normalized_clear[name]["median_clear_site_variance"] for name in methods]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), layout="constrained")
    width = 0.36
    for axis, all_values, clear_values, title in ((axes[0], raw_all_values, raw_clear_values, "Raw"),
                                                  (axes[1], norm_all_values, norm_clear_values, "Radial standardized")):
        axis.bar(x - width / 2, all_values, width, label="all fixed sites")
        axis.bar(x + width / 2, clear_values, width, label="complete-footprint clear")
        axis.axhline(1, color="black", linestyle="--", linewidth=1)
        axis.set(title=title, ylabel="primary median score variance", xticks=x, xticklabels=methods)
        axis.tick_params(axis="x", rotation=30, labelsize=8)
        axis.grid(axis="y", alpha=0.2)
        axis.legend()
    fig.suptitle("11-pixel response: effect of complete known-planet footprint exclusion")
    fig.savefig(root / "comparison.png", dpi=170)
    plt.close(fig)


def run(bundle_path: Path, raw_root: Path, normalized_root: Path, root: Path) -> None:
    """Run the footprint audit from verified local products."""
    stage.require(not root.exists(), f"output already exists: {root}")
    metadata, arrays, receipt = raw.verify_bundle(bundle_path)
    for parent in (raw_root, normalized_root):
        completion = read(parent / "complete.json")
        stage.verify([value for key, value in completion.items() if key != "status"])
    planet = raw.planet_position(arrays["baseline"].shape, metadata)
    footprint_geometry = geometry(metadata, planet)
    raw_records = read(raw_root / "records.json")
    normalized_records = read(normalized_root / "records.json")
    raw_result = read(raw_root / "results.json")
    normalized_result = read(normalized_root / "results.json")
    result = {"purpose": "audit complete candidate-footprint avoidance of the fitted KLIP planet",
              "mode": 200, "planet_row_column": list(planet),
              "planet_exclusion_radius": stage.PLANET_EXCLUSION_RADIUS,
              "geometry": footprint_geometry,
              "residual": residual_diagnostics(np.asarray(arrays["baseline"], dtype=np.float64), planet),
              "parent_primary": {"raw": [row for row in raw_result["primary"] if row["support"] == 11 and
                                           row["band_role"] == "narrow"],
                                 "normalized": normalized_result["primary"]},
              "clear_scores": {"raw": score_audit(raw_records, footprint_geometry, False),
                               "normalized": score_audit(normalized_records, footprint_geometry, True)},
              "inputs": {"bundle": stage.fingerprint(bundle_path), "bundle_receipt": receipt,
                         "raw": stage.fingerprint(raw_root / "complete.json"),
                         "normalized": stage.fingerprint(normalized_root / "complete.json")},
              "script": stage.fingerprint(Path(__file__))}
    root.mkdir(parents=True)
    stage.write_json(root / "results.json", result)
    write_report(root, result)
    completion = {"status": "complete", "results": stage.fingerprint(root / "results.json"),
                  "report": stage.fingerprint(root / "README.md"),
                  "figure": stage.fingerprint(root / "comparison.png")}
    stage.write_json(root / "complete.json", completion)
    print(root / "README.md", flush=True)


def check() -> None:
    """Check exact native-pixel footprint intersection behavior."""
    planet = (75.72, 61.49)
    stage.require(not footprint_clear((76, 69), 11, planet) and
                  footprint_clear((76, 80), 11, planet) and
                  not footprint_clear((64, 64), 47, planet),
                  "planet-footprint intersection check failed")
    print("KLIP Stage-B planet-footprint audit checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the footprint-audit command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("bundle", type=Path)
    run_parser.add_argument("raw", type=Path)
    run_parser.add_argument("normalized", type=Path)
    run_parser.add_argument("output", type=Path)
    return result


def main() -> None:
    """Dispatch the requested footprint audit."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        run(args.bundle.resolve(), args.raw.resolve(), args.normalized.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
