#!/usr/bin/env python3
"""Audit 11-pixel PSD training behind larger KLIP response footprints."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_response_47 as response47  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402


TRAINING_SUPPORT = 11
RESPONSE_SUPPORTS = (11, 31, 47)
RADII = (6.0, 7.5, 10.0, 12.0, 16.0, 20.0, 24.0)


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def verify_parent(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Verify the completed coupled-footprint preflight."""
    protocol = footprint.load_experiment(root)
    completion = read(root / "complete.json")
    state = read(root / "state.json")
    stage.require(completion["status"] == "complete" and state["status"] == "complete",
                  "coupled-footprint preflight is incomplete")
    stage.verify([completion["results"], completion["report"]])
    results = read(root / "results.json")
    stage.require(not results["selection_uses_baseline_values"] and
                  results["baseline_used_only_for_common_finite_mask"],
                  "parent preflight was not score blind")
    return protocol, results


def prepare(args: argparse.Namespace) -> None:
    """Freeze the parent preflight, decoupled support grid, and software."""
    root, parent = args.root.resolve(), args.footprint.resolve()
    stage.require(not root.exists(), f"output already exists: {root}")
    parent_protocol, _ = verify_parent(parent)
    resources = parent_protocol["resources"]
    stage.require(sorted(os.sched_getaffinity(0)) == resources["cpu_affinity"],
                  "prepare must use the parent CPU affinity")

    root.mkdir(parents=True)
    (root / "software").mkdir()
    runner = Path(__file__).resolve()
    dependencies = [runner, runner.with_name("run_klip_stage_b_footprint_preflight.py"),
                    runner.with_name("run_klip_response_47.py"),
                    runner.with_name("run_klip_response_tail_linearity.py"),
                    runner.with_name("run_klip_response_stamp_convergence.py"),
                    runner.with_name("run_klip_covariance_stage_a.py")]
    for path in dependencies:
        shutil.copy2(path, root / "software" / path.name)

    protocol = {
        "schema": 1,
        "stage": "KLIP covariance matched-filter Stage B decoupled-footprint preflight",
        "purpose": "test 11-pixel Welch training behind 11-, 31-, and 47-pixel response exclusions",
        "parent_footprint_preflight": str(parent),
        "parent_response": parent_protocol["parent_response"],
        "training_support": TRAINING_SUPPORT,
        "response_exclusion_supports": list(RESPONSE_SUPPORTS),
        "radii": list(RADII),
        "training_band_half_widths": list(footprint.BAND_WIDTHS[TRAINING_SUPPORT]),
        "training_center_step": (TRAINING_SUPPORT - 1) // 2,
        "minimum_patches_per_disjoint_half": footprint.MINIMUM_SPLIT_PATCHES,
        "candidate_sites": "reuse each support's score-blind deterministic sites from the parent preflight",
        "exclusion": "union of five response-sized candidate footprints plus the known planet",
        "resources": resources,
    }
    stage.write_json(root / "protocol.json", protocol)
    frozen = [stage.fingerprint(root / "protocol.json")]
    frozen.extend(stage.fingerprint(root / "software" / path.name) for path in dependencies)
    parent_records = [stage.fingerprint(parent / name)
                      for name in ("protocol.json", "manifest.json", "results.json", "results.md",
                                   "complete.json")]
    stage.write_json(root / "manifest.json", {"schema": 1, "frozen_records": frozen,
                                              "parent_records": parent_records})
    stage.write_json(root / "state.json", {"status": "prepared"})
    print(root, flush=True)


def load_experiment(root: Path) -> dict[str, object]:
    """Load one prepared experiment and verify its lineage."""
    protocol = read(root / "protocol.json")
    manifest = read(root / "manifest.json")
    stage.verify(manifest["frozen_records"] + manifest["parent_records"])
    parent_protocol, _ = verify_parent(Path(str(protocol["parent_footprint_preflight"])))
    stage.require(parent_protocol["resources"] == protocol["resources"], "parent resource contract changed")
    stage.require(sorted(os.sched_getaffinity(0)) == protocol["resources"]["cpu_affinity"],
                  "CPU affinity changed from the frozen protocol")
    return protocol


def analyze(root: Path, protocol: dict[str, object]) -> dict[str, object]:
    """Audit 11-pixel training geometry for all response-sized exclusions."""
    parent_root = Path(str(protocol["parent_footprint_preflight"]))
    parent = read(parent_root / "results.json")
    response_root = Path(str(protocol["parent_response"]))
    paths = response47.product_paths(response_root)
    baseline = np.asarray(fits.getdata(paths["final"]), dtype=np.float64)
    finite = np.all(np.isfinite(baseline), axis=0)
    planet = tuple(float(value) for value in parent["known_planet_row_column"])

    geometry = {}
    for response_support in RESPONSE_SUPPORTS:
        by_radius = {}
        for radius in RADII:
            parent_row = parent["geometry"][str(response_support)][str(radius)]
            sites = parent_row["selected_sites"]
            query_records = []
            for site_index, site in enumerate(sites):
                searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                            for delta_row, delta_column in footprint.SEARCH_OFFSETS]
                for search_index, query in enumerate(searches):
                    rings = footprint.training_rings(
                        finite, query, searches, TRAINING_SUPPORT, planet,
                        exclusion_support=response_support)
                    query_records.append({"site_index": site_index, "search_index": search_index,
                                          "row": query[0], "column": query[1],
                                          "bands": footprint.band_counts(rings, TRAINING_SUPPORT)})
            by_radius[str(radius)] = {
                "eligible_five_pixel_centers": parent_row["eligible_five_pixel_centers"],
                "selected_sites": sites,
                "training": footprint.geometry_summary(query_records, TRAINING_SUPPORT)
                if query_records else None,
            }
        geometry[str(response_support)] = by_radius

    result = {
        "purpose": protocol["purpose"],
        "training_support": TRAINING_SUPPORT,
        "geometry": geometry,
        "selection_uses_baseline_values": False,
        "baseline_used_only_for_common_finite_mask": True,
        "parent_response_hash": stage.fingerprint(response_root / "complete.json")["sha256"],
        "response_protocol_hash": stage.fingerprint(response_root / "protocol.json")["sha256"],
    }
    stage.write_json(root / "results.json", result)
    write_report(root, result)
    return result


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write the decoupled-footprint coverage tables."""
    lines = ["# KLIP Stage-B decoupled-footprint preflight", "",
             "All rows use 11-by-11 Welch training patches with five-pixel center spacing. Response support "
             "controls the candidate data and source-exclusion footprint.", "",
             "## Narrowest band with eight training patches in each detector half", "",
             "| Radius | 11-pixel response | 31-pixel response | 47-pixel response |",
             "| ---: | ---: | ---: | ---: |"]
    for radius in RADII:
        values = []
        for support in RESPONSE_SUPPORTS:
            training = result["geometry"][str(support)][str(radius)]["training"]
            width = training["narrowest_split_supported_half_width"] if training else None
            values.append("—" if width is None else str(width))
        lines.append(f"| {radius:g} | {values[0]} | {values[1]} | {values[2]} |")
    lines.extend(["", "## Full-range minimum counts", "",
                  "Entries give the minimum across all selected sites and five search pixels.", "",
                  "| Radius | Response support | Accepted | First half | Second half | Rank ceiling |",
                  "| ---: | ---: | ---: | ---: | ---: | ---: |"])
    widest = str(footprint.BAND_WIDTHS[TRAINING_SUPPORT][-1])
    for radius in RADII:
        for support in RESPONSE_SUPPORTS:
            training = result["geometry"][str(support)][str(radius)]["training"]
            if training is None:
                values = ["—"] * 4
            else:
                row = training["band_half_widths"][widest]
                values = [f"{row[key]['minimum']:g}"
                          for key in ("accepted", "top", "bottom", "centered_rank_ceiling")]
            lines.append(f"| {radius:g} | {support} | {values[0]} | {values[1]} | {values[2]} | {values[3]} |")
    lines.extend(["", "The calculation uses baseline finiteness only. It does not fit a PSD, covariance, or "
                  "matched-filter score."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    """Run the decoupled preflight and write completion receipts."""
    root = args.root.resolve()
    protocol = load_experiment(root)
    if (root / "complete.json").is_file():
        completion = read(root / "complete.json")
        stage.verify([completion["results"], completion["report"]])
        print(root / "results.md", flush=True)
        return
    stage.write_json(root / "state.json", {"status": "running", "pid": os.getpid()})
    try:
        analyze(root, protocol)
        completion = {"status": "complete", "results": stage.fingerprint(root / "results.json"),
                      "report": stage.fingerprint(root / "results.md")}
        stage.write_json(root / "complete.json", completion)
        stage.write_json(root / "state.json", {"status": "complete"})
        print(root / "results.md", flush=True)
    except Exception as error:
        stage.write_json(root / "state.json", {"status": "failed", "error_type": type(error).__name__,
                                               "error": str(error)})
        raise


def check() -> None:
    """Check that larger exclusions monotonically reduce training coverage."""
    finite = np.ones((128, 128), dtype=bool)
    searches = [(76, 64), (75, 64), (77, 64), (76, 63), (76, 65)]
    totals = []
    for exclusion in RESPONSE_SUPPORTS:
        rings = footprint.training_rings(finite, searches[0], searches, TRAINING_SUPPORT,
                                         (51.0, 61.0), exclusion_support=exclusion)
        totals.append(footprint.band_counts(rings, TRAINING_SUPPORT)["60"]["accepted"])
    stage.require(totals[0] >= totals[1] >= totals[2], "training coverage is not monotonic with exclusion size")
    stage.require(footprint.BAND_WIDTHS[TRAINING_SUPPORT][-1] == 60, "full training band changed")
    print("KLIP Stage-B decoupled-footprint checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the decoupled-preflight command-line parser."""
    repo = Path(__file__).resolve().parents[3]
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--footprint", type=Path,
                                default=repo / "working/roc/klip_stage_b_footprint_preflight_20260923")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    return result


def main() -> None:
    """Dispatch the requested action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    elif args.action == "prepare":
        prepare(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
