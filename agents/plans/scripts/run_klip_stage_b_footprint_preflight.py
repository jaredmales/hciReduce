#!/usr/bin/env python3
"""Audit KLIP response footprints and covariance-training geometry."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_response_47 as response47  # noqa: E402


SUPPORTS = (11, 31, 47)
RADII = (6.0, 7.5, 10.0, 12.0, 16.0, 20.0, 24.0)
PRIMARY_RADII = RADII[1:]
SEARCH_OFFSETS = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))
SITES_PER_RADIUS = 12
BAND_WIDTHS = {11: (0, 5, 10, 20, 40, 60),
               31: (0, 15, 30, 45, 60),
               47: (0, 23, 46, 69)}
MINIMUM_SPLIT_PATCHES = 8


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def verify_parent(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Verify the completed 47-pixel response campaign and its products."""
    protocol = response47.load_experiment(root)
    completion = read(root / "complete.json")
    state = read(root / "state.json")
    stage.require(completion["status"] == "complete" and state["status"] == "complete",
                  "47-pixel response campaign is incomplete")
    stage.verify([completion["response"], completion["results"], completion["report"]])
    response_receipt = read(root / "response" / "complete.json")
    stage.verify(response_receipt["products"])
    results = read(root / "results.json")
    stage.require(response_receipt["validation"]["baseline_bitwise_identical"] and
                  results["central_11_archive_replay_passed"] and
                  results["external_47_replay_passed"], "47-pixel response campaign failed a gate")
    return protocol, completion


def prepare(args: argparse.Namespace) -> None:
    """Freeze the exact-response lineage, geometry policy, and software."""
    root, parent = args.root.resolve(), args.response.resolve()
    stage.require(not root.exists(), f"output already exists: {root}")
    parent_protocol, _ = verify_parent(parent)
    resources = parent_protocol["resources"]
    stage.require(sorted(os.sched_getaffinity(0)) == resources["cpu_affinity"],
                  "prepare must use the response campaign CPU affinity")

    root.mkdir(parents=True)
    (root / "software").mkdir()
    runner = Path(__file__).resolve()
    dependencies = [runner, runner.with_name("run_klip_response_47.py"),
                    runner.with_name("run_klip_response_tail_linearity.py"),
                    runner.with_name("run_klip_response_stamp_convergence.py"),
                    runner.with_name("run_klip_covariance_stage_a.py")]
    for path in dependencies:
        shutil.copy2(path, root / "software" / path.name)

    protocol = {
        "schema": 1,
        "stage": "KLIP covariance matched-filter Stage B footprint preflight",
        "purpose": "freeze response support and covariance-training coverage before reading baseline scores",
        "parent_response": str(parent),
        "modes": stage.MODES,
        "primary_mode": 200,
        "supports": list(SUPPORTS),
        "radii": list(RADII),
        "primary_radii": list(PRIMARY_RADII),
        "radius_bin_half_width": 0.5,
        "search_offsets_row_column": [list(offset) for offset in SEARCH_OFFSETS],
        "sites_per_radius": SITES_PER_RADIUS,
        "training_geometry": {
            "angular_and_radial_step": "(support-1)/2 pixels",
            "band_half_widths": {str(key): list(value) for key, value in BAND_WIDTHS.items()},
            "minimum_patches_per_disjoint_half": MINIMUM_SPLIT_PATCHES,
            "radial_center_range": [6.0, 60.0],
            "orientation": "radial/tangential bilinear extraction",
            "exclusions": "known planet and union of five complete candidate footprints",
        },
        "known_planet": {
            "separation": stage.PLANET_SEPARATION,
            "position_angle": stage.PLANET_PA,
            "exclusion_radius": stage.PLANET_EXCLUSION_RADIUS,
        },
        "resources": resources,
    }
    stage.write_json(root / "protocol.json", protocol)
    frozen = [stage.fingerprint(root / "protocol.json")]
    frozen.extend(stage.fingerprint(root / "software" / path.name) for path in dependencies)
    parent_records = [stage.fingerprint(parent / name)
                      for name in ("protocol.json", "manifest.json", "results.json", "results.md",
                                   "complete.json", "response/complete.json")]
    response_receipt = read(parent / "response" / "complete.json")
    stage.write_json(root / "manifest.json", {"schema": 1, "frozen_records": frozen,
                                              "parent_records": parent_records,
                                              "response_products": response_receipt["products"]})
    stage.write_json(root / "state.json", {"status": "prepared"})
    print(root, flush=True)


def load_experiment(root: Path) -> dict[str, object]:
    """Load the preflight and verify frozen software and response products."""
    protocol = read(root / "protocol.json")
    manifest = read(root / "manifest.json")
    stage.verify(manifest["frozen_records"] + manifest["parent_records"] + manifest["response_products"])
    parent_protocol, _ = verify_parent(Path(str(protocol["parent_response"])))
    stage.require(parent_protocol["resources"] == protocol["resources"], "parent resource contract changed")
    stage.require(sorted(os.sched_getaffinity(0)) == protocol["resources"]["cpu_affinity"],
                  "CPU affinity changed from the frozen protocol")
    return protocol


def summarize(values: np.ndarray) -> dict[str, object]:
    """Return fixed quantiles for one finite vector."""
    array = np.asarray(values, dtype=np.float64)
    array = array[np.isfinite(array)]
    if not len(array):
        return {"count": 0, "minimum": None, "p10": None, "median": None, "p90": None, "maximum": None}
    return {"count": len(array), "minimum": float(np.min(array)),
            "p10": float(np.percentile(array, 10)), "median": float(np.median(array)),
            "p90": float(np.percentile(array, 90)), "maximum": float(np.max(array))}


def response_summary(response: np.ndarray, validity: np.ndarray, radii: np.ndarray,
                     support: int, full_energy: np.ndarray) -> tuple[np.ndarray, dict[str, object]]:
    """Summarize completeness and energy captured by one central support."""
    center = response.shape[1] // 2
    half = support // 2
    selection = (slice(None), slice(center - half, center + half + 1),
                 slice(center - half, center + half + 1))
    local_response = response[selection]
    local_validity = validity[selection]
    safe = np.where(local_validity, local_response, 0)
    energy = np.einsum("ijk,ijk->i", safe, safe, dtype=np.float64)
    fraction = energy / full_energy
    complete = np.all(local_validity, axis=(1, 2))
    border = np.zeros((support, support), dtype=bool)
    border[[0, -1], :] = True
    border[:, [0, -1]] = True
    border_energy = np.einsum("ijk,ijk->i", safe * border, safe * border, dtype=np.float64) / energy
    bins = {}
    for nominal in RADII:
        chosen = np.abs(radii - nominal) <= 0.5
        bins[str(nominal)] = {
            "responses": int(np.count_nonzero(chosen)),
            "complete_responses": int(np.count_nonzero(chosen & complete)),
            "energy_fraction_of_47": summarize(fraction[chosen]),
            "border_energy_fraction": summarize(border_energy[chosen]),
        }
    return complete, {"complete_response_count": int(np.count_nonzero(complete)), "bins": bins}


def complete_data_support(baseline: np.ndarray, coordinates: np.ndarray, support: int) -> np.ndarray:
    """Mark exact locations whose data stamp is finite in every mode."""
    half = support // 2
    height, width = baseline.shape[1:]
    finite = np.all(np.isfinite(baseline), axis=0)
    result = np.zeros(len(coordinates), dtype=bool)
    for index, (row_value, column_value, _, _) in enumerate(coordinates):
        row, column = int(row_value), int(column_value)
        if row - half < 0 or row + half >= width or column - half < 0 or column + half >= height:
            continue
        result[index] = np.all(finite[column - half:column + half + 1,
                                     row - half:row + half + 1])
    return result


def eligible_centers(coordinates: np.ndarray, radii: np.ndarray, complete: np.ndarray,
                     data_complete: np.ndarray, nominal: float,
                     planet_row: float, planet_column: float) -> list[dict[str, object]]:
    """Find common-mode five-pixel search centers using geometry only."""
    lookup = {(int(row), int(column)): index for index, (row, column, _, _) in enumerate(coordinates)}
    centers = []
    for index, (row_value, column_value, _, _) in enumerate(coordinates):
        row, column = int(row_value), int(column_value)
        if abs(radii[index] - nominal) > 0.5:
            continue
        neighbors = [lookup.get((row + delta_row, column + delta_column))
                     for delta_row, delta_column in SEARCH_OFFSETS]
        if any(source is None or not complete[int(source)] or not data_complete[int(source)]
               for source in neighbors):
            continue
        if any(math.hypot(row + delta_row - planet_row, column + delta_column - planet_column) <=
                   stage.PLANET_EXCLUSION_RADIUS for delta_row, delta_column in SEARCH_OFFSETS):
            continue
        centers.append({"row": row, "column": column,
                        "angle_radians": math.atan2(column - 63.5, row - 63.5)})
    return centers


def select_sites(centers: list[dict[str, object]]) -> list[dict[str, object]]:
    """Choose up to twelve deterministic angularly distributed centers."""
    if not centers:
        return []
    ordered = sorted(centers, key=lambda item: (float(item["angle_radians"]),
                                                int(item["row"]), int(item["column"])))
    count = min(SITES_PER_RADIUS, len(ordered))
    indices = np.floor((np.arange(count) + 0.5) * len(ordered) / count).astype(int)
    return [ordered[int(index)] for index in indices]


def training_rings(finite: np.ndarray, query: tuple[int, int], searches: list[tuple[int, int]],
                   support: int, planet: tuple[float, float]) -> dict[int, dict[str, int]]:
    """Count accepted radial/tangential patches on half-overlap rings."""
    row_query, column_query = query
    center_row = 0.5 * (finite.shape[1] - 1)
    center_column = 0.5 * (finite.shape[0] - 1)
    radius = math.hypot(row_query - center_row, column_query - center_column)
    angle = math.atan2(column_query - center_column, row_query - center_row)
    half = support // 2
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
        sample_row = (centers_row[:, None] + cosine * delta_row.ravel() -
                      sine * delta_column.ravel())
        sample_column = (centers_column[:, None] + sine * delta_row.ravel() +
                         cosine * delta_column.ravel())
        row0 = np.floor(sample_row).astype(int)
        column0 = np.floor(sample_column).astype(int)
        row_fraction = sample_row - row0
        column_fraction = sample_column - column0
        excluded = np.zeros(count, dtype=bool)
        incomplete = np.zeros(count, dtype=bool)
        used_columns = []
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
                    bad |= ((np.abs(native_row - search_row) <= half) &
                            (np.abs(native_column - search_column) <= half))
                excluded |= np.any(used & bad & ~outside, axis=1)
                incomplete |= np.any(used & (outside | ~finite[safe_column, safe_row]), axis=1)
                used_columns.append(np.where(used, native_column, np.nan))
        accepted = ~(excluded | incomplete)
        columns = np.stack(used_columns, axis=-1)[accepted]
        if len(columns):
            halves = np.where(np.nanmax(columns, axis=(1, 2)) < center_column, 0,
                              np.where(np.nanmin(columns, axis=(1, 2)) > center_column, 1, -1))
        else:
            halves = np.empty(0, dtype=int)
        rings[offset] = {"attempted": count, "accepted": int(np.count_nonzero(accepted)),
                         "top": int(np.count_nonzero(halves == 0)),
                         "bottom": int(np.count_nonzero(halves == 1)),
                         "straddling": int(np.count_nonzero(halves < 0))}
    return rings


def band_counts(rings: dict[int, dict[str, int]], support: int) -> dict[str, dict[str, object]]:
    """Aggregate ring counts into the frozen radial half-widths."""
    pixels = support * support
    result = {}
    for width in BAND_WIDTHS[support]:
        selected = [value for offset, value in rings.items() if abs(offset) <= width]
        counts = {key: sum(value[key] for value in selected)
                  for key in ("attempted", "accepted", "top", "bottom", "straddling")}
        counts["centered_rank_ceiling"] = min(pixels, max(0, counts["accepted"] - 1))
        counts["top_centered_rank_ceiling"] = min(pixels, max(0, counts["top"] - 1))
        counts["bottom_centered_rank_ceiling"] = min(pixels, max(0, counts["bottom"] - 1))
        counts["split_minimum_passed"] = min(counts["top"], counts["bottom"]) >= MINIMUM_SPLIT_PATCHES
        result[str(width)] = counts
    return result


def geometry_summary(records: list[dict[str, object]], support: int) -> dict[str, object]:
    """Summarize training counts across selected sites and all five searches."""
    widths = {}
    for width in BAND_WIDTHS[support]:
        rows = [record["bands"][str(width)] for record in records]
        widths[str(width)] = {
            key: summarize(np.asarray([row[key] for row in rows], dtype=np.float64))
            for key in ("accepted", "top", "bottom", "straddling", "centered_rank_ceiling")
        }
        widths[str(width)]["all_queries_split_minimum_passed"] = all(
            row["split_minimum_passed"] for row in rows)
    passing = [width for width in BAND_WIDTHS[support]
               if widths[str(width)]["all_queries_split_minimum_passed"]]
    return {"query_count": len(records), "band_half_widths": widths,
            "narrowest_split_supported_half_width": passing[0] if passing else None}


def analyze(root: Path, protocol: dict[str, object]) -> dict[str, object]:
    """Measure footprint energy, common search support, and training coverage."""
    parent = Path(str(protocol["parent_response"]))
    paths = response47.product_paths(parent)
    coordinates = np.asarray(fits.getdata(paths["coordinates"]), dtype=np.float64).T
    baseline = np.asarray(fits.getdata(paths["final"]), dtype=np.float64)
    center_row = 0.5 * (baseline.shape[2] - 1)
    center_column = 0.5 * (baseline.shape[1] - 1)
    radii = np.hypot(coordinates[:, 0] - center_row, coordinates[:, 1] - center_column)
    planet_row = center_row - stage.PLANET_SEPARATION * math.sin(math.radians(stage.PLANET_PA))
    planet_column = center_column + stage.PLANET_SEPARATION * math.cos(math.radians(stage.PLANET_PA))

    all_mode_complete = {support: np.ones(len(coordinates), dtype=bool) for support in SUPPORTS}
    response_results = {}
    for mode, response_path, validity_path in zip(stage.MODES, paths["responses"], paths["validities"]):
        response = np.asarray(fits.getdata(response_path, memmap=True))
        validity = np.asarray(fits.getdata(validity_path, memmap=True)) > 0.5
        safe = np.where(validity, response, 0)
        full_energy = np.einsum("ijk,ijk->i", safe, safe, dtype=np.float64)
        stage.require(np.all(full_energy > 0), f"mode {mode} contains zero-energy responses")
        supports = {}
        for support in SUPPORTS:
            complete, supports[str(support)] = response_summary(
                response, validity, radii, support, full_energy)
            all_mode_complete[support] &= complete
        response_results[str(mode)] = supports
        del response, validity, safe

    finite = np.all(np.isfinite(baseline), axis=0)
    geometry = {}
    for support in SUPPORTS:
        data_complete = complete_data_support(baseline, coordinates, support)
        support_geometry = {}
        for nominal in RADII:
            centers = eligible_centers(coordinates, radii, all_mode_complete[support], data_complete,
                                        nominal, planet_row, planet_column)
            sites = select_sites(centers)
            query_records = []
            for site_index, site in enumerate(sites):
                searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                            for delta_row, delta_column in SEARCH_OFFSETS]
                for search_index, query in enumerate(searches):
                    rings = training_rings(finite, query, searches, support, (planet_row, planet_column))
                    query_records.append({"site_index": site_index, "search_index": search_index,
                                          "row": query[0], "column": query[1],
                                          "bands": band_counts(rings, support)})
            support_geometry[str(nominal)] = {
                "eligible_five_pixel_centers": len(centers),
                "selected_sites": sites,
                "training": geometry_summary(query_records, support) if query_records else None,
            }
        geometry[str(support)] = support_geometry

    result = {
        "purpose": protocol["purpose"],
        "response": response_results,
        "geometry": geometry,
        "common_mode_response_completeness": {
            str(support): int(np.count_nonzero(all_mode_complete[support])) for support in SUPPORTS},
        "known_planet_row_column": [planet_row, planet_column],
        "selection_uses_baseline_values": False,
        "baseline_used_only_for_common_finite_mask": True,
    }
    stage.write_json(root / "results.json", result)
    write_report(root, result)
    return result


def write_report(root: Path, result: dict[str, object]) -> None:
    """Write compact mode-200 footprint and geometry tables."""
    lines = ["# KLIP Stage-B footprint preflight", "", "## Mode-200 response energy", "",
             "Values are median fractions of the complete 47-pixel response energy.", "",
             "| Radius | 11 pixels | 31 pixels | 47 pixels |", "| ---: | ---: | ---: | ---: |"]
    mode = result["response"]["200"]
    for radius in RADII:
        values = [mode[str(support)]["bins"][str(radius)]["energy_fraction_of_47"]["median"]
                  for support in SUPPORTS]
        lines.append(f"| {radius:g} | {values[0]:.4f} | {values[1]:.4f} | {values[2]:.4f} |")
    lines.extend(["", "## Common eight-mode five-pixel search geometry", "",
                  "| Radius | 11 pixels | 31 pixels | 47 pixels |",
                  "| ---: | ---: | ---: | ---: |"])
    for radius in RADII:
        values = [result["geometry"][str(support)][str(radius)]["eligible_five_pixel_centers"]
                  for support in SUPPORTS]
        lines.append(f"| {radius:g} | {values[0]} | {values[1]} | {values[2]} |")
    lines.extend(["", "## Narrowest radial band with eight patches in each disjoint half", "",
                  "Each entry is the radial half-width in pixels across all selected sites and their five search "
                  "pixels. A dash means no tested band passes.", "",
                  "| Radius | 11 pixels | 31 pixels | 47 pixels |",
                  "| ---: | ---: | ---: | ---: |"])
    for radius in RADII:
        values = []
        for support in SUPPORTS:
            training = result["geometry"][str(support)][str(radius)]["training"]
            width = training["narrowest_split_supported_half_width"] if training else None
            values.append("—" if width is None else str(width))
        lines.append(f"| {radius:g} | {values[0]} | {values[1]} | {values[2]} |")
    lines.extend(["", "The preflight uses response values only for template-energy diagnostics. Candidate and "
                  "training selection use coordinates, validity, and the baseline finite mask; no baseline scores "
                  "or covariance values are read."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    """Run the footprint preflight and write immutable completion receipts."""
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
    """Check the frozen support grids and band aggregation algebra."""
    stage.require(tuple((support - 1) // 2 for support in SUPPORTS) == (5, 15, 23),
                  "half-overlap steps changed")
    rings = {-5: {"attempted": 10, "accepted": 8, "top": 4, "bottom": 4, "straddling": 0},
             0: {"attempted": 12, "accepted": 10, "top": 5, "bottom": 4, "straddling": 1},
             5: {"attempted": 14, "accepted": 12, "top": 6, "bottom": 6, "straddling": 0}}
    bands = band_counts(rings, 11)
    stage.require(bands["0"]["accepted"] == 10 and bands["5"]["accepted"] == 30 and
                  bands["5"]["centered_rank_ceiling"] == 29 and
                  bands["5"]["split_minimum_passed"], "band aggregation changed")
    stage.require(BAND_WIDTHS[47][-1] >= 60, "47-pixel grid no longer reaches the full radial range")
    print("KLIP Stage-B footprint preflight checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the footprint-preflight command-line parser."""
    repo = Path(__file__).resolve().parents[3]
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--response", type=Path,
                                default=repo / "working/roc/klip_response_47_20260921")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    return result


def main() -> None:
    """Dispatch the requested footprint-preflight action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    elif args.action == "prepare":
        prepare(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
