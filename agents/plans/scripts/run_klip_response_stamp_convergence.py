#!/usr/bin/env python3
"""Measure KLIP response-stamp convergence at fixed geometry-only sites.

Each selected site receives one positive and one negative signal-free source
perturbation. Their full-image central difference is extracted at several
stamp sizes, so widening the diagnostic footprint does not require another
reduction. The 11-pixel extraction is also compared with the archived native
exact response at the same integer coordinate.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_klip_covariance_stage_a as stage  # noqa: E402


RADII = stage.PRIMARY_RADII
STAMP_SIZES = [11, 15, 19, 23, 31, 39, 47, 55, 63]
DECISION_STAMP_SIZES = [15, 19, 23, 31]
SITES_PER_RADIUS = 12


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def safe_radius_tag(radius: float) -> str:
    """Format a radius for a task name."""
    return format(radius, "g").replace(".", "p")


def vector(header: fits.Header, keyword: str) -> list[float]:
    """Read a comma-separated numeric FITS-header vector."""
    return [float(token) for token in str(header[keyword]).split(",") if token.strip()]


def select_sites(parent_results: dict[str, object], exact_manifest: Path) -> list[dict[str, object]]:
    """Select fixed angularly distributed sites using geometry only."""
    products = stage.exact_products(exact_manifest)
    coordinates = np.asarray(fits.getdata(products["coordinates"]), dtype=np.float64).T
    lookup = {(int(row), int(column)): index for index, (row, column, _, _) in enumerate(coordinates)}
    sites = []
    for radius in RADII:
        candidates = parent_results["response_audit"]["geometry"][str(radius)]["centers"]
        stage.require(len(candidates) >= SITES_PER_RADIUS,
                      f"radius {radius:g} has too few convergence candidates")
        ordered = sorted(candidates, key=lambda item: (float(item["angle_radians"]), item["row"], item["column"]))
        selected_indices = np.floor((np.arange(SITES_PER_RADIUS) + 0.5) * len(ordered) /
                                    SITES_PER_RADIUS).astype(int)
        stage.require(len(set(map(int, selected_indices))) == SITES_PER_RADIUS,
                      f"radius {radius:g} angular selection is not unique")
        for site_index, candidate_index in enumerate(selected_indices):
            candidate = ordered[int(candidate_index)]
            row, column = int(candidate["row"]), int(candidate["column"])
            delta_row, delta_column = row - 63.5, column - 63.5
            separation = math.hypot(delta_row, delta_column)
            position_angle = math.degrees(math.atan2(-delta_row, delta_column)) % 360
            stage.require((row, column) in lookup, "selected convergence site lacks an exact response")
            sites.append({"name": f"r{safe_radius_tag(radius)}_s{site_index:02d}",
                          "nominal_radius": radius, "site_index": site_index, "row": row, "column": column,
                          "separation": separation, "position_angle": position_angle,
                          "exact_source_index": lookup[(row, column)],
                          "geometry_selection_index": int(candidate_index)})
    stage.require(len(sites) == len(RADII) * SITES_PER_RADIUS and
                  len({site["name"] for site in sites}) == len(sites), "site selection is incomplete or duplicated")
    return sites


def verify_parent(parent: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Verify the completed canonical Stage-A experiment."""
    protocol, _ = stage.load_protocol(parent)
    completion = read(parent / "complete.json")
    stage.verify([completion["results"], completion["report"]])
    state = read(parent / "state.json")
    results = read(parent / "results.json")
    stage.require(state["status"] == "complete" and completion["status"] == "complete",
                  "parent Stage-A experiment is incomplete")
    stage.require(results["stage_a_passed"] and results["baseline"]["compatible"],
                  "parent Stage-A compatibility gates did not pass")
    stage.require(results["larger_stamp_experiment_required"],
                  "parent Stage-A result did not activate stamp convergence")
    stage.require(protocol["binary_provenance"]["hash_matches"]["klipReduce"],
                  "parent Stage-A reduction binary does not match the response archive")
    return protocol, results


def prepare(args: argparse.Namespace) -> None:
    """Freeze the convergence sites, tasks, inputs, software, and thresholds."""
    root, parent = args.root.resolve(), args.stage_a.resolve()
    stage.require(not root.exists(), f"output already exists: {root}")
    parent_protocol, parent_results = verify_parent(parent)
    exact_manifest = Path(str(parent_protocol["paths"]["exact_manifest"]))
    sites = select_sites(parent_results, exact_manifest)
    resources = parent_protocol["resources"]
    stage.require(sorted(os.sched_getaffinity(0)) == resources["cpu_affinity"],
                  "prepare must use the parent Stage-A CPU affinity")
    tasks = [{"name": site["name"] + "_" + sign, "site": site["name"], "sign": sign,
              "contrast": stage.PLANET_CONTRAST * (1 if sign == "plus" else -1)}
             for site in sites for sign in ("plus", "minus")]

    root.mkdir(parents=True)
    (root / "software").mkdir()
    runner = Path(__file__).resolve()
    stage_runner = runner.with_name("run_klip_covariance_stage_a.py")
    shutil.copy2(runner, root / "software" / runner.name)
    shutil.copy2(stage_runner, root / "software" / stage_runner.name)
    shutil.copy2(parent / "reduction.conf", root / "reduction.conf")
    shutil.copy2(parent / "inputs.txt", root / "inputs.txt")

    protocol = {
        "schema": 1,
        "purpose": "choose the smallest response footprint that resolves the Stage-A edge-energy trigger",
        "parent_stage_a": str(parent),
        "modes": stage.MODES,
        "primary_mode": 200,
        "radii": RADII,
        "sites_per_radius": SITES_PER_RADIUS,
        "site_selection": "geometry-only candidates sorted by angle; 12 equal-count bin midpoints per radius",
        "sites": sites,
        "tasks": tasks,
        "perturbation_half_contrast": stage.PLANET_CONTRAST,
        "known_planet": {"separation": stage.PLANET_SEPARATION, "position_angle": stage.PLANET_PA,
                         "contrast": stage.PLANET_CONTRAST},
        "stamp_sizes": STAMP_SIZES,
        "decision_stamp_sizes": DECISION_STAMP_SIZES,
        "reference_stamp_size": 63,
        "edge_thresholds": {"median_fraction": 0.01, "individual_fraction": 0.05},
        "response_replay_gate": {"minimum_cosine": 0.999, "projection_absolute_tolerance": 0.01},
        "resources": resources,
        "paths": {"klipreduce": parent_protocol["paths"]["klipreduce"],
                  "psf": parent_protocol["paths"]["psf"], "exact_manifest": str(exact_manifest),
                  "parent_baseline": str(parent / "current_baseline/finim.fits")},
        "expected_reductions": len(tasks),
    }
    stage.write_json(root / "protocol.json", protocol)
    records = [stage.fingerprint(root / name) for name in ("reduction.conf", "inputs.txt", "protocol.json")]
    records.extend(stage.fingerprint(root / "software" / name) for name in
                   (runner.name, stage_runner.name))
    parent_records = [stage.fingerprint(parent / name) for name in
                      ("protocol.json", "manifest.json", "results.json", "complete.json",
                       "current_baseline/finim.fits")]
    external_records = [stage.fingerprint(Path(protocol["paths"][name])) for name in
                        ("klipreduce", "psf", "exact_manifest")]
    exact_products = stage.product_paths(exact_manifest)
    external_records.extend(stage.fingerprint(path) for path in exact_products if path != exact_manifest)
    stage.write_json(root / "manifest.json", {"schema": 1, "frozen_records": records,
                                              "parent_records": parent_records,
                                              "external_records": external_records})
    stage.write_json(root / "state.json", {"status": "prepared", "completed_tasks": []})
    print(root, flush=True)


def load_experiment(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Load a convergence experiment and verify all frozen inputs."""
    protocol = read(root / "protocol.json")
    manifest = read(root / "manifest.json")
    stage.verify(manifest["frozen_records"] + manifest["parent_records"] + manifest["external_records"])
    parent_protocol, _ = verify_parent(Path(protocol["parent_stage_a"]))
    stage.require(parent_protocol["resources"] == protocol["resources"], "parent resource contract changed")
    stage.require(sorted(os.sched_getaffinity(0)) == protocol["resources"]["cpu_affinity"],
                  "CPU affinity changed from the frozen protocol")
    return protocol, manifest


def task_command(root: Path, protocol: dict[str, object], task: dict[str, object]) -> list[str]:
    """Construct one signal-free positive or negative perturbation command."""
    site = next(site for site in protocol["sites"] if site["name"] == task["site"])
    paths = protocol["paths"]
    return [str(paths["klipreduce"]), "--config", str(root / "reduction.conf"), "--input.directory=",
            "--input.fileList", str(root / "inputs.txt"), "--klip.Nmodes", ",".join(map(str, stage.MODES)),
            "--psfResponse.file=", "--psfResponse.outputModels=false", "--psfResponse.filter=false",
            "--planet.sep", str(stage.PLANET_SEPARATION), "--planet.PA", str(stage.PLANET_PA),
            "--planet.contrast", str(stage.PLANET_CONTRAST), "--fake.method", "single",
            "--fake.fileName", str(paths["psf"]), "--fake.sep", str(site["separation"]),
            "--fake.PA", str(site["position_angle"]), "--fake.contrast", str(task["contrast"]),
            "--fake.subtractPlanet=true", "--output.directory", str(root / "reductions" / task["name"]),
            "--output.fileName", "finim.fits", "--output.exactFName=true", "--showTiming=true"]


def validate_reduction(path: Path, task: dict[str, object], site: dict[str, object], protocol: dict[str, object]) -> None:
    """Validate one completed perturbation cube and its source metadata."""
    data, header = fits.getdata(path, header=True)
    stage.require(data.shape == (len(stage.MODES), 128, 128) and stage.read_modes(header) == stage.MODES,
                  f"unexpected perturbation cube: {path}")
    fake_separation, fake_pa, fake_contrast = (vector(header, key) for key in ("FAKESEP", "FAKEPA", "FAKECONT"))
    planet_separation, planet_pa, planet_contrast = (vector(header, key) for key in
                                                     ("PLANETSEP", "PLANETPA", "PLANETCONT"))
    stage.require(len(fake_separation) == len(fake_pa) == len(fake_contrast) == 1,
                  f"perturbation metadata is incomplete: {path}")
    stage.require(len(planet_separation) == len(planet_pa) == len(planet_contrast) == 1,
                  f"known-planet metadata is incomplete: {path}")
    stage.require(math.isclose(fake_separation[0], site["separation"], rel_tol=0, abs_tol=5e-4) and
                  math.isclose(fake_pa[0], site["position_angle"], rel_tol=0, abs_tol=5e-4) and
                  math.isclose(fake_contrast[0], task["contrast"], rel_tol=0, abs_tol=5e-9),
                  f"perturbation source metadata changed: {path}")
    stage.require(math.isclose(planet_separation[0], stage.PLANET_SEPARATION, rel_tol=0, abs_tol=5e-4) and
                  math.isclose(planet_pa[0], stage.PLANET_PA, rel_tol=0, abs_tol=5e-4) and
                  math.isclose(planet_contrast[0], stage.PLANET_CONTRAST, rel_tol=0, abs_tol=5e-9),
                  f"known-planet subtraction metadata changed: {path}")
    stage.require(Path(str(header["FAKEFILE"]).strip()).resolve() == Path(protocol["paths"]["psf"]).resolve(),
                  f"perturbation PSF metadata changed: {path}")


def run(args: argparse.Namespace) -> None:
    """Run or resume every frozen perturbation task and analyze on completion."""
    root = args.root.resolve()
    protocol, _ = load_experiment(root)
    if (root / "complete.json").is_file():
        completion = read(root / "complete.json")
        stage.verify([completion["results"], completion["report"], *completion["reduction_receipts"]])
        for task in protocol["tasks"]:
            receipt = read(root / "reductions" / task["name"] / "complete.json")
            stage.verify(receipt["products"])
        stage.require(read(root / "state.json")["status"] == "complete", "completion state is inconsistent")
        print(root / "results.md", flush=True)
        return
    reductions = root / "reductions"
    reductions.mkdir(exist_ok=True)
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(protocol["resources"]["openmp_threads"]), OMP_PROC_BIND="true",
                       OMP_PLACES="cores", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    sites = {site["name"]: site for site in protocol["sites"]}
    completed = []
    limit = args.maximum_tasks if args.maximum_tasks > 0 else len(protocol["tasks"])
    stage.write_json(root / "state.json", {"status": "running", "pid": os.getpid(), "completed_tasks": completed})
    try:
        for task in protocol["tasks"]:
            directory = reductions / task["name"]
            receipt = directory / "complete.json"
            if receipt.is_file():
                record = read(receipt)
                stage.verify(record["products"])
                validate_reduction(directory / "finim.fits", task, sites[task["site"]], protocol)
                completed.append(task["name"])
                continue
            if len(completed) >= limit:
                break
            stage.require(not directory.exists(), f"incomplete task exists; archive before retry: {directory}")
            directory.mkdir()
            command = task_command(root, protocol, task)
            elapsed = stage.run_command(command, directory, directory / "run.log", environment)
            validate_reduction(directory / "finim.fits", task, sites[task["site"]], protocol)
            stage.write_json(receipt, {"task": task, "elapsed_seconds": elapsed,
                                      "products": [stage.fingerprint(directory / "finim.fits")]})
            completed.append(task["name"])
            stage.write_json(root / "state.json", {"status": "running", "pid": os.getpid(),
                                                   "completed_tasks": completed})
            print(f"completed {len(completed)}/{len(protocol['tasks'])}: {task['name']}", flush=True)
        if len(completed) == len(protocol["tasks"]):
            analyze(root, protocol)
        else:
            stage.write_json(root / "state.json", {"status": "partial", "completed_tasks": completed,
                                                   "remaining_tasks": len(protocol["tasks"]) - len(completed)})
            print(f"partial run: {len(completed)}/{len(protocol['tasks'])} tasks complete", flush=True)
    except Exception as error:
        stage.write_json(root / "state.json", {"status": "failed", "completed_tasks": completed,
                                               "error_type": type(error).__name__, "error": str(error)})
        raise


def extract_stamp(image: np.ndarray, row: int, column: int, size: int) -> np.ndarray:
    """Extract an Eigen-oriented square stamp from a FITS-oriented image."""
    half = size // 2
    stamp = image[column - half:column + half + 1, row - half:row + half + 1].T
    stage.require(stamp.shape == (size, size), "response stamp extends outside the image")
    return stamp


def one_response_metrics(response: np.ndarray, size: int) -> dict[str, object]:
    """Measure support, energy, sign, and border content of one response stamp."""
    valid = np.isfinite(response)
    energy = float(np.sum(np.square(response[valid])))
    border = np.zeros((size, size), dtype=bool)
    border[[0, -1], :] = True
    border[:, [0, -1]] = True
    stage.require(energy > 0, "finite-difference response has no positive energy")
    return {"support_fraction": float(np.count_nonzero(valid)) / response.size,
            "energy": energy,
            "border_energy_fraction": float(np.sum(np.square(response[valid & border]))) / energy,
            "negative_energy_fraction": float(np.sum(np.square(response[valid & (response < 0)]))) / energy}


def replay_metrics(measured: np.ndarray, exact: np.ndarray, validity: np.ndarray) -> dict[str, object]:
    """Compare an external finite difference with the archived exact response."""
    support = validity & np.isfinite(measured) & np.isfinite(exact)
    first, second = measured[support], exact[support]
    first_energy, second_energy = float(first @ first), float(second @ second)
    stage.require(first_energy > 0 and second_energy > 0, "response replay has zero common-support energy")
    cross = float(first @ second)
    projection = cross / second_energy
    residual = math.sqrt(float(np.sum(np.square(first - projection * second))) / first_energy)
    return {"support_fraction": float(np.count_nonzero(support)) / measured.size,
            "cosine": cross / math.sqrt(first_energy * second_energy),
            "exact_to_measured_projection": projection, "best_scaled_relative_residual": residual,
            "maximum_absolute_difference": float(np.max(np.abs(first - second)))}


def summarize(values: list[float]) -> dict[str, object]:
    """Return count and fixed quantiles for finite diagnostic values."""
    array = np.asarray(values, dtype=np.float64)
    array = array[np.isfinite(array)]
    stage.require(len(array) > 0, "cannot summarize an empty diagnostic")
    return {"count": len(array), "minimum": float(np.min(array)), "median": float(np.median(array)),
            "p90": float(np.percentile(array, 90)), "maximum": float(np.max(array))}


def analyze(root: Path, protocol: dict[str, object]) -> None:
    """Analyze completed pairs, choose a converged footprint, and write receipts."""
    sites = protocol["sites"]
    exact = stage.exact_products(Path(protocol["paths"]["exact_manifest"]))
    exact_responses = [np.asarray(fits.getdata(path), dtype=np.float64) for path in exact["responses"]]
    exact_validities = [np.asarray(fits.getdata(path), dtype=np.float64) > 0.5 for path in exact["validities"]]
    baseline = np.asarray(fits.getdata(protocol["paths"]["parent_baseline"]), dtype=np.float64)
    measurements = []
    for site_number, site in enumerate(sites, start=1):
        plus = np.asarray(fits.getdata(root / "reductions" / (site["name"] + "_plus") / "finim.fits"),
                          dtype=np.float64)
        minus = np.asarray(fits.getdata(root / "reductions" / (site["name"] + "_minus") / "finim.fits"),
                           dtype=np.float64)
        response_cube = (plus - minus) / (2 * protocol["perturbation_half_contrast"])
        midpoint = 0.5 * (plus + minus)
        source = int(site["exact_source_index"])
        for mode_index, mode in enumerate(stage.MODES):
            sizes = {}
            for size in protocol["stamp_sizes"]:
                response = extract_stamp(response_cube[mode_index], site["row"], site["column"], size)
                sizes[str(size)] = one_response_metrics(response, size)
            reference_energy = sizes[str(protocol["reference_stamp_size"])]["energy"]
            for size in protocol["stamp_sizes"]:
                sizes[str(size)]["energy_fraction_of_reference"] = sizes[str(size)]["energy"] / reference_energy
            measured_11 = extract_stamp(response_cube[mode_index], site["row"], site["column"], 11)
            archived = exact_responses[mode_index][source].T
            validity = exact_validities[mode_index][source].T
            replay = replay_metrics(measured_11, archived, validity)
            midpoint_stamp = extract_stamp(midpoint[mode_index] - baseline[mode_index],
                                           site["row"], site["column"], 31)
            midpoint_values = midpoint_stamp[np.isfinite(midpoint_stamp)]
            measurements.append({"site": site["name"], "nominal_radius": site["nominal_radius"],
                                 "row": site["row"], "column": site["column"], "mode": mode,
                                 "sizes": sizes, "archived_11_replay": replay,
                                 "midpoint_baseline_rms_31": float(np.sqrt(np.mean(np.square(midpoint_values)))),
                                 "midpoint_baseline_maximum_31": float(np.max(np.abs(midpoint_values)))})
        print(f"analyzed {site_number}/{len(sites)}: {site['name']}", flush=True)

    summaries = {}
    for mode in stage.MODES:
        mode_summary = {}
        for radius in RADII:
            rows = [row for row in measurements if row["mode"] == mode and row["nominal_radius"] == radius]
            size_summary = {}
            for size in protocol["stamp_sizes"]:
                values = [row["sizes"][str(size)] for row in rows]
                size_summary[str(size)] = {
                    key: summarize([value[key] for value in values]) for key in
                    ("support_fraction", "border_energy_fraction", "negative_energy_fraction",
                     "energy_fraction_of_reference")}
            replay_values = [row["archived_11_replay"] for row in rows]
            mode_summary[str(radius)] = {"sites": len(rows), "sizes": size_summary,
                                         "archived_11_replay": {key: summarize([value[key] for value in replay_values])
                                                                 for key in ("cosine", "exact_to_measured_projection",
                                                                             "best_scaled_relative_residual",
                                                                             "maximum_absolute_difference")}}
        summaries[str(mode)] = mode_summary

    replay_rows = [row["archived_11_replay"] for row in measurements]
    minimum_cosine = min(row["cosine"] for row in replay_rows)
    maximum_projection_error = max(abs(row["exact_to_measured_projection"] - 1) for row in replay_rows)
    replay_passed = (minimum_cosine >= protocol["response_replay_gate"]["minimum_cosine"] and
                     maximum_projection_error <= protocol["response_replay_gate"]["projection_absolute_tolerance"])
    sufficient = []
    for size in protocol["decision_stamp_sizes"]:
        groups = [summaries["200"][str(radius)]["sizes"][str(size)] for radius in RADII]
        if (all(group["support_fraction"]["minimum"] == 1 for group in groups) and
                all(group["border_energy_fraction"]["median"] <= protocol["edge_thresholds"]["median_fraction"]
                for group in groups) and
                all(group["border_energy_fraction"]["maximum"] <= protocol["edge_thresholds"]["individual_fraction"]
                    for group in groups)):
            sufficient.append(size)
    result = {"response_replay_passed": replay_passed, "minimum_archived_replay_cosine": minimum_cosine,
              "maximum_archived_replay_projection_error": maximum_projection_error,
              "sufficient_stamp_sizes_mode200": sufficient,
              "selected_smallest_sufficient_stamp_size": min(sufficient) if sufficient else None,
              "summaries": summaries, "measurements": measurements}
    stage.write_json(root / "results.json", result)
    write_report(root, protocol, result)
    stage.require(replay_passed, "external 11-pixel responses do not reproduce the archived exact products")
    completion = {"status": "complete", "results": stage.fingerprint(root / "results.json"),
                  "report": stage.fingerprint(root / "results.md"),
                  "reduction_receipts": [stage.fingerprint(root / "reductions" / task["name"] / "complete.json")
                                         for task in protocol["tasks"]]}
    stage.write_json(root / "complete.json", completion)
    stage.write_json(root / "state.json", {"status": "complete",
                                           "completed_tasks": [task["name"] for task in protocol["tasks"]]})
    print(root / "results.md", flush=True)


def write_report(root: Path, protocol: dict[str, object], result: dict[str, object]) -> None:
    """Write the compact response-stamp convergence report."""
    selected = result["selected_smallest_sufficient_stamp_size"]
    lines = ["# KLIP response-stamp convergence", "",
             f"The study used {len(protocol['sites'])} fixed geometry-only sites and "
             f"{len(protocol['tasks'])} signal-free KLIP reductions. Each plus/minus pair supplies every stamp size.",
             "", "## Gates", "", "| Check | Result | Diagnostic |", "| --- | --- | --- |",
             f"| External 11-pixel response versus archive | {'pass' if result['response_replay_passed'] else 'fail'} | minimum cosine {result['minimum_archived_replay_cosine']:.6f}; maximum projection error {result['maximum_archived_replay_projection_error']:.6g} |",
             f"| Smallest footprint meeting edge thresholds | {selected if selected is not None else 'none'} | mode 200 across all six primary radii |",
             "", "## Mode-200 footprint diagnostics", "",
             "Support is median / minimum; energy entries are median / maximum across the 12 sites at that radius. "
             "Energy fractions use the 63-pixel diagnostic extraction as their reference.", "",
             "| Radius | Stamp | Support | Border energy | Energy / 63-pixel energy |", "| ---: | ---: | ---: | ---: | ---: |"]
    for radius in RADII:
        for size in protocol["stamp_sizes"]:
            summary = result["summaries"]["200"][str(radius)]["sizes"][str(size)]
            lines.append(f"| {radius:g} | {size} | {summary['support_fraction']['median']:.3f} / "
                         f"{summary['support_fraction']['minimum']:.3f} | "
                         f"{summary['border_energy_fraction']['median']:.4f} / "
                         f"{summary['border_energy_fraction']['maximum']:.4f} | "
                         f"{summary['energy_fraction_of_reference']['median']:.4f} / "
                         f"{summary['energy_fraction_of_reference']['maximum']:.4f} |")
    if selected is None:
        lines.extend(["", "No tested footprint satisfies both preregistered edge-energy thresholds at every primary radius."])
    else:
        lines.extend(["", f"The smallest tested footprint satisfying both edge-energy thresholds is {selected} by {selected} pixels."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def check() -> None:
    """Run deterministic checks of site selection and stamp-energy arithmetic."""
    image = np.zeros((41, 41), dtype=np.float64)
    yy, xx = np.indices(image.shape)
    image = np.exp(-0.5 * ((xx - 20) ** 2 + (yy - 20) ** 2) / 9)
    small = extract_stamp(image, 20, 20, 11)
    large = extract_stamp(image, 20, 20, 31)
    small_metrics, large_metrics = one_response_metrics(small, 11), one_response_metrics(large, 31)
    stage.require(small_metrics["energy"] < large_metrics["energy"] and
                  small_metrics["border_energy_fraction"] > large_metrics["border_energy_fraction"],
                  "stamp convergence arithmetic check failed")
    measured = np.arange(121, dtype=np.float64).reshape(11, 11) - 60
    replay = replay_metrics(measured, measured.copy(), np.ones((11, 11), dtype=bool))
    stage.require(np.isclose(replay["cosine"], 1) and np.isclose(replay["exact_to_measured_projection"], 1) and
                  np.isclose(replay["best_scaled_relative_residual"], 0), "response replay check failed")
    print("KLIP response-stamp convergence checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the convergence runner command-line parser."""
    repo = Path(__file__).resolve().parents[3]
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--stage-a", type=Path,
                                default=repo / "working/roc/klip_covariance_stage_a_20260921")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    run_parser.add_argument("--maximum-tasks", type=int, default=0,
                            help="debug/resume limit; zero runs every remaining task")
    return result


def main() -> None:
    """Dispatch the requested convergence action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    elif args.action == "prepare":
        prepare(args)
    else:
        stage.require(args.maximum_tasks >= 0, "maximum-tasks must be nonnegative")
        run(args)


if __name__ == "__main__":
    main()
