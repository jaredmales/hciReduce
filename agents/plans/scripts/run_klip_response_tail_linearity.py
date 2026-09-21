#!/usr/bin/env python3
"""Test whether the extended KLIP response tail is contrast-linear.

The response-stamp convergence experiment supplies the reference central
difference at the fitted-planet contrast. This experiment repeats a fixed
geometry-only subset at half and twice that contrast. It compares the three
derivatives in square shells, with the 31-to-47-pixel shell controlling whether
the 47-pixel diagnostic footprint can be promoted to a full response campaign.
"""
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
import run_klip_response_stamp_convergence as convergence  # noqa: E402


RADII = [7.5, 10.0, 12.0]
SITE_INDICES = [0, 2, 4, 6, 8, 10]
CONTRAST_SCALES = [0.5, 2.0]
REFERENCE_SCALE = 1.0
STAMP_SIZE = 63
PROMOTED_SIZE = 47
REGIONS = {
    "core_11": (None, 5),
    "middle_11_to_31": (5, 15),
    "tail_31_to_47": (15, 23),
    "outer_47_to_63": (23, 31),
    "full_31": (None, 15),
    "full_47": (None, 23),
    "full_63": (None, 31),
}


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def scale_tag(scale: float) -> str:
    """Format a contrast scale for a task name."""
    return format(scale, "g").replace(".", "p")


def verify_parent(parent: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Verify the completed response-stamp convergence experiment."""
    protocol, _ = convergence.load_experiment(parent)
    completion = read(parent / "complete.json")
    state = read(parent / "state.json")
    results = read(parent / "results.json")
    stage.require(completion["status"] == "complete" and state["status"] == "complete",
                  "parent response-stamp experiment is incomplete")
    stage.verify([completion["results"], completion["report"], *completion["reduction_receipts"]])
    for task in protocol["tasks"]:
        receipt = read(parent / "reductions" / task["name"] / "complete.json")
        stage.verify(receipt["products"])
    stage.require(results["response_replay_passed"], "parent archived-response replay failed")
    stage.require(results["selected_smallest_sufficient_stamp_size"] is None,
                  "parent unexpectedly selected a candidate footprint")

    groups = [results["summaries"][str(mode)][str(radius)]["sizes"][str(PROMOTED_SIZE)]
              for mode in stage.MODES for radius in stage.PRIMARY_RADII]
    thresholds = protocol["edge_thresholds"]
    stage.require(all(group["support_fraction"]["minimum"] == 1 for group in groups),
                  "parent 47-pixel diagnostic lacks complete support")
    stage.require(all(group["border_energy_fraction"]["median"] <= thresholds["median_fraction"]
                      for group in groups), "parent 47-pixel median edge criterion failed")
    stage.require(all(group["border_energy_fraction"]["maximum"] <= thresholds["individual_fraction"]
                      for group in groups), "parent 47-pixel individual edge criterion failed")
    return protocol, results


def selected_sites(parent_protocol: dict[str, object]) -> list[dict[str, object]]:
    """Select six fixed, alternating angular sites at each inner radius."""
    sites = [site for site in parent_protocol["sites"]
             if float(site["nominal_radius"]) in RADII and int(site["site_index"]) in SITE_INDICES]
    sites.sort(key=lambda site: (float(site["nominal_radius"]), int(site["site_index"])))
    stage.require(len(sites) == len(RADII) * len(SITE_INDICES), "linearity site selection is incomplete")
    for radius in RADII:
        radius_sites = [site for site in sites if float(site["nominal_radius"]) == radius]
        stage.require([int(site["site_index"]) for site in radius_sites] == SITE_INDICES,
                      f"linearity sites changed at radius {radius:g}")
    return sites


def prepare(args: argparse.Namespace) -> None:
    """Freeze the parent responses, new contrasts, sites, tasks, and software."""
    root, parent = args.root.resolve(), args.convergence.resolve()
    stage.require(not root.exists(), f"output already exists: {root}")
    parent_protocol, _ = verify_parent(parent)
    resources = parent_protocol["resources"]
    stage.require(sorted(os.sched_getaffinity(0)) == resources["cpu_affinity"],
                  "prepare must use the parent CPU affinity")
    sites = selected_sites(parent_protocol)
    tasks = [{"name": f"{site['name']}_h{scale_tag(scale)}_{sign}", "site": site["name"],
              "scale": scale, "sign": sign,
              "contrast": stage.PLANET_CONTRAST * scale * (1 if sign == "plus" else -1)}
             for site in sites for scale in CONTRAST_SCALES for sign in ("plus", "minus")]

    root.mkdir(parents=True)
    (root / "software").mkdir()
    runner = Path(__file__).resolve()
    convergence_runner = runner.with_name("run_klip_response_stamp_convergence.py")
    stage_runner = runner.with_name("run_klip_covariance_stage_a.py")
    for path in (runner, convergence_runner, stage_runner):
        shutil.copy2(path, root / "software" / path.name)
    shutil.copy2(parent / "reduction.conf", root / "reduction.conf")
    shutil.copy2(parent / "inputs.txt", root / "inputs.txt")

    protocol = {
        "schema": 1,
        "purpose": "test whether the 31-to-47-pixel KLIP response tail is contrast-linear",
        "parent_convergence": str(parent),
        "modes": stage.MODES,
        "primary_mode": 200,
        "radii": RADII,
        "site_indices": SITE_INDICES,
        "site_selection": "alternating even-index sites from the parent's geometry-only angular selection",
        "sites": sites,
        "tasks": tasks,
        "reference_contrast_scale": REFERENCE_SCALE,
        "new_contrast_scales": CONTRAST_SCALES,
        "perturbation_base_contrast": stage.PLANET_CONTRAST,
        "stamp_size": STAMP_SIZE,
        "promoted_size": PROMOTED_SIZE,
        "regions": REGIONS,
        "linearity_gate": {"minimum_median_cosine": 0.95,
                           "projection_absolute_tolerance": 0.10,
                           "controlling_region": "tail_31_to_47"},
        "edge_thresholds": parent_protocol["edge_thresholds"],
        "known_planet": parent_protocol["known_planet"],
        "resources": resources,
        "paths": {"klipreduce": parent_protocol["paths"]["klipreduce"],
                  "psf": parent_protocol["paths"]["psf"],
                  "parent_baseline": parent_protocol["paths"]["parent_baseline"]},
        "expected_reductions": len(tasks),
    }
    stage.write_json(root / "protocol.json", protocol)

    records = [stage.fingerprint(root / name) for name in ("reduction.conf", "inputs.txt", "protocol.json")]
    records.extend(stage.fingerprint(root / "software" / path.name)
                   for path in (runner, convergence_runner, stage_runner))
    parent_records = [stage.fingerprint(parent / name) for name in
                      ("protocol.json", "manifest.json", "results.json", "complete.json")]
    for site in sites:
        for sign in ("plus", "minus"):
            directory = parent / "reductions" / f"{site['name']}_{sign}"
            parent_records.extend([stage.fingerprint(directory / "complete.json"),
                                   stage.fingerprint(directory / "finim.fits")])
    external_records = [stage.fingerprint(Path(protocol["paths"][name]))
                        for name in ("klipreduce", "psf", "parent_baseline")]
    stage.write_json(root / "manifest.json", {"schema": 1, "frozen_records": records,
                                              "parent_records": parent_records,
                                              "external_records": external_records})
    stage.write_json(root / "state.json", {"status": "prepared", "completed_tasks": []})
    print(root, flush=True)


def load_experiment(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Load the experiment and verify every frozen input and resource setting."""
    protocol = read(root / "protocol.json")
    manifest = read(root / "manifest.json")
    stage.verify(manifest["frozen_records"] + manifest["parent_records"] + manifest["external_records"])
    parent_protocol, parent_results = verify_parent(Path(protocol["parent_convergence"]))
    stage.require(parent_protocol["resources"] == protocol["resources"], "parent resource contract changed")
    stage.require(sorted(os.sched_getaffinity(0)) == protocol["resources"]["cpu_affinity"],
                  "CPU affinity changed from the frozen protocol")
    return protocol, parent_results


def task_command(root: Path, protocol: dict[str, object], task: dict[str, object]) -> list[str]:
    """Construct one half- or double-contrast perturbation command."""
    return convergence.task_command(root, protocol, task)


def validate_reduction(path: Path, task: dict[str, object], site: dict[str, object],
                       protocol: dict[str, object]) -> None:
    """Validate one new perturbation cube and its source metadata."""
    header = fits.getheader(path)
    recorded_contrast = convergence.vector(header, "FAKECONT")
    stage.require(len(recorded_contrast) == 1 and
                  math.isclose(recorded_contrast[0], float(task["contrast"]), rel_tol=1e-6, abs_tol=1e-12),
                  f"perturbation contrast metadata changed: {path}")
    serialized_task = dict(task)
    serialized_task["contrast"] = recorded_contrast[0]
    convergence.validate_reduction(path, serialized_task, site, protocol)


def repair(args: argparse.Namespace) -> None:
    """Repair the initial metadata-tolerance failure without discarding products."""
    root = args.root.resolve()
    protocol = read(root / "protocol.json")
    manifest = read(root / "manifest.json")
    repair_path = root / "metadata_validation_repair.json"
    runner_path = root / "software" / Path(__file__).name
    if repair_path.is_file():
        stage.verify(manifest["frozen_records"] + manifest["parent_records"] + manifest["external_records"])
        print(repair_path, flush=True)
        return

    stage.verify(manifest["frozen_records"] + manifest["parent_records"] + manifest["external_records"])
    state = read(root / "state.json")
    stage.require(state["status"] == "failed" and
                  state.get("error") == f"perturbation source metadata changed: "
                  f"{root / 'reductions/r7p5_s00_h2_plus/finim.fits'}",
                  "experiment is not at the recognized initial metadata-tolerance failure")
    tasks = {task["name"]: task for task in protocol["tasks"]}
    sites = {site["name"]: site for site in protocol["sites"]}
    incomplete = [directory for directory in (root / "reductions").iterdir()
                  if directory.is_dir() and not (directory / "complete.json").is_file()]
    stage.require([directory.name for directory in incomplete] == ["r7p5_s00_h2_plus"],
                  "unexpected incomplete reduction set; preserve it and inspect before repair")
    task = tasks[incomplete[0].name]
    product = incomplete[0] / "finim.fits"
    validate_reduction(product, task, sites[task["site"]], protocol)

    old_runner = stage.fingerprint(runner_path)
    temporary_runner = runner_path.with_suffix(".py.repair")
    shutil.copy2(Path(__file__).resolve(), temporary_runner)
    temporary_runner.replace(runner_path)
    new_runner = stage.fingerprint(runner_path)
    stage.write_json(incomplete[0] / "complete.json",
                     {"task": task, "elapsed_seconds": None,
                      "recovered_after_metadata_validation_fix": True,
                      "products": [stage.fingerprint(product)]})
    repair_record = {"schema": 1,
                     "reason": "FAKECONT uses six-significant-digit stream serialization; the original "
                               "5e-9 absolute tolerance missed the doubled contrast by 7.17e-12",
                     "validation": {"relative_tolerance": 1e-6, "absolute_tolerance": 1e-12},
                     "old_runner": old_runner, "new_runner": new_runner,
                     "recovered_task": task,
                     "recovered_product": stage.fingerprint(product)}
    stage.write_json(repair_path, repair_record)
    replaced = False
    for index, record in enumerate(manifest["frozen_records"]):
        if Path(str(record["path"])).resolve() == runner_path.resolve():
            manifest["frozen_records"][index] = new_runner
            replaced = True
            break
    stage.require(replaced, "frozen runner record is missing")
    manifest["frozen_records"].append(stage.fingerprint(repair_path))
    stage.write_json(root / "manifest.json", manifest)
    completed = [task["name"] for task in protocol["tasks"]
                 if (root / "reductions" / task["name"] / "complete.json").is_file()]
    stage.write_json(root / "state.json", {"status": "repaired", "completed_tasks": completed,
                                           "repair": stage.fingerprint(repair_path)})
    stage.verify(manifest["frozen_records"] + manifest["parent_records"] + manifest["external_records"])
    print(repair_path, flush=True)


def run(args: argparse.Namespace) -> None:
    """Run or resume every frozen perturbation and analyze on completion."""
    root = args.root.resolve()
    protocol, _ = load_experiment(root)
    if (root / "complete.json").is_file():
        completion = read(root / "complete.json")
        stage.verify([completion["results"], completion["report"], *completion["reduction_receipts"]])
        for task in protocol["tasks"]:
            receipt = read(root / "reductions" / task["name"] / "complete.json")
            stage.verify(receipt["products"])
        stage.require(read(root / "state.json")["status"] == "complete",
                      "completion state is inconsistent")
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
    stage.write_json(root / "state.json", {"status": "running", "pid": os.getpid(),
                                           "completed_tasks": completed})
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
            elapsed = stage.run_command(task_command(root, protocol, task), directory,
                                        directory / "run.log", environment)
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


def region_masks(size: int = STAMP_SIZE) -> dict[str, np.ndarray]:
    """Construct centered square core, shell, and cumulative masks."""
    coordinate = np.arange(size) - size // 2
    chebyshev_radius = np.maximum(np.abs(coordinate[:, None]), np.abs(coordinate[None, :]))
    masks = {}
    for name, (inner, outer) in REGIONS.items():
        mask = chebyshev_radius <= outer
        if inner is not None:
            mask &= chebyshev_radius > inner
        masks[name] = mask
    return masks


def comparison_metrics(candidate: np.ndarray, reference: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    """Compare a candidate derivative with the reference on one region."""
    support = mask & np.isfinite(candidate) & np.isfinite(reference)
    first, second = candidate[support], reference[support]
    first_energy, second_energy = float(first @ first), float(second @ second)
    stage.require(first_energy > 0 and second_energy > 0, "response-comparison region has zero energy")
    cross = float(first @ second)
    projection = cross / second_energy
    residual = math.sqrt(float(np.sum(np.square(first - projection * second))) / first_energy)
    relative_difference = math.sqrt(float(np.sum(np.square(first - second))) / second_energy)
    return {"support_fraction": float(np.count_nonzero(support)) / np.count_nonzero(mask),
            "cosine": cross / math.sqrt(first_energy * second_energy),
            "candidate_to_reference_projection": projection,
            "best_scaled_relative_residual": residual,
            "relative_difference": relative_difference,
            "energy_ratio": first_energy / second_energy}


def summarize(values: list[float]) -> dict[str, object]:
    """Return fixed quantiles for finite values."""
    return convergence.summarize(values)


def derivative(plus: np.ndarray, minus: np.ndarray, half_contrast: float) -> np.ndarray:
    """Calculate a central-difference derivative cube."""
    stage.require(half_contrast > 0, "central-difference half contrast must be positive")
    return (plus - minus) / (2 * half_contrast)


def analyze(root: Path, protocol: dict[str, object]) -> None:
    """Compare response derivatives by contrast, decide the tail gate, and write receipts."""
    parent = Path(protocol["parent_convergence"])
    baseline = np.asarray(fits.getdata(protocol["paths"]["parent_baseline"]), dtype=np.float64)
    masks = region_masks()
    measurements = []
    for site_number, site in enumerate(protocol["sites"], start=1):
        responses = {}
        midpoints = {}
        parent_plus = np.asarray(fits.getdata(parent / "reductions" / f"{site['name']}_plus" / "finim.fits"),
                                 dtype=np.float64)
        parent_minus = np.asarray(fits.getdata(parent / "reductions" / f"{site['name']}_minus" / "finim.fits"),
                                  dtype=np.float64)
        responses[REFERENCE_SCALE] = derivative(parent_plus, parent_minus, stage.PLANET_CONTRAST)
        midpoints[REFERENCE_SCALE] = 0.5 * (parent_plus + parent_minus)
        for scale in CONTRAST_SCALES:
            prefix = root / "reductions" / f"{site['name']}_h{scale_tag(scale)}"
            plus = np.asarray(fits.getdata(str(prefix) + "_plus/finim.fits"), dtype=np.float64)
            minus = np.asarray(fits.getdata(str(prefix) + "_minus/finim.fits"), dtype=np.float64)
            responses[scale] = derivative(plus, minus, stage.PLANET_CONTRAST * scale)
            midpoints[scale] = 0.5 * (plus + minus)

        for mode_index, mode in enumerate(stage.MODES):
            stamps = {scale: convergence.extract_stamp(cube[mode_index], int(site["row"]),
                                                       int(site["column"]), STAMP_SIZE)
                      for scale, cube in responses.items()}
            reference = stamps[REFERENCE_SCALE]
            reference_energy = float(np.nansum(np.square(reference[masks["full_63"]])))
            stage.require(reference_energy > 0, "reference response has zero energy")
            region_energy = {name: float(np.nansum(np.square(reference[mask]))) / reference_energy
                             for name, mask in masks.items()}
            comparisons = {
                str(scale): {name: comparison_metrics(stamps[scale], reference, mask)
                             for name, mask in masks.items()}
                for scale in CONTRAST_SCALES
            }
            richardson = (4 * stamps[0.5] - reference) / 3
            edge = {str(scale): convergence.one_response_metrics(
                        convergence.extract_stamp(responses[scale][mode_index], int(site["row"]),
                                                  int(site["column"]), PROMOTED_SIZE), PROMOTED_SIZE)
                    for scale in (0.5, 1.0, 2.0)}
            half = PROMOTED_SIZE // 2
            richardson_47 = richardson[STAMP_SIZE // 2 - half:STAMP_SIZE // 2 + half + 1,
                                       STAMP_SIZE // 2 - half:STAMP_SIZE // 2 + half + 1]
            edge["richardson"] = convergence.one_response_metrics(richardson_47, PROMOTED_SIZE)
            midpoint_rms = {}
            for scale, midpoint in midpoints.items():
                residual = convergence.extract_stamp(midpoint[mode_index] - baseline[mode_index],
                                                     int(site["row"]), int(site["column"]), PROMOTED_SIZE)
                finite = residual[np.isfinite(residual)]
                midpoint_rms[str(scale)] = float(np.sqrt(np.mean(np.square(finite))))
            measurements.append({"site": site["name"], "nominal_radius": site["nominal_radius"],
                                 "site_index": site["site_index"], "mode": mode,
                                 "reference_region_energy_fraction": region_energy,
                                 "comparisons": comparisons, "edge_metrics_47": edge,
                                 "midpoint_baseline_rms_47": midpoint_rms})
        print(f"analyzed {site_number}/{len(protocol['sites'])}: {site['name']}", flush=True)

    summaries = {}
    for mode in stage.MODES:
        mode_summary = {}
        for radius in RADII:
            rows = [row for row in measurements
                    if row["mode"] == mode and float(row["nominal_radius"]) == radius]
            comparison_summary = {}
            for scale in CONTRAST_SCALES:
                region_summary = {}
                for name in REGIONS:
                    region_summary[name] = {
                        key: summarize([row["comparisons"][str(scale)][name][key] for row in rows])
                        for key in ("support_fraction", "cosine", "candidate_to_reference_projection",
                                    "best_scaled_relative_residual", "relative_difference", "energy_ratio")
                    }
                comparison_summary[str(scale)] = region_summary
            edge_summary = {}
            for scale in ("0.5", "1.0", "2.0", "richardson"):
                edge_summary[scale] = {
                    key: summarize([row["edge_metrics_47"][scale][key] for row in rows])
                    for key in ("support_fraction", "border_energy_fraction", "negative_energy_fraction")
                }
            mode_summary[str(radius)] = {
                "sites": len(rows),
                "reference_region_energy_fraction": {
                    name: summarize([row["reference_region_energy_fraction"][name] for row in rows])
                    for name in REGIONS
                },
                "comparisons": comparison_summary,
                "edge_metrics_47": edge_summary,
                "midpoint_baseline_rms_47": {
                    str(scale): summarize([row["midpoint_baseline_rms_47"][str(scale)] for row in rows])
                    for scale in (0.5, 1.0, 2.0)
                },
            }
        summaries[str(mode)] = mode_summary

    gate = protocol["linearity_gate"]
    controlling = gate["controlling_region"]
    gate_groups = [summaries["200"][str(radius)]["comparisons"][str(scale)][controlling]
                   for radius in RADII for scale in CONTRAST_SCALES]
    cosine_passed = all(group["cosine"]["median"] >= gate["minimum_median_cosine"]
                        for group in gate_groups)
    projection_passed = all(abs(group["candidate_to_reference_projection"]["median"] - 1) <=
                            gate["projection_absolute_tolerance"] for group in gate_groups)
    edge_groups = [summaries["200"][str(radius)]["edge_metrics_47"]["richardson"] for radius in RADII]
    edge_passed = (all(group["support_fraction"]["minimum"] == 1 for group in edge_groups) and
                   all(group["border_energy_fraction"]["median"] <=
                       protocol["edge_thresholds"]["median_fraction"] for group in edge_groups) and
                   all(group["border_energy_fraction"]["maximum"] <=
                       protocol["edge_thresholds"]["individual_fraction"] for group in edge_groups))
    tail_linearity_passed = cosine_passed and projection_passed
    promote_47 = tail_linearity_passed and edge_passed
    recorded_elapsed = [read(root / "reductions" / task["name"] / "complete.json").get("elapsed_seconds")
                        for task in protocol["tasks"]]
    elapsed = [float(value) for value in recorded_elapsed if value is not None]
    stage.require(elapsed, "no reduction durations were recorded")
    result = {"tail_linearity_passed": tail_linearity_passed,
              "cosine_gate_passed": cosine_passed,
              "projection_gate_passed": projection_passed,
              "richardson_edge_gate_passed": edge_passed,
              "promote_47_pixel_response": promote_47,
              "elapsed_seconds": {"recorded_reductions": len(elapsed), "total": float(sum(elapsed)),
                                  "median_per_reduction": float(np.median(elapsed)),
                                  "maximum_per_reduction": float(max(elapsed))},
              "summaries": summaries, "measurements": measurements}
    stage.write_json(root / "results.json", result)
    write_report(root, protocol, result)
    completion = {"status": "complete", "results": stage.fingerprint(root / "results.json"),
                  "report": stage.fingerprint(root / "results.md"),
                  "reduction_receipts": [stage.fingerprint(root / "reductions" / task["name"] / "complete.json")
                                         for task in protocol["tasks"]]}
    stage.write_json(root / "complete.json", completion)
    stage.write_json(root / "state.json", {"status": "complete",
                                           "completed_tasks": [task["name"] for task in protocol["tasks"]]})
    print(root / "results.md", flush=True)


def write_report(root: Path, protocol: dict[str, object], result: dict[str, object]) -> None:
    """Write the compact contrast-linearity result report."""
    lines = ["# KLIP response-tail contrast linearity", "",
             f"The study used {len(protocol['sites'])} fixed inner sites and "
             f"{len(protocol['tasks'])} new signal-free KLIP reductions. The completed convergence run supplies "
             "the reference derivative at unit contrast scale.", "", "## Gates", "",
             "| Check | Result |", "| --- | --- |",
             f"| Mode-200 31-to-47 shell cosine | {'pass' if result['cosine_gate_passed'] else 'fail'} |",
             f"| Mode-200 31-to-47 shell projection | {'pass' if result['projection_gate_passed'] else 'fail'} |",
             f"| Richardson-extrapolated 47-pixel edge | "
             f"{'pass' if result['richardson_edge_gate_passed'] else 'fail'} |",
             f"| Promote 47-pixel response | {'yes' if result['promote_47_pixel_response'] else 'no'} |",
             "", "## Mode-200 tail diagnostics", "",
             "Values are medians across six alternating angular sites. Projection is the candidate derivative "
             "onto the unit-scale derivative. Tail energy is the unit-scale 31-to-47 shell as a fraction of "
             "the 63-pixel response energy.", "",
             "| Radius | Scale | Tail cosine | Tail projection | Tail scaled residual | Tail energy fraction |",
             "| ---: | ---: | ---: | ---: | ---: | ---: |"]
    for radius in RADII:
        summary = result["summaries"]["200"][str(radius)]
        tail_energy = summary["reference_region_energy_fraction"]["tail_31_to_47"]["median"]
        for scale in CONTRAST_SCALES:
            comparison = summary["comparisons"][str(scale)]["tail_31_to_47"]
            lines.append(f"| {radius:g} | {scale:g} | {comparison['cosine']['median']:.4f} | "
                         f"{comparison['candidate_to_reference_projection']['median']:.4f} | "
                         f"{comparison['best_scaled_relative_residual']['median']:.4f} | "
                         f"{tail_energy:.4f} |")
    lines.extend(["", "## Runtime", "",
                  f"The {result['elapsed_seconds']['recorded_reductions']} reductions with recorded durations used "
                  f"{result['elapsed_seconds']['total'] / 60:.2f} summed task-minutes; median task time was "
                  f"{result['elapsed_seconds']['median_per_reduction']:.2f} seconds.", ""])
    if result["promote_47_pixel_response"]:
        lines.append("The outer response shell is repeatable across perturbation contrasts, and the "
                     "Richardson-extrapolated 47-pixel response clears the original edge thresholds. "
                     "A 47-pixel full-field response campaign may proceed.")
    else:
        lines.append("The 47-pixel response is not promoted. Resolve the failed tail or edge gate before a "
                     "full-field response campaign.")
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def check() -> None:
    """Run deterministic checks of masks, comparisons, and Richardson arithmetic."""
    masks = region_masks()
    shell_union = (masks["core_11"] | masks["middle_11_to_31"] |
                   masks["tail_31_to_47"] | masks["outer_47_to_63"])
    stage.require(np.array_equal(shell_union, masks["full_63"]), "response shells do not cover the stamp")
    shell_sum = sum(mask.astype(int) for mask in
                    (masks["core_11"], masks["middle_11_to_31"],
                     masks["tail_31_to_47"], masks["outer_47_to_63"]))
    stage.require(np.all(shell_sum == 1), "response shells overlap")
    reference = np.arange(STAMP_SIZE * STAMP_SIZE, dtype=np.float64).reshape(STAMP_SIZE, STAMP_SIZE) + 1
    metrics = comparison_metrics(1.1 * reference, reference, masks["tail_31_to_47"])
    stage.require(np.isclose(metrics["cosine"], 1) and
                  np.isclose(metrics["candidate_to_reference_projection"], 1.1) and
                  np.isclose(metrics["best_scaled_relative_residual"], 0, atol=1e-15),
                  "response comparison arithmetic failed")
    half_contrast = 0.25
    signal = np.ones((2, 3, 3), dtype=np.float64)
    offset = np.arange(18, dtype=np.float64).reshape(2, 3, 3)
    calculated = derivative(offset + half_contrast * signal, offset - half_contrast * signal, half_contrast)
    stage.require(np.allclose(calculated, signal), "central-difference arithmetic failed")
    serialized_double_contrast = 0.00914873
    stage.require(math.isclose(serialized_double_contrast, 2 * stage.PLANET_CONTRAST,
                               rel_tol=1e-6, abs_tol=1e-12),
                  "serialized contrast tolerance does not cover production FITS precision")
    print("KLIP response-tail linearity checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the contrast-linearity command-line parser."""
    repo = Path(__file__).resolve().parents[3]
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--convergence", type=Path,
                                default=repo / "working/roc/klip_response_stamp_convergence_20260921")
    repair_parser = subparsers.add_parser("repair")
    repair_parser.add_argument("root", type=Path)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    run_parser.add_argument("--maximum-tasks", type=int, default=0,
                            help="debug/resume limit; zero runs every remaining task")
    return result


def main() -> None:
    """Dispatch the requested contrast-linearity action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    elif args.action == "prepare":
        prepare(args)
    elif args.action == "repair":
        repair(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
