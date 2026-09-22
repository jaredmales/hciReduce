#!/usr/bin/env python3
"""Generate and validate the promoted 47-pixel exact KLIP response field."""
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
import run_klip_response_tail_linearity as linearity  # noqa: E402


STAMP_SIZE = 47
ARCHIVED_STAMP_SIZE = 11
EXPECTED_LOCATIONS = 11192
EXPECTED_TRIALS = 2 * EXPECTED_LOCATIONS
ELEMENT_BYTES = np.dtype(np.float32).itemsize + np.dtype(np.uint8).itemsize
EXPECTED_RETAINED_BYTES = EXPECTED_LOCATIONS * len(stage.MODES) * STAMP_SIZE * STAMP_SIZE * ELEMENT_BYTES


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def verify_parent(parent: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Verify the completed response-tail experiment and its promotion decision."""
    protocol, _ = linearity.load_experiment(parent)
    completion = read(parent / "complete.json")
    state = read(parent / "state.json")
    results = read(parent / "results.json")
    stage.require(completion["status"] == "complete" and state["status"] == "complete",
                  "parent response-tail experiment is incomplete")
    stage.verify([completion["results"], completion["report"], *completion["reduction_receipts"]])
    for task in protocol["tasks"]:
        receipt = read(parent / "reductions" / task["name"] / "complete.json")
        stage.verify(receipt["products"])
    stage.require(results["tail_linearity_passed"] and results["richardson_edge_gate_passed"] and
                  results["promote_47_pixel_response"], "parent did not promote the 47-pixel response")
    return protocol, results


def lineage(parent: Path, parent_protocol: dict[str, object]) -> tuple[Path, dict[str, object], dict[str, object]]:
    """Resolve and verify the convergence and Stage-A ancestors."""
    convergence_root = Path(str(parent_protocol["parent_convergence"]))
    convergence_protocol = read(convergence_root / "protocol.json")
    stage_a_root = Path(str(convergence_protocol["parent_stage_a"]))
    stage_a_protocol, _ = stage.load_protocol(stage_a_root)
    stage.require(stage_a_protocol["resources"] == parent_protocol["resources"],
                  "Stage-A and response-tail resource contracts differ")
    stage.require(stage_a_protocol["binary_provenance"]["hash_matches"]["klipReduce"],
                  "Stage-A klipReduce does not match the archived exact-response binary")
    return stage_a_root, stage_a_protocol, convergence_protocol


def prepare(args: argparse.Namespace) -> None:
    """Freeze the promoted footprint, complete response command, inputs, and lineage."""
    root, parent = args.root.resolve(), args.linearity.resolve()
    stage.require(not root.exists(), f"output already exists: {root}")
    parent_protocol, _ = verify_parent(parent)
    stage_a_root, stage_a_protocol, convergence_protocol = lineage(parent, parent_protocol)
    resources = parent_protocol["resources"]
    stage.require(sorted(os.sched_getaffinity(0)) == resources["cpu_affinity"],
                  "prepare must use the parent CPU affinity")
    archived_manifest = Path(str(stage_a_protocol["paths"]["exact_manifest"]))
    archived_products = stage.exact_products(archived_manifest)
    archived_header = archived_products["header"]
    stage.require(int(archived_header["KLIP PSF MEASUREMENT COUNT"]) == EXPECTED_LOCATIONS and
                  int(archived_header["KLIP PSF REFIT TRIAL COUNT"]) == EXPECTED_TRIALS,
                  "archived response location or trial count changed")

    root.mkdir(parents=True)
    (root / "software").mkdir()
    runner = Path(__file__).resolve()
    dependencies = [runner, runner.with_name("run_klip_response_tail_linearity.py"),
                    runner.with_name("run_klip_response_stamp_convergence.py"),
                    runner.with_name("run_klip_covariance_stage_a.py")]
    for path in dependencies:
        shutil.copy2(path, root / "software" / path.name)
    shutil.copy2(stage_a_root / "reduction.conf", root / "reduction.conf")
    shutil.copy2(stage_a_root / "inputs.txt", root / "inputs.txt")

    protocol = {
        "schema": 1,
        "purpose": "generate the promoted 47-pixel exact KLIP response at every 6-to-60-pixel search location",
        "parent_linearity": str(parent),
        "parent_convergence": str(parent_protocol["parent_convergence"]),
        "parent_stage_a": str(stage_a_root),
        "modes": stage.MODES,
        "primary_mode": 200,
        "stamp_size": STAMP_SIZE,
        "search_radii": {"minimum": 6.0, "maximum": 60.0},
        "expected_locations": EXPECTED_LOCATIONS,
        "expected_trials": EXPECTED_TRIALS,
        "expected_retained_bytes": EXPECTED_RETAINED_BYTES,
        "known_planet": {"separation": stage.PLANET_SEPARATION, "position_angle": stage.PLANET_PA,
                         "contrast": stage.PLANET_CONTRAST},
        "perturbation_half_contrast": stage.PLANET_CONTRAST,
        "response_replay_gate": {"minimum_cosine": 0.999,
                                 "projection_absolute_tolerance": 0.01},
        "resources": resources,
        "paths": {"klipreduce": stage_a_protocol["paths"]["klipreduce"],
                  "psf": stage_a_protocol["paths"]["psf"],
                  "parent_baseline": str(stage_a_root / "current_baseline/finim.fits"),
                  "archived_exact_manifest": str(archived_manifest)},
    }
    stage.write_json(root / "protocol.json", protocol)
    records = [stage.fingerprint(root / name) for name in ("reduction.conf", "inputs.txt", "protocol.json")]
    records.extend(stage.fingerprint(root / "software" / path.name) for path in dependencies)
    lineage_records = []
    for ancestor in (parent, Path(str(parent_protocol["parent_convergence"])), stage_a_root):
        for name in ("protocol.json", "manifest.json", "results.json", "complete.json"):
            path = ancestor / name
            if path.is_file():
                lineage_records.append(stage.fingerprint(path))
    external_records = [stage.fingerprint(Path(protocol["paths"][name]))
                        for name in ("klipreduce", "psf", "parent_baseline", "archived_exact_manifest")]
    external_records.extend(stage.fingerprint(path) for path in stage.product_paths(archived_manifest)
                            if path != archived_manifest)
    stage.write_json(root / "manifest.json", {"schema": 1, "frozen_records": records,
                                              "lineage_records": lineage_records,
                                              "external_records": external_records})
    stage.write_json(root / "state.json", {"status": "prepared"})
    print(root, flush=True)


def load_experiment(root: Path) -> dict[str, object]:
    """Load one prepared campaign and verify every frozen input."""
    protocol = read(root / "protocol.json")
    manifest = read(root / "manifest.json")
    stage.verify(manifest["frozen_records"] + manifest["lineage_records"] + manifest["external_records"])
    parent_protocol, _ = verify_parent(Path(protocol["parent_linearity"]))
    stage_a_root, stage_a_protocol, _ = lineage(Path(protocol["parent_linearity"]), parent_protocol)
    stage.require(str(stage_a_root) == protocol["parent_stage_a"] and
                  stage_a_protocol["resources"] == protocol["resources"], "campaign lineage changed")
    stage.require(sorted(os.sched_getaffinity(0)) == protocol["resources"]["cpu_affinity"],
                  "CPU affinity changed from the frozen protocol")
    return protocol


def response_command(root: Path, protocol: dict[str, object]) -> list[str]:
    """Construct the promoted exact-response command."""
    paths = protocol["paths"]
    return [str(paths["klipreduce"]), "--config", str(root / "reduction.conf"),
            "--input.directory=", "--input.fileList", str(root / "inputs.txt"),
            "--klip.Nmodes", ",".join(map(str, stage.MODES)),
            "--planet.sep", str(stage.PLANET_SEPARATION), "--planet.PA", str(stage.PLANET_PA),
            "--planet.contrast", str(stage.PLANET_CONTRAST), "--fake.method", "single",
            "--fake.fileName", str(paths["psf"]), "--fake.subtractPlanet=true",
            "--psfResponse.file", str(paths["psf"]), "--psfResponse.stampSize", str(STAMP_SIZE),
            "--psfResponse.sampleEveryPixel=true", "--psfResponse.method", "refitDifference",
            "--psfResponse.sampleAvoidRadius", "0", "--psfResponse.refitContrast", str(stage.PLANET_CONTRAST),
            "--psfResponse.outputModels=true", "--psfResponse.filter=false",
            "--psfResponse.filterMinGoodFract", "1", "--psfResponse.outputPrefix", "klipPSF_",
            "--output.directory", str(root / "response"), "--output.fileName", "finim.fits",
            "--output.exactFName=true", "--showTiming=true"]


def product_paths(root: Path) -> dict[str, object]:
    """Return paths for the expected promoted response product set."""
    directory = root / "response" / "finim_outputs"
    manifest = directory / "klipPSF_manifest.fits"
    return {"final": root / "response" / "finim.fits", "manifest": manifest,
            "coordinates": directory / "klipPSF_coordinates.fits",
            "responses": [directory / f"klipPSF_mode{index:03d}_pixel_response.fits"
                          for index in range(len(stage.MODES))],
            "validities": [directory / f"klipPSF_mode{index:03d}_pixel_validity.fits"
                           for index in range(len(stage.MODES))]}


def validate_products(root: Path, protocol: dict[str, object]) -> dict[str, object]:
    """Validate schema, dimensions, coordinates, headers, and baseline replay."""
    paths = product_paths(root)
    required = [paths["final"], paths["manifest"], paths["coordinates"],
                *paths["responses"], *paths["validities"]]
    stage.require(all(path.is_file() for path in required), "47-pixel response product set is incomplete")
    manifest_data, header = fits.getdata(paths["manifest"], header=True)
    stage.require(np.asarray(manifest_data).size == 1 and float(np.ravel(manifest_data)[0]) == 1,
                  "47-pixel manifest sentinel is invalid")
    stage.require(int(header["KLIP PSF COMPLETE"]) == 1 and int(header["KLIP PSF PRODUCT SCHEMA"]) == 2 and
                  str(header["KLIP PSF SPATIAL MODEL"]).strip() == "PIXEL_EXACT" and
                  str(header["KLIP PSF RESPONSE METHOD"]).strip() == "refitDifference" and
                  str(header["KLIP PSF ACCUMULATION"]).strip() == "PAIRED_FINAL_DIFFERENCE" and
                  int(header["KLIP PSF SAMPLE EVERY PIXEL"]) == 1,
                  "47-pixel response manifest contract changed")
    stage.require(stage.read_modes(header) == stage.MODES and int(header["KLIP PSF STAMP SIZE"]) == STAMP_SIZE,
                  "47-pixel response modes or stamp size changed")
    stage.require(int(header["KLIP PSF MEASUREMENT COUNT"]) == protocol["expected_locations"] and
                  int(header["KLIP PSF REFIT TRIAL COUNT"]) == protocol["expected_trials"] and
                  int(header["KLIP PSF RETAINED BYTES"]) == protocol["expected_retained_bytes"],
                  "47-pixel response counts or retained storage changed")
    stage.require(math.isclose(float(header["KLIP PSF REFIT CONTRAST"]),
                               protocol["perturbation_half_contrast"], rel_tol=1e-6, abs_tol=1e-12),
                  "47-pixel response perturbation changed")

    coordinates = np.asarray(fits.getdata(paths["coordinates"]), dtype=np.float64)
    archived = stage.exact_products(Path(protocol["paths"]["archived_exact_manifest"]))
    archived_coordinates = np.asarray(fits.getdata(archived["coordinates"]), dtype=np.float64)
    stage.require(coordinates.shape == archived_coordinates.shape == (4, protocol["expected_locations"]) and
                  np.array_equal(coordinates, archived_coordinates), "47-pixel response coordinates changed")
    for mode, response_path, validity_path in zip(stage.MODES, paths["responses"], paths["validities"]):
        response = np.asarray(fits.getdata(response_path, memmap=True))
        validity = np.asarray(fits.getdata(validity_path, memmap=True))
        expected_shape = (protocol["expected_locations"], STAMP_SIZE, STAMP_SIZE)
        stage.require(response.shape == validity.shape == expected_shape,
                      f"mode {mode} 47-pixel product shape changed")
        stage.require(np.all((validity == 0) | (validity == 1)), f"mode {mode} validity is not binary")
        support = validity > 0.5
        stage.require(np.all(np.isfinite(response[support])) and
                      np.all(support[:, STAMP_SIZE // 2, STAMP_SIZE // 2]),
                      f"mode {mode} response has invalid supported values or anchors")

    final, final_header = fits.getdata(paths["final"], header=True)
    baseline, baseline_header = fits.getdata(protocol["paths"]["parent_baseline"], header=True)
    stage.require(stage.read_modes(final_header) == stage.read_modes(baseline_header) == stage.MODES and
                  final.shape == baseline.shape and np.array_equal(np.isfinite(final), np.isfinite(baseline)),
                  "47-pixel final image and parent baseline schemas differ")
    difference = np.where(np.isfinite(final), np.asarray(final, dtype=np.float64) - baseline, 0)
    baseline_bitwise = np.array_equal(final, baseline, equal_nan=True)
    return {"baseline_bitwise_identical": baseline_bitwise,
            "baseline_maximum_absolute_difference": float(np.max(np.abs(difference))),
            "baseline_rms_difference": float(np.sqrt(np.mean(np.square(difference)))),
            "coordinates_identical_to_archive": True}


def summarize(values: np.ndarray) -> dict[str, float]:
    """Summarize one finite numeric vector."""
    stage.require(values.size > 0 and np.all(np.isfinite(values)), "cannot summarize empty or nonfinite values")
    return {"minimum": float(np.min(values)), "median": float(np.median(values)),
            "maximum": float(np.max(values))}


def comparison_vectors(candidate: np.ndarray, reference: np.ndarray,
                       support: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return per-source cosine, projection, and scaled residual vectors."""
    candidate_values = np.where(support, candidate, 0).reshape(candidate.shape[0], -1)
    reference_values = np.where(support, reference, 0).reshape(reference.shape[0], -1)
    candidate_energy = np.sum(np.square(candidate_values), axis=1)
    reference_energy = np.sum(np.square(reference_values), axis=1)
    cross = np.sum(candidate_values * reference_values, axis=1)
    valid = (candidate_energy > 0) & (reference_energy > 0)
    stage.require(np.all(valid), "response replay contains a zero-energy source")
    projection = cross / reference_energy
    cosine = cross / np.sqrt(candidate_energy * reference_energy)
    residual_energy = np.sum(np.square(candidate_values - projection[:, None] * reference_values), axis=1)
    residual = np.sqrt(residual_energy / candidate_energy)
    return cosine, projection, residual


def analyze(root: Path, protocol: dict[str, object], validation: dict[str, object]) -> dict[str, object]:
    """Compare the new central core with the archive and full stamps with external derivatives."""
    paths = product_paths(root)
    archived = stage.exact_products(Path(protocol["paths"]["archived_exact_manifest"]))
    coordinates = np.asarray(fits.getdata(paths["coordinates"]), dtype=np.float64).T
    coordinate_lookup = {(int(row), int(column)): index for index, (row, column, _, _) in enumerate(coordinates)}
    parent_protocol = read(Path(protocol["parent_linearity"]) / "protocol.json")
    convergence_root = Path(str(parent_protocol["parent_convergence"]))
    external = {}
    for site in parent_protocol["sites"]:
        plus = np.asarray(fits.getdata(convergence_root / "reductions" / f"{site['name']}_plus" / "finim.fits"),
                          dtype=np.float64)
        minus = np.asarray(fits.getdata(convergence_root / "reductions" / f"{site['name']}_minus" / "finim.fits"),
                           dtype=np.float64)
        external[site["name"]] = (plus - minus) / (2 * stage.PLANET_CONTRAST)

    modes = {}
    replay_cosines, replay_projection_errors = [], []
    half_old = ARCHIVED_STAMP_SIZE // 2
    center_new = STAMP_SIZE // 2
    central_slice = slice(center_new - half_old, center_new + half_old + 1)
    for mode_index, mode in enumerate(stage.MODES):
        response = np.asarray(fits.getdata(paths["responses"][mode_index], memmap=True), dtype=np.float64)
        validity = np.asarray(fits.getdata(paths["validities"][mode_index], memmap=True), dtype=np.float64) > 0.5
        archived_response = np.asarray(fits.getdata(archived["responses"][mode_index]), dtype=np.float64)
        archived_validity = np.asarray(fits.getdata(archived["validities"][mode_index]), dtype=np.float64) > 0.5
        central_response = response[:, central_slice, central_slice]
        central_validity = validity[:, central_slice, central_slice]
        stage.require(np.array_equal(central_validity, archived_validity),
                      f"mode {mode} central validity differs from the 11-pixel archive")
        cosine, projection, residual = comparison_vectors(central_response, archived_response, archived_validity)
        central = {"cosine": summarize(cosine),
                   "projection": summarize(projection),
                   "projection_error_maximum": float(np.max(np.abs(projection - 1))),
                   "best_scaled_relative_residual": summarize(residual),
                   "maximum_absolute_difference": float(np.max(np.abs(np.where(
                       archived_validity, central_response - archived_response, 0))))}

        site_replays = []
        for site in parent_protocol["sites"]:
            source = coordinate_lookup[(int(site["row"]), int(site["column"]))]
            measured = convergence.extract_stamp(external[site["name"]][mode_index], int(site["row"]),
                                                 int(site["column"]), STAMP_SIZE)
            metrics = convergence.replay_metrics(measured, response[source].T, validity[source].T)
            site_replays.append({"site": site["name"], "nominal_radius": site["nominal_radius"], **metrics})
            replay_cosines.append(metrics["cosine"])
            replay_projection_errors.append(abs(metrics["exact_to_measured_projection"] - 1))
        modes[str(mode)] = {"central_11_archive_replay": central,
                            "external_47_replay": {
                                key: convergence.summarize([row[key] for row in site_replays])
                                for key in ("support_fraction", "cosine", "exact_to_measured_projection",
                                            "best_scaled_relative_residual", "maximum_absolute_difference")},
                            "external_47_sites": site_replays}

    gate = protocol["response_replay_gate"]
    central_passed = all(result["central_11_archive_replay"]["cosine"]["minimum"] >= gate["minimum_cosine"] and
                         result["central_11_archive_replay"]["projection_error_maximum"] <=
                         gate["projection_absolute_tolerance"] for result in modes.values())
    external_passed = (min(replay_cosines) >= gate["minimum_cosine"] and
                       max(replay_projection_errors) <= gate["projection_absolute_tolerance"])
    result = {"product_validation": validation,
              "central_11_archive_replay_passed": central_passed,
              "external_47_replay_passed": external_passed,
              "minimum_external_47_cosine": min(replay_cosines),
              "maximum_external_47_projection_error": max(replay_projection_errors),
              "modes": modes}
    stage.write_json(root / "results.json", result)
    write_report(root, protocol, result)
    return result


def write_report(root: Path, protocol: dict[str, object], result: dict[str, object]) -> None:
    """Write the promoted response campaign report."""
    validation = result["product_validation"]
    lines = ["# KLIP 47-pixel exact-response campaign", "", "## Gates", "",
             "| Check | Result | Diagnostic |", "| --- | --- | --- |",
             f"| Signal-free baseline replay | {'pass' if validation['baseline_bitwise_identical'] else 'fail'} | "
             f"maximum difference {validation['baseline_maximum_absolute_difference']:.6g} |",
             f"| Central 11 pixels versus archive | "
             f"{'pass' if result['central_11_archive_replay_passed'] else 'fail'} | all "
             f"{protocol['expected_locations']} locations and eight modes |",
             f"| Full 47 pixels versus external derivatives | "
             f"{'pass' if result['external_47_replay_passed'] else 'fail'} | minimum cosine "
             f"{result['minimum_external_47_cosine']:.6f}; maximum projection error "
             f"{result['maximum_external_47_projection_error']:.6g} |", "", "## Per-mode replay", "",
             "| Mode | Central minimum cosine | Central maximum projection error | External minimum cosine | "
             "External maximum projection error |",
             "| ---: | ---: | ---: | ---: | ---: |"]
    for mode in stage.MODES:
        mode_result = result["modes"][str(mode)]
        central = mode_result["central_11_archive_replay"]
        external = mode_result["external_47_sites"]
        lines.append(f"| {mode} | {central['cosine']['minimum']:.6f} | "
                     f"{central['projection_error_maximum']:.6g} | "
                     f"{min(row['cosine'] for row in external):.6f} | "
                     f"{max(abs(row['exact_to_measured_projection'] - 1) for row in external):.6g} |")
    lines.extend(["", f"The response contains {protocol['expected_locations']} exact locations, "
                  f"{protocol['expected_trials']} signed refit trials, and 47-by-47 stamps in all eight modes."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def finalize(root: Path, protocol: dict[str, object], elapsed: float | None,
             recovered: bool) -> None:
    """Validate, analyze, and receipt a completed native response run."""
    validation = validate_products(root, protocol)
    paths = product_paths(root)
    product_records = [stage.fingerprint(path) for path in
                       [paths["final"], paths["manifest"], paths["coordinates"],
                        *paths["responses"], *paths["validities"]]]
    response_receipt = {"elapsed_seconds": elapsed, "recovered_after_runner_interruption": recovered,
                        "products": product_records, "validation": validation}
    stage.write_json(root / "response" / "complete.json", response_receipt)
    result = analyze(root, protocol, validation)
    stage.require(validation["baseline_bitwise_identical"], "47-pixel run changed the signal-free baseline")
    stage.require(result["central_11_archive_replay_passed"], "47-pixel central response does not match archive")
    stage.require(result["external_47_replay_passed"], "47-pixel response does not match external derivatives")
    completion = {"status": "complete", "response": stage.fingerprint(root / "response" / "complete.json"),
                  "results": stage.fingerprint(root / "results.json"),
                  "report": stage.fingerprint(root / "results.md")}
    stage.write_json(root / "complete.json", completion)
    stage.write_json(root / "state.json", {"status": "complete"})
    print(root / "results.md", flush=True)


def run(args: argparse.Namespace) -> None:
    """Run the full campaign or validate and finish an existing native product."""
    root = args.root.resolve()
    protocol = load_experiment(root)
    if (root / "complete.json").is_file():
        completion = read(root / "complete.json")
        stage.verify([completion["response"], completion["results"], completion["report"]])
        response_receipt = read(root / "response" / "complete.json")
        stage.verify(response_receipt["products"])
        stage.require(read(root / "state.json")["status"] == "complete", "completion state is inconsistent")
        print(root / "results.md", flush=True)
        return

    response_root = root / "response"
    if response_root.exists():
        try:
            validate_products(root, protocol)
        except Exception as error:
            stage.write_json(root / "state.json", {"status": "incomplete", "error_type": type(error).__name__,
                                                   "error": str(error)})
            raise RuntimeError(f"incomplete response directory exists; preserve it before retry: {response_root}") \
                from error
        finalize(root, protocol, None, True)
        return

    response_root.mkdir()
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(protocol["resources"]["openmp_threads"]), OMP_PROC_BIND="true",
                       OMP_PLACES="cores", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    stage.write_json(root / "state.json", {"status": "running", "pid": os.getpid(),
                                           "progress_log": str(response_root / "run.log")})
    try:
        elapsed = stage.run_command(response_command(root, protocol), response_root,
                                    response_root / "run.log", environment)
        finalize(root, protocol, elapsed, False)
    except Exception as error:
        stage.write_json(root / "state.json", {"status": "failed", "error_type": type(error).__name__,
                                               "error": str(error),
                                               "progress_log": str(response_root / "run.log")})
        raise


def check() -> None:
    """Check storage arithmetic, command construction, and replay algebra."""
    stage.require(EXPECTED_RETAINED_BYTES == 988925120, "47-pixel retained-byte calculation changed")
    reference = np.arange(3 * 5 * 5, dtype=np.float64).reshape(3, 5, 5) + 1
    support = np.ones_like(reference, dtype=bool)
    cosine, projection, residual = comparison_vectors(1.002 * reference, reference, support)
    stage.require(np.allclose(cosine, 1) and np.allclose(projection, 1.002) and
                  np.allclose(residual, 0, atol=1e-15), "response replay arithmetic failed")
    command = response_command(Path("/tmp/klip47"),
                               {"paths": {"klipreduce": "/bin/true", "psf": "/tmp/psf.fits"}})
    joined = " ".join(command)
    for token in ("--psfResponse.stampSize 47", "--psfResponse.sampleEveryPixel=true",
                  "--psfResponse.method refitDifference", "--psfResponse.outputModels=true",
                  "--fake.subtractPlanet=true"):
        stage.require(token in joined, f"47-pixel response command lacks {token}")
    print("KLIP 47-pixel exact-response checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the exact-response campaign command-line parser."""
    repo = Path(__file__).resolve().parents[3]
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--linearity", type=Path,
                                default=repo / "working/roc/klip_response_tail_linearity_20260921")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("root", type=Path)
    return result


def main() -> None:
    """Dispatch the requested exact-response action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    elif args.action == "prepare":
        prepare(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
