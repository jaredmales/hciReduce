#!/usr/bin/env python3
"""Freeze KLIP Stage-C sites, source levels, methods, and reduction commands."""
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
from scipy.ndimage import gaussian_filter
from scipy.signal import convolve2d

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_response_47 as response47  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402
import run_klip_stage_b_local_noise_screen as raw  # noqa: E402
import run_klip_stage_b_local_radial_normalization as normalized  # noqa: E402


RADII = tuple(stage.PRIMARY_RADII)
TARGET_SNRS = (3.0, 5.0, 7.0)
ROLE_COUNTS = {"calibration": 20, "development": 6, "validation": 6, "heldout_null": 6}
TRAINING_WIDTHS = (5, 10, 20, 40, 60)
RESPONSE_SUPPORT = 11
GAUSSIAN_FWHM = 3.6
PRIMARY_MODE = 200
METHODS = (
    "native",
    "gaussian_fwhm2p4",
    "gaussian_fwhm3p0",
    "gaussian_fwhm3p6",
    "gaussian_fwhm4p2",
    "exact_identity",
    "sparse_identity",
    "exact_identity_lpf1p8",
    "exact_identity_lpf2p7",
    "exact_isotropic_fitted_mean",
    "raw_rectangular_m0p3",
    "radial_hann_m0p1_full",
    "radial_hann_m0p1_clip0p5",
    "radial_hann_m0p1_clip0p75",
    "radial_hann_m0p1_clip1p0",
    "radial_hann_m0p1_trunc0p5",
    "radial_hann_m0p1_trunc0p75",
    "radial_hann_m0p1_trunc1p0",
)


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def radius_tag(radius: float) -> str:
    """Return one path-safe radius label."""
    return format(radius, "g").replace(".", "p")


def role_tag(role: str) -> str:
    """Return a compact immutable role label."""
    return {"calibration": "cal", "development": "dev", "validation": "val",
            "heldout_null": "hold"}[role]


def gaussian_map(image: np.ndarray, fwhm: float = GAUSSIAN_FWHM) -> np.ndarray:
    """Apply the validated 15-pixel mask-normalized Gaussian reference."""
    sigma = fwhm / math.sqrt(8 * math.log(2))
    yy, xx = np.indices((15, 15), dtype=np.float64)
    kernel = np.exp(-((xx - 7) ** 2 + (yy - 7) ** 2) / (2 * sigma**2))
    finite = np.isfinite(image)
    numerator = convolve2d(np.where(finite, image, 0), kernel, mode="same")
    denominator = convolve2d(finite.astype(np.float64), kernel, mode="same")
    result = np.full(image.shape, np.nan, dtype=np.float64)
    np.divide(numerator, denominator, out=result, where=finite & (denominator > 0))
    return result


def annular_sigma(image: np.ndarray, site: dict[str, object],
                  planet: tuple[float, float]) -> tuple[float, dict[str, object]]:
    """Interpolate a one-pixel annular deviation with planet and trial exclusions."""
    native_column, native_row = np.indices(image.shape, dtype=np.float64)
    center_row = 0.5 * (image.shape[1] - 1)
    center_column = 0.5 * (image.shape[0] - 1)
    radii = np.hypot(native_row - center_row, native_column - center_column)
    usable = np.isfinite(image)
    for source_row, source_column in (planet, (float(site["row"]), float(site["column"]))):
        usable &= (np.hypot(native_row - source_row, native_column - source_column) >
                   stage.PLANET_EXCLUSION_RADIUS + 0.5)
    centers, deviations, counts = [], [], []
    for lower in range(int(math.ceil(float(np.max(radii))))):
        selected = usable & (radii >= lower) & (radii <= lower + 1)
        values = np.asarray(image[selected], dtype=np.float64)
        if len(values) < 2:
            continue
        centers.append(lower + 0.5)
        deviations.append(float(np.std(values, ddof=1)))
        counts.append(len(values))
    stage.require(len(centers) >= 2 and all(np.isfinite(deviations)) and min(deviations) >= 0,
                  "Gaussian annular profile is incomplete")
    sigma = float(np.interp(float(site["separation"]), centers, deviations))
    stage.require(np.isfinite(sigma) and sigma > 0, "Gaussian annular deviation is not positive")
    required = [count for center, count in zip(centers, counts)
                if abs(center - float(site["separation"])) <= 1.5]
    return sigma, {"minimum_local_pixels": min(required), "profile_bins": len(centers)}


def gaussian_response_peak(response: np.ndarray) -> float:
    """Return the source-only Gaussian maximum over the fixed five-pixel search."""
    sigma = GAUSSIAN_FWHM / math.sqrt(8 * math.log(2))
    filtered = gaussian_filter(np.asarray(response, dtype=np.float64), sigma=sigma,
                               mode="constant", cval=0.0, truncate=4.0)
    half = response.shape[0] // 2
    values = [float(filtered[half + delta_row, half + delta_column])
              for delta_row, delta_column in footprint.SEARCH_OFFSETS]
    peak = max(values)
    stage.require(np.isfinite(peak) and peak > 0, "Gaussian source-only response is not positive")
    return peak


def small_sample_correction(radius: float) -> float:
    """Return the production resolution-element correction at one radius."""
    count = 2 * math.pi * radius / stage.LAMBDA_D - 1
    stage.require(count > 0, "small-sample comparison count is not positive")
    return 1 / math.sqrt(1 + 1 / count)


def angular_distance(first: dict[str, object], second: dict[str, object]) -> float:
    """Return the wrapped angular distance between two same-radius candidates."""
    difference = abs(float(first["angle_radians"]) - float(second["angle_radians"]))
    return min(difference, 2 * math.pi - difference)


def maximin(rows: list[dict[str, object]], count: int) -> list[dict[str, object]]:
    """Choose a deterministic angular maximin subset."""
    stage.require(len(rows) >= count, "insufficient candidates for angular maximin selection")
    ordered = sorted(rows, key=lambda row: (float(row["angle_radians"]), int(row["row"]),
                                            int(row["column"])))
    best, best_key = None, None
    for seed_index, seed in enumerate(ordered):
        chosen = [seed]
        remaining = [row for index, row in enumerate(ordered) if index != seed_index]
        nearest = {id(row): angular_distance(row, seed) for row in remaining}
        while len(chosen) < count:
            choice = max(remaining, key=lambda row: (nearest[id(row)],
                         -float(row["angle_radians"]), -int(row["row"]), -int(row["column"])))
            chosen.append(choice)
            remaining.remove(choice)
            for row in remaining:
                nearest[id(row)] = min(nearest[id(row)], angular_distance(row, choice))
        separations = sorted(angular_distance(first, second)
                             for index, first in enumerate(chosen) for second in chosen[index + 1:])
        key = tuple(separations)
        if best_key is None or key > best_key:
            best, best_key = list(chosen), key
    return sorted(best, key=lambda row: (float(row["angle_radians"]), int(row["row"]),
                                         int(row["column"])))


def partition_roles(candidates: list[dict[str, object]], radius: float) -> list[dict[str, object]]:
    """Assign disjoint calibration, development, validation, and held-out roles."""
    available = list(candidates)
    selected = []
    for role in ("calibration", "development", "validation", "heldout_null"):
        chosen = maximin(available, ROLE_COUNTS[role])
        for index, row in enumerate(chosen):
            row = dict(row)
            row.update({"role": role, "role_index": index,
                        "name": f"r{radius_tag(radius)}_{role_tag(role)}{index:02d}"})
            selected.append(row)
        chosen_coordinates = {(int(row["row"]), int(row["column"])) for row in chosen}
        available = [row for row in available
                     if (int(row["row"]), int(row["column"])) not in chosen_coordinates]
    stage.require(len(selected) == sum(ROLE_COUNTS.values()) and
                  len({row["name"] for row in selected}) == len(selected),
                  "site partition is incomplete")
    return selected


def response_complete(paths: dict[str, object]) -> tuple[np.ndarray, list[dict[str, object]]]:
    """Require complete 47-pixel response support in every KL mode."""
    complete = None
    diagnostics = []
    for mode, path in zip(stage.MODES, paths["validities"]):
        validity = np.asarray(fits.getdata(path, memmap=True), dtype=np.float64) > 0.5
        current = np.all(validity, axis=(1, 2))
        complete = current if complete is None else complete & current
        diagnostics.append({"mode": mode, "complete_locations": int(np.count_nonzero(current))})
    stage.require(complete is not None and np.any(complete), "no complete exact responses remain")
    return complete, diagnostics


def all_mode_data_complete(baseline: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    """Mark locations with a finite 11-pixel data stamp in every KL mode."""
    finite = np.all(np.isfinite(baseline), axis=0)
    half = RESPONSE_SUPPORT // 2
    result = np.zeros(len(coordinates), dtype=bool)
    for index, (row_value, column_value, _, _) in enumerate(coordinates):
        row, column = int(row_value), int(column_value)
        if (column - half < 0 or column + half >= finite.shape[0] or
                row - half < 0 or row + half >= finite.shape[1]):
            continue
        result[index] = bool(np.all(finite[column - half:column + half + 1,
                                             row - half:row + half + 1]))
    return result


def radial_profile_supported(finite: np.ndarray, searches: list[tuple[int, int]],
                             planet: tuple[float, float]) -> tuple[bool, int]:
    """Check strict radial-profile counts without reading science values."""
    excluded = normalized.exclusion_mask(finite.shape, searches, RESPONSE_SUPPORT, planet)
    counts = normalized.profile_counts(np.where(finite, 0.0, np.nan), excluded)
    return bool(np.all(counts >= normalized.MINIMUM_PROFILE_PIXELS)), int(np.min(counts))


def split_supported_widths(finite: np.ndarray, searches: list[tuple[int, int]],
                           planet: tuple[float, float]) -> list[int]:
    """Return fixed bands supported by every query and detector half."""
    supported = set(TRAINING_WIDTHS)
    for query in searches:
        rings = footprint.training_rings(finite, query, searches, 11, planet, exclusion_support=11)
        counts = footprint.band_counts(rings, 11)
        supported &= {width for width in TRAINING_WIDTHS
                      if counts[str(width)]["split_minimum_passed"]}
    return sorted(supported)


def candidate_geometry(baseline: np.ndarray, coordinates: np.ndarray, response_ok: np.ndarray,
                       planet: tuple[float, float]) -> tuple[dict[str, list[dict[str, object]]],
                                                             dict[str, object]]:
    """Enumerate common-support candidates without reading a score."""
    data_ok = all_mode_data_complete(baseline, coordinates)
    finite = np.all(np.isfinite(baseline), axis=0)
    lookup = {(int(row), int(column)): index
              for index, (row, column, _, _) in enumerate(coordinates)}
    center_row = 0.5 * (baseline.shape[2] - 1)
    center_column = 0.5 * (baseline.shape[1] - 1)
    candidates, diagnostics = {}, {}
    for radius in RADII:
        rows = []
        for source, (row_value, column_value, _, _) in enumerate(coordinates):
            row, column = int(row_value), int(column_value)
            actual = math.hypot(row - center_row, column - center_column)
            if abs(actual - radius) > 0.5:
                continue
            searches = [(row + delta_row, column + delta_column)
                        for delta_row, delta_column in footprint.SEARCH_OFFSETS]
            indices = [lookup.get(query) for query in searches]
            if any(index is None or not response_ok[int(index)] or not data_ok[int(index)]
                   for index in indices):
                continue
            if any(math.hypot(query[0] - planet[0], query[1] - planet[1]) <=
                   stage.PLANET_EXCLUSION_RADIUS for query in searches):
                continue
            profile_ok, minimum_profile = radial_profile_supported(finite, searches, planet)
            if not profile_ok:
                continue
            widths = split_supported_widths(finite, searches, planet)
            if not widths:
                continue
            delta_row, delta_column = row - center_row, column - center_column
            rows.append({"row": row, "column": column, "source_index": source,
                         "nominal_radius": radius, "separation": actual,
                         "position_angle": math.degrees(math.atan2(-delta_row, delta_column)) % 360,
                         "angle_radians": math.atan2(delta_column, delta_row),
                         "supported_training_widths": widths,
                         "minimum_radial_profile_pixels": minimum_profile})
        selected_width = next((width for width in TRAINING_WIDTHS
                               if sum(width in row["supported_training_widths"] for row in rows) >=
                               sum(ROLE_COUNTS.values())), None)
        stage.require(selected_width is not None,
                      f"radius {radius:g} cannot support all 38 Stage-C roles")
        supported = [row for row in rows if selected_width in row["supported_training_widths"]]
        candidates[str(radius)] = supported
        diagnostics[str(radius)] = {"common_candidates_by_training_width": {
            str(width): sum(width in row["supported_training_widths"] for row in rows)
            for width in TRAINING_WIDTHS}, "selected_training_half_width": selected_width,
            "selected_width_candidates": len(supported),
            "minimum_radial_profile_pixels": min(row["minimum_radial_profile_pixels"]
                                                  for row in supported)}
        print(f"geometry radius {radius:g}: {len(supported)} candidates at band {selected_width}",
              flush=True)
    return candidates, diagnostics


def calibrate_contrasts(baseline: np.ndarray, response_path: Path,
                        sites: list[dict[str, object]], planet: tuple[float, float]) -> dict[str, object]:
    """Freeze site-specific Gaussian source-only SNR 3/5/7 contrasts."""
    mode_index = stage.MODES.index(PRIMARY_MODE)
    gaussian = gaussian_map(np.asarray(baseline[mode_index], dtype=np.float64))
    response = np.asarray(fits.getdata(response_path, memmap=True), dtype=np.float64)
    records = {}
    for site in sites:
        if site["role"] not in {"development", "validation"}:
            continue
        sigma, profile = annular_sigma(gaussian, site, planet)
        source = int(site["source_index"])
        peak = gaussian_response_peak(response[source].T)
        correction = small_sample_correction(float(site["separation"]))
        unit_snr = peak / sigma * correction
        stage.require(np.isfinite(unit_snr) and unit_snr > 0, "invalid Gaussian source-only unit SNR")
        contrasts = [target / unit_snr for target in TARGET_SNRS]
        records[str(site["name"])] = {
            "mode": PRIMARY_MODE, "gaussian_fwhm_pixels": GAUSSIAN_FWHM,
            "annular_deviation": sigma, "small_sample_correction": correction,
            "unit_contrast_five_pixel_peak": peak, "unit_contrast_source_only_snr": unit_snr,
            "target_snrs": list(TARGET_SNRS), "contrasts": contrasts,
            "trial_and_planet_excluded_from_annular_profile": True,
            "annular_profile": profile}
    stage.require(len(records) == len(RADII) * (ROLE_COUNTS["development"] +
                                                ROLE_COUNTS["validation"]),
                  "contrast calibration lost injection sites")
    return records


def task_command(root: Path, protocol: dict[str, object], task: dict[str, object]) -> list[str]:
    """Construct one frozen signal-free positive-injection reduction command."""
    site = next(row for row in protocol["sites"] if row["name"] == task["site"])
    paths = protocol["paths"]
    return [str(paths["klipreduce"]), "--config", str(root / "reduction.conf"),
            "--input.directory=", "--input.fileList", str(root / "inputs.txt"),
            "--klip.Nmodes", ",".join(map(str, stage.MODES)),
            "--psfResponse.file=", "--psfResponse.outputModels=false", "--psfResponse.filter=false",
            "--planet.sep", str(stage.PLANET_SEPARATION), "--planet.PA", str(stage.PLANET_PA),
            "--planet.contrast", str(stage.PLANET_CONTRAST), "--fake.method", "single",
            "--fake.fileName", str(paths["psf"]), "--fake.sep", str(site["separation"]),
            "--fake.PA", str(site["position_angle"]), "--fake.contrast", str(task["contrast"]),
            "--fake.subtractPlanet=true", "--output.directory", str(root / "reductions" / task["name"]),
            "--output.fileName", "finim.fits", "--output.exactFName=true", "--showTiming=true"]


def prepare(args: argparse.Namespace) -> None:
    """Freeze the complete development design before any positive reduction."""
    root, parent = args.root.resolve(), args.response.resolve()
    stage.require(not root.exists(), f"output already exists: {root}")
    parent_protocol, parent_completion = footprint.verify_parent(parent)
    stage_a_root = Path(str(parent_protocol["parent_stage_a"]))
    stage_a_protocol, stage_a_manifest = stage.load_protocol(stage_a_root)
    stage.require(sorted(os.sched_getaffinity(0)) == parent_protocol["resources"]["cpu_affinity"],
                  "prepare must use the exact-response campaign CPU affinity")
    paths = response47.product_paths(parent)
    baseline, header = fits.getdata(paths["final"], header=True)
    baseline = np.asarray(baseline, dtype=np.float64)
    stage.require(baseline.shape == (len(stage.MODES), 128, 128) and
                  stage.read_modes(header) == stage.MODES, "unexpected response-campaign baseline")
    coordinates = np.asarray(fits.getdata(paths["coordinates"]), dtype=np.float64).T
    stage.require(np.array_equal(coordinates[:, 3], np.arange(len(coordinates))),
                  "exact response coordinate indices changed")
    complete, response_diagnostics = response_complete(paths)
    planet = raw.planet_position(baseline.shape[1:], {"known_planet": {
        "separation": stage.PLANET_SEPARATION, "position_angle": stage.PLANET_PA}})
    candidates, geometry = candidate_geometry(baseline, coordinates, complete, planet)

    sites = []
    for radius in RADII:
        width = int(geometry[str(radius)]["selected_training_half_width"])
        chosen = partition_roles(candidates[str(radius)], radius)
        for row in chosen:
            row["training_half_width"] = width
        sites.extend(chosen)
    stage.require(len(sites) == len(RADII) * sum(ROLE_COUNTS.values()),
                  "Stage-C site partition has the wrong size")
    contrasts = calibrate_contrasts(baseline, paths["responses"][stage.MODES.index(PRIMARY_MODE)],
                                    sites, planet)
    development_tasks, validation_tasks = [], []
    for site in sites:
        if site["role"] not in {"development", "validation"}:
            continue
        calibration = contrasts[site["name"]]
        destination = development_tasks if site["role"] == "development" else validation_tasks
        for level_index, (target, contrast) in enumerate(zip(TARGET_SNRS, calibration["contrasts"])):
            destination.append({"name": f"{site['name']}_snr{target:g}".replace(".", "p"),
                                "site": site["name"], "role": site["role"],
                                "level_index": level_index, "target_source_snr": target,
                                "contrast": contrast})
    stage.require(len(development_tasks) == len(validation_tasks) == 108,
                  "Stage-C task count changed")

    root.mkdir(parents=True)
    (root / "software").mkdir()
    runner = Path(__file__).resolve()
    dependencies = [runner, runner.with_name("run_klip_covariance_stage_a.py"),
                    runner.with_name("check_klip_stage_b_psd_extension.py"),
                    runner.with_name("run_klip_response_47.py"),
                    runner.with_name("run_klip_response_tail_linearity.py"),
                    runner.with_name("run_klip_response_stamp_convergence.py"),
                    runner.with_name("run_klip_stage_b_footprint_preflight.py"),
                    runner.with_name("run_klip_stage_b_local_noise_screen.py"),
                    runner.with_name("run_klip_stage_b_local_radial_normalization.py")]
    for path in dependencies:
        shutil.copy2(path, root / "software" / path.name)
    shutil.copy2(parent / "reduction.conf", root / "reduction.conf")
    shutil.copy2(parent / "inputs.txt", root / "inputs.txt")

    protocol = {
        "schema": 1, "stage": "KLIP covariance matched-filter Stage C development",
        "purpose": "freeze score-blind sites and Gaussian-SNR 3/5/7 positive injections",
        "parent_response": str(parent), "parent_stage_a": str(stage_a_root),
        "modes": stage.MODES, "primary_mode": PRIMARY_MODE, "radii": list(RADII),
        "target_source_snrs": list(TARGET_SNRS), "sites_per_radius_by_role": ROLE_COUNTS,
        "site_selection": "sequential deterministic angular maximin: calibration, development, validation, held-out null",
        "site_eligibility": "complete 47-pixel exact response in all modes for all five search pixels; finite all-mode 11-pixel data; strict radial-profile support; split Welch support; planet-clear centers",
        "response_support": RESPONSE_SUPPORT,
        "training": {"patch_support": 11, "candidate_exclusion_support": 11,
                     "candidate_specific_half_width_by_radius": {
                         str(radius): geometry[str(radius)]["selected_training_half_width"] for radius in RADII},
                     "minimum_patches_per_detector_half": footprint.MINIMUM_SPLIT_PATCHES},
        "radial_standardization": "strict 3.6-pixel leave-site-out profile; known planet and union of five complete 11-pixel candidate footprints excluded",
        "known_planet": {"separation": stage.PLANET_SEPARATION,
                         "position_angle": stage.PLANET_PA, "contrast": stage.PLANET_CONTRAST,
                         "exclusion_radius": stage.PLANET_EXCLUSION_RADIUS},
        "search_offsets_row_column": [list(value) for value in footprint.SEARCH_OFFSETS],
        "contrast_calibration": "site-specific mode-200 Gaussian-FWHM-3.6 source-only five-pixel peak divided by a one-pixel annular deviation; trial and known planet excluded; production small-sample correction",
        "methods": list(METHODS),
        "method_contract": {
            "permanent_references": ["native", "gaussian_fwhm3p6", "exact_identity",
                                     "sparse_identity", "exact_identity_lpf1p8",
                                     "exact_identity_lpf2p7", "exact_isotropic_fitted_mean"],
            "development_gaussian_widths": [2.4, 3.0, 3.6, 4.2],
            "mandatory_p4_prior": "raw 11-pixel rectangular Welch PSD, mixing 0.3, no candidate fitted mean",
            "leading_covariance": "strict-radial 11-pixel Hann Welch PSD, mixing 0.1, no candidate fitted mean",
            "leading_precision_grid": {"full_inverse": True, "clip_below_mean_fractions": [0.5, 0.75, 1.0],
                                       "truncate_below_mean_fractions": [0.5, 0.75, 1.0]},
            "positive_covariance_arms": ["refit on each positive with complete source exclusion",
                                         "baseline-frozen weights replayed unchanged"]},
        "threshold_rule": "maximum of 20 fixed calibration-site five-pixel scores per radius and method; strict exceedance; freeze before positive analysis",
        "sites": sites, "development_tasks": development_tasks,
        "validation_tasks_unopened": validation_tasks,
        "expected_development_reductions": len(development_tasks),
        "expected_future_validation_reductions": len(validation_tasks),
        "resources": parent_protocol["resources"],
        "paths": {"klipreduce": parent_protocol["paths"]["klipreduce"],
                  "psf": parent_protocol["paths"]["psf"],
                  "baseline": str(paths["final"]), "exact_manifest": str(paths["manifest"]),
                  "sparse_manifest": stage_a_protocol["paths"]["sparse_manifest"]}}
    stage.write_json(root / "protocol.json", protocol)
    stage.write_json(root / "geometry.json", {"scores_or_recovery_used": False,
                                               "response_completeness": response_diagnostics,
                                               "radii": geometry, "sites": sites})
    stage.write_json(root / "contrasts.json", contrasts)
    commands = [{"task": task["name"], "command": task_command(root, protocol, task)}
                for task in development_tasks]
    stage.write_json(root / "commands.json", commands)
    frozen = [stage.fingerprint(root / name) for name in
              ("reduction.conf", "inputs.txt", "protocol.json", "geometry.json",
               "contrasts.json", "commands.json")]
    frozen.extend(stage.fingerprint(root / "software" / path.name) for path in dependencies)
    response_receipt = read(parent / "response" / "complete.json")
    lineage = [stage.fingerprint(parent / name) for name in
               ("protocol.json", "manifest.json", "results.json", "complete.json",
                "response/complete.json")]
    lineage.extend(stage.fingerprint(stage_a_root / name) for name in
                   ("protocol.json", "manifest.json", "results.json", "complete.json"))
    external = [stage.fingerprint(Path(protocol["paths"][name]))
                for name in ("klipreduce", "psf", "baseline", "exact_manifest", "sparse_manifest")]
    stage.write_json(root / "manifest.json", {
        "schema": 1, "frozen_records": frozen, "lineage_records": lineage,
        "response_products": response_receipt["products"], "external_records": external,
        "stage_a_response_products": stage_a_manifest["response_records"],
        "parent_completion": parent_completion})
    stage.write_json(root / "state.json", {"status": "prepared", "development_reductions": 0,
                                           "validation_products_opened": False,
                                           "positive_analysis_started": False})
    print(root, flush=True)


def load_experiment(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Load a prepared Stage-C experiment and verify every frozen input."""
    protocol = read(root / "protocol.json")
    manifest = read(root / "manifest.json")
    stage.verify(manifest["frozen_records"] + manifest["lineage_records"] +
                 manifest["response_products"] + manifest["external_records"] +
                 manifest["stage_a_response_products"])
    parent_protocol, completion = footprint.verify_parent(Path(str(protocol["parent_response"])))
    stage.require(completion == manifest["parent_completion"],
                  "exact-response parent completion receipt changed")
    stage.require(parent_protocol["resources"] == protocol["resources"],
                  "exact-response resource contract changed")
    stage.require(sorted(os.sched_getaffinity(0)) == protocol["resources"]["cpu_affinity"],
                  "CPU affinity changed from the frozen protocol")
    return protocol, manifest


def check() -> None:
    """Check angular selection, Gaussian filtering, contrast algebra, and commands."""
    rows = [{"row": index, "column": 0, "angle_radians": 2 * math.pi * index / 48}
            for index in range(48)]
    first = maximin(rows, 20)
    second = maximin(rows, 20)
    stage.require([(row["row"], row["column"]) for row in first] ==
                  [(row["row"], row["column"]) for row in second] and len(first) == 20,
                  "angular maximin selection is not deterministic")
    constant = np.full((128, 128), 4.25)
    constant[:2] = np.nan
    smoothed = gaussian_map(constant)
    stage.require(np.allclose(smoothed[np.isfinite(smoothed)], 4.25, rtol=0, atol=2e-14),
                  "mask-normalized Gaussian does not preserve a constant")
    radius = 12.0
    correction = small_sample_correction(radius)
    target, sigma, peak = 5.0, 0.3, 1.2
    contrast = target * sigma / correction / peak
    stage.require(np.isclose(contrast * peak / sigma * correction, target),
                  "source-only contrast calibration algebra failed")
    protocol = {"sites": [{"name": "site", "separation": 10.2, "position_angle": 31.0}],
                "paths": {"klipreduce": "/bin/true", "psf": "/tmp/psf.fits"}}
    command = task_command(Path("/tmp/stage_c"), protocol,
                           {"name": "job", "site": "site", "contrast": 0.001})
    joined = " ".join(command)
    for token in ("--fake.subtractPlanet=true", "--fake.contrast 0.001",
                  "--klip.Nmodes 125,150,175,200,225,250,300,350",
                  "--psfResponse.outputModels=false"):
        stage.require(token in joined, f"Stage-C command lacks {token}")
    stage.require(len(METHODS) == len(set(METHODS)) and len(METHODS) == 18,
                  "Stage-C method manifest changed unexpectedly")
    print("KLIP Stage-C development preparer checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-C preparation command-line parser."""
    repo = Path(__file__).resolve().parents[3]
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("root", type=Path)
    prepare_parser.add_argument("--response", type=Path,
                                default=repo / "working/roc/klip_response_47_20260921")
    return result


def main() -> None:
    """Dispatch the Stage-C preparation action."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check()
    else:
        prepare(arguments)


if __name__ == "__main__":
    main()
