#!/usr/bin/env python3
"""Run calibrated KLIP Stage-C development reductions and matched-filter analysis."""
from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_klip_stage_b_psd_extension as extension  # noqa: E402
import prepare_klip_stage_c_development as preparation  # noqa: E402
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_response_47 as response47  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402
import run_klip_stage_b_local_noise_screen as raw  # noqa: E402
import run_klip_stage_b_local_radial_normalization as radial  # noqa: E402


METHODS = tuple(preparation.METHODS)
REFERENCE_WEIGHT_METHODS = ("exact_identity", "sparse_identity",
                            "exact_identity_lpf1p8", "exact_identity_lpf2p7")
COVARIANCE_METHODS = (
    "raw_rectangular_m0p3",
    "radial_hann_m0p1_full",
    "radial_hann_m0p1_clip0p5",
    "radial_hann_m0p1_clip0p75",
    "radial_hann_m0p1_clip1p0",
    "radial_hann_m0p1_trunc0p5",
    "radial_hann_m0p1_trunc0p75",
    "radial_hann_m0p1_trunc1p0",
)
RADIAL_POLICIES = (
    ("radial_hann_m0p1_full", "full", 0.0),
    ("radial_hann_m0p1_clip0p5", "clip", 0.5),
    ("radial_hann_m0p1_clip0p75", "clip", 0.75),
    ("radial_hann_m0p1_clip1p0", "clip", 1.0),
    ("radial_hann_m0p1_trunc0p5", "truncate", 0.5),
    ("radial_hann_m0p1_trunc0p75", "truncate", 0.75),
    ("radial_hann_m0p1_trunc1p0", "truncate", 1.0),
)
GAUSSIAN_METHODS = (("gaussian_fwhm2p4", 2.4), ("gaussian_fwhm3p0", 3.0),
                    ("gaussian_fwhm3p6", 3.6), ("gaussian_fwhm4p2", 4.2))
SUPPORT = 11
HALF = SUPPORT // 2
ANALYSIS_ROLES = ("calibration", "development")


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def radius_tag(radius: float) -> str:
    """Return a path-safe radius label."""
    return format(radius, "g").replace(".", "p")


def archive_incomplete(root: Path, directory: Path, group: str) -> Path:
    """Move one unreceipted directory to a numbered interruption archive."""
    archive = root / "interrupted" / group / directory.name
    archive.mkdir(parents=True, exist_ok=True)
    attempt = 1
    while (archive / f"attempt_{attempt:04d}").exists():
        attempt += 1
    destination = archive / f"attempt_{attempt:04d}"
    shutil.move(directory, destination)
    return destination


def environment(protocol: dict[str, object], reduction: bool) -> dict[str, str]:
    """Return the frozen single-BLAS-thread ROC execution environment."""
    result = os.environ.copy()
    result.update(OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    if reduction:
        result.update(OMP_NUM_THREADS=str(protocol["resources"]["openmp_threads"]),
                      OMP_PROC_BIND="true", OMP_PLACES="cores")
    else:
        result.update(OMP_NUM_THREADS="1")
    return result


def ensure_frozen_runner(root: Path, arguments: list[str]) -> None:
    """Copy this runner into the prepared package and re-execute the frozen copy."""
    destination = root / "software" / Path(__file__).name
    if Path(__file__).resolve() == destination.resolve():
        return
    stage.require(root.is_dir() and (root / "manifest.json").is_file(),
                  "Stage-C preparation does not exist")
    if destination.exists() and destination.read_bytes() != Path(__file__).read_bytes():
        state = read(root / "state.json")
        stage.require(state["status"] == "prepared" and not (root / "reductions").exists() and
                      not (root / "calibration" / "complete.json").exists(),
                      "development runner can change only before calibration and reductions complete")
        archive = root / "interrupted" / "runner_repairs"
        archive.mkdir(parents=True, exist_ok=True)
        attempt = 1
        while (archive / f"attempt_{attempt:04d}").exists():
            attempt += 1
        repair = archive / f"attempt_{attempt:04d}"
        repair.mkdir()
        destination.replace(repair / destination.name)
        previous_runner = stage.fingerprint(repair / destination.name)
        manifest = root / "development_manifest.json"
        previous_manifest = None
        if manifest.exists():
            manifest.replace(repair / manifest.name)
            previous_manifest = stage.fingerprint(repair / manifest.name)
        stage.write_json(repair / "repair.json", {
            "reason": "pre-calibration runner correction",
            "previous_runner": previous_runner,
            "previous_development_manifest": previous_manifest,
            "replacement": stage.fingerprint(Path(__file__)),
            "calibration_units_retained": len(list((root / "calibration" / "units").glob("*/complete.json")))
            if (root / "calibration" / "units").exists() else 0,
            "completed_null_sites_retained": len(list((root / "calibration" / "sites").glob("*/complete.json")))
            if (root / "calibration" / "sites").exists() else 0})
    if not destination.exists():
        shutil.copy2(Path(__file__), destination)
    os.execv(sys.executable, [sys.executable, str(destination), *arguments])


def enable(root: Path) -> tuple[dict[str, object], dict[str, object], Path]:
    """Verify preparation and freeze the runner plus production analyzer."""
    protocol, manifest = preparation.load_experiment(root)
    stage_a_protocol, _ = stage.load_protocol(Path(str(protocol["parent_stage_a"])))
    analyzer = Path(str(stage_a_protocol["paths"]["hcianalyze"]))
    stage.require(analyzer.is_file(), "frozen hciAnalyze executable is missing")
    path = root / "development_manifest.json"
    expected = {
        "schema": 1,
        "purpose": "freeze Stage-C calibration, reduction, and development-analysis software",
        "preparation": stage.fingerprint(root / "manifest.json"),
        "runner": stage.fingerprint(Path(__file__)),
        "hcianalyze": stage.fingerprint(analyzer),
    }
    if path.exists():
        stage.require(read(path) == expected, "development runner manifest changed")
    else:
        state = read(root / "state.json")
        stage.require(state["status"] == "prepared" and not (root / "reductions").exists(),
                      "runner can only be enabled on untouched preparation")
        stage.write_json(path, expected)
    stage.verify([expected["preparation"], expected["runner"], expected["hcianalyze"]])
    return protocol, manifest, analyzer


def image_radius(shape: tuple[int, int]) -> np.ndarray:
    """Return native-pixel radius for a FITS-order image."""
    column, row = np.indices(shape, dtype=np.float64)
    return np.hypot(row - 0.5 * (shape[1] - 1), column - 0.5 * (shape[0] - 1))


def required_bins(sites: list[dict[str, object]], radius_map: np.ndarray) -> tuple[list[int], np.ndarray]:
    """Select complete one-pixel annuli bracketing all five-pixel searches."""
    bins: set[int] = set()
    for site in sites:
        for delta_row, delta_column in footprint.SEARCH_OFFSETS:
            row = int(site["row"]) + delta_row
            column = int(site["column"]) + delta_column
            lower = math.floor(float(radius_map[column, row]) - 0.5)
            bins.update((lower, lower + 1))
    selected = np.zeros(radius_map.shape, dtype=bool)
    for lower in bins:
        selected |= (radius_map > lower) & (radius_map <= lower + 1)
    return sorted(bins), selected


def candidate_mask(query: tuple[int, int], planet: tuple[float, float]) -> np.ndarray:
    """Return the 11-pixel candidate coordinates outside the fixed planet disk."""
    row, column = query
    native_row, native_column = np.mgrid[row - HALF:row + HALF + 1,
                                         column - HALF:column + HALF + 1]
    return np.hypot(native_row - planet[0], native_column - planet[1]) > stage.PLANET_EXCLUSION_RADIUS


def stamp(image: np.ndarray, query: tuple[int, int]) -> np.ndarray:
    """Extract one Eigen-oriented 11-pixel stamp from a FITS-order image."""
    row, column = query
    result = np.asarray(image[column - HALF:column + HALF + 1,
                              row - HALF:row + HALF + 1], dtype=np.float64).T
    stage.require(result.shape == (SUPPORT, SUPPORT), "candidate stamp is incomplete")
    return result


def normalized_weight(template: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Return one masked identity-covariance weight with unit template response."""
    vector = np.asarray(template, dtype=np.float64).ravel()
    selected = np.asarray(mask, dtype=bool).ravel()
    energy = float(vector[selected] @ vector[selected])
    stage.require(np.isfinite(energy) and energy > 0, "identity template has no energy")
    weight = np.zeros(vector.shape, dtype=np.float64)
    weight[selected] = vector[selected] / energy
    return weight


def combined_samples(image: np.ndarray, rings: dict[int, dict[str, object]],
                     width: int) -> tuple[np.ndarray, list[int]]:
    """Combine the two separately supported detector-half training matrices."""
    halves = [raw.band_samples(image, rings, width, detector_half)[0]
              for detector_half in (0, 1)]
    stage.require(all(len(values) >= footprint.MINIMUM_SPLIT_PATCHES for values in halves),
                  "training band lost detector-half support")
    return np.vstack(halves), [len(values) for values in halves]


def precision_grid(model: dict[str, object], template: np.ndarray, mask: np.ndarray,
                   scale: np.ndarray, policies: tuple[tuple[str, str, float], ...]) \
        -> tuple[dict[str, np.ndarray], dict[str, float], dict[str, object]]:
    """Solve full, clipped, or truncated precision on one masked covariance."""
    template = np.asarray(template, dtype=np.float64)
    scale = np.asarray(scale, dtype=np.float64)
    mask = np.asarray(mask, dtype=bool)
    transformed = template / scale
    covariance = extension.dense_covariance(model, SUPPORT)
    selected = np.flatnonzero(mask.ravel())
    submatrix = covariance[np.ix_(selected, selected)]
    eigenvalues, eigenvectors = np.linalg.eigh(submatrix)
    target = float(model["target_variance"])
    stage.require(np.all(eigenvalues > 0) and np.isclose(np.mean(np.diag(submatrix)), target,
                  rtol=2e-12, atol=0), "PSD covariance eigensystem is invalid")
    coefficients = eigenvectors.T @ transformed.ravel()[selected]
    full_precision = 1 / eigenvalues
    full_energy = float(np.sum(full_precision * np.square(coefficients)))
    full_variance = 1 / full_energy
    weights: dict[str, np.ndarray] = {}
    sigmas: dict[str, float] = {}
    records: dict[str, object] = {}
    for name, kind, cutoff in policies:
        threshold = cutoff * target
        if kind == "full":
            values = full_precision
        elif kind == "clip":
            values = 1 / np.maximum(eigenvalues, threshold)
        else:
            values = np.where(eigenvalues >= threshold, 1 / eigenvalues, 0)
        retained = values > 0
        stage.require(np.any(retained), "precision policy removed every covariance mode")
        inverse_template = eigenvectors @ (values * coefficients)
        energy = float(transformed.ravel()[selected] @ inverse_template)
        stage.require(np.isfinite(energy) and energy > 0, "precision policy has no template support")
        transformed_weight = np.zeros(template.size, dtype=np.float64)
        transformed_weight[selected] = inverse_template / energy
        physical_weight = transformed_weight / scale.ravel()
        stage.require(np.isclose(physical_weight @ template.ravel(), 1, rtol=3e-10, atol=3e-12),
                      "precision weight lost unit physical response")
        modeled_variance = float(transformed_weight[selected] @ submatrix @ transformed_weight[selected])
        weights[name] = physical_weight
        sigmas[name] = 1 / math.sqrt(energy)
        records[name] = {
            "conditional_sigma": sigmas[name],
            "psd_sigma": math.sqrt(modeled_variance),
            "efficiency_if_psd_exact": math.sqrt(full_variance / modeled_variance),
            "retained_modes": int(np.count_nonzero(retained)),
            "modified_modes": int(np.count_nonzero(eigenvalues < threshold)) if cutoff else 0,
            "condition_number": float(eigenvalues[-1] / eigenvalues[0]),
            "minimum_over_mean": float(eigenvalues[0] / target),
            "maximum_over_mean": float(eigenvalues[-1] / target),
            "samples": None,
        }
    return weights, sigmas, {"covariance": covariance, "eigenvalues": eigenvalues,
                              "selected": selected, "policies": records}


def fit_query(image: np.ndarray, query: tuple[int, int], searches: list[tuple[int, int]],
              template: np.ndarray, validity: np.ndarray, planet: tuple[float, float],
              width: int) -> dict[str, object]:
    """Fit the two frozen PSD families and complete radial precision grid."""
    finite = np.isfinite(image)
    rings = raw.training_stencils(finite, query, searches, planet, SUPPORT)
    raw_samples, raw_counts = combined_samples(image, rings, width)
    rectangular = raw.regularize(raw.fit_periodogram_base(raw_samples, SUPPORT, "rectangular"), 0.3)

    excluded = radial.exclusion_mask(image.shape, searches, SUPPORT, planet)
    counts = radial.profile_counts(image, excluded)
    stage.require(np.all(counts >= radial.MINIMUM_PROFILE_PIXELS),
                  "strict radial profile lost its frozen support")
    scale_map, profile = radial.variance_profile(image, excluded)
    standardized = image / scale_map
    radial_samples, radial_counts = combined_samples(standardized, rings, width)
    hann = raw.regularize(raw.fit_periodogram_base(radial_samples, SUPPORT, "hann"), 0.1)

    support = np.asarray(validity, dtype=bool) & candidate_mask(query, planet)
    data = stamp(image, query)
    scale = stamp(scale_map, query)
    stage.require(support[HALF, HALF] and np.all(np.isfinite(data[support])) and
                  np.all(np.isfinite(scale[support])) and np.all(scale[support] > 0),
                  "candidate support is invalid")
    raw_weights, raw_sigmas, raw_detail = precision_grid(
        rectangular, template, support, np.ones(template.shape),
        (("raw_rectangular_m0p3", "full", 0.0),))
    radial_weights, radial_sigmas, radial_detail = precision_grid(
        hann, template, support, scale, RADIAL_POLICIES)
    for record in raw_detail["policies"].values():
        record["samples"] = len(raw_samples)
        record["split_samples"] = raw_counts
    for record in radial_detail["policies"].values():
        record["samples"] = len(radial_samples)
        record["split_samples"] = radial_counts
    return {
        "weights": raw_weights | radial_weights,
        "sigmas": raw_sigmas | radial_sigmas,
        "isotropic_mean": np.asarray(rectangular["mean"], dtype=np.float64),
        "isotropic_sigma": math.sqrt(float(rectangular["target_variance"]) /
                                      float(np.sum(np.square(template[support])))),
        "support": support,
        "scale": scale,
        "raw_model": rectangular,
        "radial_model": hann,
        "raw_detail": raw_detail,
        "radial_detail": radial_detail,
        "profile_minimum_pixels": int(np.min(counts)),
        "profile": profile,
    }


def sparse_products(manifest: Path, mode_index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load one mode of the frozen sparse radial response model."""
    prefix = manifest.name.removesuffix("manifest.fits")
    response_path = manifest.parent / (prefix + f"mode{mode_index:03d}_radial_response.fits")
    validity_path = manifest.parent / (prefix + f"mode{mode_index:03d}_radial_validity.fits")
    response, header = fits.getdata(response_path, header=True, memmap=True)
    validity = np.asarray(fits.getdata(validity_path, memmap=True), dtype=np.float64) > 0.5
    radii = np.asarray([float(value) for value in str(header["KLIP PSF SAMPLE RADII"]).split(",")])
    return np.asarray(response, dtype=np.float64), validity, radii


def reference_weights(template: np.ndarray, validity: np.ndarray, query: tuple[int, int],
                      sparse: tuple[np.ndarray, np.ndarray], planet: tuple[float, float]) \
        -> tuple[np.ndarray, np.ndarray]:
    """Build exact, sparse, and response-smoothed identity weights."""
    mask = np.asarray(validity, dtype=bool) & candidate_mask(query, planet)
    sparse_template, sparse_validity = sparse
    sparse_mask = np.asarray(sparse_validity, dtype=bool) & candidate_mask(query, planet)
    stage.require(mask[HALF, HALF] and sparse_mask[HALF, HALF], "reference source anchor is invalid")
    templates = [template, sparse_template]
    for fwhm in (1.8, 2.7):
        sigma = fwhm / math.sqrt(8 * math.log(2))
        templates.append(gaussian_filter(template, sigma=sigma, mode="constant", cval=0.0,
                                         truncate=4.0))
    masks = [mask, sparse_mask, mask, mask]
    weights = np.stack([normalized_weight(value, support)
                        for value, support in zip(templates, masks)])
    responses = np.asarray([weight @ template.ravel() for weight in weights], dtype=np.float64)
    return weights, responses


def unit_name(radius: float, mode: int) -> str:
    """Return one calibration-unit directory name."""
    return f"r{radius_tag(radius)}_mode{mode}"


def load_unit(root: Path, radius: float, mode: int) -> dict[str, np.ndarray]:
    """Load and verify one completed calibration unit."""
    directory = root / "calibration" / "units" / unit_name(radius, mode)
    stage.verify(read(directory / "complete.json")["products"])
    with np.load(directory / "models.npz", allow_pickle=False) as source:
        return {name: np.array(source[name], copy=True) for name in source.files}


def calculate_unit(root_value: str, radius_value: float, mode_index: int) -> str:
    """Calculate baseline maps and frozen site weights for one radius and mode."""
    root = Path(root_value)
    protocol = read(root / "protocol.json")
    radius_value = float(radius_value)
    mode = stage.MODES[mode_index]
    directory = root / "calibration" / "units" / unit_name(radius_value, mode)
    directory.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()

    baseline_cube, header = fits.getdata(protocol["paths"]["baseline"], header=True, memmap=True)
    image = np.asarray(baseline_cube[mode_index], dtype=np.float64)
    parent = Path(str(protocol["parent_response"]))
    paths = response47.product_paths(parent)
    coordinates = np.asarray(fits.getdata(paths["coordinates"], memmap=True), dtype=np.float64).T
    responses = np.asarray(fits.getdata(paths["responses"][mode_index], memmap=True), dtype=np.float64)
    validities = np.asarray(fits.getdata(paths["validities"][mode_index], memmap=True), dtype=np.float64) > 0.5
    sparse_cube, sparse_validity, sparse_radii = sparse_products(
        Path(str(protocol["paths"]["sparse_manifest"])), mode_index)
    center_row = 0.5 * (image.shape[1] - 1)
    center_column = 0.5 * (image.shape[0] - 1)
    planet = raw.planet_position(image.shape, {"known_planet": protocol["known_planet"]})
    radius_map = image_radius(image.shape)
    sites = [site for site in protocol["sites"]
             if float(site["nominal_radius"]) == radius_value and site["role"] in ANALYSIS_ROLES]
    bins, selected_pixels = required_bins(sites, radius_map)
    selected_sources = []
    for source, (row_value, column_value, _, _) in enumerate(coordinates):
        row, column = int(row_value), int(column_value)
        central_validity = raw.crop(np.asarray(validities[source], dtype=bool).T, SUPPORT)
        query = (row, column)
        if (selected_pixels[column, row] and np.all(central_validity) and
                candidate_mask(query, planet)[HALF, HALF]):
            selected_sources.append(source)
    stage.require(selected_sources, "calibration unit has no exact-response positions")

    count = len(selected_sources)
    positions = np.empty((count, 2), dtype=np.int16)
    source_indices = np.asarray(selected_sources, dtype=np.int32)
    reference = np.full((count, len(REFERENCE_WEIGHT_METHODS), SUPPORT * SUPPORT), np.nan)
    covariance = np.full((count, len(COVARIANCE_METHODS), SUPPORT * SUPPORT), np.nan)
    isotropic_means = np.full((count, SUPPORT * SUPPORT), np.nan)
    conditional_sigmas = np.full((count, 1 + len(COVARIANCE_METHODS)), np.nan)
    source_responses = np.full((count, len(METHODS)), np.nan)
    maps = np.full((len(METHODS), *image.shape), np.nan, dtype=np.float32)
    maps[METHODS.index("native")][selected_pixels] = image[selected_pixels]
    for name, fwhm in GAUSSIAN_METHODS:
        smoothed = preparation.gaussian_map(image, fwhm)
        maps[METHODS.index(name)][selected_pixels] = smoothed[selected_pixels]

    generic_diagnostics = []
    lookup: dict[tuple[int, int], int] = {}
    width = int(protocol["training"]["candidate_specific_half_width_by_radius"][str(radius_value)])
    for local, source in enumerate(selected_sources):
        row, column = map(int, coordinates[source, :2])
        query = (row, column)
        positions[local] = query
        lookup[query] = local
        template = raw.crop(np.asarray(responses[source], dtype=np.float64).T, SUPPORT)
        validity = raw.crop(np.asarray(validities[source], dtype=bool).T, SUPPORT)
        actual_radius = math.hypot(row - center_row, column - center_column)
        angle = math.atan2(row - center_row, column - center_column)
        sparse_template, sparse_support = stage.evaluate_sparse_response(
            sparse_cube, sparse_validity, sparse_radii, actual_radius, angle)
        weights, response_scales = reference_weights(template, validity, query,
                                                      (sparse_template, sparse_support), planet)
        reference[local] = weights
        data = stamp(image, query).ravel()
        for method_index, name in enumerate(REFERENCE_WEIGHT_METHODS):
            maps[METHODS.index(name), column, row] = weights[method_index] @ data
            source_responses[local, METHODS.index(name)] = response_scales[method_index]
        source_responses[local, METHODS.index("native")] = template[HALF, HALF]
        for name, fwhm in GAUSSIAN_METHODS:
            sigma = fwhm / math.sqrt(8 * math.log(2))
            filtered = gaussian_filter(template, sigma=sigma, mode="constant", cval=0.0, truncate=4.0)
            source_responses[local, METHODS.index(name)] = filtered[HALF, HALF]
        searches = [(row + delta_row, column + delta_column)
                    for delta_row, delta_column in footprint.SEARCH_OFFSETS]
        try:
            fitted = fit_query(image, query, searches, template, validity, planet, width)
        except RuntimeError as error:
            generic_diagnostics.append({"row": row, "column": column, "valid": False,
                                        "reason": str(error)})
            continue
        isotropic_means[local] = fitted["isotropic_mean"]
        identity_weight = weights[REFERENCE_WEIGHT_METHODS.index("exact_identity")]
        maps[METHODS.index("exact_isotropic_fitted_mean"), column, row] = (
            identity_weight @ (data - fitted["isotropic_mean"]))
        source_responses[local, METHODS.index("exact_isotropic_fitted_mean")] = (
            identity_weight @ template.ravel())
        conditional_sigmas[local, 0] = fitted["isotropic_sigma"]
        for covariance_index, name in enumerate(COVARIANCE_METHODS):
            weight = fitted["weights"][name]
            covariance[local, covariance_index] = weight
            maps[METHODS.index(name), column, row] = weight @ data
            source_responses[local, METHODS.index(name)] = weight @ template.ravel()
            conditional_sigmas[local, covariance_index + 1] = fitted["sigmas"][name]
        generic_diagnostics.append({
            "row": row, "column": column, "valid": True,
            "profile_minimum_pixels": fitted["profile_minimum_pixels"],
            "raw_samples": fitted["raw_detail"]["policies"]["raw_rectangular_m0p3"]["samples"],
            "radial_samples": fitted["radial_detail"]["policies"]["radial_hann_m0p1_full"]["samples"],
        })

    site_names = np.asarray([str(site["name"]) for site in sites])
    site_positions = np.empty((len(sites), len(footprint.SEARCH_OFFSETS), 2), dtype=np.int16)
    site_covariance = np.empty((len(sites), len(footprint.SEARCH_OFFSETS),
                                len(COVARIANCE_METHODS), SUPPORT * SUPPORT), dtype=np.float64)
    site_means = np.empty((len(sites), len(footprint.SEARCH_OFFSETS), SUPPORT * SUPPORT), dtype=np.float64)
    site_sigmas = np.empty((len(sites), len(footprint.SEARCH_OFFSETS),
                            1 + len(COVARIANCE_METHODS)), dtype=np.float64)
    site_baseline = np.empty((len(sites), len(footprint.SEARCH_OFFSETS), len(METHODS)), dtype=np.float64)
    site_diagnostics = {}
    for site_index, site in enumerate(sites):
        searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                    for delta_row, delta_column in footprint.SEARCH_OFFSETS]
        details = []
        for search_index, query in enumerate(searches):
            stage.require(query in lookup, f"site query is outside required annular maps: {site['name']}")
            local = lookup[query]
            source = int(source_indices[local])
            template = raw.crop(np.asarray(responses[source], dtype=np.float64).T, SUPPORT)
            validity = raw.crop(np.asarray(validities[source], dtype=bool).T, SUPPORT)
            fitted = fit_query(image, query, searches, template, validity, planet, width)
            data = stamp(image, query).ravel()
            site_positions[site_index, search_index] = query
            site_means[site_index, search_index] = fitted["isotropic_mean"]
            site_sigmas[site_index, search_index, 0] = fitted["isotropic_sigma"]
            site_baseline[site_index, search_index] = maps[:, query[1], query[0]]
            identity_weight = reference[local, REFERENCE_WEIGHT_METHODS.index("exact_identity")]
            site_baseline[site_index, search_index, METHODS.index("exact_isotropic_fitted_mean")] = (
                identity_weight @ (data - fitted["isotropic_mean"]))
            for covariance_index, name in enumerate(COVARIANCE_METHODS):
                weight = fitted["weights"][name]
                site_covariance[site_index, search_index, covariance_index] = weight
                site_sigmas[site_index, search_index, covariance_index + 1] = fitted["sigmas"][name]
                site_baseline[site_index, search_index, METHODS.index(name)] = weight @ data
            details.append({
                "query_row_column": list(query),
                "profile_minimum_pixels": fitted["profile_minimum_pixels"],
                "raw": fitted["raw_detail"]["policies"],
                "radial": fitted["radial_detail"]["policies"],
            })
        site_diagnostics[str(site["name"])] = details

    np.savez_compressed(directory / "models.npz", positions=positions, source_indices=source_indices,
                        reference_weights=reference, covariance_weights=covariance,
                        isotropic_means=isotropic_means, conditional_sigmas=conditional_sigmas,
                        source_responses=source_responses, site_names=site_names,
                        site_positions=site_positions, site_covariance_weights=site_covariance,
                        site_isotropic_means=site_means, site_conditional_sigmas=site_sigmas,
                        site_baseline_amplitudes=site_baseline)
    output_header = header.copy()
    output_header["HCI FILTER LABELS"] = ",".join(METHODS)
    output_header["HCI RADIAL BINS"] = ",".join(map(str, bins))
    fits.writeto(directory / "baseline_amplitudes.fits", maps, output_header)
    valid_counts = {name: int(np.count_nonzero(np.isfinite(maps[index])))
                    for index, name in enumerate(METHODS)}
    result = {"radius": radius_value, "mode": mode, "training_half_width": width,
              "required_bins": bins, "generic_positions": count,
              "valid_amplitude_pixels": valid_counts, "sites": list(site_names),
              "generic_diagnostics": generic_diagnostics,
              "site_diagnostics": site_diagnostics,
              "elapsed_seconds": time.monotonic() - started}
    stage.write_json(directory / "results.json", result)
    products = [stage.fingerprint(directory / name)
                for name in ("models.npz", "baseline_amplitudes.fits", "results.json")]
    stage.write_json(directory / "complete.json", {"status": "complete", "products": products})
    return str(directory / "complete.json")


def annular_oracle(amplitudes: np.ndarray, exclusions: list[tuple[float, float, float]],
                    supported_edge: bool = False) -> np.ndarray:
    """Reconstruct production annular SNR, optionally extending the nearest supported edge bin."""
    column, row = np.indices(amplitudes.shape)
    center_row = 0.5 * (amplitudes.shape[1] - 1)
    center_column = 0.5 * (amplitudes.shape[0] - 1)
    radii = np.hypot(row - center_row, column - center_column).astype(np.float32)
    noise = np.isfinite(amplitudes)
    for source_row, source_column, source_radius in exclusions:
        distance = np.hypot(row - source_row, column - source_column).astype(np.float32)
        noise &= distance > source_radius + 0.5
    maximum_bin = int(math.ceil(float(np.max(radii)))) + 1
    centers, means, deviations = [], [], []
    for lower in range(maximum_bin):
        values = amplitudes[noise & (radii > lower) & (radii <= lower + 1)].astype(np.float64)
        centers.append(lower + 0.5)
        means.append(float(np.mean(values)) if len(values) else np.nan)
        deviations.append(float(np.std(values, ddof=1)) if len(values) > 1 else np.nan)
    if supported_edge:
        valid = np.isfinite(means) & np.isfinite(deviations)
        stage.require(np.count_nonzero(valid) >= 2, "annular profile lacks supported interpolation bins")
        mean = np.interp(radii, np.asarray(centers)[valid], np.asarray(means)[valid])
        deviation = np.interp(radii, np.asarray(centers)[valid],
                              np.asarray(deviations)[valid]).astype(np.float32)
    else:
        mean = np.interp(radii, centers, means)
        deviation = np.interp(radii, centers, deviations).astype(np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        score = ((amplitudes - mean) / deviation).astype(np.float32)
        count = 2 * np.pi * radii / stage.LAMBDA_D - 1
        correction = np.where(count > 0, 1 / np.sqrt(1 + 1 / count), 0).astype(np.float32)
        score *= correction
    score[(radii < 0) | (radii > 60) | ~np.isfinite(amplitudes)] = np.nan
    return score


def production_snr(root: Path, analyzer: Path, protocol: dict[str, object], maps: np.ndarray,
                   header: fits.Header, site: dict[str, object], directory: Path) \
        -> tuple[np.ndarray, float, list[str]]:
    """Run hciAnalyze on arbitrary amplitude maps and verify an independent oracle."""
    directory.mkdir(parents=True, exist_ok=True)
    planet = protocol["known_planet"]
    command = [str(analyzer), "--file=amplitudes.fits", f"--lambdaD={stage.LAMBDA_D}",
               f"--planet.sep={planet['separation']},{site['separation']}",
               f"--planet.PA={planet['position_angle']},{site['position_angle']}",
               f"--planet.R={planet['exclusion_radius']},{planet['exclusion_radius']}",
               "--snr.apertureR=60", "--snr.minRad=0", "--snr.maxRad=60",
               "--filter.psfResponse=", "--filter.lpfGaussFW=0", "--filter.hpfGaussFW=0",
               "--noise.model=identity", "--noise.only=false", "--noise.outputDiagnostics=false"]
    flat = maps.reshape((-1, *maps.shape[-2:])).astype(np.float32)
    working_header = header.copy()
    working_header["HCI FILTER LABELS"] = ",".join(
        f"m{mode}_{method}" for mode in stage.MODES for method in METHODS)
    with tempfile.TemporaryDirectory(dir=directory) as temporary_value:
        temporary = Path(temporary_value)
        fits.writeto(temporary / "amplitudes.fits", flat, working_header)
        with (directory / "analysis.log").open("w", encoding="utf-8") as log:
            subprocess.run(command, cwd=temporary, env=environment(protocol, False), stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        snr, snr_header = fits.getdata(temporary / "amplitudes_snr.fits", header=True)
    snr = np.asarray(snr, dtype=np.float64).reshape(maps.shape)
    stage.require(int(snr_header["SNRMEAN"]) == 1 and int(snr_header["SNRSMALL"]) == 1,
                  "hciAnalyze annular SNR contract changed")
    planet_position = raw.planet_position(maps.shape[-2:], {"known_planet": planet})
    exclusions = [(planet_position[0], planet_position[1], float(planet["exclusion_radius"])),
                  (float(site["row"]), float(site["column"]), float(planet["exclusion_radius"]))]
    checked = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
               for delta_row, delta_column in footprint.SEARCH_OFFSETS]
    radii = image_radius(maps.shape[-2:])
    maximum_error = 0.0
    boundary_fallbacks = []
    for mode_index, mode in enumerate(stage.MODES):
        for method_index, method in enumerate(METHODS):
            expected = annular_oracle(maps[mode_index, method_index], exclusions)
            supported = annular_oracle(maps[mode_index, method_index], exclusions, supported_edge=True)
            for row, column in checked:
                observed_value = float(snr[mode_index, method_index, column, row])
                expected_value = float(expected[column, row])
                if np.isfinite(expected_value):
                    stage.require(np.isfinite(observed_value) and np.isclose(
                                  observed_value, expected_value, rtol=2e-6, atol=2e-6),
                                  f"annular oracle mismatch for {method}")
                    maximum_error = max(maximum_error, abs(observed_value - expected_value))
                    continue
                radius_value = float(radii[column, row])
                replacement = float(supported[column, row])
                stage.require(radius_value < stage.SNR_MIN_RADIUS + 0.5 and
                              observed_value == 0 and np.isfinite(replacement),
                              f"unsupported non-boundary annular SNR for {method}")
                snr[mode_index, method_index, column, row] = replacement
                boundary_fallbacks.append({
                    "mode": mode, "method": method, "row": row, "column": column,
                    "radius": radius_value, "production_value": observed_value,
                    "replacement": replacement,
                    "policy": "nearest supported one-pixel annular mean and sample deviation; exact small-sample correction at the candidate radius"})
    stage.write_json(directory / "command.json", command)
    stage.write_json(directory / "annular_verification.json", {
        "maximum_production_oracle_error": maximum_error,
        "boundary_fallbacks": boundary_fallbacks})
    return snr, maximum_error, command, boundary_fallbacks


def calibration_site(root: Path, protocol: dict[str, object], analyzer: Path,
                     radius_value: float, site: dict[str, object], header: fits.Header) -> dict[str, object]:
    """Measure one maximum-null site from the already frozen baseline maps and weights."""
    directory = root / "calibration" / "sites" / str(site["name"])
    complete = directory / "complete.json"
    if complete.exists():
        stage.verify(read(complete)["products"])
        existing = read(directory / "scores.json")
        if "annular_boundary_fallbacks" in existing:
            return existing
        archive_incomplete(root, directory, "calibration_sites")
    if directory.exists():
        archive_incomplete(root, directory, "calibration_sites")
    directory.mkdir(parents=True)
    maps = []
    amplitudes = []
    sigmas = []
    for mode in stage.MODES:
        unit = load_unit(root, radius_value, mode)
        unit_directory = root / "calibration" / "units" / unit_name(radius_value, mode)
        maps.append(np.asarray(fits.getdata(unit_directory / "baseline_amplitudes.fits"), dtype=np.float64))
        site_index = list(unit["site_names"]).index(str(site["name"]))
        site_amplitude = unit["site_baseline_amplitudes"][site_index]
        site_sigma = unit["site_conditional_sigmas"][site_index]
        amplitudes.append(site_amplitude)
        sigmas.append(site_sigma)
        for search_index, (row, column) in enumerate(unit["site_positions"][site_index]):
            maps[-1][:, int(column), int(row)] = site_amplitude[search_index]
    maps_array = np.stack(maps)
    snr, oracle_error, command, boundary_fallbacks = production_snr(
        root, analyzer, protocol, maps_array, header, site, directory)
    records = {}
    searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                for delta_row, delta_column in footprint.SEARCH_OFFSETS]
    for mode_index, mode in enumerate(stage.MODES):
        methods = {}
        for method_index, name in enumerate(METHODS):
            values = [float(snr[mode_index, method_index, column, row]) for row, column in searches]
            raw_values = [float(amplitudes[mode_index][search_index, method_index])
                          for search_index in range(len(searches))]
            stage.require(np.all(np.isfinite(values)) and np.all(np.isfinite(raw_values)),
                          f"calibration site lacks {name} support")
            entry = {"snr_pixels": values, "search_score": max(values),
                     "amplitude_pixels": raw_values}
            if name == "exact_isotropic_fitted_mean":
                sigma = np.asarray(sigmas[mode_index])[:, 0]
                entry["conditional_score_pixels"] = (np.asarray(raw_values) / sigma).tolist()
            elif name in COVARIANCE_METHODS:
                sigma = np.asarray(sigmas[mode_index])[:, 1 + COVARIANCE_METHODS.index(name)]
                entry["conditional_score_pixels"] = (np.asarray(raw_values) / sigma).tolist()
            methods[name] = entry
        records[str(mode)] = methods
    result = {"site": str(site["name"]), "radius": radius_value,
              "annular_oracle_maximum_error": oracle_error,
              "annular_boundary_fallbacks": boundary_fallbacks,
              "methods_by_mode": records}
    stage.write_json(directory / "scores.json", result)
    products = [stage.fingerprint(directory / name)
                for name in ("analysis.log", "command.json", "annular_verification.json", "scores.json")]
    stage.write_json(complete, {"status": "complete", "products": products})
    return result


def thresholds(records: list[dict[str, object]]) -> dict[str, object]:
    """Freeze method-specific maxima of the 20 calibration-site searches."""
    result = {}
    for radius_value in preparation.RADII:
        radius_records = [record for record in records if record["radius"] == radius_value]
        stage.require(len(radius_records) == preparation.ROLE_COUNTS["calibration"],
                      "calibration radius lost fixed sites")
        modes = {}
        for mode in stage.MODES:
            methods = {}
            for name in METHODS:
                values = [float(record["methods_by_mode"][str(mode)][name]["search_score"])
                          for record in radius_records]
                ordered = sorted(values, reverse=True)
                methods[name] = {"threshold": ordered[0], "site_maxima": values,
                                 "largest": ordered[0], "second_largest": ordered[1],
                                 "largest_second_gap": ordered[0] - ordered[1],
                                 "strict_exceedance": True}
            modes[str(mode)] = methods
        result[str(radius_value)] = modes
    return result


def verify_calibration(root: Path) -> dict[str, object]:
    """Verify the calibration receipt and every nested unit and site product."""
    completion = read(root / "calibration" / "complete.json")
    stage.require(completion["status"] == "complete" and completion["thresholds_frozen"] and
                  not completion["positive_products_opened"], "calibration receipt is incomplete")
    stage.verify(completion["products"])
    for receipt_path in sorted((root / "calibration" / "units").glob("*/complete.json")):
        stage.verify(read(receipt_path)["products"])
    for receipt_path in sorted((root / "calibration" / "sites").glob("*/complete.json")):
        stage.verify(read(receipt_path)["products"])
    return completion


def calibrate(root: Path, workers: int) -> None:
    """Complete and receipt all null calibration before any positive reduction."""
    protocol, _, analyzer = enable(root)
    completion = root / "calibration" / "complete.json"
    if completion.exists():
        verify_calibration(root)
        print("Stage-C calibration already complete.", flush=True)
        return
    stage.require(not (root / "reductions").exists(),
                  "calibration must finish before positive reductions exist")
    calibration_root = root / "calibration"
    calibration_root.mkdir(exist_ok=True)
    (calibration_root / "units").mkdir(exist_ok=True)
    (calibration_root / "sites").mkdir(exist_ok=True)
    pending = []
    receipts = []
    for radius_value in preparation.RADII:
        for mode_index, mode in enumerate(stage.MODES):
            directory = calibration_root / "units" / unit_name(radius_value, mode)
            receipt = directory / "complete.json"
            if receipt.exists():
                stage.verify(read(receipt)["products"])
                receipts.append(receipt)
            else:
                if directory.exists():
                    archive_incomplete(root, directory, "calibration_units")
                pending.append((radius_value, mode_index))
    if pending:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(calculate_unit, str(root), radius_value, mode_index):
                       (radius_value, stage.MODES[mode_index])
                       for radius_value, mode_index in pending}
            for completed, future in enumerate(as_completed(futures), start=1):
                receipt = Path(future.result())
                receipts.append(receipt)
                radius_value, mode = futures[future]
                print(f"calibration unit {completed}/{len(pending)}: radius {radius_value:g}, mode {mode}",
                      flush=True)
    baseline, header = fits.getdata(protocol["paths"]["baseline"], header=True, memmap=True)
    stage.require(np.asarray(baseline).shape == (len(stage.MODES), 128, 128),
                  "Stage-C baseline shape changed")
    records = []
    sites = [site for site in protocol["sites"] if site["role"] == "calibration"]
    for index, site in enumerate(sites, start=1):
        records.append(calibration_site(root, protocol, analyzer,
                                        float(site["nominal_radius"]), site, header))
        print(f"calibration site {index}/{len(sites)}: {site['name']}", flush=True)
    frozen_thresholds = thresholds(records)
    stage.write_json(calibration_root / "records.json", records)
    stage.write_json(calibration_root / "thresholds.json", frozen_thresholds)
    products = [stage.fingerprint(root / "development_manifest.json"),
                stage.fingerprint(calibration_root / "records.json"),
                stage.fingerprint(calibration_root / "thresholds.json")]
    products.extend(stage.fingerprint(path) for path in sorted(receipts))
    products.extend(stage.fingerprint(calibration_root / "sites" / str(site["name"]) / "complete.json")
                    for site in sites)
    stage.write_json(completion, {"status": "complete", "positive_products_opened": False,
                                  "thresholds_frozen": True, "products": products})
    stage.write_json(root / "state.json", {"status": "calibrated", "development_reductions": 0,
                                           "validation_products_opened": False,
                                           "positive_analysis_started": False})
    print(completion, flush=True)


def validate_reduction(path: Path) -> dict[str, object]:
    """Validate one positive KLIP cube and return its fingerprint."""
    data, header = fits.getdata(path, header=True, memmap=True)
    stage.require(np.asarray(data).shape == (len(stage.MODES), 128, 128) and
                  stage.read_modes(header) == stage.MODES, "development reduction schema changed")
    return stage.fingerprint(path)


def reduce(root: Path) -> None:
    """Run or verify all 108 frozen development reductions after calibration."""
    protocol, _, _ = enable(root)
    verify_calibration(root)
    commands = read(root / "commands.json")
    stage.require(len(commands) == protocol["expected_development_reductions"],
                  "frozen development command count changed")
    reductions = root / "reductions"
    reductions.mkdir(exist_ok=True)
    for index, record in enumerate(commands, start=1):
        name = str(record["task"])
        directory = reductions / name
        complete = directory / "complete.json"
        if complete.exists():
            stage.verify(read(complete)["products"])
        else:
            if directory.exists():
                archive_incomplete(root, directory, "reductions")
            directory.mkdir()
            elapsed = stage.run_command([str(value) for value in record["command"]], directory,
                                        directory / "run.log", environment(protocol, True))
            product = validate_reduction(directory / "finim.fits")
            stage.write_json(complete, {"status": "complete", "elapsed_seconds": elapsed,
                                        "products": [product, stage.fingerprint(directory / "run.log"),
                                                     stage.fingerprint(directory / "run.command.json")]})
        print(f"development reduction {index}/{len(commands)}: {name}", flush=True)
    stage.write_json(root / "state.json", {"status": "reductions_complete",
                                           "development_reductions": len(commands),
                                           "validation_products_opened": False,
                                           "positive_analysis_started": False})


def positive_maps(image: np.ndarray, radius_value: float, mode: int, unit: dict[str, np.ndarray],
                  baseline_unit: Path) -> np.ndarray:
    """Apply all baseline-frozen generic filters to one positive image."""
    maps = np.full((len(METHODS), *image.shape), np.nan, dtype=np.float64)
    baseline_maps = np.asarray(fits.getdata(baseline_unit / "baseline_amplitudes.fits"), dtype=np.float64)
    selected = np.any(np.isfinite(baseline_maps), axis=0)
    maps[METHODS.index("native")][selected] = image[selected]
    for name, fwhm in GAUSSIAN_METHODS:
        smoothed = preparation.gaussian_map(image, fwhm)
        maps[METHODS.index(name)][selected] = smoothed[selected]
    positions = unit["positions"]
    patches = np.stack([stamp(image, tuple(map(int, query))).ravel() for query in positions])
    reference_values = np.einsum("nkp,np->nk", unit["reference_weights"], patches)
    for reference_index, name in enumerate(REFERENCE_WEIGHT_METHODS):
        for local, (row, column) in enumerate(positions):
            maps[METHODS.index(name), int(column), int(row)] = reference_values[local, reference_index]
    identity = unit["reference_weights"][:, REFERENCE_WEIGHT_METHODS.index("exact_identity")]
    isotropic = np.einsum("np,np->n", identity, patches - unit["isotropic_means"])
    covariance_values = np.einsum("nkp,np->nk", unit["covariance_weights"], patches)
    for local, (row, column) in enumerate(positions):
        maps[METHODS.index("exact_isotropic_fitted_mean"), int(column), int(row)] = isotropic[local]
        for covariance_index, name in enumerate(COVARIANCE_METHODS):
            maps[METHODS.index(name), int(column), int(row)] = covariance_values[local, covariance_index]
    return maps


def fidelity(delta: np.ndarray, template: np.ndarray, fitted: dict[str, object]) -> dict[str, object]:
    """Compare a finite positive difference with the exact response in raw and PSD metrics."""
    mask = np.asarray(fitted["support"], dtype=bool).ravel()
    raw_delta = np.asarray(delta, dtype=np.float64).ravel()[mask]
    raw_template = np.asarray(template, dtype=np.float64).ravel()[mask]

    def metric(data: np.ndarray, model_template: np.ndarray, covariance: np.ndarray) -> dict[str, float]:
        submatrix = covariance[np.ix_(mask, mask)]
        inverse_data = np.linalg.solve(submatrix, data)
        inverse_template = np.linalg.solve(submatrix, model_template)
        cross = float(model_template @ inverse_data)
        template_energy = float(model_template @ inverse_template)
        data_energy = float(data @ inverse_data)
        projection = cross / template_energy
        cosine = cross / math.sqrt(template_energy * data_energy)
        residual = data - projection * model_template
        residual_energy = float(residual @ np.linalg.solve(submatrix, residual))
        return {"cosine": cosine, "projection_scale": projection,
                "best_scaled_relative_residual": math.sqrt(residual_energy / data_energy)}

    identity = np.eye(len(raw_delta), dtype=np.float64)
    result = {"unweighted": metric(raw_delta, raw_template, identity)}
    raw_covariance = fitted["raw_detail"]["covariance"]
    result["raw_rectangular_m0p3"] = metric(raw_delta, raw_template, raw_covariance)
    scale = np.asarray(fitted["scale"], dtype=np.float64).ravel()[mask]
    radial_covariance = fitted["radial_detail"]["covariance"]
    result["radial_hann_m0p1"] = metric(raw_delta / scale, raw_template / scale, radial_covariance)
    return result


def analyze_task(root_value: str, task_name: str) -> str:
    """Analyze one completed positive cube with frozen and per-image-refit covariance."""
    root = Path(root_value)
    protocol = read(root / "protocol.json")
    analyzer = Path(str(read(root / "development_manifest.json")["hcianalyze"]["path"]))
    task = next(value for value in protocol["development_tasks"] if value["name"] == task_name)
    site = next(value for value in protocol["sites"] if value["name"] == task["site"])
    radius_value = float(site["nominal_radius"])
    directory = root / "development_analysis" / task_name
    complete = directory / "complete.json"
    if complete.exists():
        stage.verify(read(complete)["products"])
        return str(directory / "measurements.json")
    if directory.exists():
        archive_incomplete(root, directory, "development_analysis")
    directory.mkdir(parents=True)
    started = time.monotonic()
    source = root / "reductions" / task_name / "finim.fits"
    validate_reduction(source)
    positive_cube, header = fits.getdata(source, header=True, memmap=True)
    baseline_cube = fits.getdata(protocol["paths"]["baseline"], memmap=True)
    parent = Path(str(protocol["parent_response"]))
    paths = response47.product_paths(parent)
    coordinates = np.asarray(fits.getdata(paths["coordinates"], memmap=True), dtype=np.float64).T
    coordinate_lookup = {(int(row), int(column)): index
                         for index, (row, column, _, _) in enumerate(coordinates)}
    planet = raw.planet_position((128, 128), {"known_planet": protocol["known_planet"]})
    searches = [(int(site["row"]) + delta_row, int(site["column"]) + delta_column)
                for delta_row, delta_column in footprint.SEARCH_OFFSETS]
    width = int(site["training_half_width"])
    frozen_maps, refit_maps = [], []
    baseline_amplitudes, source_responses = [], []
    fit_records, fidelity_records = {}, {}
    for mode_index, mode in enumerate(stage.MODES):
        unit = load_unit(root, radius_value, mode)
        unit_directory = root / "calibration" / "units" / unit_name(radius_value, mode)
        positive = np.asarray(positive_cube[mode_index], dtype=np.float64)
        baseline = np.asarray(baseline_cube[mode_index], dtype=np.float64)
        frozen = positive_maps(positive, radius_value, mode, unit, unit_directory)
        refit = np.array(frozen, copy=True)
        site_index = list(unit["site_names"]).index(str(site["name"]))
        baseline_amplitudes.append(unit["site_baseline_amplitudes"][site_index])
        local_source_responses = []
        responses = np.asarray(fits.getdata(paths["responses"][mode_index], memmap=True), dtype=np.float64)
        validities = np.asarray(fits.getdata(paths["validities"][mode_index], memmap=True), dtype=np.float64) > 0.5
        mode_fits = []
        for search_index, query in enumerate(searches):
            source_index = coordinate_lookup[query]
            template = raw.crop(np.asarray(responses[source_index], dtype=np.float64).T, SUPPORT)
            validity = raw.crop(np.asarray(validities[source_index], dtype=bool).T, SUPPORT)
            data = stamp(positive, query).ravel()
            baseline_data = stamp(baseline, query).ravel()
            frozen_weights = unit["site_covariance_weights"][site_index, search_index]
            identity_local = np.where(np.all(unit["positions"] == np.asarray(query), axis=1))[0]
            stage.require(len(identity_local) == 1, "site query is missing from generic unit")
            identity_weight = unit["reference_weights"][identity_local[0],
                                                       REFERENCE_WEIGHT_METHODS.index("exact_identity")]
            frozen[METHODS.index("exact_isotropic_fitted_mean"), query[1], query[0]] = (
                identity_weight @ (data - unit["site_isotropic_means"][site_index, search_index]))
            for covariance_index, name in enumerate(COVARIANCE_METHODS):
                frozen[METHODS.index(name), query[1], query[0]] = frozen_weights[covariance_index] @ data
            fitted = fit_query(positive, query, searches, template, validity, planet, width)
            refit[METHODS.index("exact_isotropic_fitted_mean"), query[1], query[0]] = (
                identity_weight @ (data - fitted["isotropic_mean"]))
            for name in COVARIANCE_METHODS:
                refit[METHODS.index(name), query[1], query[0]] = fitted["weights"][name] @ data
            response_vector = unit["source_responses"][identity_local[0]]
            local_source_responses.append(response_vector)
            mode_fits.append({"query_row_column": list(query),
                              "profile_minimum_pixels": fitted["profile_minimum_pixels"],
                              "raw": fitted["raw_detail"]["policies"],
                              "radial": fitted["radial_detail"]["policies"]})
            if search_index == 0:
                baseline_fit = fit_query(baseline, query, searches, template, validity,
                                         planet, width)
                for covariance_index, name in enumerate(COVARIANCE_METHODS):
                    stage.require(np.allclose(baseline_fit["weights"][name],
                                  frozen_weights[covariance_index], rtol=3e-10, atol=3e-12),
                                  "recomputed baseline weight differs from frozen calibration")
                delta = (stamp(positive, query) - stamp(baseline, query)) / float(task["contrast"])
                fidelity_records[str(mode)] = fidelity(delta, template, baseline_fit)
        frozen_maps.append(frozen)
        refit_maps.append(refit)
        source_responses.append(np.stack(local_source_responses))
        fit_records[str(mode)] = mode_fits
    frozen_array = np.stack(frozen_maps)
    refit_array = np.stack(refit_maps)
    frozen_snr, frozen_error, _, frozen_fallbacks = production_snr(
        root, analyzer, protocol, frozen_array, header, site, directory / "frozen")
    refit_snr, refit_error, _, refit_fallbacks = production_snr(
        root, analyzer, protocol, refit_array, header, site, directory / "refit")
    baseline_amplitudes_array = np.stack(baseline_amplitudes)
    source_responses_array = np.stack(source_responses)
    arms = {}
    for arm, maps, snr in (("frozen", frozen_array, frozen_snr), ("refit", refit_array, refit_snr)):
        mode_records = {}
        for mode_index, mode in enumerate(stage.MODES):
            method_records = {}
            for method_index, name in enumerate(METHODS):
                if arm == "refit" and name not in COVARIANCE_METHODS:
                    continue
                values = [float(snr[mode_index, method_index, column, row])
                          for row, column in searches]
                amplitudes = [float(maps[mode_index, method_index, column, row])
                              for row, column in searches]
                stage.require(np.all(np.isfinite(values)) and np.all(np.isfinite(amplitudes)),
                              f"positive task lacks {arm} {name} support")
                best = int(np.argmax(values))
                response_scale = float(source_responses_array[mode_index, 0, method_index])
                stage.require(np.isfinite(response_scale) and response_scale != 0,
                              f"invalid source response for {name}")
                center_estimate = amplitudes[0] / response_scale
                baseline_center = baseline_amplitudes_array[mode_index, 0, method_index] / response_scale
                method_records[name] = {
                    "snr_pixels": values, "search_score": values[best],
                    "center_snr": values[0], "best_search_index": best,
                    "localization_error_pixels": math.hypot(*footprint.SEARCH_OFFSETS[best]),
                    "amplitude_pixels": amplitudes, "source_response_center": response_scale,
                    "center_contrast_estimate": center_estimate,
                    "center_relative_contrast_error": (center_estimate - float(task["contrast"])) /
                                                      float(task["contrast"]),
                    "positive_minus_baseline_throughput": (center_estimate - baseline_center) /
                                                          float(task["contrast"]),
                }
            mode_records[str(mode)] = method_records
        arms[arm] = mode_records
    result = {"task": task, "site": site, "source": stage.fingerprint(source),
              "annular_oracle_maximum_error": {"frozen": frozen_error, "refit": refit_error},
              "annular_boundary_fallbacks": {"frozen": frozen_fallbacks, "refit": refit_fallbacks},
              "arms": arms, "fit_diagnostics": fit_records,
              "response_fidelity": fidelity_records,
              "elapsed_seconds": time.monotonic() - started}
    stage.write_json(directory / "measurements.json", result)
    products = [stage.fingerprint(directory / "measurements.json"),
                stage.fingerprint(directory / "frozen" / "analysis.log"),
                stage.fingerprint(directory / "frozen" / "command.json"),
                stage.fingerprint(directory / "frozen" / "annular_verification.json"),
                stage.fingerprint(directory / "refit" / "analysis.log"),
                stage.fingerprint(directory / "refit" / "command.json"),
                stage.fingerprint(directory / "refit" / "annular_verification.json")]
    stage.write_json(complete, {"status": "complete", "products": products})
    return str(directory / "measurements.json")


def summarize(root: Path, protocol: dict[str, object], measurements: list[dict[str, object]]) -> None:
    """Write complete all-mode development tables without selecting a policy."""
    frozen_thresholds = read(root / "calibration" / "thresholds.json")
    rows = []
    for mode in stage.MODES:
        for radius_value in preparation.RADII:
            for target in preparation.TARGET_SNRS:
                chosen = [record for record in measurements
                          if float(record["site"]["nominal_radius"]) == radius_value and
                          float(record["task"]["target_source_snr"]) == target]
                stage.require(len(chosen) == preparation.ROLE_COUNTS["development"],
                              "development summary lost paired sites")
                for arm in ("frozen", "refit"):
                    names = METHODS if arm == "frozen" else COVARIANCE_METHODS
                    for name in names:
                        values = [record["arms"][arm][str(mode)][name] for record in chosen]
                        threshold = float(frozen_thresholds[str(radius_value)][str(mode)][name]["threshold"])
                        rows.append({
                            "mode": mode, "radius": radius_value, "target_source_snr": target,
                            "arm": arm, "method": name, "sites": len(values),
                            "threshold": threshold,
                            "recoveries": sum(float(value["search_score"]) > threshold for value in values),
                            "mean_search_snr": float(np.mean([value["search_score"] for value in values])),
                            "median_search_snr": float(np.median([value["search_score"] for value in values])),
                            "mean_center_snr": float(np.mean([value["center_snr"] for value in values])),
                            "mean_localization_error_pixels": float(np.mean([
                                value["localization_error_pixels"] for value in values])),
                            "median_relative_contrast_error": float(np.median([
                                value["center_relative_contrast_error"] for value in values])),
                            "relative_contrast_error_scatter": float(np.std([
                                value["center_relative_contrast_error"] for value in values], ddof=1)),
                            "mean_positive_minus_baseline_throughput": float(np.mean([
                                value["positive_minus_baseline_throughput"] for value in values])),
                        })
    result = {"purpose": "KLIP Stage-C calibrated development result",
              "selection_performed": False, "validation_products_opened": False,
              "measurements": len(measurements), "summary_rows": rows,
              "maximum_annular_oracle_error": max(
                  value for record in measurements for value in record["annular_oracle_maximum_error"].values())}
    stage.write_json(root / "development_results.json", result)
    fields = list(rows[0])
    with (root / "development_results.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    primary = [row for row in rows if row["mode"] == preparation.PRIMARY_MODE and
               row["arm"] == "frozen"]
    lines = ["# KLIP Stage-C development result", "",
             "All maximum-null thresholds were frozen from the 20 calibration sites before any positive "
             "image was opened. This report contains development data only; it does not select a policy or "
             "open validation products.", "", "## Mode 200 frozen-weight comparison", "",
             "| Radius | Target | Method | Recovered / 6 | Threshold | Mean max SNR | Median max SNR | Mean center SNR | Mean throughput |",
             "| ---: | ---: | :--- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in primary:
        lines.append(f"| {row['radius']:g} | {row['target_source_snr']:g} | {row['method']} | "
                     f"{row['recoveries']} | {row['threshold']:.4f} | {row['mean_search_snr']:.4f} | "
                     f"{row['median_search_snr']:.4f} | {row['mean_center_snr']:.4f} | "
                     f"{row['mean_positive_minus_baseline_throughput']:.4f} |")
    lines.extend(["", "## Covariance refit minus frozen", "",
                  "The refit arm changes the fitted covariance at the five source-search pixels. Its annular "
                  "normalization uses the same baseline-frozen surrounding map; the complete trial disk is "
                  "excluded from that profile.", "",
                  "| Radius | Target | Method | Recovery change | Mean max-SNR change |",
                  "| ---: | ---: | :--- | ---: | ---: |"])
    frozen_lookup = {(row["radius"], row["target_source_snr"], row["method"]): row
                     for row in primary}
    for row in rows:
        if row["mode"] != preparation.PRIMARY_MODE or row["arm"] != "refit":
            continue
        parent = frozen_lookup[(row["radius"], row["target_source_snr"], row["method"])]
        lines.append(f"| {row['radius']:g} | {row['target_source_snr']:g} | {row['method']} | "
                     f"{row['recoveries'] - parent['recoveries']:+d} | "
                     f"{row['mean_search_snr'] - parent['mean_search_snr']:+.4f} |")
    (root / "development_results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    plotted = ("gaussian_fwhm3p6", "exact_identity", "raw_rectangular_m0p3",
               "radial_hann_m0p1_full", "radial_hann_m0p1_clip1p0",
               "radial_hann_m0p1_trunc0p5")
    pooled = []
    for target in preparation.TARGET_SNRS:
        for name in plotted:
            selected = [row for row in primary if row["target_source_snr"] == target and
                        row["method"] == name]
            pooled.append((target, name, sum(row["recoveries"] for row in selected),
                           float(np.mean([row["mean_search_snr"] for row in selected]))))
    figure, axes = plt.subplots(1, 2, figsize=(15, 5), layout="constrained")
    for name in plotted:
        chosen = [row for row in pooled if row[1] == name]
        axes[0].plot([row[0] for row in chosen], [row[2] for row in chosen], marker="o", label=name)
        axes[1].plot([row[0] for row in chosen], [row[3] for row in chosen], marker="o", label=name)
    axes[0].set(xlabel="Target Gaussian source-only SNR", ylabel="Recoveries across 36 sites",
                title="Mode 200 frozen-threshold recovery")
    axes[1].set(xlabel="Target Gaussian source-only SNR", ylabel="Mean five-pixel maximum SNR",
                title="Mode 200 mean detection SNR")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    figure.savefig(root / "development_mode200.png", dpi=170)
    plt.close(figure)


def analyze(root: Path, workers: int) -> None:
    """Analyze all development positives after reduction and calibration completion."""
    protocol, _, _ = enable(root)
    verify_calibration(root)
    for task in protocol["development_tasks"]:
        stage.verify(read(root / "reductions" / str(task["name"]) / "complete.json")["products"])
    analysis_root = root / "development_analysis"
    analysis_root.mkdir(exist_ok=True)
    state = read(root / "state.json")
    stage.require(not state["validation_products_opened"], "validation products were opened prematurely")
    stage.write_json(root / "state.json", {"status": "analyzing",
                                           "development_reductions": len(protocol["development_tasks"]),
                                           "validation_products_opened": False,
                                           "positive_analysis_started": True})
    tasks = [str(task["name"]) for task in protocol["development_tasks"]]
    paths = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(analyze_task, str(root), name): name for name in tasks}
        for completed, future in enumerate(as_completed(futures), start=1):
            paths.append(Path(future.result()))
            print(f"development analysis {completed}/{len(tasks)}: {futures[future]}", flush=True)
    measurements = [read(path) for path in sorted(paths)]
    summarize(root, protocol, measurements)
    products = [stage.fingerprint(root / name) for name in
                ("development_results.json", "development_results.csv",
                 "development_results.md", "development_mode200.png")]
    products.extend(stage.fingerprint(root / "development_analysis" /
                    str(task["name"]) / "complete.json") for task in protocol["development_tasks"])
    stage.write_json(root / "development_complete.json", {"status": "complete",
                                                           "validation_products_opened": False,
                                                           "products": products})
    stage.write_json(root / "state.json", {"status": "development_complete",
                                           "development_reductions": len(protocol["development_tasks"]),
                                           "validation_products_opened": False,
                                           "positive_analysis_started": True})
    print(root / "development_results.md", flush=True)


def run(root: Path, workers: int) -> None:
    """Run calibration, development reductions, and development analysis in order."""
    calibrate(root, workers)
    reduce(root)
    analyze(root, workers)


def check() -> None:
    """Verify precision, identity normalization, annular algebra, and method mappings."""
    generator = np.random.default_rng(90625)
    samples = generator.normal(size=(48, SUPPORT * SUPPORT))
    template = generator.normal(size=(SUPPORT, SUPPORT))
    model = extension.fit_extended_psd(samples, SUPPORT, "hann", 0.1)
    mask = np.ones(template.shape, dtype=bool)
    mask[0, :3] = False
    scale = np.exp(generator.normal(scale=0.1, size=template.shape))
    weights, _, detail = precision_grid(model, template, mask, scale, RADIAL_POLICIES)
    stage.require(set(weights) == {name for name, _, _ in RADIAL_POLICIES} and
                  all(np.isclose(weight @ template.ravel(), 1, rtol=3e-10, atol=3e-12)
                      for weight in weights.values()), "precision-grid check failed")
    identity = normalized_weight(template, mask)
    stage.require(np.isclose(identity @ template.ravel(), 1), "identity check failed")
    image = generator.normal(size=(128, 128))
    image[:5] = np.nan
    score = annular_oracle(image, [(70.0, 70.0, 7.0)])
    stage.require(np.any(np.isfinite(score)), "annular oracle check failed")
    radii = image_radius(image.shape)
    boundary_image = generator.normal(size=image.shape)
    boundary_image[(radii > stage.SNR_MIN_RADIUS - 1) &
                   (radii <= stage.SNR_MIN_RADIUS)] = np.nan
    unsupported = annular_oracle(boundary_image, [])
    supported = annular_oracle(boundary_image, [], supported_edge=True)
    stage.require(not np.isfinite(unsupported[59, 59]) and np.isfinite(supported[59, 59]),
                  "annular supported-edge check failed")
    stage.require(tuple(METHODS) == tuple(preparation.METHODS) and
                  set(COVARIANCE_METHODS).issubset(METHODS) and
                  set(detail["policies"]) == {name for name, _, _ in RADIAL_POLICIES},
                  "Stage-C method mapping changed")
    print("KLIP Stage-C development runner checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the Stage-C development command-line parser."""
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    for action in ("calibrate", "reduce", "analyze", "run"):
        current = subparsers.add_parser(action)
        current.add_argument("root", type=Path)
        if action in {"calibrate", "analyze", "run"}:
            current.add_argument("--workers", type=int, default=4)
    return result


def main() -> None:
    """Dispatch the requested Stage-C development phase."""
    arguments = parser().parse_args()
    if arguments.action == "check":
        check()
        return
    root = arguments.root.resolve()
    ensure_frozen_runner(root, sys.argv[1:])
    workers = int(getattr(arguments, "workers", 1))
    stage.require(workers >= 1, "worker count must be positive")
    if arguments.action == "calibrate":
        calibrate(root, workers)
    elif arguments.action == "reduce":
        reduce(root)
    elif arguments.action == "analyze":
        analyze(root, workers)
    else:
        run(root, workers)


if __name__ == "__main__":
    main()
