#!/usr/bin/env python3
"""Run the frozen Stage-A audit for the KLIP covariance matched-filter study.

The preparation phase fingerprints every reduction input, archived response
product, executable, and analysis input before any new score is calculated.
The run phase makes one current signal-free baseline, verifies compatibility
with the archived exact-response baseline, audits response support and energy,
and independently reconstructs the exact identity-filter amplitude cube and
its annular SNR map.
"""
from __future__ import annotations

import argparse
import configparser
import csv
import hashlib
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


MODES = [125, 150, 175, 200, 225, 250, 300, 350]
PRIMARY_RADII = [7.5, 10.0, 12.0, 16.0, 20.0, 24.0]
PLANET_SEPARATION = 12.387750470790344
PLANET_PA = 260.64315155951886
PLANET_CONTRAST = 0.0045743624964148452
PLANET_EXCLUSION_RADIUS = 7.0
LAMBDA_D = 3.6
SNR_MIN_RADIUS = 6.0
SNR_MAX_RADIUS = 60.0
ARCHIVED_APERTURE_RADIUS = 2.0
ARCHIVED_PLANET_RADIUS = 5.0
NEIGHBORS = [(0, 0), (-1, 0), (1, 0), (0, -1), (0, 1)]


def require(condition: bool, message: str) -> None:
    """Raise a runtime error when an experimental contract is not satisfied."""
    if not condition:
        raise RuntimeError(message)


def finite_json(value: object) -> object:
    """Convert NumPy values and replace nonfinite floats for strict JSON."""
    if isinstance(value, dict):
        return {str(key): finite_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return finite_json(value.tolist())
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def write_json(path: Path, value: object) -> None:
    """Atomically write strict, indented JSON."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(finite_json(value), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def fingerprint(path: Path) -> dict[str, object]:
    """Return a resolved path, byte count, and SHA-256 digest for one file."""
    path = path.resolve()
    require(path.is_file(), f"required file is missing: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def verify(records: list[dict[str, object]]) -> None:
    """Verify that every frozen file still has its recorded size and digest."""
    for expected in records:
        observed = fingerprint(Path(str(expected["path"])))
        require(observed == expected, f"frozen file changed: {expected['path']}")


def resolve_executable(value: str) -> Path:
    """Resolve an executable spelling to an existing absolute path."""
    match = shutil.which(value)
    require(match is not None, f"executable was not found: {value}")
    return Path(match).resolve()


def read_config_inputs(path: Path) -> list[Path]:
    """Resolve the maintained KLIP configuration's ordered input sequence."""
    config = configparser.ConfigParser(interpolation=None)
    config.optionxform = str
    config.read_string("[root]\n" + path.read_text(encoding="utf-8"))
    section = config["input"]
    file_list = section.get("fileList", "").strip()
    if file_list:
        list_path = Path(file_list).expanduser().resolve()
        paths = [Path(line.strip()).expanduser().resolve() for line in list_path.read_text().splitlines() if line.strip()]
    else:
        directory = Path(section["directory"]).expanduser().resolve()
        paths = sorted(directory.glob(section.get("prefix", "") + "*" + section.get("extension", "")))
    front = int(section.get("deleteFront", "0"))
    back = int(section.get("deleteBack", "0"))
    if back:
        paths = paths[front:-back]
    else:
        paths = paths[front:]
    require(paths and all(path.is_file() for path in paths), "KLIP input sequence is empty or incomplete")
    return paths


def product_paths(manifest: Path) -> list[Path]:
    """List a manifest and every response product belonging to its prefix."""
    prefix = manifest.name.removesuffix("manifest.fits")
    paths = sorted(manifest.parent.glob(prefix + "*.fits"))
    require(manifest in paths, f"manifest is outside its product listing: {manifest}")
    return paths


def read_modes(header: fits.Header) -> list[int]:
    """Read integer KL mode labels from a FITS header."""
    return [int(float(token)) for token in str(header["NMODES"]).split(",")]


def exact_products(manifest: Path) -> dict[str, object]:
    """Validate and resolve the schema-2 exact-response product set."""
    data, header = fits.getdata(manifest, header=True)
    require(np.asarray(data).size == 1 and float(np.ravel(data)[0]) == 1, "exact response manifest is incomplete")
    require(int(header["KLIP PSF PRODUCT SCHEMA"]) == 2, "exact response must use schema 2")
    require(str(header["KLIP PSF SPATIAL MODEL"]).strip().startswith("PIXEL_EXACT"),
            "exact response manifest does not declare PIXEL_EXACT")
    require(int(header["KLIP PSF COMPLETE"]) == 1, "exact response manifest is not complete")
    modes = read_modes(header)
    require(modes == MODES, f"unexpected exact-response modes: {modes}")
    stamp_size = int(header["KLIP PSF STAMP SIZE"])
    require(stamp_size == 11, f"unexpected exact-response stamp size: {stamp_size}")
    prefix = manifest.name.removesuffix("manifest.fits")
    coordinate_path = manifest.parent / (prefix + "coordinates.fits")
    responses = [manifest.parent / (prefix + f"mode{index:03d}_pixel_response.fits") for index in range(len(modes))]
    validities = [manifest.parent / (prefix + f"mode{index:03d}_pixel_validity.fits") for index in range(len(modes))]
    require(all(path.is_file() for path in [coordinate_path, *responses, *validities]),
            "exact response product set is incomplete")
    return {"header": header, "modes": modes, "stamp_size": stamp_size, "coordinates": coordinate_path,
            "responses": responses, "validities": validities}


def science_shape(path: Path) -> tuple[int, int, int]:
    """Read and validate the expected mode-major final-image cube shape."""
    data, header = fits.getdata(path, header=True)
    require(data.ndim == 3 and data.shape[0] == len(MODES), f"unexpected science cube shape: {data.shape}")
    require(read_modes(header) == MODES, "science cube mode ordering changed")
    return tuple(int(value) for value in data.shape)


def protocol_paths(args: argparse.Namespace) -> dict[str, Path]:
    """Resolve every external input used by Stage A."""
    archive = args.archive.resolve()
    sparse_experiment = args.sparse_experiment.resolve()
    exact_manifest = archive / "signal_free_pixel_response/finim_outputs/klipPSF_manifest.fits"
    return {
        "archive": archive,
        "archived_baseline": archive / "signal_free_pixel_response/finim.fits",
        "exact_manifest": exact_manifest,
        "archived_planet_summary": archive / "hciAnalyze/hciAnalyze_filter_summary.csv",
        "archived_provenance": archive / "provenance.txt",
        "original_science": sparse_experiment / "science_only/finim.fits",
        "sparse_manifest": sparse_experiment / "radial_ld_refit4_filter/finim_outputs/klipPSF_manifest.fits",
        "config": args.config.resolve(),
        "psf": args.psf.resolve(),
        "klipreduce": resolve_executable(args.klipreduce),
        "hcianalyze": resolve_executable(args.hcianalyze),
        "runner": Path(__file__).resolve(),
    }


def prepare(args: argparse.Namespace) -> None:
    """Freeze a Stage-A experiment before running the new baseline or analysis."""
    root = args.root.resolve()
    require(not root.exists(), f"output already exists: {root}")
    paths = protocol_paths(args)
    for name, path in paths.items():
        if name not in {"archive"}:
            require(path.is_file(), f"missing {name}: {path}")
    products = exact_products(paths["exact_manifest"])
    archived_shape = science_shape(paths["archived_baseline"])
    require(science_shape(paths["original_science"]) == archived_shape, "archived science cubes have different shapes")
    sparse_data, sparse_header = fits.getdata(paths["sparse_manifest"], header=True)
    require(np.asarray(sparse_data).size == 1 and float(np.ravel(sparse_data)[0]) == 1,
            "sparse response manifest is incomplete")
    require(read_modes(sparse_header) == MODES, "sparse response mode ordering changed")

    input_paths = read_config_inputs(paths["config"])
    root.mkdir(parents=True)
    (root / "software").mkdir()
    shutil.copy2(paths["runner"], root / "software" / paths["runner"].name)
    shutil.copy2(paths["config"], root / "reduction.conf")

    print(f"fingerprinting {len(input_paths)} KLIP inputs", flush=True)
    input_records = [fingerprint(path) for path in input_paths]
    response_paths = product_paths(paths["exact_manifest"]) + product_paths(paths["sparse_manifest"])
    fixed_paths = [paths[name] for name in ("archived_baseline", "archived_planet_summary", "archived_provenance",
                                            "original_science",
                                            "config", "psf", "klipreduce", "hcianalyze")]
    fixed_records = [fingerprint(path) for path in fixed_paths]
    response_records = [fingerprint(path) for path in response_paths]
    copied_records = [fingerprint(root / "software" / paths["runner"].name), fingerprint(root / "reduction.conf")]
    (root / "inputs.txt").write_text("".join(str(path) + "\n" for path in input_paths), encoding="utf-8")
    copied_records.append(fingerprint(root / "inputs.txt"))

    archived_hashes = {}
    for line in paths["archived_provenance"].read_text(encoding="utf-8").splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) == 2 and len(fields[0]) == 64 and all(character in "0123456789abcdef" for character in fields[0]):
            archived_hashes[Path(fields[1]).name] = fields[0]
    current_hashes = {name: fingerprint(paths[name])["sha256"] for name in ("klipreduce", "hcianalyze")}
    protocol = {
        "schema": 1,
        "stage": "KLIP covariance matched-filter Stage A",
        "modes": MODES,
        "primary_mode": 200,
        "primary_radii_pixels": PRIMARY_RADII,
        "radius_bin_half_width_pixels": 0.5,
        "five_pixel_search_offsets_row_column": NEIGHBORS,
        "known_planet": {"separation_pixels": PLANET_SEPARATION, "position_angle_degrees": PLANET_PA,
                         "contrast": PLANET_CONTRAST, "geometry_exclusion_radius_pixels": PLANET_EXCLUSION_RADIUS},
        "baseline_compatibility": {"rtol": 2e-6, "atol": 5e-7, "equal_nan": True},
        "identity_replay": {"rtol": 3e-6, "atol": 2e-5},
        "response_edge_triggers": {"median_primary_mode_bin": 0.01, "selected_location": 0.05},
        "sparse_comparison_samples_per_radius_bin": 64,
        "archived_planet_analysis": {"lambda_d": LAMBDA_D, "planet_radius": ARCHIVED_PLANET_RADIUS,
                                     "snr_min_radius": SNR_MIN_RADIUS, "snr_max_radius": SNR_MAX_RADIUS,
                                     "aperture_radius": ARCHIVED_APERTURE_RADIUS},
        "paths": {name: str(path) for name, path in paths.items() if name != "archive"},
        "archive": str(paths["archive"]),
        "science_shape_fits_order": archived_shape,
        "input_count": len(input_records),
        "resources": {"cpu_affinity": sorted(os.sched_getaffinity(0)),
                      "openmp_threads": len(os.sched_getaffinity(0)), "blas_threads": 1},
        "binary_provenance": {
            "archived": {"klipReduce": archived_hashes.get("klipReduce"),
                         "hciAnalyze": archived_hashes.get("hciAnalyze")},
            "current": {"klipReduce": current_hashes["klipreduce"],
                        "hciAnalyze": current_hashes["hcianalyze"]},
            "hash_matches": {"klipReduce": current_hashes["klipreduce"] == archived_hashes.get("klipReduce"),
                             "hciAnalyze": current_hashes["hcianalyze"] == archived_hashes.get("hciAnalyze")},
        },
    }
    write_json(root / "protocol.json", protocol)
    copied_records.append(fingerprint(root / "protocol.json"))
    manifest = {"schema": 1, "input_records": input_records, "response_records": response_records,
                "fixed_records": fixed_records, "copied_records": copied_records}
    write_json(root / "manifest.json", manifest)
    write_json(root / "state.json", {"status": "prepared", "completed_steps": []})
    print(root, flush=True)


def load_protocol(root: Path) -> tuple[dict[str, object], dict[str, object]]:
    """Load a prepared experiment and verify all immutable files."""
    protocol = json.loads((root / "protocol.json").read_text(encoding="utf-8"))
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    verify(manifest["input_records"] + manifest["response_records"] + manifest["fixed_records"] +
           manifest["copied_records"])
    return protocol, manifest


def run_command(command: list[str], directory: Path, log_path: Path, environment: dict[str, str] | None = None) -> float:
    """Run one recorded command, sending combined output to its log."""
    write_json(log_path.with_suffix(".command.json"), command)
    start = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log:
        subprocess.run(command, cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
    return time.monotonic() - start


def baseline_command(protocol: dict[str, object], root: Path) -> list[str]:
    """Build the current signal-free baseline command from the frozen protocol."""
    paths = protocol["paths"]
    return [str(paths["klipreduce"]), "--config", str(root / "reduction.conf"),
            "--input.directory=", "--input.fileList", str(root / "inputs.txt"),
            "--klip.Nmodes", ",".join(map(str, MODES)),
            "--psfResponse.file=", "--psfResponse.outputModels=false", "--psfResponse.filter=false",
            "--planet.sep", str(PLANET_SEPARATION), "--planet.PA", str(PLANET_PA),
            "--planet.contrast", str(PLANET_CONTRAST), "--fake.method", "single",
            "--fake.fileName", str(paths["psf"]), "--fake.subtractPlanet=true",
            "--output.directory", str(root / "current_baseline"), "--output.fileName", "finim.fits",
            "--output.exactFName=true", "--showTiming=true"]


def run_baseline(protocol: dict[str, object], root: Path) -> tuple[Path, dict[str, object]]:
    """Run or verify the one new signal-free baseline and compare it with the archive."""
    directory = root / "current_baseline"
    complete = directory / "complete.json"
    if complete.is_file():
        record = json.loads(complete.read_text(encoding="utf-8"))
        verify(record["products"])
        return directory / "finim.fits", record["comparison"]
    require(not directory.exists(), "incomplete current_baseline exists; archive it before retrying")
    directory.mkdir()
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS=str(protocol["resources"]["openmp_threads"]), OMP_PROC_BIND="true",
                       OMP_PLACES="cores", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    duration = run_command(baseline_command(protocol, root), directory, directory / "run.log", environment)
    current_path = directory / "finim.fits"
    archived_path = Path(str(protocol["paths"]["archived_baseline"]))
    current, current_header = fits.getdata(current_path, header=True)
    archived, archived_header = fits.getdata(archived_path, header=True)
    require(current.shape == archived.shape, "current and archived signal-free baselines have different shapes")
    require(read_modes(current_header) == read_modes(archived_header) == MODES, "baseline mode labels differ")
    finite_equal = np.array_equal(np.isfinite(current), np.isfinite(archived))
    difference = current.astype(np.float64) - archived.astype(np.float64)
    finite = np.isfinite(difference)
    close = finite_equal and np.allclose(current, archived, rtol=2e-6, atol=5e-7, equal_nan=True)
    comparison = {"compatible": bool(close), "finite_mask_equal": bool(finite_equal),
                  "maximum_absolute_difference": float(np.max(np.abs(difference[finite]))),
                  "rms_difference": float(np.sqrt(np.mean(np.square(difference[finite])))),
                  "current_header_modes": read_modes(current_header), "elapsed_seconds": duration}
    write_json(complete, {"products": [fingerprint(current_path)], "comparison": comparison})
    return current_path, comparison


def coordinates_and_radius(products: dict[str, object], image_shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    """Load exact coordinates and calculate their radii from the half-pixel image center."""
    coordinates = np.asarray(fits.getdata(products["coordinates"]), dtype=np.float64).T
    require(coordinates.ndim == 2 and coordinates.shape[1] == 4, "exact coordinate table has the wrong shape")
    require(np.array_equal(coordinates[:, 3], np.arange(len(coordinates))), "exact coordinate indices are not sequential")
    require(np.all(coordinates[:, :3] == np.trunc(coordinates[:, :3])), "exact coordinates are not integral")
    center_row = 0.5 * (image_shape[1] - 1)
    center_column = 0.5 * (image_shape[0] - 1)
    radius = np.hypot(coordinates[:, 0] - center_row, coordinates[:, 1] - center_column)
    return coordinates, radius


def summarize_values(values: np.ndarray) -> dict[str, object]:
    """Summarize a finite one-dimensional diagnostic sample."""
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not len(values):
        return {"count": 0, "median": None, "p10": None, "p90": None, "maximum": None}
    return {"count": len(values), "median": float(np.median(values)),
            "p10": float(np.percentile(values, 10)), "p90": float(np.percentile(values, 90)),
            "maximum": float(np.max(values))}


def cubic_weight(distance: float) -> float:
    """Evaluate mxlib's default negative-half cubic-convolution kernel."""
    distance = abs(distance)
    cubic = -0.5
    if distance <= 1:
        return (cubic + 2) * distance**3 - (cubic + 3) * distance**2 + 1
    if distance < 2:
        return cubic * distance**3 - 5 * cubic * distance**2 + 8 * cubic * distance - 4 * cubic
    return 0.0


def rotate_response(response: np.ndarray, validity: np.ndarray, angle: float) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce RadialPSFModel's cubic response rotation and validity rule."""
    rows, columns = response.shape
    center_row = 0.5 * (rows - 1)
    center_column = 0.5 * (columns - 1)
    cosine, sine = math.cos(angle), math.sin(angle)
    output = np.zeros_like(response, dtype=np.float64)
    output_validity = np.zeros_like(validity, dtype=bool)
    for output_column in range(columns):
        delta_column = output_column - center_column
        for output_row in range(rows):
            delta_row = output_row - center_row
            input_row = center_row + delta_row * cosine + delta_column * sine
            input_column = center_column - delta_row * sine + delta_column * cosine
            floor_row, floor_column = math.floor(input_row), math.floor(input_column)
            footprint_row, footprint_column = floor_row - 1, floor_column - 1
            row_fraction, column_fraction = input_row - floor_row, input_column - floor_column
            row_weights = [cubic_weight(value) for value in
                           (1 + row_fraction, row_fraction, 1 - row_fraction, 2 - row_fraction)]
            column_weights = [cubic_weight(value) for value in
                              (1 + column_fraction, column_fraction, 1 - column_fraction, 2 - column_fraction)]
            value, valid = 0.0, True
            for column_offset, column_weight in enumerate(column_weights):
                for row_offset, row_weight in enumerate(row_weights):
                    weight = row_weight * column_weight
                    if weight == 0:
                        continue
                    input_sample_row = footprint_row + row_offset
                    input_sample_column = footprint_column + column_offset
                    if not (0 <= input_sample_row < rows and 0 <= input_sample_column < columns):
                        continue
                    if not validity[input_sample_row, input_sample_column]:
                        valid = False
                        break
                    value += response[input_sample_row, input_sample_column] * weight
                if not valid:
                    break
            if valid and math.isfinite(value):
                output[output_row, output_column] = value
                output_validity[output_row, output_column] = True
    return output, output_validity


def evaluate_sparse_response(response_cube: np.ndarray, validity_cube: np.ndarray, radii: np.ndarray,
                             radius: float, angle: float) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate one schema-1 radial response using production interpolation conventions."""
    upper = int(np.searchsorted(radii, radius, side="left"))
    if upper == 0:
        lower = upper = 0
        upper_fraction = 0.0
    elif upper == len(radii):
        lower = upper = len(radii) - 1
        upper_fraction = 0.0
    elif radii[upper] == radius:
        lower = upper
        upper_fraction = 0.0
    else:
        lower = upper - 1
        upper_fraction = (radius - radii[lower]) / (radii[upper] - radii[lower])
    if lower == upper:
        canonical = response_cube[lower].T
        canonical_validity = validity_cube[lower].T
    else:
        canonical = ((1 - upper_fraction) * response_cube[lower] + upper_fraction * response_cube[upper]).T
        canonical_validity = (validity_cube[lower] & validity_cube[upper]).T
    return rotate_response(canonical, canonical_validity, -angle)


def response_audit(protocol: dict[str, object], root: Path) -> tuple[dict[str, object], list[np.ndarray]]:
    """Audit exact-response coordinates, support, energy, lobes, and edge content."""
    manifest = Path(str(protocol["paths"]["exact_manifest"]))
    products = exact_products(manifest)
    shape = tuple(int(value) for value in protocol["science_shape_fits_order"])
    coordinates, radius = coordinates_and_radius(products, shape[1:])
    count = len(coordinates)
    response_cubes: list[np.ndarray] = []
    validity_cubes: list[np.ndarray] = []
    mode_results: dict[str, object] = {}
    border = np.zeros((11, 11), dtype=bool)
    border[[0, -1], :] = True
    border[:, [0, -1]] = True
    for mode, response_path, validity_path in zip(MODES, products["responses"], products["validities"]):
        response = np.asarray(fits.getdata(response_path), dtype=np.float64)
        validity_raw = np.asarray(fits.getdata(validity_path), dtype=np.float64)
        require(response.shape == validity_raw.shape == (count, 11, 11), f"mode {mode} exact product shape changed")
        require(np.all((validity_raw == 0) | (validity_raw == 1)), f"mode {mode} validity is not binary")
        validity = validity_raw > 0.5
        require(np.all(np.isfinite(response[validity])), f"mode {mode} valid response contains nonfinite values")
        energy = np.sum(np.where(validity, response * response, 0), axis=(1, 2))
        require(np.all(energy > 0), f"mode {mode} contains a zero-energy response")
        edge_energy = np.sum(np.where(validity & border[None], response * response, 0), axis=(1, 2)) / energy
        negative_energy = np.sum(np.where(validity & (response < 0), response * response, 0), axis=(1, 2)) / energy
        full = np.all(validity, axis=(1, 2))
        bins = {}
        for nominal in [6.0, *PRIMARY_RADII]:
            selected = np.abs(radius - nominal) <= 0.5
            bins[str(nominal)] = {"responses": int(np.count_nonzero(selected)),
                                  "complete_responses": int(np.count_nonzero(selected & full)),
                                  "energy": summarize_values(energy[selected]),
                                  "negative_energy_fraction": summarize_values(negative_energy[selected]),
                                  "border_energy_fraction": summarize_values(edge_energy[selected])}
        mode_results[str(mode)] = {"complete_response_count": int(np.count_nonzero(full)), "bins": bins}
        response_cubes.append(response)
        validity_cubes.append(validity)

    primary_validity = validity_cubes[MODES.index(200)]
    complete = np.all(primary_validity, axis=(1, 2))
    lookup = {(int(row), int(column)): index for index, (row, column, _, _) in enumerate(coordinates)}
    planet_row = 0.5 * (shape[2] - 1) - PLANET_SEPARATION * math.sin(math.radians(PLANET_PA))
    planet_column = 0.5 * (shape[1] - 1) + PLANET_SEPARATION * math.cos(math.radians(PLANET_PA))
    geometry = {}
    for nominal in [6.0, *PRIMARY_RADII]:
        native = []
        for index, (row_value, column_value, _, _) in enumerate(coordinates):
            row, column = int(row_value), int(column_value)
            if abs(radius[index] - nominal) > 0.5:
                continue
            response_indices = [lookup.get((row + delta_row, column + delta_column))
                                for delta_row, delta_column in NEIGHBORS]
            if any(value is None or not complete[int(value)] for value in response_indices):
                continue
            if any(math.hypot(row + delta_row - planet_row, column + delta_column - planet_column) <=
                   PLANET_EXCLUSION_RADIUS for delta_row, delta_column in NEIGHBORS):
                continue
            native.append({"row": row, "column": column, "angle_radians": math.atan2(row - 63.5, column - 63.5)})
        geometry[str(nominal)] = {"eligible_five_pixel_centers": len(native), "centers": native}
    require([geometry[str(radius)]["eligible_five_pixel_centers"] for radius in PRIMARY_RADII] ==
            [42, 44, 60, 91, 113, 160], "exact-response primary geometry no longer matches the preregistration")

    edge_failures = []
    for nominal in PRIMARY_RADII:
        summary = mode_results["200"]["bins"][str(nominal)]["border_energy_fraction"]
        if summary["median"] is not None and summary["median"] > 0.01:
            edge_failures.append({"radius": nominal, "median_border_energy_fraction": summary["median"]})
    result = {"coordinate_count": count, "mode_results": mode_results, "geometry": geometry,
              "known_planet_row_column": [planet_row, planet_column],
              "primary_edge_trigger_passed": not edge_failures, "primary_edge_trigger_failures": edge_failures,
              "selected_location_edge_trigger_deferred": True}
    write_json(root / "response_audit.json", result)
    return result, response_cubes


def sparse_response_audit(protocol: dict[str, object], root: Path,
                          exact_responses: list[np.ndarray]) -> dict[str, object]:
    """Compare exact responses with deterministic samples of the sparse radial model."""
    exact = exact_products(Path(str(protocol["paths"]["exact_manifest"])))
    coordinates = np.asarray(fits.getdata(exact["coordinates"]), dtype=np.float64).T
    center = 0.5 * (int(protocol["science_shape_fits_order"][1]) - 1)
    source_radii = np.hypot(coordinates[:, 0] - center, coordinates[:, 1] - center)
    source_angles = np.arctan2(coordinates[:, 0] - center, coordinates[:, 1] - center)
    sparse_manifest = Path(str(protocol["paths"]["sparse_manifest"]))
    sparse_prefix = sparse_manifest.name.removesuffix("manifest.fits")
    results: dict[str, object] = {}
    for mode_index, mode in enumerate(MODES):
        response_path = sparse_manifest.parent / (sparse_prefix + f"mode{mode_index:03d}_radial_response.fits")
        validity_path = sparse_manifest.parent / (sparse_prefix + f"mode{mode_index:03d}_radial_validity.fits")
        sparse_cube, header = fits.getdata(response_path, header=True)
        sparse_cube = np.asarray(sparse_cube, dtype=np.float64)
        sparse_validity = np.asarray(fits.getdata(validity_path), dtype=np.float64) > 0.5
        radial_nodes = np.asarray([float(token) for token in str(header["KLIP PSF SAMPLE RADII"]).split(",")])
        require(sparse_cube.shape == sparse_validity.shape == (len(radial_nodes), 11, 11),
                f"mode {mode} sparse response shape changed")
        exact_validity = np.asarray(fits.getdata(exact["validities"][mode_index]), dtype=np.float64) > 0.5
        mode_bins = {}
        for nominal in [6.0, *PRIMARY_RADII]:
            candidates = np.flatnonzero((np.abs(source_radii - nominal) <= 0.5) &
                                        np.all(exact_validity, axis=(1, 2)))
            order = candidates[np.argsort(source_angles[candidates])]
            if len(order) > 64:
                selected_indices = np.floor(np.arange(64) * len(order) / 64).astype(int)
                order = order[selected_indices]
            cosines, projections, residuals = [], [], []
            for source in order:
                sparse_response, sparse_support = evaluate_sparse_response(
                    sparse_cube, sparse_validity, radial_nodes, float(source_radii[source]),
                    float(source_angles[source]))
                exact_response = exact_responses[mode_index][source].T
                support = sparse_support & exact_validity[source].T & np.isfinite(sparse_response) & np.isfinite(exact_response)
                sparse_values, exact_values = sparse_response[support], exact_response[support]
                sparse_energy = float(sparse_values @ sparse_values)
                exact_energy = float(exact_values @ exact_values)
                if sparse_energy <= 0 or exact_energy <= 0:
                    continue
                cross = float(sparse_values @ exact_values)
                projection = cross / sparse_energy
                cosines.append(cross / math.sqrt(sparse_energy * exact_energy))
                projections.append(projection)
                residuals.append(math.sqrt(float(np.sum(np.square(exact_values - projection * sparse_values))) /
                                           exact_energy))
            mode_bins[str(nominal)] = {"sampled_locations": len(order), "valid_comparisons": len(cosines),
                                       "cosine": summarize_values(np.asarray(cosines)),
                                       "sparse_to_exact_projection_scale": summarize_values(np.asarray(projections)),
                                       "best_scaled_relative_residual": summarize_values(np.asarray(residuals))}
        results[str(mode)] = mode_bins
    result = {"sampling": "up to 64 complete exact responses per radial bin, evenly indexed in sorted angle",
              "comparison_orientation": "sparse model projected onto exact response on common valid support",
              "modes": results}
    write_json(root / "sparse_response_audit.json", result)
    return result


def identity_amplitude_cube(science: np.ndarray, products: dict[str, object], responses: list[np.ndarray]) -> np.ndarray:
    """Independently apply the signed normalized exact-response filter."""
    coordinates = np.asarray(fits.getdata(products["coordinates"]), dtype=np.float64).T
    output = np.full(science.shape, np.nan, dtype=np.float32)
    half = int(products["stamp_size"]) // 2
    for mode_index, (response, validity_path) in enumerate(zip(responses, products["validities"])):
        validity = np.asarray(fits.getdata(validity_path), dtype=np.float64) > 0.5
        for source, (row_value, column_value, _, _) in enumerate(coordinates):
            row, column = int(row_value), int(column_value)
            if not validity[source, half, half]:
                continue
            # FITS images are (column,row), while native response stamps are
            # exposed in FITS as (column,row) and become Eigen-oriented here.
            stamp = science[mode_index, column - half:column + half + 1, row - half:row + half + 1].T
            if stamp.shape != (2 * half + 1, 2 * half + 1):
                continue
            template = response[source].T
            support = validity[source].T & np.isfinite(template) & np.isfinite(stamp)
            if np.count_nonzero(support) / support.size < 1:
                continue
            template_values = template[support]
            denominator = float(template_values @ template_values)
            if denominator > 0:
                output[mode_index, column, row] = float(template_values @ stamp[support] / denominator)
    return output


def annular_oracle(amplitude: np.ndarray, source_x: float, source_y: float,
                   source_radius: float) -> np.ndarray:
    """Reconstruct hciAnalyze's annular normalization and small-sample factor."""
    yy, xx = np.indices(amplitude.shape)
    center_x = 0.5 * (amplitude.shape[1] - 1)
    center_y = 0.5 * (amplitude.shape[0] - 1)
    radius = np.hypot(xx - center_x, yy - center_y).astype(np.float32)
    noise = np.isfinite(amplitude) & (np.hypot(xx - source_x, yy - source_y).astype(np.float32) > source_radius + 0.5)
    maximum_bin = int(math.ceil(float(np.max(radius)))) + 1
    centers, means, deviations = [], [], []
    for lower in range(maximum_bin):
        values = amplitude[noise & (radius > lower) & (radius <= lower + 1)].astype(np.float64)
        centers.append(lower + 0.5)
        means.append(float(np.mean(values)) if len(values) else np.nan)
        deviations.append(float(np.std(values, ddof=1)) if len(values) > 1 else np.nan)
    mean = np.interp(radius, centers, means)
    deviation = np.interp(radius, centers, deviations).astype(np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        score = ((amplitude - mean) / deviation).astype(np.float32)
        count = 2 * np.pi * radius / LAMBDA_D - 1
        correction = np.where(count > 0, 1 / np.sqrt(1 + 1 / count), 0).astype(np.float32)
        score *= correction
    score[(radius < SNR_MIN_RADIUS) | (radius > SNR_MAX_RADIUS) | ~np.isfinite(amplitude)] = np.nan
    return score


def analyze_command(binary: str, science: Path, response: str = "", gaussian: float = 0.0) -> list[str]:
    """Construct the archived known-planet hciAnalyze command."""
    return [binary, "--file", str(science), "--lambdaD", str(LAMBDA_D), "--planet.sep", str(PLANET_SEPARATION),
            "--planet.PA", str(PLANET_PA), "--planet.contrast", str(PLANET_CONTRAST),
            "--planet.R", str(ARCHIVED_PLANET_RADIUS), "--snr.minRad", str(SNR_MIN_RADIUS),
            "--snr.maxRad", str(SNR_MAX_RADIUS), "--snr.apertureR", str(ARCHIVED_APERTURE_RADIUS),
            "--filter.hpfGaussFW", "0", "--filter.lpfGaussFW", str(gaussian),
            "--filter.psfResponse", response, "--noise.model=identity", "--noise.outputDiagnostics=false"]


def parse_hcianalyze(path: Path) -> list[dict[str, object]]:
    """Parse hciAnalyze's whitespace result table."""
    lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    require(lines and lines[0] == ["plane", "mode", "signal", "separation", "position_angle", "contrast", "snr"],
            f"unexpected hciAnalyze result format: {path}")
    return [{"plane": int(row[0]), "mode": int(float(row[1])), "signal": int(row[2]), "snr": float(row[6])}
            for row in lines[1:]]


def exact_identity_replay(protocol: dict[str, object], root: Path, baseline: Path,
                          responses: list[np.ndarray]) -> dict[str, object]:
    """Compare the independent exact identity filter with direct hciAnalyze output."""
    directory = root / "identity_replay"
    complete = directory / "complete.json"
    if complete.is_file():
        result = json.loads(complete.read_text(encoding="utf-8"))
        verify(result["products"])
        require(result["passed"], "completed identity replay did not pass")
        return result
    require(not directory.exists(), "incomplete identity_replay exists; archive it before retrying")
    directory.mkdir()
    science, header = fits.getdata(baseline, header=True)
    products = exact_products(Path(str(protocol["paths"]["exact_manifest"])))
    amplitude = identity_amplitude_cube(np.asarray(science, dtype=np.float64), products, responses)
    fits.writeto(directory / "independent_amplitude.fits", amplitude, header)
    direct_input = directory / "science.fits"
    shutil.copy2(baseline, direct_input)
    command = analyze_command(str(protocol["paths"]["hcianalyze"]), direct_input,
                              str(protocol["paths"]["exact_manifest"]), 0.0)
    elapsed = run_command(command, directory, directory / "hciAnalyze.log")
    actual = np.asarray(fits.getdata(directory / "science_snr.fits"), dtype=np.float32)
    source_x = 0.5 * (amplitude.shape[2] - 1) - PLANET_SEPARATION * math.sin(math.radians(PLANET_PA))
    source_y = 0.5 * (amplitude.shape[1] - 1) + PLANET_SEPARATION * math.cos(math.radians(PLANET_PA))
    expected = np.stack([annular_oracle(plane, source_x, source_y, ARCHIVED_PLANET_RADIUS) for plane in amplitude])
    fits.writeto(directory / "independent_snr.fits", expected, header)
    yy, xx = np.indices(amplitude.shape[1:])
    radius = np.hypot(xx - 0.5 * (amplitude.shape[2] - 1), yy - 0.5 * (amplitude.shape[1] - 1))
    checked = np.isfinite(expected) & np.isfinite(amplitude) & (radius[None] >= SNR_MIN_RADIUS) & (radius[None] <= SNR_MAX_RADIUS)
    require(np.any(checked), "identity replay has no comparable pixels")
    difference = actual[checked] - expected[checked]
    close = np.allclose(actual[checked], expected[checked], rtol=3e-6, atol=2e-5)
    result = {"passed": bool(close), "compared_pixels": int(np.count_nonzero(checked)),
              "maximum_absolute_snr_difference": float(np.max(np.abs(difference))),
              "rms_snr_difference": float(np.sqrt(np.mean(np.square(difference)))), "elapsed_seconds": elapsed,
              "products": [fingerprint(directory / name) for name in
                           ("independent_amplitude.fits", "independent_snr.fits", "science_snr.fits")]}
    write_json(directory / "complete.json", result)
    require(close, "independent exact identity-filter replay differs from hciAnalyze")
    return result


def planet_snr_reproduction(protocol: dict[str, object], root: Path) -> dict[str, object]:
    """Rerun the four archived known-planet filter controls with current hciAnalyze."""
    directory = root / "planet_controls"
    complete = directory / "complete.json"
    if complete.is_file():
        result = json.loads(complete.read_text(encoding="utf-8"))
        verify(result["products"])
        require(result["passed"], "completed planet-control reproduction did not pass")
        return result
    require(not directory.exists(), "incomplete planet_controls exists; archive it before retrying")
    directory.mkdir()
    cases = {
        "unfiltered": ("", 0.0),
        "sparse_response": (str(protocol["paths"]["sparse_manifest"]), 0.0),
        "exact_response": (str(protocol["paths"]["exact_manifest"]), 0.0),
        "gaussian_fwhm_3p6": ("", 3.6),
    }
    rows = []
    products = []
    for name, (response, gaussian) in cases.items():
        case = directory / name
        case.mkdir()
        science = case / "science.fits"
        shutil.copy2(Path(str(protocol["paths"]["original_science"])), science)
        command = analyze_command(str(protocol["paths"]["hcianalyze"]), science, response, gaussian)
        run_command(command, case, case / "results.txt")
        parsed = parse_hcianalyze(case / "results.txt")
        require([row["mode"] for row in parsed] == MODES, f"{name} returned unexpected mode rows")
        for row in parsed:
            rows.append({"case": name, **row})
        products.extend([fingerprint(case / "results.txt"), fingerprint(case / "science_snr.fits")])

    archived = {}
    with Path(str(protocol["paths"]["archived_planet_summary"])).open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            archived[(row["case"], int(row["mode"]))] = float(row["snr"])
    comparisons = []
    for row in rows:
        reference = archived[(str(row["case"]), int(row["mode"]))]
        difference = float(row["snr"]) - reference
        comparisons.append({**row, "archived_snr": reference, "difference": difference})
    maximum = max(abs(float(row["difference"])) for row in comparisons)
    passed = all(np.isclose(float(row["snr"]), float(row["archived_snr"]), rtol=2e-5, atol=2e-5)
                 for row in comparisons)
    result = {"passed": bool(passed), "maximum_absolute_snr_difference": maximum,
              "comparisons": comparisons, "products": products}
    write_json(directory / "complete.json", result)
    require(passed, "current hciAnalyze does not reproduce the archived planet controls")
    return result


def write_report(root: Path, protocol: dict[str, object], baseline: dict[str, object], response: dict[str, object],
                 sparse: dict[str, object], replay: dict[str, object], planet: dict[str, object]) -> None:
    """Write a compact Markdown report for review before Stage B."""
    lines = ["# KLIP covariance Stage A", "",
             "Stage A verifies the archived exact response and current analysis path before covariance fitting.", "",
             "## Gates", "", "| Check | Result | Diagnostic |", "| --- | --- | --- |",
             f"| klipReduce archive hash | {'pass' if protocol['binary_provenance']['hash_matches']['klipReduce'] else 'different'} | current versus archived executable |",
             f"| hciAnalyze archive hash | {'pass' if protocol['binary_provenance']['hash_matches']['hciAnalyze'] else 'different'} | current versus archived executable |",
             f"| Current signal-free baseline | {'pass' if baseline['compatible'] else 'fail'} | max abs difference {baseline['maximum_absolute_difference']:.6g} |",
             f"| Exact-response geometry | pass | {response['coordinate_count']} native response locations |",
             f"| Independent identity-filter replay | {'pass' if replay['passed'] else 'fail'} | max SNR difference {replay['maximum_absolute_snr_difference']:.6g} |",
             f"| Archived planet controls | {'pass' if planet['passed'] else 'fail'} | max SNR difference {planet['maximum_absolute_snr_difference']:.6g} |",
             f"| 11x11 response edge-energy trigger | {'pass' if response['primary_edge_trigger_passed'] else 'triggered'} | mode 200, primary radial bins |",
             "", "## Five-pixel common-support geometry", "", "| Radius (pixels) | Eligible centers |",
             "| ---: | ---: |"]
    for radius in [6.0, *PRIMARY_RADII]:
        lines.append(f"| {radius:g} | {response['geometry'][str(radius)]['eligible_five_pixel_centers']} |")
    lines.extend(["", "## Exact versus sparse response at mode 200", "",
                  "| Radius (pixels) | Samples | Median cosine | Median projection | Median relative residual |",
                  "| ---: | ---: | ---: | ---: | ---: |"])
    for radius in [6.0, *PRIMARY_RADII]:
        summary = sparse["modes"]["200"][str(radius)]
        lines.append(f"| {radius:g} | {summary['valid_comparisons']} | {summary['cosine']['median']:.5f} | "
                     f"{summary['sparse_to_exact_projection_scale']['median']:.5f} | "
                     f"{summary['best_scaled_relative_residual']['median']:.5f} |")
    lines.extend(["", "## Known-planet SNR reproduction", "",
                  "| Mode | Unfiltered | Gaussian 3.6 | Exact response | Sparse response |", "| ---: | ---: | ---: | ---: | ---: |"])
    lookup = {(row["case"], row["mode"]): row["snr"] for row in planet["comparisons"]}
    for mode in MODES:
        lines.append(f"| {mode} | {lookup[('unfiltered', mode)]:.5f} | {lookup[('gaussian_fwhm_3p6', mode)]:.5f} | "
                     f"{lookup[('exact_response', mode)]:.5f} | {lookup[('sparse_response', mode)]:.5f} |")
    if not response["primary_edge_trigger_passed"]:
        lines.extend(["", "The preregistered response edge-energy trigger fired. A larger-stamp response experiment is required before an optimal-filter claim."])
    (root / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> None:
    """Run or resume the prepared Stage-A experiment."""
    root = args.root.resolve()
    protocol, _ = load_protocol(root)
    require(sorted(os.sched_getaffinity(0)) == protocol["resources"]["cpu_affinity"],
            "CPU affinity changed between prepare and run")
    write_json(root / "state.json", {"status": "running", "pid": os.getpid(), "completed_steps": []})
    try:
        baseline_path, baseline = run_baseline(protocol, root)
        write_json(root / "state.json", {"status": "running", "pid": os.getpid(), "completed_steps": ["baseline"]})
        response, response_cubes = response_audit(protocol, root)
        sparse = sparse_response_audit(protocol, root, response_cubes)
        write_json(root / "state.json", {"status": "running", "pid": os.getpid(),
                                         "completed_steps": ["baseline", "response_audit", "sparse_response_audit"]})
        replay_baseline = baseline_path if baseline["compatible"] else Path(str(protocol["paths"]["archived_baseline"]))
        replay = exact_identity_replay(protocol, root, replay_baseline, response_cubes)
        replay["baseline_source"] = "current" if baseline["compatible"] else "archived exact-response baseline"
        write_json(root / "state.json", {"status": "running", "pid": os.getpid(),
                                         "completed_steps": ["baseline", "response_audit", "sparse_response_audit",
                                                             "identity_replay"]})
        planet = planet_snr_reproduction(protocol, root)
        results = {"baseline": baseline, "response_audit": response, "sparse_response_audit": sparse,
                   "identity_replay": replay,
                   "planet_controls": planet,
                   "stage_a_passed": bool(baseline["compatible"] and replay["passed"] and planet["passed"]),
                   "larger_stamp_experiment_required": not response["primary_edge_trigger_passed"]}
        write_json(root / "results.json", results)
        write_report(root, protocol, baseline, response, sparse, replay, planet)
        completion = {"status": "complete", "results": fingerprint(root / "results.json"),
                      "report": fingerprint(root / "results.md")}
        write_json(root / "complete.json", completion)
        write_json(root / "state.json", {"status": "complete",
                                         "completed_steps": ["baseline", "response_audit", "sparse_response_audit",
                                                             "identity_replay", "planet_controls", "report"]})
        print(root / "results.md", flush=True)
    except Exception as error:
        write_json(root / "state.json", {"status": "failed", "error_type": type(error).__name__, "error": str(error)})
        raise


def check() -> None:
    """Run small deterministic checks of the JSON and identity-filter adapters."""
    template = np.arange(-60, 61, dtype=np.float64).reshape(11, 11)
    template[5, 5] = 7
    data = 2.5 * template + 3
    validity = np.ones((11, 11), dtype=bool)
    amplitude = float(template[validity] @ data[validity] / (template[validity] @ template[validity]))
    expected = 2.5 + 3 * float(np.sum(template)) / float(template.ravel() @ template.ravel())
    require(np.isclose(amplitude, expected, rtol=1e-14, atol=1e-14), "identity-filter arithmetic check failed")
    encoded = json.dumps(finite_json({"nan": np.nan, "integer": np.int64(3)}), allow_nan=False)
    require(encoded == '{"nan": null, "integer": 3}', "strict JSON adapter check failed")
    print("KLIP covariance Stage-A runner checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    repo = Path(__file__).resolve().parents[3]
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check", help="run deterministic adapter checks")
    for action in ("prepare", "run"):
        command = subparsers.add_parser(action)
        command.add_argument("root", type=Path, help="Stage-A experiment directory")
        if action == "prepare":
            command.add_argument("--archive", type=Path,
                                 default=repo / "working/roc/klip_signal_free_pixel_response_20260914T232006Z")
            command.add_argument("--sparse-experiment", type=Path,
                                 default=repo / "working/roc/klip_cpp_response_20260913T222913Z")
            command.add_argument("--config", type=Path,
                                 default=repo / "agents/plans/scripts/klipReduce_afLepNaco_psf_response.conf")
            command.add_argument("--psf", type=Path,
                                 default=Path("/home/jrmales/Source/mxWork/NACO/AFLep/2011-10-21/out/psf_reg_median.fits"))
            command.add_argument("--klipreduce", default="klipReduce")
            command.add_argument("--hcianalyze", default="hciAnalyze")
    return result


def main() -> None:
    """Dispatch the selected Stage-A action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    elif args.action == "prepare":
        prepare(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
