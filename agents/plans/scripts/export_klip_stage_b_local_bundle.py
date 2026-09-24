#!/usr/bin/env python3
"""Export a compact mode-200 KLIP response bundle for local Stage-B work."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from astropy.io import fits

SCRIPT_DIRECTORY = Path(__file__).resolve().parent
REPOSITORY_SCRIPTS = Path.cwd() / "agents" / "plans" / "scripts"
sys.path.insert(0, str(SCRIPT_DIRECTORY))
sys.path.insert(0, str(REPOSITORY_SCRIPTS))
import run_klip_covariance_stage_a as stage  # noqa: E402
import run_klip_response_47 as response47  # noqa: E402
import run_klip_stage_b_footprint_preflight as footprint  # noqa: E402


MODE = 200
SUPPORTS = (11, 31, 47)
RADII = (7.5, 10.0, 12.0, 16.0, 20.0, 24.0)


def read(path: Path) -> object:
    """Read one JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def selected_positions(result: dict[str, object]) -> np.ndarray:
    """Return the sorted union of fixed five-pixel query positions."""
    positions = set()
    for support in SUPPORTS:
        for radius in RADII:
            sites = result["geometry"][str(support)][str(radius)]["selected_sites"]
            for site in sites:
                row, column = int(site["row"]), int(site["column"])
                positions.update((row + delta_row, column + delta_column)
                                 for delta_row, delta_column in footprint.SEARCH_OFFSETS)
    stage.require(positions, "decoupled preflight contains no selected positions")
    return np.asarray(sorted(positions), dtype=np.int16)


def product_record(receipt: dict[str, object], path: Path) -> dict[str, object]:
    """Find one product fingerprint in the response receipt."""
    resolved = str(path.resolve())
    matches = [record for record in receipt["products"] if record["path"] == resolved]
    stage.require(len(matches) == 1, f"response receipt does not contain {resolved}")
    return matches[0]


def export(args: argparse.Namespace) -> None:
    """Write and verify one compressed local-analysis bundle."""
    response_root = args.response.resolve()
    preflight_root = args.preflight.resolve()
    output = args.output.resolve()
    stage.require(not output.exists() and not output.with_suffix(output.suffix + ".receipt.json").exists(),
                  f"output already exists: {output}")

    response_completion = read(response_root / "complete.json")
    response_results = read(response_root / "results.json")
    stage.require(response_completion["status"] == "complete" and
                  response_results["central_11_archive_replay_passed"] and
                  response_results["external_47_replay_passed"], "response campaign did not pass")
    stage.verify([response_completion["response"], response_completion["results"],
                  response_completion["report"]])
    response_receipt = read(response_root / "response" / "complete.json")

    preflight_completion = read(preflight_root / "complete.json")
    stage.require(preflight_completion["status"] == "complete", "decoupled preflight is incomplete")
    stage.verify([preflight_completion["results"], preflight_completion["report"]])
    preflight_result = read(preflight_root / "results.json")
    stage.require(not preflight_result["selection_uses_baseline_values"],
                  "decoupled site selection used baseline values")

    paths = response47.product_paths(response_root)
    mode_index = stage.MODES.index(MODE)
    required = [paths["final"], paths["coordinates"], paths["responses"][mode_index],
                paths["validities"][mode_index]]
    stage.verify([product_record(response_receipt, path) for path in required])

    coordinates = np.asarray(fits.getdata(paths["coordinates"]), dtype=np.float64).T
    lookup = {(int(row), int(column)): index
              for index, (row, column, _, _) in enumerate(coordinates)}
    positions = selected_positions(preflight_result)
    source_indices = np.asarray([lookup.get((int(row), int(column)), -1)
                                 for row, column in positions], dtype=np.int32)
    stage.require(np.all(source_indices >= 0), "selected position is absent from the exact response field")

    baseline_cube = np.asarray(fits.getdata(paths["final"], memmap=True), dtype=np.float32)
    response = np.asarray(fits.getdata(paths["responses"][mode_index], memmap=True)[source_indices],
                          dtype=np.float32)
    validity = np.asarray(fits.getdata(paths["validities"][mode_index], memmap=True)[source_indices] > 0.5,
                          dtype=np.uint8)
    stage.require(baseline_cube.shape[0] == len(stage.MODES) and
                  response.shape == validity.shape == (len(positions), 47, 47),
                  "bundle source arrays have unexpected shapes")
    stage.require(np.all(validity[:, 23, 23] == 1), "bundle contains an invalid response anchor")

    metadata = {
        "schema": 1,
        "purpose": "compact local mode-200 KLIP Stage-B noise-screen input",
        "mode": MODE,
        "mode_index": mode_index,
        "supports": list(SUPPORTS),
        "radii": list(RADII),
        "search_offsets": [list(offset) for offset in footprint.SEARCH_OFFSETS],
        "position_count": len(positions),
        "response_root": str(response_root),
        "preflight_root": str(preflight_root),
        "response_complete_sha256": stage.fingerprint(response_root / "complete.json")["sha256"],
        "preflight_complete_sha256": stage.fingerprint(preflight_root / "complete.json")["sha256"],
        "source_product_records": [product_record(response_receipt, path) for path in required],
        "known_planet": read(response_root / "protocol.json")["known_planet"],
        "selection": "union of the fixed decoupled-preflight sites and their five search pixels",
        "geometry": preflight_result["geometry"],
    }
    encoded = np.frombuffer(json.dumps(metadata, sort_keys=True).encode("utf-8"), dtype=np.uint8)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, baseline=baseline_cube[mode_index], positions=positions,
                        source_indices=source_indices, response=response, validity=validity,
                        metadata_json=encoded)

    with np.load(output, allow_pickle=False) as bundle:
        recovered = json.loads(bytes(bundle["metadata_json"]).decode("utf-8"))
        stage.require(recovered == metadata and
                      bundle["baseline"].shape == baseline_cube[mode_index].shape and
                      np.array_equal(bundle["positions"], positions) and
                      np.array_equal(bundle["source_indices"], source_indices) and
                      np.array_equal(bundle["response"], response) and
                      np.array_equal(bundle["validity"], validity),
                      "written local bundle does not replay")
    receipt = {"status": "complete", "bundle": stage.fingerprint(output),
               "position_count": len(positions), "mode": MODE,
               "source_records": metadata["source_product_records"]}
    stage.write_json(output.with_suffix(output.suffix + ".receipt.json"), receipt)
    print(json.dumps({"bundle": receipt["bundle"], "positions": len(positions)}, indent=2), flush=True)


def check() -> None:
    """Check deterministic position union behavior."""
    result = {"geometry": {str(support): {str(radius): {"selected_sites": [
        {"row": 20 + int(radius), "column": 30 + support}]} for radius in RADII}
        for support in SUPPORTS}}
    first = selected_positions(result)
    second = selected_positions(result)
    stage.require(np.array_equal(first, second) and len(first) <= len(SUPPORTS) * len(RADII) * 5,
                  "local-bundle position selection is not deterministic")
    print("KLIP Stage-B local-bundle checks passed", flush=True)


def parser() -> argparse.ArgumentParser:
    """Build the local-bundle command-line parser."""
    repo = Path.cwd()
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="action", required=True)
    subparsers.add_parser("check")
    export_parser = subparsers.add_parser("export")
    export_parser.add_argument("output", type=Path)
    export_parser.add_argument("--response", type=Path,
                               default=repo / "working/roc/klip_response_47_20260921")
    export_parser.add_argument("--preflight", type=Path,
                               default=repo / "working/roc/klip_stage_b_decoupled_preflight_20260923")
    return result


def main() -> None:
    """Dispatch the requested bundle action."""
    args = parser().parse_args()
    if args.action == "check":
        check()
    else:
        export(args)


if __name__ == "__main__":
    main()
