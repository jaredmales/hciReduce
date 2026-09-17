#!/usr/bin/env python3
"""Check the production analytic P4 kernel against independent SVD and FP64 paired refits.

Replay baseline/unit-source captures from a completed detector convergence sweep.
The comparison uses the same sampled unit direction for all three methods; paired
refits here perturb the captured FP64 arrays, without resampling a source image.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from analyze_p4_detector_convergence import capture_paths, read_capture, reference_derivative, same_geometry, write_csv
from run_p4_response_convergence import fingerprint, ratio


def main() -> None:
    """Preserve native replay products and compare only explicitly resolved/supported response columns."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep", required=True, type=Path)
    parser.add_argument("--binary", required=True, type=Path, help="p4ResponseBenchmark executable")
    parser.add_argument("--output", required=True, type=Path, help="new analysis directory")
    parser.add_argument("--amplitude", type=float, default=1e-5, help="FP64 refit half-amplitude")
    args = parser.parse_args()
    if not np.isfinite(args.amplitude) or args.amplitude <= 0:
        parser.error("amplitude must be finite and positive")
    root, binary, output = args.sweep.resolve(), args.binary.resolve(), args.output.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    settings = manifest["arguments"]
    if (not (root / "complete.json").is_file() or settings.get("precision") != "P4-D64"
            or not settings.get("capture_detectors")):
        parser.error("a completed FP64 detector-capture sweep is required")
    libraries = {}
    for line in subprocess.check_output(["ldd", str(binary)], text=True).splitlines():
        if "=>" not in line:
            continue
        name, resolved = line.split("=>", 1)
        if name.strip().startswith(("libhcireduce.", "libmxlib.", "libopenblas.", "liblapack.")):
            libraries[name.strip()] = fingerprint(Path(resolved.split(" (", 1)[0].strip()).resolve())
    output.mkdir(parents=True, exist_ok=False)
    (output / "provenance.json").write_text(json.dumps({
        "schema": 1, "sweep": fingerprint(root / "manifest.json"), "binary": fingerprint(binary),
        "libraries": libraries,
        "scripts": [fingerprint(Path(__file__).resolve().with_name(name)) for name in
                    (Path(__file__).name, "analyze_p4_detector_convergence.py", "run_p4_response_convergence.py")],
        "amplitude": args.amplitude,
        "scope": "captured FP64 inputs; source direction = unit minus baseline; no further source sampling",
        "thread_environment": {key: os.environ.get(key) for key in
                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    }, indent=2) + "\n")
    rows, summaries = [], []
    positions = [(r, a) for r in settings["radii"] for a in settings["angles"]]
    for position, (radius, angle) in enumerate(positions):
        source = root / f"position_{position:03d}"
        destination = output / source.name
        destination.mkdir()
        baseline_paths = capture_paths(source / "baseline")
        unit_paths = capture_paths(source / "unit_source")
        if baseline_paths.keys() != unit_paths.keys():
            raise ValueError("baseline and unit-source capture coordinates differ")
        norms = np.zeros((len(settings["modes"]), 3))
        counts = np.zeros(len(settings["modes"]), dtype=int)
        for name, baseline_path in sorted(baseline_paths.items()):
            baseline, unit = read_capture(baseline_path), read_capture(unit_paths[name])
            same_geometry(baseline, unit)
            metadata = baseline["metadata"]
            if not np.dtype(metadata["dtype"]).isnative or not np.dtype(unit["metadata"]["dtype"]).isnative:
                raise ValueError("native C++ replay requires native-endian capture files")
            reference, diagnostics = reference_derivative(baseline, unit)
            replay = destination / f"{name}.bin"
            command = [str(binary), str(metadata["rows"]), str(metadata["columns"]), str(metadata["rank_tolerance"]),
                       ",".join(map(str, metadata["modes"])), str(args.amplitude),
                       str(baseline_path.with_suffix(".bin")), str(unit_paths[name].with_suffix(".bin")), str(replay)]
            record = {"command": command,
                      "captures": [fingerprint(path) for path in
                                   (baseline_path, baseline_path.with_suffix(".bin"),
                                    unit_paths[name], unit_paths[name].with_suffix(".bin"))]}
            with (destination / f"{name}.log").open("w") as log:
                process = subprocess.run(command, capture_output=True, text=True, check=False)
                log.write(process.stdout + process.stderr)
            process.check_returncode()
            result = json.loads(process.stdout)
            record["result"] = result
            (destination / f"{name}.json").write_text(json.dumps(record, indent=2) + "\n")
            shape = (metadata["rows"], len(metadata["modes"]))
            data = np.fromfile(replay, dtype=np.float64)
            elements = shape[0] * shape[1]
            if data.size != 2 * elements or result["schema"] != 1 or len(result["modes"]) != shape[1]:
                raise ValueError("invalid analytic replay dimensions/schema")
            analytic = data[:elements].reshape(shape, order="F")
            paired = data[elements:].reshape(shape, order="F")
            for index, mode in enumerate(result["modes"]):
                if mode["count"] != metadata["modes"][index]:
                    raise ValueError("replay returned a different retained count")
                available = mode["status"] == "differentiable"
                if available and not np.all(np.isfinite(analytic[:, index])):
                    raise ValueError("nonfinite available analytic response")
                if not available and not np.all(np.isnan(analytic[:, index])):
                    raise ValueError("unavailable analytic response must be NaN")
                row = {"position": position, "fit": name, "mode_fraction": settings["modes"][index],
                       "retained_modes": mode["count"], "analytic_status": mode["status"],
                       "reference_resolved": diagnostics[index]["reference_resolved"],
                       "paired_supported": mode["paired_supported"],
                       "analytic_relative_error": None, "paired_relative_error": None}
                if available and diagnostics[index]["reference_resolved"] and mode["paired_supported"]:
                    if not np.all(np.isfinite(paired[:, index])):
                        raise ValueError("nonfinite supported paired response")
                    norm = np.linalg.norm(reference[:, index])
                    analytic_error = np.linalg.norm(analytic[:, index] - reference[:, index])
                    paired_error = np.linalg.norm(paired[:, index] - reference[:, index])
                    row.update(analytic_relative_error=ratio(analytic_error, norm),
                               paired_relative_error=ratio(paired_error, norm))
                    norms[index] += np.array([norm, analytic_error, paired_error])**2
                    counts[index] += 1
                rows.append(row)
        for index, fraction in enumerate(settings["modes"]):
            norm, analytic_error, paired_error = np.sqrt(norms[index])
            summaries.append({"position": position, "radius_pixels": radius, "angle_degrees": angle,
                              "mode_fraction": fraction, "compared_fits": int(counts[index]),
                              "total_fits": len(baseline_paths), "analytic_relative_error": ratio(analytic_error, norm),
                              "paired_relative_error": ratio(paired_error, norm)})
        print(f"Compared position {position}: {len(baseline_paths)} detector fits", flush=True)
    write_csv(output / "detectors.csv", rows)
    write_csv(output / "summary.csv", summaries)
    (output / "summary.json").write_text(json.dumps(summaries, indent=2, allow_nan=False) + "\n")
    (output / "complete.json").write_text(json.dumps({"detector_mode_comparisons": len(rows)}) + "\n")


if __name__ == "__main__":
    main()
