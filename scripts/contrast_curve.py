#!/usr/bin/env python3
## \file contrast_curve.py
# \brief Run KLIP or P4 injection grids and estimate contrast from recovered SNR.

"""Run isolated planet injections and estimate contrast from mean recovered SNR.

Each separation/PA pair gets one reduction and one hciAnalyze invocation. Mode
planes are averaged independently before scaling to the requested SNR. See
doc/contrast_curve.dox for examples and the statistical interpretation.
"""

from __future__ import annotations

import argparse
import bisect
import configparser
import csv
from decimal import Decimal
import fcntl
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import shutil
import statistics
import subprocess
import sys


## Result-table header used to recognize hciAnalyze measurements.
SNR_HEADER = "plane mode signal separation position_angle contrast snr"

## Reduction options controlled by the script to preserve each trial's identity.
REDUCTION_OWNED = {
    "config", "mode", "preProcess.only", "fake.sep", "fake.PA", "fake.contrast",
    "output.directory", "output.fileName", "output.exactFName",
    "p4Batch.file", "p4Batch.sourceSize", "p4Batch.keepGoing",
}
## Analysis options controlled by the script to measure the requested injection.
ANALYSIS_OWNED = {
    "config", "file", "noise.only", "planet.sep", "planet.PA", "planet.contrast",
    "fake.sep", "fake.PA", "fake.contrast",
}


def read_config(path: Path) -> tuple[str, dict[str, str]]:
    """Read an mxlib INI file, including its sectionless application options."""
    text = path.read_text()
    parser = configparser.ConfigParser(interpolation=None, strict=False, inline_comment_prefixes=("#", ";"))
    parser.optionxform = str
    parser.read_string("[__global__]\n" + text)
    values = {}
    for section in parser.sections():
        for key, value in parser.items(section):
            values[key if section == "__global__" else f"{section}.{key}"] = value
    return text, values


def overrides(items: list[str], owned: set[str]) -> dict[str, str]:
    """Parse shell-free KEY=VALUE application overrides, protecting trial identity."""
    result = {}
    for item in items:
        key, separator, value = item.removeprefix("--").partition("=")
        if not separator or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.]*", key) or key in owned:
            raise ValueError(f"invalid or script-controlled override: {item}")
        result[key] = value
    return result


def samples(minimum: float, maximum: float, step: float) -> list[float]:
    """Generate an inclusive grid without accumulated floating-point step error."""
    if not all(math.isfinite(x) for x in (minimum, maximum, step)) or step <= 0 or maximum < minimum:
        raise ValueError("sample bounds must be finite and ordered, with a positive step")
    low, high, delta = (Decimal(str(x)) for x in (minimum, maximum, step))
    count = int((high - low) // delta) + 1
    if count > 100000:
        raise ValueError("sample grid exceeds 100000 points")
    return [float(low + i * delta) for i in range(count)]


def contrasts(args: argparse.Namespace, radii: list[float]) -> list[float]:
    """Evaluate a constant contrast or linearly interpolate a separation/contrast CSV."""
    if args.contrast_table:
        with args.contrast_table.open(newline="") as stream:
            reader = csv.DictReader(stream)
            missing = [name for name in ("separation", "contrast") if name not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(
                    f"contrast table {args.contrast_table}: missing column(s) {', '.join(missing)}; "
                    f"include the comma-separated header 'separation,contrast' before the data rows "
                    f"(found {reader.fieldnames or []!r})"
                )
            points = [(float(row["separation"]), float(row["contrast"])) for row in reader]
        if not points or any(not math.isfinite(r) or not math.isfinite(c) or r < 0 or c <= 0 for r, c in points):
            raise ValueError("contrast table needs finite nonnegative separations and positive contrasts")
        points.sort()
        knots = [point[0] for point in points]
        if len(set(knots)) != len(knots) or radii[0] < knots[0] or radii[-1] > knots[-1]:
            raise ValueError("contrast table must have unique separations and cover the entire radius grid")
        result = []
        for radius in radii:
            upper = bisect.bisect_left(knots, radius)
            if knots[upper] == radius:
                result.append(points[upper][1])
            else:
                r0, c0 = points[upper - 1]
                r1, c1 = points[upper]
                result.append(c0 + (c1 - c0) * (radius - r0) / (r1 - r0))
        return result
    if args.contrast is None or not math.isfinite(args.contrast) or args.contrast <= 0:
        raise ValueError("provide --contrast with a positive value, or --contrast-table")
    return [args.contrast] * len(radii)


def executable(name: str) -> dict:
    """Resolve an executable and record its identity for safe resumption."""
    found = shutil.which(name)
    if not found:
        raise ValueError(f"executable not found: {name}")
    path = Path(found).resolve()
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def cli(values: dict) -> list[str]:
    """Serialize application options as individual arguments, without a shell."""
    return [f"--{key}={value}" for key, value in values.items()]


def make_plan(args: argparse.Namespace) -> dict:
    """Validate settings and construct a reproducible trial manifest."""
    config_path = args.config.resolve()
    reduction_text, reduction = read_config(config_path)
    selected = [name for name, key in (("p4Reduce", "p4.modeFractions"), ("klipReduce", "klip.Nmodes")) if key in reduction]
    reducer = args.reducer or (selected[0] if len(selected) == 1 else None)
    if reducer is None:
        raise ValueError("cannot infer reducer from config; specify --reducer")
    if args.p4_batch and reducer != "p4Reduce":
        raise ValueError("--p4-batch requires p4Reduce")
    if args.p4_source_size <= 0 or (not args.p4_batch and args.p4_source_size != 10):
        raise ValueError("--p4-source-size must be positive and requires --p4-batch")
    reduce_args = overrides(args.reduce_arg, REDUCTION_OWNED)
    reduction.update(reduce_args)
    if reduction.get("p4Batch.file"):
        raise ValueError("the script owns the injection list; use --p4-batch instead of p4Batch.file")
    if not reduction.get("fake.fileName"):
        raise ValueError("fake.fileName must specify the injection PSF in the config or --reduce-arg")
    analysis_text, analysis = read_config(args.analyze_config.resolve()) if args.analyze_config else ("", {})
    analyze_args = overrides(args.analyze_arg, ANALYSIS_OWNED)
    analysis.update(analyze_args)
    has_fake = any(key.startswith("fake.") for key in analysis)
    has_planet = any(key.startswith("planet.") for key in analysis)
    if has_fake and has_planet:
        raise ValueError("analysis config must not mix planet and fake groups")
    for key, value in {"lambdaD": "4.8", "filter.lpfGaussFW": "4.8" if reducer == "p4Reduce" else "2"}.items():
        if key not in analysis:
            analyze_args[key] = value
    radii = samples(args.radius_min, args.radius_max, args.radius_step)
    if radii[0] <= 0:
        raise ValueError("injection separations must be positive")
    # 0 and 360 degrees are the same trial; retain the first occurrence only.
    pas = list(dict.fromkeys(round(pa % 360, 10) for pa in samples(args.pa_min, args.pa_max, args.pa_step)))
    if len(radii) * len(pas) > 100000:
        raise ValueError("injection grid exceeds 100000 trials")
    if not math.isfinite(args.target_snr) or args.target_snr <= 0:
        raise ValueError("target SNR must be finite and positive")
    if not math.isfinite(args.min_pa_fraction) or not 0 < args.min_pa_fraction <= 1:
        raise ValueError("minimum PA fraction must lie in (0,1]")
    return {
        "schema": 1, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reducer": reducer, "reducer_binary": executable(args.reducer_bin or reducer),
        "analyzer_binary": executable(args.analyzer_bin),
        "reduction_config": {"path": str(config_path), "text": reduction_text},
        "analysis_config": {"path": str(args.analyze_config.resolve()) if args.analyze_config else None,
                            "text": analysis_text},
        "reduction_cwd": str(args.work_dir.resolve() if args.work_dir else config_path.parent),
        "analysis_cwd": str(args.analyze_config.resolve().parent if args.analyze_config else
                            (args.work_dir.resolve() if args.work_dir else config_path.parent)),
        "reduce_args": reduce_args, "analyze_args": analyze_args,
        "signal_group": "fake" if has_fake else "planet", "radii": radii, "pas": pas,
        "contrasts": contrasts(args, radii), "target_snr": args.target_snr,
        "min_pa_fraction": args.min_pa_fraction,
        "p4_batch": args.p4_batch, "p4_source_size": args.p4_source_size if args.p4_batch else None,
    }


def commands(plan: dict, root: Path, trial: dict) -> tuple[list[str], list[str]]:
    """Construct reduction and analysis commands for one isolated injection."""
    directory = root / trial["name"]
    reduce_command = [plan["reducer_binary"]["path"], f"--config={root / 'reduction.conf'}"]
    reduce_command += cli(plan["reduce_args"])
    reduce_command += cli({"mode": "normal", "preProcess.only": "false", "fake.sep": trial["separation"],
                           "fake.PA": trial["pa"], "fake.contrast": trial["contrast"],
                           "output.directory": directory, "output.fileName": "finim.fits", "output.exactFName": "true"})
    analyze_command = [plan["analyzer_binary"]["path"], f"--config={root / 'analysis.conf'}"]
    analyze_command += cli(plan["analyze_args"])
    group = plan["signal_group"]
    analyze_command += cli({"file": directory / "finim.fits", "noise.only": "false",
                            f"{group}.sep": trial["separation"], f"{group}.PA": trial["pa"],
                            f"{group}.contrast": trial["contrast"]})
    return reduce_command, analyze_command


## \ingroup programming_library
def batch_command(plan: dict, root: Path, keep_going: bool = False) -> list[str]:
    """Build the persistent P4 command, whose configured output is the uninjected baseline."""
    return [plan["reducer_binary"]["path"], f"--config={root / 'reduction.conf'}"] + cli(plan["reduce_args"]) + cli({
        "mode": "normal", "preProcess.only": "false", "p4Optimize.enabled": "false",
        "fake.sep": plan["radii"][0], "fake.PA": plan["pas"][0], "fake.contrast": 0,
        "output.directory": root, "output.fileName": "baseline.fits", "output.exactFName": "true",
        "p4Batch.file": root / "batch_trials.txt", "p4Batch.sourceSize": plan["p4_source_size"],
        "p4Batch.keepGoing": str(keep_going).lower(),
    })


## \ingroup programming_library
def batch_row(root: Path, trial: dict) -> str:
    """Serialize a single injection using the C++ batch reader's quoted-path convention."""
    path = str(root / trial["name"] / "finim.fits")
    if "\n" in path or "\r" in path:
        raise ValueError("P4 batch output paths cannot contain newline characters")
    quoted = '"' + path.replace('\\', '\\\\').replace('"', '\\"') + '"'
    return f"{trial['separation']!r} {trial['pa']!r} {trial['contrast']!r} {quoted}"


## \ingroup programming_library
def run_batch(plan: dict, root: Path, trials: list[dict], resume: bool, keep_going: bool) -> str | None:
    """Reduce pending injections in one process and recover individually published products after failure."""
    pending = []
    for trial in trials:
        directory = root / trial["name"]
        directory.mkdir(exist_ok=True)
        completed = directory / "result.json"
        if resume and completed.exists() and json.loads(completed.read_text())["status"] == "ok":
            continue
        if resume and (directory / "finim.fits").is_file() and (
                (directory / "reduction.complete").exists() or (directory / "finim.fits.complete").exists()):
            (directory / "reduction.complete").touch()
            continue
        pending.append(trial)
    if not pending:
        return None
    rows = [batch_row(root, trial) for trial in pending]
    (root / "batch_trials.txt").write_text("\n".join(rows) + "\n")
    command = batch_command(plan, root, keep_going)
    write_json(root / "batch_commands.json", {"reduction": command, "trials": [trial["name"] for trial in pending]})
    print(f"P4 batch: one baseline plus {len(pending)} injections; progress: {root / 'batch.stderr'}", flush=True)
    error = None
    try:
        run_command(command, plan["reduction_cwd"], root / "batch.stdout")
    except (OSError, RuntimeError) as exception:
        error = str(exception)
    for trial in pending:
        directory = root / trial["name"]
        if (directory / "finim.fits").is_file() and (directory / "finim.fits.complete").is_file():
            (directory / "reduction.complete").touch()
    return error


def parse_snr(text: str, trial: dict) -> list[dict]:
    """Read hciAnalyze's result table and verify it measures the requested trial."""
    rows = []
    in_table = False
    for line in text.splitlines():
        if line.strip() == SNR_HEADER:
            if in_table:
                raise ValueError("duplicate hciAnalyze result header")
            in_table = True
            continue
        if not in_table or not line.strip():
            continue
        fields = line.split()
        if len(fields) != 7:
            raise ValueError(f"malformed hciAnalyze result: {line}")
        plane, mode, signal, separation, pa, contrast, snr = fields
        if int(signal) != 0 or int(plane) < 0:
            raise ValueError("expected one injected signal per nonnegative mode plane")
        for value, expected in ((float(separation), trial["separation"]), (float(pa), trial["pa"]),
                                (float(contrast), trial["contrast"])):
            if not math.isclose(value, expected, rel_tol=1e-5, abs_tol=0):
                raise ValueError("hciAnalyze measured a different injection than requested")
        if not math.isfinite(float(mode)):
            raise ValueError("hciAnalyze output lacks finite mode labels")
        rows.append({"plane": int(plane), "mode": float(mode),
                     "snr": float(snr) if math.isfinite(float(snr)) else None})
    if not rows or sorted(row["plane"] for row in rows) != list(range(len(rows))):
        raise ValueError("hciAnalyze returned missing, duplicate, or non-contiguous mode planes")
    return sorted(rows, key=lambda row: row["plane"])


def write_json(path: Path, value: dict) -> None:
    """Atomically replace a manifest or completion record."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def run_command(command: list[str], cwd: str, log: Path) -> None:
    """Run one application and retain separate stdout and stderr logs."""
    with log.open("w") as stdout, log.with_suffix(".stderr").open("w") as stderr:
        result = subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr, check=False)
    if result.returncode:
        raise RuntimeError(f"{Path(command[0]).name} exited {result.returncode}; see {log.with_suffix('.stderr')}")


def run_trial(plan: dict, root: Path, trial: dict, resume: bool, batch_error: str | None = None) -> dict:
    """Resume completed analysis, or run the missing stages of one trial."""
    directory = root / trial["name"]
    directory.mkdir(exist_ok=True)
    completed = directory / "result.json"
    if resume and completed.exists():
        result = json.loads(completed.read_text())
        if result["status"] == "ok":
            return result
    reduce_command, analyze_command = commands(plan, root, trial)
    write_json(directory / "commands.json", {
        "reduction": str(root / "batch_commands.json") if plan.get("p4_batch") else reduce_command,
        "analysis": analyze_command,
        **({"batch_entry": batch_row(root, trial)} if plan.get("p4_batch") else {}),
    })
    result = {**trial, "status": "failed", "measurements": []}
    try:
        marker = directory / "reduction.complete"
        if plan.get("p4_batch"):
            if not (marker.exists() and (directory / "finim.fits").is_file()):
                raise RuntimeError(batch_error or f"P4 batch did not complete {trial['name']}; see {root / 'batch.stderr'}")
        elif not (resume and marker.exists() and (directory / "finim.fits").is_file()):
            marker.unlink(missing_ok=True)
            run_command(reduce_command, plan["reduction_cwd"], directory / "reduction.stdout")
            if not (directory / "finim.fits").is_file():
                raise RuntimeError("reducer did not write finim.fits")
            marker.touch()
        run_command(analyze_command, plan["analysis_cwd"], directory / "analysis.stdout")
        result["measurements"] = parse_snr((directory / "analysis.stdout").read_text(), trial)
        result["status"] = "ok"
    except (OSError, ValueError, RuntimeError) as error:
        result["error"] = str(error)
    write_json(completed, result)
    return result


def summarize(plan: dict, results: list[dict]) -> tuple[list[dict], list[dict]]:
    """Average signed finite SNRs by mode and radius before contrast rescaling."""
    labels = None
    for result in results:
        if result["status"] != "ok":
            continue
        current = [(row["plane"], row["mode"]) for row in result["measurements"]]
        if labels is not None and labels != current:
            raise ValueError("output mode planes changed between injections; refusing to mix reductions")
        labels = current
    curves, best = [], []
    expected = len(plan["pas"])
    required = math.ceil(plan["min_pa_fraction"] * expected)
    for radius, contrast in zip(plan["radii"], plan["contrasts"]):
        at_radius = [result for result in results if result["separation"] == radius]
        for plane, mode in labels or []:
            snrs = [result["measurements"][plane]["snr"] for result in at_radius
                    if result["status"] == "ok" and result["measurements"][plane]["snr"] is not None]
            mean = statistics.fmean(snrs) if snrs else math.nan
            status = "insufficient_pa" if len(snrs) < required else "nonpositive_snr" if mean <= 0 else "ok"
            curves.append({"separation": radius, "plane": plane, "mode": mode, "injected_contrast": contrast,
                           "pa_expected": expected, "pa_used": len(snrs), "mean_snr": mean,
                           "snr_stddev": statistics.stdev(snrs) if len(snrs) > 1 else math.nan,
                           "target_snr": plan["target_snr"],
                           "contrast": contrast * plan["target_snr"] / mean if status == "ok" else math.nan,
                           "status": status})
        valid = [row for row in curves if row["separation"] == radius and row["status"] == "ok"]
        if valid:
            best.append(min(valid, key=lambda row: row["contrast"]))
    return curves, best


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    """Write an explicitly headed CSV, including empty result sets."""
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(plan: dict, root: Path, results: list[dict], plot: bool) -> None:
    """Write trial measurements, per-mode curves, their best-mode envelope, and plots."""
    trials = []
    for result in results:
        for measurement in result["measurements"] or [{"plane": "", "mode": "", "snr": None}]:
            trials.append({key: result[key] for key in ("name", "separation", "pa", "contrast", "status")} |
                          measurement | {"error": result.get("error", "")})
    write_csv(root / "measurements.csv", ["name", "separation", "pa", "contrast", "status", "plane", "mode", "snr", "error"], trials)
    curves, best = summarize(plan, results)
    fields = ["separation", "plane", "mode", "injected_contrast", "pa_expected", "pa_used", "mean_snr",
              "snr_stddev", "target_snr", "contrast", "status"]
    write_csv(root / "contrast_by_mode.csv", fields, curves)
    write_csv(root / "contrast_best.csv", fields, best)
    if plot:
        import matplotlib.pyplot as plt

        figure, axes = plt.subplots(figsize=(9, 6), layout="constrained")
        for plane in sorted({row["plane"] for row in curves}):
            rows = [row for row in curves if row["plane"] == plane]
            label = "mode fraction" if plan["reducer"] == "p4Reduce" else "modes"
            axes.plot([row["separation"] for row in rows], [row["contrast"] for row in rows],
                      marker=".", alpha=0.6, linewidth=1, label=f"{label} {rows[0]['mode']:g}")
        envelope = {row["separation"]: row["contrast"] for row in best}
        axes.plot(plan["radii"], [envelope.get(radius, math.nan) for radius in plan["radii"]],
                  "ko-", linewidth=2, label="best mode at each separation")
        axes.set(xlabel="Separation [pixels]", ylabel=f"Estimated contrast at SNR = {plan['target_snr']:g}",
                 yscale="log", title=f"{plan['reducer']}: mean recovered SNR across PA")
        axes.grid(True, which="both", alpha=0.2)
        axes.legend(fontsize="small", ncols=2)
        for extension in ("png", "pdf"):
            figure.savefig(root / f"contrast_curve.{extension}", dpi=160)
        plt.close(figure)


def argument_parser() -> argparse.ArgumentParser:
    """Define the injection grid, application overrides, and execution controls."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="klipReduce or p4Reduce configuration")
    parser.add_argument("--reducer", choices=("klipReduce", "p4Reduce"), help="default: infer from configuration")
    parser.add_argument("--reducer-bin", help="path/name of the selected reduction executable")
    parser.add_argument("--p4-batch", action="store_true", help="reuse one P4 baseline and refit direct/OR-affected targets")
    parser.add_argument("--p4-source-size", type=int, default=10, help="batch source crop width; default 10, parity preserved")
    parser.add_argument("--analyzer-bin", default="hciAnalyze")
    parser.add_argument("--analyze-config", type=Path)
    parser.add_argument("--reduce-arg", action="append", default=[], metavar="KEY=VALUE")
    parser.add_argument("--analyze-arg", action="append", default=[], metavar="KEY=VALUE")
    parser.add_argument("--work-dir", type=Path, help="reduction cwd; default: reduction config directory")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--radius-min", type=float, required=True)
    parser.add_argument("--radius-max", type=float, required=True)
    parser.add_argument("--radius-step", "--radius-delta", type=float, required=True)
    parser.add_argument("--pa-min", type=float, default=0)
    parser.add_argument("--pa-max", type=float, default=330)
    parser.add_argument("--pa-step", "--pa-delta", type=float, default=30)
    brightness = parser.add_mutually_exclusive_group(required=True)
    brightness.add_argument("--contrast", type=float, help="constant injected contrast")
    brightness.add_argument("--contrast-table", type=Path, help="CSV columns: separation,contrast; linear interpolation")
    parser.add_argument("--target-snr", type=float, default=5)
    parser.add_argument("--min-pa-fraction", type=float, default=1, help="required finite PA fraction per mode; default: all")
    parser.add_argument("--resume", action="store_true", help="reuse completed stages with an identical manifest")
    parser.add_argument("--keep-going", action="store_true", help="record failed trials and continue; return failure at end")
    parser.add_argument("--dry-run", action="store_true", help="print commands without creating files or launching applications")
    parser.add_argument("--no-plot", action="store_true", help="write CSVs without requiring matplotlib")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Execute the requested grid, preserving logs and resumable completion records."""
    args = argument_parser().parse_args(argv)
    try:
        plan = make_plan(args)
        root = args.output_dir.resolve()
        trials = [{"name": f"r{ri:04d}_pa{pi:04d}", "separation": radius, "pa": pa, "contrast": contrast}
                  for ri, (radius, contrast) in enumerate(zip(plan["radii"], plan["contrasts"]))
                  for pi, pa in enumerate(plan["pas"])]
        if args.dry_run:
            print(f"{plan['reducer']}: {len(trials)} reductions and analyses; output: {root}")
            if plan["p4_batch"]:
                print(shlex.join(batch_command(plan, root, args.keep_going)))
            for trial in trials:
                if plan["p4_batch"]:
                    print("# batch trial: " + batch_row(root, trial))
                for command in commands(plan, root, trial)[1 if plan["p4_batch"] else 0:]:
                    print(shlex.join(command))
            return 0
        if not args.no_plot:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot  # Check plotting availability before starting expensive reductions.
        root.mkdir(parents=True, exist_ok=True)
        with (root / ".lock").open("a") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            manifest = root / "manifest.json"
            if args.resume:
                if not manifest.is_file() or json.loads(manifest.read_text()) != plan:
                    raise ValueError("resume requires an identical config, executable identity, script, and trial plan")
            elif any(path.name != ".lock" for path in root.iterdir()):
                raise ValueError("output directory is not empty; choose a new directory or use --resume")
            else:
                write_json(manifest, plan)
                (root / "reduction.conf").write_text(plan["reduction_config"]["text"])
                (root / "analysis.conf").write_text(plan["analysis_config"]["text"])
            batch_error = run_batch(plan, root, trials, args.resume, args.keep_going) if plan["p4_batch"] else None
            results = []
            for index, trial in enumerate(trials, 1):
                print(f"[{index}/{len(trials)}] r={trial['separation']:g} PA={trial['pa']:g} contrast={trial['contrast']:g}", flush=True)
                result = run_trial(plan, root, trial, args.resume, batch_error)
                results.append(result)
                write_outputs(plan, root, results, plot=False)
                if result["status"] != "ok":
                    print(result["error"], file=sys.stderr)
                    if not args.keep_going:
                        return 1
            write_outputs(plan, root, results, plot=not args.no_plot)
            print(f"Wrote measurements and contrast curves to {root}")
            return int(any(result["status"] != "ok" for result in results))
    except (OSError, ValueError, KeyError, configparser.Error, ImportError) as error:
        print(f"contrast_curve: {error}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Interrupted; use --resume to reuse completed stages.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
