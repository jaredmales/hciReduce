#!/usr/bin/env python3
## \file contrast_curve_checks.py
# \brief Exercise the production contrast-curve script through Catch2-driven checks.

"""Behavioral checks invoked by the Catch2 contrast-curve integration tests."""

import csv
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile


SCRIPT = Path(__file__).resolve().parents[2] / "scripts/contrast_curve.py"
SPEC = importlib.util.spec_from_file_location("contrast_curve", SCRIPT)
CURVE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CURVE)


MOCK_APP = '''#!/usr/bin/env python3
import json
from pathlib import Path
import shlex
import sys
options = dict(argument[2:].split("=", 1) for argument in sys.argv[1:])
base = Path(__file__).parent
role = Path(__file__).name
with (base / "calls.jsonl").open("a") as log:
    log.write(json.dumps({"role": role, "options": options}) + "\\n")
if role == "hciAnalyze":
    if (base / "fail").exists() and float(options["planet.PA"]) == 180:
        sys.exit(9)
    original = json.loads(Path(options["file"]).read_text())
    for key in ("sep", "PA", "contrast"):
        assert options["planet." + key] == original["fake." + key]
    print("plane mode signal separation position_angle contrast snr")
    pa = float(options["planet.PA"])
    for plane, mode in enumerate((0.1, 0.5)):
        snr = ((-2, 6), (4, 8))[plane][int(pa == 180)]
        if (base / "missing").exists() and plane == 0 and pa == 180:
            snr = float("nan")
        print(plane, mode, 0, options["planet.sep"], options["planet.PA"], options["planet.contrast"], snr)
else:
    assert options["output.exactFName"] == "true"
    assert options["preProcess.only"] == "false"
    assert options["mode"] == "normal"
    output = Path(options["output.directory"]) / options["output.fileName"]
    output.write_text(json.dumps(options))
    if "p4Batch.file" in options:
        failed = False
        for line in Path(options["p4Batch.file"]).read_text().splitlines():
            sep, pa, contrast, path = shlex.split(line)
            if (base / "batch_fail").exists() and float(pa) == 180:
                failed = True
                if options["p4Batch.keepGoing"] == "false":
                    break
                continue
            values = dict(options, **{"fake.sep": sep, "fake.PA": pa, "fake.contrast": contrast})
            Path(path).write_text(json.dumps(values))
            Path(path + ".complete").touch()
        if failed:
            sys.exit(9)
'''


def fixture(root: Path) -> list[str]:
    """Create executable application doubles and a conflicting base injection config."""
    for name in ("p4Reduce", "hciAnalyze"):
        path = root / name
        path.write_text(MOCK_APP)
        path.chmod(0o755)
    config = root / "p4.conf"
    config.write_text("mode=normal\n[p4]\nmodeFractions=0.1,0.5\n[fake]\nfileName=psf.fits\n"
                      "sep=99\nPA=99\ncontrast=1\n[output]\nfileName=do-not-use.fits\n")
    return [sys.executable, str(SCRIPT), str(config), "--reducer-bin", str(root / "p4Reduce"),
            "--analyzer-bin", str(root / "hciAnalyze"), "--radius-min=10", "--radius-max=20",
            "--radius-step=5", "--pa-min=0", "--pa-max=360", "--pa-step=180", "--contrast=0.0001",
            "--output-dir", str(root / "output with spaces"), "--no-plot"]


def run(command: list[str], expected: int = 0) -> subprocess.CompletedProcess:
    """Run the real script and include its output in any assertion failure."""
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode == expected, result.stdout + result.stderr
    return result


def read_csv(path: Path) -> list[dict]:
    """Read a generated curve or measurement table."""
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def grid_and_parser() -> None:
    """Check sampling, contrast interpolation, defaults, and strict result parsing."""
    assert CURVE.samples(0.1, 0.3, 0.1) == [0.1, 0.2, 0.3]
    for bounds in ((1, 0, 1), (0, 1, 0), (0, math.inf, 1)):
        try:
            CURVE.samples(*bounds)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid grid accepted")
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        command = fixture(root)
        table = root / "contrast.csv"
        table.write_text("separation,contrast\n10,0.0001\n20,0.0003\n")
        command.remove("--contrast=0.0001")
        command += ["--contrast-table", str(table)]
        args = CURVE.argument_parser().parse_args(command[2:])
        plan = CURVE.make_plan(args)
        assert plan["pas"] == [0, 180]
        assert math.isclose(plan["contrasts"][1], 0.0002)
        assert plan["analyze_args"] == {"lambdaD": "4.8", "filter.lpfGaussFW": "4.8"}
        args.reducer = "klipReduce"
        assert CURVE.make_plan(args)["analyze_args"]["filter.lpfGaussFW"] == "2"
        analysis = root / "analysis.conf"
        analysis.write_text("lambdaD=3\n[filter]\nlpfGaussFW=1.2\n[fake]\nR=3\n")
        args.analyze_config = analysis
        args.analyze_arg = ["filter.lpfGaussFW=2.7"]
        plan = CURVE.make_plan(args)
        assert plan["signal_group"] == "fake"
        assert plan["analyze_args"] == {"filter.lpfGaussFW": "2.7"}
        dry = run(command + ["--dry-run"])
        assert "6 reductions" in dry.stdout
        assert not (root / "output with spaces").exists()
        run(command + ["--reduce-arg=output.fileName=wrong.fits"], 1)
        run(command + ["--radius-max=25"], 1)
        for content in ("10,0.0001\n20,0.0003\n", "radius,contrast\n10,0.0001\n", ""):
            table.write_text(content)
            for extra in ([], ["--dry-run"]):
                failure = run(command + extra, 1)
                assert str(table) in failure.stderr
                assert "header 'separation,contrast' before the data rows" in failure.stderr
                assert "Traceback" not in failure.stderr
                assert not (root / "calls.jsonl").exists()
                assert not (root / "output with spaces").exists()
    trial = {"separation": 10, "pa": 0, "contrast": 0.0001}
    parsed = CURVE.parse_snr(CURVE.SNR_HEADER + "\n0 0.1 0 10 0 0.0001 nan\n", trial)
    assert parsed[0]["snr"] is None
    for text in ("", CURVE.SNR_HEADER + "\n0 0.1 0 99 0 0.0001 5\n",
                 CURVE.SNR_HEADER + "\n1 0.1 0 10 0 0.0001 5\n"):
        try:
            CURVE.parse_snr(text, trial)
        except ValueError:
            pass
        else:
            raise AssertionError("malformed analysis accepted")


def reduction_workflow() -> None:
    """Check real CLI orchestration, signed-SNR averaging, and immutable resume."""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        command = fixture(root)
        run(command)
        output = root / "output with spaces"
        rows = read_csv(output / "contrast_by_mode.csv")
        assert len(rows) == 6
        for row in rows:
            assert row["pa_used"] == row["pa_expected"] == "2"
            mean = 2 if row["plane"] == "0" else 6
            assert float(row["mean_snr"]) == mean
            assert math.isclose(float(row["contrast"]), 0.0001 * 5 / mean)
        assert all(row["plane"] == "1" for row in read_csv(output / "contrast_best.csv"))
        assert len(read_csv(output / "measurements.csv")) == 12
        calls = (root / "calls.jsonl").read_text()
        assert len(calls.splitlines()) == 12
        run(command + ["--resume"])
        assert (root / "calls.jsonl").read_text() == calls
        run(command, 1)
        run(command + ["--resume", "--analyze-arg=lambdaD=3"], 1)
        assert (root / "calls.jsonl").read_text() == calls


def failure_and_coverage() -> None:
    """Check failure persistence, reduction reuse, incomplete PA coverage, and invalid SNRs."""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        command = fixture(root)
        (root / "fail").touch()
        run(command + ["--keep-going"], 1)
        rows = read_csv(root / "output with spaces/contrast_by_mode.csv")
        assert all(row["status"] == "insufficient_pa" and math.isnan(float(row["contrast"])) for row in rows)
        (root / "fail").unlink()
        run(command + ["--resume"])
        calls = [json.loads(line) for line in (root / "calls.jsonl").read_text().splitlines()]
        assert sum(call["role"] == "p4Reduce" for call in calls) == 6
        assert sum(call["role"] == "hciAnalyze" for call in calls) == 9
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        command = fixture(root)
        (root / "missing").touch()
        run(command)
        output = root / "output with spaces"
        measurements = read_csv(output / "measurements.csv")
        assert len(measurements) == 12 and all(row["status"] == "ok" for row in measurements)
        assert sum(row["snr"] == "" for row in measurements) == 3
        curves = read_csv(output / "contrast_by_mode.csv")
        assert all(row["status"] == "insufficient_pa" and row["pa_used"] == "1"
                   for row in curves if row["plane"] == "0")
        assert all(row["status"] == "ok" and row["pa_used"] == "2"
                   for row in curves if row["plane"] == "1")
        assert all(row["plane"] == "1" for row in read_csv(output / "contrast_best.csv"))
        before = (root / "calls.jsonl").read_text()
        run(command + ["--resume"])
        assert (root / "calls.jsonl").read_text() == before
    plan = {"pas": [0, 180], "radii": [10], "contrasts": [0.0001], "target_snr": 5, "min_pa_fraction": 0.5}
    result = {"status": "ok", "separation": 10, "measurements": [{"plane": 0, "mode": 0.1, "snr": None}]}
    curves, best = CURVE.summarize(plan, [result])
    assert curves[0]["status"] == "insufficient_pa" and not best
    result["measurements"][0]["snr"] = -2
    curves, best = CURVE.summarize(plan, [result])
    assert curves[0]["status"] == "nonpositive_snr" and not best
    result["measurements"][0]["snr"] = 2
    curves, best = CURVE.summarize(plan, [result])
    assert curves[0]["pa_used"] == 1 and best


def batch_workflow() -> None:
    """Check one-process P4 dispatch, failed-batch recovery, and analysis-only resumption."""
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        command = fixture(root) + ["--p4-batch", "--p4-source-size=10"]
        dry = run(command + ["--dry-run"])
        assert dry.stdout.count("--p4Batch.file=") == 1
        assert dry.stdout.count("# batch trial:") == 6
        assert not (root / "output with spaces").exists()
        (root / "batch_fail").touch()
        run(command + ["--keep-going"], 1)
        output = root / "output with spaces"
        calls = [json.loads(line) for line in (root / "calls.jsonl").read_text().splitlines()]
        assert sum(call["role"] == "p4Reduce" for call in calls) == 1
        assert sum(call["role"] == "hciAnalyze" for call in calls) == 3
        assert all(row["status"] == "insufficient_pa" for row in read_csv(output / "contrast_by_mode.csv"))
        (root / "batch_fail").unlink()
        (root / "fail").touch()
        run(command + ["--resume", "--keep-going"], 1)
        assert len((output / "batch_trials.txt").read_text().splitlines()) == 3
        (root / "fail").unlink()
        run(command + ["--resume"])
        calls = [json.loads(line) for line in (root / "calls.jsonl").read_text().splitlines()]
        assert sum(call["role"] == "p4Reduce" for call in calls) == 2
        assert sum(call["role"] == "hciAnalyze" for call in calls) == 9
        assert all(row["status"] == "ok" for row in read_csv(output / "contrast_by_mode.csv"))
        previous_calls = (root / "calls.jsonl").read_text()
        run(command + ["--resume"])
        assert (root / "calls.jsonl").read_text() == previous_calls
        run(command + ["--resume", "--p4-source-size=12"], 1)
        run(command + ["--reducer=klipReduce", "--dry-run"], 1)


if __name__ == "__main__":
    {"grid": grid_and_parser, "workflow": reduction_workflow, "failure": failure_and_coverage,
     "batch": batch_workflow}[sys.argv[1]]()
