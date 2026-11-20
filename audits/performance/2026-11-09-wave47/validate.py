"""Replay wave-47 checks sequentially under the shared CI lock.

Usage: validate.py identity|focused|lots|suite|rust|corpora|recent|boundaries|ode-boundaries
Build/copy the release extension successfully before running checks.
Commands, exit statuses, timings and output are retained in raw/.
"""
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

AUDIT = Path(__file__).resolve().parent
ROOT = AUDIT.parents[2]
RAW = AUDIT / "raw"
PYTHON = "/home/keckj/symox-ci/venv/bin/python"
CARGO = "/home/keckj/.cargo/bin/cargo"
TARGET = "/tmp/symox-wave46-target"
EVIDENCE = Path("/tmp/symox-wave47-evidence")
ENV = dict(os.environ, PYTHONPATH=str(ROOT / "python"),
           PATH="/home/keckj/symox-ci/venv/bin:/home/keckj/.cargo/bin:" + os.environ.get("PATH", ""),
           CARGO_TARGET_DIR=TARGET, PYO3_PYTHON=PYTHON,
           PYTEST_DISABLE_PLUGIN_AUTOLOAD="1", OPENBLAS_NUM_THREADS="1",
           OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
           XLA_PYTHON_CLIENT_PREALLOCATE="false", MPLBACKEND="Agg")


def identity():
    paths = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard",
         "crates", "python", "tests", "Cargo.toml", "Cargo.lock"], cwd=ROOT
    ).decode().splitlines()
    files = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
             for p in sorted(set(paths)) if (ROOT / p).is_file()}
    return {
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip(),
        "source_sha256": hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest(),
        "so_sha256": hashlib.sha256((ROOT / "python/symox/_core.abi3.so").read_bytes()).hexdigest(),
        "files": files,
    }


def run(name, command):
    print("START", name, flush=True)
    start = time.time()
    with (RAW / (name + ".log")).open("w") as output:
        result = subprocess.run(command, cwd=ROOT, env=ENV, stdout=output, stderr=subprocess.STDOUT)
    record = {"command": list(map(str, command)), "returncode": result.returncode,
              "elapsed_seconds": round(time.time() - start, 3),
              "environment": {k: ENV[k] for k in ("PYTHONPATH", "PATH", "XLA_PYTHON_CLIENT_PREALLOCATE")}}
    extension = Path(ENV["PYTHONPATH"].split(os.pathsep)[0]) / "symox/_core.abi3.so"
    if extension.is_file():
        record["extension_sha256"] = hashlib.sha256(extension.read_bytes()).hexdigest()
    (RAW / (name + ".command.json")).write_text(json.dumps(record, indent=2) + "\n")
    print("END", name, json.dumps(record), flush=True)
    if result.returncode:
        print((RAW / (name + ".log")).read_text()[-12000:], flush=True)
        raise SystemExit(result.returncode)


def main():
    RAW.mkdir(exist_ok=True)
    EVIDENCE.mkdir(exist_ok=True)
    mode = sys.argv[1]
    if mode == "identity":
        record = identity()
        record["previous_so_sha256"] = record.pop("so_sha256")
        (RAW / "build-input.json").write_text(json.dumps(record, indent=2) + "\n")
        return
    with open("/home/keckj/symox-ci/heavy.lock", "rb") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        record = identity()
        expected = json.loads((RAW / "build-input.json").read_text())
        if record["source_sha256"] != expected["source_sha256"]:
            raise SystemExit("Source tree changed since the recorded build input; rebuild and record it again.")
        (RAW / (mode + "-identity.json")).write_text(json.dumps(record, indent=2) + "\n")
        if mode == "focused":
            run("focused", [PYTHON, "-m", "pytest", "-q", "--tb=short",
                "-o", "cache_dir=/tmp/symox-wave47-pytest-cache",
                "tests/test_series_wave47.py", "tests/test_native_matcher_wave47.py",
                "tests/test_polytools_wave47.py", "tests/test_limits_wave46.py", "tests/test_native_matcher_wave46.py",
                "tests/test_polytools_wave46.py", "tests/test_native_matcher_parity.py",
                "tests/test_algfactor_parity.py"])
        elif mode == "lots":
            run("matcher47-probe", [PYTHON,
                "audits/performance/2026-11-09-matcher47/probe_baseline.py",
                str(RAW / "matcher47-after.json")])
            run("polytools47-probe", [PYTHON,
                "audits/performance/2026-11-09-polytools47/probe.py",
                str(RAW / "polytools47-after.json")])
            run("series47-boundary", [PYTHON,
                "audits/performance/2026-11-09-series47/probe_certified_boundary.py",
                str(RAW / "series47-boundary.json")])
            run("series47-cdir", [PYTHON,
                "audits/performance/2026-11-09-series47/probe_cdir_admission.py",
                str(RAW / "series47-cdir-after.json"), "--lock-held"])
            matcher = json.loads((RAW / "matcher47-after.json").read_text())
            polytools = json.loads((RAW / "polytools47-after.json").read_text())
            if not all(row["public_agrees"] for row in matcher["rows"]):
                raise SystemExit("New matcher witness disagrees with reference; inspect rows.")
            for row in matcher["rows"]:
                reference, native = row["sympy"], row["native"]
                if reference["ok"]:
                    correct = native.get("handled") and native.get("value") == reference["value"]
                else:
                    correct = (native.get("type"), native.get("message")) == (reference["type"], reference["message"])
                if not correct:
                    raise SystemExit("Native matcher witness mismatch: " + str(row))
            if not all(row["reference"] == row["actual"] for row in polytools):
                raise SystemExit("New cancel witness disagrees with reference; inspect rows.")
            cdir = json.loads((RAW / "series47-cdir-after.json").read_text())
            for name, row in cdir["current"]["cases"].items():
                expected_engine = "sympy" if name in {"irrational_log_default", "irrational_log_positive"} else "baseline"
                expected = cdir[expected_engine]["cases"][name]["public_nseries"]
                if row["public_nseries"] != expected:
                    raise SystemExit("Series direction admission mismatch: " + name)
        elif mode == "suite":
            run("suite-final", [PYTHON, "audits/performance/2026-11-09-wave46/run_suite.py",
                str(ROOT), str(EVIDENCE / "suite-final")])
        elif mode == "rust":
            run("cargo-test", [CARGO, "test", "--offline", "--release"])
            run("clippy", [CARGO, "clippy", "--offline", "--release", "--",
                           "-A", "clippy::all", "-D", "clippy::correctness"])
        elif mode == "corpora":
            run("ode", ["sh", "audits/performance/2026-11-05-ode/scripts/run_corpus.sh",
                        str(ROOT), str(EVIDENCE / "ode")])
            run("ode-vs-wave46", [PYTHON, "audits/performance/2026-11-05-ode/scripts/score.py",
                str(EVIDENCE / "ode"), "/tmp/symox-wave47-baseline/ode"])
            solve = "audits/performance/2026-11-06-tsolve/scripts/"
            for corpus in ("tsolve", "unrad"):
                result = str(RAW / (corpus + "-symox.json"))
                run(corpus, [PYTHON, solve + "run_cases.py", "symox", result,
                             "--corpus", corpus, "--trace"])
                run(corpus + "-compare", [PYTHON, solve + "compare.py",
                    "audits/performance/2026-11-06-" + corpus + "/raw/sympy.json", result])
            matcher = "audits/performance/2026-11-06-matcher/"
            run("matcher-sympy", [PYTHON, matcher + "scripts/compare.py", "sympy",
                matcher + "corpus.json.gz", str(RAW / "matcher-sympy.json")])
            run("matcher", [PYTHON, matcher + "scripts/compare.py", "symox",
                            matcher + "corpus.json.gz", str(RAW / "matcher-after.json")])
            run("matcher-report", [PYTHON, matcher + "scripts/compare.py", "report",
                str(RAW / "matcher-sympy.json"), str(RAW / "matcher-after.json"),
                matcher + "corpus.json.gz", str(RAW / "matcher-report-after.json")])
            for family, date in (("algfactor", "2026-11-07"), ("algpoly", "2026-11-06"),
                                 ("ratint2", "2026-11-06")):
                folder = "audits/performance/" + date + "-" + family + "/"
                result = str(RAW / (family + "-symox.jsonl"))
                run(family, [PYTHON, folder + "diffcorpus.py", "symox", result])
                run(family + "-compare", [PYTHON, folder + "compare.py", folder + "raw/sympy.jsonl", result])
        elif mode == "ode-boundaries":
            run("ode-paired-boundaries", [PYTHON, str(AUDIT / "audit_ode_boundaries.py"),
                "--out", str(RAW / "ode-paired-boundaries"), "--lock-held", "--jobs", "2"])
        elif mode == "boundaries":
            before_rows = json.loads((ROOT / "audits/performance/2026-11-09-wave46/raw/liegroup-symox.json").read_text())
            after_rows = json.loads((RAW / "liegroup-symox.json").read_text())
            before = {row["id"]: row["result"] for row in before_rows}
            after = {row["id"]: row["result"] for row in after_rows}
            assert len(before) == len(before_rows) == 428 and before.keys() == after.keys()
            changes = {}
            for key, old in before.items():
                new = after[key]
                old = {k: v for k, v in old.items() if k != "seconds"}
                new = {k: v for k, v in new.items() if k != "seconds"}
                if old != new and any(row.get("error") in {"DEADLINE", "CRASH", "MemoryError"}
                                      for row in (old, new)):
                    changes[key] = {"before": old, "after": new}
            (RAW / "lie-boundary-selection.json").write_text(json.dumps(changes, indent=2) + "\n")
            ENV["PYTHONHASHSEED"] = "0"
            for label, python_path in (
                ("baseline", "/tmp/symox-wave47-baseline/python"),
                ("final", str(ROOT / "python")),
            ):
                ENV["PYTHONPATH"] = python_path
                if changes:
                    run("lie-boundary-" + label, [PYTHON, str(AUDIT / "boundary_cases.py"),
                        str(RAW / ("lie-boundary-" + label + ".json")), *sorted(changes)])
            (RAW / "lie-boundary-environment.json").write_text(json.dumps({
                "PYTHONHASHSEED": "0", "cases": sorted(changes), "case_limit_seconds": 120,
                "isolation_grace_seconds": 10,
            }, indent=2) + "\n")
        elif mode == "recent":
            baseline = ROOT / "audits/performance/2026-11-09-wave46/raw"
            for family, date in (("pde", "2026-11-07"), ("liegroup", "2026-11-06")):
                folder = "audits/performance/" + date + "-" + family + "/"
                reference = RAW / (family + "-sympy.json")
                reference.write_bytes(gzip.decompress((ROOT / folder / "sympy.json.gz").read_bytes()))
                result = RAW / (family + "-symox.json")
                report = RAW / (family + "-report.json")
                run(family, [PYTHON, folder + "scripts/corpus.py", "symox", str(result)])
                run(family + "-report", [PYTHON, folder + "scripts/corpus.py", "report",
                                         str(reference), str(result), str(report)])
                before = json.loads((baseline / (family + "-report.json")).read_text())
                after = json.loads(report.read_text())
                changes = {k: [v["verdict"], after["cases"][k]["verdict"]]
                           for k, v in before["cases"].items()
                           if v["verdict"] != after["cases"][k]["verdict"]}
                (RAW / (family + "-changes.json")).write_text(json.dumps(changes, indent=2) + "\n")
                print(family, "category changes", changes, flush=True)
            folder = "audits/performance/2026-11-06-hasmatch/"
            run("has", [PYTHON, folder + "scripts/has_corpus.py", "symox",
                        folder + "corpus.json.gz", str(RAW / "has-symox.json")])
            run("has-report", [PYTHON, folder + "scripts/has_corpus.py", "report",
                folder + "sympy-has.json", str(RAW / "has-symox.json"), str(baseline / "has-symox.json")])
            fps = "audits/performance/2026-11-07-fps/"
            sympy_root = Path(importlib.util.find_spec("sympy").origin).parent
            run("fps-formal", [PYTHON, fps + "scripts/tl.py", "symox", str(RAW / "fps-formal.json"),
                               str(sympy_root / "series/tests/test_formal.py")])
            run("fps-extra", [PYTHON, fps + "scripts/extra.py", "symox", str(RAW / "fps-extra.json")])
            for suffix in ("formal", "extra"):
                run("fps-" + suffix + "-compare", [PYTHON, fps + "scripts/compare.py",
                    str(ROOT / fps / "raw" / ("sympy_" + suffix)),
                    str(baseline / ("fps-" + suffix)), str(RAW / ("fps-" + suffix))])
            limits = "audits/performance/2026-11-06-limits2/"
            ENV["WORK"] = "/tmp/symox-wave47-limits"
            run("limits2", ["bash", limits + "scripts/run_tl.sh", "wave47", "symox", str(ROOT)])
            destination = RAW / "limits2"
            destination.mkdir(exist_ok=True)
            for result in (ROOT / limits / "raw/tl").glob("wave47_*"):
                result.rename(destination / result.name)
            run("limits2-compare", [PYTHON, limits + "scripts/compare.py",
                str(baseline / "limits2/wave46"), str(destination / "wave47")])
        else:
            raise SystemExit("unknown phase: " + mode)


if __name__ == "__main__":
    main()
