"""Probe complex nseries admission, old behavior, and whole-series routing.

The supervisor takes heavy.lock nonblockingly and never replays a corpus.
If another validation owns the lock, it exits 3 without importing engines.
"""

import argparse
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[3]
BASELINE = Path("/tmp/symox-wave47-baseline/python")
CASES = [
    ("irrational_log_I", "log(x**sqrt(13)*(1 + x**r))", "I"),
    ("irrational_log_minus_I", "log(x**sqrt(13)*(1 + x**r))", "-I"),
    ("irrational_log_default", "log(x**sqrt(13)*(1 + x**r))", "0"),
    ("irrational_log_positive", "log(x**sqrt(13)*(1 + x**r))", "1"),
    ("legacy_integer_log_I", "log(x**4*(1 + x))", "I"),
    ("legacy_abs_I", "Abs(x)", "I"),
    ("legacy_abs_log_product_I", "Abs(x)*log(1 + x)", "I"),
    ("new_nonlog_inverse_I", "1/(1 + x**r)", "I"),
]


def snapshot(module, operation):
    try:
        result = operation()
        return {"status": "value", "text": str(result), "srepr": module.srepr(result)}
    except Exception as error:
        return {"status": "error", "type": type(error).__name__, "message": str(error)}


def worker(engine):
    if engine == "sympy":
        import sympy as module
    else:
        import symox as module
    x = module.Symbol("x")
    namespace = dict(vars(module), x=x, r=module.log(3)/module.log(2))
    rows = {}
    for name, source, direction in CASES:
        expr, cdir = eval(source, namespace), eval(direction, namespace)
        public = snapshot(module, lambda: expr.nseries(x, n=3, cdir=cdir))
        if engine == "sympy":
            front = snapshot(module, lambda: expr.series(x, 0, 3, cdir=cdir))
        else:
            front = snapshot(module, lambda: module._core._series(expr, x, 0, 3, "+", None, cdir))
        rows[name] = {"source": source, "cdir": direction, "public_nseries": public,
                      "whole_series_direction_front": front}
    print(json.dumps({"engine": engine, "module_file": module.__file__, "cases": rows}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, nargs="?")
    parser.add_argument("--worker", choices=("sympy", "baseline", "current"))
    parser.add_argument("--lock-held", action="store_true",
                        help="the supervising validation process already holds heavy.lock")
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
        return 0
    if args.output is None:
        parser.error("output JSON path is required")
    with Path("/home/keckj/symox-ci/heavy.lock").open("r") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            if not args.lock_held:
                print("heavy.lock is busy; no probe ran", file=sys.stderr)
                return 3
        else:
            if args.lock_held:
                print("--lock-held requires the supervising process to own heavy.lock", file=sys.stderr)
                return 3
        results = {}
        for engine in ("sympy", "baseline", "current"):
            env = dict(os.environ)
            env["PYTHONPATH"] = str(BASELINE if engine == "baseline" else ROOT / "python")
            command = [sys.executable, str(Path(__file__).resolve()), "--worker", engine]
            completed = subprocess.run(command, cwd=ROOT, env=env, text=True,
                                       capture_output=True, timeout=30, check=True)
            results[engine] = json.loads(completed.stdout)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, indent=2) + "\n")
        print(json.dumps({"engines": list(results), "cases_per_engine": len(CASES),
                          "output": str(args.output)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
