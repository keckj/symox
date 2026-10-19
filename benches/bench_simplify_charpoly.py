#!/usr/bin/env python3
"""Focused simplify/charpoly timings, with results saved for exact comparison.

Run the same script in separate processes/environments before and after a build:
    python benches/bench_simplify_charpoly.py --lib symox --out before.json
    python benches/bench_simplify_charpoly.py --lib sympy --out sympy.json

Setup, warm-up and result printing are outside the timed region. Report both
the minimum and median of repeated calls; this measures warmed computations.
"""

import argparse
import importlib
import json
import statistics
import time
from pathlib import Path


def cases(sp):
    x, y, t = sp.symbols("x y t")
    identity = sp.sin(x)**4 - sp.cos(x)**4 + sp.cos(2*x)
    tangent = sp.tan(x)*sp.cos(x)
    rational = (x**2 - y**2)/(x - y)
    generic = sp.Matrix(6, 6, lambda i, j: sp.Symbol(f"a{i}{j}"))
    shared = sp.Matrix(4, 4, lambda i, j: x**((i + j) % 3) + y*(i - j) + (i + 1)*(j + 2) % 3)
    return {
        "simplify_identity": lambda: sp.simplify(identity),
        "simplify_tangent": lambda: sp.simplify(tangent),
        "simplify_rational": lambda: sp.simplify(rational),
        "charpoly_generic6": lambda: generic.charpoly(t),
        "charpoly_shared4": lambda: shared.charpoly(t),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lib", choices=("symox", "sympy"), default="symox")
    parser.add_argument("--repeats", type=int, default=11)
    parser.add_argument("--case", action="append")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    sp = importlib.import_module(args.lib)
    available = cases(sp)
    report = {"library": args.lib, "version": sp.__version__, "repeats": args.repeats, "cases": {}}
    for name in args.case or available:
        if name not in available:
            parser.error(f"unknown case: {name}")
        fn = available[name]
        fn()
        samples = []
        for _ in range(args.repeats):
            start = time.perf_counter()
            result = fn()
            samples.append(time.perf_counter() - start)
        expression = result.as_expr() if hasattr(result, "as_expr") else result
        report["cases"][name] = {
            "minimum": min(samples), "median": statistics.median(samples),
            "samples": samples, "result": str(expression),
        }
        print(f"{name:24s} min {min(samples)*1e3:8.3f} ms   median {statistics.median(samples)*1e3:8.3f} ms", flush=True)
    if args.out:
        args.out.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
