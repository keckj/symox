"""integrate(): symox vs SymPy on the differential test corpus."""

import pathlib
import sys
import time

import sympy
import symox

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent / "tests"))
from compare_integrate import INDEFINITE, ns  # noqa: E402


def main():
    sx, mx = ns(sympy), ns(symox)
    cases = [c for c in INDEFINITE]
    fs = [eval(c, sx) for c in cases]
    fm = [eval(c, mx) for c in cases]
    sympy.core.cache.clear_cache()
    t0 = time.perf_counter()
    for f in fs:
        sympy.integrate(f, sx["x"])
    ts = time.perf_counter() - t0
    t0 = time.perf_counter()
    for f in fm:
        symox.integrate(f, mx["x"])
    tm = time.perf_counter() - t0
    print(f"{len(cases)} integrals: sympy {ts*1e3:.0f} ms, symox {tm*1e3:.2f} ms, speedup {ts/tm:.0f}x")


if __name__ == "__main__":
    main()
