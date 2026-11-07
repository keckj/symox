"""Live parity of ``has`` / ``has_free`` / ``has_xfree`` with SymPy 1.14.

The corpus (audits/performance/2026-11-06-hasmatch/corpus.json.gz) holds
hand-written cases (sub-products and sub-sums, noncommutative operands,
powers, classes, Wild, unapplied functions, several patterns, derivatives,
integrals and other binders, relationals, Piecewise, Max/Min, sets,
matrices) and seeded random sub-collections of the arguments of sums and
products.  Every case is evaluated by SymPy and by symox in this process.
"""
import gzip
import json
import pathlib

import pytest

import sympy
import symox

ROOT = pathlib.Path(__file__).resolve().parents[1]
AUDIT = ROOT / "audits" / "performance" / "2026-11-06-hasmatch"


def _namespace(module):
    ns = dict(vars(module))
    x, y, z, w = module.symbols('x y z w')
    A, B, C = module.symbols('A B C', commutative=False)
    ns.update(x=x, y=y, z=z, w=w, A=A, B=B, C=C,
              p=module.Symbol('p', positive=True), f=module.Function('f'), g=module.Function('g'),
              a=module.Wild('a'), b=module.Wild('b'), n=module.Symbol('n', integer=True),
              M=module.MatrixSymbol('M', 2, 2), N=module.MatrixSymbol('N', 2, 2))
    return ns


def _answer(code, ns):
    try:
        return repr(eval(code, ns))
    except Exception as error:
        return 'EXC %s' % type(error).__name__


def _cases():
    return json.loads(gzip.open(AUDIT / "corpus.json.gz").read())["cases"]


def test_has_corpus_agrees_with_sympy():
    ref_ns, ns = _namespace(sympy), _namespace(symox)
    different = []
    cases = _cases()
    for case in cases:
        expected, got = _answer(case["code"], ref_ns), _answer(case["code"], ns)
        if expected != got:
            different.append((case["code"], expected, got))
    assert len(cases) >= 1200
    assert different == [], different[:20]


@pytest.mark.parametrize("code", [
    "(2*x*y*z).has(x*y)",
    "(2*x*y*z).has(2*x)",
    "(-x*y).has(x*y)",
    "(x*y*sin(x*y*z)).has(x*z)",
    "Integral(x*y*z, (x, 0, 1)).has(x*z)",
    "Matrix([[x*y*z]]).has(x*z)",
    "Max(x, y, z).has(Max(x, z))",
    "(A + B + C).has(A + C)",
    "(A*B*C).has(A*C)",
    "(x*A*B*C).has(x*B*C)",
    "x.has('x')",
    "x.has([x])",
    "(A*B*C).has_free(A*C)",
    "(x*y*z).has_free(x*z)",
    "Integral(x*y*z, (x, 0, 1)).has_free(y*z)",
    "Max(x, y, z).has_free(Max(x, z))",
])
def test_has_reference_cases(code):
    assert _answer(code, _namespace(symox)) == _answer(code, _namespace(sympy))
