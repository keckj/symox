"""All real Cardano roots must remain real through numerical evaluation."""
import pickle
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('source', ['x**3-3*x+1', 'x**3-7*x+1', 'x**3-10*x+2'])
@pytest.mark.parametrize('digits', [30, 80, 200])
def test_every_real_cubic_root_evaluates_at_requested_precision(source, digits):
    x = sx.Symbol('x')
    roots = sx.solve(sx.sympify(source), x)
    expected = sp.Poly(sp.sympify(source)).all_roots(radicals=False)
    assert len(roots) == len(expected) == 3
    actual = []
    for root in roots:
        assert root.is_real is True
        value = sx.N(root, digits)
        assert value.is_Float
        assert value._prec == sp.Float(1, digits)._prec
        restored = pickle.loads(pickle.dumps(root))
        assert restored.is_real is True
        assert sx.N(restored, digits) == value
        actual.append(value.to_sympy())
    for value, target in zip(sorted(actual), expected):
        target = sp.N(target, digits+10)
        assert abs(value-target) < abs(target)*sp.Rational(1, 10)**(digits-5)


def test_native_conjugate_proofs_do_not_import_the_reference():
    script = '''
import sys
class NoReference:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0, NoReference())
import symox as sx
x = sx.Symbol('x')
for source in ('x**3-3*x+1', 'x**3-7*x+1'):
    roots = sx.solve(sx.sympify(source), x)
    assert len(roots) == 3
    for root in roots:
        assert root.is_real is True
        assert sx.N(root, 80).is_Float
'''
    subprocess.run([sys.executable, '-c', script], check=True, timeout=30)


def test_cubic_inequality_contains_the_entire_exact_positive_set():
    x = sx.Symbol('x')
    # The endpoints are the Cardano roots of solve (the reference raises
    # NotImplementedError for this inequality). Each is an exact zero, so
    # every interval is open.
    result = sx.solve_univariate_inequality(x**3-3*x+1>0, x, relational=False)
    roots = sx.solve(x**3-3*x+1, x)
    assert isinstance(result, sx.Union) and len(result.args) == 2
    first, second = sorted(result.args, key=lambda i: float(sx.re(i.start).evalf()))
    assert first.left_open and first.right_open and second.left_open
    assert second.end == sx.oo
    assert {first.start, first.end, second.start} == set(roots)
    for endpoint, angle in [(first.start, 8), (first.end, 4), (second.start, 2)]:
        assert abs(complex(endpoint.evalf(30)) - 2*complex(sp.cos(angle*sp.pi/9).evalf(30))) < 1e-25
    # z=exp(2*pi*I/3) on its principal branch. The three endpoints are exactly
    # 2*cos(8*pi/9), 2*cos(4*pi/9), 2*cos(2*pi/9), respectively. Monotonicity
    # of cosine places them in (-2,-1), (0,1), (1,2). Exact minimal polynomials
    # and Sturm counts prove they are ALL distinct roots; the positive leading
    # coefficient then gives the full alternating sign intervals above.
    t = sp.Symbol('t')
    polynomial = sp.Poly(t**3-3*t+1, t)
    for angle, lower, upper in [(8, -2, -1), (4, 0, 1), (2, 1, 2)]:
        algebraic = 2*sp.cos(angle*sp.pi/9)
        assert sp.Poly(sp.minpoly(algebraic, t), t) == polynomial
        assert polynomial.count_roots(lower, upper) == 1
