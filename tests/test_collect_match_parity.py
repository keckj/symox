"""Live parity of ``collect`` with Wild patterns against SymPy 1.14.

``collect`` matches each factor of a term against the pattern with
``Expr.match`` in the reference; the native ``match_factors`` now goes
through the port of ``matches`` (``matches::reference_match``) instead of the
older structural matcher, which missed shifted or negated bases.
"""
import pytest

import sympy
import symox


def _cases(m):
    x, y, z, r = m.symbols('x y z r')
    a, b = m.Wild('a'), m.Wild('b')
    return [
        lambda: m.collect((x - 2)*y + z, -x + a, evaluate=False),
        lambda: m.collect((x - 2)**2*y + z, (-x + a)**2, evaluate=False),
        lambda: m.collect(y*(x - 2) + y*z + z, a - x),
        lambda: m.collect(x**r*y + x**r*z + 1, a**2, evaluate=False),
        lambda: m.collect(3*y + z*m.pi, a**(-2), evaluate=False),
        lambda: m.collect(m.exp(x)*y + m.exp(x)*z, x**a*m.exp(x*b), evaluate=False),
        lambda: m.collect(y*m.sin(x)**2 + z*m.sin(x)**2, m.sin(a)*m.sin(b), evaluate=False),
        lambda: m.collect(a*x**2 + b*x + x**2*y + x*z, x**a, evaluate=False),
        lambda: m.collect((x + 1)*y + (x + 1)*z, a + 1, evaluate=False),
        lambda: m.collect(x*y + x*z + 2*y, a*y, evaluate=False),
    ]


def _show(value):
    if isinstance(value, dict):
        return sorted((str(k), str(v)) for k, v in value.items())
    return str(value)


def _answer(case):
    try:
        return _show(case())
    except Exception as error:
        return 'EXC %s' % type(error).__name__


@pytest.mark.parametrize("index", range(len(_cases(sympy))))
def test_collect_wild_patterns_agree_with_sympy(index):
    assert _answer(_cases(symox)[index]) == _answer(_cases(sympy)[index])
