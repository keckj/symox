"""Objects forwarding attribute access to an expression, compared with SymPy."""
import pytest
import sympy as sp
import symox as sx


class Proxy:
    def __init__(self, expr):
        self.expr = expr

    def __getattr__(self, name):
        return getattr(self.expr, name)


_CASES = {
    'Add': lambda m, x, y, p: m.Add(x, p),
    'Mul': lambda m, x, y, p: m.Mul(2, p),
    'diff': lambda m, x, y, p: m.diff(p, x),
    'simplify': lambda m, x, y, p: m.simplify(p),
    'factor': lambda m, x, y, p: m.factor(p),
    'subs': lambda m, x, y, p: p.subs(x, 2),
    'free symbols': lambda m, x, y, p: sorted(map(str, p.free_symbols)),
}


def _env(lib):
    x, y = lib.symbols('x y')
    return lib, x, y, Proxy(x**2 + y)


# SymPy deprecates non-Expr arguments of Add and Mul but still accepts them.
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.parametrize('name', list(_CASES))
def test_forwarding_proxy_operations_match_sympy(name):
    assert str(_CASES[name](*_env(sx))) == str(_CASES[name](*_env(sp)))


def test_sympify_returns_a_forwarding_proxy_unchanged():
    for lib in (sx, sp):
        proxy = _env(lib)[3]
        assert lib.sympify(proxy) is proxy


def test_forwarding_marker_does_not_capture_unrelated_objects():
    class Other:
        _symox_forwarded_ = ('not an expression',)
    with pytest.raises(sx.SympifyError):
        sx.sympify(Other())
    # SymPy objects keep their own conversion path.
    assert sx.sympify(sp.Symbol('x') + 1) == sx.Symbol('x') + 1
