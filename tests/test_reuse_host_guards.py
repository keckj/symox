"""Cross-call result reuse never involves host instances of user classes.

The bounded reuse of exact arithmetic, Add construction, substitution,
expansion, numerical evaluation and NumPy kernels keys on native identities.
Host nodes dispatch user hooks on every call and own Python instances that
must be collected once unreferenced, so they always take the ordinary path.
"""
import gc
import weakref

import numpy as np
import sympy
import symox as sx


class Colored(sx.Symbol):
    pass


class Tagged(sx.Function):
    def __new__(cls, *args, tag=None):
        obj = super().__new__(cls, *args)
        obj.tag = tag
        return obj


def test_host_instances_are_not_retained_by_reuse():
    x = sx.Symbol('x')
    h = Tagged(sx.Symbol('t'), tag='kept')
    ref = weakref.ref(h)
    values = []
    for _ in range(2):
        values += [h + x, h * x, h ** 2, (h + x) * (h + x), sx.Add(2 * h, 3 * x ** 2),
                   sx.expand((h + x) ** 3)]
    assert all(isinstance(value, sx.Basic) for value in values)
    del h, values
    gc.collect()
    assert ref() is None


def test_add_reuse_keeps_host_and_plain_operands_apart():
    x = sx.Symbol('x')
    c = Colored('x')
    terms = (2 * x, 3 * x ** 2)
    plain = sx.Add(*terms)
    assert sx.Add(*terms) is plain
    hosted_terms = (2 * c, 3 * c ** 2)
    hosted = sx.Add(*hosted_terms)
    assert hosted != plain and sx.Add(*hosted_terms) == hosted
    assert all(isinstance(s, Colored) for s in hosted.free_symbols)


def _repeated_hooks(sp):
    x, y = sp.symbols('x y')
    calls = []

    class Raising(sp.Function):
        def _eval_subs(self, old, new):
            calls.append(old)
            raise ValueError('hook')

    class Power(sp.Expr):
        def _eval_power(self, e):
            return sp.Integer(8) if e == 3 else None

    class Value(sp.Function):
        def _eval_evalf(self, prec):
            return sp.Float(1.5)

    outcomes = []
    for _ in range(2):
        try:
            (Raising(x) * x + y).subs(y, 2)
            outcomes.append('no error')
        except ValueError as error:
            outcomes.append(str(error))
    p = Power(y)
    outcomes += [sp.srepr(p ** 3), sp.srepr(p ** 3), sp.srepr(p ** 2 * x)]
    outcomes += [str((Value(x) + 1).evalf()), str((Value(x) + 1).evalf())]
    return outcomes


def test_hooks_dispatch_again_on_repeated_calls():
    assert _repeated_hooks(sx) == _repeated_hooks(sympy)


def test_numpy_kernels_keep_host_functions_on_the_reference_path():
    x = sx.Symbol('x')

    plain = sx.lambdify(x, sx.sin(x) + x ** 2, 'numpy')
    values = np.linspace(0.0, 1.0, 5)
    assert np.allclose(plain(values), np.sin(values) + values ** 2)
    expr = sx.sin(x) + Colored('x') ** 2
    hosted = sx.lambdify([x, Colored('x')], expr, 'numpy')
    assert np.allclose(hosted(values, values), np.sin(values) + values ** 2)
