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


def _expansion_hooks(sp):
    x = sp.Symbol('x')

    class Doubled(sp.Function):
        def _eval_expand_double(self, **hints):
            return 2*self.args[0]

    expr = Doubled(x) + x
    return [str(sp.expand(expr, double=True)), str(sp.expand(expr, double=True)),
            str(expr.expand(double=True)), str(sp.expand((x + 1)**2))]


def test_expansion_hooks_survive_the_native_node_fast_path():
    assert _expansion_hooks(sx) == _expansion_hooks(sympy)


def test_nominal_classes_take_the_reuse_paths():
    from symox import _core
    from symox._node_types import is_native_node
    x = sx.Symbol('x')
    plain = (x, 2*x, x**2, sx.Integer(3))
    assert all(is_native_node(term) for term in plain)
    assert _core._exact_expr_key(plain) is not None
    assert _core._exact_expr_args(plain)
    hosted = (Colored('c'), 2*x)
    assert not is_native_node(hosted[0])
    assert _core._exact_expr_key(hosted) is None
    assert _core._exact_expr_key((2*Colored('c'), x)) is None
    assert not _core._exact_expr_args(hosted)


def _forwarding_wrapper(module):
    class Wrapped:
        # A host wrapper that forwards every attribute to the expression it
        # holds, like Lcapy's expressions; ``sympify`` returns it unchanged.
        def __init__(self, value):
            self.sympy = module.Integer(value)

        def __getattr__(self, attr):
            return getattr(self.sympy, attr)

    return Wrapped


def test_dense_polynomial_construction_keeps_forwarding_wrappers_on_the_general_path():
    results = []
    for module in (sympy, sx):
        Wrapped = _forwarding_wrapper(module)
        x = module.Symbol('x')
        poly = module.Poly(module.sympify([Wrapped(1), Wrapped(2), Wrapped(3)]), x)
        results.append((str(poly.as_expr()), str(poly.domain), poly.all_coeffs()))
    assert str(results[0]) == str(results[1])
