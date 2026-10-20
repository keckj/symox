"""Native nodes have SymPy's classes, singletons and class-level protocols.

Every case runs on SymPy and Symox and compares the normalized result, or
the exception type and message.
"""
import functools
import importlib
import inspect
import pickle
import warnings

import pytest
import sympy
import symox

MODULES = {'sympy': sympy, 'symox': symox}


def outcome(case, name):
    module = MODULES[name]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            result = case(module)
        except Exception as error:  # noqa: BLE001 - the error is compared
            return ('error', type(error).__name__, str(error).replace(name, 'M'))
    return ('value', repr(result).replace(name, 'M'))


def submodule(module, path):
    return importlib.import_module(module.__name__ + '.' + path)


def trig(m):
    return submodule(m, 'functions.elementary.trigonometric')


def hyperbolic(m):
    return submodule(m, 'functions.elementary.hyperbolic')


def lattice(m):
    return submodule(m, 'core.operations').LatticeOp


def application(m):
    return submodule(m, 'core.function').Application


def dispatch(m, value):
    @functools.singledispatch
    def kind(e):
        return 'other'

    @kind.register(m.sin)
    def _(e):
        return 'sin'

    @kind.register(trig(m).TrigonometricFunction)
    def _(e):
        return 'trig'

    @kind.register(m.And)
    def _(e):
        return 'And'

    @kind.register(m.Max)
    def _(e):
        return 'Max'

    @kind.register(lattice(m))
    def _(e):
        return 'Lattice'

    return kind(value(m))


def subclass(m, base):
    return type('Sub', (getattr(m, base),), {})


def xy(m):
    return m.symbols('x y')


CASES = {
    # singletons
    'zero identity': lambda m: (m.Integer(0) is m.S.Zero, m.S(1) is m.S.One, -m.S.One is m.S.NegativeOne),
    'half identity': lambda m: m.Rational(1, 2) is m.S.Half,
    'constant identity': lambda m: (m.oo is m.S.Infinity, m.nan is m.S.NaN, m.zoo is m.S.ComplexInfinity, m.pi is m.S.Pi),
    'boolean identity': lambda m: (m.true is m.S.true, m.sympify(True) is m.true, m.Eq(1, 1) is m.true),
    'singleton classes': lambda m: (type(m.S.Zero).__name__, type(m.S.Half).__name__, type(m.oo).__name__, type(m.true).__name__),
    'singleton isinstance': lambda m: (isinstance(m.S.Zero, m.Integer), isinstance(m.S.Half, m.Rational), isinstance(m.oo, m.Number)),
    'singleton pickle': lambda m: pickle.loads(pickle.dumps(m.S.Zero)) is m.S.Zero,
    # function heads
    'head is class': lambda m: (inspect.isclass(m.Abs), inspect.isclass(m.sin), repr(m.Abs), str(m.sin)),
    'head isinstance': lambda m: (isinstance(m.Abs(xy(m)[0]), m.Abs), isinstance(m.sin(xy(m)[0]), m.Abs)),
    'head func': lambda m: m.Abs(xy(m)[0]).func is m.Abs,
    'head type call': lambda m: type(m.sin(xy(m)[0]))(xy(m)[1]),
    'head subclasses': lambda m: (issubclass(m.sin, m.Function), issubclass(m.sin, m.Expr), issubclass(m.sin, m.Basic)),
    'head application': lambda m: issubclass(m.Abs, application(m)),
    'trig bases': lambda m: (issubclass(m.sin, trig(m).TrigonometricFunction),
                             isinstance(m.asin(xy(m)[0]), trig(m).InverseTrigonometricFunction)),
    'hyperbolic bases': lambda m: (issubclass(m.sinh, hyperbolic(m).HyperbolicFunction),
                                   isinstance(m.asinh(xy(m)[0]), hyperbolic(m).InverseHyperbolicFunction)),
    'head pickle': lambda m: (pickle.loads(pickle.dumps(m.Abs)) is m.Abs,
                              pickle.loads(pickle.dumps(m.Abs(xy(m)[0]))) == m.Abs(xy(m)[0])),
    'dispatch sin': lambda m: dispatch(m, lambda m: m.sin(xy(m)[0])),
    'dispatch cos': lambda m: dispatch(m, lambda m: m.cos(xy(m)[0])),
    'dispatch exp': lambda m: dispatch(m, lambda m: m.exp(xy(m)[0])),
    'nargs': lambda m: (m.sin.nargs, m.log.nargs),
    'unevaluated head': lambda m: m.Abs(-xy(m)[0], evaluate=False),
    'subclass instance': lambda m: (subclass(m, 'Abs')(xy(m)[0]), type(subclass(m, 'Abs')(xy(m)[0])).__name__,
                                    isinstance(subclass(m, 'Abs')(xy(m)[0]), m.Abs)),
    'subclass evaluates': lambda m: subclass(m, 'Abs')(-3),
    'subclass diff': lambda m: subclass(m, 'sin')(xy(m)[0]).diff(xy(m)[0]),
    'subclass zero': lambda m: subclass(m, 'sin')(0),
    'subclass trig base': lambda m: issubclass(subclass(m, 'sin'), trig(m).TrigonometricFunction),
    'subclass subs': lambda m: subclass(m, 'Abs')(xy(m)[0]).subs(xy(m)[0], -2),
    # gates and lattices
    'gate dispatch': lambda m: tuple(dispatch(m, value) for value in (
        lambda m: m.And(*xy(m)), lambda m: m.Max(*xy(m)), lambda m: m.Min(*xy(m)), lambda m: xy(m)[0])),
    'gate func': lambda m: (m.And(*xy(m)).func is m.And, (~xy(m)[0]).func is m.Not,
                            m.Max(*xy(m)).func is m.Max, m.Xor(*xy(m)).func is m.Xor),
    'gate isinstance': lambda m: (isinstance(m.And(*xy(m)), m.And), isinstance(m.And(*xy(m)), m.Or),
                                  isinstance(m.Max(*xy(m)), m.Max), isinstance(m.Max(*xy(m)), m.Min),
                                  isinstance(m.Max(*xy(m)), m.Expr), isinstance(m.And(*xy(m)), m.Expr)),
    'gate issubclass': lambda m: (issubclass(m.And, lattice(m)), issubclass(m.Max, m.Expr), issubclass(m.And, m.Expr),
                                  issubclass(type(m.Max(*xy(m))), m.Max)),
    'lattice elements': lambda m: (m.And.zero, m.And.identity, m.Or.zero, m.Max.zero, m.Min.identity,
                                   m.Max(*xy(m)).zero),
    'gate repr': lambda m: (repr(m.Max), str(m.Min), repr(m.And), str(m.Not)),
    'gate pickle': lambda m: (pickle.loads(pickle.dumps(m.And)) is m.And, pickle.loads(pickle.dumps(m.Max)) is m.Max),
    'make_args': lambda m: (m.And.make_args(m.And(*xy(m))) == frozenset(xy(m)), m.Or.make_args(m.And(*xy(m)))),
    'lattice values': lambda m: (m.Max(1, 2), m.Min(xy(m)[0], m.oo), m.Max(), m.Min(), m.Max(xy(m)[0])),
    'lattice nested': lambda m: m.Min(xy(m)[0], m.Max(xy(m)[1], m.Min(*xy(m)))),
    'lattice unevaluated': lambda m: (m.Max(xy(m)[0], xy(m)[0], evaluate=False), m.Max(1, 2, evaluate=False).args),
    'lattice not comparable': lambda m: m.Max(m.I, 1),
    'lattice subclass': lambda m: (subclass(m, 'Max')(*xy(m)), type(subclass(m, 'Max')(*xy(m))).__name__,
                                   subclass(m, 'Max')(1, 2), isinstance(subclass(m, 'Max')(*xy(m)), m.Max)),
    'gate subclass': lambda m: type(subclass(m, 'And')(*xy(m))).__name__,
    'class equality': lambda m: (type(m.Eq(*xy(m))) == m.Eq, type(m.Eq(*xy(m))) != m.Ne, type(m.sin(xy(m)[0])) == m.sin,
                                 {m.Eq: 'Eq', m.Lt: 'Lt'}.get(type(m.Lt(*xy(m)))), {m.Add: 1}.get(type(m.Add(*xy(m)))),
                                 hash(type(m.Add(*xy(m)))) == hash(m.Add), type(m.Max(*xy(m))) == m.Max),
    'class flags': lambda m: (m.Float(1).is_symbol, xy(m)[0].is_symbol, m.Dummy().is_symbol, m.pi.is_NumberSymbol,
                              m.E.is_NumberSymbol, m.Integer(2).is_NumberSymbol, m.Not(xy(m)[0]).is_Not,
                              m.And(*xy(m)).is_Not, (xy(m)[0] + 1).is_scalar, m.Integer(2).is_AlgebraicNumber,
                              m.Max(*xy(m)).is_Indexed, m.sin(xy(m)[0]).is_Poly),
    'function subclass is_number': lambda m: (subclass(m, 'Function')(1).is_number, subclass(m, 'Function')(1, 2.5).is_number,
                                              subclass(m, 'Function')(xy(m)[0]).is_number, m.Function('h')(1).is_number),
    'gate boolean function': lambda m: (isinstance(m.And(*xy(m)), m.logic.boolalg.BooleanFunction),
                                        issubclass(m.Not, m.logic.boolalg.BooleanFunction)),
}


@pytest.mark.parametrize('name', sorted(CASES))
def test_class_protocol_matches_sympy(name):
    case = CASES[name]
    assert outcome(case, 'symox') == outcome(case, 'sympy')


def test_kind_of_shared_tree_is_not_unfolded():
    x, one = symox.Symbol('x'), symox.Integer(1)
    value = symox._core._restore_arithmetic('Mul', symox.Integer(2), [(x, one)], False)
    for _ in range(80):
        value = symox._core._restore_arithmetic('Mul', one, [(value, one), (value, one)], False)
    assert value.kind is symox.core.kind.NumberKind
