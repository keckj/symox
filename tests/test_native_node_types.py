"""Nominal SymPy classes of native nodes and the conversion protocols they serve.

SageMath and Cadabra identify SymPy objects through ``type(expr)``, its
``repr``, ``__module__`` and MRO, attributes patched onto SymPy classes after
import, and the polynomial and domain-element protocols. These tests pin
that behavior without the downstream packages.
"""
import fractions
import numbers
import pickle
import sys
import types

import mpmath
import pytest

import symox as s
from symox._node_types import NodeType, native_module, node_class


@pytest.fixture
def facade(monkeypatch):
    """Simulate the routing facade: ``sympy`` names the native package."""
    module = types.ModuleType('sympy')
    module.__symox_version__ = s.__version__
    monkeypatch.setitem(sys.modules, 'sympy', module)
    return module


def _names(cls, count):
    return [base.__name__ for base in cls.__mro__[:count]]


def test_types_follow_sympy_names_and_mro():
    x, y = s.symbols('x y')
    f = s.Function('f')
    cases = [
        (x + y, 'Add', 'symox.core.add', ['Add', 'Expr', 'AssocOp', 'Basic']),
        (x*y, 'Mul', 'symox.core.mul', ['Mul', 'Expr', 'AssocOp', 'Basic']),
        (x**2, 'Pow', 'symox.core.power', ['Pow', 'Expr', 'Basic']),
        (s.Integer(3), 'Integer', 'symox.core.numbers', ['Integer', 'Rational', 'Number']),
        (s.Rational(1, 2), 'Half', 'symox.core.numbers', ['Half', 'RationalConstant', 'Rational']),
        (s.Float(1.5), 'Float', 'symox.core.numbers', ['Float', 'Number']),
        (x, 'Symbol', 'symox.core.symbol', ['Symbol', 'AtomicExpr', 'Atom', 'Expr']),
        (s.Derivative(f(x), x), 'Derivative', 'symox.core.function', ['Derivative', 'Expr', 'Basic']),
        (s.Eq(x, 1), 'Equality', 'symox.core.relational', ['Equality', 'Relational']),
        (s.sin(x), 'sin', 'symox.functions.elementary.trigonometric',
         ['sin', 'TrigonometricFunction', 'DefinedFunction', 'Function']),
    ]
    for expr, name, module, mro in cases:
        cls = type(expr)
        assert isinstance(cls, NodeType)
        assert expr.__class__ is cls
        assert cls.__name__ == name and cls.__module__ == module
        assert _names(cls, len(mro)) == mro
        assert isinstance(expr, s.Basic)
    applied = type(f(x))
    assert applied.__name__ == 'f' and repr(applied) == 'f'
    assert 'Function' in [base.__name__ for base in applied.__mro__]
    assert repr(type(x + y)) == "<class 'symox.core.add.Add'>"
    assert repr(type(s.sin(x))) == 'sin'


def test_facade_reports_sympy_modules(facade):
    x, y = s.symbols('x y')
    assert repr(type(x + y)) == "<class 'sympy.core.add.Add'>"
    assert type(x + y).__module__ == 'sympy.core.add'
    assert str(type(x)) == "<class 'sympy.core.symbol.Symbol'>"
    assert str(type(s.Integer(2))) == "<class 'sympy.core.numbers.Integer'>"
    assert repr(type(s.FiniteSet(1))) == "<class 'sympy.sets.sets.FiniteSet'>"
    assert type(s.FiniteSet(1)).__module__ == 'sympy.sets.sets'
    # The class dictionary keeps the defining module.
    assert native_module(type(x + y)) == 'symox.core.add'
    # The metaclasses report their own module as a plain string.
    assert isinstance(NodeType.__module__, str)
    assert NodeType.__module__ == 'symox._node_types'


def test_facade_does_not_change_native_behavior(facade):
    x, y = s.symbols('x y')
    f = s.Function('f')
    expr = s.sin(x) + x*y + 2
    assert sorted(map(str, expr.atoms(type(x)))) == ['x', 'y']
    assert expr.has(type(x)) and expr.count(type(x)) == 3
    assert expr.find(type(s.sin(x))) == {s.sin(x)}
    assert s.Matrix([[x, s.sin(x)], [1, 2]])*s.Matrix([1, x]) == s.Matrix([x + x*s.sin(x), 2*x + 1])
    assert s.sympify(s.FiniteSet(1, 2)) == s.FiniteSet(1, 2)
    assert s.sympify(s.EmptySet) == s.EmptySet
    assert s.sympify(s.Interval(0, 1)) == s.Interval(0, 1)
    assert s.FiniteSet(x, 1).subs(x, 2) == s.FiniteSet(1, 2)
    assert s.Tuple(x, 1).xreplace({x: 3}) == s.Tuple(3, 1)
    assert (f(x) + x).atoms(type(f(x))) == {f(x)}
    assert s.Matrix([[x]]).applyfunc(lambda t: t + 1) == s.Matrix([[x + 1]])
    assert str(s.groebner([x**2 - y, x*y - 1], x, y)) == \
        "GroebnerBasis([x - y**2, y**3 - 1], x, y, domain='ZZ', order='lex')"


def test_class_patches_reach_instances_and_follow_the_mro():
    x = s.Symbol('x')
    assert 'symox_probe' not in dir(x)
    s.Basic.symox_probe = lambda self: 'basic'
    try:
        assert x.symox_probe() == 'basic'
        assert (x + 1).symox_probe() == 'basic'
        s.Symbol.symox_probe = lambda self: 'symbol'
        # Patched after the instances exist, and the subclass wins.
        assert x.symox_probe() == 'symbol'
        assert (x + 1).symox_probe() == 'basic'
        del s.Symbol.symox_probe
        assert x.symox_probe() == 'basic'
    finally:
        del s.Basic.symox_probe
    assert not hasattr(x, 'symox_probe')


def test_sage_placeholder_is_resolved_on_the_class():
    x = s.Symbol('x')
    assert node_class('Basic').__dict__['_sage_'].__name__ == '_sage_'
    s.Symbol._sage_ = lambda self: ('sage', self.name)
    try:
        assert x._sage_() == ('sage', 'x')
        assert s.Symbol._sage_(x) == ('sage', 'x')
    finally:
        del s.Symbol._sage_
    assert x._sage_.__func__ is node_class('Basic').__dict__['_sage_']


def test_rebuild_pickle_numbers_and_subclass_checks():
    x, y = s.symbols('x y')
    assert type(x + y)(x, y) == x + y
    assert type(s.sin(x))(y) == s.sin(y)
    for cls in (type(x + y), type(s.Integer(2)), type(s.sin(x))):
        assert pickle.loads(pickle.dumps(cls)) is cls
    assert pickle.loads(pickle.dumps(x + y)) == x + y
    assert isinstance(s.Integer(3), numbers.Integral)
    assert isinstance(s.Rational(1, 3), numbers.Rational)
    assert isinstance(s.Float(1.5), numbers.Real)
    assert issubclass(type(x + y), s.Add) and issubclass(type(x), s.Basic)
    assert issubclass(type(s.Integer(2)), s.Rational) and not issubclass(type(x), s.Add)
    dispatch = {s.Add: 'add', s.Symbol: 'symbol'}
    assert [next(v for k, v in dispatch.items() if issubclass(type(e), k)) for e in (x + y, x)] == ['add', 'symbol']


def test_named_sets_keep_sympy_class_names():
    assert type(s.EmptySet).__name__ == 'EmptySet'
    assert type(s.S.Reals).__name__ == 'Reals'
    assert type(s.S.Naturals).__name__ == 'Naturals'
    assert type(s.FiniteSet(1)).__name__ == 'FiniteSet'
    assert type(s.Interval(0, 1)).__name__ == 'Interval'


def test_symbol_name_is_a_persistent_string():
    # SymPy stores ``name`` on the instance; Cadabra's bridge keeps a
    # borrowed handle to it.
    x = s.Symbol('tau')
    assert x.name is x.name
    assert s.Symbol('tau').name is x.name
    assert s.MatrixSymbol('A', 2, 2).name is s.MatrixSymbol('A', 2, 2).name
    assert s.Function('f')(x).name == 'f'


def test_derivative_and_function_protocols_for_cadabra():
    tau = s.Symbol('tau')
    r = s.Function('r')(tau)
    second = s.parsing.sympy_parser.parse_expr('diff(diff(r(tau), tau), tau)')
    assert type(second).__name__ == 'Derivative'
    assert second.expr == r and second.variables == (tau, tau)
    assert isinstance(second.variables, tuple)
    from symox.printing.str import StrPrinter
    printer = StrPrinter()
    x = s.Symbol('x')
    assert printer._as_ordered_terms(x**2 + x + 1) == [x**2, x, 1]
    assert printer._as_ordered_terms(x**2 + x + 1, order='none') == list((x**2 + x + 1).args)
    assert printer._as_ordered_terms(x**2 + x + 1, order='old') == [1, x, x**2]


def test_constant_coefficient_ode_uses_exact_characteristic_roots():
    tau, m = s.symbols('tau m')
    X = s.Function('X')
    solution = s.dsolve(s.Eq(X(tau).diff(tau, 2) + m**2*X(tau), 0), X(tau))
    C1, C2 = s.symbols('C1 C2')
    assert solution.rhs == C1*s.exp(-s.I*m*tau) + C2*s.exp(s.I*m*tau)


def test_root_intervals():
    from symox.polys.rootisolation import ComplexInterval, RealInterval
    x = s.Symbol('x')
    root = s.CRootOf(x**3 - 2, 0)
    interval = root._get_interval()
    assert isinstance(interval, RealInterval)
    assert interval.a <= interval.b and interval.dx == interval.b - interval.a
    refined = interval.refine_size(s.Rational(1, 100))
    assert refined.dx < s.Rational(1, 100)
    assert refined.a**3 <= 2 <= refined.b**3
    root._set_interval(refined)
    assert root._get_interval() is refined
    field = s.QQ
    corner = ComplexInterval((field(0), field(1)), (field(1), field(2)), None, None, None, None,
                             None, None, s.ZZ)
    assert corner.conjugate().as_tuple() == ((0, -2), (1, -1))
    assert corner.is_disjoint(corner.conjugate())


def test_float_reads_mpmath_objects_without_their_sympy_hook():
    class Foreign:
        def _mpmath_(self, prec, rounding):
            return mpmath.mpf('1.25')

        def _sympy_(self):
            raise AssertionError('Float must not call _sympy_')

    assert s.Float(Foreign()) == s.Float('1.25')
    assert s.Float(Foreign(), 30) == s.Float('1.25', 30)


def test_sympify_returns_the_hook_result_unchanged():
    marker = object()

    class Hook:
        def _sympy_(self):
            return marker

    assert s.sympify(Hook()) is marker


def test_matrix_instances_have_a_dictionary():
    matrix = s.Matrix([[1, 2]])
    matrix.cache = 'value'
    assert matrix.__dict__ == {'cache': 'value'}


def test_poly_rep_is_a_dense_view():
    from symox.polys.polyclasses import DMP
    x, y = s.symbols('x y')
    rep = s.Poly(x**2 + 2*x + 3, x).rep
    assert type(rep).__name__ == 'DMP_Python' and isinstance(rep, DMP)
    assert rep.dom == s.ZZ and rep.lev == 0
    assert rep.to_list() == [1, 2, 3]
    assert rep.terms() == [((2,), 1), ((1,), 2), ((0,), 3)]
    assert s.Poly(0, x).rep.terms() == [((0,), 0)]
    assert s.Poly(x*y + 1, x, y).rep.terms() == [((1, 1), 1), ((0, 0), 1)]
    assert pickle.loads(pickle.dumps(rep)) == rep and hash(rep) == hash(DMP([1, 2, 3], s.ZZ))
    q, r = divmod(DMP([1, 2, 0], s.ZZ), DMP([2, 0], s.ZZ))
    assert q.to_list() == [] and r.to_list() == [1, 2, 0]
    assert DMP([1, 2], s.ZZ)*DMP([1, 0], s.ZZ) == DMP([1, 2, 0], s.ZZ)


class _Ring:
    """A minimal SageMath-like ring over Python numbers."""

    def __init__(self):
        self._base = lambda c: fractions.Fraction(int(c.numerator), int(c.denominator))

    def __call__(self, value):
        if isinstance(value, str):
            return s.Symbol(value)
        return s.sympify(value)

    def base_ring(self):
        return self._base

    def zero(self):
        return s.Integer(0)


def test_domain_elements_follow_sage_conversion_protocol(monkeypatch):
    assert s.ZZ(3)._integer_(int) == 3
    assert s.QQ(4, 2)._integer_(int) == 2
    with pytest.raises(TypeError):
        s.QQ(1, 2)._integer_(int)
    rational = types.ModuleType('sage.rings.rational_field')
    rational.QQ = lambda pair: fractions.Fraction(*pair)
    monkeypatch.setitem(sys.modules, 'sage', types.ModuleType('sage'))
    monkeypatch.setitem(sys.modules, 'sage.rings', types.ModuleType('sage.rings'))
    monkeypatch.setitem(sys.modules, 'sage.rings.rational_field', rational)
    assert s.QQ(3, 6)._rational_() == fractions.Fraction(1, 2)
    x, y = s.symbols('x y')
    ring = s.QQ[x, y]
    element = ring.from_sympy(x**2*y/2 + 3)
    assert element._polynomial_(_Ring()) == x**2*y/2 + 3
    with pytest.raises(TypeError):
        s.ZZ(3)._polynomial_(_Ring())
    with pytest.raises(TypeError):
        ring.from_sympy(x)._rational_()
