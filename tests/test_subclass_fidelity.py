"""User subclasses of Symbol, Dummy and Function keep their class and state.

Every test runs against SymPy 1.14 and symox, so the expected values are the
ones the reference produces. The patterns come from EinsteinPy, Optlang,
Lcapy, Galgebra and Devito.
"""
import copy
import gc
import pickle
import weakref

import pytest
import sympy
import symox


@pytest.fixture(params=[sympy, symox], ids=['sympy', 'symox'])
def sp(request):
    return request.param


def _described_symbol(sp):
    # EinsteinPy, Galgebra: a Symbol carrying extra Python state.
    class C(sp.Symbol):
        def __new__(cls, name, descriptive_name=None, **assumptions):
            instance = super(C, cls).__new__(cls, name, **assumptions)
            instance._descriptive_name = descriptive_name
            return instance
    return C


def _variable_class(sp):
    # Optlang: a renamable Dummy whose class is named Symbol.
    class Symbol(sp.Dummy):
        def __new__(cls, name, *args, **kwargs):
            return sp.Dummy.__new__(cls, name)

        def __init__(self, name, lb=None):
            self._name = name
            super().__init__()
            self.lb = lb

        @property
        def name(self):
            return self._name

        @name.setter
        def name(self, value):
            self._name = value
    return Symbol


def _unit_step(sp):
    # Lcapy: an evaluating function with class assumptions and a rewrite rule.
    class UnitStep(sp.Function):
        is_integer = True
        is_real = True

        @classmethod
        def eval(cls, nval):
            if nval.is_Number:
                return sp.S.One if nval >= 0 else sp.S.Zero

        def _eval_rewrite_as_Heaviside(self, arg, **kwargs):
            return sp.Heaviside(arg)
    return UnitStep


def _module_level(factory, label):
    # Pickle looks classes up by qualified name in this module.
    classes = {}
    for module in (sympy, symox):
        cls = factory(module)
        cls.__qualname__ = '%s_%s' % (label, module.__name__)
        globals()[cls.__qualname__] = cls
        classes[module] = cls
    return classes.__getitem__


described_symbol = _module_level(_described_symbol, 'C')
variable_class = _module_level(_variable_class, 'Variable')
unit_step = _module_level(_unit_step, 'UnitStep')



def test_symbol_subclass_keeps_class_and_state(sp):
    C = described_symbol(sp)
    c, x = C('c', 'speed'), sp.Symbol('x')
    assert type(c) is C and c._descriptive_name == 'speed'
    assert isinstance(c, sp.Symbol) and not isinstance(x, C)
    assert c.func is C and c.args == ()
    # SymPy's constructor cache returns the same instance, so the state is
    # passed again to keep it.
    assert c != sp.Symbol('c') and C('c', 'speed') is c and hash(c) == hash(C('c', 'speed'))
    expr = c*x + 1
    assert [s for s in expr.free_symbols if type(s) is C] == [c]
    assert [a for a in (c*x).args if isinstance(a, C)] == [c]
    for rebuilt in (expr.subs(x, 2), (c + x).xreplace({x: 3}), sp.Matrix([[c, x]])[0, 0]):
        (atom,) = [s for s in rebuilt.free_symbols if s.name == 'c']
        assert type(atom) is C and atom._descriptive_name == 'speed'
    assert str(c) == 'c' and sp.latex(c) == 'c' and sp.srepr(c) == "C('c')"
    assert sp.solve(c - 3, c) == [3]


def test_symbol_subclass_pickle_and_copy(sp):
    C = described_symbol(sp)
    c, x = C('c', 'speed'), sp.Symbol('x')
    assert pickle.loads(pickle.dumps(c)) == c
    assert type(pickle.loads(pickle.dumps(c))) is C
    assert copy.copy(c) == c and copy.deepcopy(c) == c
    restored = pickle.loads(pickle.dumps(c + x))
    assert restored == c + x
    assert C in {type(a) for a in restored.args}
    reduced = c.__reduce_ex__(4)
    assert reduced[1][0] is C and reduced[1][1:] == ('c',)


def test_dummy_subclass_prints_and_rebuilds(sp):
    class D(sp.Dummy):
        pass
    d = D('d')
    assert str(d) == '_d' and sp.srepr(d).startswith("D('d'")
    assert type(d + 1 - 1) is D and d.dummy_index is not None
    assert d != D('d')


def test_renamable_dummy_named_symbol(sp):
    Variable = variable_class(sp)
    x = sp.Symbol('x')
    v = Variable('v', lb=3)
    assert v.lb == 3 and str(v + x) == 'x + v'
    v.name = 'w'
    assert str(2*v + x) == 'x + 2*w'
    assert sp.srepr(v).startswith("Symbol('w'")
    (atom,) = [s for s in (2*v + x).free_symbols if s != x]
    assert atom is v


def test_instances_held_by_expressions_stay_alive(sp):
    Variable = variable_class(sp)
    x = sp.Symbol('x')
    expr = Variable('v', lb=7) + x
    gc.collect()
    (atom,) = [s for s in expr.free_symbols if s != x]
    assert atom.lb == 7 and atom.name == 'v'
    assert [s for s in expr.subs(x, 1).free_symbols][0].lb == 7


def test_unreferenced_instance_cycles_are_collected():
    # SymPy's cache keeps such instances; the native tree references are
    # visible to the garbage collector, so symox frees them.
    Variable = variable_class(symox)
    sp = symox
    v = Variable('v')
    v.problem = 2*v + 1
    ref = weakref.ref(v)
    del v
    gc.collect()
    assert ref() is None


def test_add_mul_from_args(sp):
    C = described_symbol(sp)
    x, c = sp.Symbol('x'), C('c')
    assert sp.Add._from_args([x, c]) == x + c
    assert sp.Mul._from_args((x,)) == x
    assert sp.Add._from_args([]) == 0


def test_function_subclass_eval_and_assumptions(sp):
    UnitStep = unit_step(sp)
    n, y = sp.symbols('n y')
    u = UnitStep(n)
    assert type(u) is UnitStep and u.func is UnitStep
    assert UnitStep(3) == 1 and UnitStep(-2) == 0
    assert isinstance(u, sp.Function) and issubclass(UnitStep, sp.Function)
    assert not isinstance(u, sp.core.function.AppliedUndef)
    assert u.is_integer is True and u.is_real is True
    assert str(UnitStep) == 'UnitStep' and UnitStep.nargs == sp.FiniteSet(1)
    assert u.subs(n, 2) == 1 and type(u.subs(n, y)) is UnitStep
    assert (u*n).xreplace({n: 5}) == 5
    assert u != sp.Function('UnitStep')(n)
    with pytest.raises(TypeError, match='UnitStep takes exactly 1 argument .2 given.'):
        UnitStep(1, 2)


def test_function_subclass_rewrite(sp):
    UnitStep = unit_step(sp)
    n, x, t = sp.symbols('n x t')
    assert UnitStep(n).rewrite(sp.Heaviside) == sp.Heaviside(n)

    class tri(sp.Function):
        def _eval_rewrite_as_Heaviside(self, val, **kwargs):
            return sp.Heaviside(val)*val

    class rect(sp.Function):
        def rewrite(self, *args, **hints):
            arg = self.args[0]
            return sp.Heaviside(arg + sp.S.Half) - sp.Heaviside(arg - sp.S.Half)

    assert (2*tri(t) + sp.exp(x)).rewrite(sp.Heaviside) == 2*t*sp.Heaviside(t) + sp.exp(x)
    assert (tri(t) + sp.cos(x)).rewrite(sp.exp) == tri(t) + sp.exp(sp.I*x)/2 + sp.exp(-sp.I*x)/2
    assert tri(sp.sin(t)).rewrite(sp.exp).func is tri
    assert rect(t).rewrite() == sp.Heaviside(t + sp.S.Half) - sp.Heaviside(t - sp.S.Half)


def test_function_subclass_derivatives_and_printing(sp):
    x, y = sp.symbols('x y')

    class G(sp.Function):
        nargs = 2

        def fdiff(self, argindex=1):
            return self.args[1] if argindex == 1 else self.args[0]

        def _latex(self, printer):
            return r'\mathcal{G}(%s)' % printer._print(self.args[0])

        def _sympystr(self, printer):
            return 'GG[%s]' % printer._print(self.args[0])

    g = G(x, y)
    assert g.diff(x) == y and g.diff(y) == x
    assert sp.latex(g + 1) == r'\mathcal{G}(x) + 1'
    assert str(g + 1) == 'GG[x] + 1'
    assert G.nargs == sp.FiniteSet(2)
    assert sp.Matrix([[g]])[0].func is G

    UnitStep = unit_step(sp)
    n = sp.Symbol('n')
    assert sp.srepr(UnitStep(n)) == "UnitStep(Symbol('n'))"
    assert sp.latex(UnitStep(n)) == r'\operatorname{UnitStep}{\left(n \right)}'
    assert UnitStep(n).diff(n) == sp.Derivative(UnitStep(n), n)


def test_function_subclass_eval_is_and_arity(sp):
    x, y = sp.symbols('x y')

    class H(sp.Function):
        def _eval_is_positive(self):
            return True

    assert H(x).is_positive is True and H.nargs == sp.Naturals0
    assert H(x, y).args == (x, y)

    class K(sp.Function):
        @classmethod
        def eval(cls, a, b=1):
            if a == 0:
                return b

    assert K.nargs == sp.FiniteSet(1, 2)
    assert K(0) == 1 and K(0, 5) == 5 and K(x).subs(x, 0) == 1

    class Step(sp.Function):
        is_integer = True

        @classmethod
        def eval(cls, nval, zero=None):
            if nval.is_zero:
                return sp.S.Half if zero is None else zero

    class SubStep(Step):
        pass

    assert SubStep(x).func is SubStep and isinstance(SubStep(x), Step)
    assert SubStep(x) != Step(x) and SubStep(0) == sp.S.Half
    assert Step(sp.S(0), evaluate=False).args == (0,)


def test_function_subclass_sort_order_and_pickle(sp):
    x, n = sp.symbols('x n')

    class rect(sp.Function):
        @classmethod
        def eval(cls, val):
            if val.is_Number:
                return sp.S.Zero if abs(val) > sp.S.Half else sp.S.One

    class tri(sp.Function):
        pass

    UnitStep = unit_step(sp)
    ordered = sorted([rect(x), tri(x), sp.sin(x), x], key=sp.default_sort_key)
    assert [str(e) for e in ordered] == ['x', 'tri(x)', 'sin(x)', 'rect(x)']
    expr = x*UnitStep(n) + 1
    assert pickle.loads(pickle.dumps(expr)) == expr
    assert copy.copy(tri(x)) == tri(x) and copy.deepcopy(UnitStep(n)) == UnitStep(n)
    assert hash(UnitStep(n)) == hash(UnitStep(n))
    assert (rect(x) + 1).replace(rect, tri) == tri(x) + 1
    assert rect(x + 1).atoms(rect) == {rect(x + 1)}
    assert (rect(x) + x).has(rect)


def test_number_constructor(sp):
    assert sp.Number(-4) == -4 and type(sp.Number(-4)) is type(sp.Integer(-4))
    assert sp.Number(1, 2) == sp.Rational(1, 2)
    assert sp.Number('nan') is sp.nan and sp.Number('-inf') == -sp.oo
    assert sp.Number(sp.S(3)) == 3
    with pytest.raises(TypeError):
        sp.Number([1])
    with pytest.raises(ValueError):
        sp.Number('x')


def test_builtin_function_isinstance(sp):
    x = sp.Symbol('x')
    assert isinstance(sp.sin(x), sp.sin)
    assert not isinstance(sp.cos(x), sp.sin)
    assert isinstance(sp.Heaviside(x), (sp.Heaviside, sp.DiracDelta))
    assert not isinstance(x, sp.sin)


def test_lambdify_module_dict_overrides_functions(sp):
    x = sp.Symbol('x')

    class sincn(sp.Function):
        @classmethod
        def eval(cls, val):
            if val.is_zero:
                return sp.S.One

    f = sp.lambdify(x, sincn(x) + 1, [{'sincn': lambda a: 2*a}, 'math'])
    assert f(3) == 7
    g = sp.lambdify(x, sp.Function('h')(x)*2, [{'h': lambda a: a + 1}])
    assert g(3) == 8
    h = sp.lambdify(x, sp.sin(x), [{'sin': lambda a: 10*a}, 'math'])
    assert h(2) == 20


def test_array_atoms(sp):
    x, y = sp.symbols('x y')
    a = sp.Array([[x, 1], [y, 2]])
    assert a.atoms(sp.Symbol) == {x, y}
    assert a.atoms() == {x, y, sp.S(1), sp.S(2)}


# --- expression classes: Expr, AtomicExpr, Number, Basic and Singleton -------

def _lattice_expr(sp):
    # TorchInductor's MinMaxBase: Expr.__new__(cls, *args) then extra state.
    class Lattice(sp.Expr):
        def __new__(cls, *args):
            obj = sp.Expr.__new__(cls, *args)
            obj._argset = frozenset(args)
            return obj

        def _eval_derivative(self, s):
            return self.args[0].diff(s)
    return Lattice


def _pointer(sp):
    # Devito's AbstractObject / LocalObject / Pointer(LocalObject, Expr).
    class AbstractObject(sp.Basic):
        is_AbstractObject = True

        def __new__(cls, *args, **kwargs):
            obj = sp.Basic.__new__(cls)
            obj.__init__(*args, **kwargs)
            return obj

        def __init__(self, name, dtype=None):
            self.name = name
            self.dtype = dtype

        def _hashable_content(self):
            return (self.name, self.dtype)

        def __getnewargs_ex__(self):
            return ((self.name,), {'dtype': self.dtype})

        def __repr__(self):
            return self.name

        __str__ = __repr__

        free_symbols = property(lambda self: {self})

    class LocalObject(AbstractObject):
        pass

    class Pointer(LocalObject, sp.Expr):
        pass
    return Pointer


def _dimension(sp):
    # Devito's AbstractSymbol: an uncached Symbol carrying extra state.
    class Dimension(sp.Symbol):
        def __new__(cls, name, **assumptions):
            obj = sp.Symbol.__xnew__(cls, name, **assumptions)
            obj.spacing = 'h_' + name
            return obj
    return Dimension


lattice_expr = _module_level(_lattice_expr, 'Lattice')
pointer = _module_level(_pointer, 'Pointer')
dimension = _module_level(_dimension, 'Dimension')


def test_expr_subclass_is_basic_new(sp):
    x, y = sp.symbols('x y')

    class F(sp.Expr):
        pass

    class A(sp.AtomicExpr):
        pass

    f, a = F(x, sp.S(2)), A()
    assert f.args == (x, 2) and f.func is F and type(f) is F
    assert isinstance(f, sp.Expr) and isinstance(f, sp.Basic)
    assert not f.is_Function and not f.is_Atom and a.is_Atom
    assert str(f) == 'F(x, 2)' and str(a) == 'A()'
    assert sp.latex(f) == r'\operatorname{F}\left(x, 2\right)' and sp.latex(a) == r'\text{A}'
    assert sp.srepr(f) == "F(Symbol('x'), Integer(2))"
    assert f == F(x, sp.S(2)) and f != F(x, sp.S(3)) and hash(f) == hash(F(x, sp.S(2)))
    assert A() == A() and f.subs(x, y) == F(y, sp.S(2)) and type(f.subs(x, y)) is F
    assert f.xreplace({x: sp.S(3)}) == F(sp.S(3), sp.S(2)) and f.has(x) and f.atoms(sp.Symbol) == {x}
    assert f.free_symbols == {x} and a.free_symbols == set()
    assert not f.is_number and not a.is_number and F(sp.S(1)).is_number
    assert 2*f + f == 3*f and sp.sympify(f) is f and F().args == ()
    assert sp.Basic.__new__(F, x) == F(x)
    assert sorted([F(x), x, sp.sin(x), A()], key=sp.default_sort_key) == [A(), x, sp.sin(x), F(x)]
    assert copy.copy(f) == f and copy.deepcopy(f) == f


def test_expr_subclass_state_and_derivative(sp):
    x, y = sp.symbols('x y')
    Lattice = lattice_expr(sp)
    g = Lattice(x**2, y)
    assert g._argset == frozenset({x**2, y}) and type(g) is Lattice
    assert g.diff(x) == 2*x
    # subs rebuilds through the class, so the new instance has its state.
    h = g.subs(y, 1)
    assert type(h) is Lattice and h._argset == frozenset({x**2, 1})
    assert [arg for arg in (g + 1).args if type(arg) is Lattice][0]._argset == g._argset
    assert sp.Matrix([[g]])[0, 0] == g and (g*x).args == (x, g)
    restored = pickle.loads(pickle.dumps(g))
    assert restored == g and type(restored) is Lattice and restored._argset == g._argset


def test_atomic_expr_subclass_with_state(sp):
    # TorchInductor's SingletonInt.
    class SingletonInt(sp.AtomicExpr):
        _op_priority = 99999

        def __new__(cls, *args, coeff=None, **kwargs):
            instance = super().__new__(cls, *args, **kwargs)
            instance._val = args[0]
            instance._coeff = coeff if coeff is not None else 1
            return instance

        @property
        def free_symbols(self):
            return set()

        def __mul__(self, other):
            return SingletonInt(self._val, coeff=self._coeff*other)

        __rmul__ = __mul__

        def _sympystr(self, printer):
            return 'j%s' % self._val

    j = SingletonInt(1)
    assert str(j) == 'j1' and j.args == (1,) and j._val == 1 and j.is_Atom
    assert (j*2)._coeff == 2 and (2*j)._coeff == 2 and j.free_symbols == set()
    assert not j.is_number and j == SingletonInt(1)


def test_singleton_number_subclass(sp):
    from importlib import import_module
    Singleton = import_module(sp.__name__ + '.core.singleton').Singleton

    # TorchInductor's IntInfinity.
    class IntInfinity(sp.Number, metaclass=Singleton):
        is_integer = True
        is_commutative = True
        is_number = True
        is_extended_real = True
        is_extended_positive = True
        is_prime = False
        _op_priority = 100.0
        __slots__ = ()

        def __new__(cls):
            return sp.AtomicExpr.__new__(cls)

        def _sympystr(self, printer):
            return 'int_oo'

        def __add__(self, other):
            return self

        __radd__ = __add__

    oo = IntInfinity()
    assert oo is IntInfinity() and sp.S.IntInfinity is oo and type(oo) is IntInfinity
    assert str(oo) == 'int_oo' and oo + 1 is oo and 1 + oo is oo
    assert oo.is_integer and oo.is_Number and isinstance(oo, sp.Number)
    assert not isinstance(sp.Integer(2), IntInfinity)
    assert copy.copy(oo) is oo and copy.deepcopy(oo) is oo
    assert hash(oo) == hash(IntInfinity()) and oo.args == ()


def test_atomic_expr_eval_power(sp):
    # TorchInductor's int_oo**2 is int_oo.
    class Absorbing(sp.AtomicExpr):
        is_number = True
        is_extended_positive = True

        def _sympystr(self, printer):
            return 'w'

        def _eval_power(self, expt):
            if expt.is_extended_positive:
                return self
            return None

    w = Absorbing()
    x = sp.Symbol('x')
    assert w**2 is w and sp.Pow(w, 3) is w
    assert w**0 == 1 and w**1 is w
    assert str(w**-1) == '1/w' and str(w**x) == 'w**x'
    assert str(sp.Pow(w, 2, evaluate=False)) == 'w**2'


def test_number_subclass_arithmetic_folds_coefficients(sp):
    # TorchInductor's sympy.Add(2, int_oo) is int_oo: Add and Mul fold their
    # Number terms through the arithmetic of the Number class.
    class Absorbing(sp.Number):
        is_number = True
        is_extended_positive = True
        _op_priority = 100.0

        def __new__(cls):
            return sp.AtomicExpr.__new__(cls)

        def _sympystr(self, printer):
            return 'w'

        def __add__(self, other):
            other = sp.sympify(other)
            if isinstance(other, sp.Number):
                return other if other in (sp.oo, -sp.oo) else self
            return sp.Number.__add__(self, other)

        __radd__ = __add__

        def __mul__(self, other):
            other = sp.sympify(other)
            if isinstance(other, sp.Number):
                return sp.nan if other.is_zero else self
            return sp.Number.__mul__(self, other)

        __rmul__ = __mul__

        def __hash__(self):
            return hash('w')

        def __eq__(self, other):
            return isinstance(other, Absorbing)

        def __ne__(self, other):
            return not self == other

    w = Absorbing()
    x = sp.Symbol('x')
    assert sp.Add(2, w) is w and sp.Add(w, 2) is w and sp.Add(w, -x, x) is w
    assert sp.Mul(2, w) is w and sp.Mul(w, 3) is w
    assert sp.Add(w, sp.oo) == sp.oo and sp.Mul(w, 0) == sp.nan
    assert 2 + w is w and w * 2 is w


def test_singleton_registry_is_lazy(sp):
    from importlib import import_module
    Singleton = import_module(sp.__name__ + '.core.singleton').Singleton

    class LazySingleton(sp.Basic, metaclass=Singleton):
        pass

    assert sp.S._classes_to_install['LazySingleton'] is LazySingleton
    assert sp.S.LazySingleton is LazySingleton()
    assert 'LazySingleton' not in sp.S._classes_to_install


def test_basic_expr_mixin(sp):
    x = sp.Symbol('x')
    Pointer = pointer(sp)
    a, b = Pointer('ptr', dtype='int'), Pointer('ptr', dtype='int')
    assert type(a) is Pointer and str(a) == 'ptr' and a.args == ()
    assert a.name == 'ptr' and a.dtype == 'int'
    assert a == b and hash(a) == hash(b) and a != Pointer('ptr', dtype='float')
    assert {a: 1}[b] == 1 and (a*x).subs(x, 2) == 2*a
    restored = pickle.loads(pickle.dumps(a))
    assert restored == a and restored.dtype == 'int'


def test_symbol_xnew(sp):
    Dimension = dimension(sp)
    d = Dimension('d')
    assert type(d) is Dimension and d.spacing == 'h_d'
    assert d is not Dimension('d') and d == Dimension('d')
    with pytest.raises(TypeError):
        sp.Symbol.__xnew__(sp.Symbol, 1)


def test_public_classes_are_basic_subclasses(sp):
    assert issubclass(sp.Integer, sp.Basic) and issubclass(sp.Symbol, sp.Basic)
    assert issubclass(sp.Number, sp.Expr) or sp is symox
