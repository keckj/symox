"""Native algorithms dispatch user hooks of Expr, Function and Symbol subclasses.

Each case builds its classes against one module and returns a value; symox
must return exactly what SymPy 1.14 returns, including the exceptions that a
hook raises. Basic results are compared through srepr, so argument order and
structure count, not only printing.
"""
import pytest
import sympy
import symox


def _outcome(case, sp):
    try:
        value = case(sp)
    except Exception as exc:
        return ('raises', type(exc).__name__, str(exc))
    return _plain(value, sp)


def _plain(value, sp):
    if isinstance(value, sp.Basic):
        return sp.srepr(value)
    if isinstance(value, (list, tuple)):
        return type(value)(_plain(item, sp) for item in value)
    if isinstance(value, dict):
        return sorted((_plain(key, sp), _plain(item, sp)) for key, item in value.items())
    if isinstance(value, (set, frozenset)):
        return sorted(_plain(item, sp) for item in value)
    return value


CASES = {}


def case(fn):
    CASES[fn.__name__] = fn
    return fn


# _eval_subs

def _subs_function(sp):
    x = sp.Symbol('x')

    class F1(sp.Function):
        def _eval_subs(self, old, new):
            if old == x:
                return sp.Integer(42)
    return F1


@case
def eval_subs_function(sp):
    x = sp.Symbol('x')
    return _subs_function(sp)(x).subs(x, 1)


@case
def eval_subs_function_nested(sp):
    x, y = sp.symbols('x y')
    return (_subs_function(sp)(x) + y).subs(x, 1)


@case
def eval_subs_expr_without_pattern(sp):
    x, y = sp.symbols('x y')

    class E1(sp.Expr):
        def _eval_subs(self, old, new):
            if old == y:
                return sp.Integer(7)
    return (E1(x) * 2 + y).subs(y, 3)


@case
def eval_subs_atom(sp):
    y = sp.Symbol('y')

    class A1(sp.Symbol):
        def _eval_subs(self, old, new):
            if old == y:
                return sp.Integer(5)
    return (A1('a1') + 1).subs(y, 3)


# _eval_is_* with deduction through the assumption rules and Add/Mul

def _positive(sp):
    class P1(sp.Function):
        def _eval_is_positive(self):
            return True

        def _eval_is_real(self):
            return True
    return P1


@case
def is_positive_direct(sp):
    return _positive(sp)(sp.Symbol('x')).is_positive


@case
def is_nonnegative_deduced(sp):
    return _positive(sp)(sp.Symbol('x')).is_nonnegative


@case
def is_zero_deduced(sp):
    return _positive(sp)(sp.Symbol('x')).is_zero


@case
def add_positive(sp):
    return (_positive(sp)(sp.Symbol('x')) + 1).is_positive


@case
def mul_positive(sp):
    return (_positive(sp)(sp.Symbol('x')) * 2).is_positive


@case
def sqrt_of_square(sp):
    return sp.sqrt(_positive(sp)(sp.Symbol('x')) ** 2)


@case
def abs_of_positive(sp):
    return sp.Abs(_positive(sp)(sp.Symbol('x')))


@case
def expr_positive(sp):
    class PE(sp.Expr):
        is_commutative = True

        def _eval_is_positive(self):
            return True
    value = PE(sp.Symbol('x'))
    return value.is_positive, (value + 2).is_positive


def test_symbol_handler_deduces_facts():
    # In SymPy 1.14, Symbols with equal assumptions share one fact base, so
    # any earlier query of a plain Symbol('n') hides the subclass handler;
    # these are the answers of a fresh SymPy process.
    class PS(symox.Symbol):
        def _eval_is_integer(self):
            return True
    n = PS('n')
    assert n.is_integer is True
    assert PS('m').is_rational is True
    assert symox.cos(symox.pi * n) == (-1) ** n


# fdiff and _eval_derivative

@case
def fdiff(sp):
    x = sp.Symbol('x')

    class D1(sp.Function):
        def fdiff(self, argindex=1):
            return sp.cos(self.args[0])
    return sp.diff(D1(x ** 2), x)


@case
def function_eval_derivative(sp):
    x = sp.Symbol('x')

    class D2(sp.Function):
        def _eval_derivative(self, s):
            return sp.Integer(99)
    return sp.diff(D2(x), x), sp.diff(D2(x) * x, x)


@case
def expr_eval_derivative(sp):
    x = sp.Symbol('x')

    class D3(sp.Expr):
        def _eval_derivative(self, s):
            return sp.Integer(11)
    return sp.diff(D3(x) + x ** 2, x)


# _eval_evalf

@case
def eval_evalf(sp):
    x = sp.Symbol('x')

    class V1(sp.Function):
        def _eval_evalf(self, prec):
            return sp.Float(1.5)
    return str(V1(x).evalf()), str((V1(x) + 1).evalf()), str(sp.N(2 * V1(x)))


@case
def eval_evalf_numeric_argument(sp):
    class V2(sp.Function):
        def _eval_evalf(self, prec):
            return sp.Float(2.25)
    return str(V2(1).evalf()), str((V2(1) * 2 + sp.pi).evalf(5))


def _mpf_numbers(sp):
    from mpmath.libmp import finf, from_rational

    class Half(sp.Number):
        is_comparable = True
        is_extended_real = True

        def __new__(cls):
            return sp.Expr.__new__(cls)

        def _as_mpf_val(self, prec):
            return from_rational(1, 2, prec)

        def __hash__(self):
            return 17

    class Huge(sp.Number):
        is_comparable = True
        is_extended_real = True

        def __new__(cls):
            return sp.Expr.__new__(cls)

        def _as_mpf_val(self, prec):
            return finf

        def __hash__(self):
            return 19
    return Half(), Huge()


@case
def number_evalf_from_mpf_value(sp):
    # Number._eval_evalf: Float._new(self._as_mpf_val(prec), prec).
    half, huge = _mpf_numbers(sp)
    # (SymPy's Number arithmetic recurses for classes without it.)
    return half.evalf(2), sp.N(half, 30), huge.evalf(2), sp.N(huge), float(half), float(huge)


@case
def number_evalf_raising_mpf_value(sp):
    class Broken(sp.Number):
        def __new__(cls):
            return sp.Expr.__new__(cls)

        def _as_mpf_val(self, prec):
            raise ArithmeticError('no value')

        def __hash__(self):
            return 23
    return Broken().evalf()


# _eval_rewrite_as_*, _eval_expand_*, _eval_simplify

@case
def eval_rewrite(sp):
    x = sp.Symbol('x')

    class R1(sp.Function):
        def _eval_rewrite_as_sin(self, *args, **kwargs):
            return sp.sin(args[0])
    return R1(x).rewrite(sp.sin), (R1(x) + 1).rewrite(sp.sin)


@case
def eval_expand_func(sp):
    x = sp.Symbol('x')

    class X1(sp.Function):
        def _eval_expand_func(self, **hints):
            return self.args[0] + 1
    return sp.expand_func(X1(x)), sp.expand_func(X1(x) * 2)


@case
def eval_expand_mul(sp):
    x = sp.Symbol('x')

    class X2(sp.Function):
        def _eval_expand_mul(self, **hints):
            return sp.Integer(3)
    return (X2(x) * (x + 1)).expand()


@case
def eval_expand_custom_hint(sp):
    x = sp.Symbol('x')

    class X3(sp.Function):
        def _eval_expand_custom(self, **hints):
            return sp.Integer(17)
    return (X3(x) + 1).expand(custom=True)


@case
def eval_simplify(sp):
    x = sp.Symbol('x')

    class SM(sp.Function):
        def _eval_simplify(self, **kwargs):
            return sp.Integer(5)
    return sp.simplify(SM(x)), sp.simplify(SM(x) + sp.sin(x) ** 2 + sp.cos(x) ** 2)


# func, free_symbols, _eval_power and _op_priority overrides

@case
def func_override(sp):
    x, y = sp.symbols('x y')

    class FO(sp.Expr):
        @property
        def func(self):
            return lambda *args: sp.Integer(123)
    return (FO(x) + 1).xreplace({x: y}), (FO(x) + 1).subs(x, y)


@case
def free_symbols_override(sp):
    x, y, z = sp.symbols('x y z')

    class FS(sp.Function):
        @property
        def free_symbols(self):
            return {z}
    return FS(x).free_symbols, (FS(x) + y).free_symbols


@case
def eval_power(sp):
    x = sp.Symbol('x')

    class PW(sp.Expr):
        def _eval_power(self, e):
            return sp.Integer(8)
    return PW(x) ** 3, sp.Pow(PW(x), 2)


@case
def op_priority(sp):
    x, y = sp.symbols('x y')

    class OP(sp.Expr):
        _op_priority = 20.0

        def __rmul__(self, other):
            return 'rmul'

        def __radd__(self, other):
            return 'radd'

        def __mul__(self, other):
            return 'mul'
    return x * OP(y), x + OP(y), sp.Integer(2) * OP(y)


# Exceptions raised by hooks propagate

@case
def raising_eval(sp):
    x = sp.Symbol('x')

    class EX(sp.Function):
        @classmethod
        def eval(cls, a):
            if a == 1:
                raise ValueError('boom eval')
    return EX(x).subs(x, 1)


@case
def raising_fdiff(sp):
    x = sp.Symbol('x')

    class EX2(sp.Function):
        def fdiff(self, argindex=1):
            raise ValueError('boom fdiff')
    return sp.diff(EX2(x), x)


@case
def raising_printer(sp):
    x = sp.Symbol('x')

    class EX3(sp.Function):
        def _sympystr(self, printer):
            raise ValueError('boom str')
    return str(EX3(x) + 1)


@case
def raising_assumption(sp):
    x = sp.Symbol('x')

    class EX4(sp.Function):
        def _eval_is_positive(self):
            raise ValueError('boom assume')
    return (EX4(x) + 1).is_positive


@case
def raising_assumption_is_not_remembered(sp):
    x = sp.Symbol('x')
    calls = []

    class EX5(sp.Function):
        def _eval_is_positive(self):
            calls.append(1)
            if len(calls) == 1:
                raise ValueError('first call')
            return True

        def _eval_is_real(self):
            return True
    value = EX5(x) + 1
    try:
        value.is_positive
    except ValueError:
        pass
    return value.is_positive


# _hashable_content

@case
def hashable_content(sp):
    class HC(sp.Symbol):
        def __new__(cls, name, tag):
            obj = sp.Symbol.__new__(cls, name)
            obj.tag = tag
            return obj

        def _hashable_content(self):
            return super()._hashable_content() + (self.tag,)
    h1, h2, h3 = HC('h', 1), HC('h', 2), HC('h', 1)
    return (h1 == h2, h1 == h3, hash(h1) == hash(h2), hash(h1) == hash(h3),
            h1 + h2, h1 + h3)


@case
def hashable_content_extends_base(sp):
    x, y = sp.symbols('x y')

    class S(sp.Symbol):
        def __new__(cls, name, tag):
            obj = sp.Symbol.__xnew__(cls, name)
            obj.tag = tag
            return obj

        def __getnewargs__(self):
            return (self.name, self.tag)

        def _hashable_content(self):
            return super()._hashable_content() + (self.tag,)
    s1, s2, s3 = S('s', 1), S('s', 2), S('s', 1)
    return (s1 == s2, s1 == s3, hash(s1) == hash(s2), hash(s1) == hash(s3),
            s1 + s2, s1 + s3, s1 + y == s2 + y, s1 * y == s3 * y,
            len({s1 + 1, s2 + 1, s3 + 1}), s1._hashable_content())


@case
def raising_hashable_content(sp):
    class S(sp.Symbol):
        def __new__(cls, name):
            return sp.Symbol.__xnew__(cls, name)

        def _hashable_content(self):
            raise ValueError('boom content')
    a, b = S('s'), S('s')
    outcomes = []
    for operation in (lambda: a == b, lambda: hash(a), lambda: a + b, lambda: a == sp.Symbol('s')):
        try:
            outcomes.append(_plain(operation(), sp))
        except ValueError as exc:
            outcomes.append(str(exc))
    return outcomes


@case
def base_hashable_content(sp):
    x = sp.Symbol('x', positive=True)
    f = sp.Function('f')
    return [value._hashable_content() for value in
            (x, sp.Symbol('y'), sp.Float('1.5', 30), sp.Rational(2, 3), sp.Integer(3),
             x + 1, sp.sin(x), f(x), sp.pi, sp.Wild('w', exclude=[x]))]


@case
def args_order_uses_host_class_names(sp):
    # Basic.compare ranks type(self).__name__ through ordering_of_classes.
    x, y = sp.symbols('x y')

    class C(sp.Symbol):
        pass

    class A(sp.Symbol):
        pass

    class Zed(sp.Symbol):
        pass

    class Symbol(sp.Dummy):
        pass

    class F(sp.Function):
        pass
    c, b, a, d = C('c'), A('b'), Zed('a'), Symbol('d')
    return ([[str(arg) for arg in e.args] for e in
             (x + c, b + c + x, a + c + b + x + y, d + x + c, x * c * b * a, F(x) + c + sp.sin(x) + x)],
            sp.Add._from_args([x, c]) == x + c)


@case
def wild_subclass(sp):
    x, y = sp.symbols('x y')

    class Tag(sp.Wild):
        pass
    t = Tag('t', exclude=[x])
    found = (2 * y).match(2 * t)
    return (type(t).__name__, isinstance(t, sp.Wild), type(next(iter(found))).__name__,
            found, (2 * x).match(2 * t), t + 1, (t + 1).subs(t, 5),
            (sp.sin(x) + 2).replace(sp.sin(Tag('u')), Tag('u') ** 2),
            t == Tag('t', exclude=[x]), t == sp.Wild('t', exclude=[x]),
            t.__getnewargs__(), t.__getnewargs_ex__())


@case
def symbol_subclass_newargs(sp):
    class S(sp.Symbol):
        pass

    class D(sp.Dummy):
        pass
    return S('s').__getnewargs__(), S('s').__getnewargs_ex__(), D('d').__getnewargs__()


# is_commutative: None and False, as Mul.flatten treats them

def _noncommuting(sp):
    class CN(sp.Expr):
        is_commutative = None

    class F(sp.Function):
        pass

    class G(sp.Function):
        is_commutative = False
    return CN, F, G


@case
def commutativity_defaults(sp):
    x = sp.Symbol('x')
    CN, F, G = _noncommuting(sp)

    class CF(sp.Function):
        is_commutative = None

    class E(sp.Expr):
        pass

    class S1(sp.Symbol):
        pass

    class S2(sp.Symbol):
        is_commutative = None
    return [value.is_commutative for value in
            (CN(x), CF(x), G(x), F(x), E(x), S1('s'), S2('t'), F(CN(x)))]


NONCOMMUTING = {
    'add': lambda a, b, x, y, F, G, sp: a + 1,
    'mul_number': lambda a, b, x, y, F, G, sp: a * 2,
    'pow': lambda a, b, x, y, F, G, sp: a ** 2,
    'sin': lambda a, b, x, y, F, G, sp: sp.sin(a),
    'function': lambda a, b, x, y, F, G, sp: F(a),
    'ab': lambda a, b, x, y, F, G, sp: a * b,
    'ba': lambda a, b, x, y, F, G, sp: b * a,
    'sum': lambda a, b, x, y, F, G, sp: b + a,
    'square_product': lambda a, b, x, y, F, G, sp: (a * b) ** 2,
    'combine_adjacent': lambda a, b, x, y, F, G, sp: sp.Pow(a, 2) * a,
    'aba': lambda a, b, x, y, F, G, sp: a * b * a,
    'exponents': lambda a, b, x, y, F, G, sp: a ** 2 * a ** 3,
    'symbolic_exponents': lambda a, b, x, y, F, G, sp: a ** x * a ** y,
    'coefficients': lambda a, b, x, y, F, G, sp: 2 * a * 3 * b,
    'commutative_part': lambda a, b, x, y, F, G, sp: b * x * a * y,
    'product': lambda a, b, x, y, F, G, sp: (a + b) * (a - b),
    'expand_product': lambda a, b, x, y, F, G, sp: ((a + b) * (a - b)).expand(),
    'expand_power': lambda a, b, x, y, F, G, sp: ((a + b) ** 2).expand(),
    'expand_mul': lambda a, b, x, y, F, G, sp: sp.expand(b * (a + x)),
    'subs': lambda a, b, x, y, F, G, sp: (a * b).subs(x, y),
    'subs_factor': lambda a, b, x, y, F, G, sp: (b * a).subs(b, x),
    'xreplace': lambda a, b, x, y, F, G, sp: (b * a).xreplace({a: b}),
    'commutator': lambda a, b, x, y, F, G, sp: a * b - b * a,
    'rational_power': lambda a, b, x, y, F, G, sp: (4 * a * b) ** sp.Rational(1, 2),
    'rational_power_positive': lambda a, b, x, y, F, G, sp:
        (4 * sp.Symbol('p', positive=True) * a * b) ** sp.Rational(1, 3),
    'symbolic_power': lambda a, b, x, y, F, G, sp: (x * a * b) ** y,
    'inverse': lambda a, b, x, y, F, G, sp: ((a * b) ** -1, (2 * a * b) ** -1),
    'power_times_base': lambda a, b, x, y, F, G, sp: (a * b) ** 2 * (a * b),
    'abab': lambda a, b, x, y, F, G, sp: (a * b) * (a * b),
    'base_times_power': lambda a, b, x, y, F, G, sp: a * (a * b) ** 2,
    'negative_cube': lambda a, b, x, y, F, G, sp: (-a * b) ** 3,
    'cancel_adjacent': lambda a, b, x, y, F, G, sp: a * a ** -1,
    'no_cancel_across': lambda a, b, x, y, F, G, sp: b * a / b,
    'division': lambda a, b, x, y, F, G, sp: b / a,
    'noncommutative_function': lambda a, b, x, y, F, G, sp: (G(y) * G(x) * G(y), G(x) * x * G(y) * 2),
    'zero': lambda a, b, x, y, F, G, sp: 0 * a * b,
    'equality': lambda a, b, x, y, F, G, sp: (a * b == b * a, a * b == a * b),
    'args_cnc': lambda a, b, x, y, F, G, sp: (x * b * a).args_cnc(),
    'collect_like_terms': lambda a, b, x, y, F, G, sp: a * b + 2 * a * b,
    'term_sum': lambda a, b, x, y, F, G, sp: 2 * (a * b + b * a),
    'free_symbols': lambda a, b, x, y, F, G, sp: (b * a).free_symbols,
    'product_commutativity': lambda a, b, x, y, F, G, sp: (a * b).is_commutative,
    'mixed': lambda a, b, x, y, F, G, sp: F(x) * a * F(y) * b,
    'repeated': lambda a, b, x, y, F, G, sp: a * a * b * b,
    'exp': lambda a, b, x, y, F, G, sp: (sp.exp(a) * sp.exp(a), sp.exp(a) * sp.exp(b)),
    'derivative': lambda a, b, x, y, F, G, sp: sp.diff(F(x) * type(a)(F(x)), x),
    'derivative_order': lambda a, b, x, y, F, G, sp: sp.diff(F(x) * G(x) * F(x) ** 2, x),
    'hash': lambda a, b, x, y, F, G, sp: hash(a * b) == hash(b * a),
}


def _noncommuting_case(name):
    def run(sp):
        x, y = sp.symbols('x y')
        CN, F, G = _noncommuting(sp)
        return NONCOMMUTING[name](CN(x), CN(y), x, y, F, G, sp)
    return run


for _name in NONCOMMUTING:
    CASES['noncommuting_' + _name] = _noncommuting_case(_name)


@pytest.mark.parametrize('name', sorted(CASES))
def test_hook_matches_sympy(name):
    assert _outcome(CASES[name], symox) == _outcome(CASES[name], sympy)


def test_raised_hook_errors_are_the_original_exceptions():
    x = symox.Symbol('x')

    class Failure(Exception):
        pass

    class EX(symox.Function):
        def fdiff(self, argindex=1):
            raise Failure('kept')
    with pytest.raises(Failure, match='kept'):
        symox.diff(EX(x) * x, x)
    # The failure does not leak into the next native computation.
    assert symox.diff(x ** 2, x) == 2 * x


def test_noncommuting_order_is_kept_in_products():
    x, y = symox.symbols('x y')

    class CN(symox.Expr):
        is_commutative = None
    a, b = CN(x), CN(y)
    assert (b * a).args == (b, a)
    assert (b * a) != (a * b)
    assert (x * b * a * y).args == (x, y, b, a)
