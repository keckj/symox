"""Elementary SymPy objects completed natively: Mod/Rem, reciprocal hyperbolic
and inverse functions, sinc, frac, real_root, Id, Catalan/GoldenRatio/
TribonacciConstant, KroneckerDelta/LeviCivita/Eijk and relational aliases.

Every differential case is evaluated in a SymPy 1.14 namespace and a symox
namespace and the printed results must agree exactly.
"""

import math
import pickle

import pytest
import sympy as sp
from sympy.functions.special.tensor_functions import Eijk as sp_Eijk
import symox as sx


def namespace(lib):
    ns = dict(vars(lib))
    if lib is sp:
        ns["Eijk"] = sp_Eijk
    S = lib.Symbol
    ns.update(x=S("x"), y=S("y"), i=S("i"), j=S("j"), a=S("a"), b=S("b"),
              n=S("n", integer=True), m=S("m", integer=True),
              r=S("r", real=True), p=S("p", positive=True), q=S("q", negative=True),
              k=S("k", integer=True, nonzero=True))
    return ns


SP, SX = namespace(sp), namespace(sx)


def render(ns, case):
    try:
        value = eval(case, ns)
    except Exception as exc:  # noqa: BLE001 - error class is part of the contract
        return "ERR " + type(exc).__name__
    return value if isinstance(value, str) else str(value)


CASES = r'''
sech(0)
sech(oo)
sech(-oo)
sech(-x)
sech(I*pi)
sech(I*pi/3)
sech(asech(x))
sech(log(2))
sech(nan)
sech(zoo)
csch(0)
csch(oo)
csch(-x)
csch(I*pi/2)
csch(acsch(x))
csch(I*pi)
asec(0)
asec(1)
asec(-1)
asec(oo)
asec(zoo)
asec(2)
asec(-2)
asec(sqrt(2))
asec(2/sqrt(3))
asec(-2/sqrt(3))
asec(x)
asec(sec(x))
acsc(0)
acsc(1)
acsc(-1)
acsc(2)
acsc(-2)
acsc(sqrt(2))
acsc(oo)
acsc(zoo)
acsc(-x)
asech(0)
asech(1)
asech(-1)
asech(oo)
asech(I*oo)
asech(2)
asech(sqrt(2))
asech(2/sqrt(3))
asech(-sqrt(2))
asech(x)
acsch(0)
acsch(1)
acsch(-1)
acsch(oo)
acsch(I)
acsch(-I)
acsch(2)
acsch(-x)
acsch(zoo)
acoth(0)
acoth(1)
acoth(-1)
acoth(oo)
acoth(-oo)
acoth(I)
acoth(-x)
acoth(I*x)
acoth(zoo)
acoth(nan)
sinc(0)
sinc(pi)
sinc(2*pi)
sinc(3*pi/2)
sinc(-5*pi/2)
sinc(pi/2)
sinc(pi/3)
sinc(-x)
sinc(oo)
sinc(-oo)
sinc(nan)
sinc(zoo)
sinc(I*pi)
sinc(n*pi)
sinc(k*pi)
frac(0)
frac(Rational(7,3))
frac(-Rational(7,3))
frac(2.5)
frac(-2.5)
frac(n)
frac(zoo)
frac(nan)
frac(pi)
frac(x + n)
frac(I*Rational(7,2))
frac(3 + 2*I)
frac(x)
frac(r).is_real
frac(r).is_nonnegative
frac(r).is_finite
frac(x).is_finite
frac(n).is_zero
KroneckerDelta(1, 1)
KroneckerDelta(1, 2)
KroneckerDelta(i, j)
KroneckerDelta(j, i)
KroneckerDelta(i, i)
KroneckerDelta(i, i + 1)
KroneckerDelta(n, n + 1)
KroneckerDelta(p, 0)
KroneckerDelta(i, j).is_integer
KroneckerDelta(i, j).is_nonnegative
LeviCivita(1, 2, 3)
LeviCivita(2, 1, 3)
LeviCivita(1, 3, 2)
LeviCivita(1, 2, 2)
LeviCivita(3, 1, 2)
LeviCivita(1, 2, 3, 4)
LeviCivita(2, 1, 4, 3)
LeviCivita(i, j, k)
LeviCivita(i, j, i)
LeviCivita(1, 2, 1)
Eijk(1, 2, 3)
Mod(5, 3)
Mod(-5, 3)
Mod(5, -3)
Mod(Rational(7,2), 3)
Mod(3, 2.5)
Mod(x, 1)
Mod(n, 1)
Mod(nan, 3)
Mod(x, oo)
Mod(oo, 3)
Mod(0, x)
Mod(x, x)
Mod(-x, x)
Mod(2*n + 1, 2)
Mod(2*n, 2)
Mod(x, 2).is_nonnegative
Mod(x, -2).is_nonpositive
Mod(n, 3).is_integer
Mod(x, 3).is_integer
Mod(n, m).is_integer
Rem(5, 3)
Rem(-5, 3)
Rem(5, -3)
Rem(5.5, 2)
Rem(-7.5, 2)
Rem(x, 3)
Rem(0, x)
Rem(x, x)
Rem(n, 1)
Rem(nan, 2)
sech(x).is_real
sech(r).is_real
sech(r).is_positive
csch(p).is_positive
csch(q).is_negative
csch(r).is_real
asech(p).is_finite
acsch(r).is_real
sinc(r).is_real
sinc(r).is_finite
asec(x).is_real
diff(sech(x), x)
diff(csch(x), x)
diff(asec(x), x)
diff(acsc(x), x)
diff(asech(x), x)
diff(acsch(x), x)
diff(acoth(x), x)
diff(sinc(x), x)
diff(sinc(x**2), x)
diff(acoth(2*x), x)
N(sech(1))
N(csch(1))
N(asec(2))
N(acsc(2))
N(asech(0.5))
N(acsch(2))
N(acoth(2))
N(sinc(1))
N(asec(0.5))
N(acsc(0.5))
N(asec(-0.5))
N(acsc(-0.5))
N(asech(2))
N(asech(-0.5))
N(acsch(0.5))
N(acoth(0.5))
N(acoth(-0.5))
N(acoth(2*I))
N(asec(I))
N(acsc(2*I))
N(sech(1), 40)
N(asec(3), 40)
N(acoth(3), 40)
N(sinc(2), 40)
N(Catalan, 50)
N(GoldenRatio, 50)
N(TribonacciConstant, 50)
N(Catalan)
N(Mod(pi, 1), 30)
N(Mod(E*10, 3), 20)
N(frac(pi*10), 30)
sech(1.0)
csch(2.0)
asec(2.0)
acoth(2.0)
sinc(0.5)
frac(Float(2.75))
Mod(Float(7.5), 2)
Catalan.is_positive
Catalan.is_irrational
Catalan.is_rational
Catalan.is_algebraic
GoldenRatio.is_irrational
GoldenRatio.is_algebraic
TribonacciConstant.is_positive
TribonacciConstant.is_algebraic
GoldenRatio > 1
Catalan < 1
floor(GoldenRatio)
floor(TribonacciConstant)
floor(Catalan)
expand_func(GoldenRatio)
expand_func(TribonacciConstant)
Catalan + 1
2*GoldenRatio - 1
str(Mod(x, 2)**2)
str(2*Mod(x, 2))
str(-Mod(x, y))
str(x/Mod(y, 3))
str(-x/Mod(y, 3))
str(Mod(x, 3)/2)
str(2**Mod(x, 2))
str(Mod(x, 2)*Mod(y, 3))
str(1/Mod(x, 2))
srepr(Mod(x, 2))
srepr(sech(x))
srepr(KroneckerDelta(i, j))
srepr(Catalan)
srepr(GoldenRatio + TribonacciConstant)
latex(Mod(x, 3))
latex(2*Mod(x, 3))
latex(Mod(x, 3)**2)
latex(Mod(x + 1, 2*y))
latex(-Mod(x, y))
latex(Mod(x, 3) + 1)
latex(y*Mod(x, 3)**2)
latex(sech(x))
latex(csch(x)**2)
latex(asec(x))
latex(acoth(x)**2)
latex(sinc(x))
latex(frac(x))
latex(frac(x)**2)
latex(KroneckerDelta(i, j))
latex(KroneckerDelta(i + 1, j))
latex(LeviCivita(i, j, k))
latex(LeviCivita(i + 1, j, k))
latex(Rem(x, 3))
latex(Catalan + GoldenRatio + TribonacciConstant)
latex(asec(x), inv_trig_style='full')
latex(asech(x), inv_trig_style='full')
latex(acoth(x), inv_trig_style='power')
pretty(Mod(x, 3))
pretty(2*Mod(x, 3))
pretty(Mod(x, 3)**2)
pretty(Mod(x + 1, 2*y))
pretty(-Mod(x, y))
pretty(Mod(x, 3) + 1)
pretty(KroneckerDelta(i, j))
pretty(KroneckerDelta(i + 1, j))
pretty(KroneckerDelta(i, j), use_unicode=False)
pretty(sech(x))
pretty(Catalan + GoldenRatio + TribonacciConstant)
pretty(GoldenRatio, use_unicode=False)
sympify('Mod(x, 3)')
sympify('sech(x) + acoth(y)')
sympify('KroneckerDelta(i, j)')
sympify('Catalan + GoldenRatio + TribonacciConstant')
sympify('frac(x) + sinc(x)')
sympify('Rem(x, 3)')
sympify('LeviCivita(1, 2, 3)')
x % 3
7 % Integer(3)
Integer(7) % 3
x // 3
divmod(x, 3)
divmod(Integer(7), 3)
Integer(7) // 2.5
Mod(x, 3).subs(x, 7)
sech(x).subs(x, 0)
KroneckerDelta(i, j).subs({i: 1, j: 1})
Mod(x, 3).args
Mod(x, 3).func(y, 2)
KroneckerDelta(1, 1)**2
KroneckerDelta(i, j)**2
KroneckerDelta(i, j)**3
KroneckerDelta(i, j)*KroneckerDelta(i, j)
Mod(x, 3).free_symbols
real_root(-8, 3)
real_root(x, 3)
real_root(-8)
real_root(-32, 5)
real_root(16, 4)
Id(x)
'''.strip().splitlines()


@pytest.mark.parametrize("case", CASES)
def test_matches_sympy(case):
    assert render(SX, case) == render(SP, case)


MOD_CASES = [
    '(x**2 - 1, x + 1)',
    '(2*x, 4)',
    '(6*x, 4*y)',
    '(x*y, x)',
    '(n*m, n)',
    '(0.6*x, 0.3*y)',
    '(x + y + 2, x)',
    '(-x, -y)',
    '(-x, 3)',
    '(Mod(x, 6), 3)',
    '(Mod(x, 2), 4)',
    '(3*n*Mod(m, 4), 4)',
    '(p - 1, p)',
    '(2*n + 1, 2)',
    '(pi, 3)',
    '(I, 3)',
    '(Integer(3)**100, 7)',
    '(factorial(10), 11)',
    '(x**2, x)',
    '(2.5, 1)',
    '(5, -3)',
    '(-5, 3)',
    '(Rational(7, 2), 3)',
    '(x, oo)',
    '(2*x + 4, 2)',
    '(4*x + 2, 4)',
    '(x*1.0, 3)',
    '(-2*x, -4)',
    '(x - 3, 3)',
    '(n**2, n)',
    '(4*x, 2*x)',
    '(x, 2*x)',
    '(n*m*6, 4)',
    '(Mod(x, 3), -3)',
    '(-Mod(x, 3), 3)',
    '(-Mod(x, 3), -3)',
    '(-pi, 3)',
    '(Rational(-7, 3), Rational(2, 3))',
    '(6.0*x, 4)',
    '(x, 2.5)',
    '(2.5*x, 0.5)',
    '(x**2 + 2*x, x)',
    '(3, 2.5)',
    '(-3, 2.5)',
    '(7.5, -2)',
    '(Rational(1, 3), 0.1)',
    '(-x + y, -3)',
    '(n*Mod(n, 5), 3)',
    '(3.0, 1.5)',
    '(factorial(n), n)',
    '(factorial(n), n + 1)',
]


@pytest.mark.parametrize("args", MOD_CASES)
def test_mod_canonicalisation_table(args):
    ns_sp = dict(SP, p=sp.Symbol("p", positive=True, integer=True))
    ns_sx = dict(SX, p=sx.Symbol("p", positive=True, integer=True))
    case = "Mod" + args
    assert render(ns_sx, case) == render(ns_sp, case)


KD_SUMS = [
    "Sum(KroneckerDelta(i, j), (i, 0, n))",
    "Sum(KroneckerDelta(i, j)*x**i, (i, 0, oo))",
    "Sum(KroneckerDelta(i, 2), (i, 0, 5))",
    "Sum(KroneckerDelta(i, 7), (i, 0, 5))",
    "Sum(KroneckerDelta(i, j)*i**2, (i, a, b))",
    "Sum(KroneckerDelta(2*i, j), (i, 0, n))",
    "Sum(KroneckerDelta(i + 1, j)*x, (i, 1, n))",
    "summation(KroneckerDelta(i, j), (i, 1, 3))",
    "Sum(KroneckerDelta(i, j)*KroneckerDelta(i, a), (i, 0, n))",
]


@pytest.mark.parametrize("case", KD_SUMS)
def test_kronecker_delta_summation(case):
    i, j, n = (sp.Symbol(s) for s in "ijn")
    ns_sp = dict(SP, n=n)
    ns_sx = dict(SX, n=sx.Symbol("n"))
    assert render(ns_sx, case + ".doit()") == render(ns_sp, case + ".doit()")


def test_division_by_zero():
    x = sx.Symbol("x")
    with pytest.raises(ZeroDivisionError):
        sx.Mod(x, 0)
    with pytest.raises(ZeroDivisionError):
        sx.Rem(x, 0)
    with pytest.raises(ZeroDivisionError):
        x % 0
    with pytest.raises(ZeroDivisionError):
        sx.Integer(3) // 0
    assert sx.Mod(0, x) == 0


def test_operators():
    x = sx.Symbol("x")
    assert divmod(sx.Integer(7), 3) == (2, 1)
    assert divmod(x, 3) == (sx.floor(x / 3), sx.Mod(x, 3))
    assert 7 % sx.Integer(3) == 1 and x % 3 == sx.Mod(x, 3)
    assert 7 // sx.Integer(2) == 3 and str(x // 3) == "floor(x/3)"


OBJECTS = [
    "Mod(x, 3)", "Rem(x, 3)", "sech(x)", "csch(x)", "asec(x)", "acsc(x)",
    "asech(x)", "acsch(x)", "acoth(x)", "sinc(x)", "frac(x)",
    "KroneckerDelta(i, j)", "LeviCivita(i, j, a)", "Catalan", "GoldenRatio",
    "TribonacciConstant", "2*Mod(x, 3) + sech(x)**2 + Catalan*x",
]


@pytest.mark.parametrize("case", OBJECTS)
def test_conversion_and_pickle(case):
    e = eval(case, SX)
    s = sx.to_sympy(e)
    assert s == eval(case, SP)
    assert sx.from_sympy(s) == e
    assert pickle.loads(pickle.dumps(e)) == e
    assert sx.sympify(str(e)) == e
    assert sx.srepr(e) == sp.srepr(s)


@pytest.mark.parametrize("case,value", [
    ("Mod(x, 3)", 7.5), ("Mod(x, 3)", -7.5), ("Rem(x, 3)", -7.5),
    ("sech(x)", 1.0), ("csch(x)", 1.0), ("asec(x)", 2.0), ("acsc(x)", 2.0),
    ("asech(x)", 0.5), ("acsch(x)", 2.0), ("acoth(x)", 2.0), ("sinc(x)", 0.5),
    ("sinc(x)", 0.0), ("frac(x)", -2.25), ("Catalan*x", 1.0),
    ("GoldenRatio*x", 1.0), ("TribonacciConstant*x", 1.0),
])
def test_lambdify(case, value):
    got = sx.lambdify(SX["x"], eval(case, SX))(value)
    want = complex(sp.N(eval(case, SP).subs(SP["x"], value), 30))
    assert abs(got - want) < 1e-12


CODE = [
    ("ccode", "Mod(x, 3)"), ("ccode", "Mod(n, 3)"), ("ccode", "sech(x)"),
    ("ccode", "acoth(x)"), ("ccode", "sinc(x)"), ("ccode", "frac(x)"),
    ("ccode", "GoldenRatio*x"), ("ccode", "2*Mod(x, 3) + 1"),
    ("pycode", "Mod(x, 3)"), ("pycode", "sech(x)"), ("pycode", "acoth(x)"),
    ("pycode", "sinc(x)"), ("pycode", "frac(x)"), ("pycode", "Mod(x + 1, y)**2"),
    ("pycode", "KroneckerDelta(x, y)"), ("cxxcode", "Mod(x, 3)"),
    ("rust_code", "sech(x)"), ("rust_code", "frac(x)"),
]


@pytest.mark.parametrize("printer,case", CODE)
def test_code_printers(printer, case):
    assert getattr(sx, printer)(eval(case, SX)) == getattr(sp, printer)(eval(case, SP))


def test_constants_evalf():
    for name in ("Catalan", "GoldenRatio", "TribonacciConstant"):
        assert str(sx.N(SX[name], 60)) == str(sp.N(SP[name], 60))
        assert float(SX[name]) == float(SP[name])
    assert math.isclose(float(sx.GoldenRatio), (1 + math.sqrt(5)) / 2)


def test_relational_aliases_and_atoms():
    x, y = sx.symbols("x y")
    pairs = [(sx.Equality, sx.Eq(x, y)), (sx.Unequality, sx.Ne(x, y)),
             (sx.GreaterThan, x >= y), (sx.LessThan, x <= y),
             (sx.StrictGreaterThan, x > y), (sx.StrictLessThan, x < y)]
    for cls, rel in pairs:
        assert isinstance(rel, cls) and isinstance(rel, sx.Relational)
        assert cls(x, y) == rel
    assert not isinstance(x < y, sx.GreaterThan)
    assert sx.Rel(x, y, "<=") == (x <= y) and sx.Rel(x, y, "ne") == sx.Ne(x, y)
    assert sx.Rel(x, y) == sx.Eq(x, y)
    assert isinstance(x, sx.Atom) and isinstance(sx.Integer(2), sx.AtomicExpr)
    assert not isinstance(x + 1, sx.Atom)
    assert isinstance(sx.true, sx.Atom) and not isinstance(sx.true, sx.AtomicExpr)


def test_real_root_id_eijk():
    x = sx.Symbol("x")
    assert sx.real_root(-8, 3) == -2 and sx.real_root(-32, 5) == -2
    assert sx.real_root(sx.Integer(-8) ** sx.Rational(1, 3)) == -2
    assert sx.Id(x) == x and sx.Id(sx.Integer(3)) == 3
    assert sx.Eijk is sx.LeviCivita and sx.Eijk(2, 1, 3) == -1


def test_assumption_edges():
    x = sx.Symbol("x")
    n = sx.Symbol("n", integer=True)
    p = sx.Symbol("p", positive=True, integer=True)
    assert sx.Mod(p - 1, p) == p - 1
    assert sx.Mod(x, 2).is_nonnegative and sx.Mod(x, -2).is_nonpositive
    assert sx.frac(n) == 0 and sx.sinc(sx.pi * sx.Symbol("k", integer=True, nonzero=True)) == 0
    assert sx.KroneckerDelta(sx.Symbol("i"), sx.Symbol("j")) ** 5 == sx.KroneckerDelta(sx.Symbol("i"), sx.Symbol("j"))
