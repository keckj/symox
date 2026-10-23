"""Parity with SymPy 1.14 for a set of small gaps: multi-symbol
reduce_inequalities, the polys ``wrt`` option, integrate of vectors, exact
lambdify on integers, Add/Mul of matrix symbols, pretty/latex Sum and Limit
and the ANP printing of algebraic-field ring elements.

Every expected value was computed with the reference SymPy.
"""
import numpy as np
import pytest

from symox import *
del test  # the star import's suite runner must not be collected by pytest
from symox.vector import CoordSys3D
from symox.polys.rings import ring

x, y, z, t = symbols('x y z t')
f = x**2*y + x*y**2 + x*y
N = CoordSys3D('N')
v = x*N.i + x**2*N.j + sin(x)*N.k
A = MatrixSymbol('A', 2, 2); B = MatrixSymbol('B', 2, 2); C = MatrixSymbol('C', 2, 3)


def MD(e):
    return (e.func.__name__, isinstance(e, MatAdd), isinstance(e, MatMul), str(e),
            e.shape if hasattr(e, 'shape') else None)


K = QQ.algebraic_field(sqrt(2))
R, X = ring("X", K)
s2 = K.from_sympy(sqrt(2))
p = X**2 + X*s2 + 1
K3 = QQ.algebraic_field(sqrt(2), sqrt(3))
R3, Y, Z = ring("Y,Z", K3)
q = Y*Z*K3.from_sympy(sqrt(3)) + 2

CASES = {
    'ri1': lambda: reduce_inequalities(x*y > 1, [x]),
    'ri2': lambda: reduce_inequalities(x*y < 1, [x]),
    'ri3': lambda: reduce_inequalities(x*y >= 1, [y]),
    'ri4': lambda: reduce_inequalities(x + y > 1, [x]),
    'ri5': lambda: reduce_inequalities(x*y**2 > 1, [x]),
    'ri6': lambda: reduce_inequalities([x*y > 1, x < 3], [x]),
    'ri7': lambda: reduce_inequalities(x*y > 1, x),
    'ri9': lambda: reduce_inequalities(Symbol('p',positive=True)*x > 1, [x]),
    'fl1': lambda: factor_list(f, wrt=x),
    'fl2': lambda: factor_list(f, wrt=y),
    'fl3': lambda: factor_list(x**2*y**2 - 1, wrt=x),
    'sq1': lambda: sqf_list((x+y)**2*(x-1), wrt=x),
    'sq2': lambda: sqf_list(x**2*y**2 + 2*x*y + 1, wrt=y),
    'dv1': lambda: div(x**2 + y, x + y, wrt=x),
    'dv2': lambda: div(x**2*y + 1, x*y + 1, wrt=y),
    'fa1': lambda: factor(f, wrt=x),
    'sq3': lambda: sqf((x+y)**2*(x-1), wrt=x),
    'rem1': lambda: rem(x**2 + y, x + y, wrt=x),
    'quo1': lambda: quo(x**2 + y, x + y, wrt=x),
    'gcd1': lambda: gcd(x**2 - y**2, x - y, wrt=x),
    'pd1': lambda: pdiv(x**2 + y, x + y, wrt=x),
    'ter1': lambda: terms_gcd(x*y + x**2*y, wrt=x),
    'vi1': lambda: integrate(v, x),
    'vi2': lambda: integrate(v, (x, 0, 1)),
    'vi3': lambda: integrate(N.x*N.i, N.x),
    'vi4': lambda: Integral(v, x).doit(),
    'vi5': lambda: integrate(N.i*0 + 3*N.j, x),
    'vi6': lambda: v.integrate(x),
    'la1': lambda: (lambda r: (r.dtype.name, r.tolist()))(lambdify(x, x + 1)(np.array([1,2,3]))),
    'la2': lambda: (lambda r: (type(r).__name__, r))(lambdify(x, x + 1)(2**60)),
    'la3': lambda: (lambda r: (r.dtype.name, r.tolist()))(lambdify(x, 2*x**2)(np.array([1,2,3]))),
    'la4': lambda: (lambda r: (r.dtype.name, r.tolist()))(lambdify(x, x/2)(np.array([1,2,3]))),
    'la5': lambda: (lambda r: (type(r).__name__, r))(lambdify(x, x**2 + 3*x)(2**60)),
    'la6': lambda: (lambda r: (type(r).__name__, r))(lambdify((x,y), x*y - 1)(3, 4)),
    'la7': lambda: (lambda r: (r.dtype.name, r.tolist()))(lambdify((x,y), x*y - 1)(np.array([3]), np.array([4]))),
    'la8': lambda: (lambda r: (type(r).__name__, r))(lambdify(x, x + 1, 'math')(2**60)),
    'la9': lambda: (lambda r: (type(r).__name__, r))(lambdify(x, x/3)(2)),
    'la10': lambda: (lambda r: (type(r).__name__, r))(lambdify(x, Rational(1,2)*x)(2**60)),
    'ma1': lambda: MD(Add(A, B)),
    'ma2': lambda: MD(Mul(A, B)),
    'ma3': lambda: MD(Mul(2, A)),
    'ma4': lambda: MD(Add(A, A)),
    'ma5': lambda: MD(Mul(A, B, A)),
    'ma6': lambda: MD(Add(A, B*A)),
    'ma7': lambda: MD(Mul(x, A)),
    'ma8': lambda: MD(Add(A, ZeroMatrix(2,2))),
    'ma9': lambda: MD(Add(A.T, B)),
    'ma10': lambda: MD(Add(A, B, evaluate=False)),
    'ma12': lambda: MD(Mul(A, C)),
    'ma15': lambda: MD(Add(A, -A)),
    'ma16': lambda: MD(Mul(A, Identity(2))),
    'ma17': lambda: MD(Mul(A, B, evaluate=False)),
    'ma18': lambda: MD(Add(Mul(2, A), B)),
    'ma19': lambda: MD(Mul(A, Inverse(A))),
    'ma20': lambda: MD(Add(A, Matrix([[1,2],[3,4]]))),
    'ma21': lambda: MD(Mul(A, Matrix([[1,2],[3,4]]))),
    'pe0': lambda: str(s2*X),
    'pe1s': lambda: str(p),
    'pe1r': lambda: repr(p),
    'pe2s': lambda: str(q),
    'pe2r': lambda: repr(q),
    'pe3s': lambda: str(R(K.from_sympy(sqrt(2)/2))),
    'pe4s': lambda: str(X - K.from_sympy(1+sqrt(2))),
    'pe5': lambda: str(K.from_sympy(1+sqrt(2))),
    'pe6': lambda: repr(K.from_sympy(1+sqrt(2))),
    'pe7': lambda: str(R.from_expr(x**2*sqrt(2)+3) if False else R(3)*X*s2),
    'pe8': lambda: str((X + s2)**3),
    'pe9': lambda: repr(R),
}

EXPECTED = {'dv1': '(x - y, y**2 + y)',
 'dv2': '(x, 1 - x)',
 'fa1': 'x*y*(x + y + 1)',
 'fl1': '(1, [(x, 1), (y, 1), (x + y + 1, 1)])',
 'fl2': '(1, [(x, 1), (y, 1), (x + y + 1, 1)])',
 'fl3': '(1, [(x*y - 1, 1), (x*y + 1, 1)])',
 'gcd1': 'x - y',
 'la1': "('int64', [2, 3, 4])",
 'la10': "('float', 5.764607523034235e+17)",
 'la2': "('int', 1152921504606846977)",
 'la3': "('int64', [2, 8, 18])",
 'la4': "('float64', [0.5, 1.0, 1.5])",
 'la5': "('int', 1329227995784915876362571574100885504)",
 'la6': "('int', 11)",
 'la7': "('int64', [11])",
 'la8': "('int', 1152921504606846977)",
 'la9': "('float', 0.6666666666666666)",
 'ma1': "('MatAdd', True, False, 'A + B', (2, 2))",
 'ma10': "('MatAdd', True, False, 'A + B', (2, 2))",
 'ma12': "('MatMul', False, True, 'A*C', (2, 3))",
 'ma15': "('ZeroMatrix', False, False, '0', (2, 2))",
 'ma16': "('MatrixSymbol', False, False, 'A', (2, 2))",
 'ma17': "('MatMul', False, True, 'A*B', (2, 2))",
 'ma18': "('MatAdd', True, False, '2*A + B', (2, 2))",
 'ma19': "('Identity', False, False, 'I', (2, 2))",
 'ma2': "('MatMul', False, True, 'A*B', (2, 2))",
 'ma20': "('MatAdd', True, False, 'Matrix([\\n[1, 2],\\n[3, 4]]) + A', (2, 2))",
 'ma21': "('MatMul', False, True, 'A*Matrix([\\n[1, 2],\\n[3, 4]])', (2, 2))",
 'ma3': "('MatMul', False, True, '2*A', (2, 2))",
 'ma4': "('MatMul', False, True, '2*A', (2, 2))",
 'ma5': "('MatMul', False, True, 'A*B*A', (2, 2))",
 'ma6': "('MatAdd', True, False, 'A + B*A', (2, 2))",
 'ma7': "('MatMul', False, True, 'x*A', (2, 2))",
 'ma8': "('MatrixSymbol', False, False, 'A', (2, 2))",
 'ma9': "('MatAdd', True, False, 'B + A.T', (2, 2))",
 'pd1': '(x - y, y**2 + y)',
 'pe0': 'ANP([MPQ(1,1), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)*X',
 'pe1r': 'X**2 + ANP([MPQ(1,1), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)*X + '
         'ANP([MPQ(1,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)',
 'pe1s': 'X**2 + ANP([MPQ(1,1), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)*X + '
         'ANP([MPQ(1,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)',
 'pe2r': '-ANP([MPQ(1,2), MPQ(0,1), MPQ(-11,2), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-10,1), '
         'MPQ(0,1), MPQ(1,1)], QQ)*Y*Z + ANP([MPQ(2,1)], [MPQ(1,1), MPQ(0,1), MPQ(-10,1), '
         'MPQ(0,1), MPQ(1,1)], QQ)',
 'pe2s': '-ANP([MPQ(1,2), MPQ(0,1), MPQ(-11,2), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-10,1), '
         'MPQ(0,1), MPQ(1,1)], QQ)*Y*Z + ANP([MPQ(2,1)], [MPQ(1,1), MPQ(0,1), MPQ(-10,1), '
         'MPQ(0,1), MPQ(1,1)], QQ)',
 'pe3s': 'ANP([MPQ(1,2), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)',
 'pe4s': 'X - ANP([MPQ(-1,1), MPQ(-1,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)',
 'pe5': 'ANP([MPQ(1,1), MPQ(1,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)',
 'pe6': 'ANP([MPQ(1,1), MPQ(1,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)',
 'pe7': 'ANP([MPQ(3,1), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)*X',
 'pe8': 'X**3 + ANP([MPQ(3,1), MPQ(0,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)*X**2 + '
        'ANP([MPQ(6,1)], [MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)*X + ANP([MPQ(2,1), MPQ(0,1)], '
        '[MPQ(1,1), MPQ(0,1), MPQ(-2,1)], QQ)',
 'pe9': 'Polynomial ring in X over QQ<sqrt(2)> with lex order',
 'ps0_a': ' oo     \n'
          '____    \n'
          '\\   `   \n'
          ' \\    1 \n'
          '  \\   --\n'
          '  /    2\n'
          ' /    n \n'
          '/___,   \n'
          'n = 1   ',
 'ps0_l': '\\sum_{n=1}^{\\infty} \\frac{1}{n^{2}}',
 'ps0_u': ' ∞      \n'
          '____    \n'
          '╲       \n'
          ' ╲    1 \n'
          '  ╲   ──\n'
          '  ╱    2\n'
          ' ╱    n \n'
          '╱       \n'
          '‾‾‾‾    \n'
          'n = 1   ',
 'ps1_a': ' 10     \n___     \n\\  `    \n \\     n\n /    x \n/__,    \nn = 0   ',
 'ps1_l': '\\sum_{n=0}^{10} x^{n}',
 'ps1_u': ' 10     \n ___    \n ╲      \n  ╲    n\n  ╱   x \n ╱      \n ‾‾‾    \nn = 0   ',
 'ps2_a': '  2     3      \n'
          ' __    __      \n'
          ' \\ `   \\ `     \n'
          '  )     )   x*y\n'
          ' /_,   /_,     \n'
          'y = 1 x = 1    ',
 'ps2_l': '\\sum_{\\substack{1 \\leq x \\leq 3\\\\1 \\leq y \\leq 2}} x y',
 'ps2_u': '  2     3      \n'
          ' ___   ___     \n'
          ' ╲     ╲       \n'
          '  ╲     ╲      \n'
          '  ╱     ╱   x⋅y\n'
          ' ╱     ╱       \n'
          ' ‾‾‾   ‾‾‾     \n'
          'y = 1 x = 1    ',
 'ps3_a': '     /sin(x)\\\n lim |------|\nx->0+\\  x   /',
 'ps3_l': '\\lim_{x \\to 0^+}\\left(\\frac{\\sin{\\left(x \\right)}}{x}\\right)',
 'ps3_u': '     ⎛sin(x)⎞\n lim ⎜──────⎟\nx─→0⁺⎝  x   ⎠',
 'ps4_a': '     1\n lim -\nx->0-x',
 'ps4_l': '\\lim_{x \\to 0^-} \\frac{1}{x}',
 'ps4_u': '     1\n lim ─\nx─→0⁻x',
 'ps5_a': ' lim x\nx->oo ',
 'ps5_l': '\\lim_{x \\to \\infty} x',
 'ps5_u': 'lim x\nx─→∞ ',
 'ps6_a': '  3     n    \n'
          ' __    __    \n'
          ' \\ `   \\ `   \n'
          '  )     )   x\n'
          ' /_,   /_,   \n'
          'n = 1 x = 1  ',
 'ps6_l': '\\sum_{\\substack{1 \\leq x \\leq n\\\\1 \\leq n \\leq 3}} x',
 'ps6_u': '  3     n    \n'
          ' ___   ___   \n'
          ' ╲     ╲     \n'
          '  ╲     ╲    \n'
          '  ╱     ╱   x\n'
          ' ╱     ╱     \n'
          ' ‾‾‾   ‾‾‾   \n'
          'n = 1 x = 1  ',
 'ps7_a': '    1\nlim -\nx->0x',
 'ps7_l': '\\lim_{x \\to 0} \\frac{1}{x}',
 'ps7_u': '    1\nlim ─\nx─→0x',
 'ps8_a': '  n        \n __        \n \\ `       \n  )   x + 1\n /_,       \nx = 1      ',
 'ps8_l': '\\sum_{x=1}^{n} x + 1',
 'ps8_u': '  n        \n'
          ' ___       \n'
          ' ╲         \n'
          '  ╲        \n'
          '  ╱   x + 1\n'
          ' ╱         \n'
          ' ‾‾‾       \n'
          'x = 1      ',
 'ps9_a': '  b     \n'
          '____    \n'
          '\\   `   \n'
          ' \\     2\n'
          '  \\   n \n'
          '  /   --\n'
          ' /    2 \n'
          '/___,   \n'
          'n = a   ',
 'ps9_l': '\\sum_{n=a}^{b} \\frac{n^{2}}{2}',
 'ps9_u': ' b      \n'
          '____    \n'
          '╲       \n'
          ' ╲     2\n'
          '  ╲   n \n'
          '  ╱   ──\n'
          ' ╱    2 \n'
          '╱       \n'
          '‾‾‾‾    \n'
          'n = a   ',
 'quo1': 'x - y',
 'rem1': 'y**2 + y',
 'ri1': 'x*y > 1',
 'ri2': 'x*y < 1',
 'ri3': 'x*y >= 1',
 'ri4': 'x > 1 - y',
 'ri5': 'x*y**2 > 1',
 'ri6': '(-oo < x) & (x < 3) & (x*y > 1)',
 'ri7': 'x*y > 1',
 'ri9': 'x > 1/p',
 'sq1': '(1, [(x - 1, 1), (x + y, 2)])',
 'sq2': '(1, [(x*y + 1, 2)])',
 'sq3': '(x - 1)*(x + y)**2',
 'ter1': 'x*y*(x + 1)',
 'vi1': '(Integral(x, x))*N.i + (Integral(x**2, x))*N.j + (Integral(sin(x), x))*N.k',
 'vi2': '(Integral(x, (x, 0, 1)))*N.i + (Integral(x**2, (x, 0, 1)))*N.j + (Integral(sin(x), (x, 0, '
        '1)))*N.k',
 'vi3': 'N.x**2/2*N.i',
 'vi4': 'x**2/2*N.i + x**3/3*N.j + (-cos(x))*N.k',
 'vi5': '3*x*N.j',
 'vi6': '(Integral(x, x))*N.i + (Integral(x**2, x))*N.j + (Integral(sin(x), x))*N.k'}


def _value(func):
    try:
        r = func()
    except Exception as e:
        return "ERR %s: %s" % (type(e).__name__, str(e)[:120])
    return r if isinstance(r, str) else repr(r)


@pytest.mark.parametrize('name', sorted(CASES))
def test_probe(name):
    assert _value(CASES[name]) == EXPECTED[name]


def test_multisymbol_inequality_not_solvable():
    with pytest.raises(NotImplementedError):
        reduce_inequalities(x**2*y > 1, [x])


def test_reduce_inequalities_two_symbols_of_interest():
    with pytest.raises(NotImplementedError):
        reduce_inequalities(x*y > 1, [x, y])


def test_matrix_shape_errors():
    from symox import ShapeError
    with pytest.raises(ShapeError):
        Add(A, C)
    with pytest.raises(ShapeError):
        Mul(C, A)
    with pytest.raises(TypeError):
        A + 1


# --- the wrt/sort option of the polys functions ---------------------------
_f = x**2*y + x*y**2 + x*y
_g = x*y + y**2
_h = (x + y)**2*(x - 1)
UNARY = ['factor_list','sqf_list','factor','sqf','sqf_part','terms_gcd','degree','LC','LM','LT','monic','content','primitive','degree_list','total_degree','sqf_norm','gff_list','decompose','sturm','discriminant','trunc','real_roots','nroots','ground_roots','intervals','count_roots','apart','cancel','together','horner','poly','Poly']
BINARY = ['div','rem','quo','pdiv','prem','pquo','exquo','pexquo','gcd','lcm','cofactors','gcdex','half_gcdex','invert','resultant','subresultants','compose','reduced','groebner']


def _wrt_value(fn, w):
    F = globals()[fn]
    try:
        if fn in UNARY:
            arg = _h if fn in ('sqf_list', 'sqf', 'sqf_part', 'gff_list') else _f
            if fn in ('sturm', 'real_roots', 'nroots', 'ground_roots', 'intervals', 'count_roots',
                      'discriminant', 'decompose', 'trunc', 'apart', 'sqf_norm'):
                arg = x**2*y + 2*x + y if fn != 'trunc' else 7*x**2*y + 3*x
            if fn == 'trunc':
                r = F(arg, 3, wrt=w)
            elif fn == 'apart':
                r = F(1/((x + 1)*(x + y)), wrt=w)
            else:
                r = F(arg, wrt=w)
        elif fn == 'reduced':
            r = F(_f, [_g], wrt=w)
        elif fn == 'groebner':
            r = F([_f, _g], wrt=w)
        elif fn in ('gcdex', 'half_gcdex', 'invert'):
            r = F(x**2 + y, x + y, wrt=w)
        elif fn in ('exquo', 'pexquo'):
            r = F(_f, x*y, wrt=w)
        else:
            r = F(_f, _g, wrt=w)
        return str(r)
    except Exception as e:
        return "ERR %s" % type(e).__name__


WRT_EXPECTED = {'LC|x': '1',
 'LC|y': '1',
 'LC|y,x': '1',
 'LM|x': 'x**2*y',
 'LM|y': 'x*y**2',
 'LM|y,x': 'x*y**2',
 'LT|x': 'x**2*y',
 'LT|y': 'x*y**2',
 'LT|y,x': 'x*y**2',
 'Poly|x': "Poly(x**2*y + x*y**2 + x*y, x, y, domain='ZZ')",
 'Poly|y': "Poly(y**2*x + y*x**2 + y*x, y, x, domain='ZZ')",
 'Poly|y,x': "Poly(y**2*x + y*x**2 + y*x, y, x, domain='ZZ')",
 'apart|x': 'ERR NotImplementedError',
 'apart|y': 'ERR NotImplementedError',
 'apart|y,x': 'ERR NotImplementedError',
 'cancel|x': 'x**2*y + x*y**2 + x*y',
 'cancel|y': 'x**2*y + x*y**2 + x*y',
 'cancel|y,x': 'x**2*y + x*y**2 + x*y',
 'cofactors|x': '(y, x**2 + x*y + x, x + y)',
 'cofactors|y': '(y, x**2 + x*y + x, x + y)',
 'compose|x': 'x**2*y**3 + 2*x*y**4 + x*y**3 + x*y**2 + y**5 + y**4 + y**3',
 'compose|y': 'x**3*y**2 + x**3*y + 2*x**2*y**3 + x**2*y**2 + x**2*y + x*y**4 + x*y**2',
 'content|x': '1',
 'content|y': '1',
 'content|y,x': '1',
 'count_roots|x': 'ERR TypeError',
 'count_roots|y': 'ERR TypeError',
 'count_roots|y,x': 'ERR TypeError',
 'decompose|x': 'ERR ValueError',
 'decompose|y': 'ERR ValueError',
 'decompose|y,x': 'ERR ValueError',
 'degree_list|x': '(2, 2)',
 'degree_list|y': '(2, 2)',
 'degree_list|y,x': '(2, 2)',
 'degree|x': 'ERR TypeError',
 'degree|y': 'ERR TypeError',
 'degree|y,x': 'ERR TypeError',
 'discriminant|x': '4 - 4*y**2',
 'discriminant|y': '1',
 'discriminant|y,x': '1',
 'div|x': '(x + 1, -y**2)',
 'div|y': '(x, x*y)',
 'exquo|x': 'x + y + 1',
 'exquo|y': 'x + y + 1',
 'factor_list|x': '(1, [(x, 1), (y, 1), (x + y + 1, 1)])',
 'factor_list|y': '(1, [(x, 1), (y, 1), (x + y + 1, 1)])',
 'factor_list|y,x': '(1, [(x, 1), (y, 1), (x + y + 1, 1)])',
 'factor|x': 'x*y*(x + y + 1)',
 'factor|y': 'x*y*(x + y + 1)',
 'factor|y,x': 'x*y*(x + y + 1)',
 'gcdex|x': 'ERR ValueError',
 'gcdex|y': 'ERR ValueError',
 'gcd|x': 'y',
 'gcd|y': 'y',
 'gff_list|x': 'ERR ValueError',
 'gff_list|y': 'ERR ValueError',
 'gff_list|y,x': 'ERR ValueError',
 'groebner|x': "GroebnerBasis([x*y, y**2], x, y, domain='ZZ', order='lex')",
 'groebner|y': "GroebnerBasis([y**2, x*y], y, x, domain='ZZ', order='lex')",
 'ground_roots|x': 'ERR MultivariatePolynomialError',
 'ground_roots|y': 'ERR MultivariatePolynomialError',
 'ground_roots|y,x': 'ERR MultivariatePolynomialError',
 'half_gcdex|x': 'ERR ValueError',
 'half_gcdex|y': 'ERR ValueError',
 'horner|x': 'x*(x*y + y*(y + 1))',
 'horner|y': 'y*(x*y + x*(x + 1))',
 'horner|y,x': 'y*(x*y + x*(x + 1))',
 'intervals|x': 'ERR TypeError',
 'intervals|y': 'ERR TypeError',
 'intervals|y,x': 'ERR TypeError',
 'invert|x': 'ERR ValueError',
 'invert|y': 'ERR ValueError',
 'lcm|x': 'x**3*y + 2*x**2*y**2 + x**2*y + x*y**3 + x*y**2',
 'lcm|y': 'x**3*y + 2*x**2*y**2 + x**2*y + x*y**3 + x*y**2',
 'monic|x': 'x**2*y + x*y**2 + x*y',
 'monic|y': 'x**2*y + x*y**2 + x*y',
 'monic|y,x': 'x**2*y + x*y**2 + x*y',
 'nroots|x': 'ERR TypeError',
 'nroots|y': 'ERR TypeError',
 'nroots|y,x': 'ERR TypeError',
 'pdiv|x': '(x*y**2 + y**2, -y**4)',
 'pdiv|y': '(x, x*y)',
 'pexquo|x': 'x*y**2 + y**3 + y**2',
 'pexquo|y': 'x**3 + x**2*y + x**2',
 'poly|x': "Poly(x**2*y + x*y**2 + x*y, x, y, domain='ZZ')",
 'poly|y': "Poly(y**2*x + y*x**2 + y*x, y, x, domain='ZZ')",
 'poly|y,x': "Poly(y**2*x + y*x**2 + y*x, y, x, domain='ZZ')",
 'pquo|x': 'x*y**2 + y**2',
 'pquo|y': 'x',
 'prem|x': '-y**4',
 'prem|y': 'x*y',
 'primitive|x': '(1, x**2*y + x*y**2 + x*y)',
 'primitive|y': '(1, x**2*y + x*y**2 + x*y)',
 'primitive|y,x': '(1, x**2*y + x*y**2 + x*y)',
 'quo|x': 'x + 1',
 'quo|y': 'x',
 'real_roots|x': 'ERR TypeError',
 'real_roots|y': 'ERR TypeError',
 'real_roots|y,x': 'ERR TypeError',
 'reduced|x': '([x + 1], -y**2)',
 'reduced|y': '([x], x*y)',
 'rem|x': '-y**2',
 'rem|y': 'x*y',
 'resultant|x': '-y**4',
 'resultant|y': '0',
 'sqf_list|x': '(1, [(x - 1, 1), (x + y, 2)])',
 'sqf_list|y': '(1, [(x - 1, 1), (x + y, 2)])',
 'sqf_list|y,x': '(1, [(x - 1, 1), (x + y, 2)])',
 'sqf_norm|x': 'ERR DomainError',
 'sqf_norm|y': 'ERR DomainError',
 'sqf_norm|y,x': 'ERR DomainError',
 'sqf_part|x': 'x**2 + x*y - x - y',
 'sqf_part|y': 'x**2 + x*y - x - y',
 'sqf_part|y,x': 'x**2 + x*y - x - y',
 'sqf|x': '(x - 1)*(x + y)**2',
 'sqf|y': '(x - 1)*(x + y)**2',
 'sqf|y,x': '(x - 1)*(x + y)**2',
 'sturm|x': 'ERR ValueError',
 'sturm|y': 'ERR ValueError',
 'sturm|y,x': 'ERR ValueError',
 'subresultants|x': '[x**2*y + x*y**2 + x*y, x*y + y**2, -y**4]',
 'subresultants|y': '[x**2*y + x*y**2 + x*y, x*y + y**2, -x*y]',
 'terms_gcd|x': 'x*y*(x + y + 1)',
 'terms_gcd|y': 'x*y*(x + y + 1)',
 'terms_gcd|y,x': 'x*y*(x + y + 1)',
 'together|x': 'ERR TypeError',
 'together|y': 'ERR TypeError',
 'together|y,x': 'ERR TypeError',
 'total_degree|x': 'ERR TypeError',
 'total_degree|y': 'ERR TypeError',
 'total_degree|y,x': 'ERR TypeError',
 'trunc|x': 'x**2*y',
 'trunc|y': 'x**2*y',
 'trunc|y,x': 'x**2*y'}


@pytest.mark.parametrize('key', sorted(WRT_EXPECTED))
def test_wrt_option(key):
    fn, w = key.split('|')
    assert _wrt_value(fn, w) == WRT_EXPECTED[key]


# --- lambdify keeps integer dtypes and exact Python ints --------------------
LAMBDIFY_CASES = [
    (x, x + 1, (np.int64(3),), 'int64 int64 4'),
    (x, x**3 - 2*x, (np.array([[1, 2], [3, 4]]),), 'ndarray int64 [[-1, 4], [21, 56]]'),
    (x, x + 1, (np.array([1, 2], dtype=np.uint8),), 'ndarray uint8 [2, 3]'),
    ((x, y), x*y, (np.array([1, 2]), 2.5), 'ndarray float64 [2.5, 5.0]'),
    ((x, y), x*y + y**2, (np.array([1, 2]), np.array([3, 4])), 'ndarray int64 [12, 24]'),
    (x, 1/x, (np.array([1, 2]),), 'ndarray float64 [1.0, 0.5]'),
    (x, sin(x), (np.array([1, 2]),), 'ndarray float64 [0.8414709848078965, 0.9092974268256817]'),
    (x, x**2, (np.array([True, False]),), 'ndarray int8 [1, 0]'),
    ((x, y), x - y, (np.int32(2), 3), 'int32 int32 -1'),
    (x, Integer(5), (np.array([1, 2]),), 'int None 5'),
    (x, 2**70*x, (np.array([1]),), 'ERR OverflowError'),
    (x, x + 1, (np.array([1, 2], dtype=np.int32),), 'ndarray int32 [2, 3]'),
    (x, sqrt(2)*x, (np.array([1, 2]),), 'ndarray float64 [1.4142135623730951, 2.8284271247461903]'),
]


@pytest.mark.parametrize('args, expr, values, expected', LAMBDIFY_CASES)
def test_lambdify_numpy_dtypes(args, expr, values, expected):
    try:
        r = lambdify(args, expr)(*values)
        got = '%s %s %s' % (type(r).__name__, getattr(r, 'dtype', None),
                             r.tolist() if hasattr(r, 'tolist') else r)
    except Exception as e:
        got = 'ERR %s' % type(e).__name__
    assert got == expected


# --- pretty and latex of Sum and Limit -------------------------------------
n, a, b = symbols('n a b')
_F = Function('f')
PRINT_EXPRS = [
    Sum(x + 1, (x, 1, n)), 2*Sum(1/x, (x, 1, n)), Limit(x + 1, x, 1), Limit(sin(x)/x, x, 0)**2,
    None, Limit(exp(x), x, -oo), Sum(_F(n), (n, 1, n**2)), Sum(Rational(1, 2)**n, (n, 0, oo)),
    Limit(x*y, x, 0), Limit(-x, x, 0), Sum(x, (x, 1, n))*y, Limit(_F(x), x, a, '-'),
    Sum(x**y/(x + 1), (x, a, b)), Sum(sin(x), (x, 0, n), (n, 0, 3)),
]
PRINT_EXPECTED = {'0a': '  n          \n __          \n \\ `         \n  )   (x + 1)\n /_,         \nx = 1        ',
 '0l': '\\sum_{x=1}^{n} \\left(x + 1\\right)',
 '0u': '  n          \n'
       ' ___         \n'
       ' ╲           \n'
       '  ╲          \n'
       '  ╱   (x + 1)\n'
       ' ╱           \n'
       ' ‾‾‾         \n'
       'x = 1        ',
 '10a': '    n    \n   __    \n   \\ `   \ny*  )   x\n   /_,   \n  x = 1  ',
 '10l': 'y \\sum_{x=1}^{n} x',
 '10u': '    n    \n   ___   \n   ╲     \n    ╲    \ny⋅  ╱   x\n   ╱     \n   ‾‾‾   \n  x = 1  ',
 '11a': ' lim f(x)\nx->a-    ',
 '11l': '\\lim_{x \\to a^-} f{\\left(x \\right)}',
 '11u': ' lim f(x)\nx─→a⁻    ',
 '12a': '  b        \n'
        '____       \n'
        '\\   `      \n'
        ' \\      y  \n'
        '  \\    x   \n'
        '  /   -----\n'
        ' /    x + 1\n'
        '/___,      \n'
        'x = a      ',
 '12l': '\\sum_{x=a}^{b} \\frac{x^{y}}{x + 1}',
 '12u': ' b         \n'
        '____       \n'
        '╲          \n'
        ' ╲      y  \n'
        '  ╲    x   \n'
        '  ╱   ─────\n'
        ' ╱    x + 1\n'
        '╱          \n'
        '‾‾‾‾       \n'
        'x = a      ',
 '13a': '  3     n         \n'
        ' __    __         \n'
        ' \\ `   \\ `        \n'
        '  )     )   sin(x)\n'
        ' /_,   /_,        \n'
        'n = 0 x = 0       ',
 '13l': '\\sum_{\\substack{0 \\leq x \\leq n\\\\0 \\leq n \\leq 3}} \\sin{\\left(x \\right)}',
 '13u': '  3     n         \n'
        ' ___   ___        \n'
        ' ╲     ╲          \n'
        '  ╲     ╲         \n'
        '  ╱     ╱   sin(x)\n'
        ' ╱     ╱          \n'
        ' ‾‾‾   ‾‾‾        \n'
        'n = 0 x = 0       ',
 '1a': '   n     \n  ___    \n  \\  `   \n   \\    1\n2*  )   -\n   /    x\n  /__,   \n  x = 1  ',
 '1l': '2 \\sum_{x=1}^{n} \\frac{1}{x}',
 '1u': '   n     \n'
       '  ____   \n'
       '  ╲      \n'
       '   ╲     \n'
       '    ╲   1\n'
       '2⋅  ╱   ─\n'
       '   ╱    x\n'
       '  ╱      \n'
       '  ‾‾‾‾   \n'
       '  x = 1  ',
 '2a': ' lim (x + 1)\nx->1+       ',
 '2l': '\\lim_{x \\to 1^+}\\left(x + 1\\right)',
 '2u': ' lim (x + 1)\nx─→1⁺       ',
 '3a': '               2\n/     /sin(x)\\\\ \n| lim |------|| \n\\x->0+\\  x   // ',
 '3l': '\\left(\\lim_{x \\to 0^+}\\left(\\frac{\\sin{\\left(x \\right)}}{x}\\right)\\right)^{2}',
 '3u': '               2\n⎛     ⎛sin(x)⎞⎞ \n⎜ lim ⎜──────⎟⎟ \n⎝x─→0⁺⎝  x   ⎠⎠ ',
 '5a': '       x\n lim  e \nx->-oo  ',
 '5l': '\\lim_{x \\to -\\infty} e^{x}',
 '5u': '      x\n lim ℯ \nx─→-∞  ',
 '6a': '  2       \n n        \n __       \n \\ `      \n  )   f(n)\n /_,      \nn = 1     ',
 '6l': '\\sum_{n=1}^{n^{2}} f{\\left(n \\right)}',
 '6u': '  2       \n'
       ' n        \n'
       ' ___      \n'
       ' ╲        \n'
       '  ╲       \n'
       '  ╱   f(n)\n'
       ' ╱        \n'
       ' ‾‾‾      \n'
       'n = 1     ',
 '7a': ' oo      \n___      \n\\  `     \n \\     -n\n /    2  \n/__,     \nn = 0    ',
 '7l': '\\sum_{n=0}^{\\infty} \\left(\\frac{1}{2}\\right)^{n}',
 '7u': '  ∞      \n ___     \n ╲       \n  ╲    -n\n  ╱   2  \n ╱       \n ‾‾‾     \nn = 0    ',
 '8a': ' lim (x*y)\nx->0+     ',
 '8l': '\\lim_{x \\to 0^+}\\left(x y\\right)',
 '8u': ' lim (x⋅y)\nx─→0⁺     ',
 '9a': ' lim (-x)\nx->0+    ',
 '9l': '\\lim_{x \\to 0^+}\\left(- x\\right)',
 '9u': ' lim (-x)\nx─→0⁺    '}


@pytest.mark.parametrize('key', sorted(PRINT_EXPECTED))
def test_sum_limit_printing(key):
    e = PRINT_EXPRS[int(key[:-1])]
    mode = key[-1]
    got = latex(e) if mode == 'l' else pretty(e, use_unicode=(mode == 'u'))
    assert got == PRINT_EXPECTED[key]


PS_EXPRS = [Sum(1/n**2, (n, 1, oo)), Sum(x**n, (n, 0, 10)), Sum(x*y, (x, 1, 3), (y, 1, 2)),
            Limit(sin(x)/x, x, 0), Limit(1/x, x, 0, '-'), Limit(x, x, oo), Sum(Sum(x, (x, 1, n)), (n, 1, 3)),
            Limit(1/x, x, 0, '+-'), Sum(x, (x, 1, n)) + 1, Sum(n**2/2, (n, a, b))]


@pytest.mark.parametrize('key', sorted(k for k in EXPECTED if k.startswith('ps')))
def test_sum_limit_printing_basic(key):
    index, mode = key[2:].split('_')
    e = PS_EXPRS[int(index)]
    got = latex(e) if mode == 'l' else pretty(e, use_unicode=(mode == 'u'))
    assert got == EXPECTED[key]
