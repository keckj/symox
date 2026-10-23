"""Regression tests for the simp-rest wave: cse (opt_cse + tree_cse port),
Function.inverse, inversecombine, trigsimp(inverse=True) for asec/acsc,
_TR11/_TR56, Osborne's rule with sech/csch, TR4 on held arguments, the
_eval_rewrite_as_exp hook and kroneckersimp.  Expected values are the str()
outputs of SymPy 1.14."""
import pytest

from symox import *

# ``symox.test`` (the SymPy-style test runner) must not be collected here.
del test
from symox.simplify.cse_main import cse, tree_cse, opt_cse
from symox.simplify.cse_opts import sub_pre, sub_post
from symox.simplify._cse_diff import _forward_jacobian
from symox.simplify.fu import _TR11, _TR56, _osborne, _osbornei, TR4
from symox.simplify.simplify import inversecombine, kroneckersimp

x,y,z,w=symbols('x y z w'); d=Dummy('d')
A,B=symbols('A B',commutative=False)
i,j=symbols('i j'); K=KroneckerDelta
g=Function('g')
class f(Function):
    def inverse(self, argindex=1):
        return g
X=MatrixSymbol('X',3,3); Y=MatrixSymbol('Y',3,3)
with evaluate(False):
    held=[cos(v) for v in (0, pi/6, pi/4, pi/3, pi/2)]
cases={
 'cse_basic': lambda: cse([(x+y)**2+sqrt(x+y), (x+y)*z]),
 'cse_4020': lambda: cse([x**5+x**3, x**5+x**4]),
 'cse_4498': lambda: cse(w/(x - y) + z/(y - x), optimizations='basic'),
 'cse_nc': lambda: cse(A*B + B*A*B),
 'cse_piecewise': lambda: cse(Piecewise((x+y, x < 1), (x*y, True)) + (x+y)**2),
 'cse_8891': lambda: [cse(m) for m in (Matrix([[x+y, (x+y)**2], [x, 1]]), ImmutableMatrix([[x+y, (x+y)**2], [x, 1]]))],
 'cse_tuple': lambda: cse((x+y, (x+y)*z, Tuple(x+y, z))),
 'cse_ignore': lambda: cse([(x+y)*exp(y), (x+y)*exp(2*y)], ignore=(y,)),
 'cse_matsym': lambda: cse([(X*Y)[0,0]*(X*Y)[0,0] + (X*Y)[0,0], (X*Y)[0,0]]),
 'cse_matexpr': lambda: cse(X*Y + (X*Y).T),
 'cse_postprocess': lambda: cse([x*y+z, x*y], postprocess=lambda r, e: (r, [i+1 for i in e])),
 'cse_list_false': lambda: cse(x*y+(x*y)**2, list=False),
 'cse_none': lambda: cse([x*y+(x*y)**2, sin(x*y)], order='none'),
 'sub_pre': lambda: sub_pre(-x - y + z),
 'sub_post': lambda: sub_post(sub_pre(-x - y + z)),
 'fjac': lambda: _forward_jacobian(Matrix([x*y+sin(x), x**2*z]), Matrix([x, y, z])),
 'inverse': lambda: [asin.inverse(None), exp(x).inverse(), log(x).inverse(), erf(x).inverse(), tanh(x).inverse(), hasattr(sin(x), 'inverse')],
 'inversecombine': lambda: [inversecombine(f(g(x))), inversecombine(asin(sin(x))), inversecombine(2*log(exp(3*x))), inversecombine(exp(log(x)))],
 'simplify_inverse': lambda: [simplify(f(g(x))), simplify(f(g(x)), inverse=True), simplify(f(g(sin(x)**2 + cos(x)**2)), inverse=True), simplify(f(g(x, y)), inverse=True), simplify(2*asin(sin(3*x)), inverse=True)],
 'trigsimp_inv': lambda: [trigsimp(asec(sec(x)), inverse=True), trigsimp(acsc(csc(x)), inverse=True), trigsimp(asin(sin(x)) + acsc(csc(y)), inverse=True)],
 'TR11': lambda: [_TR11(sin(x/3)*sin(2*x)*sin(x/4)/(cos(x/6)*cos(x/8))), _TR11(sin(x/3)/(cos(x/6)))],
 'TR56': lambda: [_TR56(sin(x)**3, sin, cos, lambda x: 1 - x, 4, False), _TR56(sin(x)**10, sin, cos, lambda x: 1 - x, 10, True), _TR56(sin(x)**6, sin, cos, lambda x: 1 - x, 6, False), _TR56(sin(x)**8, sin, cos, lambda x: 1 - x, 10, True)],
 'osborne': lambda: [_osborne(e, d) for e in (sinh(x), cosh(x), tanh(x), coth(x), sech(x), csch(x))],
 'osbornei': lambda: [_osbornei(e, d) for e in (sin(d*x), cos(d*x), tan(d*x), cot(d*x), sec(d*x), csc(d*x))],
 'hyper_trig': lambda: [trigsimp(e) for e in (tanh(x)**2 + sech(x)**2, csch(x)**2 - coth(x)**2, cosh(x)*sech(x), sinh(2*x)*csch(x), sinh(x)**2*sech(x)**2 + sech(x)**2)],
 'TR4': lambda: [TR4(e) for e in held] + [TR4(sin(x) + held[1])],
 'rewrite_exp': lambda: [cos(x)._eval_rewrite_as_exp(x), tan(x)._eval_rewrite_as_exp(x), sinh(x)._eval_rewrite_as_exp(x)],
 'simplify_complex': lambda: simplify(cos(x)._eval_rewrite_as_exp(x)*tan(x)._eval_rewrite_as_exp(x)),
 'kronecker': lambda: [simplify(K(0, j).rewrite(Piecewise) * K(1, j)), simplify(K(1, i) + Piecewise((1, Eq(j, 2)), (0, True))), simplify(K(0, j) * K(1, j)), kroneckersimp(K(0, j).rewrite(Piecewise) * K(1, j)), kroneckersimp(K(1, i) + Piecewise((1, Eq(j, 2)), (0, True))), kroneckersimp(1 + K(0, j)*K(1, j))],
}

EXPECTED = {
    'cse_basic': '([(x0, x + y)], [sqrt(x0) + x0**2, x0*z])',
    'cse_4020': '([(x0, x**5)], [x**3 + x0, x**4 + x0])',
    'cse_4498': '([], [(w - z)/(x - y)])',
    'cse_nc': '([], [A*B + B*A*B])',
    'cse_piecewise': '([(x0, x + y)], [x0**2 + Piecewise((x0, x < 1), (x*y, True))])',
    'cse_8891': '[([(x0, x + y)], [Matrix([\n[x0, x0**2],\n[ x,     1]])]), ([(x0, x + y)], [Matrix([\n[x0, x0**2],\n[ x,     1]])])]',
    'cse_tuple': '([(x0, x + y)], [x0, x0*z, (x0, z)])',
    'cse_ignore': '([], [(x + y)*exp(y), (x + y)*exp(2*y)])',
    'cse_matsym': '([(x0, X[0, 0]*Y[0, 0] + X[0, 1]*Y[1, 0] + X[0, 2]*Y[2, 0])], [x0**2 + x0, x0])',
    'cse_matexpr': '([], [X*Y + Y.T*X.T])',
    'cse_postprocess': '([(x0, x*y)], [x0 + z + 1, x0 + 1])',
    'cse_list_false': '([], x**2*y**2 + x*y)',
    'cse_none': '([(x0, x*y)], [x**2*y**2 + x0, sin(x0)])',
    'sub_pre': '-(x + y - z)',
    'sub_post': '-(x + y - z)',
    'fjac': 'Matrix([[y + cos(x), x, 0], [2*x*z, 0, x**2]])',
    'inverse': '[sin, log, exp, erfinv, atanh, False]',
    'inversecombine': '[x, x, 6*x, x]',
    'simplify_inverse': '[f(g(x)), x, 1, f(g(x, y)), 6*x]',
    'trigsimp_inv': '[x, x, x + y]',
    'TR11': '[4*sin(x/8)*sin(x/6)*sin(2*x), 2*sin(x/6)]',
    'TR56': '[(1 - cos(x)**2)*sin(x), sin(x)**10, (1 - cos(x)**2)**3, (1 - cos(x)**2)**4]',
    'osborne': '[I*sin(_d*x), cos(_d*x), I*tan(_d*x), -I*cot(_d*x), sec(_d*x), -I*csc(_d*x)]',
    'osbornei': '[-I*sinh(x), cosh(x), -I*tanh(x), I*coth(x), sech(x), I*csch(x)]',
    'hyper_trig': '[1, -1, 1, 2*cosh(x), 1]',
    'TR4': '[1, sqrt(3)/2, sqrt(2)/2, 1/2, 0, sin(x) + sqrt(3)/2]',
    'rewrite_exp': '[exp(I*x)/2 + exp(-I*x)/2, I*(-exp(I*x) + exp(-I*x))/(exp(I*x) + exp(-I*x)), exp(x)/2 - exp(-x)/2]',
    'simplify_complex': 'sin(x)',
    'kronecker': '[0, KroneckerDelta(1, i) + KroneckerDelta(2, j), 0, 0, KroneckerDelta(1, i) + KroneckerDelta(2, j), 1]',
}


@pytest.mark.parametrize("name", sorted(EXPECTED))
def test_matches_sympy(name):
    assert str(cases[name]()) == EXPECTED[name]
