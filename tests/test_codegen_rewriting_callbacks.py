"""Callable queries pass the whole matched expression to their replacement."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp
from symox.codegen.rewriting import ReplaceOptim,exp2_opt,log1p_opt
from sympy.codegen.rewriting import ReplaceOptim as Reference


@pytest.mark.parametrize('source',['x**2','sin(x)**3','x**2+y**3','sin(x**2)','(x+y)**3'])
def test_callable_replacement_receives_one_expression(source):
    value=sx.sympify(source);expected=sp.sympify(source);seen=[]
    def replace(expression):
        seen.append(expression)
        return expression.base+expression.exp
    actual=ReplaceOptim(lambda node:node.is_Pow,replace)(expr=value)
    reference=Reference(lambda node:node.is_Pow,lambda node:node.base+node.exp)(expr=expected)
    assert actual.to_sympy()==reference
    assert seen and all(node.is_Pow for node in seen)


def test_function_class_query_still_passes_function_arguments():
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=ReplaceOptim(sx.sin,sx.cos)(sx.sin(x)+sx.sin(x*x))
    expected=Reference(sp.sin,sp.cos)(sp.sin(X)+sp.sin(X*X))
    assert actual.to_sympy()==expected


def test_atom_queries_are_not_skipped():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    for expression,reference in [(x,X),(x+1,X+1)]:
        actual=ReplaceOptim(lambda node:node==x,lambda node:y)(expression)
        expected=Reference(lambda node:node==X,lambda node:Y)(reference)
        assert actual.to_sympy()==expected


def test_builtin_rewrites_keep_their_correct_arity():
    x=sx.Symbol('x')
    assert str(exp2_opt(2**x))=='exp2(x)'
    assert str(log1p_opt(sx.log(1+x)))=='log1p(x)'


def test_rewrite_callbacks_do_not_import_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
from symox.codegen.rewriting import ReplaceOptim,exp2_opt
x=Symbol('x')
assert ReplaceOptim(lambda n:n.is_Pow,lambda n:n.base+n.exp)(x*x)==x+2
assert str(exp2_opt(2**x))=='exp2(x)'
'''
    subprocess.run([sys.executable,'-c',source],check=True)


# Ports of SymPy 1.14's codegen/tests/test_rewriting.py for the optimizations
# producing the numpy_nodes and scipy_nodes functions.
from symox import Pow, Rational, S, Symbol, MatrixSymbol, Q, assuming, cos, exp, log, pi, sin, sinc
from symox.codegen.cfunctions import expm1, log1p, log2
from symox.codegen.matrix_nodes import MatrixSolve
from symox.codegen.numpy_nodes import logaddexp, logaddexp2
from symox.codegen.scipy_nodes import cosm1, powm1
from symox.codegen.rewriting import (
    optimize, cosm1_opt, expm1_opt, powm1_opt, matinv_opt, logaddexp_opt, logaddexp2_opt,
    optims_numpy, optims_scipy, sinc_opts, FuncMinusOneOptim)
from symox.utilities import lambdify

try:
    import numpy
except ImportError:
    numpy = None
try:
    import scipy
except ImportError:
    scipy = None

def test_cosm1_opt():
    x = Symbol('x')

    expr1 = cos(x) - 1
    opt1 = optimize(expr1, [cosm1_opt])
    assert cosm1(x) - opt1 == 0
    assert opt1.rewrite(cos) == expr1

    expr2 = 3*cos(x) - 3
    opt2 = optimize(expr2, [cosm1_opt])
    assert 3*cosm1(x) == opt2
    assert opt2.rewrite(cos) == expr2

    expr3 = 3*cos(x) - 5
    opt3 = optimize(expr3, [cosm1_opt])
    assert 3*cosm1(x) - 2 == opt3
    assert opt3.rewrite(cos) == expr3
    cosm1_opt_non_opportunistic = FuncMinusOneOptim(cos, cosm1, opportunistic=False)
    assert expr3 == optimize(expr3, [cosm1_opt_non_opportunistic])
    assert opt1 == optimize(expr1, [cosm1_opt_non_opportunistic])
    assert opt2 == optimize(expr2, [cosm1_opt_non_opportunistic])

    expr4 = 3*cos(x) + log(x) - 3
    opt4 = optimize(expr4, [cosm1_opt])
    assert 3*cosm1(x) + log(x) == opt4
    assert opt4.rewrite(cos) == expr4

    expr5 = 3*cos(2*x) - 3
    opt5 = optimize(expr5, [cosm1_opt])
    assert 3*cosm1(2*x) == opt5
    assert opt5.rewrite(cos) == expr5

    expr6 = 2 - 2*cos(x)
    opt6 = optimize(expr6, [cosm1_opt])
    assert -2*cosm1(x) == opt6
    assert opt6.rewrite(cos) == expr6


def test_cosm1_two_cos_terms():
    x, y = map(Symbol, 'x y'.split())
    expr1 = cos(x) + cos(y) - 2
    opt1 = optimize(expr1, [cosm1_opt])
    assert opt1 == cosm1(x) + cosm1(y)


def test_expm1_cosm1_mixed():
    x = Symbol('x')
    expr1 = exp(x) + cos(x) - 2
    opt1 = optimize(expr1, [expm1_opt, cosm1_opt])
    assert opt1 == cosm1(x) + expm1(x)


def _check_num_lambdify(expr, opt, val_subs, approx_ref, lambdify_kw=None, poorness=1e10):
    """ poorness=1e10 signifies that `expr` loses precision of at least ten decimal digits. """
    num_ref = expr.subs(val_subs).evalf()
    eps = numpy.finfo(numpy.float64).eps
    assert abs(num_ref - approx_ref) < approx_ref*eps
    f1 = lambdify(list(val_subs.keys()), opt, **(lambdify_kw or {}))
    args_float = tuple(map(float, val_subs.values()))
    num_err1 = abs(f1(*args_float) - approx_ref)
    assert num_err1 < abs(num_ref*eps)
    f2 = lambdify(list(val_subs.keys()), expr, **(lambdify_kw or {}))
    num_err2 = abs(f2(*args_float) - approx_ref)
    assert num_err2 > abs(num_ref*eps*poorness)   # this only ensures that the *test* works as intended


def test_cosm1_apart():
    x = Symbol('x')

    expr1 = 1/cos(x) - 1
    opt1 = optimize(expr1, [cosm1_opt])
    assert opt1 == -cosm1(x)/cos(x)
    if scipy:
        _check_num_lambdify(expr1, opt1, {x: S(10)**-30}, 5e-61, lambdify_kw={"modules": 'scipy'})

    expr2 = 2/cos(x) - 2
    opt2 = optimize(expr2, optims_scipy)
    assert opt2 == -2*cosm1(x)/cos(x)
    if scipy:
        _check_num_lambdify(expr2, opt2, {x: S(10)**-30}, 1e-60, lambdify_kw={"modules": 'scipy'})

    expr3 = pi/cos(3*x) - pi
    opt3 = optimize(expr3, [cosm1_opt])
    assert opt3 == -pi*cosm1(3*x)/cos(3*x)
    if scipy:
        _check_num_lambdify(expr3, opt3, {x: S(10)**-30/3}, float(5e-61*pi), lambdify_kw={"modules": 'scipy'})


def test_powm1():
    args = x, y = map(Symbol, "xy")

    expr1 = x**y - 1
    opt1 = optimize(expr1, [powm1_opt])
    assert opt1 == powm1(x, y)
    for arg in args:
        assert expr1.diff(arg) == opt1.diff(arg)
    if scipy and tuple(map(int, scipy.version.version.split('.')[:3])) >= (1, 10, 0):
        subs1_a = {x: Rational(*(1.0+1e-13).as_integer_ratio()), y: pi}
        ref1_f64_a = 3.139081648208105e-13
        _check_num_lambdify(expr1, opt1, subs1_a, ref1_f64_a, lambdify_kw={"modules": 'scipy'}, poorness=10**11)

        subs1_b = {x: pi, y: Rational(*(1e-10).as_integer_ratio())}
        ref1_f64_b = 1.1447298859149205e-10
        _check_num_lambdify(expr1, opt1, subs1_b, ref1_f64_b, lambdify_kw={"modules": 'scipy'}, poorness=10**9)


def test_logaddexp_opt():
    x, y = map(Symbol, 'x y'.split())
    expr1 = log(exp(x) + exp(y))
    opt1 = optimize(expr1, [logaddexp_opt])
    assert logaddexp(x, y) - opt1 == 0
    assert logaddexp(y, x) - opt1 == 0
    assert opt1.rewrite(log) == expr1


def test_logaddexp2_opt():
    x, y = map(Symbol, 'x y'.split())
    expr1 = log(2**x + 2**y)/log(2)
    opt1 = optimize(expr1, [logaddexp2_opt])
    assert logaddexp2(x, y) - opt1 == 0
    assert logaddexp2(y, x) - opt1 == 0
    assert opt1.rewrite(log) == expr1


def test_sinc_opts():
    def check(d):
        for k, v in d.items():
            assert optimize(k, sinc_opts) == v

    x = Symbol('x')
    check({
        sin(x)/x       : sinc(x),
        sin(2*x)/(2*x) : sinc(2*x),
        sin(3*x)/x     : 3*sinc(3*x),
        x*sin(x)       : x*sin(x)
    })

    y = Symbol('y')
    check({
        sin(x*y)/(x*y)       : sinc(x*y),
        y*sin(x/y)/x         : sinc(x/y),
        sin(sin(x))/sin(x)   : sinc(sin(x)),
        sin(3*sin(x))/sin(x) : 3*sinc(3*sin(x)),
        sin(x)/y             : sin(x)/y
    })


def test_optims_numpy():
    def check(d):
        for k, v in d.items():
            assert optimize(k, optims_numpy) == v

    x = Symbol('x')
    check({
        sin(2*x)/(2*x) + exp(2*x) - 1: sinc(2*x) + expm1(2*x),
        log(x+3)/log(2) + log(x**2 + 1): log1p(x**2) + log2(x+3)
    })


def test_matsolve_top_level():
    n = Symbol('n', integer=True)
    A = MatrixSymbol('A', n, n)
    x = MatrixSymbol('x', n, 1)
    with assuming(Q.fullrank(A)):
        assert optimize(A**(-1) * x, [matinv_opt]) == MatrixSolve(A, x)
        assert isinstance(optimize(A**(-1) * x, [matinv_opt]), MatrixSolve)
    # Without the full rank assumption the product is kept.
    assert optimize(A**(-1) * x, [matinv_opt]) == A**(-1) * x


def test_rewriting_exposes_reference_names():
    import symox.codegen.rewriting as module
    import sympy.codegen.rewriting as reference
    names = {name for name in vars(reference) if not name.startswith('__')}
    missing = sorted(name for name in names if not hasattr(module, name))
    assert missing == []
    assert [type(o).__name__ for o in module.optims_scipy] == ['FuncMinusOneOptim'] * 2
    assert len(module.optims_numpy) == 9
    assert module.optims_numpy[5:7] == (module.logaddexp_opt, module.logaddexp2_opt)
