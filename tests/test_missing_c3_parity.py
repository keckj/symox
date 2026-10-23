"""Printing and code-generation parity for the c3 batch (values from SymPy 1.14)."""
import numpy as np
import pytest

import symox as S
from symox.tensor import IndexedBase, Idx
from symox.printing.latex import latex

x, y, z = S.symbols('x y z')
n, m = S.symbols('n m', integer=True)


def test_indexed_prints_through_code_printers():
    i = Idx('i', m)
    A = IndexedBase('A')
    yy = IndexedBase('y')
    assert S.fcode(A[i], source_format='free') == 'A(i)'
    assert S.jscode(yy[i] + 1) == 'y[i] + 1'


def test_indexed_contraction_loops():
    i, j = Idx('i', m), Idx('j', n)
    A, B, yy = IndexedBase('A'), IndexedBase('B'), IndexedBase('y')
    code = S.ccode(A[i, j]*B[j], assign_to=yy[i], contract=True)
    # SymPy prints the flat index as n*i + j; symox orders the factors i*n.
    assert code.startswith('for (int i=0; i<m; i++){\n   y[i] = 0;\n}\nfor (int i=0; i<m; i++){\n   for (int j=0; j<n; j++){\n      y[i] = A[')
    assert code.endswith(' + j]*B[j] + y[i];\n   }\n}')
    flat = S.ccode(A[i, j]*B[j], assign_to=yy[i], contract=False)
    assert flat in ('y[i] = A[n*i + j]*B[j];', 'y[i] = A[i*n + j]*B[j];')


def test_latex_containers_and_matrices():
    assert latex((1, x, [y, {z: 1}])) == r'\left( 1, \  x, \  \left[ y, \  \left\{ z : 1\right\}\right]\right)'
    assert latex(S.Matrix([[x, 1], [2, y]])) == r'\left[\begin{matrix}x & 1\\2 & y\end{matrix}\right]'
    assert latex(S.Matrix([[x, 1], [2, y]]), mode='inline') == r'$\left[\begin{smallmatrix}x & 1\\2 & y\end{smallmatrix}\right]$'
    assert latex(S.Float('1.5'), decimal_separator='comma') == '1{,}5'


def test_latex_long_frac_ratio():
    r = S.Symbol('r')
    assert latex(S.Integral(r, r)/2/S.pi, long_frac_ratio=2) == r'\frac{\int r\, dr}{2 \pi}'
    assert latex(S.Integral(r, r)/2/S.pi, long_frac_ratio=0) == r'\frac{1}{2 \pi} \int r\, dr'


def test_rewriting_optims_c99():
    from symox.codegen.rewriting import optimize, optims_c99
    assert str(optimize(S.log(x + 3)/S.log(2) + S.log(x**2 + 1), optims_c99)) == 'log1p(x**2) + log2(x + 3)'
    assert str(optimize(3*S.exp(2*x) - 3, optims_c99)) == '3*expm1(2*x)'
    assert str(optimize(S.log(3*x + 3), optims_c99)) == 'log1p(x) + log(3)'
    assert str(optimize(S.exp(2*x) - 1 - S.exp(-33), optims_c99)) == 'expm1(2*x) - exp(-33)'
    assert str(optimize(S.log(2*x + 3), optims_c99)) == 'log(2*x + 3)'


def test_series_approx():
    from symox.codegen.rewriting import optimize
    from symox.codegen.approximations import SeriesApprox
    expr = S.sin(x)*S.sin(y)
    approx = SeriesApprox(bounds={x: (-.1, .1), y: (S.pi - 1, S.pi + 1)}, reltol=1e-2)
    assert str(optimize(expr, [approx])) == 'x*(-y + (y - pi)**3/6 + pi)'
    with pytest.raises(ValueError):
        SeriesApprox(bounds={x: (0, 1)}, reltol=1e-2, n_point_checks=3)


def test_matrix_solve_printers():
    from symox.codegen.matrix_nodes import MatrixSolve
    from symox.printing.numpy import NumPyPrinter
    M, v = S.MatrixSymbol('M', 3, 3), S.MatrixSymbol('v', 3, 1)
    assert NumPyPrinter().doprint(MatrixSolve(M, v)) == 'numpy.linalg.solve(M, v)'
    assert S.octave_code(MatrixSolve(M, v)) == 'M \\ v'


def test_newtons_method_has_assignment():
    from symox.codegen.algorithms import newtons_method
    from symox.codegen.ast import Assignment
    dx = S.Symbol('dx')
    expr = S.cos(x) - x**3
    algo = newtons_method(expr, x, delta=dx)
    assert algo.has(Assignment(dx, -expr/expr.diff(x)))


def test_utilities_codegen_reference_format():
    from symox.utilities.codegen import codegen, Variable
    assert Variable(S.Symbol('x')).get_datatype('c') == 'double'
    assert Variable(S.Symbol('x')).get_datatype('fortran') == 'REAL*8'
    [(_, c_code), (_, header)] = codegen(('f', x + y*z), 'C89', 'test', header=False, empty=False)
    assert c_code == ('#include "test.h"\n#include <math.h>\ndouble f(double x, double y, double z) {\n'
                      '   double f_result;\n   f_result = x + y*z;\n   return f_result;\n}\n')
    assert header == '#ifndef PROJECT__TEST__H\n#define PROJECT__TEST__H\ndouble f(double x, double y, double z);\n#endif\n'
    [(_, f_code), _] = codegen(('f', x + y*z), 'F95', header=False, empty=False,
                               argument_sequence=(x, y), global_vars=(z,))
    assert f_code == ('REAL*8 function f(x, y)\nimplicit none\nREAL*8, intent(in) :: x\n'
                      'REAL*8, intent(in) :: y\nf = x + y*z\nend function\n')


def test_smtlib_rationals_and_warnings():
    from typing import Callable
    warnings = []
    text = S.smtlib_code(S.sin(x).series(x).removeO(), log_warn=warnings.append)
    assert text == '(declare-const x Real)\n(+ x (* (/ -1 6) (pow x 3)) (* (/ 1 120) (pow x 5)))'
    assert warnings == ['Could not infer type of `x`. Defaulting to float.',
                        'Non-Boolean expression `x**5/120 - x**3/6 + x` will not be asserted. Converting to SMTLib verbatim.']
    f, g = S.Function('f'), S.Function('g')
    text = S.smtlib_code(f(x) + g(x), symbol_table={g: Callable[[int], float]},
                         known_functions={f: 'existing_smtlib_fcn', S.Add: 'sum'})
    assert text == '(declare-const x Int)\n(declare-fun g (Int) Real)\n(sum (existing_smtlib_fcn x) (g x))'


def test_public_decorator():
    from symox.utilities.decorator import public
    namespace = {}
    exec('def f(): pass', namespace)
    public(namespace['f'])
    assert namespace['__all__'] == ['f']
    with pytest.raises(TypeError):
        public(1)


def test_ufuncify_helpers_and_elementwise_backends():
    from symox.utilities.autowrap import ufuncify
    helper = S.Function('helper_func')
    f = ufuncify((x, y), x**2 + y*helper(x), helpers=[('helper_func', x**3, [x])])
    assert list(f(np.array([1., 2.]), np.array([3., 4.]))) == [4., 36.]


def test_plot_implicit_compound_region():
    from symox import plot_implicit
    p = plot_implicit(S.And(y > x, y > -x), show=False)
    X, Y, Z, mode = p[0].get_data()
    assert mode == 'contourf'
    assert np.all(Z == ((Y > X) & (Y > -X)))
    with pytest.raises(NotImplementedError):
        plot_implicit(S.And(S.Eq(y, x), y > 0), show=False)
