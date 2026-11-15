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


# Code generation parity fixes (values observed with SymPy 1.14).

def test_cxx_using_fields_and_new_type_names():
    from symox.codegen import ast
    from symox.codegen.cnodes import static
    from symox.codegen.cxxnodes import using, constexpr
    u = using(ast.real, 'R')
    assert using._fields == ('type', 'alias')
    assert repr(u.alias) == "String('R')" and str(u.type) == 'real'
    assert S.cxxcode(u) == 'using R = double'
    assert S.cxxcode(using('std::vector')) == 'using std::vector'
    assert S.cxxcode(using(ast.real)) == 'using real'
    with pytest.raises(ValueError, match='C\\+\\+98 does not support type aliases'):
        S.cxxcode(u, standard='c++98')
    assert [(str(t.name), t.nbits, t.nmant, t.nexp) for t in (ast.float16, ast.float128, ast.float256)] == [
        ('float16', 16, 10, 5), ('float128', 128, 112, 15), ('float256', 256, 236, 19)]
    assert type(ast.intp).__name__ == 'IntBaseType' and str(ast.intp) == 'intp'
    assert type(static).__name__ == type(constexpr).__name__ == 'Attribute'
    assert S.ccode(ast.Declaration(ast.Variable(x, type=ast.real, attrs=[static]))) == 'double x'
    assert S.ccode(ast.Declaration(ast.Variable(x, type=ast.float128, value=1))) == 'float128 x = 1'


def test_implemented_function_evalf():
    from symox.utilities.lambdify import implemented_function
    h = implemented_function('h_c3_evalf', lambda a, b: a + 2*b)
    assert str(h(2, 3).evalf()) == '8.00000000000000'
    assert str(S.N(h(2, 3), 30)) == '8.00000000000000000000000000000'
    assert str(h(x, 3).subs(x, 2).evalf()) == '8.00000000000000'
    assert str((S.sin(h(2, 3)) + 1).evalf(5)) == '1.9894'
    assert str(h(x, 3).evalf()) == 'h_c3_evalf(x, 3)'
    assert S.lambdify(x, h(x, 3))(2) == 8
    with pytest.raises(ValueError, match='symfunc should be either a string'):
        implemented_function(3, lambda a: a)


def test_pycode_unsupported_nodes_and_c_nodes():
    from symox.codegen import ast
    from symox.codegen.cnodes import goto, Label
    from symox.printing.codeprinter import PrintMethodNotImplementedError
    loop = ast.For(x, range(3), [ast.Assignment(y, x)])
    with pytest.raises(PrintMethodNotImplementedError, match='Unsupported by .*PythonCodePrinter.*For'):
        S.pycode(loop)
    assert S.pycode(loop, strict=False).startswith('  # Not supported in Python:\n  # For\nFor(')
    assert S.ccode(goto('l')) == 'goto l'
    assert repr(goto('l')) == "goto(Label(String('l')))" and repr(Label('l')) == "Label(String('l'))"
    assert S.cxxcode(ast.Raise(ast.RuntimeError_(ast.QuotedString('e')))) == 'throw std::runtime_error("e")'
    assert S.pycode(ast.Raise(ast.RuntimeError_('e'))) == 'raise RuntimeError(e)'
    assert repr(ast.RuntimeError_('e')) == "RuntimeError_(String('e'))"


def test_fortran77_initialized_declaration_raises():
    from symox.codegen import ast
    s = S.Symbol('s')
    fd = ast.FunctionDefinition(ast.real, 'f', [ast.Variable(x, ast.real)],
                                [ast.Declaration(ast.Variable(s, ast.real, 0)), ast.Return(x + s)])
    with pytest.raises(NotImplementedError, match='F77 init./parameter statem. req. multiple lines.'):
        S.fcode(fd, standard=77)
    assert S.fcode(fd, standard=95, source_format='free') == (
        'real*8 function f(x)\nreal*8 :: x\nreal*8 :: s = 0\nf = s + x\nend function')


def test_newtons_method_function_matches_reference():
    from symox.codegen.algorithms import newtons_method, newtons_method_function
    from symox.codegen.pyutils import render_as_module
    from symox.codegen import ast
    expr = S.cos(x) - x**3
    assert S.ccode(newtons_method_function(expr, x)) == (
        'double newton(double x){\n   double d_x = INFINITY;\n'
        '   while (fabs(d_x) > 3.9999999999999999e-16*fabs(x) + 9.9999999999999998e-13) {\n'
        '      d_x = (pow(x, 3) - cos(x))/(-3*pow(x, 2) - sin(x));\n      x += d_x;\n   };\n'
        '   return x;\n}')
    source = render_as_module(newtons_method_function(expr, x))
    assert source == (
        'import math\n\ndef newton(x):\n    d_x = math.inf\n'
        '    while (abs(d_x) > 4.0e-16*abs(x) + 1.0e-12):\n'
        '        d_x = (x**3 - math.cos(x))/(-3*x**2 - math.sin(x))\n        x += d_x\n    return x')
    namespace = {}
    exec(source, namespace, namespace)
    assert abs(namespace['newton'](0.5) - 0.865474033102) < 1e-12
    assert type(newtons_method(expr, x)).__name__ == 'Scope'
    with pytest.raises(ValueError, match='Missing symbols in params: y'):
        newtons_method_function(expr + y, x)
    with pytest.raises(ValueError, match='Only Python 3 is supported.'):
        render_as_module(ast.Assignment(x, 1), standard='python2')


def test_utilities_codegen_matrix_and_indexed_arguments():
    from symox.utilities.codegen import codegen, make_routine
    A, B = S.MatrixSymbol('A', 2, 1), S.MatrixSymbol('B', 2, 2)
    assert codegen(('rot', S.Eq(A, S.Matrix([x*y, x - y]))), 'C99', 'rot', header=False, empty=False) == [
        ('rot.c', '#include "rot.h"\n#include <math.h>\nvoid rot(double x, double y, double *A) {\n'
                  '   A[0] = x*y;\n   A[1] = x - y;\n}\n'),
        ('rot.h', '#ifndef PROJECT__ROT__H\n#define PROJECT__ROT__H\nvoid rot(double x, double y, double *A);\n#endif\n')]
    assert codegen(('mi', B[0, 0] + B[1, 1]*x), 'F95', 'mi', header=False, empty=False)[0] == (
        'mi.f90', 'REAL*8 function mi(B, x)\nimplicit none\nREAL*8, intent(in) :: x\n'
                  'REAL*8, intent(in), dimension(1:2, 1:2) :: B\nmi = x*B(2, 2) + B(1, 1)\nend function\n')
    assert codegen(('iom', S.Eq(A, S.Matrix([A[0, 0] + x, A[1, 0]*y]))), 'Julia', 'iom', header=False, empty=False) == [
        ('iom.jl', 'function iom(A, x, y)\n    A = [x + A[1,1], y .* A[2,1]]\n    return A\nend\n')]
    i = Idx('i', m)
    xi, yi = IndexedBase('x'), IndexedBase('y')
    assert codegen(('l1', S.Eq(yi[i], 2*xi[i] + 1)), 'C99', 'l1', header=False, empty=False)[0] == (
        'l1.c', '#include "l1.h"\n#include <math.h>\nvoid l1(int m, double *x, double *y) {\n'
                '   for (int i=0; i<m; i++){\n      y[i] = 2*x[i] + 1;\n   }\n}\n')
    assert codegen(('l1', S.Eq(yi[i], 2*xi[i] + 1)), 'Octave', 'l1', header=False, empty=False) == [
        ('l1.m', 'function y = l1(m, x)\n  for i = 1:m\n    y(i) = 2*x(i) + 1;\n  end\nend\n')]
    routine = make_routine('io', S.Eq(yi[i], yi[i] + xi[i]))
    assert [(type(a).__name__, str(a.name), str(a.dimensions)) for a in routine.arguments] == [
        ('InputArgument', 'm', 'None'), ('InputArgument', 'x', '[(0, m - 1)]'),
        ('InOutArgument', 'y', '((0, m - 1),)')]


def test_utilities_codegen_banner_names_symox():
    from symox.utilities.codegen import codegen
    banner = codegen(('f', x), 'C99', 'f', project='p')[0][1].splitlines()[1]
    assert banner.split() == ['*', 'Code', 'generated', 'with', 'Symox', S.__version__, '*']
    assert 'SymPy' not in codegen(('f', x), 'F95', 'f')[0][1]


def test_autowrap_removes_its_own_temporary_directory(tmp_path, monkeypatch):
    import tempfile
    from symox.utilities import autowrap as aw
    made = []
    real_mkdtemp = tempfile.mkdtemp
    def mkdtemp(*args, **kwargs):
        made.append(real_mkdtemp(*args, **kwargs))
        return made[-1]
    monkeypatch.setattr(tempfile, 'mkdtemp', mkdtemp)
    aw.autowrap(x + y, backend='dummy')
    assert len(made) == 1 and not __import__('os').path.exists(made[0])
    kept = tmp_path / 'build'
    aw.autowrap(x + y, backend='dummy', tempdir=str(kept))
    assert kept.is_dir() and any(kept.iterdir())
