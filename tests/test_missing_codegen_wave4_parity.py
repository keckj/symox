"""Code generation and printer parity with SymPy 1.14 (wave 4, WP-C14).

The expected values were computed with the reference SymPy 1.14 and are
hard-coded here. ``CORRECTED`` lists the few inputs where symox keeps a
documented correction (see docs/algorithmic-parity.md).
"""
import pytest

SETUP = "\nfrom symox import (symbols, Symbol, Function, Matrix, MatrixSymbol, Piecewise, Rational, sin, cos, exp,\n    log, sqrt, floor, ceiling, Abs, pi, Eq, Lt, oo, nan, Integer, Float, IndexedBase, Idx, sign, coth)\nfrom symox import ccode, fcode, pycode, rust_code, jscode, julia_code, octave_code, rcode, glsl_code\nfrom symox.codegen.ast import (Assignment, AddAugmentedAssignment, CodeBlock, Declaration, Variable,\n    Print, FunctionPrototype, FunctionDefinition, FunctionCall, Return, While, For, real, integer,\n    float32, float64, int8, uint8, complex64, String, QuotedString, Stream, Type, value_const, none, Pointer)\nfrom symox.codegen.cfunctions import (expm1, log1p, exp2, log2, fma, log10, Sqrt, Cbrt, hypot, isnan, isinf)\nfrom symox.codegen.fnodes import (Do, Module, Program, Subroutine, SubroutineCall, use, use_rename,\n    Extent, dimension, intent_in, intent_out, array, bind_C, ArrayConstructor, GoTo, FortranReturn,\n    ImpliedDoLoop, sum_, size, isign, dsign, cmplx, merge, literal_dp, elemental, pure, assumed_extent)\nfrom symox.codegen.cnodes import PreIncrement, PostDecrement, sizeof, alignof, struct\nfrom symox.printing.fortran import FCodePrinter\nfrom symox.printing.rust import RustCodePrinter\nfrom symox.printing.pycode import PythonCodePrinter\nx, y, z, tau = symbols('x y z tau')\nn = Symbol('n', integer=True)\ni = Symbol('i', integer=True)\nr = Symbol('r', real=True)\nA = MatrixSymbol('A', 3, 1)\nmat = Matrix([x**2, Piecewise((x + 1, x > 0), (x, True)), sin(x)])\nM2 = Matrix([[x, y], [x*y, 1]])\n"


@pytest.fixture(scope='module')
def ns():
    namespace = {}
    exec(SETUP, namespace)
    return namespace


def _run(ns, code):
    try:
        value = eval(code, ns)
    except Exception as exc:
        return ['err', type(exc).__name__, str(exc)]
    return ['ok', value if isinstance(value, str) else str(value)]

# (code, expected) from SymPy 1.14; "sympy" module paths read "symox".
REFERENCE = [
    ('expm1(x).fdiff()',
     ['ok', 'exp(x)']),
    ('log1p(x).fdiff()',
     ['ok', '1/(x + 1)']),
    ('exp2(x).fdiff()',
     ['ok', 'log(2)*exp2(x)']),
    ('log2(x).fdiff()',
     ['ok', '1/(x*log(2))']),
    ('fma(x, y, z).fdiff(2)',
     ['ok', 'x']),
    ('log10(x).fdiff()',
     ['ok', '1/(x*log(10))']),
    ('Sqrt(x).fdiff()',
     ['ok', '1/(2*sqrt(x))']),
    ('Cbrt(x).fdiff()',
     ['ok', '1/(3*x**(2/3))']),
    ('hypot(x, y).fdiff(1)',
     ['ok', 'x/hypot(x, y)']),
    ('expm1(x).rewrite(exp)',
     ['ok', 'exp(x) - 1']),
    ('log1p(x).rewrite(log)',
     ['ok', 'log(x + 1)']),
    ('exp2(x).rewrite(exp)',
     ['ok', 'exp2(x)']),
    ('log2(x).rewrite(log)',
     ['ok', 'log(x)/log(2)']),
    ("fma(x, y, z).rewrite(Function('Add'))",
     ['ok', 'fma(x, y, z)']),
    ('hypot(x, y).rewrite(sqrt)',
     ['ok', 'hypot(x, y)']),
    ('expm1(Rational(1, 3)).evalf(30)',
     ['ok', '0.395612425086089528628125319603']),
    ('log1p(Rational(1, 10**20)).evalf(25)',
     ['ok', '9.999999999999999999950000e-21']),
    ('isnan(nan)',
     ['ok', 'True']),
    ('isnan(Integer(1))',
     ['ok', 'False']),
    ('isinf(oo)',
     ['ok', 'True']),
    ('isinf(Integer(1))',
     ['ok', 'False']),
    ('isnan(x)',
     ['ok', 'isnan(x)']),
    ('ccode(expm1(x))',
     ['ok', 'expm1(x)']),
    ('ccode(log1p(x) + exp2(y) + log2(z))',
     ['ok', 'exp2(y) + log1p(x) + log2(z)']),
    ('ccode(fma(x, y, z))',
     ['ok', 'fma(x, y, z)']),
    ("ccode(hypot(x, y), standard='c99')",
     ['ok', 'hypot(x, y)']),
    ('ccode(Cbrt(x))',
     ['ok', 'cbrt(x)']),
    ('ccode(isnan(x))',
     ['ok', 'isnan(x)']),
    ('pycode(log1p(x))',
     ['ok', 'math.log1p(x)']),
    ('pycode(Sqrt(x))',
     ['ok', 'math.sqrt(x)']),
    ('pycode(isinf(x))',
     ['ok', 'math.isinf(x)']),
    ('pycode(Cbrt(x))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.pycode.PythonCodePrinter'>: Cbrt\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ('fcode(sin(x)**2 + Rational(1, 3)*y)',
     ['ok', '      (1.0d0/3.0d0)*y + sin(x)**2']),
    ('fcode(Piecewise((x, x < 1), (x**2, True)), tau)',
     ['ok', '      if (x < 1) then\n         tau = x\n      else\n         tau = x**2\n      end if']),
    ("fcode(Piecewise((x, x < 1), (x**2, True)), tau, standard=95, source_format='free')",
     ['ok', 'if (x < 1) then\n   tau = x\nelse\n   tau = x**2\nend if']),
    ("fcode(floor(x) + ceiling(x), user_functions={'floor': 'FLOOR1', 'ceiling': 'CEIL'})",
     ['ok', '      CEIL(x) + FLOOR1(x)']),
    ('fcode(mat, A)',
     ['ok', '      A(1, 1) = x**2\n      if (x > 0) then\n         A(2, 1) = x + 1\n      else\n         A(2, 1) = x\n      end if\n      A(3, 1) = sin(x)']),
    ('fcode(2*x + pi, tau)',
     ['ok', '      parameter (pi = 3.1415926535897932d0)\n      tau = 2*x + pi']),
    ("fcode(sqrt(1 - x**2).series(x, n=20).removeO(), source_format='free')",
     ['ok', '-715.0d0/65536.0d0*x**18 - 429.0d0/32768.0d0*x**16 - 33.0d0/2048.0d0*x** &\n      14 - 21.0d0/1024.0d0*x**12 - 7.0d0/256.0d0*x**10 - 5.0d0/128.0d0* &\n      x**8 - 1.0d0/16.0d0*x**6 - 1.0d0/8.0d0*x**4 - 1.0d0/2.0d0*x**2 + &\n      1']),
    ('fcode(exp(-x))',
     ['ok', '      exp(-x)']),
    ('fcode(coth(x), standard=95)',
     ['ok', '      ((exp(x) + exp(-x))/(exp(x) - exp(-x)))']),
    ("fcode(sign(x), standard=95, source_format='free')",
     ['ok', 'merge(0d0, dsign(1d0, x), x == 0d0)']),
    ("fcode(sign(r), standard=95, source_format='free')",
     ['ok', 'merge(cmplx(0d0, 0d0), r/abs(r), abs(r) == 0d0)']),
    ('fcode(sign(x))',
     ['err', 'NotImplementedError', 'merge requires Fortran 95 or newer']),
    ("fcode(1/n, source_format='free')",
     ['ok', '1d0/n']),
    ('fcode(oo)',
     ['ok', '      (huge(0d0) + 1)']),
    ("fcode(Function('g')(x))",
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.fortran.FCodePrinter'>: g\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ("fcode(Function('g')(x), strict=False)",
     ['ok', 'C     Not supported in Fortran:\nC     g\n      g(x)']),
    ('fcode(x, standard=11)',
     ['err', 'ValueError', 'Unknown Fortran standard: 11']),
    ("fcode(isign(1, n), source_format='free')",
     ['ok', 'isign(1, n)']),
    ("fcode(merge(x, y, Eq(x, 1)), source_format='free')",
     ['err', 'NotImplementedError', 'merge requires Fortran 95 or newer']),
    ('fcode(merge(x, y, Eq(x, 1)))',
     ['err', 'NotImplementedError', 'merge requires Fortran 95 or newer']),
    ("fcode(cmplx(x, y), source_format='free')",
     ['ok', 'cmplx(x, y)']),
    ("fcode(literal_dp(1)/x, source_format='free')",
     ['ok', '1.0d0/x']),
    ('repr(literal_dp(1))',
     ['ok', '1.00000000000000']),
    ("fcode(Do([AddAugmentedAssignment(x, 1/n), Print([n, x])], n, 1, 10), source_format='free')",
     ['ok', 'do n = 1, 10\n   x = x + 1d0/n\n   print *, n, x\nend do']),
    ("fcode(Do([Print([i])], i, 1, 9, 2), source_format='free')",
     ['ok', 'do i = 1, 9, 2\n   print *, i\nend do']),
    ("fcode(Module('signallib', ['implicit none'], []), source_format='free')",
     ['ok', 'module signallib\nimplicit none\n\ncontains\n\n\nend module']),
    ("fcode(Program('myprog', [Print([x])]), source_format='free')",
     ['ok', 'program myprog\n   print *, x\nend program']),
    ("fcode(use('signallib'), source_format='free')",
     ['ok', 'use signallib']),
    ("fcode(use('signallib', [('metric', 'snr')]), source_format='free')",
     ['ok', 'use signallib, metric => snr']),
    ("fcode(use('signallib', only=['snr', 'convolution2d']), source_format='free')",
     ['ok', 'use signallib, only: snr, convolution2d']),
    ("fcode(use('signallib', only=['snr', use_rename('thingy', 'convolution2d')]), source_format='free')",
     ['ok', 'use signallib, only: snr, thingy => convolution2d']),
    ("fcode(Extent(-3, 3), source_format='free')",
     ['ok', '-3:3']),
    ("fcode(Extent(), source_format='free')",
     ['ok', ':']),
    ("fcode(Variable(x, real, attrs=[dimension(3, 3)]).as_Declaration(), source_format='free')",
     ['err', 'NotImplementedError', 'F77 init./parameter statem. req. multiple lines.']),
    ("fcode(Variable(x, real, attrs=[intent_in, dimension(assumed_extent)]).as_Declaration(), source_format='free')",
     ['err', 'NotImplementedError', 'F77 init./parameter statem. req. multiple lines.']),
    ("fcode(array(x, (3,), 'in'), source_format='free')",
     ['err', 'ValueError', 'Could not deduce type from expr.']),
    ("fcode(ArrayConstructor([1, 2, 3]), source_format='free')",
     ['ok', '(/1, 2, 3/)']),
    ("fcode(ArrayConstructor([ImpliedDoLoop(i**2, i, 1, 3)]), source_format='free')",
     ['ok', '(/(i**2, i = 1, 3)/)']),
    ("fcode(GoTo([10, 20], n), source_format='free')",
     ['ok', 'go to (10, 20), n']),
    ("fcode(FortranReturn(x), source_format='free')",
     ['ok', 'return x']),
    ("fcode(SubroutineCall('rout', [x, y]), source_format='free')",
     ['ok', 'call rout(x, y)']),
    ("fcode(FunctionDefinition(real, 'f', [Variable(x, real)], [Return(x)], attrs=[bind_C('f')]), source_format='free', standard=2003)",
     ['ok', 'real*8 function f(x) bind(C, name="f")\nreal*8 :: x\nf = x\nend function']),
    ("fcode(While(x > 0, [Assignment(x, x - 1)]), source_format='free')",
     ['ok', 'do while (x > 0)\n   x = x - 1\nend do']),
    ("fcode(For(i, range(3), [Print([i])]), source_format='free')",
     ['ok', 'do i = 0, 2, 1\n   print *, i\nend do']),
    ('fcode(Print([x, y], "fmt"), source_format=\'free\')',
     ['ok', 'write(*, fmt=""fmt"", advance="no"), x, y']),
    ("fcode(Declaration(Variable(x, real, 3.0, attrs=[value_const])), source_format='free')",
     ['err', 'NotImplementedError', 'F77 init./parameter statem. req. multiple lines.']),
    ('fcode(Declaration(Variable(x, real, 3.0, attrs=[value_const])))',
     ['err', 'NotImplementedError', 'F77 init./parameter statem. req. multiple lines.']),
    ("fcode(sum_(Symbol('a'))*2, source_format='free')",
     ['ok', '2*sum(a)']),
    ("FCodePrinter({'source_format': 'free'}).doprint(size(Symbol('a')))",
     ['ok', 'size(a)']),
    ('rust_code(mat, A)',
     ['ok', 'A = [x.powi(2), if (x > 0.0) {\n    x + 1\n} else {\n    x\n}, x.sin()];']),
    ('rust_code(Piecewise((x, x < 1), (x**2, True)), tau)',
     ['ok', 'tau = if (x < 1.0) {\n    x\n} else {\n    x.powi(2)\n};']),
    ("rust_code(Abs(x) + ceiling(x), user_functions={'ceiling': 'CEIL'})",
     ['ok', 'x.abs() + x.ceil()']),
    ('rust_code(floor(x) + ceiling(x) + Abs(x))',
     ['ok', 'x.abs() + x.ceil() + x.floor()']),
    ('rust_code(x**3 + sqrt(x) + 1/x)',
     ['ok', 'x.sqrt() + x.powi(3) + x.recip()']),
    ('rust_code(x**2 + sin(x) + 1)',
     ['ok', 'x.powi(2) + x.sin() + 1']),
    ('rust_code(isnan(x))',
     ['ok', 'isnan(x)']),
    ('rust_code(Eq(x, y))',
     ['ok', 'x == y']),
    ("RustCodePrinter().doprint(sin(x), 'y')",
     ['ok', 'y = x.sin();']),
    ('jscode(mat, A)',
     ['ok', 'A[0] = Math.pow(x, 2);\nif (x > 0) {\n   A[1] = x + 1;\n}\nelse {\n   A[1] = x;\n}\nA[2] = Math.sin(x);']),
    ('jscode(Piecewise((x, x < 1), (x**2, True)), tau)',
     ['ok', 'if (x < 1) {\n   tau = x;\n}\nelse {\n   tau = Math.pow(x, 2);\n}']),
    ('julia_code(mat, A)',
     ['ok', 'A = [x .^ 2, ((x > 0) ? (x + 1) : (x)), sin(x)]']),
    ('julia_code(M2)',
     ['ok', '[     x y;\nx .* y 1]']),
    ('julia_code(M2, inline=False)',
     ['ok', '[     x y;\nx .* y 1]']),
    ("julia_code(x**2*y + 1, assign_to='s')",
     ['ok', 's = x .^ 2 .* y + 1']),
    ('octave_code(mat, A)',
     ['ok', 'A = [x.^2; ((x > 0).*(x + 1) + (~(x > 0)).*(x)); sin(x)];']),
    ('octave_code(M2)',
     ['ok', '[x y; x.*y 1]']),
    ('octave_code(Piecewise((x, x < 1), (x**2, True)), assign_to=tau)',
     ['ok', 'tau = ((x < 1).*(x) + (~(x < 1)).*(x.^2));']),
    ('rcode(mat, A)',
     ['ok', 'A[0] = x^2;\nA[1] = ifelse(x > 0,x + 1,x);\nA[2] = sin(x);']),
    ('rcode(Piecewise((x, x < 1), (x**2, True)), tau)',
     ['ok', 'tau = ifelse(x < 1,x,x^2);']),
    ('glsl_code(mat, A)',
     ['ok', 'A[0][0] = pow(x, 2.0);\nif (x > 0) {\n   A[1][0] = x + 1;\n}\nelse {\n   A[1][0] = x;\n}\nA[2][0] = sin(x);']),
    ('glsl_code(M2)',
     ['ok', 'mat2(  x, y, x*y, 1)']),
    ('glsl_code(Matrix([1, 2, 3, 4]))',
     ['ok', 'vec4(1, 2, 3, 4)']),
    ('glsl_code(M2, mat_nested=True)',
     ['ok', 'mat2(  x, y, x*y, 1)']),
    ('glsl_code(Matrix([1, 2, 3, 4]), glsl_types=False)',
     ['ok', 'float[4](1, 2, 3, 4)']),
    ('ccode(PreIncrement(x))',
     ['ok', '++(x)']),
    ('ccode(PostDecrement(x))',
     ['ok', '(x)--']),
    ("ccode(sizeof('unsigned int'))",
     ['ok', 'sizeof(unsigned int)']),
    ("ccode(alignof('int'))",
     ['ok', 'alignof(int)']),
    ('ccode(Print([x, y], "coordinate: %12.5g %12.5g\\n"))',
     ['ok', 'printf("coordinate: %12.5g %12.5g\n", x, y)']),
    ('pycode(Print([x, y], "coordinate: %12.5g %12.5g\\n"))',
     ['ok', 'print("coordinate: %12.5g %12.5g\n" % (x, y), end="")']),
    ('pycode(Print([x, y]))',
     ['ok', 'print(x, y)']),
    ("pycode(FunctionDefinition(real, 'f', [Variable(x, real)], [Return(x**2)]))",
     ['ok', 'def f(x):\n    return x**2']),
    ('pycode(While(x > 0, [Assignment(x, x - 1)]))',
     ['ok', 'while (x > 0):\n    x = x - 1']),
    ("ccode(FunctionPrototype(real, 'f', [Variable(x, real), Variable(n, integer)]))",
     ['ok', 'double f(double x, int n)']),
    ("ccode(FunctionDefinition(real, 'f', [Variable(x, real)], [Return(x + 1)]))",
     ['ok', 'double f(double x){\n   return x + 1;\n}']),
    ('ccode(Declaration(Variable(x, float32, 3, attrs=[value_const])))',
     ['ok', 'const float x = 3']),
    ('ccode(Variable(x, float64, pointer=True).as_Declaration() if False else Declaration(Pointer(x, real)))',
     ['ok', 'double * x']),
    ('str(Variable.deduced(n).type)',
     ['ok', 'integer']),
    ('str(Variable.deduced(x).type)',
     ['err', 'ValueError', 'Could not deduce type from expr.']),
    ('str(Variable.deduced(r).type)',
     ['ok', 'real']),
    ('str(Variable.deduced(x, 3.0).type)',
     ['ok', 'real']),
    ('str(Type.from_expr(Integer(3)))',
     ['ok', 'integer']),
    ('str(Type.from_expr(Float(3)))',
     ['ok', 'real']),
    ('str(Type.from_expr(x))',
     ['err', 'ValueError', 'Could not deduce type from expr.']),
    ('int8.cast_check(Integer(127))',
     ['ok', '127']),
    ('int8.cast_check(Integer(128))',
     ['err', 'ValueError', 'Value is too big: 128 > 127']),
    ('uint8.cast_check(Integer(-1))',
     ['err', 'ValueError', 'Value is too small: -1 < 0']),
    ("float32.cast_check(Float('1e40'))",
     ['err', 'ValueError', 'Value is too big: 9999999999733643722579435254086528663552 > 340282346638528859811704183484516925440']),
    ('float32.cast_check(Integer(3))',
     ['ok', '3.00000000']),
    ("complex64.cast_check(1 + 2*Symbol('I'))",
     ['err', 'ValueError', 'string-float not recognized: 2.0*re(I) + 1.0']),
    ('float64.cast_check(Rational(1, 3))',
     ['ok', '0.33333333333333333']),
    ('CodeBlock(Assignment(x, 1), Assignment(y, x + 1)).left_hand_sides',
     ['ok', '(x, y)']),
    ('CodeBlock(Assignment(x, 1), Assignment(y, x + 1)).right_hand_sides',
     ['ok', '(1, x + 1)']),
    ('CodeBlock(Assignment(y, x + 1), Assignment(x, 1)).topological_sort([Assignment(y, x + 1), Assignment(x, 1)]).left_hand_sides',
     ['ok', '(x, y)']),
    ('ccode(CodeBlock(Assignment(x, 1), Assignment(z, (x + y)**2), Assignment(tau, (x + y)**2 + 1)).cse())',
     ['ok', 'x = 1;\nx0 = pow(x + y, 2);\nz = x0;\ntau = x0 + 1;']),
    ('str(Print([x]).format_string)',
     ['ok', 'NoneToken()']),
]

# Same exception type as SymPy; the message names a symox class.
SAME_ERROR = [
    ('jscode(M2)', 'PrintMethodNotImplementedError'),
    ('rcode(M2)', 'PrintMethodNotImplementedError'),
]

# Documented corrections: (code, sympy result, symox result).
CORRECTED = [
    ("fcode(Do([], i, 1, 3, concurrent=True), source_format='free')",
     ['ok', 'do concurrent i = 1, 3\n\nend do'],
     ['ok', 'do concurrent (i = 1:3)\n\nend do']),
    ("fcode(GoTo([10]), source_format='free')",
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'sympy.printing.fortran.FCodePrinter'>: <class 'sympy.codegen.ast.NoneToken'>\nSet the printer option 'strict' to False in order to generate partially printed code."],
     ['ok', 'go to 10']),
    ("fcode(Subroutine('rout', [Variable(x, real, attrs=[intent_out])], [Assignment(x, 1)]), source_format='free')",
     ['err', 'NotImplementedError', 'F77 init./parameter statem. req. multiple lines.'],
     ['ok', 'subroutine rout(x)\nreal*8 x\nx = 1\nend subroutine']),
    ("fcode(FunctionDefinition(real, 'f', [Variable(x, real, attrs=[intent_in])], [Return(x**2)], attrs=[elemental]), source_format='free')",
     ['err', 'NotImplementedError', 'F77 init./parameter statem. req. multiple lines.'],
     ['ok', 'elemental real*8 function f(x)\nreal*8 x\nf = x**2\nend function']),
    ("fcode(FunctionPrototype(real, 'f', [Variable(x, real)]), source_format='free')",
     ['err', 'NotImplementedError', 'F77 init./parameter statem. req. multiple lines.'],
     ['ok', 'interface\n   real*8 function f(x)\n   real*8 x\n   end function\nend interface']),
]


@pytest.mark.parametrize('code,expected', REFERENCE)
def test_reference_outputs(ns, code, expected):
    assert _run(ns, code) == expected


@pytest.mark.parametrize('code,error', SAME_ERROR)
def test_reference_errors(ns, code, error):
    result = _run(ns, code)
    assert result[:2] == ['err', error]


@pytest.mark.parametrize('code,reference,expected', CORRECTED)
def test_documented_corrections(ns, code, reference, expected):
    assert _run(ns, code) == expected


def test_isnan_substituted_evalf_is_false(ns):
    # SymPy 1.14 raises AttributeError ('BooleanFalse' object has no
    # attribute 'evalf'); symox returns the evaluated predicate.
    from symox import N, Rational, Symbol, false
    from symox.codegen.cfunctions import isnan, isinf
    x = Symbol('x')
    assert N(isnan(x).subs({x: Rational(3, 7)})) == false
    assert N(isinf(x).subs({x: Rational(3, 7)})) == false


def test_pycode_unsupported_cfunction_error_type(ns):
    from symox.printing.codeprinter import PrintMethodNotImplementedError
    with pytest.raises(PrintMethodNotImplementedError, match='Unsupported by .*PythonCodePrinter.*: Cbrt'):
        eval('pycode(Cbrt(x))', ns)


# codegen.numpy_nodes, scipy_nodes, pynodes, their printers and the
# rewriting optimizations producing them; (code, expected) from SymPy 1.14.
NODES_SETUP = "\nfrom symox import (Pow, symbols, Symbol, Float, Rational, Integer, I, pi, oo, sin, cos, exp, log, sqrt, Max, Min,\n    Derivative, Tuple, srepr, latex, N, pycode, ccode, cxxcode, simplify, MatrixSymbol, Q, assuming)\nfrom symox.codegen.numpy_nodes import logaddexp, logaddexp2, amin, amax, minimum, maximum\nfrom symox.codegen.scipy_nodes import cosm1, powm1\nfrom symox.codegen.pynodes import List, NumExprEvaluate\nfrom symox.codegen.ast import Assignment, Return, CodeBlock\nfrom symox.codegen.cfunctions import log2, expm1, log1p\nfrom symox.codegen import rewriting as rw\nfrom symox.codegen.rewriting import optimize, optims_c99, optims_numpy, optims_scipy\nfrom symox.printing.pycode import PythonCodePrinter, MpmathPrinter\nfrom symox.printing.numpy import NumPyPrinter, SciPyPrinter\nfrom symox.printing.lambdarepr import NumExprPrinter, LambdaPrinter\nfrom symox.printing.c import C99CodePrinter\nfrom symox.printing.cxx import CXX17CodePrinter\nx, y, z = symbols('x y z')\n"

NODES_REFERENCE = [
    ('logaddexp(y, x)',
     ['ok', 'logaddexp(x, y)']),
    ('logaddexp2(y, x)',
     ['ok', 'logaddexp2(x, y)']),
    ('logaddexp(Float(1.5), Float(2.5))',
     ['ok', '2.81326168751822']),
    ("logaddexp(Float('0.5', 30), 1)",
     ['ok', 'logaddexp(0.5, 1)']),
    ('logaddexp2(Float(1), Float(3))',
     ['ok', '3.32192809488736']),
    ('cosm1(Float(0.5))',
     ['ok', '-0.122417438109627']),
    ('powm1(Float(2), Float(0.5))',
     ['ok', '0.414213562373095']),
    ('cosm1(0)',
     ['ok', 'cosm1(0)']),
    ('powm1(x, 0)',
     ['ok', 'powm1(x, 0)']),
    ('minimum(x, y)',
     ['ok', 'minimum(x, y)']),
    ('maximum(Float(1), Float(2))',
     ['ok', 'maximum(1.0, 2.0)']),
    ('logaddexp(x, y).args',
     ['ok', '(x, y)']),
    ('logaddexp(x, y, evaluate=False)',
     ['err', 'TypeError', "logaddexp.__new__() got an unexpected keyword argument 'evaluate'"]),
    ('amin(x)',
     ['ok', 'amin(x)']),
    ('amax(x, axis=1)',
     ['ok', 'amax(x, axis=1)']),
    ('amin(x).axis',
     ['ok', 'NoneToken()']),
    ('amax(x, 2).axis',
     ['ok', '2']),
    ('amin([1, 2])',
     ['ok', 'amin([1, 2])']),
    ('amin(x) == amin(x)',
     ['ok', 'True']),
    ('List(1, 2, x)',
     ['ok', '[1, 2, x]']),
    ('List(1, 2, x).func',
     ['ok', "<class 'symox.codegen.pynodes.List'>"]),
    ('type(List(x).subs(x, 2)).__name__',
     ['ok', 'List']),
    ('List(x, y)[1]',
     ['ok', 'y']),
    ('len(List(x, y))',
     ['ok', '2']),
    ('NumExprEvaluate(x + 1)',
     ['ok', 'NumExprEvaluate(x + 1)']),
    ('NumExprEvaluate(x).expr',
     ['ok', 'x']),
    ('logaddexp(x, y).diff(x)',
     ['ok', '1/(exp(-x + y) + 1)']),
    ('logaddexp(x, y).diff(y)',
     ['ok', '1/(exp(x - y) + 1)']),
    ('logaddexp2(x, y).diff(x)',
     ['ok', '1/(2**(-x + y) + 1)']),
    ('logaddexp2(x, y).diff(y)',
     ['ok', '1/(2**(x - y) + 1)']),
    ('cosm1(x).diff(x)',
     ['ok', '-sin(x)']),
    ('powm1(x, y).diff(x)',
     ['ok', 'x**y*y/x']),
    ('powm1(x, y).diff(y)',
     ['ok', 'x**y*log(x)']),
    ('logaddexp(x, y).fdiff(3)',
     ['err', 'ArgumentIndexError', 'Invalid operation with argument number 3 for Function logaddexp(x, y)']),
    ('cosm1(x).fdiff(2)',
     ['err', 'ArgumentIndexError', 'Invalid operation with argument number 2 for Function cosm1(x)']),
    ('logaddexp(x**2, x).diff(x)',
     ['ok', '2*x/(exp(-x**2 + x) + 1) + 1/(exp(x**2 - x) + 1)']),
    ('logaddexp(x, y).rewrite(log)',
     ['ok', 'log(exp(x) + exp(y))']),
    ('logaddexp2(x, y).rewrite(log)',
     ['ok', 'log(2**x + 2**y)/log(2)']),
    ('cosm1(x).rewrite(cos)',
     ['ok', 'cos(x) - 1']),
    ('powm1(x, y).rewrite(Pow)',
     ['ok', 'x**y - 1']),
    ('maximum(x, y).rewrite(Max)',
     ['ok', 'Max(x, y)']),
    ('minimum(x, y, z).rewrite(Min)',
     ['ok', 'Min(x, y, z)']),
    ('logaddexp(log(2), log(3)).evalf()',
     ['ok', '1.60943791243410']),
    ('N(logaddexp(1, 2), 30)',
     ['ok', '2.31326168751822283404899549497']),
    ('cosm1(Rational(1, 3)).evalf(20)',
     ['ok', '-0.055043053685262335612']),
    ('powm1(2, Rational(1, 2)).evalf()',
     ['ok', '0.414213562373095']),
    ('logaddexp(log(2), log(3)).simplify()',
     ['ok', 'log(5)']),
    ('logaddexp2(1, 1).simplify()',
     ['ok', 'log(4)/log(2)']),
    ('cosm1(x).simplify()',
     ['ok', 'cos(x) - 1']),
    ('powm1(x, y).simplify()',
     ['ok', 'x**y - 1']),
    ('simplify(cosm1(x) + 1)',
     ['ok', 'cosm1(x) + 1']),
    ('simplify(powm1(x, y) + 1)',
     ['ok', 'powm1(x, y) + 1']),
    ('logaddexp(x, y).is_real',
     ['ok', 'None']),
    ('cosm1(x).is_number',
     ['ok', 'False']),
    ('srepr(logaddexp(x, y))',
     ['ok', "logaddexp(Symbol('x'), Symbol('y'))"]),
    ('srepr(cosm1(x))',
     ['ok', "cosm1(Symbol('x'))"]),
    ('srepr(amax(x, axis=1))',
     ['ok', "amax(Symbol('x'), axis=Integer(1))"]),
    ('srepr(List(x, 1))',
     ['ok', "List(Symbol('x'), Integer(1))"]),
    ('str(amin(x))',
     ['ok', 'amin(x)']),
    ('str(NumExprEvaluate(x + 1))',
     ['ok', 'NumExprEvaluate(x + 1)']),
    ('latex(logaddexp(x, y))',
     ['ok', '\\operatorname{logaddexp}{\\left(x,y \\right)}']),
    ('latex(cosm1(x))',
     ['ok', '\\operatorname{cosm}_{1}{\\left(x \\right)}']),
    ('latex(powm1(x, y))',
     ['ok', '\\operatorname{powm}_{1}{\\left(x,y \\right)}']),
    ('latex(maximum(x, y))',
     ['ok', '\\operatorname{maximum}{\\left(x,y \\right)}']),
    ('pycode(logaddexp(x, y))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.pycode.PythonCodePrinter'>: logaddexp\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ('pycode(cosm1(x))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.pycode.PythonCodePrinter'>: cosm1\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ('pycode(cosm1(x), strict=False)',
     ['ok', '  # Not supported in Python:\n  # cosm1\ncosm1(x)']),
    ('pycode(logaddexp(x, y) + cosm1(x), strict=False)',
     ['ok', '  # Not supported in Python:\n  # cosm1\n  # logaddexp\ncosm1(x) + logaddexp(x, y)']),
    ('pycode(List(x, 1))',
     ['ok', '[x, 1]']),
    ('pycode(amin(x))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.pycode.PythonCodePrinter'>: <class 'symox.codegen.numpy_nodes.amin'>\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ('pycode(amin(x), strict=False)',
     ['ok', "  # Not supported in Python:\n  # amin\namin(Symbol('x'))"]),
    ('PythonCodePrinter().doprint(cosm1(x))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.pycode.PythonCodePrinter'>: cosm1\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ("PythonCodePrinter({'strict': False}).doprint(powm1(x, y))",
     ['ok', '  # Not supported in Python:\n  # powm1\npowm1(x, y)']),
    ('MpmathPrinter().doprint(logaddexp(x, y))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.pycode.MpmathPrinter'>: logaddexp\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ("MpmathPrinter({'strict': False}).doprint(logaddexp(x, y))",
     ['ok', '  # Not supported in Python with mpmath:\n  # logaddexp\nlogaddexp(x, y)']),
    ('NumPyPrinter().doprint(logaddexp(x, y))',
     ['ok', 'numpy.logaddexp(x, y)']),
    ('NumPyPrinter().doprint(logaddexp2(x, y))',
     ['ok', 'numpy.logaddexp2(x, y)']),
    ('NumPyPrinter().doprint(minimum(x, y, z))',
     ['ok', 'functools.reduce(numpy.minimum, [x, y, z])']),
    ('NumPyPrinter().doprint(maximum(x))',
     ['ok', 'x']),
    ('NumPyPrinter().doprint(maximum())',
     ['err', 'NotImplementedError', 'Need at least one argument for numpy.maximum']),
    ('NumPyPrinter().doprint(amin(x, axis=0))',
     ['ok', 'numpy.amin(x, axis=0)']),
    ('NumPyPrinter().doprint(amax(x))',
     ['ok', 'numpy.amax(x, axis=None)']),
    ('NumPyPrinter().doprint(cosm1(x))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.numpy.NumPyPrinter'>: cosm1\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ("NumPyPrinter({'strict': False}).doprint(cosm1(x))",
     ['ok', '  # Not supported in Python with numpy:\n  # cosm1\ncosm1(x)']),
    ('NumPyPrinter().doprint(Max(x, y, z))',
     ['ok', 'functools.reduce(numpy.maximum, [x, y, z])']),
    ('NumPyPrinter().doprint(Min(x, y))',
     ['ok', 'functools.reduce(numpy.minimum, [x, y])']),
    ('SciPyPrinter().doprint(cosm1(x))',
     ['ok', 'scipy.special.cosm1(x)']),
    ('SciPyPrinter().doprint(powm1(x, y))',
     ['ok', 'scipy.special.powm1(x, y)']),
    ('SciPyPrinter().doprint(logaddexp(x, y))',
     ['ok', 'numpy.logaddexp(x, y)']),
    ('SciPyPrinter().doprint(pi)',
     ['ok', 'scipy.constants.pi']),
    ('SciPyPrinter().doprint(Derivative(cosm1(x), x))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.numpy.SciPyPrinter'>: <class 'symox.core.function.Derivative'>\nPrinter has no method: _print_Derivative_cosm1\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ("SciPyPrinter({'strict': False}).doprint(Derivative(cosm1(x), x))",
     ['ok', '  # Not supported in Python with SciPy and NumPy:\n  # Derivative\nDerivative(cosm1(x), x)']),
    ('LambdaPrinter().doprint(cosm1(x))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.lambdarepr.LambdaPrinter'>: cosm1\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ('NumExprPrinter().doprint(x + 1)',
     ['ok', "numexpr.evaluate('x + 1', truediv=True)"]),
    ('NumExprPrinter().doprint(sin(x) + I)',
     ['ok', "numexpr.evaluate('sin(x) + 1j', truediv=True)"]),
    ('NumExprPrinter().doprint(cosm1(x))',
     ['err', 'TypeError', "numexpr does not support function 'cosm1'"]),
    ('NumExprPrinter().doprint(Return(x + 1))',
     ['ok', "return numexpr.evaluate('x + 1', truediv=True)"]),
    ('NumExprPrinter().doprint([x])',
     ['err', 'TypeError', 'numexpr cannot be used with list']),
    ('NumExprPrinter().doprint(NumExprEvaluate(x**2))',
     ['ok', "numexpr.evaluate('x**2', truediv=True)"]),
    ('NumExprPrinter().doprint(Assignment(y, x + 1))',
     ['ok', "y = numexpr.evaluate('x + 1', truediv=True)"]),
    ('ccode(logaddexp(x, y))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.c.C99CodePrinter'>: logaddexp\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ('ccode(cosm1(x), strict=False)',
     ['ok', '/* Not supported in C: */\n/* cosm1 */\ncosm1(x)']),
    ('cxxcode(powm1(x, y))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.cxx.CXX11CodePrinter'>: powm1\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ('ccode(amin(x))',
     ['err', 'PrintMethodNotImplementedError', "Unsupported by <class 'symox.printing.c.C99CodePrinter'>: <class 'symox.codegen.numpy_nodes.amin'>\nSet the printer option 'strict' to False in order to generate partially printed code."]),
    ("C99CodePrinter({'strict': False}).doprint(logaddexp(x, y))",
     ['ok', '/* Not supported in C: */\n/* logaddexp */\nlogaddexp(x, y)']),
    ("CXX17CodePrinter({'strict': False}).doprint(amax(x))",
     ['ok', "/* Not supported in C++: */\n/* amax */\namax(Symbol('x'))"]),
    ('optimize(cos(x) - 1, [rw.cosm1_opt])',
     ['ok', 'cosm1(x)']),
    ('optimize(3*cos(x) - 5, [rw.cosm1_opt])',
     ['ok', '3*cosm1(x) - 2']),
    ('optimize(2 - 2*cos(x), [rw.cosm1_opt])',
     ['ok', '-2*cosm1(x)']),
    ('optimize(cos(x) + cos(y) - 2, [rw.cosm1_opt])',
     ['ok', 'cosm1(x) + cosm1(y)']),
    ('optimize(exp(x) + cos(x) - 2, [rw.expm1_opt, rw.cosm1_opt])',
     ['ok', 'cosm1(x) + expm1(x)']),
    ('optimize(1/cos(x) - 1, [rw.cosm1_opt])',
     ['ok', '-cosm1(x)/cos(x)']),
    ('optimize(2/cos(x) - 2, optims_scipy)',
     ['ok', '-2*cosm1(x)/cos(x)']),
    ('optimize(x**y - 1, [rw.powm1_opt])',
     ['ok', 'powm1(x, y)']),
    ('optimize(2**x - 1, [rw.powm1_opt])',
     ['ok', 'powm1(2, x)']),
    ('optimize(log(exp(x) + exp(y)), optims_numpy)',
     ['ok', 'logaddexp(x, y)']),
    ('optimize(log(2**x + 2**y), optims_numpy)',
     ['ok', 'log(exp2(x) + exp2(y))']),
    ('optimize(log(exp(x) + exp(y)), [rw.logaddexp_opt])',
     ['ok', 'logaddexp(x, y)']),
    ('optimize(log(2**x + 2**y), [rw.logaddexp2_opt])',
     ['ok', 'log(2)*logaddexp2(x, y)']),
    ('optimize(log(exp(2*x) + exp(3)), [rw.logaddexp_opt])',
     ['ok', 'logaddexp(3, 2*x)']),
    ('optimize(sin(x)/x, optims_numpy)',
     ['ok', 'sinc(x)']),
    ('optimize(3*exp(2*x) - 3, optims_c99)',
     ['ok', '3*expm1(2*x)']),
    ('optimize(log(x + 3)/log(2) + log(x**2 + 1), optims_c99)',
     ['ok', 'log1p(x**2) + log2(x + 3)']),
    ('rw.logsumexp_2terms_opt(log(exp(x) + exp(y)))',
     ['ok', 'log1p(exp(Min(x, y))) + Max(x, y)']),
    ('type(rw.logaddexp_opt).__name__',
     ['ok', 'ReplaceOptim']),
    ('type(rw.cosm1_opt).__name__',
     ['ok', 'FuncMinusOneOptim']),
    ('rw.cosm1_opt.func_m_1.__name__',
     ['ok', 'cosm1']),
    ('rw.powm1_opt.func.__name__',
     ['ok', 'Pow']),
    ('[type(o).__name__ for o in optims_scipy]',
     ['ok', "['FuncMinusOneOptim', 'FuncMinusOneOptim']"]),
    ('len(optims_numpy)',
     ['ok', '9']),
    ('rw.create_expand_pow_optimization(3)(x**5 + x**3)',
     ['ok', 'x**5 + x*x*x']),
    ('rw.FuncMinusOneOptim(cos, cosm1, opportunistic=False)(3*cos(x) - 5)',
     ['ok', '3*cos(x) - 5']),
]


@pytest.fixture(scope='module')
def nodes_ns():
    namespace = {}
    exec(NODES_SETUP, namespace)
    return namespace


@pytest.mark.parametrize('code,expected', NODES_REFERENCE)
def test_nodes_reference_outputs(nodes_ns, code, expected):
    assert _run(nodes_ns, code) == expected


def test_nodes_lambdify_numpy_and_scipy():
    import math
    np = pytest.importorskip('numpy')
    pytest.importorskip('scipy')
    from symox import lambdify, symbols
    from symox.codegen.numpy_nodes import logaddexp, logaddexp2, minimum, maximum
    from symox.codegen.scipy_nodes import cosm1, powm1
    x, y = symbols('x y')
    for modules in ('numpy', 'scipy'):
        assert lambdify((x, y), logaddexp(x, y), modules)(1.0, 2.0) == pytest.approx(math.log(math.exp(1) + math.exp(2)))
        assert lambdify((x, y), logaddexp2(x, y), modules)(1.0, 2.0) == pytest.approx(math.log2(6))
    assert lambdify(x, cosm1(x), 'scipy')(1e-3) == pytest.approx(-4.999999583333347e-07, rel=1e-12)
    assert lambdify((x, y), powm1(x, y), 'scipy')(1.5, 1e-10) == pytest.approx(4.054651081899185e-11, rel=1e-12)
    a, b = np.array([1., 5., 3.]), np.array([4., 2., 6.])
    assert lambdify((x, y), minimum(x, y), 'numpy')(a, b).tolist() == [1., 2., 3.]
    assert lambdify((x, y), maximum(x, y), 'numpy')(a, b).tolist() == [4., 5., 6.]


def test_numpy_max_of_three_arrays_keeps_its_inputs():
    # numpy.maximum(x, y, z) would take z as the output array.
    np = pytest.importorskip('numpy')
    from symox import Max, Min, lambdify, symbols
    from symox.printing.numpy import NumPyPrinter
    x, y, z = symbols('x y z')
    assert NumPyPrinter().doprint(Max(x, y, z)) == 'functools.reduce(numpy.maximum, [x, y, z])'
    a, b, c = np.array([1., 5., 3.]), np.array([4., 2., 6.]), np.array([0., 9., -1.])
    for F, expected in ((Max, [4., 9., 6.]), (Min, [0., 2., -1.])):
        args = (a.copy(), b.copy(), c.copy())
        assert lambdify((x, y, z), F(x, y, z), 'numpy')(*args).tolist() == expected
        assert [v.tolist() for v in args] == [a.tolist(), b.tolist(), c.tolist()]


def _codegen_lambdify_args(n, arrays):
    np = pytest.importorskip('numpy')
    if arrays:
        return [np.array([0.5, 1.5, 2.5]) * (i + 1) - i for i in range(n)]
    return [0.5 + 0.25 * i for i in range(n)]


@pytest.mark.parametrize('n', [2, 3, 5])
@pytest.mark.parametrize('modules', [None, 'numpy', 'scipy'])
@pytest.mark.parametrize('arrays', [False, True])
def test_codegen_minimum_maximum_lambdify_keeps_inputs(n, modules, arrays):
    # SymPy 1.14: NumPyPrinter._print_minimum folds numpy.minimum with
    # functools.reduce; numpy.minimum(x, y, z) would write into z (out).
    np = pytest.importorskip('numpy')
    if modules != 'numpy':
        pytest.importorskip('scipy')
    import inspect
    from symox import Max, Min, lambdify, symbols
    from symox.codegen.numpy_nodes import maximum, minimum
    syms = symbols('x0:%d' % n)
    for F, ufunc in ((minimum, np.minimum), (maximum, np.maximum), (Min, np.minimum), (Max, np.maximum)):
        values = _codegen_lambdify_args(n, arrays)
        kept = [np.copy(v) for v in values]
        f = lambdify(syms, F(*syms), modules)
        result = f(*values)
        assert np.array_equal(result, ufunc.reduce(np.array(kept)))
        assert all(np.array_equal(v, k) for v, k in zip(values, kept)), 'inputs modified'
        if F in (minimum, maximum):
            assert inspect.getsource(f).splitlines()[-1].strip() == \
                'return reduce(%s, [%s])' % (F.__name__, ', '.join(map(str, syms)))
            if not arrays:
                assert type(result) is np.float64


@pytest.mark.parametrize('modules', ['math', 'mpmath'])
def test_codegen_minimum_maximum_lambdify_python_printers(modules):
    # The math and mpmath printers print minimum/maximum under their own
    # name: the generated function raises NameError when called.
    from symox import lambdify, symbols
    from symox.codegen.numpy_nodes import maximum, minimum
    x, y, z = symbols('x y z')
    for F in (minimum, maximum):
        for args in ((x, y), (x, y, z)):
            f = lambdify(args, F(*args), modules)
            with pytest.raises(NameError, match="name '%s' is not defined" % F.__name__):
                f(*[0.5] * len(args))


def test_codegen_amin_amax_lambdify():
    np = pytest.importorskip('numpy')
    pytest.importorskip('scipy')
    import inspect
    from symox import Tuple, lambdify, sin, symbols
    from symox.codegen.numpy_nodes import amax, amin
    from symox.printing.codeprinter import PrintMethodNotImplementedError
    x = symbols('x')
    a = np.arange(6.).reshape(2, 3)
    kept = a.copy()
    for modules in (None, 'numpy', 'scipy'):
        for options in ({}, {'cse': True}, {'dummify': True}):
            assert lambdify(x, amin(x), modules, **options)(a) == 0.0
            assert lambdify(x, amax(x), modules, **options)(a) == 5.0
            assert lambdify(x, amin(x, axis=0), modules, **options)(a).tolist() == [0.0, 1.0, 2.0]
            assert lambdify(x, amax(x**2, axis=1), modules, **options)(a).tolist() == [4.0, 25.0]
            assert lambdify(x, amin([x, 2*x]), modules, **options)(a) == 0.0
            assert lambdify(x, amax(x, axis=Tuple(0, 1)), modules, **options)(a) == 5.0
    assert np.array_equal(a, kept)
    assert inspect.getsource(lambdify(x, amin(x), 'numpy')).splitlines()[-1].strip() == 'return amin(x, axis=None)'
    assert inspect.getsource(lambdify(x, amax(sin(x), axis=0), 'scipy')).splitlines()[-1].strip() == \
        'return amax(sin(x), axis=0)'
    for modules in ('math', 'mpmath'):
        for options in ({}, {'cse': True}):
            with pytest.raises(PrintMethodNotImplementedError, match='Unsupported by'):
                lambdify(x, amin(x), modules, **options)
            with pytest.raises(PrintMethodNotImplementedError, match='Unsupported by'):
                lambdify(x, amax(x, axis=0), modules, **options)


def test_codegen_nodes_lambdify_modules():
    # SymPy 1.14 lambdify of the numpy_nodes, scipy_nodes and cfunctions
    # classes: the printers' module functions compute, the names they do
    # not know raise NameError when the generated function is called.
    import math
    np = pytest.importorskip('numpy')
    pytest.importorskip('scipy')
    mpmath = pytest.importorskip('mpmath')
    from scipy import special
    from symox import lambdify, symbols
    from symox.codegen import cfunctions as cf
    from symox.codegen.numpy_nodes import logaddexp, logaddexp2
    from symox.codegen.scipy_nodes import cosm1, powm1
    x, y, z = symbols('x y z')
    nodes = {
        'logaddexp': ((x, y), logaddexp(x, y), lambda a, b: np.logaddexp(a, b)),
        'logaddexp2': ((x, y), logaddexp2(x, y), lambda a, b: np.logaddexp2(a, b)),
        'cosm1': ((x,), cosm1(x), lambda a: special.cosm1(a)),
        'powm1': ((x, y), powm1(x, y), lambda a, b: special.powm1(a, b)),
        'expm1': ((x,), cf.expm1(x), lambda a: np.expm1(a)),
        'log1p': ((x,), cf.log1p(x), lambda a: np.log1p(a)),
        'log2': ((x,), cf.log2(x), lambda a: np.log2(a)),
        'exp2': ((x,), cf.exp2(x), lambda a: np.exp2(a)),
        'log10': ((x,), cf.log10(x), lambda a: np.log10(a)),
        'fma': ((x, y, z), cf.fma(x, y, z), None),
        'hypot': ((x, y), cf.hypot(x, y), lambda a, b: np.hypot(a, b)),
        'Sqrt': ((x,), cf.Sqrt(x), lambda a: np.sqrt(a)),
        'Cbrt': ((x,), cf.Cbrt(x), None),
    }
    # Names each printer leaves undefined (NameError at call time).
    undefined = {
        'numpy': {'cosm1', 'powm1', 'fma', 'Cbrt'},
        'scipy': {'fma', 'Cbrt'},
        'math': {'logaddexp', 'logaddexp2', 'cosm1', 'powm1', 'exp2', 'fma', 'Cbrt'},
        'mpmath': {'logaddexp', 'logaddexp2', 'cosm1', 'exp2', 'fma', 'Cbrt'},
    }
    for modules, names in undefined.items():
        for name, (args, expr, reference) in nodes.items():
            scalars = [0.5 + 0.25 * i for i in range(len(args))]
            for options in ({}, {'cse': True}):
                f = lambdify(args, expr, modules, **options)
                if name in names:
                    with pytest.raises(NameError, match="name '%s' is not defined" % name):
                        f(*scalars)
                    continue
                value = f(*scalars)
                assert float(value) == pytest.approx(float(reference(*scalars)), rel=1e-14)
                if modules == 'mpmath' and not options:
                    assert isinstance(value, mpmath.mpf)
                if modules in ('numpy', 'scipy') and not options:
                    assert type(value) is np.float64
                    arrays = [np.array([0.5, 1.5, 2.5]) + i for i in range(len(args))]
                    kept = [a.copy() for a in arrays]
                    assert np.allclose(f(*arrays), reference(*kept), rtol=1e-14)
                    assert all(np.array_equal(a, k) for a, k in zip(arrays, kept))
    assert math.isclose(lambdify((x, y), powm1(x, y), 'mpmath')(0.5, 0.75), -0.405396442498, rel_tol=1e-10)


def test_derivative_printing_matches_codeprinter():
    from symox import Derivative, Function, Symbol, sin, zeta
    from symox.printing.codeprinter import PrintMethodNotImplementedError
    from symox.printing.pycode import MpmathPrinter, PythonCodePrinter
    x = Symbol('x')
    assert MpmathPrinter().doprint(Derivative(zeta(x), (x, 3))) == 'mpmath.zeta(x, derivative=3)'
    with pytest.raises(PrintMethodNotImplementedError, match='Printer has no method: _print_Derivative_sin'):
        PythonCodePrinter().doprint(Derivative(sin(x), x))
    assert PythonCodePrinter({'strict': False}).doprint(Derivative(sin(x), x)) == \
        '  # Not supported in Python:\n  # Derivative\nDerivative(sin(x), x)'
    with pytest.raises(ValueError, match='only supports functions with symbols as arguments'):
        PythonCodePrinter().doprint(Derivative(sin(2*x), x))
