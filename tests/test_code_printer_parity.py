"""SymPy's C and C++ code printers.

The printers are ports of SymPy's codeprinter.py, c.py and cxx.py. Every
case prints on SymPy and Symox and compares the code, or the exception type
and message.
"""
import warnings

import pytest
import sympy
import symox

MODULES = {'sympy': sympy, 'symox': symox}

EXPRESSIONS = [
    'm.sin(x)**2+m.sqrt(y)/3+m.Rational(1,3)', 'm.Max(x,y)+m.pi', "m.Float('0.1',30)*x", 'x**-1',
    'x**m.Rational(3,2)', 'm.exp(x)+m.log(y)', 'm.Abs(x)', 'm.ceiling(x)+m.floor(y)',
    'm.Piecewise((x,x>0),(y,True))', 'm.Piecewise((x,x>0),(y,x<-1),(z,True))', 'm.Min(x,y,z)', 'm.atan2(x,y)',
    'm.gamma(x)', 'm.erf(x)', 'm.E**x', 'm.sinh(x)', 'm.asin(x)', 'x**2', 'x**3', 'x**y', 'm.Mod(x,y)',
    'm.And(x>0,y<1)', 'm.Or(x>0,m.Not(y<1))', 'm.Eq(x,y)', 'm.Ne(x,1)', '2*x*y/3', '-x/y', '1/(x*y)', 'm.I*x',
    'm.oo', '-m.oo', 'm.nan', 'm.zoo', 'm.Float(1.5)', "m.Float('1e-30')", 'm.Float(2)**x', 'm.sign(x)',
    "m.Function('f')(x)", "m.Function('sin')(x)", 'm.Heaviside(x)', 'm.log(x,2)', 'm.log(x,10)', 'm.exp(x)-1',
    'm.log(1+x)', 'm.Integer(10)**40', 'm.Rational(1,2)', 'm.S.Half*x', 'x**m.Rational(1,3)', 'm.cbrt(x)',
    'm.re(x)', "m.Symbol('if')", 'm.loggamma(x)', 'm.sec(x)', 'm.cot(x)', 'm.acosh(x)',
    "m.MatrixSymbol('A',2,2)[0,1]", "m.MatrixSymbol('A',2,3)[1,2]+x", 'm.Xor(x>0,y>0)', 'm.ITE(x>0,y>0,z>0)',
    'm.Lambda(x,x+1)', 'm.Tuple(x,y)', 'x>=y', 'm.sqrt(2)', 'm.exp(1)', 'm.Pow(x,-2,evaluate=False)', 'm.sinc(x)',
    'm.Float(3.25)*x**2', "m.Matrix([[x,y],[1,2]])"
]

C_SETTINGS = [
    {}, {'precision': 5}, {'standard': 'c89'}, {'standard': 'c11'}, {'assign_to': 'w'},
    {'user_functions': {'f': 'ff'}}, {'contract': False}, {'human': False}, {'full_prec': True},
    {'allow_unknown_functions': True}, {'math_macros': {}}, {'dereference': ['x']}, {'strict': False},
]
CXX_SETTINGS = [{}, {'standard': 'c++98'}, {'standard': 'c++17'}, {'precision': 3}]

CASES = ([('ccode', e, s) for e in EXPRESSIONS for s in C_SETTINGS]
         + [('cxxcode', e, s) for e in EXPRESSIONS for s in CXX_SETTINGS])


def outcome(name, function, expression, settings):
    m = MODULES[name]
    x, y, z = m.symbols('x y z')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            value = eval(expression, {'m': m, 'x': x, 'y': y, 'z': z})
            return ('value', repr(getattr(m, function)(value, **settings)))
        except Exception as error:  # noqa: BLE001 - the error is compared
            message = str(error).replace(name, 'M')
            # The matrix classes are not nominal: only their name differs.
            message = message.replace('M._core.Matrix', 'M.matrices.immutable.ImmutableDenseMatrix')
            return ('error', type(error).__name__, message)


@pytest.mark.parametrize('function, expression, settings', CASES,
                         ids=[f'{f}-{e}-{sorted(s)}' for f, e, s in CASES])
def test_code_matches_sympy(function, expression, settings):
    assert outcome('symox', function, expression, settings) == outcome('sympy', function, expression, settings)


def printer_subclass(m):
    """A downstream printer customizing floats, types and precedence, as
    Devito's does."""
    from importlib import import_module
    c = import_module(m.__name__ + '.printing.c')
    precedence = import_module(m.__name__ + '.printing.precedence')
    codegen = import_module(m.__name__ + '.codegen.ast')

    class Printer(c.C99CodePrinter):
        type_mappings = {**c.C99CodePrinter.type_mappings, codegen.real: 'float'}

        def _print_Float(self, expr):
            return super()._print_Float(expr) + 'F'

        def _print_Pow(self, expr):
            if expr.exp == -1:
                return '1.0F/' + self.parenthesize(expr.base, precedence.PRECEDENCE['Pow'])
            return super()._print_Pow(expr)

    return Printer


SUBCLASS_CASES = ['x + m.Float(1.5)', '1/(x+y)', 'm.sqrt(x)*m.Float(0.25)', 'm.Float(2.0)**x']


@pytest.mark.parametrize('expression', SUBCLASS_CASES)
def test_printer_subclass_matches_sympy(expression):
    def result(name):
        m = MODULES[name]
        x, y = m.symbols('x y')
        return printer_subclass(m)().doprint(eval(expression, {'m': m, 'x': x, 'y': y}))
    assert result('symox') == result('sympy')


def test_piecewise_assignment_and_comment_statements():
    from symox.codegen.ast import Assignment, CodeBlock, Comment
    from symox.printing.c import C99CodePrinter
    x, y, w = symox.symbols('x y w')
    code = C99CodePrinter().doprint(Assignment(w, symox.Piecewise((x, x > 0), (y, True))))
    assert code == 'if (x > 0) {\n   w = x;\n}\nelse {\n   w = y;\n}'
    assert C99CodePrinter().doprint(CodeBlock(Comment('step'), Assignment(w, x))) == '/* step */;\nw = x;'


def test_assignment_checks_matrix_structure():
    from symox.codegen.ast import Assignment
    x = symox.Symbol('x')
    with pytest.raises(ValueError, match='Dimensions of lhs and rhs do not align'):
        Assignment(symox.MatrixSymbol('A', 3, 1), symox.ones(2, 2))
    with pytest.raises(ValueError, match='Cannot assign a matrix to a scalar'):
        Assignment(x, symox.ones(2, 2))
    with pytest.raises(ValueError, match='Cannot assign a scalar to a matrix'):
        Assignment(symox.MatrixSymbol('A', 2, 2), x)


def test_codegen_tokens_follow_sympy_protocols():
    from symox.codegen.ast import Variable, While, none, real
    x = symox.Symbol('x')
    assert none == None  # noqa: E711 - SymPy's NoneToken equals None
    assert Variable(x, real).value == None  # noqa: E711
    loop = While(x < 1, [x])
    assert set(loop.kwargs()) == {'condition', 'body'}
    assert loop.kwargs(exclude=('body',), apply=str) == {'condition': 'x < 1'}


def test_matrix_element_parts():
    a = symox.MatrixSymbol('A', 2, 3)
    element = a[1, 2]
    assert (element.parent, element.i, element.j) == (a, 1, 2)
    assert not hasattr(symox.Symbol('x'), 'parent')
    with pytest.raises(AttributeError, match="'Symbol' object has no attribute 'i'"):
        symox.Symbol('x').i


def test_float_mpf_matches_sympy():
    for value in ['1.5', '0.1', '-2.5e-7', '0', '1e20', '3', '1e-400', '-7e100000']:
        for precision in (53, 100):
            assert symox.Float(value, precision=precision)._mpf_ == sympy.Float(value, precision=precision)._mpf_
    # Float('inf') is oo, which has no mpf value.
    assert [hasattr(m.Float(v), '_mpf_') for m in (symox, sympy) for v in ('inf', 'nan')] == [False] * 4
    with pytest.raises(AttributeError):
        symox.Integer(2)._mpf_


def test_strings_never_equal_expressions():
    x = symox.Symbol('x')
    assert (x == 'x') is False and (x != 'x') is True
    assert (symox.Integer(1) == '1') is False
    assert (sympy.Symbol('x') == 'x') is False


def test_ite_rewrites_as_piecewise():
    def rewritten(m):
        x, y, z = m.symbols('x y z')
        return str(m.ITE(x > 0, y > 0, z > 0).rewrite(m.Piecewise))
    assert rewritten(symox) == rewritten(sympy)


def test_float_strings_match_mpmath():
    import random
    from mpmath.libmp import from_man_exp, prec_to_dps, to_str
    from symox.printing import _mpf_str
    generator = random.Random(3)
    cases = [(value, dps, True, None, None) for value in [(0, 0, 0, 0), (0, 0, -456, -2), (1, 0, -789, -3),
                                                        (0, 0, -123, -1)] for dps in (0, 1, 15)]
    for _ in range(3000):
        man = (generator.getrandbits(generator.randint(1, 300)) | 1) * generator.choice([1, -1])
        exp = generator.choice([generator.randint(-80, 80), generator.randint(-5000, 5000),
                                generator.randint(-10**9, 10**9)])
        cases.append((from_man_exp(man, exp), generator.choice([0, 1, 2, 15, 17, 30, 100]),
                      generator.choice([True, False]), generator.choice([None, -3, 0]),
                      generator.choice([None, 0, 5, 40])))
    for value, dps, strip, low, high in cases:
        expected = to_str(value, dps, strip_zeros=strip, min_fixed=low, max_fixed=high)
        assert _mpf_str.to_str(value, dps, strip_zeros=strip, min_fixed=low, max_fixed=high) == expected
    assert all(_mpf_str.prec_to_dps(n) == prec_to_dps(n) for n in range(2000))


def test_float_printer_settings_match_sympy():
    def printed(m):
        value = m.Float('1.23456789012345678901234567890', 30)
        return [m.sstr(value), m.sstr(value, full_prec=False), m.sstr(value, min=-2, max=2),
                m.sstr(m.Float('1e-7') * m.Symbol('x')), m.sstr(m.Float('123456.5'), max=3),
                m.ccode(m.Float('0.1', 40)), m.ccode(m.Float('1e-300'), precision=3)]
    assert printed(symox) == printed(sympy)


def test_dot_product_code():
    a, b = symox.MatrixSymbol('a', 1, 3), symox.MatrixSymbol('b', 3, 1)
    # A Symox extension: the reference raises for any dot product.
    assert symox.ccode(symox.DotProduct(a, b)) == 'a[0]*b[0] + a[1]*b[1] + a[2]*b[2]'
    n = symox.Symbol('n')
    symbolic = symox.DotProduct(symox.MatrixSymbol('A', n, 1), symox.MatrixSymbol('B', n, 1))
    with pytest.raises(symox.printing.codeprinter.PrintMethodNotImplementedError):
        symox.ccode(symbolic)
    assert symox.ccode(symbolic, strict=False) == '/* Not supported in C: */\n/* DotProduct */\nDotProduct(A, B)'


def test_augmented_assignment_operators():
    def printed(m):
        from importlib import import_module
        ast = import_module(m.__name__ + '.codegen.ast')
        x, y = m.symbols('x y')
        return [ast.Assignment(x, y).op, ast.AddAugmentedAssignment(x, y).op, ast.aug_assign(x, '%', y).op,
                m.ccode(ast.MulAugmentedAssignment(x, y)), m.cxxcode(ast.DivAugmentedAssignment(x, y))]
    assert printed(symox) == printed(sympy)
