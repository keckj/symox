"""Error classes and messages matching the reference for invalid input.

Expected classes and texts were recorded from SymPy 1.14 for each example.
Module paths in a few messages name symox instead of sympy.
"""
import pytest

CASES = [
    ('BlockMatrix',
     'from symox import BlockMatrix, ones\ndat=[[ones(3,2), ones(3,3)*2],[ones(2,3)*3, ones(2,2)*4]]\nR=BlockMatrix(dat)',
     'ValueError', 'exact', '\nAlthough this matrix is comprised of blocks, the blocks do not fill\nthe matrix in a size-symmetric fashion. To create a full matrix from\nthese arguments, pass them directly to Matrix.'),
    ('ConditionSet',
     "from symox import ConditionSet, S, Symbol\nx=Symbol('x')\nR=ConditionSet(x + 1, x + 1 < 1, S.Integers)",
     'ValueError', 'exact', '`x + 1` is not symbol-like'),
    ('Domain.exquo',
     'from symox import ZZ\nR=ZZ.exquo(ZZ(5), ZZ(2))',
     'ExactQuotientFailed', 'exact', '2 does not divide 5 in ZZ'),
    ('Ellipse.reflect',
     'from symox import Ellipse, Point, Line\nR=Ellipse(Point(3, 4), 1, 3).reflect(Line(Point(0, -4), Point(5, 0)))',
     'NotImplementedError', 'exact', '\nGeneral Ellipse is not supported but the equation of the reflected\nEllipse is given by the zeros of: f(x, y) = (9*x/41 + 40*y/41 +\n37/41)**2 + (40*x/123 - 3*y/41 - 364/123)**2 - 1'),
    ('Relational.bool',
     "from symox import Symbol\nx=Symbol('x')\nR=(x < 1 != x < 1)",
     'TypeError', 'prefix', 'cannot determine truth value of Relational'),
    ('ITE',
     'from symox import ITE\nR=ITE(True, [], ())',
     'TypeError', 'exact', 'expecting bool or Boolean, not `[]`.'),
    ('ImmutableDenseMatrix',
     'from symox import ImmutableDenseMatrix\nm=ImmutableDenseMatrix([[1,2],[3,4]])\nm[0,0]=42',
     'TypeError', 'prefix', 'Cannot set values of'),
    ('ImmutableSparseMatrix',
     'from symox import ImmutableSparseMatrix\nm=ImmutableSparseMatrix([[1,2],[3,4]])\nm[0,0]=42',
     'TypeError', 'exact', "Cannot set values of <class 'symox._sparse.ImmutableSparseMatrix'>"),
    ('Poly.exquo',
     "from symox import Poly, Symbol\nx=Symbol('x')\nR=Poly(x**2 + 1, x).exquo(Poly(2*x - 4, x))",
     'ExactQuotientFailed', 'exact', '2*x - 4 does not divide x**2 + 1'),
    ('Poly.pexquo',
     "from symox import Poly, Symbol\nx=Symbol('x')\nR=Poly(x**2 + 1, x).pexquo(Poly(2*x - 4, x))",
     'ExactQuotientFailed', 'exact', '2*x - 4 does not divide x**2 + 1'),
    ('all_roots',
     "from symox import all_roots, sqrt, expand, Symbol\nx=Symbol('x')\np = expand((x - sqrt(2))*(x - sqrt(3)))\nR=all_roots(p)",
     'NotImplementedError', 'exact', 'sorted roots not supported over EX'),
    ('real_roots',
     "from symox import real_roots, sqrt, expand, Symbol\nx=Symbol('x')\np = expand((x - sqrt(2))*(x - sqrt(3)))\nR=real_roots(p)",
     'NotImplementedError', 'exact', 'sorted roots not supported over EX'),
    ('banded',
     'from symox import banded, ones\nR=banded({0: (2,)*5, 1: (ones(2),)*3})',
     'ValueError', 'exact', 'There is a collision at (1, 1) for 1 and 2.'),
    ('Permutation.call',
     'from symox.combinatorics import Permutation\nR=Permutation(1, 2)(3)',
     'IndexError', 'exact', 'list index out of range'),
    ('Permutation.resize',
     'from symox.combinatorics import Permutation\np=Permutation(0,1,2)\np=p.resize(4)\nR=p.resize(2)',
     'ValueError', 'exact', 'The permutation cannot be resized to 2 because the cycle (0, 1, 2) may break.'),
    ('CantSympify',
     'from symox import sympify\nfrom symox.core.sympify import CantSympify\nclass Something(dict, CantSympify):\n    pass\nR=sympify(Something())',
     'SympifyError', 'exact', 'SympifyError: {}'),
    ('padded_key',
     "from symox.crypto.crypto import padded_key\nR=padded_key('RSA', 'ARTIST')",
     'ValueError', 'exact', 'duplicate characters in symbols: T'),
    ('diff',
     "from symox import diff, sin, Symbol\nx,y=Symbol('x'),Symbol('y')\nR=diff(sin(x*y))",
     'ValueError', 'exact', '\nSince there is more than one variable in the expression, the\nvariable(s) of differentiation must be supplied to differentiate\nsin(x*y)'),
    ('exquo',
     "from symox import exquo, Symbol\nx=Symbol('x')\nR=exquo(x**2 + 1, 2*x - 4)",
     'ExactQuotientFailed', 'exact', '2*x - 4 does not divide x**2 + 1'),
    ('pexquo',
     "from symox import pexquo, Symbol\nx=Symbol('x')\nR=pexquo(x**2 + 1, 2*x - 4)",
     'ExactQuotientFailed', 'exact', '2*x - 4 does not divide x**2 + 1'),
    ('integrate',
     "from symox import integrate, Symbol\nx,y=Symbol('x'),Symbol('y')\nR=integrate(x*y)",
     'ValueError', 'exact', 'specify dummy variables for x*y'),
    ('linear_eq_to_matrix',
     "from symox import linear_eq_to_matrix, Symbol\nx,y=Symbol('x'),Symbol('y')\neqns = [(x**2 - 3*x)/(x - 3) - 3, y**2 - 3*y - y*(y - 4) + x - 4]\nR=linear_eq_to_matrix(eqns, [x, y])",
     'NonlinearError', 'exact', 'nonlinear term: 1/(x - 3)'),
    ('linsolve1',
     "from symox import linsolve, Symbol\nx=Symbol('x')\nR=linsolve([x*(1/x - 1)], x)",
     'NonlinearError', 'exact', 'nonlinear term: 1/x'),
    ('linsolve2',
     "from symox import linsolve, Symbol\nx,y=Symbol('x'),Symbol('y')\nR=linsolve([x*(y + 1)], x, y)",
     'NonlinearError', 'exact', '\nnonlinear cross-term: x*(y + 1)'),
    ('linsolve3',
     "from symox import linsolve, Symbol\nx=Symbol('x')\nR=linsolve([x**2 - 1], x)",
     'NonlinearError', 'exact', 'nonlinear term: x**2'),
    ('nsolve',
     "from symox import nsolve, cos, cosh, Symbol\nx=Symbol('x')\nf = cos(x)*cosh(x) - 1\nR=nsolve(f, 3.14*100)",
     'ValueError', 'prefix', 'Could not find root within given tolerance'),
    ('aesara_function',
     'from symox.abc import x\nfrom symox.printing.aesaracode import aesara_function\nR=aesara_function([x], [x**2 - 1], scalar=True)',
     'ImportError', 'exact', 'Aesara is required for aesara_function'),
    ('series',
     "from symox import series, tan, Symbol, oo\nx=Symbol('x')\nf=tan(x)\nR=series(f, x, 2, oo, '+')",
     'TypeError', 'exact', 'Cannot convert oo to int'),
    ('sift',
     'from symox import sift, exp, pi, Tuple\nunknown = exp(1) - pi\nargs = Tuple(1, pi, unknown)\nR=sift(args, lambda x: x.is_rational, binary=True)',
     'ValueError', 'exact', 'keyfunc gave non-binary output'),
    ('solve_poly_system',
     "from symox import solve_poly_system, Symbol\nx,y=Symbol('x'),Symbol('y')\nR=solve_poly_system([x**5 - x + y**3, y**2 - 1], x, y, strict=True)",
     'UnsolvableFactorError', 'exact', '\nStrict mode: some factors cannot be solved in radicals, so a complete\nlist of solutions cannot be returned. Call roots with strict=False to\nget solutions expressible in radicals (if there are any).'),
    ('lpmax',
     "from symox.solvers.simplex import lpmax\nfrom symox import Symbol\nx=Symbol('x')\nR=lpmax(x, [x <= -1, x >= 0])",
     'InfeasibleLPError', 'exact', '\nInconsistent/False constraint'),
    ('lpmin',
     "from symox.solvers.simplex import lpmin\nfrom symox import Symbol\nx=Symbol('x')\nR=lpmin(x, [x <= 3])",
     'UnboundedLPError', 'exact', '\nObjective function can assume arbitrarily large values!'),
    ('linear_coeffs',
     "from symox.solvers.solveset import linear_coeffs\nfrom symox import Symbol\nx,y=Symbol('x'),Symbol('y')\nR=linear_coeffs(x*(y + 1), x, y)",
     'NonlinearError', 'exact', '\nnonlinear cross-term: x*(y + 1)'),
    ('sympify',
     "from symox import sympify\nR=sympify('x***2')",
     'SympifyError', 'exact', "Sympify of expression 'could not parse 'x***2'' failed, because of exception being raised:\nSyntaxError: invalid syntax (<string>, line 1)"),
    ('raises',
     'from symox.testing.pytest import raises\nR=raises(ZeroDivisionError, lambda: 1/2)',
     'Failed', 'exact', 'DID NOT RAISE'),
    ('warns',
     'from symox.testing.pytest import warns\nwith warns(UserWarning):\n    pass',
     'Failed', 'exact', "Failed: DID NOT WARN. No warnings of type <class 'UserWarning'> was emitted. The list of emitted warnings is: []."),
    ('warns_deprecated_sympy',
     'from symox.testing.pytest import warns_deprecated_sympy\nwith warns_deprecated_sympy():\n    pass',
     'Failed', 'exact', "Failed: DID NOT WARN. No warnings of type <class 'symox.utilities.exceptions.SymPyDeprecationWarning'> was emitted. The list of emitted warnings is: []."),
]


@pytest.mark.parametrize('name,code,cls,kind,message', CASES, ids=[c[0] for c in CASES])
def test_error_class_and_message(name, code, cls, kind, message):
    ns = {}
    with pytest.raises(BaseException) as info:
        exec(code, ns)
    assert type(info.value).__name__ == cls
    text = str(info.value)
    if kind == 'exact':
        assert text == message
    else:
        assert text.startswith(message)


def test_valid_results_kept():
    from symox import ITE, Symbol, series, tan, oo, integrate, diff, sin
    from symox.combinatorics import Permutation
    x, y = Symbol('x'), Symbol('y')
    assert ITE(x > 0, 1, 0) == (x > 0)
    assert ITE(True, x, y) == x
    assert Permutation(1, 2)(2) == 1
    assert series(5, x, 0, oo) == 5
    assert integrate(5*x) == 5*x**2/2
    assert diff(sin(x)) == diff(sin(x), x)
