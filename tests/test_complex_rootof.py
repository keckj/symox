"""Indexed complex roots (``CRootOf``) compared with SymPy.

Every case is evaluated in both libraries from the same source text. The
printed value, or the exception class and message, must agree: ordering of
the roots, the arguments, the printed and ``srepr`` forms, the assumptions,
exact isolation intervals, the rational approximations and the
arbitrary-precision values.
"""
import pytest
import symox as sm

POLYS = ['x**5-x-1', '2*x**3-x+5', 'x**3/2-x+5', '(x**2-2)*(x**3-2)', '4*x**5-2*x+1',
         'x**10-2*x+3', 'x**7-3*x+1', 'x**4+1', 'x**6+x+1', 'x**5+x**2+1',
         '(x**3-x+1)*(x**3-x-1)', 'x**4-10*x**2+1']

CASES = [
    # construction, arguments and printing
    "CRootOf(x**5-x-1,0)", "CRootOf(x**5-x-1,1)", "CRootOf(x**5-x-1,4)",
    "CRootOf(x**5-x-1,0).args", "srepr(CRootOf(x**5-x-1,3))", "CRootOf(x**5-x-1,-1)",
    "CRootOf(x**5-x-1,x,-5)", "CRootOf(x**2-2,0)", "rootof(x**2-2,0)", "rootof(x**3-2,1)",
    "CRootOf(x**3-2,1,radicals=True)", "CRootOf(x**2+1,0)", "CRootOf(x-3,0)",
    "CRootOf(x**5-x-1,0).func", "CRootOf(x**5-x-1,0).poly", "CRootOf(x**5-x-1,0).index",
    "CRootOf(x**5-x-1,0).free_symbols", "CRootOf(x**5-x-1,0).is_number",
    "CRootOf(x**5-x-1,0) == CRootOf(x**5-x-1,0)", "CRootOf(x**5-x-1,0) == CRootOf(x**5-x-1,1)",
    "str(1 + 2*CRootOf(x**5-x-1,1))", "str(CRootOf(x**5-x-1,0)**2 + CRootOf(x**5-x-1,0))",
    "latex(CRootOf(x**5-x-1,0))", "CRootOf(x**5-x-1,0).as_poly()",
    # errors
    "CRootOf(x**5-x-1,5)", "CRootOf(x**5-x-1,-6)", "CRootOf(x*y-1,x,0)", "CRootOf(7,0)",
    "CRootOf(x**5+y*x-1,x,0)", "CRootOf(x**5-x-1,S(1)/2)", "Poly(x*CRootOf(x**5-x-1,0)+1,x)",
    # assumptions and conjugation
    "CRootOf(x**5-x-1,0).is_real", "CRootOf(x**5-x-1,0).is_positive", "CRootOf(x**5-x-1,1).is_real",
    "CRootOf(x**5-x-1,1).is_imaginary", "CRootOf(x**4+1,0).is_imaginary",
    "CRootOf(x**5-x-1,1).conjugate()", "CRootOf(x**5-x-1,3).conjugate()",
    "CRootOf(x**7-3*x+1,0).is_negative", "CRootOf(x**6+x+1,0).is_algebraic",
    # isolation intervals and approximations
    "CRootOf(x**10-2*x+3,9)._get_interval()", "CRootOf(x**3-x**2-x-1,0)._get_interval()",
    "CRootOf(x**5-x-1,0).eval_rational(n=10)", "CRootOf(x**5-x-1,1).eval_rational(n=10)",
    "CRootOf(x**5-x-1,1).eval_rational(dx=S(1)/100, dy=S(1)/100)",
    "CRootOf(x**5-x-1,0)._eval_evalf(20)", "CRootOf(x**5-x-1,1).evalf(30)",
    "CRootOf(x**5-x-1,0).evalf(50)", "CRootOf(x**5-x-1,3).n(5)", "N(CRootOf(x**7-3*x+1,3), 40)",
    "CRootOf(x**7-3*x+1,3).evalf(20)", "CRootOf(x**4+1,2).evalf()",
    "(1 + 2*CRootOf(x**5-x-1,1)).evalf(20)", "(solve(x**6+x+1,x)[0]+1).n()",
    # solving and calculus with root coefficients
    "solve(x**6+x+1,x)", "solve(x-CRootOf(x**5-x-1,0),x)", "solve(CRootOf(x**5-x-1,0)*x-1,x)",
    "integrate(CRootOf(x**5-x-1,0),x)", "integrate(CRootOf(x**5-x-1,0)*x,x)",
    "CRootOf(x**5-x-1,0).diff(x)", "(x*CRootOf(x**5-x-1,0)).diff(x)",
    "CRootOf(x**5-x-1,0).has(x)", "CRootOf(x**5-x-1,0).atoms(Symbol)",
    "Poly(x**5-x-1).intervals(all=True)", "Poly(x**4+1).intervals(all=True, eps=S(1)/10)",
]
CASES += ["Poly(%s).all_roots()" % p for p in POLYS]
CASES += ["Poly(%s).real_roots()" % p for p in POLYS]
CASES += ["Poly((x**2-2)**2*(x**3-x-1)).all_roots(multiple=False)",
          "Poly((x**2-2)**2*(x**3-x-1)).real_roots(multiple=False)",
          "real_roots((x-1)**2*(x**5-x-1))", "Poly(x**5-x-1).root(3)",
          "[r.evalf(20) for r in Poly(x**4-10*x**2+1).all_roots(radicals=False)]",
          "[r.evalf(25) for r in Poly(x**6+x+1).all_roots()]"]


def run(module, source):
    # Assumption queries and printing refine the cached isolation intervals
    # as a side effect in the reference; start each case from empty caches.
    module.CRootOf.clear_cache()
    namespace = dict(vars(module))
    x, y = module.symbols('x y')
    namespace.update(x=x, y=y)
    try:
        value = eval(source, namespace)
    except Exception as error:
        return '!%s: %s' % (type(error).__name__, error)
    text = str(value)
    return text.replace('sympy.polys.rootoftools', 'symox.polys.rootoftools')


@pytest.mark.parametrize('source', CASES)
def test_complex_rootof_matches_reference(source):
    import sympy
    assert run(sm, source) == run(sympy, source)


def test_root_float_rounds_once_at_each_precision():
    import sympy
    x = sm.Symbol('x')
    root = sm.CRootOf(x**7 - 3*x + 1, 4)
    oracle = sympy.CRootOf(sympy.Symbol('x')**7 - 3*sympy.Symbol('x') + 1, 4)
    for digits in (2, 5, 15, 16, 33, 60, 120):
        assert str(root.evalf(digits)) == str(oracle.evalf(digits))
