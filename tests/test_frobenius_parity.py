"""dsolve(..., hint='2nd_power_series_regular') against SymPy 1.14.

The reference strings below were produced by SymPy 1.14.0.  Where SymPy
returns a Frobenius result, Symox must print the same series, order and form.
Where SymPy stops early, Symox only adds a documented extension (the
logarithmic second solution), and the reference part is kept unchanged.  Every
basis element Symox returns is certified by substituting it into the ODE with
SymPy: the residual vanishes through the computed order.
"""
import pytest
import sympy as sp
import symox as sx


def namespace(module):
    x = module.Symbol('x')
    f = module.Function('f')
    F = f(x)
    d1, d2 = F.diff(x), F.diff(x, 2)

    def hyp(a, b, c):
        return x*(1 - x)*d2 + (c - (a + b + 1)*x)*d1 - a*b*F

    return vars(module) | {
        'x': x, 'f': f, 'F': F, 'd1': d1, 'd2': d2, 'y': x - 1,
        'h': module.Rational(1, 2), 'hyp': hyp,
    }


def solve(source, options, **extra):
    names = namespace(sx)
    return sx.dsolve(eval(source, names), names['F'],
                     hint='2nd_power_series_regular', **options, **extra)


EXACT = [
    ('x**2*d2 + x*d1 + (x**2 - h**2)*F', {}, 'Eq(f(x), C2*sqrt(x)*(x**4/120 - x**2/6 + 1) + C1*(x**4/24 - x**2/2 + 1)/sqrt(x) + O(x**6))'),
    ('x**2*d2 + x*d1 + (x**2 - Rational(1, 9))*F', {}, 'Eq(f(x), C2*x**(1/3)*(9*x**4/896 - 3*x**2/16 + 1) + C1*(9*x**4/320 - 3*x**2/8 + 1)/x**(1/3) + O(x**6))'),
    ('x**2*d2 + x*d1 + (x**2 - Rational(9, 4))*F', {}, 'Eq(f(x), C2*x**(3/2)*(1 - x**2/10) + C1*(x**6/144 - x**4/8 + x**2/2 + 1)/x**(3/2) + O(x**6))'),
    ('x*d2 + (Rational(3, 2) - x)*d1 + 3*F', {}, 'Eq(f(x), C2*(-8*x**3/105 + 4*x**2/5 - 2*x + 1) + C1*(x**5/1080 + x**4/24 - 7*x**3/6 + 35*x**2/6 - 7*x + 1)/sqrt(x) + O(x**6))'),
    ('x**2*d2 + x*d1 - F', {}, 'Eq(f(x), C2*x + C1/x + O(x**6))'),
    ('x**2*d2 - 2*x*d1 + 2*F', {}, 'Eq(f(x), C2*x**2 + C1*x + O(x**6))'),
    ('x**2*d2 - 11*x*d1 + 36*F', {}, 'Eq(f(x), O(1))'),
    ('x*d2 + 2*d1 + x*F', {}, 'Eq(f(x), C2*(x**4/120 - x**2/6 + 1) + C1*(-x**6/720 + x**4/24 - x**2/2 + 1)/x + O(x**6))'),
    ('x*d2 + 2*d1 - x*F', {}, 'Eq(f(x), C2*(x**4/120 + x**2/6 + 1) + C1*(x**6/720 + x**4/24 + x**2/2 + 1)/x + O(x**6))'),
    ('x**2*d2 + x*d1 + (x**2 - h**2)*F', {'n': 4}, 'Eq(f(x), C2*sqrt(x)*(1 - x**2/6) + C1*(1 - x**2/2)/sqrt(x) + O(x**4))'),
    ('x**2*d2 + x*d1 + (x**2 - h**2)*F', {'n': 10}, 'Eq(f(x), C2*sqrt(x)*(x**8/362880 - x**6/5040 + x**4/120 - x**2/6 + 1) + C1*(x**8/40320 - x**6/720 + x**4/24 - x**2/2 + 1)/sqrt(x) + O(x**10))'),
    ('x*d2 + 2*d1 + x*F', {'n': 8}, 'Eq(f(x), C2*(-x**6/5040 + x**4/120 - x**2/6 + 1) + C1*(x**8/40320 - x**6/720 + x**4/24 - x**2/2 + 1)/x + O(x**8))'),
    ('x*d2 + 2*d1 + x*F', {'n': 2}, 'Eq(f(x), C2*(1 - x**2/2)/x + C1 + O(x**2))'),
    ('y**2*d2 + y*d1 + (y**2 - h**2)*F', {'x0': 1}, 'Eq(f(x), C2*sqrt(x - 1)*((x - 1)**4/120 - (x - 1)**2/6 + 1) + C1*((x - 1)**4/24 - (x - 1)**2/2 + 1)/sqrt(x - 1) + O(x**6))'),
    ('y*d2 + 2*d1 + y*F', {'x0': 1}, 'Eq(f(x), C2*((x - 1)**4/120 - (x - 1)**2/6 + 1) + C1*(-(x - 1)**6/720 + (x - 1)**4/24 - (x - 1)**2/2 + 1)/(x - 1) + O(x**6))'),
    ('(x + 1)**2*d2 + (x + 1)*d1 - Rational(1, 4)*F', {'x0': -1}, 'Eq(f(x), C2*sqrt(x + 1) + C1/sqrt(x + 1) + O(x**6))'),
    ('y*d2 + 2*d1 + y*F', {'x0': 1, 'n': 4}, 'Eq(f(x), C2*(1 - (x - 1)**2/6) + C1*((x - 1)**4/24 - (x - 1)**2/2 + 1)/(x - 1) + O(x**4))'),
]

ERRORS = [
    ('hyp(1, 1, 2)', {}, IndexError, 'list index out of range'),
    ('x**2*d2 + x*d1 + F', {}, ValueError, 'Expected Expr or iterable but got None'),
    ('x**2*d2 + x*d1 + (x**2 + 1)*F', {}, ValueError, 'Expected Expr or iterable but got None'),
    ("x**2*d2 + x*d1 + (x**2 - Symbol('nu')**2)*F", {}, ValueError, 'Expected Expr or iterable but got None'),
    ('x**2*d2 + x*d1 + x**2*F - x', {}, ValueError, 'ODE x**2*f(x) + x**2*Derivative(f(x), (x, 2)) + x*Derivative(f(x), x) - x does not match hint 2nd_power_series_regular'),
]

LOGARITHMIC = [
    ('x**2*d2 + x*d1 + x**2*F', {}, 'Eq(f(x), C1*(x**4/64 - x**2/4 + 1) + O(x**6))'),
    ('x**2*d2 + x*d1 + (x**2 - 1)*F', {}, 'Eq(f(x), C1*x*(x**4/192 - x**2/8 + 1) + O(x**6))'),
    ('x**2*d2 + x*d1 + (x**2 - 4)*F', {}, 'Eq(f(x), C1*x**2*(1 - x**2/12) + O(x**6))'),
    ('x*d2 + (1 - x)*d1 + 2*F', {}, 'Eq(f(x), C1*(x**2/2 - 2*x + 1) + O(x**6))'),
    ('x*d2 + (3 - x)*d1 + 2*F', {}, 'Eq(f(x), C1*(x**2/12 - 2*x/3 + 1) + O(x**6))'),
    ('x**2*d2 + 3*x*d1 + F', {}, 'Eq(f(x), C1/x + O(x**6))'),
    ('x*d2 + F', {}, 'Eq(f(x), C1*x*(x**4/2880 - x**3/144 + x**2/12 - x/2 + 1) + O(x**6))'),
    ('x**2*d2 - x*d1 + (1 - x)*F', {}, 'Eq(f(x), C1*x*(x**4/576 + x**3/36 + x**2/4 + x + 1) + O(x**6))'),
    ('x**2*d2 + x*d1 + (x**2 - 1)*F', {'n': 3}, 'Eq(f(x), C1*x + O(x**3))'),
    ('x**2*d2 + x*d1 + x**2*F', {'n': 10}, 'Eq(f(x), C1*(x**8/147456 - x**6/2304 + x**4/64 - x**2/4 + 1) + O(x**10))'),
]

CORRECTED = [
    # The reference cuts 2F1(1/2, 1; 3/2; x) = 1 + x/3 + ... down to its
    # constant term, which the residual check cannot reject, and numbers the
    # constants differently because that series sorts differently.
    ('hyp(h, 1, Rational(3, 2))', {}, 'Eq(f(x), C2/sqrt(x) + C1 + O(x**6))', 'Eq(f(x), C2*(x**5/11 + x**4/9 + x**3/7 + x**2/5 + x/3 + 1) + C1/sqrt(x) + O(x**6))', True),
    ('hyp(1, 2, h)', {}, 'Eq(f(x), C2*(-128*x**5/14175 + 32*x**4/315 - 32*x**3/45 + 8*x**2/3 - 4*x + 1) + C1*sqrt(x)*(32*x**4/2835 - 32*x**3/315 + 8*x**2/15 - 4*x/3 + 1) + O(x**6))', 'Eq(f(x), C2*(512*x**5/21 + 128*x**4/7 + 64*x**3/5 + 8*x**2 + 4*x + 1) + C1*sqrt(x)*(1155*x**4/128 + 105*x**3/16 + 35*x**2/8 + 5*x/2 + 1) + O(x**6))', False),
    ('x*(1 - x)*d2 + (Rational(1, 3) - 2*x)*d1 - Rational(1, 4)*F', {'x0': 1}, 'Eq(f(x), C2*(-1671168*x - 57*(x - 1)**5 + 1020*(x - 1)**4 - 15232*(x - 1)**3 + 182784*(x - 1)**2 + 12812288)/56853135360 + C1*(-27869184*x - 1813*(x - 1)**5 + 28860*(x - 1)**4 - 374400*(x - 1)**3 + 3773952*(x - 1)**2 + 362299392)/(56853135360*(x - 1)**(2/3)) + O(x**6))', 'Eq(f(x), C2*(-18382848*x - 2066715*(x - 1)**5 + 2891700*(x - 1)**4 - 4406400*(x - 1)**3 + 7755264*(x - 1)**2 + 140935168)/8129998356480 + C1*(-362299392*x - 18498601*(x - 1)**5 + 27275820*(x - 1)**4 - 45302400*(x - 1)**3 + 94348800*(x - 1)**2 + 4709892096)/(8129998356480*(x - 1)**(2/3)) + O(x**6))', False),
]

RATIONAL_P = 'Eq(f(x), C1*sqrt(x)*(983*x**4/49152 - 359*x**3/9216 + 13*x**2/192 - 3*x/8 + 1) + O(x**6))'


def _exponents(expr, t):
    """Powers of t with a nonzero coefficient in the expanded expr."""
    groups = {}
    expr = sp.expand(expr)
    for term in sp.Add.make_args(expr) if expr != 0 else ():
        coeff, dependent = term.as_independent(t, as_Add=False)
        dependent = sp.powsimp(dependent)
        base, power = dependent.as_base_exp() if dependent != 1 else (t, sp.S.Zero)
        assert base == t, dependent
        groups[power] = groups.get(power, 0) + coeff
    return [p for p, c in groups.items() if sp.simplify(c) != 0]


def certify(source, options, solution):
    """Substitute each basis element into the ODE with SymPy.

    log(x - x0) is held as a symbol and t = x - x0.  With the ODE multiplied
    by t**(2 - k), k the order of its leading coefficient at x0, the residual
    of a truncated Frobenius series must start strictly above the highest
    power of t that the series contains.  Returns the number of nonzero basis
    elements.
    """
    names = namespace(sp)
    x, F, d2 = names['x'], names['F'], names['d2']
    x0 = options.get('x0', 0)
    t, L = sp.symbols('t L')
    ode = sp.expand(sp.numer(sp.together(eval(source, names))))
    k = min(m[0] for m in sp.Poly(ode.coeff(d2).subs(x, t + x0), t).monoms())
    body = sp.sympify(str(solution.rhs), locals={'f': names['f']}).removeO()
    count = 0
    for constant in sp.symbols('C1 C2'):
        y = body.diff(constant)
        if y == 0:
            continue
        count += 1
        yt = y.subs(sp.log(x - x0), L).subs(x, t + x0)
        # An irrational leading power t**r is factored out of the series.
        r = sum((b.exp for b in sp.Mul.make_args(yt) if b.is_Pow and b.base == t
                 and not b.exp.is_Rational), sp.S.Zero)
        if r:
            yt = t**r*sp.Poly(yt/t**r, t, L, extension=True).as_expr()
            y = yt.subs(t, x - x0).subs(L, sp.log(x - x0))
        residual = ode.subs(F, y).doit().subs(sp.log(x - x0), L).subs(x, t + x0)
        residual = sp.Poly(sp.expand(residual/t**r), L)
        series = sp.Poly(sp.expand(yt/t**r), L)
        lowest = min([min(_exponents(c, t) or [sp.oo]) for c in residual.all_coeffs()])
        top = min(max(_exponents(c, t)) for c in series.all_coeffs() if c != 0)
        assert lowest + 2 - k > top, (constant, lowest, top)
    return count


@pytest.mark.parametrize('source, options, reference', EXACT)
def test_matches_reference(source, options, reference):
    solution = solve(source, options)
    assert str(solution) == reference
    certify(source, options, solution)


@pytest.mark.parametrize('source, options, error, message', ERRORS)
def test_raises_like_reference(source, options, error, message):
    with pytest.raises(error) as info:
        solve(source, options)
    assert str(info.value) == message


@pytest.mark.parametrize('source, options, reference', LOGARITHMIC)
def test_logarithmic_second_solution_extends_reference(source, options, reference):
    # SymPy 1.14 returns only the series of the larger indicial root when the
    # roots coincide or the recurrence breaks at their integer difference.
    # Symox appends C2 times the logarithmic solution; without it the output
    # is the reference.
    solution = solve(source, options)
    C2 = sx.Symbol('C2')
    assert str(solution.subs(C2, 0)) == reference
    extension = solution.rhs.coeff(C2)
    assert extension.has(sx.log) and not extension.has(sx.Symbol('C1'))
    assert certify(source, options, solution) == 2


@pytest.mark.parametrize('source, options, reference, expected, truncated', CORRECTED)
def test_non_polynomial_coefficients_are_expanded_in_x(source, options, reference, expected,
                                                       truncated):
    # SymPy 1.14 calls series(ser, x=x0, n=n+1) on a non-polynomial p or q,
    # which leaves it unexpanded, and reads wrong coefficients from it.  With
    # only that expansion corrected, SymPy's own pipeline prints `expected`,
    # and the reference is not a solution beyond its first terms.
    solution = solve(source, options)
    assert str(solution) == expected
    assert certify(source, options, solution) == 2
    if not truncated:
        with pytest.raises(AssertionError):
            certify(source, options, sx.sympify(reference))


def test_non_polynomial_coefficients_with_resonant_roots():
    source = 'x**2*d2 + x/(1 - x)*d1 - Rational(1, 4)*F/(1 + x)'
    solution = solve(source, {})
    assert str(solution.subs(sx.Symbol('C2'), 0)) == RATIONAL_P
    assert solution.rhs.coeff(sx.Symbol('C2')).has(sx.log)
    assert certify(source, {}, solution) == 2


def test_irrational_indicial_roots():
    # The raw series is the reference's (simplify=False).  The default path
    # raises "Not sure of sign of 6 - x0" in SymPy 1.14 from Order.contains
    # during constantsimp; Symox decides that comparison and returns the
    # certified basis.
    source = 'x**2*d2 + x*d1 + (x - 2)*F'
    names = namespace(sp)
    reference = sp.dsolve(eval(source, names), names['F'],
                          hint='2nd_power_series_regular', simplify=False)
    assert str(solve(source, {}, simplify=False)) == str(reference)
    with pytest.raises(NotImplementedError):
        sp.dsolve(eval(source, names), names['F'], hint='2nd_power_series_regular')
    assert certify(source, {}, solve(source, {})) == 2


@pytest.mark.parametrize('source', [
    # cse neither counts nor enters Order terms
    'cse((x + 1)*y + (x + 1)*z + O(x**4))',
    'cse([C1*(x**2 + 1) + O(x**3), C2*(x**2 + 1)])',
    # Mul subs keeps the residual coefficient as one factor
    '(Rational(3, 7)*C1*(x + 5)).subs(3*C1, C1)',
    '((3*a*b + c/x)/7).subs(3*a, w)',
    '(Rational(10, 7)*C2*(x + 2)/(x - 1)**Rational(2, 3) + Rational(3, 7)*C1*(x + 5)'
    ' + O(x**6)).subs(3*C1, C1)',
    # terms_gcd(expand=False) keeps shifted generators unexpanded
    'terms_gcd(C1*(x - 1)**2/2 + C1*(x - 1)/4, clear=False, deep=True, expand=False)',
    'terms_gcd(C1*x**2*(x**2/12 - 2*x/3 + 1)/(x - 1)**2 + C2*x/(x - 1),'
    ' clear=False, deep=True, expand=False)',
])
def test_constant_simplification_helpers(source):
    def run(module):
        names = vars(module) | dict(zip('xyzwabc', module.symbols('x y z w a b c')))
        names |= dict(zip(('C1', 'C2'), module.symbols('C1 C2')))
        return str(eval(source, names))
    assert run(sx) == run(sp)
