"""Reference ``solveset`` steps that the native solver does not take.

``_solve_modular``/``_invert_modular`` follow SymPy 1.14's
``sympy.solvers.solveset`` for ``A - Mod(B, C) = 0`` over the integers, and
``_solve_radical`` follows its ``unrad`` step for equations with radicals
(every candidate is checked with ``checksol``; an undecided candidate stays
in a ConditionSet).
"""


def _term_factors(f):
    from .. import Add, Mul
    for add_arg in Add.make_args(f):
        yield from Mul.make_args(add_arg)


def _is_modular(f, symbol):
    from .. import Mod
    if not f.has(Mod):
        return False
    modterms = list(f.atoms(Mod))
    return (len(modterms) == 1 and
            modterms[0].args[0].has(symbol) and
            modterms[0].args[1].is_integer and
            any(isinstance(term, Mod) for term in _term_factors(f)))


def _invert_modular(modterm, rhs, n, symbol):
    from .. import (S, ImageSet, Lambda, Mod, discrete_log, divisors,
                    igcd, invert, nthroot_mod, totient)
    a, m = modterm.args
    if rhs.is_integer is False:
        return symbol, S.EmptySet
    if rhs.is_real is False or any(term.is_real is False
                                   for term in _term_factors(a)):
        return modterm, rhs
    if abs(rhs) >= abs(m):
        return symbol, S.EmptySet
    if a == symbol:
        return symbol, ImageSet(Lambda(n, m*n + rhs), S.Integers)
    if a.is_Add:
        g, h = a.as_independent(symbol)
        if g is not S.Zero:
            x_indep_term = rhs - Mod(g, m)
            return _invert_modular(Mod(h, m), Mod(x_indep_term, m), n, symbol)
    if a.is_Mul:
        g, h = a.as_independent(symbol)
        if g is not S.One:
            x_indep_term = rhs*invert(g, m)
            return _invert_modular(Mod(h, m), Mod(x_indep_term, m), n, symbol)
    if a.is_Pow:
        base, expo = a.args
        if expo.has(symbol) and not base.has(symbol):
            if not m.is_Integer and rhs.is_Integer and a.base.is_Integer:
                return modterm, rhs
            mdiv = m.p // igcd(m.p, rhs.p)
            try:
                remainder = discrete_log(mdiv, rhs.p, a.base.p)
            except ValueError:
                return modterm, rhs
            period = totient(m)
            for p in divisors(period):
                if pow(a.base, p, m / igcd(m.p, a.base.p)) == 1:
                    period = p
                    break
            return expo, ImageSet(Lambda(n, period*n + remainder), S.Naturals0)
        elif base.has(symbol) and not expo.has(symbol):
            try:
                remainder_list = nthroot_mod(rhs, expo, m, all_roots=True)
                if remainder_list == []:
                    return symbol, S.EmptySet
            except (ValueError, NotImplementedError):
                return modterm, rhs
            g_n = S.EmptySet
            for rem in remainder_list:
                g_n += ImageSet(Lambda(n, m*n + rem), S.Integers)
            return base, g_n
    return modterm, rhs


def _normalize(f, symbol):
    """The coefficient normalization at the start of ``_solveset``."""
    from .. import S
    if f.is_Add:
        a, h = f.as_independent(symbol)
        m, h = h.as_independent(symbol, as_Add=False)
        if m not in {S.ComplexInfinity, S.Zero, S.Infinity, S.NegativeInfinity}:
            f = a/m + h
    return f


def solve_modular(f, symbol, domain, solver):
    """SymPy's ``_solve_modular`` (``solver`` is the recursive solveset)."""
    from .. import S, ConditionSet, Dummy, Eq, FiniteSet, ImageSet, Lambda, Mod
    unsolved_result = ConditionSet(symbol, Eq(f, 0), domain)
    modterm = list(f.atoms(Mod))[0]
    rhs = -S.One*(f.subs(modterm, S.Zero))
    if f.as_coefficients_dict()[modterm].is_negative:
        rhs *= -S.One
    if not domain.is_subset(S.Integers):
        return unsolved_result
    if rhs.has(symbol):
        return unsolved_result
    n = Dummy('n', integer=True)
    f_x, g_n = _invert_modular(modterm, rhs, n, symbol)
    if f_x == modterm and g_n == rhs:
        return unsolved_result
    if f_x == symbol:
        if domain is not S.Integers:
            return domain.intersect(g_n)
        return g_n
    if isinstance(g_n, ImageSet):
        lamda_expr = g_n.lamda.expr
        lamda_vars = g_n.lamda.variables
        base_sets = g_n.base_sets
        sol_set = solver(f_x - lamda_expr, symbol, S.Integers)
        if isinstance(sol_set, FiniteSet):
            tmp_sol = S.EmptySet
            for sol in sol_set:
                tmp_sol += ImageSet(Lambda(lamda_vars, sol), *base_sets)
            sol_set = tmp_sol
        else:
            sol_set = ImageSet(Lambda(lamda_vars, sol_set), *base_sets)
        return domain.intersect(sol_set)
    return unsolved_result


def modular_step(f, symbol, domain, solver):
    """Returns the reference result for a modular equation, or None."""
    from .. import Expr, Equality, sympify
    f = sympify(f)
    if isinstance(f, Equality):
        f = f.lhs - f.rhs
    if not isinstance(f, Expr) or f.is_Mul or not f.has(symbol):
        return None
    f = _normalize(f, symbol)
    if not _is_modular(f, symbol):
        return None
    return solve_modular(f, symbol, domain, solver)


def _has_radical(f, symbol):
    for p in f.atoms(__import__('symox').Pow):
        if p.base.has(symbol) and p.exp.is_Rational and not p.exp.is_Integer:
            return True
    return False


def _invert_real_radical(indep, dep, symbol, domain, solver):
    """``_invert_real`` of ``dep = -indep`` for ``dep = g**(num/den)`` with
    an even ``den`` and an odd ``num``: the right side must be nonnegative
    and ``g = (-indep)**(den/num)``; the new equation is solved again.
    The odd-denominator branch is not taken: the reference uses the real
    root there, which is not the principal power (``x**(1/3) + 2`` has no
    real solution). Returns None to keep the native result."""
    from .. import S
    if domain.is_subset(S.Reals) is not True or not dep.is_Pow:
        return None
    base, expo = dep.args
    if expo.has(symbol) or not expo.is_Rational or expo.is_Integer:
        return None
    num, den = expo.as_numer_denom()
    if not (den % 2 == 0 and num % 2 == 1):
        return None
    rhs = -indep
    if rhs.is_negative:
        return S.EmptySet
    if rhs.is_nonnegative is not True:
        return None
    return solver(base - rhs**(1/expo), symbol, domain)


def radical_step(f, symbol, domain, solver):
    """SymPy's ``_solve_radical`` for an equation that the native solver
    leaves as a ConditionSet. Returns None to keep the native result."""
    from .. import (I, S, Complement, ConditionSet, Dummy, Eq, Expr,
                    Equality, FiniteSet, Lambda, Union, imageset, sympify)
    from .._solver_checks import checksol
    from ._algebraic import denoms, unrad
    f = sympify(f)
    if isinstance(f, Equality):
        f = f.lhs - f.rhs
    if not isinstance(f, Expr) or not getattr(symbol, 'is_Symbol', False):
        return None
    f = _normalize(f, symbol)
    if not _has_radical(f, symbol):
        return None
    # ``_invert`` only moves the symbol-free terms of a sum, so the
    # reference reaches ``unrad`` with ``equation == f`` exactly when the
    # symbol-dependent part is itself a sum. A single radical or power is
    # inverted first (and may stay a ConditionSet there).
    indep, dep = f.as_independent(symbol, as_Add=True)
    if dep.is_Mul:
        # the reference inverters divide by a symbol-free coefficient and
        # solve again; a product without one reaches ``unrad`` unchanged
        g, h = dep.as_independent(symbol, as_Add=False)
        if g is not S.One:
            if g.is_zero is not False or g.is_finite is not True:
                return None
            return solver(h + indep/g, symbol, domain)
    elif not dep.is_Add:
        return _invert_real_radical(indep, dep, symbol, domain, solver)
    try:
        u = unrad(f, symbol)
    except (NotImplementedError, ValueError):
        return None
    if not u:
        return None
    eq, cov = u
    if not cov:
        result = solver(eq, symbol, domain) - \
            Union(*[solver(g, symbol, domain) for g in denoms(f, symbol)])
    else:
        y, yeq = cov
        if not solver(y - I, y, domain):
            yreal = Dummy('yreal', real=True)
            yeq = yeq.xreplace({y: yreal})
            eq = eq.xreplace({y: yreal})
            y = yreal
        g_y_s = solver(yeq, symbol, domain)
        f_y_sols = solver(eq, y, domain)
        result = Union(*[imageset(Lambda(y, g_y), f_y_sols) for g_y in g_y_s])

    def check_finiteset(solutions):
        f_set = []
        c_set = []
        for s in solutions:
            if checksol(f, symbol, s):
                f_set.append(s)
            else:
                c_set.append(s)
        return FiniteSet(*f_set) + ConditionSet(symbol, Eq(f, 0), FiniteSet(*c_set))

    def check_set(solutions):
        if solutions is S.EmptySet:
            return solutions
        elif isinstance(solutions, ConditionSet):
            return solutions
        elif isinstance(solutions, FiniteSet):
            return check_finiteset(solutions)
        elif isinstance(solutions, Complement):
            A, B = solutions.args
            return Complement(check_set(A), B)
        elif isinstance(solutions, Union):
            return Union(*[check_set(s) for s in solutions.args])
        return solutions

    return check_set(result)


def _negative_base_pows(f, symbol):
    """The powers ``b**g(x)`` of ``f`` with a constant base that is not
    known to be positive (or zero): SymPy's ``_invert_real`` inverts them
    only for a negative base, with ``integer_log``, and otherwise keeps a
    ConditionSet (``c**x = 4`` has the solution 2 for ``c = -2``, which the
    principal logarithm misses)."""
    from .. import Pow
    return [p for p in f.atoms(Pow) if p.exp.has(symbol) and not p.base.has(symbol)
            and p.base.is_positive is not True and p.base.is_zero is not True]


def _finite_with_finite_vars(f, domain):
    """SymPy's ``_is_finite_with_finite_vars``."""
    from .. import S, Dummy
    def assumptions(s):
        A = dict(s.assumptions0)
        A.setdefault('finite', A.get('finite', True))
        if domain.is_subset(S.Reals):
            A.setdefault('real', True)
        else:
            A.setdefault('complex', True)
        return A
    reps = {s: Dummy(**assumptions(s)) for s in f.free_symbols}
    return f.xreplace(reps).is_finite


def _rational_log(a, m):
    """The integer ``k`` with ``a**k == m`` for positive Rationals ``a != 1``
    and ``m``, or None."""
    import math
    from .. import Rational
    if m == 1:
        return 0
    la = math.log(a.p) - math.log(a.q)
    lm = math.log(m.p) - math.log(m.q)
    k0 = round(lm/la)
    for k in (k0, k0 - 1, k0 + 1):
        if Rational(a)**k == m:
            return k
    return None


def _negative_base_exponents(base, rhs):
    """The real exponents ``e`` with ``base**e == rhs`` for a negative
    constant ``base``: ``('ints', [k])`` with the integer solutions,
    ``('parity', j)`` for ``base == -1`` (``e`` in ``2*Integers + j``), or
    None when it is not decided.

    For a real ``e``, ``(-a)**e = a**e*exp(I*pi*e)`` is real only for an
    integer ``e``, and then it is ``(-1)**e*a**e``. SymPy's ``_invert_real``
    finds the nonnegative ``e`` with ``integer_log(rhs, base)`` and raises
    ValueError for other right-hand sides."""
    from .. import S, Abs, Float, cancel, expand_log, log
    if rhs.free_symbols or base.free_symbols or rhs.has(S.NaN, S.ComplexInfinity):
        return None
    if rhs.is_zero:
        return ('ints', [])
    if rhs.is_zero is not False or rhs.is_finite is not True:
        return None
    if rhs.is_extended_real is not True:
        return None
    if base.atoms(Float) or rhs.atoms(Float):
        return None
    if base == -1:
        if rhs == 1:
            return ('parity', 0)
        if rhs == -1:
            return ('parity', 1)
        if rhs.is_Rational:
            return ('ints', [])
        return None
    a = -base
    if a.is_Rational and rhs.is_Rational:
        k = _rational_log(a, Abs(rhs))
    else:
        k = cancel(expand_log(log(Abs(rhs)), force=True)/expand_log(log(a), force=True))
        if not k.is_Rational:
            return None
        if not k.is_Integer:
            return ('ints', [])
        k = int(k)
        if (a**k - Abs(rhs)).equals(0) is not True:
            return None
    if k is None:
        return ('ints', [])
    sign = rhs.is_positive
    if sign is None:
        return None
    return ('ints', [k] if (k % 2 == 0) == sign else [])


def negative_base_step(f, symbol, domain, solver):
    """``_invert_real`` of an equation whose symbol sits in the exponent of a
    power with a negative constant base (``(-2)**x = 4``) on a real domain.

    The symbol-free terms and factors are moved to the right side and
    ``exp``/``log`` are inverted (as ``_invert_real`` does) until the power
    is reached; its exponent is then solved for the exact integer powers
    found by ``_negative_base_exponents``. When the inversion does not
    reach such a power, the native result is checked by
    ``negative_base_guard``. Returns None to keep the native result."""
    from .. import (S, ConditionSet, Dummy, Eq, Equality, Expr, Lambda,
                    Union, exp, imageset, log, sympify)
    from .. import Set
    if not isinstance(domain, Set):
        return None
    if not getattr(symbol, 'is_Symbol', False) or domain.is_subset(S.Reals) is not True:
        return None
    f0 = sympify(f)
    if isinstance(f0, Equality):
        f0 = f0.lhs - f0.rhs
    if not isinstance(f0, Expr) or not _negative_base_pows(f0, symbol):
        return None
    unsolved = ConditionSet(symbol, Eq(f0, 0), domain)
    if f0.is_Mul and all(_finite_with_finite_vars(m, domain) for m in f0.args):
        # ``_solveset`` solves the factors of a product of finite terms
        return Union(*[solver(m, symbol, domain) for m in f0.args])
    lhs, ys = f0, [S.Zero]
    while True:
        if isinstance(lhs, exp):
            lhs, ys = lhs.args[0], [log(y) for y in ys]
        elif isinstance(lhs, log) and len(lhs.args) == 1:
            lhs, ys = lhs.args[0], [exp(y) for y in ys]
        elif lhs.is_Add and lhs.as_independent(symbol)[0] != 0:
            g, lhs = lhs.as_independent(symbol)
            ys = [y - g for y in ys]
        elif lhs.is_Mul and lhs.as_independent(symbol)[0] != 1:
            g, lhs = lhs.as_independent(symbol)
            if g.is_zero is not False or g.is_finite is not True:
                return None
            ys = [y/g for y in ys]
        else:
            break
    if not (lhs.is_Pow and not lhs.base.has(symbol) and lhs.exp.has(symbol)):
        return None
    if not lhs.base.is_negative:
        return unsolved if lhs.base.is_positive is not True and lhs.base.is_zero is not True else None
    expo = lhs.exp
    parts = []
    for y in ys:
        found = _negative_base_exponents(lhs.base, y)
        if found is None:
            return unsolved
        kind, value = found
        if kind == 'ints':
            parts.extend(solver(expo - k, symbol, domain) for k in value)
            continue
        c = expo.as_independent(symbol, as_Add=True)[0]
        coeff = (expo - c)/symbol
        if coeff.has(symbol) or coeff.is_extended_real is not True or coeff.is_zero is not False \
                or c.is_extended_real is not True:
            return unsolved
        n = Dummy('n', integer=True)
        lattice = imageset(Lambda(n, (2*n + value - c)/coeff), S.Integers)
        parts.append(lattice if domain == S.Reals else domain.intersect(lattice))
    return Union(*parts) if parts else S.EmptySet


def negative_base_guard(f, symbol, domain, result):
    """The native solver inverts ``b**g(x)`` with the principal logarithm,
    which misses the integer exponents of a negative base on a real domain
    (``(-2)**x*2**x = 16`` has the solution 2). An EmptySet obtained that
    way is not a proof that there is no solution, and a finite set may
    miss solutions (``(-2)**x = 2**x`` holds for every even ``x``), so the
    native result becomes a ConditionSet. ``negative_base_step`` already
    split the products and inverted the single powers; for the remaining
    equations SymPy either raises ValueError (``integer_log``) or keeps a
    ConditionSet."""
    from .. import S, ConditionSet, Eq, Equality, Expr, sympify
    from .. import Set
    if not isinstance(domain, Set):
        return result
    if not getattr(symbol, 'is_Symbol', False) or isinstance(result, ConditionSet):
        return result
    if domain.is_subset(S.Reals) is not True:
        return result
    f0 = sympify(f)
    if isinstance(f0, Equality):
        f0 = f0.lhs - f0.rhs
    if not isinstance(f0, Expr) or not _negative_base_pows(f0, symbol):
        return result
    return ConditionSet(symbol, Eq(f0, 0), domain)
