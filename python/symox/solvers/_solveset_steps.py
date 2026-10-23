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
