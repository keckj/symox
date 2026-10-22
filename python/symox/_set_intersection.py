"""Port of the reference ``intersection_sets(ImageSet, Set)`` handler.

The native intersection keeps an ``ImageSet`` unevaluated unless its map is
an exact linear lattice over the integers or the other set is a bounded
numeric interval. This module adds the remaining cases of the reference
handler of SymPy 1.14 (``sympy/sets/handlers/intersection.py``): the
diophantine equation ``f(n) = g(m)`` against ``Integers`` or another integer
image set, the real part of the map against ``Reals``, and the inverted
bounds of an ``Interval``.
"""


def _imageset_intersection(self, other):
    from symox import (S, Dummy, FiniteSet, Lambda, Union, ConditionSet, Eq,
                       Interval, ImageSet, Range, Mul, imageset, Intersection)
    from symox import expand_complex
    lamda = self.lamda
    if len(lamda.variables) > 1 or tuple(lamda.signature) != tuple(lamda.variables):
        return None
    base_set = self.base_sets[0]

    if base_set == S.Integers:
        gm = None
        if isinstance(other, ImageSet) and other.base_sets == (S.Integers,):
            gm = other.lamda.expr
            var = other.lamda.variables[0]
            m = Dummy('m')
            gm = gm.subs(var, m)
        elif other == S.Integers:
            m = gm = Dummy('m')
        if gm is not None:
            from symox import diophantine
            fn = lamda.expr
            n = lamda.variables[0]
            try:
                solns = list(diophantine(fn - gm, syms=(n, m), permute=True))
            except (TypeError, NotImplementedError, ValueError):
                return None
            if len(solns) == 0:
                return S.EmptySet
            if any(s.free_symbols for tupl in solns for s in tupl):
                if len(solns) == 1:
                    soln, _ = solns[0]
                    (t,) = soln.free_symbols
                    expr = fn.subs(n, soln.subs(t, n)).expand()
                    return imageset(Lambda(n, expr), S.Integers)
                return None
            return FiniteSet(*(fn.subs(n, s[0]) for s in solns))

    if other == S.Reals:
        from symox.solvers.solvers import denoms
        from symox import solve_linear, numer

        def _solution_union(exprs, sym):
            sols = []
            for i in exprs:
                x, xis = solve_linear(i, 0, [sym])
                if x == sym:
                    sols.append(FiniteSet(xis))
                else:
                    sols.append(ConditionSet(sym, Eq(i, 0)))
            return Union(*sols)

        f = lamda.expr
        n = lamda.variables[0]
        n_ = Dummy(n.name, real=True)
        f_ = f.subs(n, n_)
        re, im = f_.as_real_imag()
        im = expand_complex(im)
        re = re.subs(n_, n)
        im = im.subs(n_, n)
        ifree = im.free_symbols
        lam = Lambda(n, re)
        if im.is_zero:
            pass
        elif im.is_zero is False:
            return S.EmptySet
        elif ifree != {n}:
            return None
        else:
            base_set &= _solution_union(Mul.make_args(numer(im)), n)
        base_set -= _solution_union(denoms(f), n)
        return imageset(lam, base_set)

    if isinstance(other, Interval):
        from symox.solvers.solveset import invert_real, invert_complex
        f = lamda.expr
        n = lamda.variables[0]

        def _through_real_part():
            # the reference intersects the real zeros of f with the interval,
            # which drops real values such as 1 in ImageSet(n*I + 1); take the
            # real part of the image set first instead
            if f.is_real or other.is_subset(S.Reals) is not True:
                return None
            real_part = _imageset_intersection(self, S.Reals)
            if real_part is None or real_part == self:
                return None
            return Intersection(real_part, other)

        new_inf, new_sup = None, None
        new_lopen, new_ropen = other.left_open, other.right_open
        inverter = invert_real if f.is_real else invert_complex
        g1, h1 = inverter(f, other.inf, n)
        g2, h2 = inverter(f, other.sup, n)
        if not all(isinstance(i, FiniteSet) for i in (h1, h2)):
            return _through_real_part()
        if g1 == n and len(h1) == 1:
            new_inf = h1.args[0]
        if g2 == n and len(h2) == 1:
            new_sup = h2.args[0]
        if any(i is None for i in (new_sup, new_inf)):
            return _through_real_part()
        if not all(i.is_real for i in (new_sup, new_inf)):
            return _through_real_part()
        if new_inf > new_sup:
            new_inf, new_sup = new_sup, new_inf
        new_interval = Interval(new_inf, new_sup, new_lopen, new_ropen)
        range_set = base_set.intersect(new_interval)
        if range_set == S.EmptySet:
            return S.EmptySet
        if isinstance(range_set, Range) and range_set.size is not S.Infinity:
            range_set = FiniteSet(*list(range_set))
        if range_set is not None:
            return imageset(Lambda(n, f), range_set)
    return None


_ACTIVE = []


def refine_intersection(result):
    """Apply the reference image-set rules to an unevaluated intersection."""
    from symox import ImageSet, Intersection, S
    if not isinstance(result, Intersection):
        return result
    args = list(result.args)
    if not any(isinstance(a, ImageSet) or a == S.Complexes for a in args):
        return result
    key = result
    if key in _ACTIVE:
        return result
    _ACTIVE.append(key)
    try:
        if S.Complexes in args:
            # reference `intersection_sets(ComplexRegion, Set)` for the
            # rectangular `Complexes`: a real set meets its real axis
            others = [a for a in args if a != S.Complexes]
            for k, t in enumerate(others):
                if t.is_subset(S.Reals) is True:
                    others[k] = Intersection(S.Reals, t)
                    return Intersection(*others)
        for i, s in enumerate(args):
            if not isinstance(s, ImageSet):
                continue
            for j, t in enumerate(args):
                if i == j:
                    continue
                try:
                    new = _imageset_intersection(s, t)
                except (TypeError, NotImplementedError, ValueError, AttributeError):
                    new = None
                if new is None:
                    continue
                rest = [a for k, a in enumerate(args) if k not in (i, j)]
                if not rest:
                    return new
                return Intersection(new, *rest)
    finally:
        _ACTIVE.pop()
    return result
