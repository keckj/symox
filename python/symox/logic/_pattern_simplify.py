"""Pattern-based simplification of ``And`` and ``Or``.

Port of SymPy 1.14 ``And._eval_simplify`` (``Eq`` substitution),
``Or._eval_simplify`` (``Le``/``Ge`` split into ``Eq`` and strict) and
``_apply_patternbased_simplification`` with the two- and three-term
pattern tables of ``sympy/logic/boolalg.py``. The ``Xor`` table is
ported but ``Xor._eval_simplify`` is not wired yet.
"""
from collections import defaultdict
from functools import lru_cache
from itertools import combinations, permutations, product

from symox import S, Tuple, Add, Not, And, Or, ITE, true, false
from symox.utilities.iterables import sift
from symox.utilities.sorting import ordered

cacheit = lru_cache(maxsize=None)


def _canonical(cond):
    from symox import Relational
    reps = {r: r.canonical for r in cond.atoms(Relational)}
    return cond.xreplace(reps)


def _is_inequality(r):
    return r.rel_op in ("<", "<=", ">", ">=")


def _is_rel(r):
    return getattr(r, "is_Relational", False)


def _is_eq(r):
    return _is_rel(r) and r.rel_op == "=="


def _apply_patternbased_simplification(rv, patterns, measure,
                                       dominatingvalue,
                                       replacementvalue=None,
                                       threeterm_patterns=None):
    """
    Replace patterns of Relational

    Parameters
    ==========

    rv : Expr
        Boolean expression

    patterns : tuple
        Tuple of tuples, with (pattern to simplify, simplified pattern) with
        two terms.

    measure : function
        Simplification measure.

    dominatingvalue : Boolean or ``None``
        The dominating value for the function of consideration.
        For example, for :py:class:`~.And` ``S.false`` is dominating.
        As soon as one expression is ``S.false`` in :py:class:`~.And`,
        the whole expression is ``S.false``.

    replacementvalue : Boolean or ``None``, optional
        The resulting value for the whole expression if one argument
        evaluates to ``dominatingvalue``.
        For example, for :py:class:`~.Nand` ``S.false`` is dominating, but
        in this case the resulting value is ``S.true``. Default is ``None``.
        If ``replacementvalue`` is ``None`` and ``dominatingvalue`` is not
        ``None``, ``replacementvalue = dominatingvalue``.

    threeterm_patterns : tuple, optional
        Tuple of tuples, with (pattern to simplify, simplified pattern) with
        three terms.

    """
    pass

    if replacementvalue is None and dominatingvalue is not None:
        replacementvalue = dominatingvalue
    # Use replacement patterns for Relationals
    Rel, nonRel = sift(rv.args, _is_rel,
                       binary=True)
    if len(Rel) <= 1:
        return rv
    Rel, nonRealRel = sift(Rel, lambda i: not any(s.is_real is False
                                                  for s in i.free_symbols),
                           binary=True)
    Rel = [i.canonical for i in Rel]

    if threeterm_patterns and len(Rel) >= 3:
        Rel = _apply_patternbased_threeterm_simplification(Rel,
                            threeterm_patterns, rv.func, dominatingvalue,
                            replacementvalue, measure)

    Rel = _apply_patternbased_twoterm_simplification(Rel, patterns,
                    rv.func, dominatingvalue, replacementvalue, measure)

    rv = rv.func(*([_canonical(i) for i in ordered(Rel)]
                 + nonRel + nonRealRel))
    return rv


def _apply_patternbased_twoterm_simplification(Rel, patterns, func,
                                               dominatingvalue,
                                               replacementvalue,
                                               measure):
    """ Apply pattern-based two-term simplification."""
    from symox import Min, Max
    pass
    changed = True
    while changed and len(Rel) >= 2:
        changed = False
        # Use only < or <=
        Rel = [r.reversed if (r.rel_op in (">", ">=")) else r for r in Rel]
        # Sort based on ordered
        Rel = list(ordered(Rel))
        # Eq and Ne must be tested reversed as well
        rtmp = [(r, ) if _is_inequality(r) else (r, r.reversed) for r in Rel]
        # Create a list of possible replacements
        results = []
        # Try all combinations of possibly reversed relational
        for ((i, pi), (j, pj)) in combinations(enumerate(rtmp), 2):
            for pattern, simp in patterns:
                res = []
                for p1, p2 in product(pi, pj):
                    # use SymPy matching
                    oldexpr = Tuple(p1, p2)
                    tmpres = oldexpr.match(pattern)
                    if tmpres:
                        res.append((tmpres, oldexpr))

                if res:
                    for tmpres, oldexpr in res:
                        # we have a matching, compute replacement
                        np = simp.xreplace(tmpres)
                        if np == dominatingvalue:
                            # if dominatingvalue, the whole expression
                            # will be replacementvalue
                            return [replacementvalue]
                        # add replacement
                        if not (getattr(np, "func", None) == ITE) and not np.has(Min, Max):
                            # We only want to use ITE and Min/Max replacements if
                            # they simplify to a relational
                            costsaving = measure(func(*oldexpr.args)) - measure(np)
                            if costsaving > 0:
                                results.append((costsaving, ([i, j], np)))
        if results:
            # Sort results based on complexity
            results = sorted(results,
                                           key=lambda pair: pair[0], reverse=True)
            # Replace the one providing most simplification
            replacement = results[0][1]
            idx, newrel = replacement
            idx.sort()
            # Remove the old relationals
            for index in reversed(idx):
                del Rel[index]
            if dominatingvalue is None or newrel != Not(dominatingvalue):
                # Insert the new one (no need to insert a value that will
                # not affect the result)
                if newrel.func == func:
                    for a in newrel.args:
                        Rel.append(a)
                else:
                    Rel.append(newrel)
            # We did change something so try again
            changed = True
    return Rel


def _apply_patternbased_threeterm_simplification(Rel, patterns, func,
                                                 dominatingvalue,
                                                 replacementvalue,
                                                 measure):
    """ Apply pattern-based three-term simplification."""
    from symox import Min, Max
    pass
    changed = True
    while changed and len(Rel) >= 3:
        changed = False
        # Use only > or >=
        Rel = [r.reversed if (r.rel_op in ("<", "<=")) else r for r in Rel]
        # Sort based on ordered
        Rel = list(ordered(Rel))
        # Create a list of possible replacements
        results = []
        # Eq and Ne must be tested reversed as well
        rtmp = [(r, ) if _is_inequality(r) else (r, r.reversed) for r in Rel]
        # Try all combinations of possibly reversed relational
        for ((i, pi), (j, pj), (k, pk)) in permutations(enumerate(rtmp), 3):
            for pattern, simp in patterns:
                res = []
                for p1, p2, p3 in product(pi, pj, pk):
                    # use SymPy matching
                    oldexpr = Tuple(p1, p2, p3)
                    tmpres = oldexpr.match(pattern)
                    if tmpres:
                        res.append((tmpres, oldexpr))

                if res:
                    for tmpres, oldexpr in res:
                        # we have a matching, compute replacement
                        np = simp.xreplace(tmpres)
                        if np == dominatingvalue:
                            # if dominatingvalue, the whole expression
                            # will be replacementvalue
                            return [replacementvalue]
                        # add replacement
                        if not (getattr(np, "func", None) == ITE) and not np.has(Min, Max):
                            # We only want to use ITE and Min/Max replacements if
                            # they simplify to a relational
                            costsaving = measure(func(*oldexpr.args)) - measure(np)
                            if costsaving > 0:
                                results.append((costsaving, ([i, j, k], np)))
        if results:
            # Sort results based on complexity
            results = sorted(results,
                                           key=lambda pair: pair[0], reverse=True)
            # Replace the one providing most simplification
            replacement = results[0][1]
            idx, newrel = replacement
            idx.sort()
            # Remove the old relationals
            for index in reversed(idx):
                del Rel[index]
            if dominatingvalue is None or newrel != Not(dominatingvalue):
                # Insert the new one (no need to insert a value that will
                # not affect the result)
                if newrel.func == func:
                    for a in newrel.args:
                        Rel.append(a)
                else:
                    Rel.append(newrel)
            # We did change something so try again
            changed = True
    return Rel


@cacheit
def _simplify_patterns_and():
    """ Two-term patterns for And."""

    from symox import Wild
    from symox import Eq, Ne, Ge, Gt, Le, Lt
    from symox import Abs
    from symox import Min, Max
    a = Wild('a')
    b = Wild('b')
    c = Wild('c')
    # Relationals patterns should be in alphabetical order
    # (pattern1, pattern2, simplified)
    # Do not use Ge, Gt
    _matchers_and = ((Tuple(Eq(a, b), Lt(a, b)), false),
                     #(Tuple(Eq(a, b), Lt(b, a)), S.false),
                     #(Tuple(Le(b, a), Lt(a, b)), S.false),
                     #(Tuple(Lt(b, a), Le(a, b)), S.false),
                     (Tuple(Lt(b, a), Lt(a, b)), false),
                     (Tuple(Eq(a, b), Le(b, a)), Eq(a, b)),
                     #(Tuple(Eq(a, b), Le(a, b)), Eq(a, b)),
                     #(Tuple(Le(b, a), Lt(b, a)), Gt(a, b)),
                     (Tuple(Le(b, a), Le(a, b)), Eq(a, b)),
                     #(Tuple(Le(b, a), Ne(a, b)), Gt(a, b)),
                     #(Tuple(Lt(b, a), Ne(a, b)), Gt(a, b)),
                     (Tuple(Le(a, b), Lt(a, b)), Lt(a, b)),
                     (Tuple(Le(a, b), Ne(a, b)), Lt(a, b)),
                     (Tuple(Lt(a, b), Ne(a, b)), Lt(a, b)),
                     # Sign
                     (Tuple(Eq(a, b), Eq(a, -b)), And(Eq(a, S.Zero), Eq(b, S.Zero))),
                     # Min/Max/ITE
                     (Tuple(Le(b, a), Le(c, a)), Ge(a, Max(b, c))),
                     (Tuple(Le(b, a), Lt(c, a)), ITE(b > c, Ge(a, b), Gt(a, c))),
                     (Tuple(Lt(b, a), Lt(c, a)), Gt(a, Max(b, c))),
                     (Tuple(Le(a, b), Le(a, c)), Le(a, Min(b, c))),
                     (Tuple(Le(a, b), Lt(a, c)), ITE(b < c, Le(a, b), Lt(a, c))),
                     (Tuple(Lt(a, b), Lt(a, c)), Lt(a, Min(b, c))),
                     (Tuple(Le(a, b), Le(c, a)), ITE(Eq(b, c), Eq(a, b), ITE(b < c, false, And(Le(a, b), Ge(a, c))))),
                     (Tuple(Le(c, a), Le(a, b)), ITE(Eq(b, c), Eq(a, b), ITE(b < c, false, And(Le(a, b), Ge(a, c))))),
                     (Tuple(Lt(a, b), Lt(c, a)), ITE(b < c, false, And(Lt(a, b), Gt(a, c)))),
                     (Tuple(Lt(c, a), Lt(a, b)), ITE(b < c, false, And(Lt(a, b), Gt(a, c)))),
                     (Tuple(Le(a, b), Lt(c, a)), ITE(b <= c, false, And(Le(a, b), Gt(a, c)))),
                     (Tuple(Le(c, a), Lt(a, b)), ITE(b <= c, false, And(Lt(a, b), Ge(a, c)))),
                     (Tuple(Eq(a, b), Eq(a, c)), ITE(Eq(b, c), Eq(a, b), false)),
                     (Tuple(Lt(a, b), Lt(-b, a)), ITE(b > 0, Lt(Abs(a), b), false)),
                     (Tuple(Le(a, b), Le(-b, a)), ITE(b >= 0, Le(Abs(a), b), false)),
                     )
    return _matchers_and


@cacheit
def _simplify_patterns_and3():
    """ Three-term patterns for And."""

    from symox import Wild
    from symox import Eq, Ge, Gt

    a = Wild('a')
    b = Wild('b')
    c = Wild('c')
    # Relationals patterns should be in alphabetical order
    # (pattern1, pattern2, pattern3, simplified)
    # Do not use Le, Lt
    _matchers_and = ((Tuple(Ge(a, b), Ge(b, c), Gt(c, a)), false),
                     (Tuple(Ge(a, b), Gt(b, c), Gt(c, a)), false),
                     (Tuple(Gt(a, b), Gt(b, c), Gt(c, a)), false),
                     # (Tuple(Ge(c, a), Gt(a, b), Gt(b, c)), S.false),
                     # Lower bound relations
                     # Commented out combinations that does not simplify
                     (Tuple(Ge(a, b), Ge(a, c), Ge(b, c)), And(Ge(a, b), Ge(b, c))),
                     (Tuple(Ge(a, b), Ge(a, c), Gt(b, c)), And(Ge(a, b), Gt(b, c))),
                     # (Tuple(Ge(a, b), Gt(a, c), Ge(b, c)), And(Ge(a, b), Ge(b, c))),
                     (Tuple(Ge(a, b), Gt(a, c), Gt(b, c)), And(Ge(a, b), Gt(b, c))),
                     # (Tuple(Gt(a, b), Ge(a, c), Ge(b, c)), And(Gt(a, b), Ge(b, c))),
                     (Tuple(Ge(a, c), Gt(a, b), Gt(b, c)), And(Gt(a, b), Gt(b, c))),
                     (Tuple(Ge(b, c), Gt(a, b), Gt(a, c)), And(Gt(a, b), Ge(b, c))),
                     (Tuple(Gt(a, b), Gt(a, c), Gt(b, c)), And(Gt(a, b), Gt(b, c))),
                     # Upper bound relations
                     # Commented out combinations that does not simplify
                     (Tuple(Ge(b, a), Ge(c, a), Ge(b, c)), And(Ge(c, a), Ge(b, c))),
                     (Tuple(Ge(b, a), Ge(c, a), Gt(b, c)), And(Ge(c, a), Gt(b, c))),
                     # (Tuple(Ge(b, a), Gt(c, a), Ge(b, c)), And(Gt(c, a), Ge(b, c))),
                     (Tuple(Ge(b, a), Gt(c, a), Gt(b, c)), And(Gt(c, a), Gt(b, c))),
                     # (Tuple(Gt(b, a), Ge(c, a), Ge(b, c)), And(Ge(c, a), Ge(b, c))),
                     (Tuple(Ge(c, a), Gt(b, a), Gt(b, c)), And(Ge(c, a), Gt(b, c))),
                     (Tuple(Ge(b, c), Gt(b, a), Gt(c, a)), And(Gt(c, a), Ge(b, c))),
                     (Tuple(Gt(b, a), Gt(c, a), Gt(b, c)), And(Gt(c, a), Gt(b, c))),
                     # Circular relation
                     (Tuple(Ge(a, b), Ge(b, c), Ge(c, a)), And(Eq(a, b), Eq(b, c))),
                     )
    return _matchers_and


@cacheit
def _simplify_patterns_or():
    """ Two-term patterns for Or."""

    from symox import Wild
    from symox import Eq, Ne, Ge, Gt, Le, Lt
    from symox import Abs
    from symox import Min, Max
    a = Wild('a')
    b = Wild('b')
    c = Wild('c')
    # Relationals patterns should be in alphabetical order
    # (pattern1, pattern2, simplified)
    # Do not use Ge, Gt
    _matchers_or = ((Tuple(Le(b, a), Le(a, b)), true),
                    #(Tuple(Le(b, a), Lt(a, b)), true),
                    (Tuple(Le(b, a), Ne(a, b)), true),
                    #(Tuple(Le(a, b), Lt(b, a)), true),
                    #(Tuple(Le(a, b), Ne(a, b)), true),
                    #(Tuple(Eq(a, b), Le(b, a)), Ge(a, b)),
                    #(Tuple(Eq(a, b), Lt(b, a)), Ge(a, b)),
                    (Tuple(Eq(a, b), Le(a, b)), Le(a, b)),
                    (Tuple(Eq(a, b), Lt(a, b)), Le(a, b)),
                    #(Tuple(Le(b, a), Lt(b, a)), Ge(a, b)),
                    (Tuple(Lt(b, a), Lt(a, b)), Ne(a, b)),
                    (Tuple(Lt(b, a), Ne(a, b)), Ne(a, b)),
                    (Tuple(Le(a, b), Lt(a, b)), Le(a, b)),
                    #(Tuple(Lt(a, b), Ne(a, b)), Ne(a, b)),
                    (Tuple(Eq(a, b), Ne(a, c)), ITE(Eq(b, c), true, Ne(a, c))),
                    (Tuple(Ne(a, b), Ne(a, c)), ITE(Eq(b, c), Ne(a, b), true)),
                    # Min/Max/ITE
                    (Tuple(Le(b, a), Le(c, a)), Ge(a, Min(b, c))),
                    #(Tuple(Ge(b, a), Ge(c, a)), Ge(Min(b, c), a)),
                    (Tuple(Le(b, a), Lt(c, a)), ITE(b > c, Lt(c, a), Le(b, a))),
                    (Tuple(Lt(b, a), Lt(c, a)), Gt(a, Min(b, c))),
                    #(Tuple(Gt(b, a), Gt(c, a)), Gt(Min(b, c), a)),
                    (Tuple(Le(a, b), Le(a, c)), Le(a, Max(b, c))),
                    #(Tuple(Le(b, a), Le(c, a)), Le(Max(b, c), a)),
                    (Tuple(Le(a, b), Lt(a, c)), ITE(b >= c, Le(a, b), Lt(a, c))),
                    (Tuple(Lt(a, b), Lt(a, c)), Lt(a, Max(b, c))),
                    #(Tuple(Lt(b, a), Lt(c, a)), Lt(Max(b, c), a)),
                    (Tuple(Le(a, b), Le(c, a)), ITE(b >= c, true, Or(Le(a, b), Ge(a, c)))),
                    (Tuple(Le(c, a), Le(a, b)), ITE(b >= c, true, Or(Le(a, b), Ge(a, c)))),
                    (Tuple(Lt(a, b), Lt(c, a)), ITE(b > c, true, Or(Lt(a, b), Gt(a, c)))),
                    (Tuple(Lt(c, a), Lt(a, b)), ITE(b > c, true, Or(Lt(a, b), Gt(a, c)))),
                    (Tuple(Le(a, b), Lt(c, a)), ITE(b >= c, true, Or(Le(a, b), Gt(a, c)))),
                    (Tuple(Le(c, a), Lt(a, b)), ITE(b >= c, true, Or(Lt(a, b), Ge(a, c)))),
                    (Tuple(Lt(b, a), Lt(a, -b)), ITE(b >= 0, Gt(Abs(a), b), true)),
                    (Tuple(Le(b, a), Le(a, -b)), ITE(b > 0, Ge(Abs(a), b), true)),
                    )
    return _matchers_or


@cacheit
def _simplify_patterns_xor():
    """ Two-term patterns for Xor."""

    from symox import Min, Max
    from symox import Wild
    from symox import Eq, Ne, Ge, Gt, Le, Lt
    a = Wild('a')
    b = Wild('b')
    c = Wild('c')
    # Relationals patterns should be in alphabetical order
    # (pattern1, pattern2, simplified)
    # Do not use Ge, Gt
    _matchers_xor = (#(Tuple(Le(b, a), Lt(a, b)), true),
                     #(Tuple(Lt(b, a), Le(a, b)), true),
                     #(Tuple(Eq(a, b), Le(b, a)), Gt(a, b)),
                     #(Tuple(Eq(a, b), Lt(b, a)), Ge(a, b)),
                     (Tuple(Eq(a, b), Le(a, b)), Lt(a, b)),
                     (Tuple(Eq(a, b), Lt(a, b)), Le(a, b)),
                     (Tuple(Le(a, b), Lt(a, b)), Eq(a, b)),
                     (Tuple(Le(a, b), Le(b, a)), Ne(a, b)),
                     (Tuple(Le(b, a), Ne(a, b)), Le(a, b)),
                     # (Tuple(Lt(b, a), Lt(a, b)), Ne(a, b)),
                     (Tuple(Lt(b, a), Ne(a, b)), Lt(a, b)),
                     # (Tuple(Le(a, b), Lt(a, b)), Eq(a, b)),
                     # (Tuple(Le(a, b), Ne(a, b)), Ge(a, b)),
                     # (Tuple(Lt(a, b), Ne(a, b)), Gt(a, b)),
                     # Min/Max/ITE
                     (Tuple(Le(b, a), Le(c, a)),
                      And(Ge(a, Min(b, c)), Lt(a, Max(b, c)))),
                     (Tuple(Le(b, a), Lt(c, a)),
                      ITE(b > c, And(Gt(a, c), Lt(a, b)),
                          And(Ge(a, b), Le(a, c)))),
                     (Tuple(Lt(b, a), Lt(c, a)),
                      And(Gt(a, Min(b, c)), Le(a, Max(b, c)))),
                     (Tuple(Le(a, b), Le(a, c)),
                      And(Le(a, Max(b, c)), Gt(a, Min(b, c)))),
                     (Tuple(Le(a, b), Lt(a, c)),
                      ITE(b < c, And(Lt(a, c), Gt(a, b)),
                          And(Le(a, b), Ge(a, c)))),
                     (Tuple(Lt(a, b), Lt(a, c)),
                      And(Lt(a, Max(b, c)), Ge(a, Min(b, c)))),
                     )
    return _matchers_xor




def and_simplify(rv, measure, ratio):
    """The part of ``And._eval_simplify`` after the standard simplify."""
    from symox.solvers.solveset import linear_coeffs
    Rel, nonRel = sift(rv.args, _is_rel, binary=True)
    if not Rel:
        return rv
    eqs, other = sift(Rel, _is_eq, binary=True)
    if eqs:
        reps = {}
        sifted = sift(ordered([
            (i.free_symbols, i) for i in eqs]),
            lambda x: len(x[0]))
        eqs = []
        nonlineqs = []
        while 1 in sifted:
            for free, e in sifted.pop(1):
                x = free.pop()
                if (e.lhs != x or x in e.rhs.free_symbols) and x not in reps:
                    try:
                        m, b = linear_coeffs(
                            Add(e.lhs, -e.rhs, evaluate=False), x)
                        enew = e.func(x, -b/m)
                        if measure(enew) <= ratio*measure(e):
                            e = enew
                        else:
                            eqs.append(e)
                            continue
                    except ValueError:
                        pass
                if x in reps:
                    eqs.append(e.subs(x, reps[x]))
                elif e.lhs == x and x not in e.rhs.free_symbols:
                    reps[x] = e.rhs
                    eqs.append(e)
                else:
                    # x is not yet identified, but may be later
                    nonlineqs.append(e)
            resifted = defaultdict(list)
            for k in sifted:
                for f, e in sifted[k]:
                    e = e.xreplace(reps)
                    f = e.free_symbols
                    resifted[len(f)].append((f, e))
            sifted = resifted
        for k in sifted:
            eqs.extend([e for f, e in sifted[k]])
        nonlineqs = [ei.subs(reps) for ei in nonlineqs]
        other = [ei.subs(reps) for ei in other]
        rv = rv.func(*([i.canonical for i in (eqs + nonlineqs + other)] + nonRel))
    if rv.func != And:
        return rv
    patterns = _simplify_patterns_and()
    threeterm_patterns = _simplify_patterns_and3()
    return _apply_patternbased_simplification(rv, patterns,
                                              measure, false,
                                              threeterm_patterns=threeterm_patterns)


def or_split_lege(expr):
    """``Or._eval_simplify``: ``a <= b`` becomes ``Eq(a, b) | (a < b)``."""
    from symox import Relational, Eq
    lege = [r for r in expr.atoms(Relational) if r.rel_op in ("<=", ">=")]
    if not lege:
        return None
    reps = {i: Or(Eq(i.lhs, i.rhs), i.strict) for i in lege}
    return expr.xreplace(reps)


def or_simplify(rv, measure):
    """The pattern part of ``Or._eval_simplify``."""
    return _apply_patternbased_simplification(rv, _simplify_patterns_or(),
                                              measure, true)
