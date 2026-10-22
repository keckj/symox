"""The substitution method of the reference ``solveset.substitution``.

Each equation is solved with ``solveset_real`` and ``solveset_complex`` for
one unknown at a time, and the solution is substituted in the remaining
equations. Periodic solutions are kept as ``ImageSet`` families: the lambda
expression is carried through the substitution, and the ``ImageSet`` is
restored in the final tuples. This follows SymPy 1.14 step by step, so
no family of solutions is dropped.
"""


def _return_conditionset(eqs, symbols):
    from .. import ConditionSet, Tuple, And, Eq, S
    eqs = (Eq(lhs, 0) for lhs in eqs)
    return ConditionSet(Tuple(*symbols), And(*eqs), S.Complexes**len(symbols))


def _intersection_args(s):
    """The two arguments of a binary Intersection in the reference order:
    the domain (Reals, an Interval, ...) first and the candidate solutions
    (FiniteSet, ImageSet) second. The native canonical order can differ."""
    from .. import FiniteSet, ImageSet, Union
    first, second = s.args[0], s.args[-1]
    if len(s.args) == 2 and isinstance(first, (FiniteSet, ImageSet, Union)) \
            and not isinstance(second, (FiniteSet, ImageSet, Union)):
        first, second = second, first
    if len(s.args) > 2:
        from .. import Intersection
        second = Intersection(*s.args[1:])
    return first, second


def _checksol(*args, **kwargs):
    from ..solvers.solvers import checksol
    return checksol(*args, **kwargs)


def substitution(system, symbols, result=None, known_symbols=(), exclude=(), all_symbols=None):
    from .. import (S, Dummy, Lambda, ImageSet, FiniteSet, ConditionSet, Complement,
                    Intersection, Union, Interval, oo, Expr, Abs, Eq, default_sort_key,
                    ordered, solveset_real, solveset_complex, sympify, simplify, Set)
    from ..utilities.iterables import is_sequence

    if not system:
        return S.EmptySet
    system = list(system)
    for i, e in enumerate(system):
        e = sympify(e)
        if isinstance(e, Eq) or getattr(e, 'is_Equality', False):
            e = e.lhs - e.rhs
        system[i] = e
    if not symbols:
        raise ValueError('Symbols must be given, for which solution of the system is to be found.')
    if not is_sequence(symbols):
        raise TypeError('symbols should be given as a sequence, e.g. a list.'
                        'Not type %s: %s' % (type(symbols), symbols))
    symbols = list(symbols)
    if not getattr(symbols[0], 'is_Symbol', False):
        raise ValueError('Iterable of symbols must be given as second argument, not type %s: %s'
                         % (type(symbols[0]), symbols[0]))
    if all_symbols is None:
        all_symbols = symbols
    all_symbols = list(all_symbols)
    known_symbols = list(known_symbols)
    exclude = [sympify(d) for d in exclude]

    # The reference mutates the dictionaries of ``result`` while it solves
    # (both passes share them); work on copies of the caller's data.
    old_result = [dict(r) for r in ([{}] if result is None else result)]
    complements = {}
    intersections = {}
    total_conditionset = -1
    total_solveset_call = -1

    def _unsolved_syms(eq, sort=False):
        free = eq.free_symbols
        unsolved = (free - set(known_symbols)) & set(all_symbols)
        if sort:
            unsolved = list(unsolved)
            unsolved.sort(key=default_sort_key)
        return unsolved

    eqs_in_better_order = list(ordered(system, lambda _: len(_unsolved_syms(_))))

    def add_intersection_complement(result, intersection_dict, complement_dict):
        final_result = []
        for res in result:
            res_copy = res
            for key_res, value_res in res.items():
                intersect_set, complement_set = None, None
                for key_sym, value_sym in intersection_dict.items():
                    if key_sym == key_res:
                        intersect_set = value_sym
                for key_sym, value_sym in complement_dict.items():
                    if key_sym == key_res:
                        complement_set = value_sym
                if intersect_set or complement_set:
                    new_value = FiniteSet(value_res)
                    if intersect_set and intersect_set != S.Complexes:
                        new_value = Intersection(new_value, intersect_set)
                    if complement_set:
                        new_value = Complement(new_value, complement_set)
                    if new_value is S.EmptySet or new_value == S.EmptySet:
                        res_copy = None
                        break
                    elif new_value.is_FiniteSet and len(new_value) == 1:
                        res_copy[key_res] = set(new_value).pop()
                    else:
                        res_copy[key_res] = new_value
            if res_copy is not None:
                final_result.append(res_copy)
        return final_result

    def _extract_main_soln(sym, sol, soln_imageset):
        if isinstance(sol, ConditionSet):
            sol = sol.base_set
        if isinstance(sol, Complement):
            complements[sym] = sol.args[1]
            sol = sol.args[0]
        if isinstance(sol, Union):
            sol_args = sol.args
            sol = S.EmptySet
            for sol_arg2 in sol_args:
                if isinstance(sol_arg2, FiniteSet):
                    sol += sol_arg2
                else:
                    sol += FiniteSet(sol_arg2)
        if isinstance(sol, Intersection):
            domain, sol = _intersection_args(sol)
            if domain not in (S.Reals, S.Complexes):
                intersections[sym] = domain
        if isinstance(sol, ImageSet):
            soln_imagest = sol
            expr2 = sol.lamda.expr
            sol = FiniteSet(expr2)
            soln_imageset[expr2] = soln_imagest
        if not isinstance(sol, FiniteSet):
            sol = FiniteSet(sol)
        return sol, soln_imageset

    def _check_exclude(rnew, imgset_yes):
        rnew_ = rnew
        if imgset_yes:
            rnew_copy = rnew.copy()
            dummy_n = imgset_yes[0]
            for key_res, value_res in rnew_copy.items():
                rnew_copy[key_res] = value_res.subs(dummy_n, 0)
            rnew_ = rnew_copy
        try:
            satisfy_exclude = any(_checksol(d, rnew_) for d in exclude)
        except TypeError:
            satisfy_exclude = None
        return satisfy_exclude

    def _restore_imgset(rnew, original_imageset, newresult):
        restore_sym = set(rnew.keys()) & set(original_imageset.keys())
        for key_sym in restore_sym:
            rnew[key_sym] = original_imageset[key_sym]
        if rnew not in newresult:
            newresult.append(rnew)

    def _append_eq(eq, result, res, delete_soln, n=None):
        u = Dummy('u')
        if n:
            eq = eq.subs(n, 0)
        satisfy = eq if eq in (True, False, S.true, S.false) else _checksol(u, u, eq, minimal=True)
        if satisfy is False or satisfy is S.false:
            delete_soln = True
            res = {}
        else:
            result.append(res)
        return result, res, delete_soln

    def _append_new_soln(rnew, sym, sol, imgset_yes, soln_imageset, original_imageset, newresult, eq=None):
        satisfy_exclude = _check_exclude(rnew, imgset_yes)
        delete_soln = False
        if not satisfy_exclude:
            local_n = None
            if imgset_yes:
                local_n = imgset_yes[0]
                base = imgset_yes[1]
                if sym and sol:
                    dummy_list = list(sol.atoms(Dummy))
                    dummy_zip = [(d, local_n) for d in dummy_list]
                    lam = Lambda(local_n, sol.subs(dummy_zip))
                    rnew[sym] = ImageSet(lam, base)
                if eq is not None:
                    newresult, rnew, delete_soln = _append_eq(eq, newresult, rnew, delete_soln, local_n)
            elif eq is not None:
                newresult, rnew, delete_soln = _append_eq(eq, newresult, rnew, delete_soln)
            elif sol in soln_imageset.keys():
                rnew[sym] = soln_imageset[sol]
                _restore_imgset(rnew, original_imageset, newresult)
            else:
                newresult.append(rnew)
        elif satisfy_exclude:
            delete_soln = True
            rnew = {}
        _restore_imgset(rnew, original_imageset, newresult)
        return newresult, delete_soln

    def _solve_using_known_values(result, solver):
        soln_imageset = {}
        total_solvest_call = 0
        total_conditionst = 0
        for index, eq in enumerate(eqs_in_better_order):
            newresult = []
            imgset_yes = False
            for res in list(result):
                original_imageset = {}
                got_symbol = set()
                for k, v in res.items():
                    if isinstance(v, ImageSet):
                        res[k] = v.lamda.expr
                        original_imageset[k] = v
                        dummy_n = v.lamda.expr.atoms(Dummy).pop()
                        (base,) = v.base_sets
                        imgset_yes = (dummy_n, base)
                eq2 = eq.subs(res).expand()
                if imgset_yes and not eq2.has(imgset_yes[0]):
                    imgset_yes = False
                unsolved_syms = _unsolved_syms(eq2, sort=True)
                if not unsolved_syms:
                    if res:
                        newresult, delete_res = _append_new_soln(
                            res, None, None, imgset_yes, soln_imageset,
                            original_imageset, newresult, eq2)
                        if delete_res and res in result:
                            result.remove(res)
                    continue
                depen1, depen2 = eq2.as_independent(*unsolved_syms)
                if (depen1.has(Abs) or depen2.has(Abs)) and solver is solveset_complex:
                    continue
                soln_imageset = {}
                for sym in unsolved_syms:
                    not_solvable = False
                    try:
                        soln = solver(eq2, sym)
                        total_solvest_call += 1
                        soln_new = S.EmptySet
                        if isinstance(soln, Complement):
                            complements[sym] = soln.args[1]
                            soln = soln.args[0]
                        if isinstance(soln, Intersection):
                            domain, rest = _intersection_args(soln)
                            if domain != Interval(-oo, oo):
                                intersections[sym] = domain
                            soln_new += rest
                        soln = soln_new if soln_new else soln
                        if index > 0 and solver is solveset_real:
                            if not isinstance(soln, (ImageSet, ConditionSet)):
                                soln += solveset_complex(eq2, sym)
                    except (NotImplementedError, ValueError):
                        continue
                    if isinstance(soln, ConditionSet):
                        if soln.base_set in (S.Reals, S.Complexes):
                            soln = S.EmptySet
                            not_solvable = True
                            total_conditionst += 1
                        else:
                            soln = soln.base_set
                    if soln is not S.EmptySet and soln != S.EmptySet:
                        soln, soln_imageset = _extract_main_soln(sym, soln, soln_imageset)
                    for sol in soln:
                        sol, soln_imageset = _extract_main_soln(sym, sol, soln_imageset)
                        sol = set(sol).pop()
                        free = sol.free_symbols
                        if got_symbol and any(ss in free for ss in got_symbol):
                            continue
                        rnew = res.copy()
                        for k, v in res.items():
                            if isinstance(v, Expr) and isinstance(sol, Expr):
                                rnew[k] = v.subs(sym, sol)
                        if sol in soln_imageset.keys():
                            imgst = soln_imageset[sol]
                            rnew[sym] = imgst.lamda(*[0 for i in range(len(imgst.lamda.variables))])
                        else:
                            rnew[sym] = sol
                        newresult, delete_res = _append_new_soln(
                            rnew, sym, sol, imgset_yes, soln_imageset,
                            original_imageset, newresult)
                        if delete_res and res in result:
                            result.remove(res)
                    if not not_solvable:
                        got_symbol.add(sym)
            if newresult:
                result = newresult
        return result, total_solvest_call, total_conditionst

    new_result_real, solve_call1, cnd_call1 = _solve_using_known_values(old_result, solveset_real)
    new_result_complex, solve_call2, cnd_call2 = _solve_using_known_values(old_result, solveset_complex)
    total_conditionset += (cnd_call1 + cnd_call2)
    total_solveset_call += (solve_call1 + solve_call2)
    if total_conditionset == total_solveset_call and total_solveset_call != -1:
        return _return_conditionset(eqs_in_better_order, all_symbols)

    def _same(a, b):
        if isinstance(a, int) and isinstance(b, int):
            return True
        if isinstance(a, ImageSet) and isinstance(b, ImageSet):
            # Basic.dummy_eq: equal up to the renaming of the bound variables.
            va, vb = a.lamda.variables, b.lamda.variables
            return (len(va) == len(vb) and a.base_sets == b.base_sets
                    and a.lamda.expr.xreplace(dict(zip(va, vb))) == b.lamda.expr)
        try:
            return a.dummy_eq(b)
        except AttributeError:
            return a == b

    filtered_complex = []
    for i in list(new_result_complex):
        for j in list(new_result_real):
            if i.keys() != j.keys():
                continue
            if all(_same(a, b) for a, b in zip(i.values(), j.values())):
                break
        else:
            filtered_complex.append(i)
    result = new_result_real + filtered_complex

    def _same_row(a, b):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)

    def _plain(value):
        return value.lamda.expr if isinstance(value, ImageSet) else value

    def _constrains(row, unsolved):
        # An equation that still depends on an unknown left free in ``row``
        # is a constraint the substitution did not solve (solveset returned
        # a ConditionSet): the unknown is not a free parameter.
        values = {k: _plain(v) for k, v in row.items() if k not in unsolved}
        for eq in system:
            residual = eq.subs(values)
            if residual.free_symbols & set(unsolved) and simplify(residual) != 0:
                return True
        return False

    def _covers(row, finite):
        # Whether the parametric ``row`` contains the point ``finite``.
        free = {k: finite[k] for k, v in row.items() if v == k}
        for k, v in row.items():
            if k in free or isinstance(v, Set) or isinstance(finite[k], Set):
                if k not in free and not _same(v, finite[k]):
                    return False
                continue
            if simplify(v.subs(free) - finite[k]) != 0:
                return False
        return True

    result_all_variables = []
    result_infinite = []
    unresolved = False
    for res in result:
        if not res:
            continue
        if len(res) < len(all_symbols):
            solved_symbols = res.keys()
            unsolved = [s for s in all_symbols if s not in solved_symbols]
            for unsolved_sym in unsolved:
                res[unsolved_sym] = unsolved_sym
            if _constrains(res, unsolved):
                if cnd_call1 + cnd_call2:
                    unresolved = True
                else:
                    # Every solveset call was decided: this row only comes
                    # from a pass in which an equation had no solution
                    # (over the reals) and it is not a solution.
                    continue
            result_infinite.append(res)
        if not any(_same_row(res, other) for other in result_all_variables):
            result_all_variables.append(res)
    if unresolved:
        # The reference returns the row with the unknown left free, which
        # claims every value of that unknown, although solveset could not
        # solve an equation for it; return the unsolved system.
        return _return_conditionset(eqs_in_better_order, all_symbols)
    if result_infinite:
        # The reference keeps only the parametric rows; also keep a finite
        # row unless a parametric row contains it.
        kept = []
        for res in result_infinite:
            if not any(_same_row(res, other) for other in kept):
                kept.append(res)
        for res in result_all_variables:
            if any(res is other for other in result_infinite):
                continue
            if not any(_covers(other, res) for other in result_infinite):
                kept.append(res)
        result_all_variables = kept
    if intersections or complements:
        result_all_variables = add_intersection_complement(result_all_variables, intersections, complements)
    out = S.EmptySet
    for r in result_all_variables:
        out += FiniteSet(tuple(r[symb] for symb in all_symbols))
    return out


def _simple_dens(f, symbols):
    from ..solvers.solvers import denoms
    dens = set()
    for d in denoms(f, symbols):
        if d.is_Pow and d.exp.is_Number:
            if d.exp.is_zero:
                continue
            d = d.base
        dens.add(d)
    return dens


def _separate_poly_nonpoly(system, symbols):
    from .. import simplify, Eq
    from ..solvers.solvers import unrad
    polys, polys_expr, nonpolys, unrad_changed = [], [], [], []
    denominators = set()
    poly = None
    for eq in system:
        denominators.update(_simple_dens(eq, symbols))
        if isinstance(eq, Eq) or getattr(eq, 'is_Equality', False):
            eq = eq.lhs - eq.rhs
        without_radicals = unrad(simplify(eq), *symbols)
        if without_radicals:
            unrad_changed.append(eq)
            eq_unrad, cov = without_radicals
            if not cov:
                eq = eq_unrad
        if getattr(eq, 'is_Relational', False) is False and hasattr(eq, 'as_numer_denom'):
            eq = eq.as_numer_denom()[0]
            poly = eq.as_poly(*symbols, extension=True)
        elif simplify(eq).is_number:
            continue
        if poly is not None:
            polys.append(poly)
            polys_expr.append(poly.as_expr())
        else:
            nonpolys.append(eq)
    return polys, polys_expr, nonpolys, denominators, unrad_changed


def _is_polynomial_system(system, symbols):
    """Whether the reference ``nonlinsolve`` sees only polynomial equations,
    in which case it never reaches ``substitution``."""
    for eq in system:
        if getattr(eq, 'is_Equality', False):
            eq = eq.lhs - eq.rhs
        numerator = eq.as_numer_denom()[0]
        if not numerator.is_polynomial(*symbols):
            return False
    return True


def _handle_poly(polys, symbols):
    from .. import groebner, nsimplify, Poly, nfloat, solve_poly_system
    no_information = [{}]
    inexact = any(not p.domain.is_Exact for p in polys)
    if inexact:
        polys = [Poly(nsimplify(p.as_expr(), rational=True), *symbols) for p in polys]
    basis = groebner([p.as_expr() for p in polys], *symbols, order='grevlex')
    if 1 in list(basis):
        return [], []
    if basis.is_zero_dimensional:
        basis = basis.fglm('lex')
        basis = list(basis)
        if inexact:
            basis = [nfloat(p) for p in basis]
        try:
            result = solve_poly_system(basis, *symbols, strict=True)
        except NotImplementedError:
            return no_information, list(basis)
        if result is None:
            return [], []
        return [dict(zip(symbols, res)) for res in result], []
    poly_eqs = list(groebner([p.as_expr() for p in polys], *symbols, order='lex'))
    if inexact:
        poly_eqs = [nfloat(p) for p in poly_eqs]
    return no_information, poly_eqs


def nonlinsolve_substitution(system, symbols):
    """The reference ``nonlinsolve`` branch for systems with at least one
    non-polynomial equation: the polynomial part is solved by Groebner
    elimination, the rest by ``substitution`` on each polynomial
    solution."""
    from .. import FiniteSet, Set
    from ..solvers.solvers import checksol
    polys, polys_expr, nonpolys, denominators, unrad_changed = _separate_poly_nonpoly(system, symbols)
    poly_eqs = []
    poly_sol = [{}]
    if polys:
        poly_sol, poly_eqs = _handle_poly(polys, symbols)
        if poly_sol and poly_sol[0]:
            poly_syms = set().union(*(eq.free_symbols for eq in polys))
            unrad_syms = set().union(*(eq.free_symbols for eq in unrad_changed))
            if unrad_syms == poly_syms and unrad_changed:
                poly_sol = [sol for sol in poly_sol if checksol(unrad_changed, sol)]
    remaining = poly_eqs + nonpolys

    def to_tuple(sol):
        return tuple(sol[s] for s in symbols)
    if not remaining:
        return FiniteSet(*map(to_tuple, poly_sol))
    subs_res = substitution(remaining, symbols, result=poly_sol, exclude=denominators)
    if not isinstance(subs_res, FiniteSet):
        return subs_res
    if unrad_changed:
        result = [dict(zip(symbols, sol)) for sol in subs_res.args]
        correct_sols = [sol for sol in result if any(isinstance(v, Set) for v in sol)
                        or checksol(unrad_changed, sol) is not False]
        return FiniteSet(*map(to_tuple, correct_sols))
    return subs_res
