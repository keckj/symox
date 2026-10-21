"""Helpers for the ODE solvers adapted from SymPy 1.14; see LICENSE-SYMPY."""
from collections import defaultdict

from symox import Add, Mul, S, Pow, Derivative, Equality
from symox.core.function import AppliedUndef, expand_mul, expand_multinomial
from symox.polys.solvers import PolyNonlinearError


def _mexpand(expr, recursive=False):
    if expr is None:
        return None
    was = None
    while was != expr:
        was, expr = expr, expand_mul(expand_multinomial(expr))
        if not recursive:
            break
    return expr


def _lin_eq2dict(a, symset):
    """Return the symbol-free part of a and its coefficient dictionary."""
    if a in symset:
        return S.Zero, {a: S.One}
    elif a.is_Add:
        terms_list = defaultdict(list)
        coeff_list = []
        for ai in a.args:
            ci, ti = _lin_eq2dict(ai, symset)
            coeff_list.append(ci)
            for mij, cij in ti.items():
                terms_list[mij].append(cij)
        coeff = Add(*coeff_list)
        terms = {sym: Add(*coeffs) for sym, coeffs in terms_list.items()}
        return coeff, terms
    elif a.is_Mul:
        terms = terms_coeff = None
        coeff_list = []
        for ai in a.args:
            ci, ti = _lin_eq2dict(ai, symset)
            if not ti:
                coeff_list.append(ci)
            elif terms is None:
                terms = ti
                terms_coeff = ci
            else:
                raise PolyNonlinearError('nonlinear cross-term: %s' % a)
        coeff = Mul(*coeff_list)
        if terms is None:
            return coeff, {}
        terms = {sym: coeff * c for sym, c in terms.items()}
        return coeff * terms_coeff, terms
    elif not a.has_xfree(symset):
        return a, {}
    else:
        raise PolyNonlinearError('nonlinear term: %s' % a)


def _preprocess(expr, func=None, hint='_Integral'):
    """See SymPy 1.14 sympy.solvers.deutils."""
    if isinstance(expr, Pow):
        # if f(x)**p=0 then f(x)=0 (p>0)
        if (expr.exp).is_positive:
            expr = expr.base
    derivs = expr.atoms(Derivative)
    if not func:
        funcs = set().union(*[d.atoms(AppliedUndef) for d in derivs])
        if len(funcs) != 1:
            raise ValueError('The function cannot be '
                'automatically detected for %s.' % expr)
        func = funcs.pop()
    fvars = set(func.args)
    if hint is None:
        return expr, func
    reps = [(d, d.doit()) for d in derivs if not hint.endswith('_Integral') or
            d.has(func) or set(d.variables) & fvars]
    eq = expr.subs(reps)
    return eq, func


def _desolve(eq, func=None, hint="default", ics=None, simplify=True, *, prep=True, **kwargs):
    """See SymPy 1.14 sympy.solvers.deutils."""
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs

    # preprocess the equation and find func if not given
    if prep or func is None:
        eq, func = _preprocess(eq, func)
        prep = False

    # type is an argument passed by the solve functions in ode and pde.py
    # that identifies whether the function caller is an ordinary
    # or partial differential equation. Accordingly corresponding
    # changes are made in the function.
    type = kwargs.get('type', None)
    xi = kwargs.get('xi')
    eta = kwargs.get('eta')
    x0 = kwargs.get('x0', 0)
    terms = kwargs.get('n')

    if type == 'ode':
        from .ode import classify_ode, allhints
        classifier = classify_ode
        string = 'ODE '
        dummy = ''

    elif type == 'pde':
        from symox.solvers.pde import classify_pde, allhints
        classifier = classify_pde
        string = 'PDE '
        dummy = 'p'

    # Magic that should only be used internally.  Prevents classify_ode from
    # being called more than it needs to be by passing its results through
    # recursive calls.
    if kwargs.get('classify', True):
        hints = classifier(eq, func, dict=True, ics=ics, xi=xi, eta=eta,
        n=terms, x0=x0, hint=hint, prep=prep)

    else:
        # Here is what all this means:
        #
        # hint:    The hint method given to _desolve() by the user.
        # hints:   The dictionary of hints that match the DE, along with other
        #          information (including the internal pass-through magic).
        # default: The default hint to return, the first hint from allhints
        #          that matches the hint; obtained from classify_ode().
        # match:   Dictionary containing the match dictionary for each hint
        #          (the parts of the DE for solving).  When going through the
        #          hints in "all", this holds the match string for the current
        #          hint.
        # order:   The order of the DE, as determined by ode_order().
        hints = kwargs.get('hint',
                           {'default': hint,
                            hint: kwargs['match'],
                            'order': kwargs['order']})
    if not hints['default']:
        # classify_ode will set hints['default'] to None if no hints match.
        if hint not in allhints and hint != 'default':
            raise ValueError("Hint not recognized: " + hint)
        elif hint not in hints['ordered_hints'] and hint != 'default':
            raise ValueError(string + str(eq) + " does not match hint " + hint)
        # If dsolve can't solve the purely algebraic equation then dsolve will raise
        # ValueError
        elif hints['order'] == 0:
            raise ValueError(
                str(eq) + " is not a solvable differential equation in " + str(func))
        else:
            raise NotImplementedError(dummy + "solve" + ": Cannot solve " + str(eq))
    if hint == 'default':
        return _desolve(eq, func, ics=ics, hint=hints['default'], simplify=simplify,
                      prep=prep, x0=x0, classify=False, order=hints['order'],
                      match=hints[hints['default']], xi=xi, eta=eta, n=terms, type=type)
    elif hint in ('all', 'all_Integral', 'best'):
        retdict = {}
        gethints = set(hints) - {'order', 'default', 'ordered_hints'}
        if hint == 'all_Integral':
            for i in hints:
                if i.endswith('_Integral'):
                    gethints.remove(i.removesuffix('_Integral'))
            # special cases
            for k in ["1st_homogeneous_coeff_best", "1st_power_series",
                "lie_group", "2nd_power_series_ordinary", "2nd_power_series_regular"]:
                if k in gethints:
                    gethints.remove(k)
        for i in gethints:
            sol = _desolve(eq, func, ics=ics, hint=i, x0=x0, simplify=simplify, prep=prep,
                classify=False, n=terms, order=hints['order'], match=hints[i], type=type)
            retdict[i] = sol
        retdict['all'] = True
        retdict['eq'] = eq
        return retdict
    elif hint not in allhints:  # and hint not in ('default', 'ordered_hints'):
        raise ValueError("Hint not recognized: " + hint)
    elif hint not in hints:
        raise ValueError(string + str(eq) + " does not match hint " + hint)
    else:
        # Key added to identify the hint needed to solve the equation
        hints['hint'] = hint
    hints.update({'func': func, 'eq': eq})
    return hints
