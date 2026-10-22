"""Reference handling of noncommutative input in the rational simplifiers.

The native rational algorithms (polynomial conversion, cancellation,
factorisation, candidate search) assume commutative generators.  When an
expression contains a noncommutative Symbol or undefined function, the
public entry points below follow SymPy 1.14 instead: ``_mask_nc``, the
``gcd_terms`` mask, ``factor_nc``, ``together`` and the PolynomialError
branches of ``cancel``, ``factor``, ``apart`` and ``ratsimp`` are ported
from the reference, and commutative subexpressions still go through the
native code.

Stages of ``simplify`` without a noncommutative reference algorithm
(trigsimp, logcombine, hyperexpand, ...) run on a *linear mask*: the
noncommutative part of every term is replaced by a commutative Dummy, the
native function runs on the masked expression, and the result is restored
only when it is still linear in the Dummies with no Dummy in a
denominator.  Otherwise the stage is skipped, so no factor is ever
commuted.
"""
from . import _core

NC_POLY_MESSAGE = 'non-commutative expressions are not supported'

_native = {}


def is_nc(e):
    """True when ``e`` contains a noncommutative Symbol or function."""
    if isinstance(e, _core.Expr):
        return _core._has_noncommutative_atoms(e)
    if getattr(e, 'is_Relational', False):
        return is_nc(e.lhs) or is_nc(e.rhs)
    return False


def _any_nc(seq):
    return any(is_nc(e) for e in seq)


def _polynomial_error():
    from . import PolynomialError
    return PolynomialError(NC_POLY_MESSAGE)


def _sympify(e):
    from . import sympify
    return sympify(e)


# --- generic helpers ---------------------------------------------------------

def _bottom_up(rv, F, atoms=False, nonbasic=False):
    args = getattr(rv, 'args', None)
    if args is not None:
        if args:
            args = tuple([_bottom_up(a, F, atoms, nonbasic) for a in args])
            if args != rv.args:
                rv = rv.func(*args)
            rv = F(rv)
        elif atoms:
            rv = F(rv)
    else:
        if nonbasic:
            try:
                rv = F(rv)
            except TypeError:
                pass
    return rv


def _linear_mask(expr):
    """Mask the noncommutative part of each term with a commutative Dummy."""
    from . import Add, Mul, Dummy
    reps, back, terms = {}, {}, []
    for t in Add.make_args(expr):
        if t.is_commutative:
            terms.append(t)
            continue
        c, nc = t.args_cnc()
        part = Mul(*nc)
        d = reps.get(part)
        if d is None:
            d = reps[part] = Dummy('nc')
            back[d] = part
        terms.append(Mul(*c)*d)
    return Add(*terms), back


def _linear_in(expr, dummies):
    from . import Add, Mul, expand_mul
    n, d = expr.as_numer_denom()
    if d.free_symbols & dummies:
        return False
    for t in Add.make_args(expand_mul(n)):
        k = 0
        for f in Mul.make_args(t):
            b, e = f.as_base_exp()
            if b in dummies:
                if e != 1:
                    return False
                k += 1
            elif f.free_symbols & dummies:
                return False
        if k > 1:
            return False
    return True


def masked(fn, expr, *args, **kwargs):
    """Apply a commutative transformation to ``expr`` without commuting."""
    if not is_nc(expr):
        return fn(expr, *args, **kwargs)
    m, back = _linear_mask(expr)
    if is_nc(m):
        return expr
    result = fn(m, *args, **kwargs)
    if not isinstance(result, _core.Expr) or not _linear_in(result, set(back)):
        return expr
    return result.xreplace(back)


# --- core.exprtools ----------------------------------------------------------

def _mask_nc(eq, name=None):
    """SymPy's ``_mask_nc``."""
    from . import Dummy as _Dummy, preorder_traversal, default_sort_key
    name = name or 'mask'

    def numbered_names():
        i = 0
        while True:
            yield name + str(i)
            i += 1

    names = numbered_names()

    def Dummy(*args, **kwargs):
        return _Dummy(next(names), *args, **kwargs)

    expr = eq
    if expr.is_commutative:
        return eq, {}, []

    rep = []
    nc_obj = set()
    nc_syms = set()
    pot = preorder_traversal(expr, keys=default_sort_key)
    for a in pot:
        if any(a == r[0] for r in rep):
            pot.skip()
        elif not a.is_commutative:
            if a.is_symbol:
                nc_syms.add(a)
                pot.skip()
            elif not (a.is_Add or a.is_Mul or a.is_Pow):
                nc_obj.add(a)
                pot.skip()

    if len(nc_obj) == 1 and not nc_syms:
        rep.append((nc_obj.pop(), Dummy()))
    elif len(nc_syms) == 1 and not nc_obj:
        rep.append((nc_syms.pop(), Dummy()))

    nc_obj = sorted(nc_obj, key=default_sort_key)
    for n in nc_obj:
        nc = Dummy(commutative=False)
        rep.append((n, nc))
        nc_syms.add(nc)
    expr = expr.subs(rep)

    nc_syms = list(nc_syms)
    nc_syms.sort(key=default_sort_key)
    return expr, {v: k for k, v in rep}, nc_syms


def gcd_terms(terms, isprimitive=False, clear=True, fraction=True):
    """SymPy's ``gcd_terms`` for noncommutative input."""
    from . import Add, Mul, Dummy, Dict
    from .core.mul import _keep_coeff
    from ._simptools import gcd_terms as public_gcd_terms

    def mask(terms):
        args = [(a, []) if a.is_commutative else a.args_cnc() for a in terms]
        reps = []
        for i, (c, nc) in enumerate(args):
            if nc:
                nc = Mul(*nc)
                d = Dummy()
                reps.append((d, nc))
                c.append(d)
                args[i] = Mul(*c)
            else:
                args[i] = c
        return args, dict(reps)

    if isinstance(terms, _core.Expr) and terms.is_Add or isinstance(terms, (list, tuple, set)):
        if isinstance(terms, _core.Expr):
            terms = list(terms.args)
        else:
            terms = [_sympify(t) for t in terms]
        terms, reps = mask(terms)
        rv = _core._st_gcd_terms_list(terms, bool(isprimitive), bool(clear), bool(fraction))
        return rv.xreplace(reps)

    if not hasattr(terms, 'args') or not terms.args:
        return terms

    if terms.is_Mul:
        c, args = terms.as_coeff_mul()
        return _keep_coeff(c, Mul(*[public_gcd_terms(i, isprimitive, clear, fraction)
            for i in args]), clear=clear)

    def handle(a):
        if not isinstance(a, _core.Expr):
            if hasattr(a, 'args'):
                if not a.args:
                    return a
                return a.func(*[handle(i) for i in a.args])
            return type(a)([handle(i) for i in a])
        return public_gcd_terms(a, isprimitive, clear, fraction)

    if isinstance(terms, Dict):
        return Dict(*[(k, handle(v)) for k, v in terms.args])
    return terms.func(*[handle(i) for i in terms.args])


def factor_terms(expr, radical=False, clear=False, fraction=False, sign=True):
    """SymPy's ``factor_terms`` for noncommutative input."""
    from . import Add, Dummy, Sum, Integral
    from .core.mul import _keep_coeff, Mul

    def do(expr):
        if not is_nc(expr):
            return _core._st_factor_terms(expr, bool(radical), bool(clear), bool(fraction), bool(sign))
        if not expr.args:
            return expr
        if expr.is_Pow or expr.is_Function or not hasattr(expr, 'args_cnc'):
            args = expr.args
            newargs = tuple([do(i) for i in args])
            if newargs == args:
                return expr
            return expr.func(*newargs)
        if isinstance(expr, (Sum, Integral)):
            return _core._st_factor_terms(expr, bool(radical), bool(clear), bool(fraction), bool(sign))

        cont, p = expr.as_content_primitive(radical=radical, clear=clear)
        if p.is_Add:
            list_args = [do(a) for a in Add.make_args(p)]
            if not any(a.as_coeff_Mul()[0].extract_multiplicatively(-1) is None
                       for a in list_args):
                cont = -cont
                list_args = [-a for a in list_args]
            special = {}
            for i, a in enumerate(list_args):
                b, e = a.as_base_exp()
                if e.is_Mul and e != Mul(*e.args):
                    list_args[i] = Dummy()
                    special[list_args[i]] = a
            p = Add._from_args(list_args)
            p = gcd_terms(p, isprimitive=True, clear=clear, fraction=fraction).xreplace(special)
        elif p.args:
            p = p.func(*[do(a) for a in p.args])
        return _keep_coeff(cont, p, clear=clear, sign=sign)

    return do(_sympify(expr))


def factor_nc(expr):
    """SymPy's ``factor_nc``."""
    from . import Add, Mul, S, Dummy, gcd, expand_power_exp, default_sort_key, powsimp
    from .core.mul import _keep_coeff
    from .utilities.iterables import common_prefix, common_suffix, variations
    from . import factor

    expr = _sympify(expr)
    if not isinstance(expr, _core.Expr) or not expr.args:
        return expr
    if not expr.is_Add:
        return expr.func(*[factor_nc(a) for a in expr.args])
    expr = expr.func(*[expand_power_exp(i) for i in expr.args])

    expr, rep, nc_symbols = _mask_nc(expr)

    if rep:
        return factor(expr).subs(rep)
    args = [a.args_cnc() for a in Add.make_args(expr)]
    c = g = l = r = S.One
    hit = False
    for i, a in enumerate(args):
        if i == 0:
            c = Mul._from_args(a[0])
        elif a[0]:
            c = gcd(c, Mul._from_args(a[0]))
        else:
            c = S.One
    if c is not S.One:
        hit = True
        c, g = c.as_coeff_Mul()
        if g is not S.One:
            for i, (cc, _) in enumerate(args):
                cc = list(Mul.make_args(Mul._from_args(list(cc))/g))
                args[i][0] = cc
        for i, (cc, _) in enumerate(args):
            if cc:
                cc[0] = cc[0]/c
            else:
                cc = [1/c]
            args[i][0] = cc
    for i, a in enumerate(args):
        if i == 0:
            n = a[1][:]
        else:
            n = common_prefix(n, a[1])
        if not n:
            if not args[0][1]:
                break
            b, e = args[0][1][0].as_base_exp()
            ok = False
            if e.is_Integer:
                for t in args:
                    if not t[1]:
                        break
                    bt, et = t[1][0].as_base_exp()
                    if et.is_Integer and bt == b:
                        e = min(e, et)
                    else:
                        break
                else:
                    ok = hit = True
                    l = b**e
                    il = b**-e
                    for _ in args:
                        _[1][0] = il*_[1][0]
                    break
            if not ok:
                break
    else:
        hit = True
        lenn = len(n)
        l = Mul(*n)
        for _ in args:
            _[1] = _[1][lenn:]
    for i, a in enumerate(args):
        if i == 0:
            n = a[1][:]
        else:
            n = common_suffix(n, a[1])
        if not n:
            if not args[0][1]:
                break
            b, e = args[0][1][-1].as_base_exp()
            ok = False
            if e.is_Integer:
                for t in args:
                    if not t[1]:
                        break
                    bt, et = t[1][-1].as_base_exp()
                    if et.is_Integer and bt == b:
                        e = min(e, et)
                    else:
                        break
                else:
                    ok = hit = True
                    r = b**e
                    il = b**-e
                    for _ in args:
                        _[1][-1] = _[1][-1]*il
                    break
            if not ok:
                break
    else:
        hit = True
        lenn = len(n)
        r = Mul(*n)
        for _ in args:
            _[1] = _[1][:len(_[1]) - lenn]
    if hit:
        mid = Add(*[Mul(*cc)*Mul(*nc) for cc, nc in args])
    else:
        mid = expr

    rep1 = [(n, Dummy()) for n in sorted(nc_symbols, key=default_sort_key)]
    unrep1 = [(v, k) for k, v in rep1]
    unrep1.reverse()
    new_mid, r2, _ = _mask_nc(mid.subs(rep1))
    new_mid = powsimp(factor(new_mid))

    new_mid = new_mid.subs(r2).subs(unrep1)

    if new_mid.is_Pow:
        return _keep_coeff(c, g*l*new_mid*r)

    if new_mid.is_Mul:
        def _pemexpand(expr):
            return expr.expand(deep=True, mul=True, power_exp=True,
                power_base=False, basic=False, multinomial=True, log=False)
        cfac = []
        ncfac = []
        for f in new_mid.args:
            if f.is_commutative:
                cfac.append(f)
            else:
                b, e = f.as_base_exp()
                if e.is_Integer:
                    ncfac.extend([b]*e)
                else:
                    ncfac.append(f)
        pre_mid = g*Mul(*cfac)*l
        target = _pemexpand(expr/c)
        for s in variations(ncfac, len(ncfac)):
            ok = pre_mid*Mul(*s)*r
            if _pemexpand(ok) == target:
                return _keep_coeff(c, ok)

    return _keep_coeff(c, g*l*mid*r)


# --- polys -------------------------------------------------------------------

def together(expr, deep=False, fraction=True):
    """SymPy's ``together`` for noncommutative input."""
    from . import Add

    def _together(expr):
        if not is_nc(expr):
            return _native['together'](expr, deep=deep, fraction=fraction) if hasattr(expr, 'args') else expr
        if not expr.args or (expr.is_Function and not deep):
            return expr
        elif expr.is_Add:
            return gcd_terms(list(map(_together, Add.make_args(expr))), fraction=fraction)
        elif expr.is_Pow:
            base = _together(expr.base)
            exp = _together(expr.exp) if deep else expr.exp
            return expr.func(base, exp)
        return expr.func(*[_together(arg) for arg in expr.args])

    return _together(_sympify(expr))


def cancel(f, *gens, _signsimp=True, **args):
    """SymPy's ``cancel`` PolynomialError branch for noncommutative input."""
    from . import signsimp, Piecewise, preorder_traversal
    from .utilities.iterables import sift
    from . import cancel as public_cancel

    f = _sympify(f)
    if _signsimp:
        f = signsimp(f)
    if f.is_Number or f.is_Relational or not isinstance(f, _core.Expr):
        return f
    f = factor_terms(f, radical=True)
    if not is_nc(f):
        return _native['cancel'](f, *gens, **args)
    if f.is_Add or f.is_Mul:
        c, nc = sift(f.args, lambda x: x.is_commutative is True and not x.has(Piecewise), binary=True)
        nc = [public_cancel(i) for i in nc]
        return f.func(public_cancel(f.func(*c)), *nc)
    reps = []
    pot = preorder_traversal(f)
    next(pot)
    for e in pot:
        if not isinstance(e, _core.Expr) or e.is_Boolean:
            continue
        try:
            reps.append((e, public_cancel(e)))
            pot.skip()
        except NotImplementedError:
            pass
    return f.xreplace(dict(reps))


def factor(f, *gens, deep=False, **args):
    """SymPy's ``factor``: noncommutative input goes to ``factor_nc``."""
    from . import Mul, Add
    from . import factor as public_factor
    f = _sympify(f)
    if deep:
        def _try_factor(expr):
            fac = public_factor(expr, *gens, **args)
            if fac.is_Mul or fac.is_Pow:
                return fac
            return expr

        f = _bottom_up(f, _try_factor)
        partials = {}
        for p in f.atoms(Mul, Add):
            fac = public_factor(p, *gens, **args)
            if (fac.is_Mul or fac.is_Pow) and fac != p:
                partials[p] = fac
        return f.xreplace(partials)
    return factor_nc(f)


def apart(f, x=None, full=False, **options):
    """SymPy's ``apart`` PolynomialError branch for noncommutative input."""
    from . import preorder_traversal
    from . import apart as public_apart
    f = _sympify(f)
    if not f.args:
        return f
    if f.is_Mul:
        c, nc = f.args_cnc(split_1=False)
        nc = f.func(*nc)
        if c:
            c = public_apart(f.func._from_args(c), x=x, full=full, **options)
            return c*nc
        return nc
    if f.is_Add:
        c = []
        nc = []
        for i in f.args:
            if i.is_commutative:
                c.append(i)
            else:
                try:
                    nc.append(public_apart(i, x=x, full=full, **options))
                except NotImplementedError:
                    nc.append(i)
        return public_apart(f.func(*c), x=x, full=full, **options) + f.func(*nc)
    reps = []
    pot = preorder_traversal(f)
    next(pot)
    for e in pot:
        try:
            reps.append((e, public_apart(e, x=x, full=full, **options)))
            pot.skip()
        except NotImplementedError:
            pass
    return f.xreplace(dict(reps))


def ratsimp(expr):
    """SymPy's ``ratsimp``: ``reduced`` rejects noncommutative input."""
    from . import cancel as public_cancel
    f, g = public_cancel(expr).as_numer_denom()
    if is_nc(f) or is_nc(g):
        raise _polynomial_error()
    return _core._st_ratsimp(f/g)


_GAMMA_LIKE = ('gamma', 'factorial', 'factorial2', 'binomial', 'RisingFactorial', 'FallingFactorial',
               'subfactorial', 'lowergamma', 'uppergamma', 'beta', 'loggamma', 'polygamma', 'digamma',
               'trigamma', 'multigamma')


def _has_gamma_like(expr):
    from . import preorder_traversal
    return any(type(a).__name__ in _GAMMA_LIKE or getattr(getattr(a, 'func', None), '__name__', '') in _GAMMA_LIKE
               for a in preorder_traversal(expr))


def combsimp(expr, native):
    """``combsimp``/``gammasimp``: without gamma-like factors the reference
    reduces to ``factor`` (``factor_nc``) of the expression."""
    from . import factor as public_factor
    if not _has_gamma_like(expr):
        return public_factor(expr)
    return masked(native, expr)


def collect_check(expr, syms):
    """Raise like SymPy when a collected symbol is noncommutative."""
    if not isinstance(syms, (list, tuple, set, frozenset)):
        syms = [syms]
    if any(getattr(_sympify(s), 'is_commutative', True) is False for s in syms):
        from . import Add, Mul
        for product in Add.make_args(expr):
            for s in syms:
                if _sympify(s).is_commutative is False and product.has(s):
                    raise AttributeError("Can not collect noncommutative symbol")


# --- simplify ----------------------------------------------------------------

def nc_simplify(expr, deep=True):
    """SymPy's ``nc_simplify`` for scalar expressions."""
    from . import Add, Mul, Pow, Symbol

    def _overlaps(args):
        m = [[([1, 0] if a == args[0] else [0]) for a in args[1:]]]
        for i in range(1, len(args)):
            overlaps = []
            j = 0
            for j in range(len(args) - i - 1):
                overlap = []
                for v in m[i-1][j+1]:
                    if j + i + 1 + v < len(args) and args[i] == args[j+i+1+v]:
                        overlap.append(v + 1)
                overlap += [0]
                overlaps.append(overlap)
            m.append(overlaps)
        return m

    def _reduce_inverses(_args):
        inv_tot = 0
        inverses = []
        args = []
        for arg in _args:
            if arg.is_Pow and arg.args[1].is_extended_negative:
                inverses = [arg**-1] + inverses
                inv_tot += 1
            else:
                if len(inverses) == 1:
                    args.append(inverses[0]**-1)
                elif len(inverses) > 1:
                    args.append(Pow(Mul(*inverses), -1))
                    inv_tot -= len(inverses) - 1
                inverses = []
                args.append(arg)
        if inverses:
            args.append(Pow(Mul(*inverses), -1))
            inv_tot -= len(inverses) - 1
        return inv_tot, tuple(args)

    def get_score(s):
        if s.is_Pow:
            return get_score(s.args[0])
        elif s.is_Add or s.is_Mul:
            return sum(get_score(a) for a in s.args)
        return 1

    def compare(s, alt_s):
        if s != alt_s and get_score(alt_s) < get_score(s):
            return alt_s
        return s

    if not isinstance(expr, _core.Expr) or not (expr.is_Add or expr.is_Mul or expr.is_Pow) or expr.is_commutative:
        return expr
    args = expr.args[:]
    if expr.is_Pow:
        if deep:
            return Pow(nc_simplify(args[0]), args[1]).doit()
        else:
            return expr
    elif expr.is_Add:
        return Add(*[nc_simplify(a, deep=deep) for a in args]).doit()
    else:
        c_args, args = expr.args_cnc()
        com_coeff = Mul(*c_args)
        if com_coeff != 1:
            return com_coeff*nc_simplify(expr/com_coeff, deep=deep)

    inv_tot, args = _reduce_inverses(args)
    invert = False
    if inv_tot > len(args)/2:
        invert = True
        args = [a**-1 for a in args[::-1]]

    if deep:
        args = tuple(nc_simplify(a) for a in args)

    m = _overlaps(args)
    simps = {}
    post = 1
    pre = 1
    max_simp_coeff = 0
    simp = None
    for i in range(1, len(args)):
        simp_coeff = 0
        l = 0
        p = 0
        if i < len(args) - 1:
            rep = m[i][0]
        start = i
        end = i+1
        if i == len(args)-1 or rep == [0]:
            if args[i].is_Pow and not isinstance(args[i].args[0], Symbol):
                subterm = args[i].args[0].args
                l = len(subterm)
                if args[i-l:i] == subterm:
                    p += 1
                    start -= l
                if args[i+1:i+1+l] == subterm:
                    p += 1
                    end += l
            if p:
                p += args[i].args[1]
            else:
                continue
        else:
            l = rep[0]
            start -= l - 1
            subterm = args[start:end]
            p = 2
            end += l

        if subterm in simps and simps[subterm] >= start:
            continue

        while end < len(args):
            if l in m[end-1][0]:
                p += 1
                end += l
            elif args[end].is_Pow and args[end].args[0].args == subterm:
                p += args[end].args[1]
                end += 1
            else:
                break

        pre_exp = 0
        pre_arg = 1
        if start - l >= 0 and args[start-l+1:start] == subterm[1:]:
            if subterm[0].is_Pow:
                pre_arg = subterm[0].args[0]
                exp = subterm[0].args[1]
            else:
                pre_arg = subterm[0]
                exp = 1
            if args[start-l].is_Pow and args[start-l].args[0] == pre_arg:
                pre_exp = args[start-l].args[1] - exp
                start -= l
                p += 1
            elif args[start-l] == pre_arg:
                pre_exp = 1 - exp
                start -= l
                p += 1

        post_exp = 0
        post_arg = 1
        if end + l - 1 < len(args) and args[end:end+l-1] == subterm[:-1]:
            if subterm[-1].is_Pow:
                post_arg = subterm[-1].args[0]
                exp = subterm[-1].args[1]
            else:
                post_arg = subterm[-1]
                exp = 1
            if args[end+l-1].is_Pow and args[end+l-1].args[0] == post_arg:
                post_exp = args[end+l-1].args[1] - exp
                end += l
                p += 1
            elif args[end+l-1] == post_arg:
                post_exp = 1 - exp
                end += l
                p += 1

        if post_exp and exp % 2 == 0 and start > 0:
            exp = exp/2
            _pre_exp = 1
            _post_exp = 1
            if args[start-1].is_Pow and args[start-1].args[0] == post_arg:
                _post_exp = post_exp + exp
                _pre_exp = args[start-1].args[1] - exp
            elif args[start-1] == post_arg:
                _post_exp = post_exp + exp
                _pre_exp = 1 - exp
            if _pre_exp == 0 or _post_exp == 0:
                if not pre_exp:
                    start -= 1
                post_exp = _post_exp
                pre_exp = _pre_exp
                pre_arg = post_arg
                subterm = (post_arg**exp,) + subterm[:-1] + (post_arg**exp,)

        simp_coeff += end-start

        if post_exp:
            simp_coeff -= 1
        if pre_exp:
            simp_coeff -= 1

        simps[subterm] = end

        if simp_coeff > max_simp_coeff:
            max_simp_coeff = simp_coeff
            simp = (start, Mul(*subterm), p, end, l)
            pre = pre_arg**pre_exp
            post = post_arg**post_exp

    if simp:
        subterm = Pow(nc_simplify(simp[1], deep=deep), simp[2])
        pre = nc_simplify(Mul(*args[:simp[0]])*pre, deep=deep)
        post = post*nc_simplify(Mul(*args[simp[3]:]), deep=deep)
        simp = pre*subterm*post
        if pre != 1 or post != 1:
            simp = nc_simplify(simp, deep=False)
    else:
        simp = Mul(*args)

    if invert:
        simp = Pow(simp, -1)

    f_expr = factor_nc(expr)
    if f_expr != expr:
        alt_simp = nc_simplify(f_expr, deep=deep)
        simp = compare(simp, alt_simp)
    return simp


def _has_class(expr, *modules_names):
    from importlib import import_module
    classes = []
    for mod, name in modules_names:
        try:
            classes.append(getattr(import_module(mod), name))
        except (ImportError, AttributeError):
            pass
    return bool(classes) and expr.has(*classes)


_TRIG_HYPER = (('symox.functions.elementary.trigonometric', 'TrigonometricFunction'),
               ('symox.functions.elementary.hyperbolic', 'HyperbolicFunction'))
_COMBINATORIAL = (('symox.functions.combinatorial.factorials', 'CombinatorialFunction'),
                  ('symox', 'gamma'))


def simplify(expr, ratio=1.7, measure=None, rational=False, inverse=False, doit=True, **kwargs):
    """SymPy's ``simplify`` pipeline for noncommutative input.

    Stages with a noncommutative reference algorithm (``nc_simplify``,
    ``cancel``, ``together``, ``factor_terms``, ``powsimp``) follow it;
    the commutative-only stages run through the linear mask."""
    from . import (S, Mul, Add, Float, Piecewise, Abs, sign, exp, log, nsimplify, nfloat, signsimp,
                   powsimp, expand_mul, expand_power_exp, hyperexpand, trigsimp, expand_log, logcombine,
                   exptrigsimp, radsimp, fraction, piecewise_fold, count_ops, Function, KroneckerDelta,
                   kroneckersimp, besselsimp, Sum, Integral, Product, sympify)
    from . import simplify as public_simplify, combsimp as public_combsimp
    from . import cancel as public_cancel, together as public_together, factor_terms as public_factor_terms
    from importlib import import_module
    ssimp = import_module('symox.simplify.simplify')
    from . import expand_multinomial

    def _mexpand(e):
        return expand_mul(expand_multinomial(e))
    if measure is None:
        measure = count_ops

    def shorter(*choices):
        if not any(s != choices[0] for s in choices[1:]):
            return choices[0]
        return min(choices, key=measure)

    def done(e):
        rv = e.doit() if doit else e
        return rv

    options = dict(ratio=ratio, measure=measure, rational=rational, inverse=inverse, doit=doit)
    expr = sympify(expr, rational=rational)
    original_expr = expr = signsimp(expr)
    if not expr.args:
        return expr
    if inverse and expr.has(Function):
        expr = ssimp.inversecombine(expr)
        if not expr.args:
            return expr

    def handled(x):
        return x.is_Add or x.is_Mul or x.is_Pow or x.func == exp

    def deep(e):
        # Expr.replace(query, value, simultaneous=False): bottom-up.
        if not getattr(e, 'args', None):
            return e
        args = tuple(deep(a) for a in e.args)
        if args != e.args:
            e = e.func(*args)
        if isinstance(e, _core.Expr) and e.args and not handled(e):
            e = e.func(*[public_simplify(i, **options) for i in e.args])
        return e

    expr = deep(expr)
    if not isinstance(expr, _core.Expr) or not handled(expr):
        return done(expr)

    if not expr.is_commutative:
        expr = nc_simplify(expr)

    floats = False
    if rational is not False and expr.has(Float):
        floats = True
        expr = nsimplify(expr, rational=True)

    expr = _bottom_up(expr, lambda w: getattr(w, 'normal', lambda: w)())
    expr = Mul(*powsimp(expr).as_content_primitive())
    _e = public_cancel(expr)
    expr1 = shorter(_e, public_cancel(_mexpand(_e)))
    expr2 = shorter(public_together(expr, deep=True), public_together(expr1, deep=True))

    if ratio is S.Infinity:
        expr = expr2
    else:
        expr = shorter(expr2, expr1, expr)
    if not isinstance(expr, _core.Expr):
        return expr

    expr = public_factor_terms(expr, sign=False)

    if expr.has(sign):
        expr = expr.rewrite(Abs)

    if expr.has(Piecewise):
        expr = done(piecewise_fold(expr))
        if expr.has(Piecewise):
            return shorter(expr, public_factor_terms(expr))

    expr = masked(hyperexpand, expr)
    if expr.has(KroneckerDelta):
        expr = masked(kroneckersimp, expr)
    if _has_class(expr, ('symox.functions.special.bessel', 'BesselBase')):
        expr = masked(besselsimp, expr)
    if _has_class(expr, *_TRIG_HYPER):
        expr = masked(trigsimp, expr, deep=True)
    if expr.has(log):
        expr = shorter(masked(expand_log, expr, deep=True), masked(logcombine, expr))
    if _has_class(expr, *_COMBINATORIAL):
        expr = public_combsimp(expr)
    if expr.has(Sum) or expr.has(Product):
        expr = masked(lambda e: public_simplify(e, **options), expr)
    if expr.has(Integral):
        expr = expr.xreplace({i: public_factor_terms(i) for i in expr.atoms(Integral)})

    short = shorter(powsimp(expr, combine='exp', deep=True), powsimp(expr), expr)
    short = shorter(short, public_cancel(short))
    short = shorter(short, public_factor_terms(short), expand_power_exp(expand_mul(short)))
    if _has_class(short, *_TRIG_HYPER) or short.has(exp):
        short = masked(exptrigsimp, short)

    def hollow(e):
        if not getattr(e, 'args', None):
            return e
        args = tuple(hollow(a) for a in e.args)
        if args != e.args:
            e = e.func(*args)
        if (e.is_Mul and len(e.args) == 2 and e.args[0].is_Number and e.args[1].is_Add
                and e.is_commutative):
            return Mul(*e.args)
        return e

    expr = hollow(short)

    numer, denom = expr.as_numer_denom()
    if denom.is_Add:
        n, d = fraction(radsimp(1/denom, symbolic=False, max_terms=1))
        if n is not S.One:
            expr = (numer*n).expand()/d

    if expr.could_extract_minus_sign():
        n, d = fraction(expr)
        if d != 0:
            expr = signsimp(-n/(-d))

    if measure(expr) > ratio*measure(original_expr):
        expr = original_expr

    if floats and rational is None:
        expr = nfloat(expr, exponent=False)

    return done(expr)



# --- solve -------------------------------------------------------------------

def _nc_coefficients(expr):
    """Split ``expr`` into ``{nc word: commutative coefficient}``.

    Distinct words of the free algebra are linearly independent, so ``expr``
    vanishes for commutative unknowns exactly when every coefficient does."""
    from . import Add, Mul, expand
    out = {}
    for t in Add.make_args(expand(expr)):
        c, nc = t.args_cnc()
        word = Mul(*nc)
        out[word] = out.get(word, 0) + Mul(*c)
    return {k: v for k, v in out.items() if v != 0}


def _equation_expr(e):
    e = _sympify(e)
    if getattr(e, 'is_Equality', False):
        return e.lhs - e.rhs
    return e


def _occurrences(word, x):
    """Positions of ``x`` among the factors of the nc word ``word``."""
    from . import Mul
    factors = list(Mul.make_args(word))
    hits = []
    for i, f in enumerate(factors):
        if f == x:
            hits.append(i)
        elif f.has(x):
            b, e = f.as_base_exp()
            if b == x and e.is_Integer and e > 1:
                raise _polynomial_error()
            raise NotImplementedError('solving for %s inside %s is not implemented' % (x, f))
    return factors, hits


def _solve_nc_symbol(expr, x):
    """Solve ``expr = 0`` for the noncommutative unknown ``x`` when it is a
    one-sided linear equation ``L*x = r`` or ``x*R = r``."""
    from . import Add, Mul
    from .utilities.iterables import sift
    left, right, rest = [], [], []
    for t in Add.make_args(_core_expand(expr)):
        c, nc = t.args_cnc()
        if not t.has(x):
            rest.append(t)
            continue
        if any(f.has(x) for f in c):
            raise NotImplementedError('solving for %s inside %s is not implemented' % (x, t))
        factors, hits = _occurrences(Mul(*nc), x)
        if len(hits) != 1:
            raise _polynomial_error()
        i = hits[0]
        left.append((Mul(*c), Mul(*factors[:i]), Mul(*factors[i + 1:])))
    rhs = -Add(*rest)
    if len(left) == 1:
        c, l, r = left[0]
        return l**-1*rhs*r**-1/c
    if all(r == 1 for _, _, r in left):
        coeff = Add(*[c*l for c, l, _ in left])
        if coeff == 0:
            return None
        return coeff**-1*rhs if not coeff.is_commutative else rhs/coeff
    if all(l == 1 for _, l, _ in left):
        coeff = Add(*[c*r for c, _, r in left])
        if coeff == 0:
            return None
        return rhs*coeff**-1 if not coeff.is_commutative else rhs/coeff
    raise NotImplementedError('two-sided (Sylvester) equations in a noncommutative unknown are not implemented')


def _core_expand(e):
    from . import expand
    return expand(e)


def solve(native, f, *symbols, **flags):
    """``solve`` for noncommutative equations.

    A commutative unknown must take a commutative value: each nc equation is
    replaced by the commutative equations ``coefficient(word) = 0`` of its
    distinct noncommutative words. A single noncommutative unknown is solved
    from a one-sided linear equation without commuting factors."""
    from . import Symbol
    many = isinstance(f, (list, tuple, set))
    eqs = [_equation_expr(e) for e in (f if many else [f])]
    if len(symbols) == 1 and isinstance(symbols[0], (list, tuple, set)):
        syms = list(symbols[0])
    else:
        syms = list(symbols)
    if not syms:
        from . import ordered
        syms = list(ordered(set().union(*[e.free_symbols for e in eqs])))
    nc_syms = [s for s in syms if s.is_commutative is False]
    if nc_syms:
        if many or len(syms) != 1:
            raise NotImplementedError('systems in noncommutative unknowns are not implemented')
        x = syms[0]
        sol = _solve_nc_symbol(eqs[0], x)
        if sol is None:
            raise NotImplementedError('the coefficient of %s vanishes identically' % x)
        if flags.get('dict'):
            return [{x: sol}]
        if flags.get('set'):
            return [x], {(sol,)}
        return [sol]
    # Commutative unknowns: coefficient equations of the free algebra.
    new = []
    for e in eqs:
        if is_nc(e):
            if any(w.has(*syms) for w in _nc_coefficients(e) if w != 1):
                raise NotImplementedError('a commutative unknown inside a noncommutative factor is not supported')
            new.extend(_nc_coefficients(e).values())
        else:
            new.append(e)
    if many or len(syms) > 1:
        return native(new, *symbols, **flags)
    if not new:
        return native(0, *symbols, **flags)
    first, others = new[0], new[1:]
    result = native(first, *symbols, **flags)
    if not others:
        return result
    return _filter_solutions(result, others, syms)


def _filter_solutions(result, others, syms):
    from . import simplify

    def ok(sub):
        return all(simplify(o.xreplace(sub)) == 0 for o in others)

    def as_sub(v):
        if isinstance(v, dict):
            return v
        if isinstance(v, tuple):
            return dict(zip(syms, v))
        return {syms[0]: v}

    if isinstance(result, list):
        return [r for r in result if ok(as_sub(r))]
    if isinstance(result, dict):
        return result if ok(result) else []
    if isinstance(result, tuple) and len(result) == 2:
        names, values = result
        return names, {v for v in values if ok(dict(zip(names, v)))}
    return result


# --- installation --------------------------------------------------------------

def install(cancel_, factor_, apart_, together_, collect_):
    """Wrap the native polynomial entry points with the noncommutative branch."""
    import functools
    _native.update(cancel=cancel_, factor=factor_, apart=apart_, together=together_, collect=collect_)

    @functools.wraps(cancel_)
    def cancel_entry(f, *gens, **args):
        if is_nc(f):
            return cancel(f, *gens, **args)
        return cancel_(f, *gens, **args)

    @functools.wraps(factor_)
    def factor_entry(f, *gens, **args):
        # The reference factors lists and tuples entrywise (Beam3D doctest).
        if isinstance(f, (list, tuple)):
            return type(f)(factor_entry(e, *gens, **args) for e in f)
        if is_nc(f):
            return factor(f, *gens, **args)
        return factor_(f, *gens, **args)

    @functools.wraps(apart_)
    def apart_entry(f, x=None, full=False, **options):
        if isinstance(f, (list, tuple)):
            return type(f)(apart_entry(e, x, full=full, **options) for e in f)
        if is_nc(f):
            return apart(f, x=x, full=full, **options)
        return apart_(f, x, full=full, **options)

    @functools.wraps(together_)
    def together_entry(expr, deep=False, fraction=True, **kwargs):
        if isinstance(expr, (list, tuple)):
            return type(expr)(together_entry(e, deep=deep, fraction=fraction, **kwargs) for e in expr)
        if is_nc(expr):
            return together(expr, deep=deep, fraction=fraction)
        if fraction is True and not kwargs:
            return together_(expr, deep=deep)
        return together_(expr, deep=deep, fraction=fraction, **kwargs)

    @functools.wraps(collect_)
    def collect_entry(expr, syms, *args, **kwargs):
        if is_nc(expr):
            collect_check(expr, syms)
        return collect_(expr, syms, *args, **kwargs)

    return cancel_entry, factor_entry, apart_entry, together_entry, collect_entry
