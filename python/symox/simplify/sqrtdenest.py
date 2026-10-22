"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'sqrtdenest')

import symox as S

class SqrtdenestStopIteration(Exception):pass

def is_sqrt(expr):
    """Return True if expr is a sqrt, otherwise False."""
    return bool(getattr(expr, 'is_Pow', False) and expr.exp.is_Rational
                and abs(expr.exp) == S.S.Half)


def sqrt_depth(p):
    """Return the maximum depth of any square root argument of p."""
    p = S.sympify(p)
    if p == S.I:
        return 1
    if p.is_Atom:
        return 0
    if p.is_Add or p.is_Mul:
        return max(sqrt_depth(x) for x in p.args)
    if is_sqrt(p):
        return sqrt_depth(p.base) + 1
    return 0


def is_algebraic(p):
    """Return True if p is comprised of only Rationals or square roots
    of Rationals and algebraic operations."""
    p = S.sympify(p)
    if p.is_Rational:
        return True
    elif p.is_Atom:
        return False
    elif is_sqrt(p) or p.is_Pow and p.exp.is_Integer:
        return is_algebraic(p.base)
    elif p.is_Add or p.is_Mul:
        return all(is_algebraic(x) for x in p.args)
    return False


def _mexpand(expr, recursive=False):
    if expr is None:
        return None
    was = None
    while was != expr:
        was, expr = expr, S.expand_mul(S.expand_multinomial(expr))
        if not recursive:
            break
    return expr


def _sqrt_match(p):
    """Return [a, b, r] for p.match(a + b*sqrt(r)) where, in addition to
    matching, sqrt(r) also has then maximal sqrt_depth among addends of p
    (SymPy 1.14 algorithm)."""
    from .radsimp import split_surds
    from ..core.sorting import default_sort_key
    p = _mexpand(S.sympify(p))
    if p.is_Number:
        res = (p, S.S.Zero, S.S.Zero)
    elif p.is_Add:
        pargs = sorted(p.args, key=default_sort_key)
        sqargs = [x**2 for x in pargs]
        if all(sq.is_Rational and sq.is_positive for sq in sqargs):
            r, b, a = split_surds(p)
            return [a, b, r]
        v = [(sqrt_depth(x), x, i) for i, x in enumerate(pargs)]
        nmax = max(v, key=default_sort_key)
        if nmax[0] == 0:
            res = []
        else:
            depth, _, i = nmax
            r = pargs.pop(i)
            v.pop(i)
            b = S.S.One
            if r.is_Mul:
                bv = []
                rv = []
                for x in r.args:
                    if sqrt_depth(x) < depth:
                        bv.append(x)
                    else:
                        rv.append(x)
                b = S.Mul._from_args(bv)
                r = S.Mul._from_args(rv)
            a1 = []
            b1 = [b]
            for x in v:
                if x[0] < depth:
                    a1.append(x[1])
                else:
                    x1 = x[1]
                    if x1 == r:
                        b1.append(1)
                    elif x1.is_Mul:
                        x1args = list(x1.args)
                        if r in x1args:
                            x1args.remove(r)
                            b1.append(S.Mul(*x1args))
                        else:
                            a1.append(x[1])
                    else:
                        a1.append(x[1])
            a = S.Add(*a1)
            b = S.Add(*b1)
            res = (a, b, r**2)
    else:
        b, r = p.as_coeff_Mul()
        if is_sqrt(r):
            res = (S.S.Zero, b, r**2)
        else:
            res = []
    return list(res)


def sqrt_biquadratic_denest(expr, a, b, r, d2):
    """Denest expr = sqrt(a + b*sqrt(r)) with a, b, r linear combinations
    of square roots of positive rationals on the rationals (SymPy 1.14)."""
    from .radsimp import radsimp, rad_rationalize
    if r <= 0 or d2 < 0 or not b or sqrt_depth(expr.base) < 2:
        return None
    for x in (a, b, r):
        for y in x.args:
            y2 = y**2
            if not y2.is_Integer or not y2.is_positive:
                return None
    sqd = _mexpand(S.sqrtdenest(S.sqrt(radsimp(d2))))
    if sqrt_depth(sqd) > 1:
        return None
    x1, x2 = [a/2 + sqd/2, a/2 - sqd/2]
    for x in (x1, x2):
        A = S.sqrtdenest(S.sqrt(x))
        if sqrt_depth(A) > 1:
            continue
        Bn, Bd = rad_rationalize(b, _mexpand(2*A))
        B = Bn/Bd
        z = A + B*S.sqrt(r)
        if z < 0:
            z = -z
        return _mexpand(z)
    return None

__all__ += ['_mexpand','_sqrt_match','SqrtdenestStopIteration','is_algebraic','is_sqrt','sqrt_biquadratic_denest','sqrt_depth']
