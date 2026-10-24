"""One-sided limits onto a branch cut.

The reference computes these through each function's leading term with a
direction (``cdir``); a plain substitution gives the principal value, which
is the limit from the other side only.  ``limit(log(I*x - 1), x, 0, '-')``
is ``-I*pi``, not ``log(-1) = I*pi``.

``branch_side_rewrite`` finds the nodes whose argument tends to a point of
the cut from the non-principal side (decided numerically on the approach
path) and replaces them:

- ``log(a)`` by ``log(-a) + s*I*pi`` and ``a**r`` by ``(-a)**r*exp(s*I*pi*r)``
  (``s`` the side of ``a``), identities that hold on that whole side, so any
  surrounding expression keeps its limit;
- an inverse trigonometric node by its one-sided value (``pi - asin(2)``),
  only when that makes the expression free of the variable.

Every replacement is checked numerically against the function on the
approach path, so a wrong side never leaks a value.
"""
from .. import (I, Pow, Rational, S, acos, acot, acsc, asec, asin, atan, im, log, oo, pi, re, zoo, nan,
                preorder_traversal)

_H = Rational(1, 10**30)
_TOL = 1e-12

_REAL_CUT = (asin, acos, asec, acsc)
_IMAG_CUT = (atan, acot)


def _close(u, v):
    try:
        d = complex(u) - complex(v)
    except (TypeError, ValueError):
        return False
    return abs(d) <= _TOL * (1 + abs(complex(v)))


def _num(e, prec=40):
    try:
        v = e.evalf(prec)
        return complex(v)
    except (TypeError, ValueError, AttributeError, ZeroDivisionError):
        return None


def _real(e):
    return e.is_extended_real is True


def _on_cut(f, w):
    """Whether the finite number ``w`` lies on the branch cut of ``f``; the
    axis ('re' or 'im') of the displacement that selects the side."""
    if isinstance(f, (log, Pow)):
        return 'im' if _real(w) and w.is_negative is True else None
    if isinstance(f, (asin, acos)):
        if _real(w) and ((w > 1) == True or (w < -1) == True):
            return 'im'
        return None
    if isinstance(f, (asec, acsc)):
        if _real(w) and w.is_zero is False and ((w < 1) == True and (w > -1) == True):
            return 'im'
        return None
    if isinstance(f, (atan, acot)):
        t = -I*w
        if not _real(t):
            return None
        if isinstance(f, atan) and (abs(t) > 1) == True:
            return 're'
        if isinstance(f, acot) and (abs(t) < 1) == True:
            return 're'
    return None


def _candidates(f, w):
    if isinstance(f, asin):
        return [pi - asin(w), -pi - asin(w)]
    if isinstance(f, acos):
        return [-acos(w), 2*pi - acos(w)]
    if isinstance(f, asec):
        return [-asec(w), 2*pi - asec(w)]
    if isinstance(f, acsc):
        return [pi - acsc(w), -pi - acsc(w)]
    if isinstance(f, atan):
        return [atan(w) - pi, atan(w) + pi]
    if isinstance(f, acot):
        return [acot(w) + pi, acot(w) - pi]
    return []


def _lower_side_at_infinity(e, z, z0):
    """At ``z0 = +-oo``: rewrite each ``log(a)`` and ``a**r`` whose argument
    approaches the negative real axis from below (``im(a) < 0`` on the whole
    tail of the path) by ``log(-a) - I*pi`` and ``(-a)**r*(-1)**(-r)``.
    These are identities on the lower half-plane, and the new arguments stay
    away from the cut; ``log(-x - I)`` at ``oo`` is ``log(x) - I*pi``, not the
    principal ``log(x) + I*pi`` a substitution of the leading term gives."""
    sgn = 1 if z0 == oo else -1
    pts = [sgn*S(10)**k for k in (10, 20, 30)]
    mapping = {}
    for f in preorder_traversal(e):
        if isinstance(f, log):
            if len(f.args) != 1 or not f.has(z):
                continue
            a = f.args[0]
        elif isinstance(f, Pow):
            r = f.exp
            if not (f.base.has(z) and not r.has(z) and r.is_number and r.is_integer is not True and _real(r)):
                continue
            a = f.base
        else:
            continue
        ok = True
        for pt in pts:
            v = _num(a.subs(z, pt), 60)
            if v is None or not (v.real < 0 and v.imag < 0 and abs(v.imag) <= 1e-6*abs(v)):
                ok = False
                break
        if not ok:
            continue
        if isinstance(f, log):
            new = log(-a) - I*pi
        else:
            new = Pow(-a, f.exp)*(-1)**(-f.exp)
        path = a.subs(z, pts[-1])
        exact = _num(Pow(path, f.exp) if isinstance(f, Pow) else log(path), 30)
        if exact is not None and _close(_num(new.subs(z, pts[-1]), 30), exact):
            mapping[f] = new
    if not mapping:
        return None
    return e.xreplace(mapping)


def _symbolic_side(a, z, z0, dir, w, axis):
    """Side of the cut for ``a`` on ``z = z0 +- h`` (``h -> 0+``) from the
    leading term of ``im`` (or ``re``) of ``a - w``, with a path point far
    enough inside that side for the numerical checks."""
    from .. import Dummy, expand
    h = Dummy('h', positive=True)
    try:
        d = expand(a.subs(z, z0 + h if dir == '+' else z0 - h) - w)
        part = expand(im(d) if axis == 'im' else re(d))
        if part.has(im, re) or part.is_zero is not False:
            return None, None
        lt = part.as_leading_term(h)
        c, p = lt.as_coeff_exponent(h)
        if c.has(h) or not c.is_number or not p.is_number or p.is_positive is not True:
            return None, None
        cv = _num(c, 30)
        if cv is None or abs(cv.imag) > _TOL*abs(cv) or cv.real == 0:
            return None, None
        # a step small enough for the leading term to dominate, while
        # the displacement stays well above the evaluation precision
        k = max(4, min(30, int(40/(float(p) + 1))))
        hp = Rational(1, 10**k)
        return (1 if cv.real > 0 else -1), (z0 + hp if dir == '+' else z0 - hp)
    except Exception:
        return None, None


def branch_side_rewrite(e, z, z0, dir):
    if not getattr(z, 'is_Symbol', False):
        return None
    if z0 == oo or z0 == -oo:
        return _lower_side_at_infinity(e, z, z0)
    if dir not in ('+', '-'):
        return None
    try:
        if z0.has(oo, -oo, zoo, nan) or z0.is_finite is not True:
            return None
    except AttributeError:
        return None
    nodes = []
    for node in preorder_traversal(e):
        if isinstance(node, (log,) + _REAL_CUT + _IMAG_CUT):
            if len(node.args) == 1 and node.has(z):
                nodes.append(node)
        elif isinstance(node, Pow):
            r = node.exp
            if node.base.has(z) and not r.has(z) and r.is_number and r.is_integer is not True and _real(r):
                nodes.append(node)
    if not nodes:
        return None
    zp = z0 + (_H if dir == '+' else -_H)
    mapping = {}
    inverse = {}
    for f in nodes:
        a = f.base if isinstance(f, Pow) else f.args[0]
        try:
            w = a.subs(z, z0)
        except Exception:
            continue
        if not w.is_number or w.has(oo, -oo, zoo, nan):
            continue
        axis = _on_cut(f, w)
        if axis is None:
            continue
        an = _num(a.subs(z, zp), 60)
        wn = _num(w, 60)
        if an is None or wn is None:
            continue
        d = an - wn
        side = d.imag if axis == 'im' else d.real
        zq = zp
        if abs(side) <= 1e-50:
            # The displacement off the axis is of higher order than the
            # step (`-I*x**2 + x - 3`): read its sign from the leading term
            # of the exact displacement on the approach path.
            s, zq = _symbolic_side(a, z, z0, dir, w, axis)
            if s is None:
                continue
        else:
            s = 1 if side > 0 else -1
        a_path = a.subs(z, zq)
        exact = _num(Pow(a_path, f.exp) if isinstance(f, Pow) else f.func(a_path), 30)
        principal = _num(Pow(w, f.exp) if isinstance(f, Pow) else f.func(w), 30)
        if exact is None or principal is None or _close(exact, principal):
            continue
        if isinstance(f, log):
            new = log(-a) + s*I*pi
            if _close(_num(new.subs(z, zq), 30), exact):
                mapping[f] = new
            continue
        if isinstance(f, Pow):
            r = f.exp
            phase = (-1)**r if s > 0 else -(-1)**(1 - r)
            new = Pow(-a, r)*phase
            if _close(_num(new.subs(z, zq), 30), exact):
                mapping[f] = new
            continue
        for c in _candidates(f, w):
            if _close(_num(c, 30), exact):
                inverse[f] = c
                break
    if inverse:
        trial = e.xreplace({**mapping, **inverse})
        if not trial.has(z) and not trial.has(zoo, nan):
            return trial
    if not mapping:
        return None
    return e.xreplace(mapping)
