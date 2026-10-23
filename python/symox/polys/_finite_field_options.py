"""``domain=GF(p)`` for the expression-level polynomial functions.

The reference ``Options`` treats ``domain=FF(p)`` (or the string
``'GF(p)'``/``'FF(p)'``) and ``modulus=p`` as the same option, so the
functions below rewrite the domain into the modulus. ``cancel`` over a
finite field is computed here with Poly arithmetic modulo p: the reference
``Poly.cancel`` divides out the gcd and makes the denominator monic.
"""
import functools
import re

_NAMES = ('factor', 'factor_list', 'sqf', 'sqf_list', 'gcd', 'lcm', 'cancel',
          'div', 'quo', 'rem', 'gcdex', 'invert')
_PATTERN = re.compile(r'^\s*(?:GF|FF)\s*\(\s*(\d+)\s*\)\s*$')


def _finite_modulus(domain):
    if isinstance(domain, str):
        match = _PATTERN.match(domain)
        return int(match.group(1)) if match else None
    if getattr(domain, 'is_FiniteField', False):
        return int(domain.characteristic())
    return None


def _rewrite(args):
    if 'domain' in args:
        p = _finite_modulus(args['domain'])
        if p is not None:
            if 'modulus' in args:
                return args
            args = dict(args)
            del args['domain']
            args['modulus'] = p
    return args


_BINARY = ('gcd', 'lcm', 'div', 'quo', 'rem', 'gcdex', 'invert')


def _relabel_roots(name, fargs, args):
    """Root objects whose defining polynomial uses an inferred generator.

    SymPy reads ``CRootOf(x**3 - x - 1, 0)`` in ``gcd(x**3 - x - 1, x -
    CRootOf(x**3 - x - 1, 0), extension=True)`` as an algebraic coefficient
    when the generators are inferred: only ``Poly`` with explicit generators
    rejects it. The backend checks the generators structurally, so each such
    root is rewritten over a fresh symbol and restored in the result.
    Returns the new arguments and the restoring map, or None.
    """
    from .._core import Expr
    count = 2 if name in _BINARY else 1
    if len(fargs) != count or 'gens' in args or args.get('extension') is None:
        return None
    if not all(isinstance(a, Expr) for a in fargs):
        return None
    from ..polys.rootoftools import ComplexRootOf
    from .. import Dummy
    roots = set()
    for a in fargs:
        roots |= a.atoms(ComplexRootOf)
    if not roots:
        return None
    free = set()
    for a in fargs:
        free |= a.free_symbols
    forward, back = {}, {}
    for r in roots:
        symbols = r.expr.free_symbols
        if not symbols & free:
            continue
        t = Dummy('t')
        new = ComplexRootOf(r.expr.xreplace({s: t for s in symbols}), r.index)
        forward[r] = new
        back[new] = r
    if not forward:
        return None
    return tuple(a.xreplace(forward) for a in fargs), back


def _restore_roots(value, back):
    from .._core import Expr
    if isinstance(value, Expr):
        return value.xreplace(back)
    if isinstance(value, (tuple, list)):
        return type(value)(_restore_roots(v, back) for v in value)
    return value


def _cancel_mod(native, f, gens, args):
    from .. import Poly, sympify, together, fraction, Tuple
    p = args['modulus']
    f = sympify(f)
    if isinstance(f, (tuple, Tuple)) or not f.free_symbols or f.is_Relational:
        return native(f, *gens, **args)
    numer, denom = fraction(together(f))
    allgens = gens or tuple(sorted((numer.free_symbols | denom.free_symbols), key=str))
    if len(allgens) != 1:
        return native(f, *gens, **args)
    try:
        P = Poly(numer, *allgens, modulus=p)
        Q = Poly(denom, *allgens, modulus=p)
    except Exception:
        return native(f, *gens, **args)
    if Q.is_zero:
        return native(f, *gens, **args)
    G = P.gcd(Q)
    P, Q = P.quo(G), Q.quo(G)
    lc = Q.LC()
    inv = pow(int(lc) % p, -1, p)
    P, Q = P.mul_ground(inv), Q.monic()
    return P.as_expr() / Q.as_expr()


def install(namespace):
    from ._modular_functions import NAMES as ring_names, selected, call
    for name in ring_names:
        native = namespace.get(name)
        if native is None or getattr(native, '_symox_ring_options', False):
            continue

        def make_ring(native, name):
            @functools.wraps(native)
            def wrapper(*fargs, **args):
                dom = selected(fargs, args)
                if dom is not None:
                    return call(name, dom, fargs, args)
                return native(*fargs, **args)
            wrapper._symox_ring_options = True
            return wrapper
        namespace[name] = make_ring(native, name)
    for name in _NAMES:
        native = namespace.get(name)
        if native is None or getattr(native, '_symox_ff_options', False):
            continue

        def make(native, name):
            @functools.wraps(native)
            def wrapper(*fargs, **args):
                args = _rewrite(args)
                if name == 'cancel' and 'modulus' in args and fargs and \
                        args.get('symmetric', True) is not False and \
                        set(args) <= {'modulus', 'symmetric'}:
                    return _cancel_mod(native, fargs[0], fargs[1:], args)
                relabeled = _relabel_roots(name, fargs, args)
                if relabeled is not None:
                    return _restore_roots(native(*relabeled[0], **args), relabeled[1])
                return native(*fargs, **args)
            wrapper._symox_ff_options = True
            return wrapper
        namespace[name] = make(native, name)
