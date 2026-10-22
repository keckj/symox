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
                return native(*fargs, **args)
            wrapper._symox_ff_options = True
            return wrapper
        namespace[name] = make(native, name)
