"""Reference parametric solution of a*x**2 + b*y**2 + c*z**2 = 0.

Private port of SymPy 1.14 ``diophantine`` for a homogeneous ternary
quadratic in normal form (``sympy/solvers/diophantine/diophantine.py``:
``HomogeneousTernaryQuadraticNormal.solve``, ``sqf_normal``, ``descent``,
``ldescent``, ``gaussian_reduce``, ``holzer``, ``reconstruct`` and
``_parametrize_ternary_quadratic``).  It is used only by the reference
``ImplicitRegion._regular_point_ellipse`` so that the rational point (and so
the rational parametrization) of a conic is the reference one.

The reference two-parameter family is in general not the complete integer
solution set, which is why the public ``diophantine`` does not return it.
"""
from math import gcd as _gcd
from fractions import Fraction
from functools import reduce

from ..core import Integer, Rational, expand, Symbol
from ..ntheory import sqrt_mod, factorint, divisors, integer_nthroot


def _sign(n):
    return (n > 0) - (n < 0)


def _igcd(*args):
    return reduce(_gcd, (int(a) for a in args))


def _ilcm(*args):
    out = 1
    for a in args:
        a = int(a)
        out = out*a//_gcd(out, a) if a else 0
    return abs(out)


def _gcdext(a, b):
    # sympy.external.ntheory.gcdext
    if not a or not b:
        g = abs(a) or abs(b)
        if not g:
            return (0, 0, 0)
        return (g, a // g, b // g)
    x_sign, a = (-1, -a) if a < 0 else (1, a)
    y_sign, b = (-1, -b) if b < 0 else (1, b)
    x, r = 1, 0
    y, s = 0, 1
    while b:
        q, c = divmod(a, b)
        a, b = b, c
        x, r = r, x - q*r
        y, s = s, y - q*s
    return (a, x * x_sign, y * y_sign)


def _remove_gcd_int(*x):
    g = _igcd(*x)
    if g in (0, 1):
        return tuple(x)
    return tuple(i//g for i in x)


def _rational_pq(a, b):
    return _remove_gcd_int(_sign(b)*a, abs(b))


def _nint_or_floor(p, q):
    w, r = divmod(p, q)
    if abs(r) <= abs(q)//2:
        return w
    return w + 1


def _base_solution_linear(c, a, b):
    a, b, c = _remove_gcd_int(a, b, c)
    if c == 0:
        return (0, 0)
    d, x0, y0 = _gcdext(abs(a), abs(b))
    x0 *= _sign(a)
    y0 *= _sign(b)
    if c % d:
        return (None, None)
    return (c*x0, c*y0)


def _square_factor(a):
    out = 1
    for p, e in factorint(int(a)).items():
        out *= int(p)**(int(e)//2)
    return out


def _sqrt_mod(a, p):
    r = sqrt_mod(int(a), int(p))
    return None if r is None else int(r)


def sqf_normal(a, b, c):
    ABC = _remove_gcd_int(a, b, c)
    sq = tuple(_square_factor(i) for i in ABC)
    sqf = A, B, C = tuple(i//j**2 for i, j in zip(ABC, sq))
    # the reference divides exactly (Rational), not with floor division
    A, B, C = Fraction(A), Fraction(B), Fraction(C)
    pc = _igcd(A, B)
    A /= pc
    B /= pc
    pa = _igcd(B, C)
    B /= pa
    C /= pa
    pb = _igcd(A, C)
    A /= pb
    B /= pb
    A *= pa
    B *= pb
    C *= pc
    if any(v.denominator != 1 for v in (A, B, C)):
        raise NotImplementedError('non-integral square-free normal form')
    return (sq, sqf, (int(A), int(B), int(C)))


def reconstruct(A, B, z):
    for p, e in factorint(_igcd(A, B)).items():
        if e != 1:
            raise ValueError('a and b should be square-free')
        z *= int(p)
    return z


def ldescent(A, B):
    if A == 0 or B == 0:
        raise ValueError("A and B must be non-zero integers")
    if abs(A) > abs(B):
        w, y, x = ldescent(B, A)
        return w, x, y
    if A == 1:
        return (1, 1, 0)
    if B == 1:
        return (1, 0, 1)
    if B == -1:
        return
    r = _sqrt_mod(A, B)
    if r is None:
        return
    Q = (r**2 - A) // B
    if Q == 0:
        return r, -1, 0
    for i in divisors(Q):
        i = int(i)
        d, _exact = integer_nthroot(abs(Q) // i, 2)
        if _exact:
            B_0 = _sign(Q)*i
            W, X, Y = ldescent(A, B_0)
            return _remove_gcd_int(-A*X + r*W, r*X - W, Y*B_0*int(d))


def gaussian_reduce(w, a, b):
    a = abs(a)

    def _dot(u, v):
        return u[0]*v[0] + a*u[1]*v[1]

    u = (b, 0)
    v = (w, 1) if b*w >= 0 else (-w, -1)
    if b**2 < w**2 + a:
        u, v = v, u
    while _dot(u, u) > (dv := _dot(v, v)):
        k = _dot(u, v) // dv
        u, v = v, (u[0] - k*v[0], u[1] - k*v[1])
    c = (v[0] - u[0], v[1] - u[1])
    if _dot(c, c) <= _dot(u, u) <= 2*_dot(u, v):
        return c
    return u


def descent(A, B):
    if abs(A) > abs(B):
        x, y, z = descent(B, A)
        return x, z, y
    if B == 1:
        return (1, 0, 1)
    if A == 1:
        return (1, 1, 0)
    if B == -A:
        return (0, 1, 1)
    if B == A:
        x, z, y = descent(-1, A)
        return (A*y, z, x)
    w = _sqrt_mod(A, B)
    x_0, z_0 = gaussian_reduce(w, A, B)
    t = (x_0**2 - A*z_0**2) // B
    t_2 = _square_factor(t)
    t_1 = t // t_2**2
    x_1, z_1, y_1 = descent(A, t_1)
    return _remove_gcd_int(x_0*x_1 + A*z_0*z_1, z_0*x_1 + x_0*z_1, t_1*t_2*y_1)


def holzer(x, y, z, a, b, c):
    if c % 2:
        k = 2*c
    else:
        k = c//2
    small = a*b*c
    step = 0
    while True:
        t1, t2, t3 = a*x**2, b*y**2, c*z**2
        if t1 + t2 != t3:
            if step == 0:
                raise ValueError('bad starting solution')
            break
        x_0, y_0, z_0 = x, y, z
        if max(t1, t2, t3) <= small:
            break
        u, v = _base_solution_linear(k, y_0, -x_0)
        if u is None or v is None:
            break
        p, q = -(a*u*x_0 + b*v*y_0), c*z_0
        r = Rational(p, q)
        if c % 2 == 0:
            w = _nint_or_floor(p, q)
            assert abs(w - r) <= Rational(1, 2)
        else:
            w = p//q
            if (a*u + b*v + c*w) % 2:
                w += 1
            assert abs(w - r) <= 1
        A = (a*u**2 + b*v**2 + c*w**2)
        B = (a*u*x_0 + b*v*y_0 + c*w*z_0)
        nx, ny, nz = x_0*A - 2*u*B, y_0*A - 2*v*B, z_0*A - 2*w*B
        assert nx % k == 0 and ny % k == 0 and nz % k == 0
        x, y, z = nx//k, ny//k, nz//k
        step += 1
    return (int(x_0), int(y_0), int(z_0))


def _solve_normal(a, b, c):
    """``HomogeneousTernaryQuadraticNormal(a*x**2 + b*y**2 + c*z**2).solve()``."""
    (sqf_of_a, sqf_of_b, sqf_of_c), (a_1, b_1, c_1), (a_2, b_2, c_2) = sqf_normal(a, b, c)
    A = -a_2*c_2
    B = -b_2*c_2
    if A < 0 and B < 0:
        return None
    if (_sqrt_mod(-b_2*c_2, a_2) is None or _sqrt_mod(-c_2*a_2, b_2) is None
            or _sqrt_mod(-a_2*b_2, c_2) is None):
        return None
    z_0, x_0, y_0 = descent(A, B)
    z_0, q = _rational_pq(z_0, abs(c_2))
    x_0 *= q
    y_0 *= q
    x_0, y_0, z_0 = _remove_gcd_int(x_0, y_0, z_0)
    if _sign(a) == _sign(b):
        x_0, y_0, z_0 = holzer(x_0, y_0, z_0, abs(a_2), abs(b_2), abs(c_2))
    elif _sign(a) == _sign(c):
        x_0, z_0, y_0 = holzer(x_0, z_0, y_0, abs(a_2), abs(c_2), abs(b_2))
    else:
        y_0, z_0, x_0 = holzer(y_0, z_0, x_0, abs(b_2), abs(c_2), abs(a_2))
    x_0 = reconstruct(b_1, c_1, x_0)
    y_0 = reconstruct(a_1, c_1, y_0)
    z_0 = reconstruct(a_1, b_1, z_0)
    sq_lcm = _ilcm(sqf_of_a, sqf_of_b, sqf_of_c)
    x_0 = abs(x_0*sq_lcm // sqf_of_a)
    y_0 = abs(y_0*sq_lcm // sqf_of_b)
    z_0 = abs(z_0*sq_lcm // sqf_of_c)
    return _remove_gcd_int(x_0, y_0, z_0)


def _parametrize(solution, coeffs):
    # _parametrize_ternary_quadratic for a*x**2 + b*y**2 + c*z**2
    x_0, y_0, z_0 = solution
    if list(solution).count(0) >= 2:
        return None
    if x_0 == 0:
        a, b, c = coeffs
        y_p, x_p, z_p = _parametrize((y_0, x_0, z_0), (b, a, c))
        return x_p, y_p, z_p
    a, b, c = coeffs
    r, p, q = (Symbol(n, integer=True) for n in ('r', 'p', 'q'))
    eq_1 = expand(a*(r*x_0)**2 + b*(r*y_0 + p)**2 + c*(r*z_0 + q)**2)
    A, B = eq_1.as_independent(r, as_Add=True)
    x = expand(A*x_0)
    y = expand(A*y_0 - expand(B/r*p))
    z = expand(A*z_0 - expand(B/r*q))
    contents = [Integer(e.as_content_primitive()[0]) for e in (x, y, z) if e != 0]
    g = _igcd(*contents) if len(contents) >= 2 else 1
    if g > 1:
        x, y, z = (expand(e/g) for e in (x, y, z))
    return x, y, z


def ternary_normal_diophantine(a, b, c):
    """The set the reference ``diophantine(a*x**2 + b*y**2 + c*z**2)`` returns.

    ``a``, ``b`` and ``c`` are nonzero integers.  The result is a set holding
    the reference two-parameter tuple (parameters ``p`` and ``q``), or the
    trivial solution ``(0, 0, 0)`` when no nontrivial one exists.
    """
    a, b, c = int(a), int(b), int(c)
    if not (a > 0 and b and c):
        # the only use (conic regular points) has a > 0; the reference
        # sign handling of a negative leading coefficient is not ported
        raise ValueError('a positive x**2 coefficient and nonzero b, c are required')
    # factor_list: unit content
    g = _igcd(a, b, c)
    a, b, c = a//g, b//g, c//g
    # signsimp(base, evaluate=False).as_coeff_Mul() of the diophantine driver
    x, y, z = (Symbol(n) for n in 'xyz')
    if (a*x**2 + b*y**2 + c*z**2).could_extract_minus_sign():
        a, b, c = -a, -b, -c
    sol = _solve_normal(a, b, c)
    param = None if sol is None else _parametrize(sol, (a, b, c))
    if param is None:
        zero = Integer(0)
        return {(zero, zero, zero)}
    return {tuple(param)}
