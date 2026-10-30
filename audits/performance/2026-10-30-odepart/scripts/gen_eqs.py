"""Equations for the odesolve differential tests (shared by SymPy and symox).

`cases(lib)` returns a list of (name, eq_string, hint) where eq_string is
evaluated in the namespace of `lib` (sympy or symox) with f, x, a, b, k, p
defined. The list is deterministic (fixed seed).
"""
import random


def _poly_from_roots(rng):
    """Integer characteristic polynomial coefficients (low degree first)."""
    kinds = rng.choice(['int', 'rat', 'sqrt', 'cplx', 'mixed', 'rep', 'coef'])
    if kinds == 'coef':
        n = rng.randint(1, 5)
        c = [rng.randint(-4, 4) for _ in range(n)] + [rng.choice([1, 1, 2, -1, 3])]
        return c
    # build as a product of factors with integer coefficients
    factors = []
    for _ in range(rng.randint(1, 3)):
        t = rng.choice(['lin', 'lin', 'quad_real', 'quad_cplx', 'quad_irr']) if kinds == 'mixed' else {
            'int': 'lin', 'rat': 'rlin', 'sqrt': 'quad_irr', 'cplx': 'quad_cplx', 'rep': 'lin'}[kinds]
        if t == 'lin':
            factors.append([rng.randint(-3, 3), 1])
        elif t == 'rlin':
            factors.append([rng.randint(-3, 3), rng.randint(1, 3)])
        elif t == 'quad_irr':
            factors.append([-rng.choice([2, 3, 5, 6, 7]), 0, 1])
        elif t == 'quad_real':
            factors.append([rng.randint(-3, 3), rng.randint(-5, 5), 1])
        else:
            b = rng.randint(-3, 3)
            factors.append([b*b + rng.randint(1, 9), 2*b, 1])
        if kinds == 'rep' and rng.random() < 0.6:
            factors.append(list(factors[-1]))
    c = [1]
    for g in factors:
        out = [0]*(len(c) + len(g) - 1)
        for i, u in enumerate(c):
            for j, v in enumerate(g):
                out[i + j] += u*v
        c = out
    if len(c) > 7:
        c = c[:7]
        if c[-1] == 0:
            c[-1] = 1
    return c


def _ode(c):
    terms = []
    for i, v in enumerate(c):
        if v == 0:
            continue
        d = 'f(x)' if i == 0 else 'f(x).diff(x, %d)' % i
        terms.append('(%s)*%s' % (v, d))
    return ' + '.join(terms)


def _euler(c):
    terms = []
    for i, v in enumerate(c):
        if v == 0:
            continue
        d = 'f(x)' if i == 0 else 'x**%d*f(x).diff(x, %d)' % (i, i)
        terms.append('(%s)*%s' % (v, d))
    return ' + '.join(terms)


RHS = ['exp(x)', 'exp(2*x)', 'x**2', 'sin(x)', 'cos(2*x)', 'x*exp(-x)', 'exp(x)*sin(x)', '3', 'x']
RHS_VOP = ['exp(x)', 'x', 'sin(x)', '1/x', 'exp(x)/x', 'log(x)']
SYMBOLIC = [
    'f(x).diff(x, 2) + a*f(x)', 'f(x).diff(x, 2) - a**2*f(x)', 'f(x).diff(x, 2) + a*f(x).diff(x) + b*f(x)',
    'f(x).diff(x, 2) + k**2*f(x)', 'f(x).diff(x, 2) - k**2*f(x)', 'f(x).diff(x, 2) + 2*k*f(x).diff(x) + k**2*f(x)',
    'f(x).diff(x, 2) + p*f(x)', 'f(x).diff(x, 2) - p*f(x)', 'a*f(x).diff(x, 3) + b*f(x)',
    'f(x).diff(x, 3) - a*f(x)', 'f(x).diff(x, 4) - a**4*f(x)', 'f(x).diff(x, 2) + (a + b)*f(x).diff(x) + a*b*f(x)',
    'f(x).diff(x, 2) + I*f(x)', 'f(x).diff(x, 2) + sqrt(2)*f(x)', 'f(x).diff(x, 2) + 2.0*f(x)',
    'f(x).diff(x, 2) + pi*f(x)', 'f(x).diff(x, 3) + I*f(x).diff(x) + f(x)',
    'f(x).diff(x, 5) + 10*f(x).diff(x) - 2*f(x)', 'f(x).diff(x, 5) - f(x).diff(x) - f(x)',
    'f(x).diff(x, 5) + f(x).diff(x, 4) + 2*f(x) - 3*f(x).diff(x)', '2*f(x).diff(x, 5) + 20*f(x).diff(x) - 4*f(x)',
    '-f(x).diff(x, 5) - 10*f(x).diff(x) + 2*f(x)', 'f(x).diff(x, 6) + f(x).diff(x) + 1*f(x)',
    '(f(x).diff(x, 5) + 10*f(x).diff(x) - 2*f(x))*(1) ', 'f(x).diff(x, 6) - 3*f(x).diff(x,5) + f(x)',
    'f(x).diff(x, 3) - 3*f(x).diff(x) + 1*f(x)', 'f(x).diff(x, 4) + f(x).diff(x) + f(x)',
    'f(x).diff(x, 3) + f(x).diff(x) + f(x)',
]
EULER_SYM = ['x**2*f(x).diff(x, 2) + a*f(x)', 'x**2*f(x).diff(x, 2) + x*f(x).diff(x) - a**2*f(x)',
             'x**3*f(x).diff(x, 3) + 2*f(x)', 'x**3*f(x).diff(x, 3) - x*f(x).diff(x) + 5*f(x)',
             'x**6*f(x).diff(x, 6) - x*f(x).diff(x) + f(x)', 'x**2*f(x).diff(x, 2) + x*f(x).diff(x) + f(x)',
             'x**2*f(x).diff(x, 2) - 3*x*f(x).diff(x) + 13*f(x)', 'x**4*f(x).diff(x, 4) + f(x)',
             'Derivative(x*f(x), x, x, x)', 'x**2*f(x).diff(x, 2) + 2.5*f(x)']


def cases(n_const=160, n_euler=80, n_undet=60, n_vop=30, seed=20261030):
    rng = random.Random(seed)
    out = []
    for k in range(n_const):
        c = _poly_from_roots(rng)
        while len(c) < 2 or all(v == 0 for v in c[1:]):
            c = _poly_from_roots(rng)
        out.append(('const_%03d' % k, _ode(c), 'nth_linear_constant_coeff_homogeneous'))
    for k, e in enumerate(SYMBOLIC):
        out.append(('csym_%02d' % k, e, 'nth_linear_constant_coeff_homogeneous'))
    for k in range(n_euler):
        c = _poly_from_roots(rng)
        while len(c) < 2 or all(v == 0 for v in c[1:]):
            c = _poly_from_roots(rng)
        out.append(('euler_%03d' % k, _euler(c), 'nth_linear_euler_eq_homogeneous'))
    for k, e in enumerate(EULER_SYM):
        out.append(('esym_%02d' % k, e, 'nth_linear_euler_eq_homogeneous'))
    for k in range(n_undet):
        c = _poly_from_roots(rng)
        while len(c) < 2 or len(c) > 5 or all(v == 0 for v in c[1:]):
            c = _poly_from_roots(rng)
        out.append(('undet_%03d' % k, _ode(c) + ' - (%s)' % rng.choice(RHS),
                    'nth_linear_constant_coeff_undetermined_coefficients'))
    for k in range(n_vop):
        c = _poly_from_roots(rng)
        while len(c) < 2 or len(c) > 3 or all(v == 0 for v in c[1:]):
            c = _poly_from_roots(rng)
        out.append(('vop_%03d' % k, _ode(c) + ' - (%s)' % rng.choice(RHS_VOP),
                    'nth_linear_constant_coeff_variation_of_parameters'))
    return out
