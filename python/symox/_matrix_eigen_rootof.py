"""Eigenvalues and eigenvectors of rational matrices whose characteristic
polynomial has irreducible factors of degree > 4.

The native kernel stops at quartic factors.  These are ports of SymPy 1.14's
``_eigenvals_dict``/``_eigenvals_list`` (``roots`` then
``PurePoly.all_roots``, i.e. ``CRootOf``) and of ``_eigenvects_DOM``
(``dom_eigenvects`` over ``QQ[lambda]/(minpoly)`` followed by
``dom_eigenvects_to_sympy``) used as the fallback for that case.
"""
from fractions import Fraction

UNSUPPORTED = "irreducible factors of degree > 4"

eigenvals_error_message = \
    "It is not always possible to express the eigenvalues of a matrix " + \
    "of size 5x5 or higher in radicals. " + \
    "We have CRootOf, but domains other than the rationals are not " + \
    "currently supported. " + \
    "If there are no symbols in the matrix, " + \
    "it should still be possible to compute numeric approximations " + \
    "of the eigenvalues using " + \
    "M.evalf().eigenvals() or M.charpoly().nroots()."


def is_unsupported(exc):
    return UNSUPPORTED in str(exc)


def _is_rational_matrix(M):
    return all(getattr(x, 'is_Rational', False) for x in M)


def eigenvals(M, error_when_incomplete=True, multiple=False, simplify=False,
              **flags):
    """``_eigenvals_list``/``_eigenvals_dict`` of ``sympy/matrices/eigen.py``
    (without the final ``simplify`` step, applied by the caller)."""
    from types import FunctionType
    from . import roots, MatrixError
    iblocks = M.strongly_connected_components()
    is_dom = _is_rational_matrix(M)
    all_list = []
    all_dict = {}
    for b in iblocks:
        if is_dom and len(b) == 1:
            val = M[b[0], b[0]]
            all_list.append(val)
            all_dict[val] = all_dict.get(val, 0) + 1
            continue
        block = M.extract(b, b)
        if isinstance(simplify, FunctionType):
            charpoly = block.charpoly(simplify=simplify)
        else:
            charpoly = block.charpoly()
        eigs = roots(charpoly, multiple=multiple, **flags)
        count = len(eigs) if multiple else sum(eigs.values())
        if count != block.rows:
            try:
                eigs = charpoly.all_roots(multiple=multiple)
                if not multiple:
                    eigs = dict(eigs)
            except NotImplementedError:
                if error_when_incomplete:
                    raise MatrixError(eigenvals_error_message)
                eigs = [] if multiple else {}
        if multiple:
            all_list += list(eigs)
        else:
            for k, v in eigs.items():
                all_dict[k] = all_dict.get(k, 0) + v
    return all_list if multiple else all_dict


# --- arithmetic in QQ[l]/(minpoly), coefficients low -> high degree ---------

def _trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def _pmul(a, b):
    if not a or not b:
        return []
    r = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return _trim(r)


def _psub(a, b):
    n = max(len(a), len(b))
    r = [(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
         for i in range(n)]
    return _trim([Fraction(x) for x in r])


def _pdivmod(a, m):
    a = list(a)
    q = [Fraction(0)] * max(len(a) - len(m) + 1, 0)
    lc = m[-1]
    while len(a) >= len(m) and a:
        c = a[-1] / lc
        k = len(a) - len(m)
        q[k] = c
        for i, y in enumerate(m):
            a[i + k] -= c * y
        a.pop()
        _trim(a)
    return _trim(q), a


class _Field:
    def __init__(self, minpoly):
        self.m = minpoly

    def red(self, a):
        return _pdivmod(a, self.m)[1] if len(a) >= len(self.m) else _trim(list(a))

    def mul(self, a, b):
        return self.red(_pmul(a, b))

    def sub(self, a, b):
        return _psub(a, b)

    def inv(self, a):
        # extended Euclid: s*a = 1 mod m
        r0, r1 = list(self.m), list(a)
        s0, s1 = [], [Fraction(1)]
        while r1:
            q, r = _pdivmod(r0, r1)
            r0, r1 = r1, r
            s0, s1 = s1, _psub(s0, _pmul(q, s1))
        # r0 is a nonzero constant since m is irreducible
        c = r0[0]
        return self.red([x / c for x in s0])


def _nullspace(F, A):
    """Nullspace of ``A`` from its reduced row echelon form (unique), as
    ``DomainMatrix.nullspace(divide_last=True)``."""
    rows = [list(r) for r in A]
    nr = len(rows)
    nc = len(rows[0]) if nr else 0
    pivots = []
    i = 0
    for j in range(nc):
        if i >= nr:
            break
        p = next((k for k in range(i, nr) if rows[k][j]), None)
        if p is None:
            continue
        rows[i], rows[p] = rows[p], rows[i]
        iv = F.inv(rows[i][j])
        rows[i] = [F.mul(iv, x) for x in rows[i]]
        for k in range(nr):
            if k != i and rows[k][j]:
                c = rows[k][j]
                rows[k] = [F.sub(x, F.mul(c, y)) for x, y in zip(rows[k], rows[i])]
        pivots.append(j)
        i += 1
    basis = []
    for j in range(nc):
        if j in pivots:
            continue
        v = [[] for _ in range(nc)]
        v[j] = [Fraction(1)]
        for r, pj in enumerate(pivots):
            v[pj] = F.sub([], rows[r][j])
        basis.append(v)
    return basis


def _frac(x):
    return Fraction(int(x.p), int(x.q))


def _to_rational(c):
    from . import Rational
    return Rational(c.numerator, c.denominator)


def eigenvects(M, simplify=True, primitive=False):
    """``_eigenvects`` for a rational matrix through ``_eigenvects_DOM``."""
    from types import FunctionType
    from . import (Dummy, Poly, QQ, roots, CRootOf, factor_list, gcd,
                   default_sort_key)
    from . import simplify as _simplify
    from . import Matrix
    n = M.rows
    A = [[_frac(M[i, j]) for j in range(n)] for i in range(n)]
    l = Dummy('lambda')
    charpoly = Poly(M.charpoly().as_expr(), M.charpoly().gens[0], domain=QQ)
    _, factors = factor_list(charpoly)
    result = []
    algebraic = []
    for base, mult in factors:
        coeffs = [_frac(c) for c in base.all_coeffs()]
        if len(coeffs) == 2:
            e = -coeffs[1] / coeffs[0]
            F = _Field([-e, Fraction(1)])
            AA = [[F.red([A[i][j]] if A[i][j] else []) for j in range(n)]
                  for i in range(n)]
            for i in range(n):
                AA[i][i] = F.sub(AA[i][i], [e])
            basis = _nullspace(F, AA)
            vects = [Matrix([_to_rational(x[0]) if x else 0 for x in v])
                     for v in basis]
            result.append((_to_rational(e), mult, vects))
        else:
            low = list(reversed(coeffs))
            F = _Field(low)
            AA = [[[A[i][j]] if A[i][j] else [] for j in range(n)]
                  for i in range(n)]
            for i in range(n):
                AA[i][i] = F.sub(AA[i][i], [Fraction(0), Fraction(1)])
            basis = _nullspace(F, AA)
            algebraic.append((low, mult, basis))
    for low, mult, basis in algebraic:
        minpoly = Poly([_to_rational(c) for c in reversed(low)], l, domain=QQ)
        vects = [[Poly([_to_rational(c) for c in reversed(x)] or [0], l,
                       domain=QQ).as_expr() for x in v] for v in basis]
        degree = minpoly.degree()
        mexpr = minpoly.as_expr()
        eigs = roots(mexpr, l)
        if len(eigs) != degree:
            eigs = [CRootOf(mexpr, l, idx) for idx in range(degree)]
        for ev in eigs:
            result.append((ev, mult, [Matrix([x.subs(l, ev) for x in v])
                                      for v in vects]))
    result = sorted(result, key=lambda x: default_sort_key(x[0]))
    if primitive:
        if isinstance(simplify, FunctionType):
            simpfunc = simplify
        else:
            simpfunc = _simplify if simplify else (lambda x: x)

        def denom_clean(vs):
            return [(v / gcd(list(v))).applyfunc(simpfunc) for v in vs]
        result = [(val, mult, denom_clean(es)) for val, mult, es in result]
    return result
