"""Hermite normal forms of DomainMatrix over ZZ.

A port of ``_hermite_normal_form`` and ``_hermite_normal_form_modulo_D``
from sympy/polys/matrices/normalforms.py (Cohen, Algorithms 2.4.5 and
2.4.8). The entries are processed as exact Python integers.
"""
from collections import defaultdict

from ..domains import ZZ, QQ
from ..matrices.domainmatrix import DomainMatrix
from ..matrices.exceptions import DMDomainError, DMShapeError


def _igcdex(a, b):
    # The reference ``ZZ.gcdex`` for Python integers.
    if not a and not b:
        return 0, 1, 0
    if not a:
        return 0, b // abs(b), abs(b)
    if not b:
        return a // abs(a), 0, abs(a)
    x_sign = -1 if a < 0 else 1
    y_sign = -1 if b < 0 else 1
    a, b = abs(a), abs(b)
    x, y, r, s = 1, 0, 0, 1
    while b:
        c, q = a % b, a // b
        a, b, r, s, x, y = b, c, x - q * r, y - q * s, r, s
    return x * x_sign, y * y_sign, a


def _gcdex(a, b):
    x, y, g = _igcdex(a, b)
    if a != 0 and b % a == 0:
        y = 0
        x = -1 if a < 0 else 1
    return x, y, g


def _add_columns(m, i, j, a, b, c, d):
    # replace m[:, i] by a*m[:, i] + b*m[:, j]
    # and m[:, j] by c*m[:, i] + d*m[:, j]
    for k in range(len(m)):
        e = m[k][i]
        m[k][i] = a * e + b * m[k][j]
        m[k][j] = c * e + d * m[k][j]


def _symmetric_residue(a, m):
    return a - m if a > m // 2 else a


def _int_rows(A):
    return [[int(v) for v in row] for row in A.to_list()]


def _hermite_normal_form(A):
    if not A.domain.is_ZZ:
        raise DMDomainError('Matrix must be over domain ZZ.')
    m, n = A.shape
    A = _int_rows(A)
    k = n
    for i in range(m - 1, -1, -1):
        if k == 0:
            break
        k -= 1
        for j in range(k - 1, -1, -1):
            if A[i][j] != 0:
                u, v, d = _gcdex(A[i][k], A[i][j])
                r, s = A[i][k] // d, A[i][j] // d
                _add_columns(A, k, j, u, v, -s, r)
        b = A[i][k]
        if b < 0:
            _add_columns(A, k, k, -1, 0, -1, 0)
            b = -b
        if b == 0:
            k += 1
        else:
            for j in range(k + 1, n):
                q = A[i][j] // b
                _add_columns(A, j, k, 1, -q, 0, 1)
    rows = [row[k:] for row in A]
    return DomainMatrix([[ZZ(v) for v in row] for row in rows], (m, n - k), ZZ).to_dense()


def _hermite_normal_form_modulo_D(A, D):
    if not A.domain.is_ZZ:
        raise DMDomainError('Matrix must be over domain ZZ.')
    if isinstance(D, bool) or not (isinstance(D, int) or ZZ.of_type(D)) or D < 1:
        raise DMDomainError('Modulus D must be positive element of domain ZZ.')
    D = int(D)

    def add_columns_mod_R(m, R, i, j, a, b, c, d):
        for k in range(len(m)):
            e = m[k][i]
            m[k][i] = _symmetric_residue((a * e + b * m[k][j]) % R, R)
            m[k][j] = _symmetric_residue((c * e + d * m[k][j]) % R, R)

    W = defaultdict(dict)
    m, n = A.shape
    if n < m:
        raise DMShapeError('Matrix must have at least as many columns as rows.')
    A = _int_rows(A)
    k = n
    R = D
    for i in range(m - 1, -1, -1):
        k -= 1
        for j in range(k - 1, -1, -1):
            if A[i][j] != 0:
                u, v, d = _gcdex(A[i][k], A[i][j])
                r, s = A[i][k] // d, A[i][j] // d
                add_columns_mod_R(A, R, k, j, u, v, -s, r)
        b = A[i][k]
        if b == 0:
            A[i][k] = b = R
        u, v, d = _gcdex(b, R)
        for ii in range(m):
            W[ii][i] = u * A[ii][k] % R
        if W[i][i] == 0:
            W[i][i] = R
        for j in range(i + 1, m):
            q = W[i][j] // W[i][i]
            _add_columns(W, j, i, 1, -q, 0, 1)
        R //= d
    rows = [[ZZ(W[i].get(j, 0)) for j in range(m)] for i in range(m)]
    return DomainMatrix(rows, (m, m), ZZ).to_dense()


def hermite_normal_form(A, *, D=None, check_rank=False):
    if not A.domain.is_ZZ:
        raise DMDomainError('Matrix must be over domain ZZ.')
    if D is not None and (not check_rank or A.convert_to(QQ).rank() == A.shape[0]):
        return _hermite_normal_form_modulo_D(A, D)
    return _hermite_normal_form(A)
