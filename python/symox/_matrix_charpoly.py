"""Characteristic-polynomial coefficient callbacks over native domains."""


def _reference_route(M):
    """None when the native kernel already follows SymPy's ``_charpoly``,
    else ``('EXRAW', dM)`` or ``('inexact', dM)`` for the ``to_DM`` domain."""
    if not M.is_square or not M.rows or all(
            getattr(e, 'is_Rational', False) for e in M):
        return None
    from . import _core
    if any(isinstance(e, _core.Expr) and not e.is_commutative for e in M):
        return None
    dM = M.to_DM()
    K = dM.domain
    if K.is_EXRAW:
        return ('EXRAW', dM)
    if not K.is_Exact and not (K.is_RR or K.is_CC):
        return ('inexact', dM)
    return None


_native = None


def install(matrix):
    global _native
    native = _native = matrix.charpoly

    def _default(self, x):
        route = _reference_route(self)
        if route is None:
            return native(self, x)
        if route[0] == 'EXRAW':
            # SymPy's ``_charpoly`` default ``simplify=_simplify``
            # runs ``simplify`` on each EXRAW Berkowitz coefficient:
            # ``[[sin(x)**2 + cos(x)**2, 1], [1, 1]]`` gives
            # ``_x**2 - 2*_x``.
            from . import simplify as _simplify
            return charpoly(self, x, _simplify)
        # ``PurePoly(cp, x, domain=K)`` over ``RR[a]``/``CC[a]``.
        # The native Berkowitz coefficients are the same values; only
        # the domain of the native polynomial is the bare ground field.
        from . import PurePoly
        poly = native(self, x)
        return PurePoly(poly.all_coeffs(), poly.gen, domain=route[1].domain)

    def charpoly(self, x=None, simplify=None):
        if simplify is None:
            from . import evaluate
            with evaluate(True):
                return _default(self, x)
        from . import Add, PurePoly, construct_domain, sympify
        from ._core import Expr
        if not callable(simplify) or isinstance(simplify, Expr) and not simplify.is_Lambda:
            raise TypeError('simplify must be callable or None')
        # The native constructor handles generator identity and coefficient
        # collisions. DomainMatrix retains unsimplified EXRAW coefficients
        # until the caller has seen them, as required by this callback API.
        generator = native(self, x).gen
        domain_matrix = self.to_DM()
        coefficients = [sympify(simplify(domain_matrix.domain.to_sympy(value)))
                        for value in domain_matrix.charpoly()]
        if not all(isinstance(value,Expr) for value in coefficients):
            raise TypeError('simplify must return scalar expressions')
        while len(coefficients) > 1 and coefficients[0] == 0:
            coefficients.pop(0)
        domain, _ = construct_domain(coefficients)
        degree = len(coefficients)-1
        expression = Add(*(coefficient*generator**(degree-index)
                           for index, coefficient in enumerate(coefficients)))
        return PurePoly._from_coeffs(coefficients, generator, str(domain), expression)

    matrix.charpoly = charpoly
