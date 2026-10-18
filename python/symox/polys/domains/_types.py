"""Constructor facades for native Domain descriptors and coefficients."""
import operator

from ... import _core


class _DomainType(type):
    def __instancecheck__(cls, value):
        return isinstance(value, _core.Domain) and bool(getattr(value, cls._flag, False))


class IntegerRing(metaclass=_DomainType):
    _flag = 'is_ZZ'

    def __new__(cls):
        return _core.ZZ


class RationalField(metaclass=_DomainType):
    _flag = 'is_QQ'

    def __new__(cls):
        return _core.QQ


class FiniteField(metaclass=_DomainType):
    _flag = 'is_FiniteField'

    def __new__(cls, mod, dom=None, symmetric=True):
        if dom not in (None, _core.ZZ):
            raise NotImplementedError('finite fields require the integer ground ring')
        return _core.GF(operator.index(mod), symmetric=symmetric)


def _numerical(prec, dps, tol, complex):
    if tol is not None:
        raise NotImplementedError('custom numerical domain tolerance is not implemented')
    if prec is not None:
        prec = operator.index(prec)
        if not 2 <= prec <= 40000:
            raise ValueError('bit precision must lie between 2 and 40000')
    if dps is not None:
        dps = operator.index(dps)
        if not 1 <= dps <= 10000:
            raise ValueError('decimal precision must lie between 1 and 10000')
    return _core._numeric_field(prec=prec, dps=dps, complex=complex)


class RealField(metaclass=_DomainType):
    _flag = 'is_RealField'

    def __new__(cls, prec=None, dps=None, tol=None):
        return _numerical(prec, dps, tol, False)


class ComplexField(metaclass=_DomainType):
    _flag = 'is_ComplexField'

    def __new__(cls, prec=None, dps=None, tol=None):
        return _numerical(prec, dps, tol, True)


class GaussianIntegerRing(metaclass=_DomainType):
    _flag = 'is_GaussianRing'

    def __new__(cls):
        return _core.ZZ_I


class GaussianRationalField(metaclass=_DomainType):
    _flag = 'is_GaussianField'

    def __new__(cls):
        return _core.QQ_I


class ExpressionDomain(metaclass=_DomainType):
    _flag = 'is_EX'

    def __new__(cls):
        return _core.EX


class ExpressionRawDomain(metaclass=_DomainType):
    _flag = 'is_EXRAW'

    def __new__(cls):
        return _core.EXRAW


def _generators(gens):
    from ... import symbols
    if len(gens) == 1 and isinstance(gens[0], (tuple, list)):
        return tuple(gens[0])
    if len(gens) == 1 and isinstance(gens[0], str):
        parsed = symbols(gens[0])
        return parsed if isinstance(parsed, tuple) else (parsed,)
    return gens


class PolynomialRing(metaclass=_DomainType):
    _flag = 'is_PolynomialRing'

    def __new__(cls, dom, *gens, order='lex'):
        return dom.poly_ring(*_generators(gens), order=str(order))


class FractionField(metaclass=_DomainType):
    _flag = 'is_FractionField'

    def __new__(cls, dom, *gens, order='lex'):
        return dom.frac_field(*_generators(gens), order=str(order))


class AlgebraicField(metaclass=_DomainType):
    _flag = 'is_AlgebraicField'

    def __new__(cls, dom, *extensions, alias=None):
        return dom.algebraic_field(*extensions, alias=alias)


class GlobalPolynomialRing(PolynomialRing):
    """Legacy construction syntax using the native polynomial representation."""


class GeneralizedPolynomialRing(PolynomialRing):
    """Supported global orders use the same native polynomial representation."""


__all__ = ['IntegerRing', 'RationalField', 'FiniteField', 'RealField', 'ComplexField',
           'GaussianIntegerRing', 'GaussianRationalField', 'ExpressionDomain',
           'ExpressionRawDomain', 'PolynomialRing', 'FractionField', 'AlgebraicField',
           'GlobalPolynomialRing', 'GeneralizedPolynomialRing']
