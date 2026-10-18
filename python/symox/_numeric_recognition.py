"""Rational conversion of native floats without a machine-precision detour."""
from fractions import Fraction
import re

from . import _core

MAX_PRECISION=100000
MAX_DECIMAL_EXPONENT=100000


def _decimal_fraction(value):
    if value._prec>MAX_PRECISION:
        raise NotImplementedError('Float rationalization exceeds 100000 bits of precision')
    text=str(value)
    exponent=re.search(r'[eE]([+-]?\d+)$',text)
    if exponent and abs(int(exponent.group(1)))>MAX_DECIMAL_EXPONENT:
        raise NotImplementedError('Float rationalization exceeds decimal exponent magnitude 100000')
    return Fraction(text)


def _exact_fraction(value):
    from . import RealField
    _decimal_fraction(value)  # Check exponent/precision before allocating an exact denominator.
    domain=RealField(prec=value._prec)
    numerator,denominator=domain.to_rational(domain.convert(value),limit=False)
    return Fraction(numerator,denominator)


def _tolerance(value):
    from . import sympify
    value=sympify(value)
    if value.is_finite is not True or value.is_positive is not True or value.is_real is not True:
        raise ValueError('tolerance must be finite and positive')
    if value.is_Rational:return Fraction(int(value.p),int(value.q))
    if value.is_Float:return _exact_fraction(value)
    return _exact_fraction(value.evalf(30))


def _rationalize(expression, tolerance, conversion):
    from . import Float, Rational
    replacements={}
    for atom in expression.atoms(Float):
        if tolerance is not None:
            value=_exact_fraction(atom)
            if tolerance<1:
                bound=-(-tolerance.denominator//tolerance.numerator)
                result=value.limit_denominator(bound)
            else:
                rounded=tolerance*round(value/tolerance)
                result=rounded.limit_denominator(int(tolerance))
        elif conversion=='exact':
            result=_exact_fraction(atom)
        else:
            result=_decimal_fraction(atom)
            # Recognize small fractions at the displayed precision. Retain
            # decimal data if a nearby small denominator is not justified.
            candidate=result.limit_denominator(10000)
            digits=max(1,min(15,int(round(atom._prec*0.3010299956639812-1))))
            if ((candidate!=0 or atom.is_zero is True)
                    and abs(result-candidate)<=Fraction(1,10**digits)):
                result=candidate
        replacements[atom]=Rational(result.numerator,result.denominator)
    return expression.xreplace(replacements) if replacements else expression


def nsimplify(expr, constants=(), tolerance=None, full=False, rational=None,
              rational_conversion='base10'):
    from . import sympify
    if rational_conversion not in ('base10','exact'):
        raise ValueError("rational_conversion must be 'base10' or 'exact'")
    if rational is not None and rational is not True and rational is not False:
        raise TypeError('rational must be True, False or None')
    tol=None if tolerance is None else _tolerance(tolerance)
    expression=sympify(expr)
    if rational is True:
        return _rationalize(expression,tol,rational_conversion)
    if tol is not None:
        raise NotImplementedError('Explicit tolerance for irrational recognition is not implemented')
    if full:
        raise NotImplementedError('Full irrational recognition search is not implemented')
    return _core.nsimplify(expression,constants=constants)
