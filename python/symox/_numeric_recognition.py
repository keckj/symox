"""Rational conversion of native floats without a machine-precision detour."""
from fractions import Fraction
from decimal import Decimal
import math

from . import _core

def _decimal_fraction(value):
    # Decimal's exact string constructor and integer-ratio protocol avoid
    # Python's decimal integer-string limit without changing process globals.
    decimal=Decimal(str(value))
    # The Float itself has compact storage even for enormous exponents. An
    # explicit rational has exponent-sized integer storage. Preserve the old
    # native numeric-engine allocation range while removing the smaller
    # rationalization-specific precision and exponent restrictions.
    if abs(decimal.adjusted())>1000000:
        raise NotImplementedError('Exact rational output exceeds the current conversion resource budget')
    return Fraction(decimal)


def _exact_fraction(value):
    from . import RealField
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


def _limit_denominator(self, max_denominator=1000000):
    """Choose the nearest bounded-denominator fraction with exact arithmetic."""
    from fractions import Fraction
    from . import Rational
    if not self.is_Rational:
        raise AttributeError('limit_denominator is defined only for rational numbers')
    fraction = Fraction(self.p, self.q)
    return Rational(fraction.limit_denominator(Fraction(int(max_denominator))))


def nsimplify(expr, constants=(), tolerance=None, full=False, rational=None,
              rational_conversion='base10'):
    from . import sympify
    if rational_conversion not in ('base10','exact'):
        raise ValueError("rational_conversion must be 'base10' or 'exact'")
    if rational is not None and rational is not True and rational is not False:
        raise TypeError('rational must be True, False or None')
    tol=None if tolerance is None else _tolerance(tolerance)
    expression=sympify(expr)
    if getattr(expression,'is_Matrix',False):
        matrix_constants=tuple(constants)
        return expression.applyfunc(lambda entry: nsimplify(entry,constants=matrix_constants,
            tolerance=tolerance,full=full,rational=rational,rational_conversion=rational_conversion))
    if rational is True or expression.free_symbols:
        return _rationalize(expression,tol,rational_conversion)
    # Numeric recognition requires the reference evalf protocol. A container
    # without it must not become an unrelated held scalar through _symox_.
    if not hasattr(expression,'evalf'):
        raise AttributeError("'%s' object has no attribute 'evalf'" % type(expression).__name__)
    if tol is not None:
        raise NotImplementedError('Explicit tolerance for irrational recognition is not implemented')
    if full:
        raise NotImplementedError('Full irrational recognition search is not implemented')
    constants=tuple(constants)
    # The original recognition kernel handles IEEE floats. Use it only to
    # propose a candidate for an MP atom, and validate that candidate at the
    # original precision before changing any value. Double under/overflow
    # must never turn a nonzero native MP value into zero or infinity.
    from . import Float, Rational, Abs, S
    replacements={}
    for atom in expression.atoms(Float):
        approximate=float(atom)
        if not math.isfinite(approximate) or approximate==0 and atom.is_zero is not True:
            continue
        proposal=_core.nsimplify(_core.Float(approximate),constants=constants)
        if proposal.is_Float:
            continue
        digits=max(30,math.ceil(atom._prec*math.log10(2))+10)
        original=atom.evalf(digits)
        candidate=proposal.evalf(digits)
        tolerance_digits=max(1,min(15,round(atom._prec*math.log10(2)-1)))
        threshold=Rational(1,10**tolerance_digits)*Abs(original)
        if (Abs(candidate-original)<=threshold)==S.true:
            replacements[atom]=proposal
    expression=expression.xreplace(replacements) if replacements else expression
    return _core.nsimplify(expression,constants=constants)
