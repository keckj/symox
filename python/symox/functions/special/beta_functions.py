from ... import beta,betainc,betainc_regularized
__all__ = ['beta','betainc','betainc_regularized']


def betainc_mpmath_fix(a,b,x1,x2,reg=0):
    import mpmath
    return mpmath.mpf(0) if x1==x2 else mpmath.betainc(a,b,x1,x2,regularized=bool(reg))
__all__.append('betainc_mpmath_fix')
