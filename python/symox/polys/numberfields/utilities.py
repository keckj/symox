"""Certified real algebraic root isolation using native Sturm intervals."""


def isolate(alg,eps=None,fast=False):
    from ... import sympify,minpoly
    alg=sympify(alg)
    if getattr(alg,"is_AlgebraicNumber",False):alg=alg.as_expr()
    if eps is not None:
        eps=sympify(eps)
        if not eps.is_Rational or eps<=0:raise ValueError('eps must be a positive rational')
    if alg.is_Rational:return alg,alg
    if alg.is_real is not True:raise NotImplementedError('isolation requires a real algebraic number')
    polynomial=minpoly(alg,polys=True)
    for (a,b),multiplicity in polynomial.intervals(eps=eps):
        if (alg-a).is_positive is True and (b-alg).is_positive is True:
            return a,b
    raise NotImplementedError('the algebraic embedding could not be certified against the isolating intervals')
