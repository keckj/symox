"""Exact Galois groups of irreducible rational polynomials of degree <= 4.

Cubic groups are determined by their discriminant. Quartics use the cubic
resolvent and, when needed, factorization over the discriminant quadratic
field. All polynomial operations and algebraic-field arithmetic are native.
"""
import operator


def _poly_galois_group(f, by_name=False, max_tries=30, randomize=False):
    from ... import Poly, QQ, ZZ, sqrt
    from ...combinatorics.galois import (
        S1TransitiveSubgroups as G1, S2TransitiveSubgroups as G2,
        S3TransitiveSubgroups as G3, S4TransitiveSubgroups as G4,
    )
    max_tries = operator.index(max_tries)
    if max_tries < 1:
        raise ValueError('max_tries must be positive')
    # The supported algorithm is deterministic; randomize is accepted for
    # signature compatibility and does not change its mathematical result.
    try:
        domain = f.domain
    except NotImplementedError:
        domain = None
    if not f.is_univariate or domain not in (ZZ, QQ) or not f.is_irreducible:
        raise ValueError('Polynomial must be irreducible and univariate over ZZ or QQ.')
    degree = int(f.degree())
    if degree < 1:
        raise ValueError('Constant polynomial has no Galois group.')
    if degree > 6:
        raise ValueError('Only polynomials up to degree 6 are supported.')
    if degree > 4:
        raise NotImplementedError('Galois groups of degree 5 and 6 are not implemented')
    if degree == 1:
        name, alternating = G1.S1, True
    elif degree == 2:
        name, alternating = G2.S2, False
    else:
        discriminant = f.discriminant()
        alternating = bool(sqrt(discriminant).is_Rational)
        if degree == 3:
            name = G3.A3 if alternating else G3.S3
        else:
            _, a, b, c, d = f.monic().all_coeffs()
            x = f.gens[0]
            resolvent = Poly(x**3 - b*x**2 + (a*c-4*d)*x + 4*b*d-a*a*d-c*c, x, domain=QQ)
            degrees = sorted(int(p.degree()) for p, multiplicity in resolvent.factor_list()[1]
                             for _ in range(multiplicity))
            if degrees == [3]:
                name = G4.A4 if alternating else G4.S4
            elif degrees == [1, 1, 1]:
                name = G4.V
            elif degrees == [1, 2]:
                # Over Q(sqrt(discriminant)), C4 has an intransitive even
                # subgroup C2, whereas D4 has the transitive subgroup V4.
                field = QQ.algebraic_field(sqrt(discriminant))
                factors = f.set_domain(field).factor_list()[1]
                name = G4.C4 if len(factors) > 1 else G4.D4
            else:
                raise ArithmeticError('Unexpected quartic resolvent factorization')
    return (name if by_name else name.get_perm_group()), alternating


def galois_group(f, *gens, by_name=False, max_tries=30, randomize=False, **args):
    """Return the Galois group and whether it lies in the alternating group.

    Irreducible univariate polynomials over ZZ and QQ through degree four are
    supported. Degrees five and six raise NotImplementedError.
    """
    from ... import Poly
    f = Poly(f, *gens, **args)
    return _poly_galois_group(f, by_name=by_name, max_tries=max_tries, randomize=randomize)
