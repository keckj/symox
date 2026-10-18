"""Three-dimensional isotropic oscillator in spherical coordinates."""
import symox as S


def E_nl(n,l,hw):return (2*S.sympify(n)+S.sympify(l)+S.Rational(3,2))*S.sympify(hw)


def R_nl(n,l,nu,r):
    n,l,nu,r = map(S.sympify,(n,l,nu,r))
    if any(v.is_number and (v<0 or v.is_integer is False) for v in (n,l)):raise ValueError('oscillator quantum numbers must be nonnegative integers')
    normalization = S.sqrt(2*S.factorial(n)*(2*nu)**(l+S.Rational(3,2))/S.gamma(n+l+S.Rational(3,2)))
    return normalization*r**l*S.exp(-nu*r*r)*S.assoc_laguerre(n,l+S.Rational(1,2),2*nu*r*r)


__all__ = ['E_nl','R_nl']
