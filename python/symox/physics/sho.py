"""Three-dimensional isotropic oscillator in spherical coordinates."""
import symox as S


def E_nl(n,l,hw):return (2*S.sympify(n)+S.sympify(l)+S.Rational(3,2))*S.sympify(hw)


def R_nl(n,l,nu,r):
    n,l,nu,r = map(S.sympify,(n,l,nu,r))
    if any(v.is_number and (v<0 or v.is_integer is False) for v in (n,l)):raise ValueError('oscillator quantum numbers must be nonnegative integers')
    n = n+1
    C = S.sqrt((2*nu)**(l+S.Rational(3,2))*2**(n+l+1)*S.factorial(n-1)/(S.sqrt(S.pi)*S.factorial2(2*n+2*l-1)))
    return C*r**l*S.exp(-nu*r**2)*S.assoc_laguerre(n-1,l+S.S.Half,2*nu*r**2)


__all__ = ['E_nl','R_nl']
