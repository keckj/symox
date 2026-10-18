"""One-dimensional quantum oscillator eigenstates and energies."""
import symox as S
hbar = S.Symbol('hbar',positive=True)


def psi_n(n,x,m,omega):
    n,x,m,omega = map(S.sympify,(n,x,m,omega))
    if n.is_number and (n.is_integer is False or n<0):raise ValueError('oscillator level must be a nonnegative integer')
    nu = m*omega/hbar
    return (nu/S.pi)**S.Rational(1,4)*S.hermite(n,S.sqrt(nu)*x)*S.exp(-nu*x*x/2)/S.sqrt(2**n*S.factorial(n))


def E_n(n,omega):return hbar*S.sympify(omega)*(S.sympify(n)+S.Rational(1,2))

def coherent_state(n,alpha):
    n,alpha = map(S.sympify,(n,alpha));return S.exp(-S.Abs(alpha)**2/2)*alpha**n/S.sqrt(S.factorial(n))


__all__ = ['E_n','coherent_state','psi_n']
