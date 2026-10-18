"""Hydrogenic energies and normalized atomic-unit wavefunctions."""
import symox as S


def _validate(n,l=0,m=0):
    if n.is_number and (n.is_integer is False or n<1):raise ValueError('principal quantum number must be a positive integer')
    if l.is_number and (l.is_integer is False or l<0):raise ValueError('orbital quantum number must be a nonnegative integer')
    if n.is_number and l.is_number and l>=n:raise ValueError('orbital quantum number must be less than principal quantum number')
    if m.is_number and (m.is_integer is False or l.is_number and abs(m)>l):raise ValueError('magnetic quantum number is outside its range')


def E_nl(n,Z=1):
    n,Z = map(S.sympify,(n,Z));_validate(n,S.Integer(0),S.Integer(0));return -Z**2/(2*n*n)


def R_nl(n,l,r,Z=1):
    n,l,r,Z = map(S.sympify,(n,l,r,Z));_validate(n,l,S.Integer(0))
    rho = 2*Z*r/n
    normalization = S.sqrt((2*Z/n)**3*S.factorial(n-l-1)/(2*n*S.factorial(n+l)))
    return normalization*rho**l*S.assoc_laguerre(n-l-1,2*l+1,rho)*S.exp(-rho/2)


def Psi_nlm(n,l,m,r,phi,theta,Z=1):
    n,l,m,r,phi,theta,Z = map(S.sympify,(n,l,m,r,phi,theta,Z));_validate(n,l,m)
    return R_nl(n,l,r,Z)*S.expand_func(S.Ynm(l,m,theta,phi))


def E_nl_dirac(n,l,spin_up=True,Z=1,c=S.Float('137.035999037')):
    n,l,Z,c = map(S.sympify,(n,l,Z,c));_validate(n,l,S.Integer(0))
    if l==0 and not spin_up:raise ValueError('l=0 requires spin_up=True')
    kappa = -(l+1) if spin_up else -l
    effective = n+kappa+S.sqrt(kappa*kappa-Z*Z/(c*c))
    return c*c*(1/S.sqrt(1+Z*Z/(c*c*effective*effective))-1)


__all__ = ['E_nl','E_nl_dirac','Psi_nlm','R_nl']
