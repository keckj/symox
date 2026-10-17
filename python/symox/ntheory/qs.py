"""Self-initializing quadratic sieve (``sympy.ntheory.qs``)."""

from .. import _core


class SievePolynomial:
    """Sieve polynomial ``(a*x + b)**2 - N`` with ``u = a*x + b``."""

    def __init__(self, a, b, N):
        self.a = a
        self.b = b
        self.a2 = a**2
        self.ab = 2*a*b
        self.b2 = b**2 - N

    def eval_u(self, x):
        return self.a*x + self.b

    def eval_v(self, x):
        return (self.a2*x + self.ab)*x + self.b2


class FactorBaseElem:
    """Factor base element of the quadratic sieve."""

    def __init__(self, prime, tmem_p, log_p):
        self.prime = prime
        self.tmem_p = tmem_p
        self.log_p = log_p
        self.soln1 = None
        self.soln2 = None
        self.b_ainv = None


def qs(N, prime_bound, M, ERROR_TERM=25, seed=1234):
    """Set of factors of N found by the quadratic sieve."""
    return set(qs_factor(N, prime_bound, M, ERROR_TERM, seed))


def qs_factor(N, prime_bound, M, ERROR_TERM=25, seed=1234):
    """Factorization dict of N found by the quadratic sieve."""
    if N < 2:
        raise ValueError("N should be greater than 1")
    return dict(_core.nt_qs_factor(int(N), int(prime_bound), int(M),
                                   float(ERROR_TERM), int(seed)))
