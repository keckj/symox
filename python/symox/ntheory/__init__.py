"""Number theory (``sympy.ntheory``).

Heavy algorithms (sieving, primality, factoring, residues, discrete
logarithms, partitions) run in the native core; these modules reproduce
SymPy's argument handling, return types and orderings.
"""

from .generate import (nextprime, prevprime, prime, primepi, primerange,
                       randprime, Sieve, sieve, primorial, cycle_length,
                       composite, compositepi)
from .primetest import (isprime, is_gaussian_prime, is_mersenne_prime,
                        is_fermat_pseudoprime, is_euler_pseudoprime,
                        is_euler_jacobi_pseudoprime, mr, is_lucas_prp,
                        is_strong_lucas_prp, is_extra_strong_lucas_prp,
                        proth_test)
from .factor_ import (divisors, proper_divisors, factorint, multiplicity,
                      perfect_power, pollard_pm1, factor_cache, pollard_rho,
                      primefactors, totient, divisor_count,
                      proper_divisor_count, divisor_sigma, factorrat,
                      reduced_totient, primenu, primeomega,
                      mersenne_prime_exponent, is_perfect, is_abundant,
                      is_deficient, is_amicable, is_carmichael, abundance,
                      dra, drm, multiplicity_in_factorial)
from .partitions_ import npartitions
from .residue_ntheory import (is_primitive_root, is_quad_residue,
                              legendre_symbol, jacobi_symbol, n_order,
                              sqrt_mod, quadratic_residues, primitive_root,
                              nthroot_mod, is_nthpow_residue, sqrt_mod_iter,
                              mobius, discrete_log, quadratic_congruence,
                              polynomial_congruence, binomial_mod)
from .multinomial import (binomial_coefficients, binomial_coefficients_list,
                          multinomial_coefficients)
from .continued_fraction import (continued_fraction_periodic,
                                 continued_fraction_iterator,
                                 continued_fraction_reduce,
                                 continued_fraction_convergents,
                                 continued_fraction)
from .digits import digits, count_digits, is_palindromic
from .egyptian_fraction import egyptian_fraction
from .ecm import ecm
from .qs import qs, qs_factor
from .bbp_pi import pi_hex_digits
from .modular import crt, crt1, crt2, solve_congruence, symmetric_residue
from ._intfunc import (igcd, ilcm, igcdex, integer_log, num_digits,
                       trailing, mod_inverse, integer_nthroot, isqrt)

__all__ = [
    'nextprime', 'prevprime', 'prime', 'primepi', 'primerange', 'randprime',
    'Sieve', 'sieve', 'primorial', 'cycle_length', 'composite', 'compositepi',
    'isprime', 'is_gaussian_prime', 'is_mersenne_prime',
    'divisors', 'proper_divisors', 'factorint', 'multiplicity',
    'perfect_power', 'pollard_pm1', 'factor_cache', 'pollard_rho',
    'primefactors', 'totient', 'divisor_count', 'proper_divisor_count',
    'divisor_sigma', 'factorrat', 'reduced_totient', 'primenu', 'primeomega',
    'mersenne_prime_exponent', 'is_perfect', 'is_abundant', 'is_deficient',
    'is_amicable', 'is_carmichael', 'abundance', 'dra', 'drm',
    'multiplicity_in_factorial',
    'npartitions',
    'is_primitive_root', 'is_quad_residue', 'legendre_symbol',
    'jacobi_symbol', 'n_order', 'sqrt_mod', 'quadratic_residues',
    'primitive_root', 'nthroot_mod', 'is_nthpow_residue', 'sqrt_mod_iter',
    'mobius', 'discrete_log', 'quadratic_congruence', 'polynomial_congruence',
    'binomial_coefficients', 'binomial_coefficients_list',
    'multinomial_coefficients',
    'continued_fraction_periodic', 'continued_fraction_iterator',
    'continued_fraction_reduce', 'continued_fraction_convergents',
    'continued_fraction',
    'digits', 'count_digits', 'is_palindromic',
    'egyptian_fraction', 'ecm', 'qs', 'qs_factor',
]
