"""SymPy 1.14 parity of the ntheory submodule algorithms (probable prime
tests, ECM, quadratic sieve, BBP, polynomial congruences, binomial_mod,
Egyptian fraction helpers, smoothness_p) and of sympy.discrete
(transforms, convolutions, recurrences)."""
import pytest
import sympy as sp
import symox as sx

from _ntheory_compare import NS, check, run

for _lib, _ns in NS.items():
    exec(
        'from {0}.ntheory.primetest import (is_fermat_pseudoprime, '
        'is_euler_pseudoprime, is_euler_jacobi_pseudoprime, mr, '
        'is_lucas_prp, is_strong_lucas_prp, is_extra_strong_lucas_prp, '
        'proth_test, is_square)\n'
        'from {0}.ntheory.residue_ntheory import (polynomial_congruence, '
        'binomial_mod)\n'
        'from {0}.ntheory.factor_ import smoothness_p\n'
        'from {0}.ntheory.bbp_pi import pi_hex_digits\n'
        'from {0}.ntheory.egyptian_fraction import (egypt_greedy, '
        'egypt_graham_jewett, egypt_takenouchi, egypt_golomb, '
        'egypt_harmonic)\n'
        'from {0}.ntheory.ecm import ecm, Point\n'
        'from {0}.ntheory.qs import qs, qs_factor, SievePolynomial\n'
        'from {0}.discrete.recurrences import linrec, linrec_coeffs\n'
        'from {0}.discrete.convolutions import (convolution_fft, '
        'convolution_ntt, convolution_fwht, convolution_subset, '
        'convolution_int)\n'.format(_lib), _ns)

PRIMETEST = r"""
[n for n in range(2, 2000) if is_fermat_pseudoprime(n, 2) and not isprime(n)]
[n for n in range(2, 3000) if is_euler_pseudoprime(n, 2) and not isprime(n)]
[n for n in range(3, 3000, 2) if is_euler_jacobi_pseudoprime(n, 3) and not isprime(n)]
is_fermat_pseudoprime(2, 1)
is_fermat_pseudoprime(9, 1)
is_fermat_pseudoprime(8, 1)
is_fermat_pseudoprime(1, 2)
is_fermat_pseudoprime(10, 3)
is_fermat_pseudoprime(2, 2)
is_fermat_pseudoprime(15, 3)
is_fermat_pseudoprime(0, 2)
is_fermat_pseudoprime(7, 0)
is_fermat_pseudoprime(7.5, 2)
is_fermat_pseudoprime(2**89 - 1, 3)
is_euler_pseudoprime(1, 2)
is_euler_pseudoprime(2, 1)
is_euler_pseudoprime(15, 1)
is_euler_pseudoprime(10, 3)
is_euler_pseudoprime(15, 3)
is_euler_pseudoprime(0, 3)
is_euler_pseudoprime(7, 0)
is_euler_pseudoprime(341, 2)
is_euler_jacobi_pseudoprime(561, 2)
is_euler_jacobi_pseudoprime(1, 5)
is_euler_jacobi_pseudoprime(15, 5)
is_euler_jacobi_pseudoprime(15, 1)
is_euler_jacobi_pseudoprime(4, 3)
is_euler_jacobi_pseudoprime(2**127 - 1, 7)
mr(1373653, [2, 3])
mr(1373653, [5])
mr(25326001, [2, 3, 5])
mr(3215031751, [2, 3, 5, 7])
mr(2, [2])
mr(1, [2])
mr(4, [3])
mr(17, [17, 34, 1, 0])
mr(2**61 - 1, [2, 3, 5, 7, 11, 13])
[n for n in range(2, 20000) if is_lucas_prp(n) and not isprime(n)]
[n for n in range(2, 20000) if is_strong_lucas_prp(n) and not isprime(n)]
[n for n in range(2, 20000) if is_extra_strong_lucas_prp(n) and not isprime(n)]
is_lucas_prp(1)
is_lucas_prp(2)
is_lucas_prp(4)
is_lucas_prp(25)
is_lucas_prp(2**89 - 1)
is_strong_lucas_prp(-7)
is_strong_lucas_prp(49)
is_strong_lucas_prp(2**107 - 1)
is_extra_strong_lucas_prp(2)
is_extra_strong_lucas_prp(9)
is_extra_strong_lucas_prp(2**127 - 1)
is_extra_strong_lucas_prp(2**127 + 1)
[n for n in range(3, 2000) if not (n - 1) % 2 and (n - 1) >> ((n - 1) & -(n - 1)).bit_length() - 1 < (n - 1) & -(n - 1) and proth_test(n)]
proth_test(97)
proth_test(3)
proth_test(5)
proth_test(7)
proth_test(2)
proth_test(13)
proth_test(25)
proth_test(57)
proth_test(81)
proth_test(3 * 2**20 + 1)
proth_test(5 * 2**100 + 1)
proth_test(9 * 2**50 + 1)
proth_test(Integer(97))
"""

FACTOR = r"""
ecm(398883434337287)
ecm(46167045131415113)
ecm(64211816600515193)
ecm(168541512131094651323)
ecm(631211032315670776841)
ecm(4132846513818654136451)
ecm(4516511326451341281684513)
ecm(3146531246531241245132451321)
ecm(4269021180054189416198169786894227)
ecm(7060005655815754299976961394452809, B1=100000, B2=1000000)
ecm(7060005655815754299976961394452809, B1=100000, B2=1000000, seed=3)
ecm(2**61 - 1)
ecm(2**4 * 3**5 * 101)
ecm(1009**2 * 1013**3 * 100003**2 * 1000003)
ecm(1)
ecm(6)
ecm(7**10)
ecm(-15)
ecm(10, B1=11)
ecm(10, B2=11)
ecm(2.0)
ecm(Integer(398883434337287))
sorted(qs_factor(10009202107, 100, 10000).items())
sorted(qs_factor(211107295182713951054568361, 1000, 10000).items())
sorted(qs_factor(5915587277*3267000013, 1000, 10000).items())
sorted(qs_factor(5915587277*3267000013*1009*1013, 1000, 10000).items())
sorted(qs_factor(5915587277*3267000013*1009, 1000, 10000).items())
sorted(qs_factor(2**5 * 5915587277*3267000013, 1000, 10000).items())
sorted(qs_factor(15, 100, 100).items())
sorted(qs_factor(10403, 100, 100).items())
sorted(qs_factor(7**4, 100, 100).items())
sorted(qs_factor(101, 100, 100).items())
qs(5915587277*3267000013, 1000, 10000)
qs(10009202107, 100, 10000, seed=7)
qs_factor(1, 10, 10)
qs_factor(-5, 10, 10)
smoothness_p(10431, m=1)
smoothness_p(10431)
smoothness_p(10431, power=1)
smoothness_p(Integer(10431))
smoothness_p(21477639576571, visual=1)
smoothness_p(21477639576571, visual=0)
smoothness_p(2**3 * 3**2 * 7, m=2, power=1, visual=True)
smoothness_p(smoothness_p(10431, visual=True))
smoothness_p(smoothness_p(10431, visual=True), visual=False)
smoothness_p(smoothness_p(10431, visual=True), visual=True)
smoothness_p(smoothness_p(10431, m=1))
smoothness_p(Mul(3, 7, evaluate=False) * 11)
smoothness_p(1)
"""

RESIDUE = r"""
polynomial_congruence(x**7 - x + 3, 7)
polynomial_congruence(x**5 + 1, 2**3 * 3**2)
polynomial_congruence(x**3 - x**2 + 2*x + 3, 24)
polynomial_congruence(x**2 + 1, 5**4 * 13)
polynomial_congruence(3*x - 2, 7)
polynomial_congruence(6*x - 4, 10)
polynomial_congruence(2*x**2 + 3*x + 1, 15)
polynomial_congruence(x**3 + 1, 9)
polynomial_congruence(x**4 - 16, 1001)
polynomial_congruence(x**8 + x**3 + 7, 5**3 * 7**2 * 11)
polynomial_congruence(x**5 + 2*x**4 + 3, 1000003)
polynomial_congruence(x**3 + 5*x + 1, 10007 * 10009)
polynomial_congruence(x**7 + 3*x + 11, 2**10)
polynomial_congruence(x**4 + x + 1, 2**5 * 3**4 * 5**2)
polynomial_congruence(-x**4 + 2*x + 7, 1)
polynomial_congruence(x**3, 64)
polynomial_congruence(x**5 - x, 30)
polynomial_congruence(x**4 + 1, 2**61 - 1)
polynomial_congruence(x/2 + 1, 7)
polynomial_congruence(x*y + 1, 7)
polynomial_congruence(sin(x) + 1, 7)
polynomial_congruence(1/x + 1, 7)
binomial_mod(10, 3, 7)
binomial_mod(1000, 500, 1000000007)
binomial_mod(10**18, 10**9, 3**5 * 7**2 * 11)
binomial_mod(123456789, 98765, 2**10 * 5**4 * 13**3 * 101)
binomial_mod(52, 5, 2**20)
binomial_mod(10**30, 10**15 + 7, 9973**2 * 2**7)
[binomial_mod(n, m, 360) for n in range(15) for m in range(n + 2)]
binomial_mod(-3, 1, 7)
binomial_mod(5, -1, 7)
binomial_mod(3, 5, 7)
binomial_mod(10, 3, 0)
binomial_mod(10, 3, 1)
binomial_mod(10, 3, 13)
"""

MISC = r"""
pi_hex_digits(0)
pi_hex_digits(1)
pi_hex_digits(10)
pi_hex_digits(91, 10)
pi_hex_digits(1000, 30)
pi_hex_digits(0, 0)
pi_hex_digits(5, 1)
pi_hex_digits(-1)
pi_hex_digits(1, -1)
pi_hex_digits(1.5)
pi_hex_digits(Integer(3), Integer(5))
egypt_greedy(4, 13)
egypt_greedy(3, 7)
egypt_greedy(1, 9)
egypt_greedy(17, 89)
egypt_graham_jewett(4, 13)
egypt_graham_jewett(3, 7)
egypt_takenouchi(4, 13)
egypt_takenouchi(3, 7)
egypt_takenouchi(3, 8)
egypt_golomb(4, 13)
egypt_golomb(17, 89)
egypt_harmonic(Rational(5, 2))
egypt_harmonic(Rational(3, 2))
egypt_harmonic(Rational(-1, 2))
egypt_harmonic(Rational(1, 3))
egypt_harmonic(Integer(2))
"""

DISCRETE = r"""
ntt([1, 2, 3, 4], 3*2**8 + 1)
intt(ntt([1, 2, 3, 4], 3*2**8 + 1), 3*2**8 + 1)
ntt([1, 2, 3, 4, 5], 3*2**8 + 1)
ntt([], 7)
ntt([5], 7)
ntt([-1, 10**30, 7], 998244353)
intt([1, 2, 3], 998244353)
ntt([1, 2, 3], 9)
ntt([1, 2, 3, 4, 5], 7)
ntt([Integer(3), 4], 5)
ntt([1.5, 2], 5)
ntt(3, 5)
fwht([4, 2, 2, 0, 0, 2, -2, 0])
ifwht([8, 0, 8, 0, 8, 8, 0, 0])
fwht([1, 2, 3])
ifwht([1, 2, 3])
fwht([Rational(1, 2), 3])
fwht([x, y, 1])
fwht([])
fwht([7])
fwht(5)
mobius_transform([1, 2, 3])
mobius_transform([1, 2, 3], subset=False)
inverse_mobius_transform([1, 3, 4, 6])
inverse_mobius_transform([6, 2, 3, 0], subset=False)
mobius_transform([x, y, 2, Rational(1, 3)])
mobius_transform([10**40, -3, 5, 7, 11])
inverse_mobius_transform([1])
convolution([1, 2, 3], [4, 5, 6])
convolution([1, 2, 3], [4, 5, 6], cycle=2)
convolution([1, 2, 3], [4, 5, 6], prime=19*2**10 + 1)
convolution([1, 2, 3], [4, 5, 6], prime=19*2**10 + 1, cycle=3)
convolution([1, 2, 3], [4, 5, 6], dyadic=True)
convolution([1, 2, 3], [4, 5, 6], subset=True)
convolution([Rational(1, 2), 2], [Rational(2, 3), 5])
convolution([x, 2], [y, 1])
convolution([1, 2], [3, 4], dps=10)
convolution([1, 2], [3, 4], prime=7, dyadic=True)
convolution([1, 2], [3, 4], cycle=-1)
convolution([], [1])
convolution_fft([1, 2], [3, 4])
convolution_ntt([1, 2], [3, 4], 17)
convolution_fwht([1, 2], [3, 4, 5])
convolution_subset([1, 2], [3, 4, 5])
convolution_int([1, -2, 3], [-4, 5])
covering_product([1, 2, 3], [4, 5, 6])
intersecting_product([1, 2, 3], [4, 5, 6])
covering_product([x, 1], [y, 2, 3])
intersecting_product([], [1])
linrec([1, 1], [0, 1], 10)
linrec([1, 1], [0, 1], 1000)
linrec([2, -1, 3], [1, 0, 2], 7)
linrec([2, -1, 3], [1, 0, 2], 0)
linrec([2, -1, 3], [1], 5)
linrec([1, 1], [x, y], 10)
linrec([x, 1], [1, 1], 4)
linrec([Rational(1, 2), 1], [1, 1], 6)
linrec([], [1], 5)
linrec([1, 1], [1, 2, 3], 5)
linrec([1, 1], [0, 1], -1)
linrec([1, 1], [0, 1], 2.5)
linrec(1, [0, 1], 3)
linrec([1, 1], 1, 3)
linrec([0, 0, 1], [1, 2, 3], 10**6)
linrec_coeffs([1, 1], 10)
linrec_coeffs([1, 1], 1)
linrec_coeffs([2, -1, 3], 9)
linrec_coeffs([x, 1], 4)
"""


def _cases(block):
    return [s for s in block.strip().splitlines() if s.strip()]


@pytest.mark.parametrize('src', _cases(PRIMETEST))
def test_primetest(src):
    check(src)


@pytest.mark.parametrize('src', _cases(FACTOR))
def test_factor(src):
    check(src)


@pytest.mark.parametrize('src', _cases(RESIDUE))
def test_residue(src):
    check(src)


@pytest.mark.parametrize('src', _cases(MISC))
def test_misc(src):
    check(src)


@pytest.mark.parametrize('src', _cases(DISCRETE))
def test_discrete(src):
    check(src)


def test_ecm_point_arithmetic():
    for lib in (sp, sx):
        Point = lib.ntheory.ecm.__class__  # noqa: F841
    from sympy.ntheory.ecm import Point as SP
    from symox.ntheory.ecm import Point as XP
    for (x, z, a24, m) in [(11, 16, 7, 29), (3, 5, 2, 101),
                           (10**20 + 3, 7, 12345, 10**21 + 117)]:
        sp_p, xp_p = SP(x, z, a24, m), XP(x, z, a24, m)
        sp_q, xp_q = sp_p.double(), xp_p.double()
        assert (sp_q.x_cord, sp_q.z_cord) == (xp_q.x_cord, xp_q.z_cord)
        sp_r, xp_r = sp_q.add(sp_p, sp_p), xp_q.add(xp_p, xp_p)
        assert (sp_r.x_cord, sp_r.z_cord) == (xp_r.x_cord, xp_r.z_cord)
        for k in (1, 2, 3, 10, 97, 2**70 + 5):
            a, b = sp_p.mont_ladder(k), xp_p.mont_ladder(k)
            assert (a.x_cord, a.z_cord, a.a_24, a.mod) == \
                (b.x_cord, b.z_cord, b.a_24, b.mod)
        assert xp_p.mont_ladder(3) == xp_r
        assert (sp_p.mont_ladder(3) == sp_r) is True
        assert (xp_p == XP(x, z, a24 + 1, m)) is False


def test_qs_helpers():
    from sympy.ntheory.qs import SievePolynomial as SPoly
    from symox.ntheory.qs import SievePolynomial as XPoly, FactorBaseElem
    a, b = SPoly(7, 3, 101), XPoly(7, 3, 101)
    for t in (-3, 0, 5):
        assert (a.eval_u(t), a.eval_v(t)) == (b.eval_u(t), b.eval_v(t))
    f = FactorBaseElem(7, 3, 1993)
    assert (f.prime, f.tmem_p, f.log_p, f.soln1) == (7, 3, 1993, None)


def test_ecm_qs_results_are_sets_and_dicts():
    assert type(sx.ntheory.ecm(398883434337287)) is set
    assert type(sx.ntheory.qs(10403, 100, 100)) is set
    assert type(sx.ntheory.qs_factor(10403, 100, 100)) is dict
    assert run('symox', 'ecm(3*5*7*1000003*1000033)') == \
        run('sympy', 'ecm(3*5*7*1000003*1000033)')


def test_linrec_large_native():
    import sympy.discrete.recurrences as spr
    import symox.discrete.recurrences as sxr
    c, b = [3, -1, 4, 1, -5], [2, 7, 1, 8, 2]
    for n in (5, 77, 1234):
        assert int(sxr.linrec(c, b, n)) == int(spr.linrec(c, b, n))
    assert sxr.linrec(c, b, 77).is_Integer is True
