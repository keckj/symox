"""Homogeneous symbolic MRV powers with explicit parameter-domain proofs."""
import pytest
import sympy as sp
import symox as sx


def issue(module, k, a=1, b=0, numerator=None):
    n = module.Symbol('n', positive=True)
    f = (n+a)**k / ((n+a)**(k+1) - (n+b)**(k+1))
    if numerator is not None:
        f *= numerator
    return f, n


@pytest.mark.parametrize('facts', ({'positive': True}, {'nonnegative': True}))
@pytest.mark.parametrize('a,b', ((1,0), (3,1), (-2,-5), (1,4)))
def test_positive_and_nonnegative_symbolic_power_shifted_family(facts,a,b):
    k = sx.Symbol('k', **facts)
    f,n = issue(sx,k,a,b)
    result = sx.limit(f,n,sx.oo)
    expected = 1 / ((a-b)*(k+1))
    assert not result.has(sx.Limit) and sx.cancel(result-expected) == 0
    k_ref = sp.Symbol('k', **facts)
    f_ref,n_ref = issue(sp,k_ref,a,b)
    assert sp.simplify(sp.limit(f_ref,n_ref,sp.oo)-1/((a-b)*(k_ref+1))) == 0


@pytest.mark.parametrize('scale', ('2','sqrt(2)','I','1+I'))
def test_certified_exact_real_and_complex_coefficients(scale):
    k = sx.Symbol('k', positive=True)
    coefficient = sx.sympify(scale)
    f,n = issue(sx,k,numerator=coefficient)
    result = sx.limit(f,n,sx.oo)
    assert not result.has(sx.Limit)
    assert sx.cancel(result-coefficient/(k+1)) == 0


@pytest.mark.parametrize('facts', ({}, {'real': True}, {'negative': True}, {'complex': True}))
@pytest.mark.parametrize('cancel_leading_factor', (False,True))
def test_unknown_parameter_pole_is_not_hidden_by_leading_coefficient_cancellation(facts,cancel_leading_factor):
    k = sx.Symbol('k', **facts)
    f,n = issue(sx,k,numerator=(k+1) if cancel_leading_factor else None)
    result = sx.limit(f,n,sx.oo)
    assert isinstance(result,sx.Limit)
    # Both source expressions are undefined for k=-1. A generic finite
    # answer1 after series cancellation would conceal this specialization.
    assert (k+1).subs(k,-1) == 0


@pytest.mark.parametrize('value', ('-3','-2','-1/2','0','1/3','1','2'))
def test_concrete_exponents_preserve_exact_reference_results(value):
    k = sx.sympify(value); k_ref = sp.sympify(value)
    f,n = issue(sx,k)
    f_ref,n_ref = issue(sp,k_ref)
    result = sx.limit(f,n,sx.oo)
    expected = sp.limit(f_ref,n_ref,sp.oo)
    assert not result.has(sx.Limit)
    assert sx.simplify(result-sx.sympify(str(expected))) == 0


@pytest.mark.parametrize('scale', (False,True))
def test_degenerate_exponent_does_not_become_a_finite_continuation(scale):
    f,n = issue(sx,sx.Integer(-1),numerator=sx.Integer(0) if scale else None)
    result = sx.limit(f,n,sx.oo)
    assert result.is_finite is not True
    assert result != sx.Integer(1)


def test_unknown_shift_difference_retains_undecidable_domain():
    k = sx.Symbol('k',positive=True)
    a,b = sx.symbols('a b',positive=True)
    f,n = issue(sx,k,a,b,numerator=(a-b)*(k+1))
    result = sx.limit(f,n,sx.oo)
    assert isinstance(result,sx.Limit)


def test_certified_positive_shift_difference():
    k,a = sx.symbols('k a',positive=True)
    f,n = issue(sx,k,2*a,a)
    result = sx.limit(f,n,sx.oo)
    assert not result.has(sx.Limit) and sx.cancel(result-1/(a*(k+1))) == 0


@pytest.mark.parametrize('value', ('I','1+I'))
def test_complex_exponents_have_an_analytic_certificate(value):
    k = sx.sympify(value)
    f,n = issue(sx,k)
    result = sx.limit(f,n,sx.oo)
    assert not result.has(sx.Limit) and sx.simplify(result-1/(k+1)) == 0
    # SymPy1.14's direct limit raises TypeError (non-real series order), so
    # prove the analytic first-order denominator coefficient separately.
    w = sp.Symbol('w',positive=True)
    kr = sp.sympify(value)
    assert sp.diff((1+w)**(kr+1),w).subs(w,0) == kr+1
    assert ((1+w)**kr).subs(w,0) == 1


@pytest.mark.parametrize('base,coefficient', (('-1-I/n','-I'), ('-1+I/n','I'),
                                            ('-1-1/n','I'), ('1+I+1/n','sqrt(1+I)')))
def test_principal_branch_sides_and_off_cut_coefficients(base,coefficient):
    k = sx.Symbol('k',positive=True)
    f,n = issue(sx,k)
    factor = sx.sqrt(sx.sympify(base,locals={'n':n}))
    expected = sx.sympify(coefficient)/(k+1)
    result = sx.limit(factor*f,n,sx.oo)
    assert not result.has(sx.Limit) and sx.simplify(result-expected) == 0


def test_certified_negative_parameter_shift():
    p = sx.Symbol('p',positive=True)
    f,n = issue(sx,-p-1)
    result = sx.limit(f,n,sx.oo)
    assert not result.has(sx.Limit) and sx.cancel(result+1/p) == 0


@pytest.mark.parametrize('precision', (1,53,139,140,211))
def test_inexact_scalar_coefficients_keep_their_stored_binary_payload(precision):
    k = sx.Symbol('k',positive=True)
    coefficient = sx.Float('0.123456789012345678901234567890123456789',precision=precision)
    f,n = issue(sx,k,numerator=coefficient)
    result = sx.limit(f,n,sx.oo)
    assert not result.has(sx.Limit)
    assert sx.simplify(result-coefficient/(k+1)) == 0
    atoms = list(result.atoms(sx.Float))
    assert len(atoms) == 1
    assert atoms[0]._mpf_ == coefficient._mpf_ and atoms[0]._prec == coefficient._prec
