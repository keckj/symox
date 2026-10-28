"""Native scalar RR/CC Gauss-Jordan: precision, pivots and API parity."""
import signal

import pytest
import sympy as sp
import symox as sx
from mpmath.libmp import dps_to_prec, prec_to_dps


def exact(value):
    if isinstance(value, sx.Basic):
        value = sx.to_sympy(value)
    value = sp.sympify(value)
    return value.xreplace({f: sp.Rational(f) for f in value.atoms(sp.Float)})


def exact_matrix(matrix):
    return sp.Matrix(matrix.rows, matrix.cols, lambda i, j: exact(matrix[i, j]))


def matrix(module, n, precision, complex_, kind='Matrix'):
    def entry(i, j):
        real = module.Float(n + 2, precision=precision) if i == j else module.Rational((i + 2)*(j + 1)%7 - 3, 16)
        return real + (module.I*module.Rational((i + j)%5 - 2, 32) if complex_ else 0)
    return getattr(module, kind)(n, n, entry)


def assert_small_exact(value, bits):
    value = sp.expand_complex(value)
    squared = sp.re(value)**2 + sp.im(value)**2
    assert squared < sp.Rational(1, 2)**(2*bits)


@pytest.mark.parametrize('kind', ['Matrix', 'ImmutableMatrix', 'SparseMatrix', 'ImmutableSparseMatrix'])
@pytest.mark.parametrize('precision', [53, 80, 139, 140, 211])
@pytest.mark.parametrize('complex_', [False, True])
def test_inverse_preserves_matrix_class_selected_precision_and_stored_binary_input(kind, precision, complex_):
    original = matrix(sx, 3, precision, complex_, kind)
    snapshot = original.copy()
    actual = original.inv()
    reference = matrix(sp, 3, precision, complex_, kind).inv()
    assert type(actual) is type(original)
    assert original == snapshot
    assert actual.shape == reference.shape
    # This bounds every identity entry over the exact stored input values,
    # independently of any reference evaluator or Float comparison.
    residual = exact_matrix(original)*exact_matrix(actual) - sp.eye(3)
    for value in residual:
        assert_small_exact(value, precision - 10)
    for value, expected in zip(actual, reference):
        assert value == 0 or value.has(sx.Float)
        presentation = dps_to_prec(prec_to_dps(precision))
        assert {f._prec for f in value.atoms(sx.Float)} <= {presentation}
        assert all(f._mpf_[3] <= min(precision, presentation) for f in value.atoms(sx.Float))
        # These well-conditioned cases take the same floating operations
        # and pivot choices: compare their actual binary coefficients too.
        assert exact(value) == exact(expected)


@pytest.mark.parametrize('complex_', [False, True])
@pytest.mark.parametrize('precision', [53, 140])
def test_partial_pivoting_and_exact_zero_entries(complex_, precision):
    f = sx.Float('2.0', precision=precision)
    original = sx.Matrix([[0, f*sx.I if complex_ else f], [2*f, 0]])
    result = original.inv()
    assert result[0, 0] == result[1, 1] == sx.Integer(0)
    assert exact_matrix(original)*exact_matrix(result) == sp.eye(2)
    assert result[0, 1] == sx.Float('0.25', precision=precision)
    assert result[1, 0] == sx.Float('0.5', precision=precision)*(-sx.I if complex_ else 1)


@pytest.mark.parametrize('n', [10, 20])
@pytest.mark.parametrize('precision', [53, 140])
@pytest.mark.parametrize('complex_', [False, True])
def test_benchmark_size_inverses_reconstruct_over_exact_binary_rationals(n, precision, complex_):
    original = matrix(sx, n, precision, complex_)
    result = original.inv()
    # Gaussian rational DomainMatrix multiplication keeps this independent
    # certificate efficient without sampling or performing Float products.
    product = (exact_matrix(original).to_DM()*exact_matrix(result).to_DM()).to_Matrix()
    for value in product - sp.eye(n):
        assert_small_exact(value, precision - 10)


@pytest.mark.parametrize('precision', [53, 140])
def test_mixed_float_precisions_do_not_demote_to_double(precision):
    high = sx.Float('1.2345678901234567890123456789', precision=precision)
    low = sx.Float('0.25')
    original = sx.Matrix([[high, low], [sx.Rational(1, 7), sx.Rational(3, 2)]])
    assert original._inexact_domain() == (False, precision)
    result = original.inv()
    for value in exact_matrix(original)*exact_matrix(result) - sp.eye(2):
        assert_small_exact(value, precision - 8)
    assert {f._prec for value in result for f in value.atoms(sx.Float)} == {precision}


@pytest.mark.parametrize('precision', [139, 211])
def test_complex_real_and_imaginary_outputs_keep_reference_display_precision(precision):
    original = sx.Matrix([[2*sx.I*sx.Float(1, precision=precision), 0],
                          [0, sx.Float(4, precision=precision)]])
    assert original._inexact_domain() == (True, precision)
    result = original.inv()
    presentation = dps_to_prec(prec_to_dps(precision))
    assert result[0, 1] == result[1, 0] == sx.Integer(0)
    assert result[0, 0] == -sx.I*sx.Float('0.5', precision=presentation)
    assert result[1, 1] == sx.Float('0.25', precision=presentation)
    assert exact_matrix(original)*exact_matrix(result) == sp.eye(2)


@pytest.mark.parametrize('precision', [1, 2, 3])
def test_lowest_domain_precisions_remain_supported(precision):
    actual = sx.Matrix([[sx.Float(2, precision=precision)]]).inv()[0]
    expected = sp.Matrix([[sp.Float(2, precision=precision)]]).inv()[0]
    assert actual._prec == expected._prec == 7
    assert exact(actual) == exact(expected) == sp.Rational(1, 2)


@pytest.mark.parametrize('complex_', [False, True])
def test_nearly_singular_high_precision_input_is_not_mistaken_for_singular(complex_):
    precision = 140
    one = sx.Float(1, precision=precision)
    delta = sx.Float(sx.Rational(1, 2)**100, precision=precision)
    scale = 1 + sx.I if complex_ else 1
    original = sx.Matrix([[one*scale, one*scale], [one*scale, (one + delta)*scale]])
    inverse = original.inv()
    # All values here are dyadic: even the ill-conditioned reconstruction is exact.
    assert (exact_matrix(original)*exact_matrix(inverse)).applyfunc(sp.expand) == sp.eye(2)


@pytest.mark.parametrize('precision', [53, 140])
@pytest.mark.parametrize('complex_', [False, True])
def test_singular_numeric_matrices_raise_the_reference_exception(precision, complex_):
    for module in (sx, sp):
        a = module.Float('1.25', precision=precision)*(1 + module.I if complex_ else 1)
        original = module.Matrix([[a, 2*a], [2*a, 4*a]])
        exception = sx.NonInvertibleMatrixError if module is sx else sp.matrices.exceptions.NonInvertibleMatrixError
        with pytest.raises(exception):
            original.inv()


@pytest.mark.parametrize('source', ['oo', '-oo', 'zoo', 'nan', 'sqrt(2)', 'Symbol("x")'])
def test_nonfinite_and_nonscalar_domains_retain_the_existing_route(source):
    value = eval(source, vars(sx))
    original = sx.Matrix([[sx.Float('1.25'), 0], [0, value]])
    assert original._inexact_domain() is None
    with pytest.raises(NotImplementedError, match='finite rational or Float'):
        original._inverse_inexact(True, 140)
    if source in ('oo', '-oo'):
        # Existing composite-domain limitation, outside the finite kernel:
        # SymPy returns a zero reciprocal of the infinite diagonal entry.
        with pytest.raises(ValueError, match='cannot build a polynomial'):
            original.inv()
    elif source in ('zoo', 'nan'):
        with pytest.raises(sx.polys.polyerrors.CoercionFailed):
            original.inv()
    else:
        for value in exact_matrix(original)*exact_matrix(original.inv()) - sp.eye(2):
            assert_small_exact(sp.simplify(value), 45)


def test_native_inverse_does_not_enter_mpmath_or_python_gauss_jordan(monkeypatch):
    import mpmath
    from symox import _matrix_inverse
    original = matrix(sx, 4, 140, True)

    def forbidden(*args, **kwargs):
        raise AssertionError('native inversion entered the former Python numeric path')

    monkeypatch.setattr(mpmath, 'workprec', forbidden)
    monkeypatch.setattr(mpmath, 'mpf', forbidden)
    monkeypatch.setattr(mpmath, 'mpc', forbidden)
    assert _matrix_inverse._inexact_domain(original) == (True, 140)
    result = original.inv()
    assert result.shape == (4, 4)


@pytest.mark.skipif(not hasattr(signal, 'setitimer'), reason='requires POSIX interval timer')
def test_pending_interruption_propagates_from_native_elimination():
    original = matrix(sx, 60, 140, True)
    sentinel = KeyboardInterrupt('matrix inversion interrupted')

    def interrupt(signum, frame):
        raise sentinel

    previous = signal.signal(signal.SIGALRM, interrupt)
    try:
        signal.setitimer(signal.ITIMER_REAL, 0.001)
        with pytest.raises(KeyboardInterrupt) as caught:
            original._inverse_inexact(True, 140)
        assert caught.value is sentinel
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)
