"""Native conjugate-coefficient preparation preserves exact principal roots."""
import pytest
import sympy as sp
import symox as sx

from symox.integrals.heurisch import _cartesian_numeric_coefficients


z, w, a = sx.symbols('z w a')


def native(polynomial, generators=(z, w)):
    return sx._core._heurisch_cartesian_numeric_coefficients(polynomial, generators)


@pytest.mark.parametrize('source,expected', [
    ('z+sqrt(Rational(1,2)+sqrt(3)*I/2)', 'z+sqrt(3)/2+I/2'),
    ('z+sqrt(Rational(1,2)-sqrt(3)*I/2)', 'z+sqrt(3)/2-I/2'),
    ('z+sqrt(-3+4*I)', 'z+1+2*I'),
    ('z+sqrt(-3-4*I)', 'z+1-2*I'),
    ('z**2+(2+I)*sqrt(Rational(1,2)+sqrt(3)*I/2)*z+3*I',
     'z**2+z*(-Rational(1,2)+sqrt(3))+I*z*(sqrt(3)/2+1)+3*I'),
    ('z*w*sqrt(Rational(1,2)+sqrt(3)*I/2)+(I+1)*w+4',
     'sqrt(3)*w*z/2+I*w*z/2+w+I*w+4'),
])
def test_native_matches_preport_exact_coefficients_and_principal_branch(source, expected):
    # These exact forms also match the Python helper in the pre-port native
    # validation snapshot 345253. The independent principal value is checked
    # with SymPy, including conjugates on both sides of the negative real axis.
    namespace = dict(vars(sx), z=z, w=w)
    polynomial = eval(source, namespace)
    actual = native(polynomial)
    assert actual == eval(expected, namespace)
    assert _cartesian_numeric_coefficients(polynomial, (z, w)) == actual
    reference = eval(source, dict(vars(sp), z=sp.Symbol('z'), w=sp.Symbol('w')))
    assert sp.simplify(sp.expand_complex(sx.to_sympy(actual)-reference)) == 0


def test_nested_principal_radicals_keep_an_exact_algebraic_cartesian_form():
    actual = native(z+sx.sqrt(2+sx.sqrt(3+4*sx.I)))
    real = sx.sqrt(2+sx.sqrt(17)/2)
    assert actual == z+real+sx.I/(2*real)
    root = sx.to_sympy(actual-z)
    # Positive real parts select both principal square roots uniquely.
    # This avoids relying on a polar half-angle simplifier to prove equality.
    inner = 2+sp.I
    assert sp.expand(inner**2) == 3+4*sp.I
    assert sp.re(inner).is_positive is True
    assert sp.simplify(sp.expand(root**2)-(2+inner)) == 0
    assert sx.to_sympy(real).is_positive is True


@pytest.mark.parametrize('source', [
    'z+f(1)+I', 'z+sqrt(a+I)', 'z+(a+I)*w',
    '1/z+sqrt(1+I)', 'sin(z)+I',
])
def test_opaque_parameter_and_nonpolynomial_inputs_preserve_the_original(source):
    polynomial = eval(source, dict(vars(sx), z=z, w=w, a=a, f=sx.Function('f')))
    actual = native(polynomial)
    assert actual == polynomial
    assert actual.free_symbols == polynomial.free_symbols


@pytest.mark.parametrize('digits', [40, 85])
def test_floating_coefficients_keep_their_declared_precision(digits):
    coefficient = sx.sqrt(sx.Float('0.5', digits)+sx.I*sx.Float('0.5', digits))
    actual = sx.to_sympy(native(z+coefficient)-z)
    reference = sp.sqrt(sp.Float('0.5', digits)+sp.I*sp.Float('0.5', digits)).evalf(digits+20)
    assert abs((actual-reference).evalf(digits+20)) < sp.Rational(1, 10**(digits-3))
    assert all(value._prec >= sx.Float('0.5', digits)._prec for value in actual.atoms(sp.Float))


@pytest.mark.parametrize('sign', [-1, 1])
def test_arbitrarily_small_exact_imaginary_parts_use_the_correct_principal_side(sign):
    epsilon = sx.Rational(1, 10**120)
    actual = native(z+sx.sqrt(-1+sign*sx.I*epsilon))-z
    assert not actual.has(sx.Float)
    reference = sp.sqrt(-1+sign*sp.I*sp.Rational(1, 10**120)).evalf(170)
    assert abs((sx.to_sympy(actual)-reference).evalf(170)) < sp.Rational(1, 10**150)
    assert sign*sp.im(sx.to_sympy(actual).evalf(170)) > 0


def test_noncommutative_polynomials_are_not_reordered_and_invalid_generators_keep_errors():
    from symox.polys.polyerrors import GeneratorsError

    u = sx.Symbol('u', commutative=False)
    polynomial = z+u+sx.sqrt(sx.Rational(1, 2)+sx.sqrt(3)*sx.I/2)
    assert native(polynomial, (z,)) == polynomial
    with pytest.raises(GeneratorsError, match='non-commutative generators'):
        native(polynomial, (u,))
    with pytest.raises(GeneratorsError, match='duplicated generators'):
        native(polynomial, (z, z))


@pytest.mark.parametrize('source', ['z+1', 'z+f(1)+I', 'z+sqrt(a+I)', 'z+Float("1.234",80)'])
def test_untouched_poly_inputs_keep_the_original_object_and_domain(source):
    expression = eval(source, dict(vars(sx), z=z, a=a, f=sx.Function('f')))
    polynomial = sx.Poly(expression, z)
    assert native(polynomial, (z,)) is polynomial
    assert _cartesian_numeric_coefficients(polynomial, (z,)) is polynomial


@pytest.mark.parametrize('coefficient', [sx.I, sx.sqrt(sx.Rational(1, 2)+sx.I*sx.sqrt(3)/2)])
def test_processed_poly_inputs_return_an_expression_even_when_the_value_is_unchanged(coefficient):
    polynomial = sx.Poly(z+coefficient, z)
    actual = native(polynomial, (z,))
    assert not isinstance(actual, sx.Poly)
    assert actual == _cartesian_numeric_coefficients(polynomial, (z,))


@pytest.mark.parametrize('options', [{}, {'heurisch': True}])
def test_global_floor_primitive_and_integration_options_keep_the_preport_result(options):
    x = sx.Symbol('x')
    actual = sx.integrate(1/(1+3*sx.cos(x)**2), x, **options)
    expected = (sx.atan(2*sx.tan(x/2)-sx.sqrt(3))/2
                + sx.atan(2*sx.tan(x/2)+sx.sqrt(3))/2
                + sx.pi*sx.floor((x/2-sx.pi/2)/sx.pi))
    assert actual == expected
