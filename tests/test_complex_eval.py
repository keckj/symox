"""Native re/im/arg/Abs/sign/conjugate evaluators against SymPy."""
import pytest
import sympy as sp
import symox as sx


def environment(module):
    names = ['I', 'pi', 'E', 'oo', 'zoo', 'nan', 'Rational', 'Float', 'sqrt',
             'sin', 'cos', 'tan', 'cot', 'sec', 'csc', 'sinh', 'cosh', 'tanh', 'coth',
             'exp', 'log', 'Abs', 'arg', 'sign', 're', 'im', 'conjugate',
             'asin', 'acos', 'atan', 'atan2', 'asinh', 'acosh', 'atanh',
             'floor', 'ceiling', 'gamma', 'erf', 'Max', 'Min', 'Heaviside',
             'transpose', 'adjoint', 'Matrix']
    result = {name: getattr(module, name) for name in names}
    for name, flags in [('x', {}), ('y', {}), ('r', {'real': True}), ('t', {'real': True}),
                        ('p', {'positive': True}), ('n', {'negative': True}),
                        ('j', {'imaginary': True}), ('k', {'integer': True}),
                        ('z', {'zero': True}), ('c', {'complex': True}),
                        ('q', {'nonnegative': True})]:
        result[name] = module.Symbol(name, **flags)
    result['f'] = module.Function('f')
    return result


def evaluate(module, source):
    return eval(source, {'__builtins__': {}}, environment(module))


CASES = [
    # re
    're(x)', 're(r)', 're(j)', 're(I*x)', 're(2+3*I)', 're(x+I*y)', 're(r+I*t)', 're(x*y)',
    're(I*r)', 're(re(x))', 're(im(x))', 're(conjugate(x))', 're(exp(I*r))', 're(oo)',
    're(zoo)', 're(nan)', 're(f(x))', 're(sqrt(-2))', 're(sqrt(2)*I)', 're(log(2))',
    're(Abs(x))', 're(arg(x))', 're(j*r)', 're(c)', 're(2*x)', 're(x/2)', 're(I*j)', 're(k)',
    're(sqrt(-2+I))', 're(exp(2+I))', 're(sin(1+I))', 're(log(1+I))', 're((1+I)**3)',
    're(Float(2.5)+3*I)', 're(2*I*x)', 're(x+2*I*y)', 're(x*I*r)', 're(I*x*y)', 're(r+I*x)',
    're(x**2)', 're(r**2)', 're(sqrt(x))', 're(sqrt(r))', 're(sqrt(p))', 're(sqrt(n))',
    're(exp(x))', 're(exp(r))', 're(exp(I*x))', 're(exp(x+I*y))', 're(log(x))', 're(log(p))',
    're(log(n))', 're(log(x+I))', 're(sin(x))', 're(sin(I*r))', 're(cos(r+I*t))',
    're(1/(r+I*t))', 're((r+I*t)**2)', 're(x*conjugate(x))', 're(Abs(x)+I*Abs(y))',
    're(re(x)+I*im(x))', 're(f(x)+I)', 're(r*f(x))', 're(2*x+3*y+I*z)', 're(sign(x))',
    're(x/y)', 're(1/x)', 're(I**x)', 're(x**I)', 're(atan(x))', 're(Rational(1,2)+I/3)',
    're(sqrt(2)*I+sqrt(3))', 're(pi*I)', 're(E*I + 1)', 're(j + r)', 're(j*j)', 're(j**2)',
    're(j**3)', 're(k*I)', 're(x*y*I)', 're(re(x)*I)', 're(conjugate(x)+I)',
    # im
    'im(x)', 'im(r)', 'im(j)', 'im(I*x)', 'im(2+3*I)', 'im(x+I*y)', 'im(r+I*t)', 'im(x*y)',
    'im(I*r)', 'im(re(x))', 'im(im(x))', 'im(conjugate(x))', 'im(oo)', 'im(zoo)', 'im(nan)',
    'im(I*j)', 'im(2*j)', 'im(2*x)', 'im(sqrt(-2+I))', 'im(exp(2+I))', 'im(log(1+I))',
    'im((1+I)**Rational(1,3))', 'im(2*I*x)', 'im(x+2*I*y)', 'im(x*I*r)', 'im(I*x*y)',
    'im(r+I*x)', 'im(x**2)', 'im(r**2)', 'im(sqrt(x))', 'im(exp(x))', 'im(exp(I*x))',
    'im(exp(I*r))', 'im(exp(x+I*y))', 'im(log(x))', 'im(log(n))', 'im(log(x+I))', 'im(sin(x))',
    'im(sin(I*r))', 'im(cos(r+I*t))', 'im(1/(r+I*t))', 'im((r+I*t)**2)', 'im(x*conjugate(x))',
    'im(re(x)+I*im(x))', 'im(f(x))', 'im(f(x)+I)', 'im(r*f(x))', 'im(2*x+3*y+I*z)', 'im(x/y)',
    'im(1/x)', 'im(Rational(1,2)+I/3)', 'im(sqrt(2)*I+sqrt(3))', 'im(pi*I)', 'im(j + r)',
    'im(j*r)', 'im(j**2)', 'im(j**3)', 'im(k*I)', 'im(x*y*I)', 'im(Abs(x)*I)', 'im(re(x)*I)',
    'im(conjugate(x)*I)',
    # conjugate
    'conjugate(x)', 'conjugate(r)', 'conjugate(j)', 'conjugate(I)', 'conjugate(2+3*I)',
    'conjugate(x+I*y)', 'conjugate(x*y)', 'conjugate(exp(x))', 'conjugate(sin(x))',
    'conjugate(log(x))', 'conjugate(log(p))', 'conjugate(sqrt(x))', 'conjugate(x**2)',
    'conjugate(x**y)', 'conjugate(conjugate(x))', 'conjugate(Abs(x))', 'conjugate(re(x))',
    'conjugate(f(x))', 'conjugate(oo)', 'conjugate(zoo)', 'conjugate(nan)', 'conjugate(pi)',
    'conjugate(asin(x))', 'conjugate(atan(x))', 'conjugate(sinh(x))', 'conjugate(exp(I*x))',
    'conjugate(I*x)', 'conjugate(1/x)', 'conjugate(x**Rational(1,2))', 'conjugate(x**k)',
    'conjugate(x**r)', 'conjugate(p**x)', 'conjugate(r**x)', 'conjugate(sin(r+I*t))',
    'conjugate(im(x))', 'conjugate(arg(x))', 'conjugate(sign(x))', 'conjugate(atan2(x,y))',
    'conjugate(Abs(x)+I)', 'conjugate(-oo)', 'conjugate(j*r)', 'conjugate(x*j)',
    'conjugate(gamma(x))', 'conjugate(log(x+I))',
    # arg
    'arg(x)', 'arg(p)', 'arg(n)', 'arg(r)', 'arg(j)', 'arg(I)', 'arg(-I)', 'arg(1+I)',
    'arg(-1-I)', 'arg(sqrt(2)*(1+I))', 'arg(exp(I*r))', 'arg(2*x)', 'arg(p*x)', 'arg(n*x)',
    'arg(I*p)', 'arg(0)', 'arg(z)', 'arg(oo)', 'arg(-oo)', 'arg(zoo)', 'arg(nan)',
    'arg(arg(x))', 'arg(conjugate(x))', 'arg(Abs(x))', 'arg(sqrt(x))', 'arg(exp(x))',
    'arg(exp(I*pi/3))', 'arg(x+I*y)', 'arg(Rational(-1,2))', 'arg(Float(-2))', 'arg(1/x)',
    'arg(sqrt(-2+I))', 'arg(1.5+2.0*I)', 'arg(-2*x-2*y)', 'arg(Float(0))', 'arg(Float(2))',
    'arg(exp(5*I))', 'arg(exp(2+3*I))', 'arg(1+sqrt(2)*I)', 'arg(f(x))', 'arg(-x)', 'arg(I*x)',
    'arg(-p*x)', 'arg(sqrt(2)*x)', 'arg(x**2)', 'arg(pi*I)', 'arg(2-3*I)', 'arg(-1+I)',
    'arg(-j)', 'arg(I*r)', 'arg(sqrt(-2))', 'arg(-sqrt(2))', 'arg(E)', 'arg(-pi)',
    'arg(arg(arg(x)))', 'arg(x*y)', 'arg(2*p*x*y)', 'arg(n*I)', 'arg(3*I)', 'arg(-3*I)',
    'arg((1+I)**2)', 'arg(Rational(1,2)+I*sqrt(3)/2)',
    # Abs
    'Abs(x)', 'Abs(r)', 'Abs(p)', 'Abs(n)', 'Abs(j)', 'Abs(I)', 'Abs(3+4*I)', 'Abs(x+I*y)',
    'Abs(r+I*t)', 'Abs(x*y)', 'Abs(exp(x))', 'Abs(exp(I*r))', 'Abs(exp(r))',
    'Abs(conjugate(x))', 'Abs(re(x))', 'Abs(Abs(x))', 'Abs(arg(x))', 'Abs(x**2)', 'Abs(r**2)',
    'Abs(sqrt(x))', 'Abs(x**y)', 'Abs(1/x)', 'Abs(oo)', 'Abs(-oo)', 'Abs(zoo)', 'Abs(nan)',
    'Abs(I*oo)', 'Abs(sign(x))', 'Abs(f(x))', 'Abs(k)', 'Abs(-k)', 'Abs(sqrt(-2+I))',
    'Abs(Float(-2.5))', 'Abs(1.0+2.0*I)', 'Abs(sqrt(2)+I)', 'Abs(2**x)', 'Abs(x**-1)',
    'Abs(I**x)', 'Abs((-2)**x)', 'Abs(2*x)', 'Abs(-x)', 'Abs(p*x)', 'Abs(x/y)', 'Abs(1/(1+I))',
    'Abs(1/(x+I))', 'Abs(r**3)', 'Abs(r**-3)', 'Abs(n**3)', 'Abs(p**x)', 'Abs(n**x)',
    'Abs(p**(x+I))', 'Abs(r+I*oo)', 'Abs(I*r)', 'Abs(I*p)', 'Abs(-j)', 'Abs(I*j)',
    'Abs(2*I*x)', 'Abs(exp(I*x))', 'Abs(x - I)', 'Abs(conjugate(x)*x)', 'Abs(1/sqrt(x))',
    'Abs((x+I*y)**2)', 'Abs(re(x)+I*im(x))', 'Abs(sin(x))', 'Abs(sin(r))', 'Abs(log(p))',
    'Abs(r*t)', 'Abs(k*I)', 'Abs(-1-I)', 'Abs(sqrt(3)+sqrt(2)*I)', 'Abs(2+I*sqrt(5))',
    'Abs(x)**2', 'Abs(r)**2', 'Abs(x)*Abs(y)',
    # sign
    'sign(x)', 'sign(1+I)', 'sign(I*r)', 'sign(j)', 'sign(conjugate(x))', 'sign(Abs(x))',
    'sign(I*p)', 'sign(-I*p)', 'sign(2*I*x)', 'sign(-3*x)', 'sign(p*x)', 'sign(n*x)',
    'sign(j*p)', 'sign(sqrt(-2))', 'sign(1-sqrt(2))', 'sign(I*(1-sqrt(2)))', 'sign(sign(x))',
    'sign(nan)', 'sign(z)', 'sign(r*t)', 'sign(x*y)', 'sign(sqrt(1-sqrt(2)))', 'sign(oo)',
    'sign(-oo)', 'sign(zoo)', 'sign(I*oo)', 'sign(k)', 'sign(exp(x))', 'sign(exp(r))',
]

EXTRA = [
    'transpose(x)', 're(Matrix([[x, I*y]]))', 'im(Matrix([[x, I*y]]))',
    'conjugate(Matrix([[x, I*y]]))', 'Abs(Matrix([[x]]))',
    'Matrix([[x, 1+I],[exp(I*r), y]]).as_real_imag()',
    'arg(x).is_real', 'arg(x).is_finite', 'conjugate(x).is_real', 'conjugate(r).is_real',
    're(c).is_finite', 'arg(p).is_zero', 'Abs(x).is_zero', 'Abs(p).is_positive',
    're(x).diff(x)', 'im(x).diff(x)', 'Abs(x).diff(x)', 'Abs(r).diff(r)',
    'arg(x).diff(x)', 'conjugate(x).diff(x)',
]

# Known differences, kept strict so a fix is noticed:
# - symox treats re(x) as finite, so Abs(x+oo) evaluates to oo;
# - im(sqrt(n)) is the equal value sqrt(-n) because of Mul canonicalization;
# - -I*oo prints as (-1)*oo*I.
KNOWN = ['Abs(x+oo)', 'im(sqrt(n))', 'conjugate(oo*I)']


@pytest.mark.parametrize('source', CASES + EXTRA)
def test_evaluator_matches_sympy(source):
    actual = evaluate(sx, source)
    expected = evaluate(sp, source)
    assert str(actual) == str(expected)
    if isinstance(expected, sp.Expr) and not isinstance(expected, sp.MatrixBase):
        assert sx.to_sympy(actual) == expected


@pytest.mark.parametrize('source', KNOWN)
@pytest.mark.xfail(strict=True, reason='documented difference')
def test_known_differences(source):
    assert str(evaluate(sx, source)) == str(evaluate(sp, source))


@pytest.mark.parametrize('k', [-5, -3, -2, 2, 3, 4, 5])
def test_abs_power_of_real_matches_sympy(k):
    r, R = sx.Symbol('r', real=True), sp.Symbol('r', real=True)
    assert sx.to_sympy(sx.Abs(r)**k) == sp.Abs(R)**k


@pytest.mark.parametrize('source', ['re(x)', 'im(x)', 'arg(x)', 'Abs(x)', 'sign(x)',
                                    'conjugate(x)', 'arg(1+I)', 'Abs(3+4*I)', 'sign(2*I*x)',
                                    're(x*I*r)', 'conjugate(exp(I*x))'])
def test_evaluated_values_round_trip(source):
    value = evaluate(sx, source)
    assert sx.sympify(value.srepr()) == value
    assert sx.sympify(sx.to_sympy(value)) == value


@pytest.mark.parametrize('point', [1+2j, -3+0.5j, -2-1j, 0.25-4j, -1.5+0j])
def test_arg_and_abs_numeric_values(point):
    import cmath
    value = sx.Float(point.real) + sx.I*sx.Float(point.imag)
    assert abs(complex(sx.arg(value).evalf()) - cmath.phase(point)) < 1e-12
    assert abs(complex(sx.Abs(value).evalf()) - abs(point)) < 1e-12
    assert abs(complex(sx.re(value).evalf()) - point.real) < 1e-12
    assert abs(complex(sx.im(value).evalf()) - point.imag) < 1e-12
