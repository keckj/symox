"""Output forms of sympy.stats documentation examples, hard-coded from SymPy 1.14.

These callables used to give an equal result in another form (support
Piecewise layout, merged exponentials, unsimplified binomials). The
distributions now build the same expressions as the reference.
"""
from symox import (Matrix, MatrixSymbol, S, Symbol, Trace, Determinant, exp, gamma,
                   pi, pretty, sqrt, symbols)
from symox.stats import (Arcsin, Benini, Dagum, Erlang, ExGaussian, Exponential, FlorySchulz,
                         Frechet, Gamma, GammaInverse, Geometric, Gompertz, Kumaraswamy, Lomax,
                         MatrixGamma, Nakagami, NormalGamma, ShiftedGompertz, Trapezoidal,
                         Uniform, UniformSum, Wishart, cdf, density, P, E)
from symox.stats.joint_rv_types import GeneralizedMultivariateLogGammaOmega

z = Symbol('z')


def _pretty(expr):
    return '\n'.join(line.rstrip() for line in pretty(expr, use_unicode=False).splitlines())


def test_cdf_piecewise_layouts():
    a, b = Symbol('a', real=True), Symbol('b', real=True)
    assert str(cdf(Arcsin('x', a, b))(z)) == \
        'Piecewise((0, a > z), (2*asin(sqrt((-a + z)/(-a + b)))/pi, b >= z), (1, True))'
    a, b = Symbol('a', negative=True), Symbol('b', positive=True)
    assert str(cdf(Uniform('x', a, b))(z)) == \
        'Piecewise((0, a > z), ((-a + z)/(-a + b), b >= z), (1, True))'
    p, a, b = symbols('p a b', positive=True)
    assert str(cdf(Dagum('x', p, a, b))(z)) == 'Piecewise(((1 + (z/b)**(-a))**(-p), z >= 0), (0, True))'
    assert str(cdf(Kumaraswamy('x', a, b))(z)) == \
        'Piecewise((0, z < 0), (1 - (1 - z**a)**b, z <= 1), (1, True))'
    assert str(cdf(GammaInverse('x', a, b))(z)) == \
        'Piecewise((uppergamma(a, b/z)/gamma(a), z > 0), (0, True))'
    l = Symbol('lambda', positive=True)
    assert str(cdf(Exponential('x', l))(z)) == 'Piecewise((1 - exp(-lambda*z), z >= 0), (0, True))'
    s, m = Symbol('s', positive=True), Symbol('m', real=True)
    assert str(cdf(Frechet('x', a, s, m))(z)) == \
        'Piecewise((exp(-1/((-m + z)/s)**a), m <= z), (0, True))'
    al, x = symbols('a, l', positive=True), Symbol('x')
    assert str(cdf(Lomax('X', *al))(x)) == 'Piecewise((1 - 1/(1 + x/l)**a, x >= 0), (0, True))'
    mu, omega = symbols('mu omega', positive=True)
    assert str(cdf(Nakagami('x', mu, omega))(z)) == \
        'Piecewise((lowergamma(mu, mu*z**2/omega)/gamma(mu), z > 0), (0, True))'
    alpha, beta, sigma = symbols('alpha beta sigma', positive=True)
    assert str(cdf(Benini('x', alpha, beta, sigma))(z)) == \
        'Piecewise((1 - exp(-alpha*log(z/sigma) - beta*log(z/sigma)**2), sigma <= z), (0, True))'
    n = Symbol('n', integer=True)
    assert str(cdf(UniformSum('x', n))(z)) == ('Piecewise((0, z < 0), (Sum((-1)**_k*(-_k + z)**n*'
        'binomial(n, _k), (_k, 0, floor(z)))/factorial(n), n >= z), (1, True))')


def test_gamma_family_cdfs():
    k, l = Symbol('k', integer=True, positive=True), Symbol('l', positive=True)
    assert _pretty(cdf(Erlang('x', k, l))(z)) == '\n'.join([
        '/lowergamma(k, l*z)', '|------------------  for z > 0', '<     Gamma(k)', '|',
        '\\        0           otherwise'])
    k, theta = symbols('k theta', positive=True)
    X = Gamma('x', k, theta)
    assert str(cdf(X)(z)) == 'Piecewise((lowergamma(k, z/theta)/gamma(k), z > 0), (0, True))'
    # With integration options the reference integrates the density.
    assert _pretty(cdf(X, meijerg=True)(z)) == '\n'.join([
        '/            /     z  \\', '|k*lowergamma|k, -----|', '|            \\   theta/',
        '<----------------------  for z >= 0', '|     Gamma(k + 1)', '|',
        '\\          0             otherwise'])


def test_exgaussian_cdf():
    mean, std, rate = Symbol('mu'), Symbol('sigma', positive=True), Symbol('lamda', positive=True)
    assert str(cdf(ExGaussian('x', mean, std, rate))(z)) == (
        '-(erf(sqrt(2)*(-lamda**2*sigma**2 + lamda*(-mu + z))/(2*lamda*sigma))/2 + 1/2)*'
        'exp(lamda**2*sigma**2/2 - lamda*(-mu + z)) + erf(sqrt(2)*(-mu + z)/(2*sigma))/2 + 1/2')


def test_density_forms():
    b, eta = symbols('b eta', positive=True)
    assert str(density(Gompertz('x', b, eta))(z)) == 'b*eta*exp(eta)*exp(b*z)*exp(-eta*exp(b*z))'
    x = Symbol('x')
    assert str(density(ShiftedGompertz('x', b, eta))(x)) == \
        'b*(eta*(1 - exp(-b*x)) + 1)*exp(-b*x)*exp(-eta*exp(-b*x))'
    assert str(density(Geometric('x', S.One/5))(z)) == '(4/5)**(z - 1)/5'
    assert str(density(FlorySchulz('x', S.One/5))(z)) == '(4/5)**(z - 1)*z/25'
    y = Symbol('y')
    assert str(density(NormalGamma('x', 0, 1, 2, 3))(y, z)) == \
        '9*sqrt(2)*z**(3/2)*exp(-3*z)*exp(-y**2*z/2)/(2*sqrt(pi))'
    alpha, beta, sigma = symbols('alpha beta sigma', positive=True)
    assert _pretty(density(Benini('x', alpha, beta, sigma))(z)) == '\n'.join([
        '/                  /  z  \\\\             /  z  \\           2/  z  \\',
        '|        2*beta*log|-----||  - alpha*log|-----| - beta*log |-----|',
        '|alpha             \\sigma/|             \\sigma/            \\sigma/',
        '|----- + -----------------|*e',
        '\\  z             z        /'])
    mu, omega = symbols('mu omega', positive=True)
    assert _pretty(density(Nakagami('x', mu, omega))(z)) == '\n'.join([
        '                                2', '                           -mu*z',
        '                           -------', '    mu      -mu  2*mu - 1   omega',
        '2*mu  *omega   *z        *e', '----------------------------------',
        '            Gamma(mu)'])


def test_trapezoidal_density():
    a, b, c, d = symbols('a b c d')
    assert _pretty(density(Trapezoidal('x', a, b, c, d))(z)) == '\n'.join([
        '/       -2*a + 2*z', '|-------------------------  for And(a <= z, b > z)',
        '|(-a + b)*(-a - b + c + d)', '|', '|           2',
        '|     --------------        for And(b <= z, c > z)', '<     -a - b + c + d', '|',
        '|        2*d - 2*z', '|-------------------------  for And(d >= z, c <= z)',
        '|(-c + d)*(-a - b + c + d)', '|', '\\            0                     otherwise'])
    # Degenerate shapes (triangles) keep a valid density.
    T = Trapezoidal('t', 0, 0, 1, 2)
    assert density(T)(S.Half) == S(2)/3
    assert P(T < 1) == S(2)/3
    assert E(Trapezoidal('u', 0, 1, 2, 3)) == S(3)/2


def test_generalized_multivariate_log_gamma_omega():
    omega = Matrix([[1, S.Half, S.Half], [S.Half, 1, S.Half], [S.Half, S.Half, 1]])
    G = GeneralizedMultivariateLogGammaOmega('G', omega, 1, [1, 1, 1], [1, 1, 1])
    y = symbols('y_1:4', positive=True)
    assert str(density(G)(y[0], y[1], y[2])) == (
        'sqrt(2)*Sum((1 - sqrt(2)/2)**n*exp((n + 1)*(y_1 + y_2 + y_3) - exp(y_1) - exp(y_2) - '
        'exp(y_3))/gamma(n + 1)**3, (n, 0, oo))/2')


def test_matrix_gamma_and_wishart_densities():
    a, b = symbols('a b', positive=True)
    M = MatrixGamma('M', a, b, [[2, 1], [1, 2]])
    assert str(density(M)([[1, 0], [0, 1]]).doit()) == \
        'exp(-4/(3*b))/(3**a*sqrt(pi)*b**(2*a)*gamma(a)*gamma(a - 1/2))'
    X = MatrixSymbol('X', 2, 2)
    expected = (exp(Trace(Matrix([[-S(2)/3, S(1)/3], [S(1)/3, -S(2)/3]])*X)/b)*Determinant(X)**(a - S(3)/2)
                /(3**a*sqrt(pi)*b**(2*a)*gamma(a)*gamma(a - S.Half)))
    assert density(M)(X).doit() == expected
    n = Symbol('n', positive=True)
    W = Wishart('W', n, [[2, 1], [1, 2]])
    expected = (exp(Trace(Matrix([[-S(1)/3, S(1)/6], [S(1)/6, -S(1)/3]])*X))*Determinant(X)**(n/2 - S(3)/2)
                /(2**n*3**(n/2)*sqrt(pi)*gamma(n/2)*gamma(n/2 - S.Half)))
    assert density(W)(X).doit() == expected
