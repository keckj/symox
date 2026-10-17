"""Symbolic products: product(), held Product nodes and Product.is_convergent."""
import pickle

import pytest
import sympy as sp
import symox as sm


def namespace(m):
    names = "product Product Rational sqrt sin cos exp I oo pi Function symbols".split()
    d = {name: getattr(m, name) for name in names}
    k, n, m_, x, y, a, b, j = m.symbols('k n m x y a b j')
    d.update(k=k, n=n, m=m_, x=x, y=y, a=a, b=b, j=j, f=m.Function('f'),
             p=m.Symbol('p', positive=True, integer=True))
    return d


def outcome(thunk):
    try:
        return str(thunk())
    except Exception as error:
        return 'EXC ' + type(error).__name__


PRODUCTS = [
    'product(k,(k,1,n))',
    'product(k**2+1,(k,1,n))',
    'product(2**k,(k,1,n))',
    'Product(2**k,(k,1,n)).doit()',
    'product(x**k,(k,1,n))',
    'product(exp(k),(k,1,n))',
    'product(k,(k,a,b))',
    'product(1+1/k,(k,1,n))',
    'product(k,(k,3,1))',
    'product(0,(k,1,n))',
    'product(x**2,(k,1,n))',
    'product(k,(k,1,n),(n,1,m))',
    'product(k*f(k),(k,1,n))',
    'product(k*x+1,(k,1,n))',
    'product(k-1,(k,1,n))',
    'product(x**(1/k),(k,1,n))',
    'product(k**k,(k,1,n))',
    'product(x*2**k+2**k,(k,1,n))',
    'product(k,(k,1,10))',
    'product(1/k,(k,1,n))',
    'product(k/(k+1),(k,1,n))',
    'product(2*k-1,(k,1,n))',
    'product((k+1)/(k+2),(k,0,n))',
    'product(k**2-1,(k,2,n))',
    'product(x,(k,a,b))',
    'product(1,(k,1,n))',
    'product(exp(2*k+1),(k,1,n))',
    'product(k,(k,n,n))',
    'product(k,(k,1,0))',
    'product(k,(k,1,-1))',
    'product(k,(k,1,-3))',
    'product(sin(k),(k,1,3))',
    'product(f(k),(k,1,3))',
    'product(f(k),(k,1,n))',
    'product(k*(k+1),(k,1,n))',
    'product(cos(x)*k,(k,1,n))',
    'product(k**2*(k+1),(k,1,n))',
    'product(1-1/k**2,(k,2,n))',
    'product(4*k**2,(k,1,n))',
    'product(k+x,(k,1,n))',
    'product(k+y,(k,0,n-1))',
    'product(-k,(k,1,n))',
    'product(I*k,(k,1,n))',
    'product(Rational(1,2)*k,(k,1,n))',
    'product(k,(k,1,oo))',
    'product(2,(k,1,oo))',
    'product(k,(k,0,n))',
    'product(k**Rational(1,2),(k,1,n))',
    'product(sqrt(k),(k,1,n))',
    'product(1+x**k,(k,1,3))',
    'product(x**k*y,(k,1,n))',
    'product(2**k*3**k,(k,1,n))',
    'product(k**2+k,(k,1,n))',
    'product(k,(k,1,n)).subs(n,5)',
    'Product(k,(k,1,n)).subs(n,5).doit()',
    'Product(k*x,(k,1,n)).doit()',
    'Product(Product(k*j,(k,1,n)),(j,1,m)).doit()',
    'Product(k,(k,1,n),(n,1,3)).doit()',
    'product(exp(x)*k,(k,1,n))',
    'product(1/(k**2+1),(k,1,n))',
    'product((k**2+1)**2,(k,1,n))',
    'product(k+a,(k,1,n))',
    'product(k,(k,1,p))',
    'product(k,(k,1,p+2))',
    'product(k-2,(k,1,n))',
    'product(x**(k**2),(k,1,n))',
    'product(2**(1/k),(k,1,n))',
    'product(k,(k,-2,2))',
    'product(k+Rational(1,2),(k,0,n))',
    'product(4**k,(k,1,n))',
    'product(Rational(1,2)**k,(k,1,n))',
    'product(2**k*x,(k,1,n))',
    'product(2**(k+1),(k,1,n))',
    'product(2**(2*k),(k,1,n))',
    'product(exp(k)*exp(x),(k,1,n))',
    'product(2**k/4,(k,1,n))',
    'product(3*k/4,(k,1,n))',
    'product(2**k*3**k*x,(k,1,n))',
    'product(p**k,(k,1,n))',
    'product(2**(k/2),(k,1,n))',
    'product(2**(-k),(k,1,n))',
    'product((k+1)/k,(k,1,n))',
    'product(1-x/k,(k,1,n))',
    'product(x/k,(k,1,n))',
    'product(k/x,(k,1,n))',
    'product(k**2/(k+1)**2,(k,1,n))',
    'product((2*k+2)/(k+3),(k,1,n))',
    'product(k*(x+1),(k,1,n))',
    'product(x**k*y**k,(k,1,n))',
    'product(sqrt(2)*k,(k,1,n))',
    'product(k*exp(k),(k,1,n))',
    'product(1/(2**k),(k,1,n))',
    'product(exp(1/k),(k,1,n))',
    'product(k*2**k,(k,1,n))',
    'product(2**k/k,(k,1,n))',
]


@pytest.mark.parametrize('source', PRODUCTS)
def test_product_matches_sympy(source):
    assert outcome(lambda: eval(source, namespace(sm))) == outcome(lambda: eval(source, namespace(sp)))


def node_cases(m):
    k, n, x, y, a, b, j = m.symbols('k n x y a b j')
    f = m.Function('f')
    P = m.Product
    return [
     ("str", lambda: str(P(k,(k,1,n)))),
     ("srepr", lambda: m.srepr(P(k**2,(k,1,n)))),
     ("srepr2", lambda: m.srepr(P(x*k,(k,1,n),(n,1,y)))),
     ("latex1", lambda: m.latex(P(k,(k,1,n)))),
     ("latex2", lambda: m.latex(P(k+x,(k,1,n)))),
     ("latex3", lambda: m.latex(P(k*j,(k,1,n),(j,1,y)))),
     ("latex4", lambda: m.latex(P(k,(k,1,n))**2)),
     ("latex5", lambda: m.latex(x*P(k**2,(k,1,n)))),
     ("pretty1", lambda: m.pretty(P(k,(k,1,n)), use_unicode=True)),
     ("pretty2", lambda: m.pretty(P(k,(k,1,n)), use_unicode=False)),
     ("pretty3", lambda: m.pretty(P(1/k,(k,1,n)), use_unicode=True)),
     ("pretty4", lambda: m.pretty(P(k**2+x,(k,1,n),(n,1,y)), use_unicode=True)),
     ("pretty5", lambda: m.pretty(x*P(k,(k,1,n)), use_unicode=True)),
     ("pretty6", lambda: m.pretty(P(f(k),(k,a,b))**2, use_unicode=True)),
     ("free", lambda: sorted(map(str,P(k*x,(k,1,n)).free_symbols))),
     ("free2", lambda: sorted(map(str,P(k*x,(k,1,n),(n,1,y)).free_symbols))),
     ("subs", lambda: str(P(k*x,(k,1,n)).subs(x,2))),
     ("subs_k", lambda: str(P(k*x,(k,1,n)).subs(k,2))),
     ("subs_n", lambda: str(P(k*x,(k,1,n)).subs(n,3))),
     ("subs_n_doit", lambda: str(P(k*x,(k,1,n)).subs(n,3).doit())),
     ("xreplace", lambda: str(P(k*x,(k,1,n)).xreplace({k:y}))),
     ("doit", lambda: str(P(k,(k,1,n)).doit())),
     ("doit_nested", lambda: str(P(P(k,(k,1,j)),(j,1,n)).doit())),
     ("doit_unev", lambda: str(P(f(k),(k,1,n)).doit())),
     ("diff_free", lambda: str(P(k,(k,1,n)).diff(x))),
     ("diff_n", lambda: str(P(k,(k,1,n)).diff(n))),
     ("diff_x_doit", lambda: str(m.simplify(P(x+k,(k,1,3)).diff(x).doit()))),
     ("bound", lambda: str(P(k,(k,1,n)).bound_symbols)),
     ("variables", lambda: str(P(k,(k,1,n),(n,1,3)).variables)),
     ("limits", lambda: str(P(k,(k,1,n)).limits)),
     ("function", lambda: str(P(k*x,(k,1,n)).function)),
     ("func", lambda: str(P(k,(k,1,n)).func(*P(k,(k,1,n)).args))),
     ("isinst", lambda: (isinstance(P(k,(k,1,n)),P), isinstance(m.Sum(k,(k,1,n)),P), isinstance(k,P))),
     ("eq", lambda: P(k,(k,1,n))==P(j,(j,1,n))),
     ("hash_eq", lambda: P(k,(k,1,n))==P(k,(k,1,n))),
     ("sympify", lambda: str(m.sympify("Product(k**2, (k, 1, n))"))),
     ("sympify_doit", lambda: str(m.sympify("Product(k**2, (k, 1, n))").doit())),
     ("is_number", lambda: (P(k,(k,1,3)).is_number, P(k,(k,1,n)).is_number)),
     ("evalf", lambda: str(P(k,(k,1,5)).evalf())),
     ("n", lambda: str(P(1+1/k,(k,1,10)).n())),
     ("mul", lambda: str(2*P(k,(k,1,n))*P(k,(k,1,n)))),
     ("noargs", lambda: str(P(k,(k,1,n)).args)),
     ("product_func", lambda: str(m.product(k,(k,1,n)))),
     ("has", lambda: (P(k,(k,1,n)).has(k), P(k,(k,1,n)).has(n))),
     ("as_dummy", lambda: str(P(k*x,(k,1,n)).subs(x,k))),
     ("doit_deep", lambda: str(P(m.Integral(x,(x,0,k)),(k,1,n)).doit())),
     ("doit_nodeep", lambda: str(P(m.Integral(x,(x,0,k)),(k,1,n)).doit(deep=False))),
     ("expand", lambda: str(m.expand(P(k*(x+1),(k,1,n))))),
     ("srepr_nested", lambda: m.srepr(P(P(k,(k,1,j)),(j,1,n)))),
     ("reverse", lambda: str(P(k,(k,5,2)).doit())),
     ("reverse_sym", lambda: str(P(1/k,(k,n+2,n)).doit())),
     ("float", lambda: str(m.product(k*m.Float(0.5),(k,1,n)))),
     ("nan", lambda: str(m.product(m.nan,(k,1,n)))),
     ("zoo", lambda: str(m.product(m.oo,(k,1,3)))),
     ("posint", lambda: str(m.product(k,(k,1,m.Symbol('p',positive=True,integer=True))))),
    ]


@pytest.mark.parametrize('index', range(len(node_cases(sm))))
def test_product_node_matches_sympy(index):
    name, actual = node_cases(sm)[index]
    _, expected = node_cases(sp)[index]
    assert outcome(actual) == outcome(expected), name


CONVERGENCE = [
    "P(1+1/k**2,(k,1,oo))", "P(1+1/k,(k,1,oo))", "P(1-1/k**2,(k,2,oo))", "P(exp(1/k**2),(k,1,oo))",
    "P(1/k,(k,1,oo))", "P(2,(k,1,oo))", "P(1,(k,1,oo))", "P(k,(k,1,oo))", "P(Rational(1,2),(k,1,oo))",
    "P(1+1/k**3,(k,1,oo))", "P(k/(k+1),(k,1,oo))", "P(1+1/2**k,(k,1,oo))", "P(1+1/k**2,(k,-oo,-1))",
    "P(1+1/k**2,(k,-oo,oo))", "P(1+1/k,(k,-oo,-1))", "P(1+1/factorial(k),(k,1,oo))", "P(1+2**k,(k,1,oo))",
    "P(2**k,(k,1,oo))", "P(Rational(1,2)**k,(k,1,oo))", "P(1+k/2**k,(k,1,oo))", "P(1-1/2**k,(k,1,oo))",
    "P(1+3**k/factorial(k),(k,1,oo))", "P(1+factorial(k)/3**k,(k,1,oo))", "P(1+2**(-k)*k**3,(k,1,oo))",
    "P(k,(k,1,5))", "P(1+x/k**2,(k,1,oo))",
    "S(2**-k,(k,1,oo))", "S(1/factorial(k),(k,1,oo))", "S(2**k,(k,1,oo))", "S(1/k**2,(k,-oo,-1))",
    "S(k/2**k,(k,1,oo))", "S(1/k**2,(k,-oo,oo))", "S(factorial(k)/3**k,(k,1,oo))", "S(1/k,(k,-oo,-1))",
    "S(k**2/2**k,(k,1,oo))", "S(1/k**2,(k,1,oo))", "S(1/k,(k,1,oo))", "S(k/(k**3+1),(k,1,oo))",
]


@pytest.mark.parametrize('source', CONVERGENCE)
def test_is_convergent_matches_sympy(source):
    def run(m):
        k, x = m.symbols('k x')
        names = dict(P=m.Product, S=m.Sum, k=k, x=x, oo=m.oo, exp=m.exp, Rational=m.Rational,
                     factorial=m.factorial)
        return outcome(lambda: eval(source, names).is_convergent())
    assert run(sm) == run(sp)


def test_unsupported_convergence_raises_not_implemented():
    k, x = sm.symbols('k x')
    with pytest.raises(NotImplementedError):
        sm.Sum(1/k**x, (k, 1, sm.oo)).is_convergent()
    with pytest.raises(NotImplementedError):
        sm.Product(1 + 1/k**x, (k, 1, sm.oo)).is_convergent()


# Non-rational terms: divergence, p-series and logarithmic comparison,
# alternating, ratio, Raabe and root tests, bounded factors.
SERIES = """S((-1)**k/k,(k,1,oo))|S((-1)**k/k**2,(k,1,oo))|S((-1)**k,(k,1,oo))|S(1/sqrt(k),(k,1,oo))
S(1/k**Rational(3,2),(k,1,oo))|S(1/(k*log(k)),(k,2,oo))|S(1/(k*log(k)**2),(k,2,oo))|S(log(k)/k**2,(k,1,oo))
S(log(k)/k,(k,1,oo))|S(sin(k)/k**2,(k,1,oo))|S(cos(k)/k**2,(k,1,oo))|S(sin(1/k),(k,1,oo))
S(sin(1/k)**2,(k,1,oo))|S(exp(-k),(k,1,oo))|S(exp(-k**2),(k,1,oo))|S(k*exp(-k),(k,1,oo))
S(1/factorial(k),(k,0,oo))|S(k**k/factorial(k),(k,1,oo))|S(factorial(k)/k**k,(k,1,oo))|S(1/(k**2+1),(k,-oo,oo))
S(3**k/2**k,(k,0,oo))|S((-1)**(k+1)/sqrt(k),(k,1,oo))|S(1/log(k),(k,2,oo))|S(log(1+1/k),(k,1,oo))
S(log(1+1/k**2),(k,1,oo))|S(atan(k)/k**2,(k,1,oo))|S(1/(k*log(k)*log(log(k))),(k,3,oo))|S(1/k**k,(k,1,oo))
S(binomial(2*k,k)/4**k,(k,1,oo))|S(cos(pi*k)/k,(k,1,oo))|S(sin(pi*k),(k,1,oo))|S(exp(1/k)-1,(k,1,oo))
S(sqrt(k+1)-sqrt(k),(k,1,oo))|S(k**3/3**k,(k,1,oo))|S(tan(1/k),(k,1,oo))|S(1/k**2+1/2**k,(k,1,oo))
S(2**k/k**2,(k,1,oo))|S(exp(k)/k,(k,1,oo))|S(5,(k,1,oo))|S(0,(k,1,oo))|S(1/(k**2-2),(k,2,oo))
S(1/(k+1)/log(k+1)**2,(k,1,oo))|S((-2)**k/k,(k,1,oo))|S((-1)**k*k/(k**2+1),(k,1,oo))|S(1/(2**k+k),(k,1,oo))
S(k/2**k,(k,-oo,-1))|S(gamma(k)/gamma(k+Rational(1,2)),(k,1,oo))|S(1/sqrt(k*(k+1)),(k,1,oo))
S(asin(1/k),(k,1,oo))|S(1/(k**2*log(k)),(k,2,oo))|S(Piecewise((1/k**2,k>5),(1,True)),(k,1,oo))
S(1/(k*log(k)*log(log(k))**2),(k,3,oo))|S(3**k*factorial(k)/k**k,(k,1,oo))|S(1/(k*sqrt(log(k))),(k,2,oo))
P(cos(pi/k),(k,1,oo))|P(cos(1/k),(k,1,oo))|P(exp(-k**2),(k,1,oo))|P(exp(1/k),(k,1,oo))
P(1+exp(-k),(k,1,oo))|P(1+sin(1/k**2),(k,1,oo))|P(sqrt(1+1/k**2),(k,1,oo))|P(1+1/(2**k+k),(k,1,oo))"""
SERIES = SERIES.replace('\n', '|').split('|')

ABSOLUTE = """S(1/k**2,(k,1,oo))|S((-1)**k/k,(k,1,oo))|S((-1)**k/k**2,(k,1,oo))|S((-1)**k/factorial(k),(k,0,oo))
S(1/factorial(k),(k,0,oo))|S(factorial(k)/k**k,(k,1,oo))|S(k**k/factorial(k),(k,1,oo))
S(binomial(2*k,k)/4**k,(k,1,oo))|S((-2)**k/k,(k,1,oo))|S((-1)**(k+1)/sqrt(k),(k,1,oo))
S(cos(pi*k)/k**2,(k,1,oo))|S(1/(k**2+1),(k,-oo,oo))|S(exp(-k),(k,1,oo))|S(1/k,(k,1,10))"""
ABSOLUTE = ABSOLUTE.replace('\n', '|').split('|')


def series_namespace(m):
    names = "Sum Product oo sqrt log sin cos exp factorial Rational atan binomial pi tan gamma asin Piecewise"
    d = {name: getattr(m, name) for name in names.split()}
    d.update(S=m.Sum, P=m.Product, k=m.Symbol('k'))
    return d


@pytest.mark.parametrize('source', SERIES)
def test_series_convergence_matches_sympy(source):
    def run(m):
        return outcome(lambda: eval(source, series_namespace(m)).is_convergent())
    assert run(sm) == run(sp)


@pytest.mark.parametrize('source', ABSOLUTE)
def test_absolute_convergence_matches_sympy(source):
    def run(m):
        return outcome(lambda: eval(source, series_namespace(m)).is_absolutely_convergent())
    assert run(sm) == run(sp)


@pytest.mark.parametrize('source,method', [
    # SymPy answers False (its comparison through O() drops the cosine
    # cancellation); 1 - cos(1/k) ~ 1/(2*k**2) converges
    ("S(1-cos(1/k),(k,1,oo))", 'is_convergent'),
    # SymPy raises: |sin(k)|/k**2 <= 1/k**2
    ("S(Abs(sin(k))/k**2,(k,1,oo))", 'is_convergent'),
    ("S(sin(k)/k**2,(k,1,oo))", 'is_absolutely_convergent'),
])
def test_convergence_beyond_sympy(source, method):
    names = series_namespace(sm)
    names['Abs'] = sm.Abs
    assert getattr(eval(source, names), method)() == sm.true


def test_absolutely_convergent_is_sum_only():
    k = sm.Symbol('k')
    with pytest.raises(AttributeError):
        sm.Product(1 + 1/k**2, (k, 1, sm.oo)).is_absolutely_convergent()


def test_abs_of_integer_powers_matches_sympy():
    def cases(m):
        n = m.Symbol('n', integer=True)
        x = m.Symbol('x', real=True)
        return [str(e) for e in (m.Abs((-1)**n), m.Abs((-1)**n*x), m.Abs(x**(2*n)),
                                 m.Abs((-1)**(n + 1)/m.sqrt(n**2 + 1)), m.Abs((-2)**n))]
    assert cases(sm) == cases(sp)


LIMIT_METHODS = """S(k, (k, 1, n)).reverse_order(0)|S(k, (k, 1, n)).reverse_order(k)|P(k, (k, 1, n)).reverse_order(0)
S(k*j, (k, 1, n), (j, 1, m)).reverse_order(0, 1)|P(k*j, (k, 1, n), (j, 1, m)).reverse_order(0, 1)
P(k*j, (k, 1, n), (j, 1, m)).reverse_order(j)|S(k, (k, 1, n)).change_index(k, k + 1)
S(k, (k, 1, n)).change_index(k, -k)|S(k, (k, 1, n)).change_index(k, -k + 3, j)
S(k**2, (k, 0, n)).change_index(k, 2*k)|S(k, (k, 1, n)).change_index(k, k**2)|P(k, (k, 1, n)).change_index(k, k - 2)
S(k*j, (k, 1, n), (j, 1, m)).reorder((k, j))|S(k*j, (k, 1, n), (j, 1, m)).reorder((0, 1))
S(k*j*x, (k, 1, n), (j, 1, m), (x, 1, 3)).reorder((k, x), (j, k))|S(k*j, (k, 1, n), (j, 1, m)).reorder_limit(0, 1)
S(k*j, (k, 1, j), (j, 1, m)).reorder_limit(0, 1)|S(k*j, (k, 1, j), (j, 1, m)).reorder((k, j))
S(k*j, (k, 1, n), (j, 1, m)).index(j)|S(k*j, (k, 1, n), (j, 1, m)).index(x)|P(k, (k, 1, n)).index(k)
S(k, (k, 1, n)).has_finite_limits|S(k, (k, 1, oo)).has_finite_limits|S(k, (k, -oo, 3)).has_finite_limits
S(k, (k, 1, 5)).has_finite_limits|P(k, (k, 1, n)).has_finite_limits|I(x, (x, 1, 2)).has_finite_limits
I(x, x).has_finite_limits|S(k, (k, 1, n)).has_reversed_limits|S(k, (k, 5, 1)).has_reversed_limits
S(k, (k, 1, 5)).has_reversed_limits|S(k, (k, 1, oo)).has_reversed_limits|S(k, (k, p, p - 3)).has_reversed_limits
P(k, (k, 1, -oo)).has_reversed_limits|I(x, (x, 2, 1)).has_reversed_limits|I(x, x).has_reversed_limits
S(k, (k, 1, 0)).has_empty_sequence|S(k, (k, 1, 1)).has_empty_sequence|S(k, (k, n + 1, n)).has_empty_sequence
S(k, (k, 1, n)).has_empty_sequence|P(k, (k, 1, 0)).has_empty_sequence|S(k*j, (k, 1, 0), (j, 1, n)).has_empty_sequence
P(k**2, (k, 1, n)).term|S(k, (k, 1, n)).change_index(k, k + 1).doit()|S(k, (k, 1, n)).reverse_order(0).doit()
P(k, (k, 1, 5)).reverse_order(0).doit()|P(k, (k, 1, n)).reverse_order(0).doit()"""
LIMIT_METHODS = LIMIT_METHODS.replace('\n', '|').split('|')


@pytest.mark.parametrize('source', LIMIT_METHODS)
def test_limit_methods_match_sympy(source):
    def run(m):
        k, j, n, m_, x = m.symbols('k j n m x')
        names = dict(S=m.Sum, P=m.Product, I=m.Integral, k=k, j=j, n=n, m=m_, x=x, oo=m.oo,
                     p=m.Symbol('p', positive=True))
        try:
            return str(eval(source, names))
        except Exception as error:
            return 'EXC %s %s' % (type(error).__name__, error)
    assert run(sm) == run(sp)


def test_reorder_error_is_not_implemented_error():
    k, j, m = sm.symbols('k j m')
    assert issubclass(sm.ReorderError, NotImplementedError)
    with pytest.raises(sm.ReorderError):
        sm.Sum(k*j, (k, 1, j), (j, 1, m)).reorder_limit(0, 1)


ASSUMPTIONS = [
    ("P(kp,(kp,1,3))", "is_positive"), ("P(kp,(kp,1,n))", "is_positive"), ("P(kp,(kp,1,n))", "is_zero"),
    ("P(kp,(kp,1,n))", "is_nonnegative"), ("P(kq,(kq,1,n))", "is_nonnegative"), ("P(kq,(kq,1,n))", "is_positive"),
    ("P(k,(k,1,n))", "is_positive"), ("P(0,(k,1,n))", "is_zero"), ("P(k,(k,1,0))", "is_positive"),
    ("P(kp,(kp,1,oo))", "is_positive"), ("P(kr,(kr,1,4))", "is_real"), ("P(k,(k,1,n))", "is_commutative"),
    ("P(k,(k,1,n))", "is_number"), ("P(k,(k,1,3))", "is_number"),
]


@pytest.mark.parametrize('source,fact', ASSUMPTIONS)
def test_product_assumptions_match_sympy(source, fact):
    def run(m):
        k, n = m.symbols('k n')
        names = dict(P=m.Product, k=k, n=n, oo=m.oo, kp=m.Symbol('kp', positive=True, integer=True),
                     kr=m.Symbol("kr", real=True), kq=m.Symbol("kq", nonnegative=True, integer=True))
        return getattr(eval(source, names), fact)
    assert run(sm) == run(sp)


@pytest.mark.parametrize('source', [
    "Product(nan,(k,1,n)).doit()", "Product(oo,(k,1,n)).doit()", "Product(0,(k,1,n)).doit()",
    "Product(k,(k,1,oo)).doit()", "Product(1+1/k,(k,1,oo)).doit()", "Product(Float('1.5'),(k,1,3)).doit()",
    "Product(k*Float('0.5'),(k,1,4)).doit()", "Product(k,(k,1,5)).evalf(30)", "product(x,(k,1,0))",
    "Product(k,(k,1,n)).subs(n,oo)", "Product(k*x,(k,1,n)).diff(n)", "Product(k,(k,1,n)).doit().subs(n,5)",
    "product(k,(k,1,zoo))", "product(nan,(k,1,3))", "product(2,(k,1,-oo))",
])
def test_product_special_values_match_sympy(source):
    def run(m):
        k, n, x = m.symbols('k n x')
        names = dict(Product=m.Product, product=m.product, k=k, n=n, x=x, oo=m.oo, nan=m.nan, zoo=m.zoo,
                     Float=m.Float)
        return outcome(lambda: eval(source, names))
    assert run(sm) == run(sp)


def test_product_derivative_uses_a_dummy_index():
    k, n, x = sm.symbols('k n x')
    derivative = sm.Product(x*k, (k, 1, n)).diff(x)
    assert isinstance(derivative, sm.Sum)
    (index, low, high), = derivative.limits
    assert (low, high) == (1, n) and index.is_Dummy
    assert derivative.function == index*sm.Product(k*x, (k, 1, index - 1))*sm.Product(k*x, (k, index + 1, n))


@pytest.mark.parametrize('limits', [(1,), 'k'])
def test_malformed_product_limits_raise_value_error(limits):
    k = sm.Symbol('k')
    with pytest.raises(ValueError):
        sm.Product(k, (k, *limits) if isinstance(limits, tuple) else k)


@pytest.mark.parametrize('source', [
    "Product(k**2,(k,1,n))", "Product(k*x,(k,1,n),(n,1,y))", "Product(Product(k,(k,1,j)),(j,1,n))",
    "2*Product(k,(k,a,b))**2", "Product(f(k),(k,1,n))",
])
def test_product_roundtrips(source):
    def build(m):
        k, n, x, y, j, a, b = m.symbols('k n x y j a b')
        return eval(source, dict(Product=m.Product, k=k, n=n, x=x, y=y, j=j, a=a, b=b, f=m.Function('f')))
    actual, expected = build(sm), build(sp)
    assert str(actual) == str(expected)
    assert actual.srepr() == sp.srepr(expected)
    assert sm.to_sympy(actual) == expected and sm.from_sympy(expected) == actual
    assert sm.sympify(str(actual)) == actual and sm.sympify(actual.srepr()) == actual
    restored = pickle.loads(pickle.dumps(actual))
    assert restored == actual and restored.args == actual.args
    assert sm.default_sort_key(actual) == sp.default_sort_key(expected)
