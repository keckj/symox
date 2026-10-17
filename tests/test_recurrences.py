"""Linear recurrences: rsolve, rsolve_poly, rsolve_ratio and rsolve_hyper."""
import pickle

import pytest
import sympy as sp
import symox as sm
from sympy.solvers.recurr import rsolve_hyper, rsolve_poly, rsolve_ratio

SOLVERS = {sp: dict(rsolve=sp.rsolve, rsolve_poly=rsolve_poly, rsolve_ratio=rsolve_ratio, rsolve_hyper=rsolve_hyper),
           sm: dict(rsolve=sm.rsolve, rsolve_poly=sm.rsolve_poly, rsolve_ratio=sm.rsolve_ratio,
                    rsolve_hyper=sm.rsolve_hyper)}
FUNCTIONS = "Eq Function Rational Symbol symbols sqrt sin cos exp factorial binomial oo".split()


def namespace(m, integer):
    d = {name: getattr(m, name) for name in FUNCTIONS}
    d.update(SOLVERS[m])
    if integer:
        n, k = m.symbols('n k', integer=True)
        d.update(a=m.Symbol('a'), A=m.Function('a'), xf=m.Function('x'), f=m.Function('f'), m=m.Symbol('m'),
                 rt=lambda c, p: sum(c[i]*p.subs(n, n + i) for i in range(len(c))))
    else:
        n, k = m.symbols('n k')
        d.update(x=m.Symbol('x'), a=m.Function('a'))
    d.update(n=n, k=k, b=m.Symbol('b'), y=m.Function('y'))
    return d


def outcome(source, m, integer):
    try:
        return str(eval(source, namespace(m, integer)))
    except Exception as error:
        return 'EXC %s: %s' % (type(error).__name__, error)


GENERIC = [
    'rsolve(y(n+1)-2*y(n),y(n))',
    'rsolve(y(n+1)-2*y(n),y(n),{y(0):3})',
    'rsolve(y(n+2)-y(n+1)-y(n),y(n))',
    'rsolve(y(n+2)-y(n+1)-y(n),y(n),{y(0):0,y(1):1})',
    'rsolve(y(n+1)-y(n)-1,y(n))',
    'rsolve(y(n+1)-y(n)-n,y(n))',
    'rsolve(y(n+1)-y(n)-2**n,y(n))',
    'rsolve(y(n+1)-2*y(n)-1,y(n),{y(0):0})',
    'rsolve(y(n+2)-4*y(n+1)+4*y(n),y(n))',
    'rsolve(y(n+2)+y(n),y(n))',
    'rsolve(y(n+1)-(n+1)*y(n),y(n))',
    'rsolve(y(n+1)-(n+1)*y(n),y(n),{y(0):1})',
    'rsolve(y(n)-y(n-1)-n,y(n))',
    'rsolve(y(n+2)-5*y(n+1)+6*y(n),y(n),[1,2])',
    'rsolve(y(n+3)-y(n),y(n))',
    'rsolve(y(n+1)-n*y(n),y(n))',
    'rsolve(Eq(y(n+1),3*y(n)+2**n),y(n))',
    'rsolve(y(n+1)-y(n)-1/n,y(n))',
    'rsolve(y(n+1)-y(n)*x,y(n))',
    'rsolve(y(n+2)-2*y(n+1)+y(n),y(n))',
    'rsolve(y(n+2)-2*y(n+1)+y(n)-1,y(n))',
    'rsolve(n*y(n+1)-(n+1)*y(n),y(n))',
    'rsolve((n+1)*y(n+1)-y(n),y(n))',
    'rsolve(y(n+1)+y(n),y(n),{y(0):1})',
    'rsolve(y(n+1)-y(n)-n*2**n,y(n))',
    'rsolve(y(n+2)-y(n)-1,y(n))',
    'rsolve(y(n)-2*y(n-1)-3*y(n-2),y(n),{y(0):1,y(1):1})',
    'rsolve(y(n+1)-y(n)-sin(n),y(n))',
    'rsolve(y(n+1)**2-y(n),y(n))',
    'rsolve_poly([-1,1],n,n)',
    'rsolve_poly([-1,1],1,n)',
    'rsolve_poly([1,-2,1],n,n)',
    'rsolve_poly([-n-1,n],1,n)',
    'rsolve_ratio([-2*n**3+n**2+2*n-1, 2*n**3+n**2-6*n, -2*n**3-11*n**2-18*n-9, 2*n**3+13*n**2+22*n+8], 0, n)',
    'rsolve_hyper([-1,1],2**n,n)',
    'rsolve_hyper([-2,1],0,n)',
    'rsolve_hyper([-1,-1,1],0,n)',
    'rsolve_hyper([-n-1,1],0,n)',
    'rsolve_hyper([-1,1],n*2**n,n)',
    'rsolve_hyper([2*n*(n+1),-n**2-3*n+2,n-1],0,n)',
    'rsolve(y(n+2)-(n+1)*y(n),y(n))',
    'rsolve(y(n+1)-2*y(n)-n,y(n),{y(1):1})',
    'rsolve(y(n+2)+4*y(n),y(n))',
    'rsolve(y(n+1)-y(n)-factorial(n),y(n))',
    'rsolve(y(n+1)-y(n)-binomial(n,2),y(n))',
]

INTEGER = [
    'rsolve_poly([-1,-1,1],0,n)',
    'rsolve_poly([-1,-1,1],1,n)',
    'rsolve_poly([-1,n+1],n,n)',
    'rsolve_poly([-1,1],n,n)',
    'rsolve_poly([-n-1,n],1,n)',
    'rsolve_poly([-4*n-2,1],4*n+1,n)',
    'rsolve_poly([-1,1],n**5+n**3,n)',
    'rsolve_poly([-1,1],n,n,symbols=True)',
    'rsolve_poly([-1,1],n,n,3)',
    'rsolve_poly([1,-2,1],0,n,symbols=True)',
    'rsolve_poly([-1,1],1/n,n)',
    'rsolve_poly([n**2-n+12,1],rt([n**2-n+12,1],n**4),n)',
    'rsolve_poly([-1,1,1,-1,1],rt([-1,1,1,-1,1],n+n**2),n)',
    'rsolve_poly([-n,1],rt([-n,1],n**3),n)',
    'rsolve_poly([-2,-1,1],rt([-2,-1,1],27*n+52*n**2-3*n**3+12*n**4-52*n**5),n)',
    'rsolve_ratio([-2*n**3+n**2+2*n-1,2*n**3+n**2-6*n,-2*n**3-11*n**2-18*n-9,2*n**3+13*n**2+22*n+8],0,n,symbols=True)',
    'rsolve_ratio([-n,n+2],n,n)',
    'rsolve_ratio([-1,1],n,n)',
    'rsolve_ratio([-1,n+1],0,n)',
    'rsolve_ratio([n,-(n+1)],0,n)',
    'rsolve_hyper([-1,-1,1],0,n)',
    'rsolve_hyper([n**2-2,-2*n-1,1],0,n)',
    'rsolve_hyper([n**2-k,-2*n-1,1],0,n)',
    'rsolve_hyper([n+2,-(2*n+3)*(17*n**2+51*n+39),n+1],0,n)',
    'rsolve_hyper([-n-1,-1,1],0,n)',
    'rsolve_hyper([-1,1],n,n)',
    'rsolve_hyper([-1,1],1+n,n)',
    'rsolve_hyper([-1,1],3*(n+n**2),n)',
    'rsolve_hyper([-a,1],0,n)',
    'rsolve_hyper([-a,0,1],0,n)',
    'rsolve_hyper([1,1,1],0,n)',
    'rsolve_hyper([1,-2*n/a-2/a,1],0,n)',
    'rsolve_hyper([n**2-n+12,1],n*(n**2-n+12)+n+1,n)',
    'rsolve_hyper([-1,1],2**n,n,symbols=True)',
    'rsolve_hyper([-1,1],sin(n),n)',
    'rsolve_hyper([-2,1],rt([-2,1],n**2),n)',
    'rsolve_hyper([-2,-1,1],rt([-2,-1,1],n**3),n)',
    'rsolve_hyper([-n,1],rt([-n,1],n+1),n)',
    'rsolve(y(n+2)-y(n+1)-y(n),y(n),[0,5])',
    'rsolve(y(n+2)-y(n+1)-y(n),y(n),{0:0,1:5})',
    'rsolve(y(n)-y(n-1)-y(n-2),y(n),[0,5])',
    'rsolve(Eq(y(n),y(n-1)+y(n-2)),y(n),[0,5])',
    'rsolve((n-1)*y(n+2)-(n**2+3*n-2)*y(n+1)+2*n*(n+1)*y(n),y(n))',
    'rsolve((n-1)*y(n+2)-(n**2+3*n-2)*y(n+1)+2*n*(n+1)*y(n),y(n),[0,3])',
    'rsolve((n-1)*y(n+2)-(n**2+3*n-2)*y(n+1)+2*n*(n+1)*y(n),y(n),[])',
    'rsolve((n-1)*y(n+2)-(n**2+3*n-2)*y(n+1)+2*n*(n+1)*y(n),y(n),{})',
    'rsolve(y(n)-y(n-1)-2,y(n),{y(0):0})',
    'rsolve(y(n)-y(n-1)-2,y(n),{y(0):1})',
    'rsolve(y(n)-y(n-1)-2,y(n),{y(0):0,y(1):1})',
    'rsolve(3*y(n-1)-y(n)-1,y(n),{y(0):0})',
    'rsolve(3*y(n-1)-y(n)-1,y(n),{y(0):2})',
    'rsolve(y(n)-1/n*y(n-1),y(n))',
    'rsolve(y(n)-1/n*y(n-1)-1,y(n))',
    'rsolve(2*y(n-1)+(1-n)*y(n)/n,y(n),{y(1):1})',
    'rsolve(2*y(n-1)+(1-n)*y(n)/n,y(n),{y(1):3})',
    'rsolve((n-1)*(n-2)*y(n+2)-(n+1)*(n+2)*y(n),y(n),{y(3):6,y(4):24})',
    'rsolve((n-1)*(n-2)*y(n+2)-(n+1)*(n+2)*y(n),y(n),{y(3):6,y(4):-24})',
    'rsolve(Eq(y(n+1),a*y(n)),y(n),{y(1):a})',
    'rsolve(y(n)+a*(y(n+1)+y(n-1))/2,y(n))',
    'rsolve((k+1)*y(k),y(k))',
    'rsolve(y(n)+y(n+1)+2**n+3**n,y(n))',
    'rsolve(y(n)-y(k+1),y(n))',
    'rsolve(y(n)-y(n+1),xf(n))',
    'rsolve(y(n)-xf(n+1),y(n))',
    'rsolve(y(n)-sqrt(n)*y(n+1),y(n))',
    'rsolve(y(n)-y(n+1),y(n),{xf(0):0})',
    'rsolve(y(n)+y(n+1)+2**n+cos(n),y(n))',
    'rsolve(y(n+2)-y(n+1)+y(n)/4,y(n))',
    'rsolve(y(n+2)-y(n+1)+y(n)/4,y(n),{y(0):0,y(1):1})',
    'rsolve(y(n+3)-y(n+2)-y(n+1)+y(n),y(n))',
    'rsolve(y(n+3)+3*y(n+2)+3*y(n+1)+y(n),y(n))',
    'rsolve(y(n)-2*y(n-3)+5*y(n-2)-4*y(n-1),y(n),[1,3,8])',
    'rsolve(y(n+3)-3*y(n+1)+2*y(n),y(n),{y(1):0,y(2):8,y(3):-2})',
    'rsolve(A(n+3)-A(n+2)-A(n+1)+A(n),A(n))',
    'rsolve(A(n)-2*A(n-1)-n,A(n),{A(0):1})',
    'rsolve(y(n)-y(n-1)-2*y(n-2)-2*n,y(n))',
    'rsolve(y(n)-y(n-1)-2*y(n-2)-2*n,y(n),{y(0):-1,y(1):1})',
    'rsolve(-2*y(n)+y(n+1)+n-1,y(n))',
    'rsolve(Eq(f(n),2*f(n-1)+n),f(n))',
    'rsolve(Eq(f(n+1),2*f(n)+n**2+1),f(n))',
    'rsolve(Eq(f(n+1),2*f(n)+n**2+1),f(n),{f(1):0})',
    'rsolve(Eq(f(n),2*f(n-1)+3*n**2),f(n))',
    'rsolve(Eq(f(n),2*f(n-1)+n**2),f(n))',
    'rsolve(Eq(f(n),2*f(n-1)+n),f(n),{f(0):1})',
    'rsolve(y(n+1)-2*y(n)-n*3**n,y(n))',
    'rsolve(y(n+2)-y(n)-2**n,y(n))',
    'rsolve(y(n+1)-y(n)-2**n-3**n,y(n))',
    'rsolve(y(n+1)-m*y(n),y(n))',
    'rsolve(y(n+1)-m*y(n),y(n),{y(0):1})',
    'rsolve(y(n+1)-y(n)-m,y(n))',
    'rsolve(y(n+1)-y(n),y(n))',
    'rsolve(y(n),y(n))',
    'rsolve(2*y(n),y(n))',
    'rsolve(y(n+1)-y(n)-2**(n**2),y(n))',
    'rsolve(y(n+1)-sin(n)*y(n),y(n))',
    'rsolve(y(n+1)-y(n)-1,y(n),{y(0):a})',
    'rsolve(y(n+1)-y(n)-1,y(n),{n:0})',
    'rsolve(y(n+1)-y(n)-1,y(n),{y(0):oo})',
    'rsolve(y(n+1)-(n+1)*y(n),y(n),{y(0):1})',
    'rsolve(y(n+1)-y(n)/2,y(n))',
    'rsolve(y(n+1)-y(n)/2,y(n),{y(0):4})',
    'rsolve(2*y(n+1)-y(n)-1,y(n))',
    'rsolve(y(n+1)+y(n)-n,y(n))',
    'rsolve(y(n+1)-y(n)-Rational(1,2),y(n))',
    'rsolve(y(n+1)-y(n)-1.5,y(n))',
    'rsolve(y(n+1)-2.0*y(n),y(n))',
    'rsolve(y(n+2)-2*y(n+1)+2*y(n),y(n))',
    'rsolve(y(n+2)-6*y(n+1)+9*y(n),y(n))',
    'rsolve(y(n+1)-y(n)-factorial(n)*n,y(n))',
    'rsolve(y(n+1)-(n+2)*y(n)-factorial(n+1),y(n))',
    'rsolve(y(n+2)-3*y(n+1)+2*y(n)-1,y(n))',
    'rsolve(y(n+1)-y(n)-1/(n*(n+1)),y(n))',
    'rsolve(n*y(n+1)-y(n),y(n))',
]


@pytest.mark.parametrize('source', GENERIC)
def test_generic_index_matches_sympy(source):
    assert outcome(source, sm, False) == outcome(source, sp, False)


@pytest.mark.parametrize('source', INTEGER)
def test_integer_index_matches_sympy(source):
    assert outcome(source, sm, True) == outcome(source, sp, True)


@pytest.mark.parametrize('source', [
    "rsolve(y(n+1)-2*y(n),y(n),{y(0):Rational(1,3)})", "rsolve(y(n+1)-y(n)/2,y(n),{y(0):0})",
    "rsolve(y(n+1)-2*y(n),y(n),{0:oo})", "rsolve(y(n+2)-y(n),y(n),(1,2))",
    "rsolve(y(n+1)-p*y(n),y(n))", "rsolve(y(p+1)-y(p)-p,y(p))", "rsolve(y(n+1)-2*y(n)-q,y(n),{y(0):q})",
    "rsolve(y(n+1)-y(n)-2.5,y(n))", "rsolve(y(n+1)-y(n)-n*2,y(n),{y(0):Rational(1,2)})",
    "rsolve(y(n+1)-y(n),y(n),{y(0):sqrt(2)})", "rsolve(0*y(n)+1,y(n))", "rsolve(y(n+1)-y(n)-1/n**2,y(n))",
    "rsolve_hyper([-1,1],0,n,symbols=True)", "rsolve_poly([-1,1],n,n,symbols=True)",
    "rsolve_ratio([-1,n+1],0,n,symbols=True)", "rsolve_poly([1],n,n)", "rsolve_hyper([0,1],0,n)",
])
def test_assumptions_inits_and_options_match_sympy(source):
    def run(m):
        n = m.Symbol('n')
        names = dict(namespace(m, False), n=n, p=m.Symbol('p', positive=True, integer=True),
                     q=m.Symbol('q', real=True))
        try:
            return str(eval(source, names))
        except Exception as error:
            return 'EXC %s: %s' % (type(error).__name__, error)
    assert run(sm) == run(sp)


def test_symbols_option_returns_the_free_constants():
    n = sm.Symbol('n')
    solution, constants = sm.rsolve_hyper([-1, -1, 1], 0, n, symbols=True)
    assert [str(c) for c in constants] == ['C0', 'C1']
    assert solution.free_symbols == {n, *constants}
    assert sm.rsolve_poly([-1, 1], sm.sin(n), n) is None


def test_solutions_print_and_roundtrip_like_sympy():
    n, N = sm.Symbol('n'), sp.Symbol('n')
    y, Y = sm.Function('y'), sp.Function('y')
    actual = sm.rsolve(y(n + 1) - 2*y(n) - n, y(n), {y(0): 1})
    expected = sp.rsolve(Y(N + 1) - 2*Y(N) - N, Y(N), {Y(0): 1})
    assert str(actual) == str(expected) and actual.srepr() == sp.srepr(expected)
    assert sm.to_sympy(actual) == expected and sm.from_sympy(expected) == actual
    assert pickle.loads(pickle.dumps(actual)) == actual
    fibonacci = sm.rsolve(y(n + 2) - y(n + 1) - y(n), y(n), {y(0): 0, y(1): 1})
    assert [fibonacci.subs(n, i).simplify() for i in range(8)] == [0, 1, 1, 2, 3, 5, 8, 13]


def test_documented_deviations():
    n = sm.Symbol('n')
    y = sm.Function('y')
    # SymPy returns the bare ansatz C0 + C1*n for this inconsistent system
    assert sm.rsolve_poly([-1, n + 1], n**2, n) is None
    # SymPy crashes with AttributeError while gathering the terms
    assert str(sm.rsolve(y(n + 1) - y(n) - (n + 1)*2**n, y(n))) == '2**n*(n - 1) + C0'
    # SymPy treats the nonlinear product as two linear terms
    with pytest.raises(ValueError, match="expected, got 'y"):
        sm.rsolve(y(n + 1) - y(n)*y(n - 1), y(n))


@pytest.mark.parametrize('source,error', [
    ("rsolve(y(n)-y(n+k),y(n))", TypeError), ("rsolve(y(n)-y(k+1),y(n))", ValueError),
    ("rsolve(y(n)-y(n+1),x(n))", ValueError),
    ("rsolve(y(n)-sqrt(n)*y(n+1),y(n))", ValueError), ("rsolve(y(n)-y(n+1),y(n),{x(0):0})", ValueError),
    ("rsolve(y(n)+y(n+1)+2**n+cos(n),y(n))", ValueError), ("rsolve(y(n+1)-y(n)-1.5,y(n))", ValueError),
])
def test_invalid_input_raises_like_sympy(source, error):
    for m in (sm, sp):
        names = dict(namespace(m, True), x=m.Function('x'))
        with pytest.raises(error):
            eval(source, names)


HYPER_TERMS = ['factorial(k)', 'binomial(2*k, k)/4**k', '2**k*k', 'factorial(k)/k**k', 'sin(k)',
               'rf(k, 3)*3**k', '5', '1/(k**2 + 1)', 'gamma(k + 1)/gamma(k + Rational(1, 2))']


@pytest.mark.parametrize('source', HYPER_TERMS)
def test_hypersimp_matches_sympy(source):
    def run(m):
        k = m.Symbol('k')
        names = {n: getattr(m, n) for n in 'factorial binomial sin rf gamma Rational'.split()}
        names['k'] = k
        return m.hypersimp(eval(source, names), k)
    ours, theirs = run(sm), run(sp)
    if theirs is None:
        assert ours is None
    else:
        # same ratio; SymPy's simplify(ratio=oo) shape may differ
        assert sp.cancel(sp.sympify(str(ours)) - theirs) == 0


@pytest.mark.parametrize('pair', [('factorial(k)', 'factorial(k + 2)'), ('factorial(k)', '2**k*k'),
                                  ('factorial(k)', 'factorial(2*k)'), ('2**k', '3**k*k'),
                                  ('binomial(k, 2)', 'k**2')])
def test_hypersimilar_matches_sympy(pair):
    def run(m):
        k = m.Symbol('k')
        names = dict(factorial=m.factorial, binomial=m.binomial, k=k)
        return bool(m.hypersimilar(eval(pair[0], names), eval(pair[1], names), k))
    assert run(sm) == run(sp)
