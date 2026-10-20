"""Deterministic public-API workloads, shared by both benchmark backends.

Inputs are constructed outside the timed region. Sizes are fixed, not selected
after seeing which library wins. These examples are not a global parity proof.
"""
from dataclasses import asdict, dataclass
import importlib
import random


@dataclass(frozen=True)
class Case:
    id: str
    family: str
    operation: str
    size: int = 0
    validation: str = "exact"
    digits: int = 0

    def record(self):
        return asdict(self)


def catalog():
    result = []

    def add(family, operation, sizes=(0,), validation="exact", digits=0):
        for size in sizes:
            result.append(Case(f"{family}.{operation}" + (f".{size}" if size else ""),
                               family, operation, size, validation, digits))

    add("construction", "incremental_sum", (100, 500))
    add("construction", "bulk_add", (100, 1000))
    add("parsing", "sympify", (50, 300))
    add("expansion", "multinomial", (8, 12, 16))
    add("expansion", "univariate", (30, 100))
    add("differentiation", "composite", (1, 4, 8))
    add("substitution", "subs", (8, 12))
    add("substitution", "xreplace", (100, 500))
    add("simplification", "rational", (5, 12))
    add("simplification", "trig_identity", (4, 8))
    add("simplification", "tangent")
    add("simplification", "trigsimp", (2, 4, 8))
    add("simplification", "powsimp", (20, 100))
    add("simplification", "radsimp", (3, 5))
    add("polynomials", "factor_cyclotomic", (20, 50, 100))
    add("polynomials", "factor_wilkinson", (10, 20))
    add("polynomials", "factor_multivariate", (2, 4))
    add("polynomials", "factor_modular", (20, 50))
    add("polynomials", "gcd", (10, 30))
    add("polynomials", "cancel_sum", (5, 15))
    add("polynomials", "apart", (4, 8))
    add("polynomials", "resultant", (3, 6))
    add("polynomials", "roots")
    add("polynomials", "nroots", (10, 30), "numeric_roots", 30)
    add("polynomials", "groebner", (3, 4))
    add("polynomials", "radical_division")
    add("matrices", "det_integer", (8, 20))
    add("matrices", "det_generic", (3, 4))
    add("matrices", "inverse_rational", (4, 10, 20))
    add("matrices", "inverse_generic", (3,))
    add("matrices", "rref_integer", (10, 30))
    add("matrices", "rref_deficient", (12,))
    add("matrices", "nullspace", (12,), "nullspace")
    add("matrices", "linear_solve", (8, 20))
    add("matrices", "charpoly_integer", (6, 12))
    add("matrices", "charpoly_generic", (3, 4, 5, 6))
    add("matrices", "charpoly_shared", (4,))
    add("matrices", "eigenvalues")
    add("matrices", "multiply_symbolic", (6, 12))
    add("matrices", "power_symbolic", (4,))
    add("matrices", "jacobian", (8, 16))
    add("integration", "polynomial", (30, 100), "antiderivative")
    add("integration", "exponential_trig", validation="antiderivative")
    add("integration", "rational", (4, 8), "antiderivative")
    add("integration", "gamma_moment", (10, 50))
    add("integration", "dirac_affine")
    add("series", "transcendental", (10, 20, 40))
    add("series", "rational", (10, 30))
    add("limits", "removable")
    add("limits", "at_infinity")
    add("summation", "geometric", (100, 10000))
    add("solvers", "polynomial", (4, 7), "unordered")
    add("solvers", "linear_system", (3, 8), "unordered")
    add("solvers", "nonlinear_system", validation="unordered")
    add("solvers", "solveset_real")
    add("solvers", "nsolve", validation="numeric", digits=50)
    add("ode", "linear_first", validation="ode")
    add("ode", "harmonic", validation="ode")
    add("numerics", "constant", (50, 100, 300), "numeric")
    add("numerics", "tiny_sine", (50, 100), "numeric")
    add("numerics", "large_argument", (50, 100), "numeric")
    add("numerics", "complex_log", (50, 100), "numeric")
    add("assumptions", "tiny_positive")
    add("assumptions", "trig_zero")
    add("number_theory", "factorint")
    add("number_theory", "perfect_power")
    add("number_theory", "primepi", (100000, 1000000))
    add("combinatorics", "permutation", (100, 1000))
    add("combinatorics", "symmetric_group", (8, 12))
    add("combinatorics", "partition", (50, 200))
    add("statistics", "normal_moment")
    add("statistics", "binomial_variance")
    add("lambdify", "compile_and_evaluate", (20, 100), "array")
    add("lambdify", "evaluate", (20, 100), "array")
    add("printing", "str", (5, 8), "text")
    add("printing", "latex", (5, 8), "text")
    assert len({c.id for c in result}) == len(result)
    return result


def integer_rows(n, columns=None, seed=41):
    rng = random.Random(seed)
    return [[rng.randint(-8, 8) for _ in range(columns or n)] for _ in range(n)]


def build(case, m):
    """Return a nullary operation; constructing this closure is not timed."""
    op, n = case.operation, case.size
    x, y, z, t = m.symbols("x y z t")
    if op == "incremental_sum":
        def run():
            value = 0
            for i in range(n):
                value += i*x**i
            return value
        return run
    if op == "bulk_add":
        terms = [i*x**i for i in range(n)]
        return lambda: m.Add(*terms)
    if op == "sympify":
        source = " + ".join(f"{i}*x**{i}*sin(y)/(z+{i})" for i in range(1, n+1))
        return lambda: m.sympify(source)
    if op in ("multinomial", "univariate"):
        expression = (x+y+z+1)**n if op == "multinomial" else (x+1)**n
        return lambda: m.expand(expression)
    if op == "composite":
        expression = m.sin(x*y)**4*m.exp(x*x)*m.log(x+y)/(x*x+y*y+1)
        return lambda: m.diff(expression, x, n)
    if op == "subs":
        expression = m.expand((x+y+z+1)**n)
        return lambda: expression.subs({x:y+1,z:2})
    if op == "xreplace":
        expression = m.Add(*[m.sin(x+i)*m.cos(y+i) for i in range(n)])
        return lambda: expression.xreplace({x:y+1,y:z})
    if op == "rational" and case.family == "simplification":
        expression = (x**n-y**n)/(x-y)
        return lambda: m.simplify(expression)
    if op == "trig_identity":
        expression = m.sin(x)**n-m.cos(x)**n+m.cos(2*x)
        return lambda: m.simplify(expression)
    if op == "tangent":
        expression = m.tan(x)*m.cos(x)
        return lambda: m.simplify(expression)
    if op == "trigsimp":
        expression = (m.sin(x)**2+m.cos(x)**2)**n + m.tan(x)*m.cos(x)
        return lambda: m.trigsimp(expression)
    if op == "powsimp":
        symbols = m.symbols(f"a0:{n}", positive=True)
        expression = m.Mul(*[m.sqrt(a)*m.sqrt(a*x) for a in symbols])
        return lambda: m.powsimp(expression)
    if op == "radsimp":
        expression = 1/sum(m.sqrt(p) for p in (2,3,5,7,11)[:n])
        return lambda: m.radsimp(expression)
    if op == "factor_cyclotomic":
        expression = x**n-1
        return lambda: m.factor(expression)
    if op == "factor_wilkinson":
        expression = m.expand(m.prod(x-i for i in range(1,n+1)))
        return lambda: m.factor(expression)
    if op == "factor_multivariate":
        expression = m.expand((x**n*y+3*z-1)*(x*y*z+y*y-2)*(x-y+2*z))
        return lambda: m.factor(expression)
    if op == "factor_modular":
        expression = x**n-1
        return lambda: m.factor_list(expression,modulus=5)
    if op == "gcd":
        common=(x**n+3*x**(n//2)-x+5)**2
        a=m.expand(common*(x**n-7*x+1));b=m.expand(common*(x**(n-1)+2*x-9))
        return lambda:m.gcd(a,b)
    if op == "cancel_sum":
        expression=sum(1/(x+i) for i in range(1,n+1))
        return lambda:m.cancel(expression)
    if op == "apart":
        expression=(x*x+1)/m.expand(m.prod((x+i)**(1+i%2) for i in range(1,n+1)))
        return lambda:m.apart(expression,x)
    if op == "resultant":
        a=x**n+y*x**(n//2)+y*y*x+1;b=x**(n-1)-y**3*x+2
        return lambda:m.resultant(a,b,x)
    if op == "roots":
        expression=m.expand((x-1)*(x+2)*(x-3)*(x*x-2)*(x*x+x+1))
        return lambda:m.roots(expression,x)
    if op == "nroots":
        expression=x**n-3*x**(n//3)+x*x-1
        return lambda:m.nroots(expression,n=case.digits,maxsteps=200)
    if op == "groebner":
        variables=m.symbols(f"g0:{n}")
        equations=[sum(m.prod(variables[(i+j)%n] for j in range(k)) for i in range(n)) for k in range(1,n)]
        equations.append(m.prod(variables)-1)
        return lambda:m.groebner(equations,*variables,order="grevlex")
    if op == "radical_division":
        a=m.Poly(x*x+2*x+3,x);b=m.Poly((m.sqrt(2)+m.sqrt(3)+1)*x+m.sqrt(2),x)
        return lambda:a.div(b)
    if case.family == "matrices":
        if op in ("det_integer","rref_integer","linear_solve","charpoly_integer"):
            matrix=m.Matrix(integer_rows(n))
            if op=="det_integer":return lambda:matrix.det()
            if op=="rref_integer":return lambda:matrix.rref()
            if op=="linear_solve":
                rhs=m.Matrix(list(range(n)))
                return lambda:matrix.LUsolve(rhs)
            return lambda:matrix.charpoly(t)
        if op in ("det_generic","inverse_generic","charpoly_generic"):
            matrix=m.Matrix(n,n,m.symbols(f"a0:{n*n}"))
            if op=="det_generic":return lambda:matrix.det()
            if op=="inverse_generic":return lambda:matrix.inv()
            return lambda:matrix.charpoly(t)
        if op=="inverse_rational":
            rng=random.Random(42)
            matrix=m.Matrix([[m.Rational(rng.randint(-8,8),rng.randint(1,9)) for _ in range(n)] for _ in range(n)])
            return lambda:matrix.inv()
        if op=="rref_deficient":
            a=integer_rows(n+6,n,43);b=integer_rows(n,n+12,44)
            matrix=m.Matrix([[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n+12)] for i in range(n+6)])
            return lambda:matrix.rref()
        if op=="nullspace":
            matrix=m.Matrix(integer_rows(n,n+6,45))
            return lambda:matrix.nullspace()
        if op=="charpoly_shared":
            matrix=m.Matrix(n,n,lambda i,j:x**((i+j)%3)+y*(i-j)+(i+1)*(j+2)%3)
            return lambda:matrix.charpoly(t)
        if op=="eigenvalues":
            matrix=m.Matrix([[6,-2,2],[-2,3,-1],[2,-1,3]])
            return lambda:matrix.eigenvals()
        if op=="multiply_symbolic":
            matrix=m.Matrix(n,n,lambda i,j:x**(i%3)+y*j-i)
            return lambda:matrix*matrix
        if op=="power_symbolic":
            matrix=m.Matrix(n,n,lambda i,j:x+i-j if (i+j)%2 else 1)
            return lambda:matrix**6
        if op=="jacobian":
            variables=m.symbols(f"q0:{n}")
            vector=m.Matrix([m.sin(variables[i])*variables[(i+1)%n]**2+m.exp(variables[i]*variables[(i+3)%n]) for i in range(n)])
            return lambda:vector.jacobian(variables)
    if case.family=="integration":
        if op=="polynomial":
            expression=m.expand((x+1)**n)
            return lambda:m.integrate(expression,x)
        if op=="exponential_trig":
            expression=m.exp(2*x)*m.cos(3*x)
            return lambda:m.integrate(expression,x)
        if op=="rational":
            expression=(x*x+1)/m.expand(m.prod(x+i for i in range(1,n+1)))
            return lambda:m.integrate(expression,x)
        if op=="gamma_moment":return lambda:m.integrate(x**n*m.exp(-x),(x,0,m.oo))
        if op=="dirac_affine":
            expression=(x*x+1)*m.DiracDelta(2*x-2)
            return lambda:m.integrate(expression,(x,-2,3))
    if case.family=="series":
        expression=m.exp(m.sin(x))*m.cos(x) if op=="transcendental" else (1+x)/(x*(1-x-x*x))
        return lambda:expression.series(x,0,n)
    if case.family=="limits":
        if op=="removable":return lambda:m.limit((m.exp(x)-1-x)/x**2,x,0)
        return lambda:m.limit(m.sqrt(x*x+x)-x,x,m.oo)
    if case.family=="summation":
        k=m.Symbol("k",integer=True)
        return lambda:m.summation(2**k,(k,0,n))
    if case.family=="solvers":
        if op=="polynomial":
            expression=m.prod(x-i for i in range(1,n+1))
            return lambda:m.solve(expression,x)
        if op=="linear_system":
            variables=m.symbols(f"v0:{n}");rows=integer_rows(n,seed=47)
            equations=[sum(rows[i][j]*variables[j] for j in range(n))-i for i in range(n)]
            return lambda:m.solve(equations,variables,dict=True)
        if op=="nonlinear_system":return lambda:m.solve([x*x+y*y-5,x-y-1],(x,y),dict=True)
        if op=="solveset_real":return lambda:m.solveset(m.sin(x)-m.Rational(1,2),x,domain=m.S.Reals)
        return lambda:m.nsolve(m.cos(x)-x,x,1,prec=case.digits)
    if case.family=="ode":
        f=m.Function("f");fx=f(x)
        equation=m.diff(fx,x)-fx if op=="linear_first" else m.diff(fx,x,2)+fx
        return lambda:m.dsolve(equation,fx)
    if case.family=="numerics":
        if op=="constant":expression=m.exp(m.pi)*m.sqrt(2)+m.log(3)
        elif op=="tiny_sine":expression=m.sin(m.Rational(1,10**1000))
        elif op=="large_argument":expression=m.sin(m.Integer(10)**100)
        else:expression=m.log(-2+m.I)
        return lambda:expression.evalf(n)
    if case.family=="assumptions":
        expression=m.sin(m.Rational(1,10**1000)) if op=="tiny_positive" else m.sin(3)**2+m.cos(3)**2-1
        return (lambda:expression.is_positive) if op=="tiny_positive" else (lambda:expression.is_zero)
    if case.family=="number_theory":
        if op=="factorint":return lambda:m.factorint(2**64+1)
        if op=="perfect_power":return lambda:m.perfect_power((2**127-1)**12)
        return lambda:m.primepi(n)
    if case.family=="combinatorics":
        package=importlib.import_module(m.__name__+".combinatorics")
        if op=="permutation":
            a=package.Permutation(list(range(1,n))+[0]);b=package.Permutation(list(reversed(range(n))))
            return lambda:(a*b)**17
        if op=="symmetric_group":
            group=package.SymmetricGroup(n)
            return lambda:group.order()
        return lambda:m.partition(n)
    if case.family=="statistics":
        stats=importlib.import_module(m.__name__+".stats")
        if op=="normal_moment":
            mu=m.Symbol("mu",real=True);sigma=m.Symbol("sigma",positive=True);X=stats.Normal("X",mu,sigma)
            return lambda:stats.E(X**4)
        p=m.Symbol("p",positive=True);B=stats.Binomial("B",12,p)
        return lambda:stats.variance(B)
    if case.family=="lambdify":
        import numpy as np
        expression=m.Add(*[m.sin(x/(i+1))*m.cos(y/(i+1)) for i in range(n)])
        values=np.linspace(0,1,8192);other=np.linspace(1,2,8192)
        if op=="evaluate":
            function=m.lambdify((x,y),expression,"numpy")
            return lambda:function(values,other)
        return lambda:m.lambdify((x,y),expression,"numpy")(values,other)
    if case.family=="printing":
        expression=m.expand((x+y+z+1)**n)
        return (lambda:str(expression)) if op=="str" else (lambda:m.latex(expression))
    raise ValueError(f"Unknown workload: {case.id}")
