"""Differential corpus of ``has`` / ``has_free`` / ``has_xfree``.

usage:
  has_corpus.py harvest OUT.json.gz          # build the cases (needs SymPy)
  has_corpus.py sympy CORPUS.json.gz OUT.json  # SymPy 1.14 answers
  has_corpus.py symox CORPUS.json.gz OUT.json  # symox answers (no SymPy import)
  has_corpus.py report SYMPY.json SYMOX.json [BEFORE.json]

A case is ``{"id", "src", "code"}``: ``code`` is evaluated in the namespace
built by ``namespace(module)``.  An answer is ``repr`` of the value, or
``EXC <class name>``.  Verdicts:

  SAME      identical answer;
  WRONG     both answered booleans that differ (a false result);
  DIFF      an exception on one side only, or different exception classes;
  DECLINED  symox raised NotImplementedError where SymPy answered.
"""
import gzip
import json
import random
import sys


def namespace(module):
    ns = dict(vars(module))
    x, y, z, w = module.symbols('x y z w')
    A, B, C = module.symbols('A B C', commutative=False)
    p = module.Symbol('p', positive=True)
    f, g = module.Function('f'), module.Function('g')
    a, b = module.Wild('a'), module.Wild('b')
    n = module.Symbol('n', integer=True)
    M = module.MatrixSymbol('M', 2, 2)
    N = module.MatrixSymbol('N', 2, 2)
    ns.update(x=x, y=y, z=z, w=w, A=A, B=B, C=C, p=p, f=f, g=g, a=a, b=b, n=n, M=M, N=N)
    return ns


HAND = r'''
# sub-products
(2*x*y*z).has(x*y)
(2*x*y*z).has(2*x)
(2*x*y*z).has(2*y*z)
(2*x*y*z).has(4*x)
(2*x*y).has(-2*x)
(-x*y).has(x*y)
(-x*y*z).has(-x*y)
(-2*x*y).has(-x)
(-2*x*y).has(2*x)
(-2*x*y).has(-2*x)
(x*y*z).has(x*y)
(x*y*z).has(-x*y)
(x*y*z).has(y*z*x)
(x**2*y).has(x*y)
(x**2*y*z).has(x**2*z)
(x**2*y*z).has(x*z)
(x*y*sin(x*y*z)).has(x*z)
exp(x*y*z).has(x*z)
(x**2+x*y*z).has(x*z)
(sqrt(2)*x*y).has(sqrt(2)*x)
(sqrt(2)*x*y).has(2*x)
(I*x*y).has(I*x)
(I*x*y).has(-I*x)
(pi*x*y).has(pi*y)
(E*x*y).has(E*x)
(Rational(1,2)*x*y).has(x/2)
(x*y/z).has(x/z)
(x*y/z).has(y/z)
(x*y/z).has(1/z)
(x*y/z).has(x*y)
(x/(y*z)).has(1/y)
(x/(y*z)).has(1/(y*z))
(x*y*exp(z)).has(x*exp(z))
(x*y*f(z)).has(y*f(z))
(x*y*f(z)).has(x*f)
(x*sin(y)*cos(z)).has(sin(y)*cos(z))
(2.0*x*y).has(2.0*x)
(2.0*x*y).has(2*x)
(oo*x*y).has(oo*x)
(x*y*z*w).has(x*w)
(x*y*z*w).has(x*y*z*w*2)
(3*x*y).has(3*x*y)
(3*x*y).has(3)
(3*x*y).has(1)
(-x).has(-1)
(x/2).has(2)
(x/2).has(S.Half)
(1/x).has(1)
(1/x).has(-1)
(x**-2).has(-2)
(x - 1).has(-1)
(x*y).has(x*y*z)
# sub-sums
(x+y+1).has(x+1)
(x+y+2).has(x+1)
(x+y+z).has(x+z)
(x+y+z).has(y+z)
(x+y+z).has(x+y+z+w)
(x-y+z).has(x-y)
(x-y+z).has(x+y)
(x-y+z).has(-y+z)
(x+2*y+z).has(2*y+z)
(x+2*y+z).has(y+z)
(2*x+2*y).has(x+y)
(x+y+1).has(x+y)
(x+y+1).has(1)
(x+y).has(x+1)
(sin(x+y+z)).has(x+z)
((x+y+z)**2).has(x+y)
(x*(y+z+w)).has(y+w)
(x + y*z + w).has(x + w)
(x + y*z + w).has(y*z + w)
(x + y*z + w).has(z + w)
(x + y + I).has(x + I)
(x + y + pi).has(y + pi)
(x + 1.5 + y).has(x + 1.5)
(x + y + oo).has(x + oo)
# non commutative
(A*B*C).has(A*C)
(A*B*C).has(A*B)
(A*B*C).has(B*C)
(A*B*C).has(B*A)
(x*A*B*C).has(x*B*C)
(x*A*B*C).has(A*B)
(x*y*A*B).has(x*A)
(x*y*A*B).has(x*B)
(x*y*A*B).has(y*A*B)
(2*A*B).has(2*A)
(A*B*A).has(A*A)
(A*B*A).has(B*A)
(A**2*B).has(A*B)
(A+B+C).has(A+C)
(A+B+C).has(A+B)
(A+B+C).has(B+C)
(A+B+x).has(A+x)
(A+B+x).has(B+x)
(x*A + B + y*C).has(B + y*C)
(x*A + B + y*C).has(x*A + y*C)
(A*B + C).has(A*B)
(A*B + C).has(B*A)
(A*B - B*A).has(A*B)
(f(A)*B).has(f(A))
# powers
(x**2).has(x)
(x**2).has(2)
(x**y).has(y)
(x**(y+z+1)).has(y+1)
(x**(2*y*z)).has(2*y)
(exp(x*y*z)).has(exp(x))
(exp(x+y)).has(exp(x))
(exp(x)*exp(y)).has(exp(x))
(sqrt(x*y)).has(x*y)
(sqrt(x*y)).has(S.Half)
((x*y)**z).has(x*y)
(2**x).has(2)
(x**2*y**2).has(x**2)
(x**2*y**2).has((x*y)**2)
# functions and classes
f(x).has(f)
f(x).has(g)
f(x).has(Function('f'))
f(x, y).has(f)
g(f(x)).has(f)
sin(x).has(sin)
sin(x).has(cos)
sin(x).has(Function)
sin(x).has(Basic)
sin(x).has(Expr)
x.has(Symbol)
x.has(Dummy)
x.has(Wild)
(x+1).has(Integer)
(x+1).has(Rational)
(x+1).has(Number)
(x+S.Half).has(Rational)
(x+S.Half).has(Integer)
(x+1.5).has(Float)
(x*y).has(Mul)
(x*y).has(Add)
(x*y+1).has(Mul)
(x**2).has(Pow)
(1/x).has(Pow)
sqrt(x).has(Pow)
exp(x).has(Pow)
exp(x).has(exp)
(x+pi).has(pi)
(x+pi).has(NumberSymbol)
(x+E).has(exp)
(x+E).has(E)
f(x).has(AppliedUndef)
f(x).has(UndefinedFunction)
sin(x).has(AppliedUndef)
sin(f(x)).has(AppliedUndef)
Derivative(f(x),x).has(Derivative)
Integral(x,x).has(Integral)
Integral(x,x).has(Derivative)
# Wild
x.has(a)
(x+a).has(a)
(x+a).has(b)
(x*a).has(x*a)
(x*y*a).has(x*a)
a.has(a)
a.has(Wild)
(a+b).has(Wild)
sin(a).has(sin)
# multiple arguments
(x+y).has(x, z)
(x+y).has(z, w)
(x+y).has()
(x+y).has(z, x+y)
(x*y*z).has(w, x*z)
(x+y+z).has(w, x+z)
(x+y).has(sin, x)
sin(x).has(cos, sin)
(x+y).has(Mul, z)
(x*y).has(Mul, z)
x.has(None)
x.has(None, x)
x.has('x')
x.has(1)
(x+1).has(1)
(x+1).has(1.0)
x.has([x])
x.has((x,))
(x+1).has(True)
# derivatives
Derivative(f(x),x).has(x)
Derivative(f(x),x).has(1)
Derivative(f(x),x).has(f)
Derivative(f(x),(x,2)).has(2)
Derivative(f(x),(x,2)).has(1)
Derivative(f(x,y),x,y).has(y)
Derivative(f(x,y),x,y).has(1)
Derivative(x*y*f(x),x).has(x*f(x))
Derivative(f(x),x).has(Tuple(x,1))
Derivative(f(x),x).has(Derivative(f(x),x))
(2*Derivative(f(x),x)*x).has(2*Derivative(f(x),x))
# integrals
Integral(x,(x,0,1)).has(0)
Integral(x,(x,0,1)).has(1)
Integral(x,(x,0,1)).has(x)
Integral(x*y*z,(x,0,1)).has(x*z)
Integral(f(x),x).has(f)
Integral(x,(x,0,y)).has(y)
Integral(x,(x,0,y)).has(Tuple(x,0,y))
Integral(x+y+z,x).has(x+z)
Sum(x*y*n,(n,0,10)).has(x*n)
Sum(x*y*n,(n,0,10)).has(10)
Sum(n,(n,0,y)).has(y)
Product(n*x*y,(n,1,y)).has(n*y)
Limit(x*y*z,x,0).has(y*z)
Limit(x*y*z,x,0).has(x*z)
Subs(f(x)*y*z, x, 0).has(y*z)
Subs(f(x)*y*z, x, 0).has(f(x)*z)
Lambda(x, x*y*z).has(y*z)
Lambda(x, x*y*z).has(x*z)
# relationals, piecewise, logic
Eq(x*y*z,1).has(x*z)
Eq(x+y+z,1).has(x+z)
(x*y*z > 1).has(x*z)
Ne(x,y).has(x)
Piecewise((x*y*z,x>0),(1,True)).has(x*z)
Piecewise((x*y*z,x>0),(1,True)).has(x>0)
Piecewise((x*y*z,x>0),(1,True)).has(True)
Piecewise((x,x>0),(y,True)).has(Piecewise)
Piecewise((x+y+z,x>0),(1,True)).has(x+z)
And(x>0,y>0,z>0).has(And(x>0,z>0))
And(x>0,y>0,z>0).has(x>0)
Or(x>0,y>0,z>0).has(Or(x>0,z>0))
Or(x>0,y>0,z>0).has(And(x>0,z>0))
Max(x,y,z).has(Max(x,z))
Min(x,y,z).has(Min(x,z))
Max(x,y,z).has(Min(x,z))
Max(x,y,z,1).has(Max(x,1))
Max(x*y*z,w).has(x*z)
Union(Interval(0,1),Interval(2,3),Interval(4,5)).has(Union(Interval(0,1),Interval(4,5)))
Union(Interval(0,1),Interval(2,3),Interval(4,5)).has(Interval(2,3))
Union(Interval(0,1),Interval(2,3),Interval(4,5)).has(3)
Interval(x*y*z, 1).has(x*z)
FiniteSet(x*y*z, 1).has(x*z)
FiniteSet(x,y,z).has(FiniteSet(x,z))
Tuple(x*y*z, 1).has(x*z)
Tuple(x,y,z).has(Tuple(x,z))
# matrices
Matrix([[x*y*z]]).has(x*z)
Matrix([[x*y*z, 1],[x+y+z, 0]]).has(x+z)
Matrix([[x, 1],[y, 0]]).has(z)
Matrix([[x, 1],[y, 0]]).has(0)
Matrix([[x, 1],[y, 0]]).has(Symbol)
ImmutableMatrix([[x*y*z]]).has(x*z)
ImmutableMatrix([[x, 1],[y, 0]]).has(y)
M.has(M)
(M*N).has(M)
(M+N).has(N)
(2*M*N).has(2*M)
(M*N).has(N*M)
(M+N+Identity(2)).has(M+Identity(2))
(M*N*M).has(M*M)
(M*N*M).has(N*M)
M[0,1].has(M)
M[0,1].has(1)
(x*M[0,1]*y).has(x*M[0,1])
Trace(M*N).has(M)
Determinant(M).has(M)
Transpose(M).has(M)
# mixed
(x*y*z + sin(x*y*z)).has(x*z)
log(x*y*z).has(x*z)
(x*y*z).has(x*z, y)
(x*y*z).has(x*z*w)
(x*y*z*w).has(x*z*w)
cos(2*x*y).has(2*x)
cos(2*x*y).has(x*y)
(x*y + z).has(x*y + z)
(x*y + z).has(x + z)
((x+y)*(x+z)).has((x+y)*(x+z))
((x+y)*(x+z)*w).has((x+y)*w)
((x+y)*(x+z)*w).has(x+y)
(x*(x+y+z)).has(x+z)
(y*exp(x*y*z)+1).has(exp(x*y*z)+1)
Abs(x*y*z).has(x*z)
conjugate(x*y*z).has(x*z)
re(x*y*z).has(x*z)
O(x**2*y).has(x**2)
O(x*y*z, x).has(y*z)
(x*y*z + O(x**4)).has(x*z)
binomial(x*y*z, 2).has(x*z)
factorial(x*y*z).has(x*z)
gamma(2*x*y).has(2*x)
hyper([x*y*z],[1],x).has(x*z)
meijerg([[x*y*z],[]],[[],[]],x).has(x*z)
KroneckerDelta(x*y*z, 0).has(x*z)
Heaviside(x*y*z).has(x*z)
DiracDelta(x*y*z).has(x*z)
UnevaluatedExpr(x*y*z).has(x*z)
(UnevaluatedExpr(x)*y*z).has(y*z)
Poly(x*y + 1, x).has(x)
Poly(x*y + 1, x).has(y)
CRootOf(x**5 - x + 1, 0).has(x)
CRootOf(x**5 - x + 1, 0).has(5)
RootSum(x**3 + x + 1, Lambda(x, log(x))).has(x)
# has_free
(x*y*z).has_free(x*y)
(x*y*z).has_free(x*z)
(2*x*y*z).has_free(x*y)
(x + y + z).has_free(x+z)
Integral(x*y, (x, 0, 1)).has_free(x)
Integral(x*y, (x, 0, 1)).has_free(y)
Integral(x*y, (x, 0, 1)).has_free(x*y)
Integral(x*y*z, (x, 0, 1)).has_free(y*z)
Integral(x, (x, 0, y)).has_free(y)
Sum(n*x, (n, 0, 3)).has_free(n)
Sum(n*x, (n, 0, 3)).has_free(x)
Lambda(x, x*y).has_free(x)
Lambda(x, x*y).has_free(y)
Subs(f(x)*y, x, 0).has_free(x)
Subs(f(x)*y, x, 0).has_free(y)
Derivative(f(x),x).has_free(x)
(x*y).has_free(Symbol)
f(x).has_free(f)
f(x).has_free(x)
(A*B*C).has_free(A*C)
(A*B*C).has_free(B*C)
(x+y).has_free()
(x+y).has_free([x])
(x+y).has_free(x, z)
(x+y).has_free(z)
Integral(x*y,(x,0,1)).has_free(Symbol)
# has_xfree
(x*y).has_xfree({x})
(x*y).has_xfree({z})
(x*y).has_xfree({x*y})
(x*y*z).has_xfree({x*y})
Integral(x*y,(x,0,1)).has_xfree({x})
Integral(x*y,(x,0,1)).has_xfree({y})
(x*y).has_xfree([x])
(x*y).has_xfree(x)
# lattice and boolean operators
(x & y & z).has(x & z)
(x | y | z).has(x | z)
(x | y | z).has(x & z)
(x & y & (x > 1)).has(x & (x > 1))
(x & y & (x > 1)).has(y & (x > 1))
(x & (x > 1) & (y > 1)).has((x > 1) & (y > 1))
Or(x, y, Eq(x, 1)).has(Or(x, Eq(x, 1)))
And(x, y, Not(z)).has(And(x, Not(z)))
And(x, y, Not(z)).has(And(y, Not(z)))
Union(Interval(0,1),Interval(2,3),Interval(4,5)).has(Union(Interval(0,1),Interval(2,3)))
Union(Interval(0,1),Interval(2,3),Interval(4,5)).has(Union(Interval(2,3),Interval(4,5)))
Union(Interval(0,1),FiniteSet(7),Interval(4,5)).has(Union(Interval(0,1),FiniteSet(7)))
Union(Interval(0,1),FiniteSet(7),Interval(4,5)).has(Union(Interval(0,1),Interval(4,5)))
Max(x, y, sin(z)).has(Max(x, sin(z)))
Min(x, y, sin(z), 2).has(Min(2, sin(z)))
Max(x, Min(y, z, w)).has(Min(y, w))
(x & y & z).has_free(x & z)
Max(x, y, z).has_free(Max(x, z))
x.has_free({x, y})
x.has_free(y, [])
Piecewise((FiniteSet(1, 2), x > 3), (4, True)).has_free(FiniteSet(1, 2))
Piecewise((1, x > 3), (4, True)).has_free(FiniteSet(1, 2))
x.has_free(Interval(0, 1))
x.has_free(ImmutableMatrix([[x]]))
x.has_free(Matrix([[x]]))
'''


ATOMS = ['x', 'y', 'z', 'w', '2', '3', '-1', 'S.Half', 'pi', 'I', 'sin(x)', 'exp(y)', 'f(x)', 'sqrt(2)']
NC_ATOMS = ['A', 'B', 'C']


def _random_expr(rng, depth, nc=False):
    atoms = ATOMS + (NC_ATOMS if nc else [])
    if depth == 0 or rng.random() < 0.25:
        return rng.choice(atoms)
    kind = rng.choice(['add', 'mul', 'mul', 'add', 'pow', 'fun'])
    if kind in ('add', 'mul'):
        k = rng.randint(2, 4)
        parts = [_random_expr(rng, depth - 1, nc) for _ in range(k)]
        return '(' + (' + ' if kind == 'add' else '*').join(parts) + ')'
    if kind == 'pow':
        return '(%s)**(%s)' % (_random_expr(rng, depth - 1, False), rng.choice(['2', '3', '-1', 'y', 'S.Half']))
    return '%s(%s)' % (rng.choice(['sin', 'cos', 'exp', 'log', 'f', 'g']), _random_expr(rng, depth - 1, nc))


def harvest(out):
    import sympy
    ns = namespace(sympy)
    rng = random.Random(20261106)
    cases = []
    src = 'hand'
    for line in HAND.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('#'):
            src = line[1:].strip()
            continue
        cases.append({'src': src, 'code': line})
    # Random expressions; patterns are random sub-collections of the
    # arguments of their Add/Mul nodes (true positives under the reference
    # subset rule), perturbed variants (likely negatives) and other random
    # expressions.
    seen = set()
    while len([c for c in cases if c['src'].startswith('random')]) < 900:
        nc = rng.random() < 0.2
        code = _random_expr(rng, 3, nc)
        try:
            e = eval(code, ns)
        except Exception:
            continue
        nodes = [node for node in sympy.preorder_traversal(e) if isinstance(node, (sympy.Add, sympy.Mul))]
        patterns = []
        for node in nodes[:4]:
            args = list(node.args)
            if len(args) < 2:
                continue
            k = rng.randint(1, len(args))
            sub = rng.sample(args, k) if not nc else args[:k]
            patterns.append(str(node.func(*sub)))
            other = rng.choice(ATOMS)
            patterns.append(str(node.func(*(sub + [eval(other, ns)]))))
        patterns.append(_random_expr(rng, 2, nc))
        for pattern in patterns:
            key = (str(e), pattern)
            if key in seen:
                continue
            seen.add(key)
            cases.append({'src': 'random-nc' if nc else 'random', 'code': '(%s).has(%s)' % (e, pattern)})
            if rng.random() < 0.15:
                cases.append({'src': 'random-free', 'code': '(%s).has_free(%s)' % (e, pattern)})
    for i, case in enumerate(cases):
        case['id'] = i
    with gzip.open(out, 'wt') as fh:
        json.dump({'cases': cases}, fh)
    print(len(cases), 'cases')


def _answer(code, ns):
    try:
        value = eval(code, ns)
    except Exception as error:
        return 'EXC %s' % type(error).__name__
    return repr(value)


def replay(module_name, corpus, out):
    if module_name == 'symox':
        class NoReference:
            def find_spec(self, fullname, path=None, target=None):
                if fullname == 'sympy' or fullname.startswith('sympy.'):
                    raise ImportError('reference forbidden')
        sys.meta_path.insert(0, NoReference())
    module = __import__(module_name)
    ns = namespace(module)
    cases = json.loads(gzip.open(corpus).read())['cases']
    rows = {case['id']: _answer(case['code'], ns) for case in cases}
    json.dump({'cases': cases, 'answers': rows}, open(out, 'w'))


def verdict(ref, got):
    if ref == got:
        return 'SAME'
    if got == 'EXC NotImplementedError' and not ref.startswith('EXC'):
        return 'DECLINED'
    if ref in ('True', 'False') and got in ('True', 'False'):
        return 'WRONG'
    return 'DIFF'


def report(sympy_out, symox_out, before_out=None):
    ref = json.load(open(sympy_out))
    got = json.load(open(symox_out))
    cases = ref['cases']
    counts, by_src, bad = {}, {}, []
    for case in cases:
        key = str(case['id'])
        v = verdict(ref['answers'][key], got['answers'][key])
        counts[v] = counts.get(v, 0) + 1
        by_src.setdefault(case['src'], {}).setdefault(v, 0)
        by_src[case['src']][v] += 1
        if v != 'SAME':
            bad.append((v, case['code'], ref['answers'][key], got['answers'][key]))
    print(json.dumps({'total': len(cases), 'counts': counts}, sort_keys=True))
    for src, c in by_src.items():
        print('  %-20s %s' % (src, json.dumps(c, sort_keys=True)))
    for row in bad:
        print('%-8s %s  sympy=%s symox=%s' % row)
    if before_out:
        before = json.load(open(before_out))
        changed = [(case['code'], before['answers'][str(case['id'])], got['answers'][str(case['id'])])
                   for case in cases if before['answers'][str(case['id'])] != got['answers'][str(case['id'])]]
        print('changed vs before:', len(changed))
        for row in changed:
            print('  %s  before=%s after=%s' % row)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'harvest':
        harvest(sys.argv[2])
    elif mode in ('sympy', 'symox'):
        replay(mode, sys.argv[2], sys.argv[3])
    else:
        report(*sys.argv[2:])
