"""Factor-system decomposition, triangular solving and integer equations."""
from itertools import product
from math import gcd,isqrt


def factor_system(eqs,gens=(),**kwargs):
    from . import Poly,sympify,factor_list,groebner
    equations = [sympify(e.lhs-e.rhs if getattr(e,'is_Equality',False) else e) for e in eqs]
    gens = tuple(gens) or tuple(sorted(set().union(*(e.free_symbols for e in equations)),key=str))
    if not gens:return [[]] if all(e == 0 for e in equations) else []
    choices = []
    for equation in equations:
        if equation == 0:continue
        if not equation.free_symbols:return []
        coefficient,factors = Poly(equation,*gens,**kwargs).factor_list()
        choices.append(tuple(f.as_expr() for f,m in factors))
    if not choices:return [[]]
    systems = []
    for selected in product(*choices):
        selected = list(dict.fromkeys(selected))
        basis = list(groebner(selected,*gens))
        if any(not f.free_symbols and f != 0 for f in basis):continue
        canonical = frozenset(basis)
        if any(existing <= canonical for existing in systems):continue
        systems = [existing for existing in systems if not canonical < existing]
        systems.append(canonical)
    return [sorted(s,key=str) for s in sorted(systems,key=lambda s:(len(s),str(sorted(s,key=str))))]


def solve_triangulated(polys,*gens,**args):
    from . import solve_poly_system
    # The native solver performs lexicographic elimination and rejects partial
    # root lists. Extension-domain solving is explicitly outside this adapter.
    extension = args.pop('extension',False)
    if extension:raise NotImplementedError('triangulated algebraic extension towers are not yet supported')
    return solve_poly_system(polys,*gens,**args)


def substitution(system,symbols,result=None,known_symbols=(),exclude=(),all_symbols=None):
    from . import sympify,solve,FiniteSet,Tuple,EmptySet,simplify
    symbols = tuple(symbols)
    if not symbols or any(not s.is_Symbol for s in symbols):raise ValueError('supply a sequence of symbols')
    all_symbols = symbols if all_symbols is None else tuple(all_symbols)
    equations = [sympify(e.lhs-e.rhs if getattr(e,'is_Equality',False) else e) for e in system]
    if not equations:return EmptySet
    rows = []
    for known in ([{}] if result is None or result == [] else result):
        residual = [e.subs(known) for e in equations]
        solved = solve(residual,*symbols,dict=True)
        for branch in solved:
            values = dict(known);values.update(branch)
            # Fully resolve chained substitutions without an unbounded loop.
            for i in range(len(values)):
                updated = {k:v.subs(values,simultaneous=True) for k,v in values.items()}
                if updated == values:break
                values = updated
            if any(simplify(sympify(e).subs(values)) == 0 for e in exclude):continue
            if any(simplify(e.subs(values)).is_zero is False for e in equations):continue
            rows.append(Tuple(*(values.get(s,s) for s in all_symbols)))
    return FiniteSet(*rows) if rows else EmptySet


def _egcd(a,b):
    old_r,r,old_s,s,old_t,t = a,b,1,0,0,1
    while r:
        q = old_r//r
        old_r,r,old_s,s,old_t,t = r,old_r-q*r,s,old_s-q*s,t,old_t-q*t
    return (old_s,old_t,old_r) if old_r >= 0 else (-old_s,-old_t,-old_r)


def _linear(coefficients,target,parameters):
    # Unimodular column reduction: A*U = (gcd(A),0,...,0).
    n = len(coefficients)
    U = [[int(i == j) for j in range(n)] for i in range(n)]
    row = list(coefficients)
    for j in range(1,n):
        a,b = row[0],row[j]
        if a == b == 0:continue
        s,t,g = _egcd(a,b)
        first,other = [U[i][0] for i in range(n)],[U[i][j] for i in range(n)]
        for i in range(n):
            U[i][0] = s*first[i]+t*other[i]
            U[i][j] = -b//g*first[i]+a//g*other[i]
        row[0],row[j] = g,0
    g = row[0]
    if g == 0:return {tuple(parameters[:n])} if target == 0 else set()
    if target % g:return set()
    return {tuple(U[i][0]*(target//g)+sum(U[i][j]*parameters[j-1] for j in range(1,n)) for i in range(n))}


def _pell_unit(D):
    a0 = isqrt(D)
    if a0*a0 == D:raise ValueError('Pell discriminant must not be a square')
    m,d,a = 0,1,a0
    p0,p1,q0,q1 = 1,a,0,1
    while p1*p1-D*q1*q1 != 1:
        m = d*a-m;d = (D-m*m)//d;a = (a0+m)//d
        p0,p1,q0,q1 = p1,a*p1+p0,q1,a*q1+q0
    return p1,q1


def diophantine(eq,param=None,syms=None,permute=False):
    from . import sympify,Symbol,Poly,Add,expand,divisors,sqrt,Integer
    eq = sympify(eq.lhs-eq.rhs if getattr(eq,'is_Equality',False) else eq)
    variables = tuple(sorted(eq.free_symbols,key=str)) if syms is None else tuple(syms)
    if not variables:return {()} if eq == 0 else set()
    if any(not v.is_Symbol for v in variables):raise TypeError('Diophantine variables must be symbols')
    if eq.free_symbols-set(variables):raise ValueError('all variables must be included in syms')
    p = Poly(eq,*variables)
    if any(c.is_Rational is not True for c in p.coeffs()):raise TypeError('rational polynomial coefficients are required')
    _,p = p.clear_denoms(convert=True)
    name = 't' if param is None else str(param)
    parameters = [Symbol('%s_%s' % (name,i),integer=True) for i in range(max(2,len(variables)))]
    if int(p.total_degree()) <= 1:
        coeffs = [int(p.as_expr().diff(v)) for v in variables]
        target = -int(p.as_expr().subs(dict.fromkeys(variables,0)))
        return _linear(coeffs,target,parameters)
    if len(variables) == 1:
        return {(r,) for r in p.ground_roots() if r.is_Integer}
    factors = p.factor_list()[1]
    if len(factors) > 1:
        result = set()
        for factor,m in factors:
            present = tuple(v for v in variables if factor.as_expr().has(v))
            for row in diophantine(factor.as_expr(),param=param,syms=present,permute=permute):
                mapping = dict(zip(present,row))
                result.add(tuple(mapping.get(v,parameters[variables.index(v)]) for v in variables))
        return result
    if len(variables) == 2:
        x,y = variables
        expression = expand(p.as_expr())
        a,b,c = (int(p.coeff_monomial(term)) for term in (x*x,x*y,y*y))
        d,e,f = (int(p.coeff_monomial(term)) for term in (x,y,Integer(1)))
        if int(p.total_degree()) == 2 and b == 0 and a > 0 and c > 0:
            # Complete the squares to derive finite, proven bounds.
            bound = d*d*c+e*e*a-4*a*c*f
            if bound < 0:return set()
            radius = isqrt(bound//c)
            lo,hi = (-radius-d)//(2*a)-1,(radius-d)//(2*a)+1
            result = set()
            for xv in range(lo,hi+1):
                discriminant = e*e-4*c*(a*xv*xv+d*xv+f)
                if discriminant < 0:continue
                root = isqrt(discriminant)
                if root*root != discriminant:continue
                for numerator in {-e-root,-e+root}:
                    if numerator%(2*c) == 0:result.add((Integer(xv),Integer(numerator//(2*c))))
            return result
        if int(p.total_degree()) == 2 and a == c == d == e == 0 and b:
            if -f%b:return set()
            target = -f//b
            if target == 0:return {(Integer(0),parameters[0]),(parameters[0],Integer(0))}
            return {(Integer(sign*k),Integer(target//(sign*k))) for k in divisors(abs(target)) for sign in (-1,1)}
        if int(p.total_degree()) == 2 and a == 1 and b == d == e == 0 and c < 0 and f == -1:
            D = -c
            if isqrt(D)**2 == D:
                return {(Integer(1),Integer(0)),(Integer(-1),Integer(0))}
            u,v = _pell_unit(D);t = parameters[0]
            unit = u+v*sqrt(D)
            real = (unit**t+unit**(-t))/2
            imag = (unit**t-unit**(-t))/(2*sqrt(D))
            return {(sign*real,imag) for sign in (-1,1)}
    if len(variables) == 3 and int(p.total_degree()) == 2:
        x,y,z = variables
        if expand(p.as_expr()-(x*x+y*y-z*z)) == 0:
            m,n,k = parameters[0],parameters[1],Symbol(name+'_2',integer=True)
            first,second,third = k*(m*m-n*n),2*k*m*n,k*(m*m+n*n)
            return {(first,second,third),(second,first,third),(first,second,-third),(second,first,-third)}
    raise NotImplementedError('this integer polynomial family is not implemented; no incomplete solution set is returned')
