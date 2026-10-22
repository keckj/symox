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
    """SymPy's Gianni-Mora triangulated solver. Without ``extension`` only
    the roots in the ground domain are returned, as in the reference."""
    from . import groebner, default_sort_key
    from .polys.polyoptions import build_options
    from .utilities.iterables import postfixes
    opt = build_options(gens, args)
    if opt.get('extension', False):
        # The reference lifts each root into an algebraic field built from
        # CRootOf generators; that domain tower is not available natively.
        raise NotImplementedError('triangulated algebraic extension towers are not yet supported')
    G = list(groebner(polys, gens, polys=True))[::-1]
    domain = opt.get('domain')
    if domain is not None:
        G = [g.set_domain(domain) for g in G]

    def _solve_univariate(f):
        return list(f.ground_roots().keys())

    f, G = G[0].ltrim(-1), G[1:]
    solutions = {(zero,) for zero in _solve_univariate(f)}
    for var, vars in zip(reversed(gens[:-1]), postfixes(gens[1:])):
        _solutions = set()
        for values in solutions:
            H, mapping = [], dict(zip(vars, values))
            for g in G:
                if g.has_only_gens(var, *vars) and g.degree(var) != 0:
                    h = g.ltrim(var).eval(mapping)
                    if g.degree(var) == h.degree():
                        H.append(h)
            p = min(H, key=lambda h: h.degree())
            for zero in _solve_univariate(p):
                _solutions.add((zero,) + values)
        solutions = _solutions
    return sorted(solutions, key=default_sort_key)


def substitution(system,symbols,result=None,known_symbols=(),exclude=(),all_symbols=None):
    # The reference substitution method, which keeps periodic solutions as
    # ImageSet families.
    from .solvers._substitution import substitution as reference
    return reference(system,symbols,result,known_symbols,exclude,all_symbols)


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


def _igcdex(a,b):
    # reference `igcdex`
    if not a and not b:return (0,1,0)
    if not a:return (0,b//abs(b),abs(b))
    if not b:return (a//abs(a),0,abs(a))
    xs,ys = (-1 if a < 0 else 1),(-1 if b < 0 else 1)
    a,b = abs(a),abs(b)
    x,y,r,s = 1,0,0,1
    while b:
        c,q = a % b,a//b
        a,b,r,s,x,y = b,c,x-q*r,y-q*s,r,s
    return (x*xs,y*ys,a)


def _base_solution_linear(c,a,b,t):
    # reference `base_solution_linear`: (numeric part, parameter part) pairs
    g = gcd(gcd(a,b),c)
    if g > 1:a,b,c = a//g,b//g,c//g
    if b < 0:t = -t
    if c == 0:return (0,b*t),(0,-a*t)
    x0,y0,d = _igcdex(abs(a),abs(b))
    x0 *= (a > 0)-(a < 0)
    y0 *= (b > 0)-(b < 0)
    if c % d:return None
    return (c*x0,b*t),(c*y0,-a*t)


def _linear_reference(A,c,params):
    # reference `Linear.solve` for two or more variables; `None` when there
    # is no solution
    from . import Add,Integer
    A = list(A)
    B = []
    if len(A) > 2:
        B.append(gcd(A[-2],A[-1]))
        A[-2] = A[-2]//B[0]
        A[-1] = A[-1]//B[0]
        for i in range(len(A)-3,0,-1):
            g = gcd(B[0],A[i])
            B[0] = B[0]//g
            A[i] = A[i]//g
            B.insert(0,g)
    B.append(A[-1])
    solutions = []
    c = Integer(c)
    for Ai,Bi in zip(A,B):
        tot_x,tot_y = [],[]
        for arg in Add.make_args(c):
            if arg.is_Integer:
                k,p,pnew = int(arg),Integer(1),params[0]
            else:
                k,p = arg.as_coeff_Mul()
                k = int(k)
                pnew = params[params.index(p)+1]
            sol = _base_solution_linear(k,Ai,Bi,pnew)
            if sol is None:
                if p == 1:return None
                continue
            (nx,px),(ny,py) = sol
            tot_x.append(nx*p+px)
            tot_y.append(ny*p+py)
        solutions.append(Add(*tot_x))
        c = Add(*tot_y)
    solutions.append(c)
    return tuple(solutions)


def _pell_unit(D):
    a0 = isqrt(D)
    if a0*a0 == D:raise ValueError('Pell discriminant must not be a square')
    m,d,a = 0,1,a0
    p0,p1,q0,q1 = 1,a,0,1
    while p1*p1-D*q1*q1 != 1:
        m = d*a-m;d = (D-m*m)//d;a = (a0+m)//d
        p0,p1,q0,q1 = p1,a*p1+p0,q1,a*q1+q0
    return p1,q1


def _canonical_signs(solution):
    # reference `DiophantineSolutionSet.add`: no `-t` unless `t` is present
    from . import sympify
    solution = [sympify(s) for s in solution]
    args = set(solution)
    for i in range(len(solution)):
        value = solution[i]
        if (-value).is_Symbol and -value not in args:
            solution = [s.subs(-value,value) for s in solution]
    return tuple(solution)


def _diop_parabolic(expr,var,name):
    """Reference `BinaryQuadratic.solve`, parabolic case `B**2 = 4*A*C`.

    `expr` is an irreducible integer quadratic in `var` (sorted as in the
    reference); the solutions are returned in the order of `var`.
    """
    from . import Poly,Symbol,Integer,S,solveset_real,Rational
    x,y = var
    P = Poly(expr,x,y)
    A,B,C,D,E,F = (int(P.coeff_monomial(mo)) for mo in (x**2,x*y,y**2,x,y,1))
    g0 = gcd(gcd(gcd(A,B),gcd(C,D)),gcd(E,F))
    if g0 > 1:
        A,B,C,D,E,F = (v//g0 for v in (A,B,C,D,E,F))
    if A == 0:
        return {(s[1],s[0]) for s in _diop_parabolic(expr,(y,x),name)}
    t = Symbol(name,integer=True)
    sgn = lambda v: (v > 0)-(v < 0)
    g = sgn(A)*gcd(A,C)
    a,c = A//g,C//g
    e = sgn(B)*sgn(A)
    sqa,sqc = isqrt(a),isqrt(c)
    _c = e*sqc*D-sqa*E
    result = set()
    if not _c:
        z = Symbol('z',real=True)
        roots = solveset_real(sqa*g*z**2+D*z+sqa*F,z).intersect(S.Integers)
        for root in roots:
            for ans in diophantine(sqa*x+e*sqc*y-root,syms=(x,y)):
                result.add(_canonical_signs(ans))
        return result
    for u in range(0,abs(_c)):
        if (sqa*g*u*u+D*u+sqa*F)%_c == 0 and (e*sqc*g*u*u+E*u+e*sqc*F)%_c == 0:
            sx = -e*sqc*g*_c*t**2-(E+2*e*sqc*g*u)*t-Integer((e*sqc*g*u*u+E*u+e*sqc*F)//_c)
            sy = sqa*g*_c*t**2+(D+2*sqa*g*u)*t+Integer((sqa*g*u*u+D*u+sqa*F)//_c)
            result.add(_canonical_signs((sx,sy)))
    return result


def diophantine(eq,param=None,syms=None,permute=False):
    from . import sympify,Symbol,Poly,Add,expand,divisors,sqrt,Integer,factor_terms,default_sort_key
    eq = sympify(eq.lhs-eq.rhs if getattr(eq,'is_Equality',False) else eq)
    variables = tuple(sorted(eq.free_symbols,key=str)) if syms is None else tuple(syms)
    if not variables:return {()} if eq == 0 else set()
    if any(not v.is_Symbol for v in variables):raise TypeError('Diophantine variables must be symbols')
    if eq.free_symbols-set(variables):raise ValueError('all variables must be included in syms')
    # reference: the numerator, with the factors free of the variables
    # removed (`2*pi*n - pi*m` is `2*n - m`)
    num,den = eq.as_numer_denom()
    if den.is_number and num.has(*variables):
        reduced = factor_terms(num).as_independent(*variables,as_Add=False)[1]
        if reduced.free_symbols <= set(variables):
            eq = reduced
    p = Poly(eq,*variables)
    if any(c.is_Rational is not True for c in p.coeffs()):raise TypeError('rational polynomial coefficients are required')
    _,p = p.clear_denoms(convert=True)
    name = 't' if param is None else str(param)
    parameters = [Symbol('%s_%s' % (name,i),integer=True) for i in range(max(2,len(variables)))]
    if int(p.total_degree()) <= 1:
        present = set(p.as_expr().free_symbols)
        if len(variables) > 1 and present == set(variables):
            # reference `Linear.solve`, in the default sorted order of the
            # variables and with the sign normalization of `factor_list`
            # and `signsimp`
            order = tuple(sorted(variables,key=default_sort_key))
            base = Poly(p.as_expr(),*order)
            base = base.quo_ground(base.content()) if base.content() != 0 else base
            if base.LC() < 0:base = -base
            base = base.as_expr()
            if base.could_extract_minus_sign():base = -base
            coeffs = [int(base.diff(v)) for v in order]
            target = -int(base.subs(dict.fromkeys(order,0)))
            row = _linear_reference(coeffs,target,[Symbol('%s_%s' % (name,i),integer=True) for i in range(len(order))])
            if row is None:return set()
            row = _canonical_signs(row)
            return {tuple(row[order.index(v)] for v in variables)}
        coeffs = [int(p.as_expr().diff(v)) for v in variables]
        target = -int(p.as_expr().subs(dict.fromkeys(variables,0)))
        return _linear(coeffs,target,parameters)
    if len(variables) == 1:
        return {(r,) for r in p.ground_roots() if r.is_Integer}
    factors = p.factor_list()[1]
    if len(factors) > 1 or (factors and factors[0][1] > 1):
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
        if int(p.total_degree()) == 2 and b*b == 4*a*c:
            # reference parabolic case, in the default sorted order of the
            # variables and with the sign normalization of `factor_list`
            # and `signsimp`
            order = tuple(sorted(variables,key=default_sort_key))
            base = Poly(expression,*order)
            if base.LC() < 0:base = -base
            base = base.as_expr()
            if base.could_extract_minus_sign():base = -base
            name = 't' if param is None else str(param)
            rows = _diop_parabolic(base,order,name)
            return {tuple(row[order.index(v)] for v in variables) for row in rows}
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
