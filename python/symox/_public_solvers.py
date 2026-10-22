"""Factor-system decomposition, triangular solving and integer equations."""
from itertools import product
from math import gcd,isqrt


def factor_system(eqs, gens=(), **kwargs):
    """Factorise a system of polynomial equations into the irreducible
    generic subsystems (SymPy 1.14 ``factor_system``)."""
    systems = _factor_system_poly_from_expr(eqs, gens, **kwargs)
    systems_generic = [sys for sys in systems if not _is_degenerate(sys)]
    return [[p.as_expr() for p in system] for system in systems_generic]


def _is_degenerate(system):
    return any(p.is_ground for p in system)


def factor_system_bool(eqs, gens=(), **kwargs):
    """``factor_system_cond`` as a Boolean: an Or of Ands of equations."""
    from . import Or, And, Eq
    systems = factor_system_cond(eqs, gens, **kwargs)
    return Or(*[And(*[Eq(eq, 0) for eq in sys]) for sys in systems])


def factor_system_cond(eqs, gens=(), **kwargs):
    """Like ``factor_system`` but keeps the subsystems that only hold
    under a condition on the parameters."""
    systems_poly = _factor_system_poly_from_expr(eqs, gens, **kwargs)
    return [[p.as_expr() for p in system] for system in systems_poly]


def _factor_system_poly_from_expr(eqs, gens, **kwargs):
    from . import Dummy
    from .polys.polytools import parallel_poly_from_expr
    from .polys.polyerrors import GeneratorsNeeded, PolificationFailed
    try:
        polys, opts = parallel_poly_from_expr(eqs, *gens, **kwargs)
        only_numbers = False
    except (GeneratorsNeeded, PolificationFailed):
        _u = Dummy('u')
        polys, opts = parallel_poly_from_expr(eqs, [_u], **kwargs)
        only_numbers = True
    if only_numbers:
        return [[]] if all(p == 0 for p in polys) else []
    return factor_system_poly(polys)


def factor_system_poly(polys):
    """Factorise a system of ``Poly`` equations (SymPy 1.14 algorithm)."""
    from . import Poly
    from .core.exprtools import factor_terms
    from .polys.polytools import sqf_part
    from .polys.polyerrors import DomainError
    if not all(isinstance(poly, Poly) for poly in polys):
        raise TypeError("polys should be a list of Poly instances")
    if not polys:
        return [[]]
    base_domain = polys[0].domain
    base_gens = polys[0].gens
    if not all(poly.domain == base_domain and poly.gens == base_gens for poly in polys[1:]):
        raise DomainError("All polynomials must have the same domain and generators")
    factor_sets = []
    for poly in polys:
        constant, factors_mult = poly.factor_list()
        if constant.is_zero is True:
            continue
        elif constant.is_zero is False:
            if not factors_mult:
                return []
            factor_sets.append([f for f, _ in factors_mult])
        else:
            constant = sqf_part(factor_terms(constant).as_coeff_Mul()[1])
            constp = Poly(constant, base_gens, domain=base_domain)
            factors = [f for f, _ in factors_mult]
            factors.append(constp)
            factor_sets.append(factors)
    if not factor_sets:
        return [[]]
    return _sort_systems(_factor_sets(factor_sets))


def _factor_sets(eqs):
    if not eqs:
        return {frozenset()}
    current_set = min(eqs, key=len)
    other_sets = [s for s in eqs if s is not current_set]
    stack = [(factor, [s for s in other_sets if factor not in s], {factor})
             for factor in current_set]
    result = set()
    while stack:
        factor, remaining_sets, current_solution = stack.pop()
        if not remaining_sets:
            result.add(frozenset(current_solution))
            continue
        next_set = min(remaining_sets, key=len)
        next_remaining = [s for s in remaining_sets if s is not next_set]
        for next_factor in next_set:
            valid_remaining = [s for s in next_remaining if next_factor not in s]
            stack.append((next_factor, valid_remaining, current_solution | {next_factor}))
    return {s1 for s1 in result if not any(s1 > s2 for s2 in result)}


def _sort_systems(systems):
    systems_list = [sorted(s, key=_poly_sort_key, reverse=True) for s in systems]
    return sorted(systems_list, key=_sys_sort_key, reverse=True)


def _poly_sort_key(poly):
    if poly.domain.is_FF:
        from .polys.domains import ZZ
        poly = poly.set_domain(ZZ)
    domain = poly.domain
    return poly.degree_list(), _rep_key(domain, poly.rep.to_list())


def _rep_key(domain, rep):
    if isinstance(rep, list):
        return [_rep_key(domain, c) for c in rep]
    return _coeff_key(domain, rep)


def _coeff_key(domain, c):
    # The comparison keys of the reference domain elements: PolyElement
    # orders by (number of terms, terms), FracElement by (denominator,
    # numerator), ANP by its coefficient list.
    from fractions import Fraction
    from . import Poly, fraction
    if domain.is_PolynomialRing:
        terms = Poly(domain.to_sympy(c), *domain.symbols, domain=domain.dom).terms()
        terms = [(m, v) for m, v in terms if v != 0]
        return (len(terms), terms)
    if domain.is_FractionField:
        numer, denom = fraction(domain.to_sympy(c))
        def key(e):
            terms = [(m, v) for m, v in Poly(e, *domain.symbols, domain=domain.dom).terms() if v != 0]
            return (len(terms), terms)
        return (key(denom), key(numer))
    if domain.is_AlgebraicField:
        return [Fraction(int(q.numerator), int(q.denominator)) for q in c.to_list()]
    return domain.to_sympy(c)


def _sys_sort_key(sys):
    return list(zip(*map(_poly_sort_key, sys)))


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


def _pow_rep_recursive(n_i, k, n_remaining, terms, p):
    # SymPy 1.14 ``pow_rep_recursive`` (BSD-3-Clause, see LICENSE-SYMPY).
    from .ntheory.primetest import is_square  # noqa: F401
    from . import integer_nthroot
    if n_i <= 0 or k <= 0:
        return
    if n_remaining < k:
        return
    if k * pow(n_i, p) < n_remaining:
        return
    if k == 0 and n_remaining == 0:
        yield tuple(terms)
    elif k == 1:
        next_term, exact = integer_nthroot(n_remaining, p)
        if exact and next_term <= n_i:
            yield tuple(terms + [next_term])
        return
    else:
        if n_i >= 1 and k > 0:
            for next_term in range(1, n_i + 1):
                residual = n_remaining - pow(next_term, p)
                if residual < 0:
                    break
                yield from _pow_rep_recursive(next_term, k - 1, residual, terms + [next_term], p)


def _even_power_representation(n, p, k):
    """SymPy's ``power_representation(n, p, k)`` without zeros, for an even
    exponent ``p > 2`` and ``k >= 2`` terms."""
    from . import integer_nthroot, perfect_power
    if n <= 0:
        return
    if k == 2:
        be = perfect_power(n)
        if be and be[1] % p == 0:
            return  # Fermat: a**n + b**n = c**n has no solution for n > 2
    if n >= k:
        a = integer_nthroot(n - (k - 1), p)[0]
        for t in _pow_rep_recursive(a, k, n, [], p):
            yield tuple(reversed(t))


def _diop_sum_of_even_powers(p, variables, permute):
    """SymPy 1.14's ``general_sum_of_even_powers`` family
    ``x_1**e + ... + x_k**e = n`` with an even ``e > 3``; None when the
    polynomial is not of that form."""
    from . import Integer
    from .utilities.iterables import signed_permutations
    degree = int(p.total_degree())
    if len(variables) < 2 or degree <= 3 or degree % 2:
        return None
    if p.LC() < 0:
        p = -p
    constant = 0
    powers = set()
    for monomial, coeff in p.terms():
        if not any(monomial):
            constant = int(coeff)
            continue
        if coeff != 1 or sorted(monomial)[-1] != degree or sum(monomial) != degree:
            return None
        powers.add(monomial.index(degree))
    if len(powers) != len(variables):
        return None
    if not constant:
        # the reference's trivial solution of a sum of even powers
        return {(Integer(0),)*len(variables)}
    n = -constant
    sign = [-1 if v.is_nonpositive else 1 for v in variables]
    negs = sign.count(-1) != 0
    result = set()
    for t in _even_power_representation(n, degree, len(variables)):
        row = tuple(Integer(sign[i]*j) if negs else Integer(j) for i, j in enumerate(t))
        if permute:
            result.update(signed_permutations(row))
        else:
            result.add(row)
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
    even = _diop_sum_of_even_powers(factors[0][0] if factors else p,variables,permute)
    if even is not None:
        return even
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
