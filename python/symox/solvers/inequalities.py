"""Real polynomial, rational and absolute-value inequalities using native sets."""
from functools import cmp_to_key

__all__=['solve_poly_inequality','solve_poly_inequalities',
         'solve_rational_inequalities','reduce_rational_inequalities','reduce_abs_inequality','reduce_abs_inequalities',
         'solve_univariate_inequality','reduce_inequalities']


def _relation(expression,relation):
    from .. import Relational
    if relation not in ('==','!=','<','<=','>','>='):
        raise ValueError("'%s' is not a valid relation" % relation)
    return Relational(expression,0,relation)


def _polynomial(poly):
    from .. import Poly
    if not isinstance(poly,Poly):raise ValueError('poly must be a Poly instance')
    if len(poly.gens)!=1:raise ValueError('a univariate polynomial is required')
    if any(not coefficient.is_Rational for coefficient in poly.all_coeffs()):
        raise NotImplementedError('Polynomial inequalities require rational coefficients')
    if poly.modulus is not None:raise NotImplementedError('Finite fields have no ordered real inequality domain')
    return poly


def _compare(left,right):
    from .. import S
    if left==right:return 0
    if (left<right)==S.true:return -1
    if (left>right)==S.true:return 1
    raise NotImplementedError('Real algebraic endpoints could not be ordered')


def _intervals(domain):
    from .. import S,FiniteSet,Interval,oo
    if domain==S.EmptySet:return []
    if domain==S.Reals:return [Interval(-oo,oo)]
    parts=[]
    for component in domain.args if domain.is_Union else (domain,):
        if component.is_FiniteSet:
            parts.extend(FiniteSet(point) for point in sorted(component.args,key=cmp_to_key(_compare)))
        elif isinstance(component,Interval):parts.append(component)
        else:raise NotImplementedError('Polynomial inequality did not produce finite real intervals')
    def endpoint(part):return part.args[0] if part.is_FiniteSet else part.start
    return sorted(parts,key=cmp_to_key(lambda left,right:_compare(endpoint(left),endpoint(right))))


def solve_poly_inequality(poly,rel):
    """Return real interval/singleton components for one rational polynomial."""
    from .. import S,solve_univariate_inequality
    poly=_polynomial(poly)
    condition=_relation(poly.as_expr(),rel)
    if poly.degree()<=0:
        if condition==S.true:return [S.Reals]
        if condition==S.false:return [S.EmptySet]
        raise NotImplementedError('The constant inequality could not be decided')
    return _intervals(solve_univariate_inequality(condition,poly.gens[0],relational=False))


def solve_poly_inequalities(polys):
    from .. import Union
    return Union(*(part for polynomial,relation in polys
                   for part in solve_poly_inequality(polynomial,relation)))


def solve_rational_inequalities(eqs):
    """Union of conjunctions, excluding every original denominator zero."""
    from .. import S,Union
    result=S.EmptySet
    count=0
    for equations in eqs:
        domain=S.Reals;present=False
        for (numerator,denominator),relation in equations:
            present=True;count+=1
            if count>4096:raise NotImplementedError('Rational inequality systems are limited to 4096 conditions')
            numerator,denominator=_polynomial(numerator),_polynomial(denominator)
            if numerator.gens!=denominator.gens:raise ValueError('Rational inequalities require one common generator')
            allowed=Union(*solve_poly_inequality(numerator*denominator,relation))
            poles=Union(*solve_poly_inequality(denominator,'=='))
            domain=domain.intersect(allowed)-poles
            if domain==S.EmptySet:break
        if present:result=result.union(domain)
    return result


def _rational_domain(expression,relation,variable):
    from .. import Poly
    numerator,denominator=expression.as_numer_denom()
    return solve_rational_inequalities([[((Poly(numerator,variable),Poly(denominator,variable)),relation)]])


def reduce_rational_inequalities(exprs,gen,relational=True):
    """Reduce disjunctions of rational conditions, equations or expression pairs."""
    from .. import Float,Poly,Symbol,S,nfloat,solve_univariate_inequality,sympify,together
    gen=sympify(gen)
    if not isinstance(gen,Symbol):raise TypeError('The inequality variable must be a symbol')
    result=S.EmptySet;exact=True;count=0
    for expressions in exprs:
        domain=S.Reals;present=False
        for expression in expressions:
            present=True;count+=1
            if count>4096:raise NotImplementedError('Rational inequality systems are limited to 4096 conditions')
            if isinstance(expression,(tuple,list)):
                expression,relation=expression
                expression=sympify(expression)
            else:
                expression=sympify(expression)
                if expression.is_Relational:
                    relation=expression.rel_op
                    expression=expression.lhs-expression.rhs
                else:relation='=='
            if expression==S.true:expression=S.Zero;relation='=='
            elif expression==S.false:expression=S.One;relation='=='
            condition=_relation(expression,relation)
            numerator,denominator=together(expression).as_numer_denom()
            p,q=Poly(numerator,gen),Poly(denominator,gen)
            if expression.atoms(Float):
                p,q=p.to_exact(),q.to_exact();exact=False
            if all(coefficient.is_Rational for poly in (p,q) for coefficient in poly.all_coeffs()):
                allowed=solve_rational_inequalities([[((p,q),relation)]])
            else:
                allowed=solve_univariate_inequality(condition,gen,relational=False)
            domain=domain.intersect(allowed)
            if domain==S.EmptySet:break
        if present:result=result.union(domain)
    if not exact:result=nfloat(result)
    return result.as_relational(gen) if relational else result


def _abs_domain(expr,rel,gen):
    from .. import Symbol, Dummy, S, Union, sympify
    from .._core import Expr
    gen=sympify(gen);expression=sympify(expr)
    if not isinstance(gen,Symbol):raise TypeError('The inequality variable must be a symbol')
    if gen.is_extended_real is False:raise TypeError('Absolute-value inequalities require a real variable')
    if not isinstance(expression,Expr):raise TypeError('A scalar expression is required')
    _relation(expression,rel)
    variable=Dummy('abs_inequality_x',real=True)
    expression=expression.xreplace({gen:variable})
    memo={};owners={};nodes=0

    def visit(node,depth=0):
        nonlocal nodes
        key=node._dag_id
        if key in memo:return memo[key]
        nodes+=1
        if nodes>10000 or depth>128:raise NotImplementedError('Absolute-value inequality graph is too large')
        owners[key]=node
        if node.is_Add or node.is_Mul:
            addition=node.is_Add
            rows=[(S.Zero if addition else S.One,())]
            for argument in node.args:
                children=visit(argument,depth+1)
                if len(rows)*len(children)>1024:raise NotImplementedError('Absolute-value inequalities are limited to 1024 sign branches')
                rows=[(left+right if addition else left*right,conditions+extra)
                      for left,conditions in rows for right,extra in children]
        elif node.is_Pow:
            base,power=node.args
            if not power.is_Integer:raise ValueError('Only integer powers are supported in absolute-value inequalities')
            rows=[(value**power,conditions) for value,conditions in visit(base,depth+1)]
        elif node.is_Function and not node.is_UndefinedFunction and node.name=='Abs':
            children=visit(node.args[0],depth+1)
            if len(children)>512:raise NotImplementedError('Absolute-value inequalities are limited to 1024 sign branches')
            rows=[branch for value,conditions in children for branch in
                  [(value,conditions+((value,'>='),)),(-value,conditions+((value,'<'),))]]
        else:rows=[(node,())]
        memo[key]=rows
        return rows

    result=S.EmptySet;cache={}
    def allowed(value,relation):
        key=(value,relation)
        if key not in cache:cache[key]=_rational_domain(value,relation,variable)
        return cache[key]
    for value,conditions in visit(expression):
        domain=allowed(value,rel)
        for value,relation in conditions:
            if domain==S.EmptySet:break
            domain=domain.intersect(allowed(value,relation))
        result=result.union(domain)
    return result


def reduce_abs_inequality(expr,rel,gen):
    """Reduce nested real absolute values by exact rational sign branches."""
    return _abs_domain(expr,rel,gen).as_relational(gen)


def reduce_abs_inequalities(exprs,gen):
    from .. import And
    return And(*(reduce_abs_inequality(expression,relation,gen)
                 for expression,relation in exprs))


def _reduce_abs_system(inequalities,symbols):
    from .. import Abs, Eq, S, Symbol, solve_univariate_inequality, sympify
    from .._core import Expr
    conditions=[sympify(condition) for condition in inequalities] if isinstance(inequalities,(list,tuple,set)) else [sympify(inequalities)]
    if not any(condition.has(Abs) for condition in conditions):
        raise NotImplementedError('This is not an absolute-value inequality system')
    if symbols is None:
        variables=set().union(*(condition.free_symbols for condition in conditions))
    else:variables=[symbols] if isinstance(symbols,Expr) else list(symbols)
    if len(variables)!=1:raise NotImplementedError('Absolute-value inequality reduction requires one variable')
    variable=next(iter(variables))
    if not isinstance(variable,Symbol):raise TypeError('The inequality variable must be a symbol')
    result=S.Reals
    for condition in conditions:
        if condition in (S.true,S.false):
            if condition==S.false:return S.false
            continue
        if not condition.is_Relational:condition=Eq(condition,0)
        if condition.is_Relational and condition.has(Abs):
            allowed=_abs_domain(condition.lhs-condition.rhs,condition.rel_op,variable)
        else:allowed=solve_univariate_inequality(condition,variable,relational=False)
        result=result.intersect(allowed)
        if result==S.EmptySet:break
    return result.as_relational(variable)


def solve_univariate_inequality(*args,**kwargs):
    from .. import solve_univariate_inequality as native
    return native(*args,**kwargs)


def reduce_inequalities(*args,**kwargs):
    from .. import reduce_inequalities as native
    return native(*args,**kwargs)


def _solve_inequality(ie,s,linear=False):
    """Port of SymPy's ``_solve_inequality``: isolate ``s`` on the left when
    that is safe, dividing only by a factor of known sign (or nonzero for
    Eq/Ne), and return the conditions the inversion needs."""
    from .. import (And,Eq,Poly,PolynomialError,S,factor_terms,
                    solve_univariate_inequality)
    from .solvers import denoms
    if s not in ie.free_symbols:
        return ie
    if ie.rhs==s:
        ie=ie.reversed
    if ie.lhs==s and s not in ie.rhs.free_symbols:
        return ie

    def classify(ie,s,i):
        try:
            v=ie.subs(s,i)
            if v is S.NaN:
                return v
            elif v not in (True,False):
                return
            return v
        except TypeError:
            return S.NaN

    rv=None
    oo=S.Infinity
    expr=ie.lhs-ie.rhs
    try:
        p=Poly(expr,s)
        if p.degree()==0:
            rv=ie.func(p.as_expr(),0)
        elif not linear and p.degree()>1:
            raise NotImplementedError
    except (PolynomialError,NotImplementedError):
        if not linear:
            try:
                rv=reduce_rational_inequalities([[ie]],s)
            except PolynomialError:
                rv=solve_univariate_inequality(ie,s)
            okoo=classify(ie,s,oo)
            if okoo is S.true and classify(rv,s,oo) is S.false:
                rv=rv.subs(s<oo,True)
            oknoo=classify(ie,s,-oo)
            if oknoo is S.true and classify(rv,s,-oo) is S.false:
                rv=rv.subs(-oo<s,True)
                rv=rv.subs(s>-oo,True)
            if rv is S.true:
                rv=(s<=oo) if okoo is S.true else (s<oo)
                if oknoo is not S.true:
                    rv=And(-oo<s,rv)
        else:
            p=Poly(expr)

    conds=[]
    if rv is None:
        e=p.as_expr()
        rhs=0
        b,ax=e.as_independent(s,as_Add=True)
        e-=b
        rhs-=b
        ef=factor_terms(e)
        a,e=ef.as_independent(s,as_Add=False)
        if (a.is_zero!=False or
                a.is_negative==a.is_positive is None and
                ie.rel_op not in ('!=','==')):
            e=ef
            a=S.One
        rhs/=a
        if a.is_positive:
            rv=ie.func(e,rhs)
        else:
            rv=ie.reversed.func(e,rhs)
        beginning_denoms=denoms(ie.lhs)|denoms(ie.rhs)
        current_denoms=denoms(rv)
        for d in beginning_denoms-current_denoms:
            c=_solve_inequality(Eq(d,0),s,linear=linear)
            if isinstance(c,Eq) and c.lhs==s:
                if classify(rv,s,c.rhs) is S.true:
                    conds.append(~c)
        for i in (-oo,oo):
            if (classify(rv,s,i) is S.true and
                    classify(ie,s,i) is not S.true):
                conds.append(s<i if i is oo else i<s)

    conds.append(rv)
    return And(*conds)


def _reduce_multisymbol(inequalities,symbols,native):
    """SymPy's ``reduce_inequalities`` when some inequality has more than one
    symbol, or the univariate inequalities are in several symbols: an
    inequality with exactly one symbol of interest is isolated by
    ``_solve_inequality``; the univariate ones keep the native reduction,
    grouped by generator. Returns None when no inequality needs this path."""
    from .. import And,Dummy,Relational,S,Symbol,sympify
    gens=set().union(*[i.free_symbols for i in inequalities])
    symbols=(set(symbols) or gens)&gens
    multi=False
    single=set()
    for i in inequalities:
        if isinstance(i,Relational):
            found=(i.lhs-i.rhs).atoms(Symbol)
            if len(found)>1:
                multi=True
            single|=found
    # univariate inequalities in different symbols are reduced per symbol
    if not multi and len(single)<2:
        return None
    if any(i.is_extended_real is False for i in symbols):
        raise TypeError('inequalities cannot contain symbols that are not real.')
    recast={i:Dummy(i.name,extended_real=True) for i in gens if i.is_extended_real is None}
    back={v:k for k,v in recast.items()}
    rsymbols={i.xreplace(recast) for i in symbols}
    keep=[]
    for i in inequalities:
        i=i.xreplace(recast)
        if isinstance(i,Relational):
            i=i.func(i.lhs.as_expr()-i.rhs.as_expr(),0)
        elif i not in (True,False):
            i=Relational(i,0,'==')
        if i==True:
            continue
        elif i==False:
            return S.false
        if i.lhs.is_number:
            raise NotImplementedError('could not determine truth value of %s' % i)
        keep.append(i)
    groups={};other=[]
    for i in keep:
        expr=i.lhs
        found=expr.atoms(Symbol)
        if len(found)==1:
            gen=found.pop()
            groups.setdefault(gen,[]).append(i.xreplace(back))
            continue
        common=expr.free_symbols&rsymbols
        if len(common)!=1:
            raise NotImplementedError('\ninequality has more than one symbol of interest.')
        other.append(_solve_inequality(i,common.pop()))
    parts=[native(group,[gen.xreplace(back)]) for gen,group in groups.items()]
    return And(*(parts+[o.xreplace(back) for o in other]))
