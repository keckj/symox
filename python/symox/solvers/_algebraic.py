"""Native denominator inspection, small determinants and radical elimination."""
import itertools
import math


def denoms(eq,*symbols):
    from .. import sympify,Mul,Expr
    from ._differential import nodes
    if len(symbols)==1 and isinstance(symbols[0],(list,tuple,set,frozenset)):
        symbols=tuple(symbols[0])
    symbols=set(map(sympify,symbols))
    expressions=list(eq) if isinstance(eq,(list,tuple,set,frozenset)) else [eq]
    result=set()
    for expression in expressions:
        for node in nodes(sympify(expression)):
            if not isinstance(node,Expr):continue
            denominator=node.as_numer_denom()[1]
            if denominator!=1:
                result.update(Mul.make_args(denominator))
    return result if not symbols else {value for value in result if value.free_symbols&symbols}


def recast_to_symbols(eqs,symbols):
    from .. import Dummy,sympify,Expr,default_sort_key
    if isinstance(eqs,Expr) or isinstance(symbols,Expr):
        raise ValueError('Equations and symbols must both be iterable')
    equations=list(eqs);symbols=list(map(sympify,symbols))
    replacements={}
    for value in sorted(set(symbols),key=default_sort_key):
        if not value.is_Symbol:
            replacements[value]=Dummy('X'+str(len(replacements)))
    changed=[sympify(value).subs(replacements,simultaneous=True) for value in equations]
    return changed,[replacements.get(value,value) for value in symbols],{v:k for k,v in replacements.items()}


def _square(M):
    from .. import Matrix,NonSquareMatrixError
    M=Matrix(M)
    if M.rows!=M.cols:raise NonSquareMatrixError('A square matrix is required')
    return M


def det_perm(M):
    """Leibniz determinant, with a native elimination fallback above order eight."""
    from .. import S,Add,Mul
    M=_square(M);n=M.rows
    if n>8:return M.det()
    terms=[]
    for permutation in itertools.permutations(range(n)):
        term=Mul(*(M[i,j] for i,j in enumerate(permutation)))
        if term==0:continue
        odd=sum(permutation[i]>permutation[j] for i in range(n) for j in range(i+1,n))%2
        terms.append(-term if odd else term)
    return Add(*terms)


def det_minor(M):
    """Laplace expansion with memoized column subsets and distributed terms."""
    from .. import S,Add
    M=_square(M);n=M.rows
    if n>12:return M.det()
    cache={():S.One}
    def minor(columns):
        if columns in cache:return cache[columns]
        row=n-len(columns)
        terms=[]
        for offset,column in enumerate(columns):
            coefficient=M[row,column]*(-1 if offset%2 else 1)
            if coefficient==0:continue
            value=minor(columns[:offset]+columns[offset+1:])
            terms.extend(coefficient*term for term in Add.make_args(value))
        value=Add(*terms)
        cache[columns]=value
        return value
    return minor(tuple(range(n)))


def det_quick(M,method=None):
    M=_square(M)
    if method is not None:return M.det(method=method)
    if any(value.free_symbols for value in M):
        return det_perm(M) if M.rows<8 and all(value.free_symbols for value in M) else det_minor(M)
    return M.det()


def inv_quick(M):
    """Exact cofactor inverse for small symbolic matrices; native kernels otherwise."""
    from .. import Matrix,S,NonInvertibleMatrixError,cancel
    M=_square(M);n=M.rows
    if n>6 or all(value.is_Number for value in M):return M.inv()
    determinant=det_quick(M)
    if determinant==0:raise NonInvertibleMatrixError('Matrix determinant is zero')
    if n==0:return Matrix(0,0,[])
    return Matrix(n,n,[cancel((-1)**(i+j)*det_quick(M.minor_submatrix(j,i))/determinant)
                       for i in range(n) for j in range(n)])


def minsolve_linear_system(system,*symbols,**flags):
    """Find a particular solution with the fewest nonzero unknowns.

    Exact search enumerates supports; the quick path greedily removes columns
    while retaining consistency. Nontrivial homogeneous solutions are retained
    when the reference linear-system convention permits free variables.
    """
    from .. import Matrix,S,solve_linear_system,cancel
    from .._linear import _symbols
    quick=flags.pop('quick',False)
    if flags:raise NotImplementedError('Unknown minimal linear-system options')
    matrix=Matrix(system);symbols=_symbols(symbols)
    if matrix.cols!=len(symbols)+1:
        raise ValueError('Augmented columns must match the symbols plus one')
    general=solve_linear_system(matrix,*symbols)
    if general is None or not general or all(value==0 for value in general.values()):return general
    n=len(symbols)
    if n>64:raise ValueError('Minimal linear solving accepts at most 64 unknowns')
    homogeneous=all(value==0 for value in matrix[:,-1:])
    def candidate(support):
        columns=list(support)+[n]
        submatrix=matrix[:,columns]
        chosen=[symbols[i] for i in support]
        solved=solve_linear_system(submatrix,*chosen)
        if solved is None:return None
        free=[symbol for symbol in chosen if symbol not in solved]
        assignments={symbol:S.Zero for symbol in free}
        if homogeneous and free:assignments[free[0]]=S.One
        values={symbol:S.Zero for symbol in symbols}
        values.update(assignments)
        values.update({symbol:cancel(value.subs(assignments,simultaneous=True)) for symbol,value in solved.items()})
        if homogeneous and not any(value!=0 for value in values.values()):return None
        return values
    support=list(range(n));best=candidate(support)
    for column in reversed(range(n)):
        reduced=[i for i in support if i!=column]
        found=candidate(reduced)
        if found is not None:support,best=reduced,found
    if quick:return best
    if n>20:raise NotImplementedError('Exact sparse support search is limited to 20 unknowns; use quick=True')
    for size in range(1,sum(value!=0 for value in best.values())):
        for subset in itertools.combinations(range(n),size):
            found=candidate(subset)
            if found is not None:return found
    return best


def _radical_normalize(expression):
    from .. import factor_terms,Poly,Mul
    expression=factor_terms(expression.expand(),clear=True)
    if expression.is_Mul:
        expression=Mul(*(value for value in expression.args if not value.is_number))
    if expression.free_symbols:
        try:
            polynomial=Poly(expression,*sorted(expression.free_symbols,key=str))
            _,polynomial=polynomial.primitive()
            expression=polynomial.as_expr()
            if polynomial.LC().is_negative is True:expression=-expression
        except (ValueError,TypeError,NotImplementedError):pass
    return expression


def unrad(eq,*syms,**flags):
    """Eliminate rational-power radicals with exact polynomial resultants.

    The returned equation has a superset of the original roots. A shared-base
    change of variable retains a low-degree polynomial when possible. Root
    candidates still require substitution into the original equation.
    """
    from .. import (sympify,Expr,Equality,Dummy,S,resultant,Poly,cancel,solve,
                    powsimp)
    from ._differential import nodes
    take=flags.pop('take',None)
    if flags:raise NotImplementedError('Internal unrad recursion options are not supported')
    expression=sympify(eq)
    if not isinstance(expression,Expr):return None
    if isinstance(expression,Equality):expression=expression.lhs-expression.rhs
    if len(syms)==1 and isinstance(syms[0],(list,tuple,set,frozenset)):syms=tuple(syms[0])
    symbols=set(map(sympify,syms)) if syms else expression.free_symbols
    expression=powsimp(expression).as_numer_denom()[0].expand()
    radicals=list(dict.fromkeys(node for node in nodes(expression) if node.is_Pow and node.exp.is_Rational
              and node.exp.q!=1 and node.base.free_symbols&symbols and (take is None or take(node))))
    if not radicals:return None
    if len(radicals)>8 or any(node.exp.q>64 for node in radicals):
        raise NotImplementedError('Radical elimination exceeds its generator or degree bound')
    # Multiple rational powers of one base share a root coordinate. This
    # preserves the covariant relation rather than treating their branches as
    # independent algebraic roots, which can lose useful information.
    bases={node.base for node in radicals}
    if len(radicals)>1 and len(bases)==1:
        base=next(iter(bases));degree=math.lcm(*(int(node.exp.q) for node in radicals))
        if degree>64:raise NotImplementedError('Shared radical degree exceeds 64')
        coordinate=Dummy('p')
        replacement={node:coordinate**int(node.exp*degree) for node in radicals}
        transformed=expression.xreplace(replacement)
        if transformed.free_symbols&symbols:
            variables=base.free_symbols&symbols
            if len(variables)==1:
                variable=next(iter(variables))
                inverse=solve(base-coordinate**degree,variable)
                if len(inverse)==1:transformed=transformed.subs(variable,inverse[0])
        if not transformed.free_symbols&symbols:
            try:Poly(transformed,coordinate)
            except (ValueError,TypeError,NotImplementedError) as error:
                raise NotImplementedError('The radical coordinate does not give a polynomial') from error
            return _radical_normalize(transformed),[coordinate,coordinate**degree-base]
    if math.prod(int(node.exp.q) for node in radicals)>256:
        raise NotImplementedError('Combined radical elimination degree exceeds 256')
    # Build definitions from inner radicals outward, then eliminate in reverse
    # dependency order. Resultants remove every algebraic root coordinate.
    ordered=sorted(radicals,key=lambda node:sum(node.has(other) for other in radicals))
    replacements={};definitions=[]
    for node in ordered:
        coordinate=Dummy('radical')
        base=node.base.xreplace(replacements)
        power=int(node.exp.p);degree=int(node.exp.q)
        definition=coordinate**degree-base**power
        replacements[node]=coordinate
        definitions.append((coordinate,definition))
    transformed=expression.xreplace(replacements)
    for coordinate,definition in reversed(definitions):
        if not transformed.has(coordinate):continue
        try:
            transformed=resultant(transformed.as_numer_denom()[0],definition,coordinate)
        except (ValueError,TypeError,NotImplementedError) as error:
            raise NotImplementedError('Radical elimination requires polynomial root coordinates') from error
        transformed=cancel(transformed).as_numer_denom()[0]
    if any(node.is_Pow and node.exp.is_Rational and node.exp.q!=1 and node.base.free_symbols&symbols
           and (take is None or take(node))
           for node in nodes(transformed)):
        raise NotImplementedError('The result still contains radicals in the requested variables')
    return _radical_normalize(transformed),[]


__all__=['denoms','recast_to_symbols','det_perm','det_minor','det_quick','inv_quick',
         'minsolve_linear_system','unrad']
