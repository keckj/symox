"""Native denominator inspection, small determinants and radical elimination."""
import itertools


def denoms(eq,*symbols):
    from .. import sympify, Mul, fraction
    from .._core import Expr
    from ._differential import nodes
    if len(symbols)==1 and isinstance(symbols[0],(list,tuple,set,frozenset)):
        symbols=tuple(symbols[0])
    symbols=set(map(sympify,symbols))
    expressions=list(eq) if isinstance(eq,(list,tuple,set,frozenset)) else [eq]
    result=set()
    for expression in expressions:
        for node in nodes(sympify(expression)):
            if not isinstance(node,Expr):continue
            # the reference `denom`: `fraction(node)[1]`, which keeps a sum
            # over its own denominator (its terms are visited on their own)
            denominator=fraction(node)[1]
            if denominator!=1:
                result.update(Mul.make_args(denominator))
    return result if not symbols else {value for value in result if value.free_symbols&symbols}


def recast_subs(value,replacements):
    """``value.subs(replacements)`` as the reference: a Number target also
    replaces its negative (``Number._eval_subs``), so ``x - 2`` becomes
    ``x - X0`` when 2 is replaced by X0."""
    if not hasattr(value,'subs'):
        return value
    numbers={}
    for old,new in replacements.items():
        if old.is_Number:
            numbers[old]=new
            numbers.setdefault(-old,-new)
    others={old:new for old,new in replacements.items() if not old.is_Number}
    if numbers:
        value=value.xreplace(numbers)
    return value.subs(others) if others else value


def recast_to_symbols(eqs,symbols,keep=None):
    """SymPy's recast_to_symbols: the non-Symbol unknowns, in canonical
    order, become Dummy('X0'), Dummy('X1'), ...; ``keep`` exempts the
    unknowns a caller handles itself."""
    from .. import Dummy, sympify, ordered
    from .._core import Expr
    if isinstance(eqs,Expr) or isinstance(symbols,Expr):
        raise ValueError('Equations and symbols must both be iterable')
    equations=list(eqs);symbols=list(map(sympify,symbols))
    replacements={}
    for value in ordered(symbols):
        if not value.is_Symbol and value not in replacements and not (keep and keep(value)):
            replacements[value]=Dummy('X%d' % len(replacements))
    changed=[recast_subs(sympify(value),replacements) for value in equations]
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


def minsolve_linear_system(system, *symbols, **flags):
    """SymPy's ``minsolve_linear_system``: a particular solution with as few
    nonzero unknowns as possible (``quick=True`` uses the greedy heuristic)."""
    from .. import Matrix, S, solve_linear_system, ordered, default_sort_key, solve
    quick = flags.get('quick', False)
    system = Matrix(system)
    s0 = solve_linear_system(system, *symbols)
    if not s0 or all(v == 0 for v in s0.values()):
        return s0
    if quick:
        s = solve_linear_system(system, *symbols)

        def update(determined, solution):
            delete = []
            for k, v in solution.items():
                solution[k] = v.subs(determined)
                if not solution[k].free_symbols:
                    delete.append(k)
                    determined[k] = solution[k]
            for k in delete:
                del solution[k]
        determined = {}
        update(determined, s)
        while s:
            k = max((k for k in s.values()),
                    key=lambda x: (len(x.free_symbols), default_sort_key(x)))
            kfree = k.free_symbols
            x = next(reversed(list(ordered(kfree))))
            if len(kfree) != 1:
                determined[x] = S.Zero
            else:
                val = solve(k, x, check=False)[0]
                if not val and not any(v.subs(x, val) for v in s.values()):
                    determined[x] = S.One
                else:
                    determined[x] = val
            update(determined, s)
        return determined
    N = len(symbols)
    bestsol = minsolve_linear_system(system, *symbols, quick=True)
    n0 = len([x for x in bestsol.values() if x != 0])
    # SymPy stops at two nonzero values (range(n0 - 1, 1, -1)), so a
    # solution with a single nonzero value is never found although its
    # docstring promises the minimal number; the search includes n = 1.
    for n in range(n0 - 1, 0, -1):
        thissol = None
        for nonzeros in itertools.combinations(range(N), n):
            subm = Matrix([system.col(i).T for i in nonzeros] + [system.col(-1).T]).T
            s = solve_linear_system(subm, *[symbols[i] for i in nonzeros])
            if s and not all(v == 0 for v in s.values()):
                subs = [(symbols[v], S.One) for v in nonzeros]
                for k, v in s.items():
                    s[k] = v.subs(subs)
                for sym in symbols:
                    if sym not in s:
                        if symbols.index(sym) in nonzeros:
                            s[sym] = S.One
                        else:
                            s[sym] = S.Zero
                thissol = s
                break
        if thissol is None:
            break
        bestsol = thissol
    return bestsol


def _radical_normalize(expression):
    """SymPy's ``unrad._canonical``: clear the content, drop numeric factors
    and factor powers, then make each factor's sign canonical."""
    from .. import factor_terms,Mul
    expression=factor_terms(expression.expand(),clear=True)
    if expression.is_Mul:
        expression=Mul(*(value.base if value.is_Pow else value for value in expression.args if not value.is_number))
    factors=list(Mul.make_args(expression))
    changed=False
    for index,value in enumerate(factors):
        if value.could_extract_minus_sign():
            factors[index]=-value;changed=True
    if changed:expression=Mul(*factors,evaluate=False)
    return expression


def _unrad_resultant(eq,*syms,**flags):
    """Eliminate rational-power radicals with exact polynomial resultants.

    The established fallback of ``solve``/``solveset`` for the equations
    where SymPy's ``unrad`` raises NotImplementedError (it is not the public
    ``unrad``, whose results are SymPy's).

    The returned equation has a superset of the original roots. A shared-base
    change of variable retains a low-degree polynomial when possible. Root
    candidates still require substitution into the original equation.
    """
    from .. import sympify, Equality, Dummy, S, resultant, Poly, cancel, solve, powsimp
    from .._core import Expr
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
        base=next(iter(bases));degree=__import__('math').lcm(*(int(node.exp.q) for node in radicals))
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
    if __import__('math').prod(int(node.exp.q) for node in radicals)>256:
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


def unrad(eq, *syms, **flags):
    """SymPy's ``unrad``: remove the radicals with symbolic arguments of an
    equation, returning ``None`` or ``(eq, cov)`` where ``cov`` is empty or
    ``[p, relation]`` for a change of variable. The algorithm is native
    (``symox::solve::unrad``); the ``_vsolve`` calls and a ``_take``
    predicate are Python callbacks."""
    from .. import Equality, sympify
    from .._core import Expr, _unrad
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs  # XXX legacy Eq as Eqn support
    elif not isinstance(eq, Expr):
        return None
    reverse, take = bool(flags.get('_reverse', False)), flags.get('_take')
    syms = [sympify(s) for s in syms]
    if not any(k in flags for k in ('cov', 'n', 'rpt')):
        return _unrad(eq, syms, reverse, take)
    # the reference recursion flags (other flags are ignored, as there):
    # ``cov`` is the caller's list, updated in place
    from .._core import _unrad_flags
    cov, n, rpt = flags.get('cov', []), flags.get('n'), flags.get('rpt', 0)
    if not isinstance(cov, (list, tuple)) or len(cov) not in (0, 2) or \
            not (n is None or isinstance(n, int)) or not isinstance(rpt, int) or \
            isinstance(n, bool) or isinstance(rpt, bool) or (n is not None and n < 0) or rpt < 0:
        raise NotImplementedError('unrad flags cov=%r, n=%r, rpt=%r are not supported' % (cov, n, rpt))
    result, state = _unrad_flags(eq, syms, [sympify(c) for c in cov], n, rpt, reverse, take)
    if list(state) != list(cov):
        if not isinstance(cov, list):
            raise TypeError("'%s' object does not support item assignment" % type(cov).__name__)
        cov[:] = state
    return result



__all__=['denoms','recast_to_symbols','det_perm','det_minor','det_quick','inv_quick',
         'minsolve_linear_system','unrad']
