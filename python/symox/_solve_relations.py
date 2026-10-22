"""Route relational solve inputs through the native inequality engine."""
from functools import wraps
import sys


def install(namespace):
    native=namespace['solve']

    @wraps(native)
    def solve(f, *symbols, **flags):
        from . import Tuple, reduce_inequalities, sympify
        sequences=(list,tuple,set,frozenset,Tuple)
        if hasattr(f,'__next__'):
            f=list(f)
        if (len(symbols)==1 and isinstance(symbols[0],sequences) and not symbols[0]
                and not flags.get('set')):
            # SymPy solves for no unknown when given an empty collection
            return []
        equations=f if isinstance(f,sequences) else (f,)
        # SymPy sympifies every equation before inspecting it, so foreign
        # objects (Sage expressions) only reach symox methods
        from .core.sympify import SympifyError
        prepared=[];relational=False
        for equation in equations:
            try:
                value=sympify(equation)
            except SympifyError:
                value=equation
            prepared.append(value)
            if getattr(value,'is_Relational',False) and value.rel_op!='==':
                relational=True
        f=prepared if isinstance(f,sequences) else prepared[0]
        if relational:
            targets=symbols[0] if len(symbols)==1 else symbols
            return reduce_inequalities(prepared,targets)
        # SymPy's solve rewrites the hyperbolic functions of the unknowns in
        # terms of exp before solving (sinh(x) - 1 has two solutions).
        # SymPy splits a Matrix equation into its entries (and returns a
        # dict); the rewrites below see those scalar entries.
        matrix=any(getattr(value,'is_Matrix',False) for value in prepared)
        scalars=[entry for value in prepared
                 for entry in (value if getattr(value,'is_Matrix',False) else (value,))]
        if matrix:
            f=prepared=scalars
        sequence=isinstance(f,sequences)
        unknowns=[s for group in symbols
                  for s in (group if isinstance(group,sequences) else (group,))]
        if unknowns:
            from . import exp
            from .functions.elementary.hyperbolic import HyperbolicFunction
            hyperbolic=lambda w:isinstance(w,HyperbolicFunction) and w.has_free(*unknowns)
            if any(hasattr(value,'replace') and value.has(HyperbolicFunction)
                   and any(hyperbolic(w) for w in value.atoms(HyperbolicFunction))
                   for value in scalars):
                rewritten=[value.replace(hyperbolic,lambda w:w.rewrite(exp))
                           if hasattr(value,'replace') else value for value in scalars]
                target=rewritten if sequence else rewritten[0]
                return solve(target,*symbols,**flags)
        # SymPy's solve rationalizes the Floats of the equations unless
        # rational=False, and with rational=None (the default) it returns
        # the solutions through nfloat(..., exponent=False).
        rational=flags.get('rational',None)
        if rational is not False:
            from . import Float,nsimplify
            converted=[]
            floats=False
            unknown_set=set(unknowns)
            # SymPy swaps the unknowns that are not symbols (f(x),
            # Derivative(f(x), x)) for Dummy symbols before this test
            applied=[u for u in unknowns if not getattr(u,'is_Symbol',False)]
            for value in scalars:
                if hasattr(value,'has') and value.has(Float):
                    # SymPy looks for Floats after dropping the equations
                    # with free symbols but no unknown that are not constant
                    if (not unknown_set or value.free_symbols&unknown_set
                            or (applied and value.has(*applied))
                            or value.is_number or value.is_constant()):
                        floats=True
                    value=nsimplify(value,rational=True)
                converted.append(value)
            if floats:
                inner={k:v for k,v in flags.items() if k!='rational'}
                target=converted if sequence else converted[0]
                result=solve(target,*symbols,**inner)
                return _restore_floats(result) if rational is None else result
            flags={k:v for k,v in flags.items() if k!='rational'}
        if not isinstance(f,sequences) and len(symbols)==1:
            linear=_linear_solution(prepared[0],symbols[0],flags)
            if linear is not None:
                if not flags.get('simplify',True):
                    return _linear_result(linear[0],linear[1],flags)
                try:
                    solutions=native(f,*symbols,**flags)
                except NotImplementedError:
                    # the native solver does not mask held derivatives
                    # and integrals of the unknown
                    from . import simplify
                    return _linear_result(linear[0],simplify(linear[1]),flags)
                return _simplified_linear(solutions,linear,flags)
        try:
            solutions=native(f,*symbols,**flags)
        except NotImplementedError:
            if isinstance(f,sequences) or len(symbols)!=1:
                raise
            values=_masked_solution(solve,prepared[0],symbols[0],flags)
            if values is None:
                values=_inverted_solution(solve,prepared[0],symbols[0],flags)
            if values is None:
                values=_unrad_solution(solve,prepared[0],symbols[0],flags)
            if values is None:
                raise
            return values
        if (not isinstance(f,sequences) and len(symbols)==1
                and not isinstance(symbols[0],sequences)):
            return _simplified_solutions(solutions,prepared[0],symbols[0],flags)
        if rational is False and isinstance(f,sequences) and len(scalars)>1 and not flags.get('set'):
            return _float_components(solutions,scalars,unknowns)
        return solutions

    namespace['solve']=solve
    # Keep public import paths coherent without replacing the private binding.
    for name,module in tuple(sys.modules.items()):
        if (name.startswith('symox.') and name!='symox._core' and module is not None
                and module.__dict__.get('solve') is native):
            module.solve=solve


def _float_components(result,equations,unknowns):
    """With rational=False, the reference solves each connected component of
    a polynomial system over the domain of its coefficients: the values of a
    component with a Float coefficient come from RR (solve_lin_sys or
    groebner) and are Floats, as nfloat(..., exponent=False) gives them."""
    from . import Float,nfloat,ordered
    if not result or not isinstance(result,(list,dict)):
        return result
    expressions=[]
    for equation in equations:
        if getattr(equation,'is_Equality',False):
            equation=equation.lhs-equation.rhs
        if hasattr(equation,'free_symbols'):
            expressions.append(equation)
    if not unknowns:
        unknowns=list(ordered(set().union(*(e.free_symbols for e in expressions))))
    known=set(unknowns)
    components=[]
    for expression in expressions:
        symbols=expression.free_symbols&known
        if not symbols:
            continue
        joined=[c for c in components if c[0]&symbols]
        for c in joined:
            components.remove(c)
            symbols=symbols|c[0]
        group=[e for c in joined for e in c[1]]+[expression]
        components.append((symbols,group))
    floating=set()
    for symbols,group in components:
        if (any(e.has(Float) for e in group)
                and all(e.is_polynomial(*symbols) for e in group)):
            floating|=symbols
    if not floating:
        return result
    convert=lambda key,value:nfloat(value,exponent=False) if key in floating else value
    def row(item):
        if isinstance(item,dict):
            return {k:convert(k,v) for k,v in item.items()}
        if isinstance(item,tuple) and len(item)==len(unknowns):
            return tuple(convert(k,v) for k,v in zip(unknowns,item))
        return item
    if isinstance(result,dict):
        return row(result)
    return [row(item) for item in result]


def _restore_floats(result):
    """SymPy's nfloat(solution, exponent=False) after a rationalized solve,
    then the removal of the duplicates it reveals and the canonical order."""
    from . import default_sort_key,nfloat
    def convert(item):
        if isinstance(item,dict):
            return {k:nfloat(v,exponent=False) for k,v in item.items()}
        if isinstance(item,tuple):
            return tuple(nfloat(v,exponent=False) for v in item)
        return nfloat(item,exponent=False)
    if isinstance(result,tuple) and len(result)==2 and isinstance(result[1],set):
        return result[0],{convert(row) for row in result[1]}
    if isinstance(result,dict):
        return convert(result)
    if not isinstance(result,list):
        return result
    found=[]
    for item in result:
        item=convert(item)
        if item not in found:
            found.append(item)
    if all(isinstance(item,dict) for item in found) or all(hasattr(item,'free_symbols') for item in found):
        found.sort(key=default_sort_key)
    return found


# force only affects _tsolve and rational only Float arguments, which the
# linear step leaves to the native solver; neither changes solve_linear.
_LINEAR_FLAGS=frozenset(('dict','simplify','check','force','rational','warn'))


def _linear_solution(equation,symbol,flags):
    """SymPy's solve_linear for one equation in one unknown, or None.

    SymPy's _solve returns simplify(solve_linear(...)) whenever the
    numerator is linear in the unknown, before any polynomial solving; the
    unsimplified value is expand_mul(-n(0)/n'(0))."""
    from . import Dummy, Expr, Float, expand_mul, sympify
    if set(flags)-_LINEAR_FLAGS:
        return None
    try:
        equation=sympify(equation)
        symbol=sympify(symbol)
    except Exception:
        return None
    if getattr(equation,'is_Equality',False):
        equation=equation.lhs-equation.rhs
    if not isinstance(equation,Expr) or not isinstance(symbol,Expr):
        return None
    if equation.has(Float) or not equation.has(symbol):
        return None
    if symbol.is_Symbol:
        xi=symbol
        if any(v is not None for v in symbol.assumptions0.values()) and symbol.assumptions0 != {'commutative': True}:
            return None
    elif getattr(symbol,'is_Function',False) or getattr(symbol,'is_Derivative',False):
        # SymPy's recast_to_symbols uses subs, so exp(-x) becomes 1/xi
        # when solving for exp(x).
        xi=Dummy('xi')
        equation=equation.subs(symbol,xi)
        if equation.has(symbol):
            return None
    else:
        return None
    # SymPy masks the derivatives and integrals of the unknown, which gives
    # an implicit solution in terms of them.
    from . import Derivative, Integral
    masked={}
    for node in equation.atoms(Derivative,Integral):
        if node.has(xi):
            masked[node]=Dummy()
    if masked:
        if any(o is not k and o.has(k) for k in masked for o in masked):
            return None
        equation=equation.xreplace(masked)
        if equation.has(Derivative,Integral) and any(
                n.has(xi) for n in equation.atoms(Derivative,Integral)):
            return None
    n,_=equation.as_numer_denom()
    if n.is_Add is not True:
        return None
    derivative=n.diff(xi)
    if not derivative or derivative.has(xi) or derivative.has(Float):
        return None
    if not derivative.free_symbols and not derivative.is_number:
        return None
    value=expand_mul(-1/derivative*n.subs(xi,0))
    if masked:
        value=value.xreplace({v:k for k,v in masked.items()})
    if symbol is not xi:
        value=value.xreplace({xi:symbol})
    return symbol,value


def _simplified_solutions(solutions,equation,symbol,flags):
    """SymPy's _solve simplifies every solution of one equation and drops
    the duplicates this reveals, then solve sorts them canonically.

    The unrad fallback, which SymPy reaches for several radicals in the
    unknown, leaves its solutions unsimplified unless asked."""
    if not flags.get('simplify',True) or flags.get('set') or not isinstance(solutions,list):
        return solutions
    if 'simplify' not in flags and _radical_generators(equation,symbol)>1:
        return solutions
    from . import default_sort_key,simplify
    as_dict=bool(solutions) and all(isinstance(s,dict) for s in solutions)
    if not as_dict and not all(hasattr(s,'free_symbols') for s in solutions):
        return solutions
    kept=set() if 'simplify' in flags else _high_order_roots(equation,symbol,solutions,as_dict)
    if kept is None:
        return solutions
    found=[]
    for solution in solutions:
        try:
            if (solution.get(symbol) if as_dict else solution) in kept:
                value=solution
            else:
                value={k:simplify(v) for k,v in solution.items()} if as_dict else simplify(solution)
        except Exception:
            return solutions
        if value not in found:
            found.append(value)
    if as_dict:
        return sorted(found,key=default_sort_key)
    return sorted(found,key=default_sort_key)


def _high_order_roots(equation,symbol,solutions,as_dict):
    """The solutions SymPy's _solve leaves unsimplified when the simplify
    flag is not given: more than two roots of a polynomial in one generator
    of the unknown ("high-order results are quite long"). None means all of
    them; the roots of each factor of a product are decided separately."""
    from . import sympify
    try:
        equation=sympify(equation)
        if getattr(equation,'is_Equality',False):
            equation=equation.lhs-equation.rhs
        if not equation.is_Mul:
            return None if _polynomial_roots(equation,symbol,len(solutions)) else set()
        kept=set()
        for factor in equation.args:
            if not factor.has(symbol):
                continue
            roots=_polynomial_roots(factor,symbol)
            if roots:
                kept.update(roots)
        if not kept:
            return set()
        values=[s.get(symbol) for s in solutions] if as_dict else list(solutions)
        return None if all(v in kept for v in values) else kept
    except Exception:
        return set()


def _polynomial_roots(f,symbol,found=None):
    """The roots SymPy's _solve gets from the polynomial branch of f when
    there are more than two of them, otherwise None.

    found is the number of solutions already known for f; more than two
    means the branch had more than two roots, so they are not recomputed."""
    from . import Add, Dummy, Poly, roots, solve, solve_linear
    f_num,_=solve_linear(f,symbols=[symbol])
    if f_num.is_zero or f_num.is_Symbol:
        return None
    poly=None
    if not f_num.is_Add:
        add_args=[i for i in f_num.atoms(Add) if symbol in i.free_symbols]
        if len(add_args)==1:
            spart=add_args[0].as_independent(symbol)[1].as_base_exp()[0]
            if spart==symbol:
                try:
                    poly=Poly(f_num,spart)
                except Exception:
                    poly=None
    if poly is None:
        poly=Poly(f_num)
    gens=[g for g in poly.gens if g.has(symbol)]
    if len(gens)!=1:
        return None
    poly=Poly(f_num,gens[0],composite=True)
    if poly.gen==symbol and (poly.degree()<=2 or found is not None and found>2):
        return None if poly.degree()<=2 else True
    soln=roots(poly)
    if sum(soln.values())<poly.degree():
        try:
            soln=poly.all_roots()
        except NotImplementedError:
            pass
    soln=list(soln)
    if poly.gen!=symbol:
        # The reference inverts the generator with tsolve=True, so the
        # inner solve skips the polynomial branch; without the guard every
        # inversion would invert its own generator again.
        if getattr(_INVERTING,'active',False):
            return None
        t=Dummy('t')
        _INVERTING.active=True
        try:
            inverse=solve(poly.gen-t,symbol)
        finally:
            _INVERTING.active=False
        soln={i.subs(t,s) for i in inverse for s in soln}
    return set(soln) if len(soln)>2 else None


_INVERTING=__import__('threading').local()


def _radical_generators(equation,symbol):
    from . import Pow,sympify
    try:
        equation=sympify(equation)
        powers=equation.atoms(Pow)
    except Exception:
        return 0
    return len({p for p in powers if p.has(symbol) and p.exp.is_Rational and not p.exp.is_Integer})


def _linear_result(symbol,value,flags):
    return [{symbol:value}] if flags.get('dict') else [value]


def _simplified_linear(solutions,linear,flags):
    from . import simplify
    symbol,value=linear
    found=[s.get(symbol) for s in solutions] if flags.get('dict') else list(solutions)
    if len(found)!=1 or found[0] is None:
        return solutions
    return _linear_result(symbol,simplify(value),flags)


def _user_inverse(node):
    """inverse() of a one-argument user Function subclass (a host node).

    SymPy's _invert rewrites f(u) = g as u = f.inverse()(g) for any such
    function; symox's native solver only knows its built-in functions."""
    from ._host import _HostFunction
    if not isinstance(node,_HostFunction) or len(getattr(node,'args',()))!=1:
        return None
    inverse=getattr(node,'inverse',None)
    if inverse is None or not callable(inverse):
        return None
    try:
        value=inverse()
    except Exception:
        return None
    return value if callable(value) else None


def _inverse_nodes(equation,symbol):
    from . import Function
    nodes=[n for n in equation.atoms(Function) if n.has(symbol) and _user_inverse(n) is not None]
    return [n for n in nodes if not any(o is not n and o.has(n) for o in nodes)]


def _inverted_solution(solve,equation,symbol,flags):
    """SymPy's factor, generator and _invert steps for user inverse() functions.

    A product is solved factor by factor; otherwise the single invertible
    generator g(u) is solved for as a symbol t and u = g.inverse()(t) is
    solved for the unknown. Returns None when this does not apply."""
    from . import Expr, sympify
    try:
        equation=sympify(equation)
        symbol=sympify(symbol)
    except Exception:
        return None
    if getattr(equation,'is_Equality',False):
        equation=equation.lhs-equation.rhs
    if not isinstance(equation,Expr) or not isinstance(symbol,Expr):
        return None
    if not _inverse_nodes(equation,symbol):
        return None
    inner={k:v for k,v in flags.items() if k!='dict'}
    values=_inverted_values(solve,equation,symbol,inner)
    if values is None:
        return None
    # SymPy's solve orders its solutions canonically
    from . import default_sort_key
    values=sorted(values,key=default_sort_key)
    return [{symbol:v} for v in values] if flags.get('dict') else values


def _inverted_values(solve,equation,symbol,flags):
    from . import Dummy
    values=[]
    def extend(found):
        for value in found:
            if value not in values:
                values.append(value)
    if equation.is_Mul:
        for factor in equation.args:
            if factor.has(symbol):
                extend(solve(factor,symbol,**flags))
        return values
    nodes=_inverse_nodes(equation,symbol)
    t=Dummy('t')
    reduced=equation.xreplace({nodes[0]:t}) if len(nodes)==1 else None
    if reduced is None or reduced.has(symbol):
        # SymPy's _tsolve then factors a sum and solves the factors.
        if equation.is_Add:
            from . import factor,powdenest
            try:
                factored=factor(powdenest(equation))
            except Exception:
                return None
            # a Float content factors out as 1.0*equation: no progress
            if factored.is_Mul and equation not in factored.args:
                return _inverted_values(solve,factored,symbol,flags)
        return None
    node=nodes[0]
    inverse=_user_inverse(node)
    argument=node.args[0]
    for root in solve(reduced,t,**flags):
        target=inverse(root)
        extend([target] if argument==symbol else solve(argument-target,symbol,**flags))
    return values


def _unrad_solution(solve,equation,symbol,flags):
    """SymPy's unrad fallback of _solve for one equation in one symbol.

    The radicals are eliminated from the numerator, the radical-free
    equation (or its change of variable) is solved, and the candidates are
    checked: a candidate is rejected when it makes a denominator zero or the
    equation provably nonzero, and kept when the check is undecidable.
    Returns None when this does not apply."""
    from . import Expr, Symbol, checksol, default_sort_key, simplify, sympify, together
    from .solvers._algebraic import denoms, unrad
    if set(flags)-{'dict','check','simplify','force','rational','warn'}:
        return None
    try:
        equation=sympify(equation)
        symbol=sympify(symbol)
    except Exception:
        return None
    if getattr(equation,'is_Equality',False):
        equation=equation.lhs-equation.rhs
    if not isinstance(equation,Expr) or not isinstance(symbol,Symbol):
        return None
    try:
        numerator=together(equation).as_numer_denom()[0]
        reduced=unrad(numerator,symbol)
    except (ValueError,NotImplementedError,TypeError):
        return None
    if not reduced:
        return None
    inner={k:v for k,v in flags.items() if k!='dict'}
    polynomial,cov=reduced
    # The resultant elimination can return a power of SymPy's radical-free
    # equation; its square-free part has the same roots.
    from . import sqf_part
    try:
        polynomial=sqf_part(polynomial)
    except Exception:
        pass
    try:
        if cov:
            coordinate,relation=cov
            inverse=solve(relation,symbol,**inner)
            if not inverse:
                return None
            values=[inverse[0].subs(coordinate,root) for root in solve(polynomial,coordinate,**inner)]
        else:
            values=solve(polynomial,symbol,**inner)
    except NotImplementedError:
        return None
    found=[]
    for value in values:
        if value not in found:
            found.append(value)
    if flags.get('simplify',False) and 'simplify' in flags:
        simplified=[]
        for value in found:
            value=simplify(value)
            if value not in simplified:
                simplified.append(value)
        found=simplified
    if flags.get('check',True):
        options={k:v for k,v in flags.items() if k in ('force','rational')}
        try:
            dens=denoms(equation,symbol)
            found=[v for v in found
                   if not any(checksol(d,{symbol:v},**options) is True for d in dens)
                   and checksol(equation,{symbol:v},**options) is not False]
        except NotImplementedError:
            pass
    found=sorted(found,key=default_sort_key)
    return [{symbol:v} for v in found] if flags.get('dict') else found


def _masked_solution(solve,equation,symbol,flags):
    """SymPy's solve masks the derivatives, integrals, sums, products and
    substitutions that contain the unknown and solves for the rest, which
    gives an implicit solution in terms of them, or no solution when the
    unknown only appears inside them. Returns None when nothing is masked."""
    from . import Derivative, Dummy, Expr, Integral, Product, Subs, Sum, sympify
    if flags.get('set'):
        return None
    try:
        equation=sympify(equation)
        symbol=sympify(symbol)
    except Exception:
        return None
    if getattr(equation,'is_Equality',False):
        equation=equation.lhs-equation.rhs
    if not isinstance(equation,Expr) or not isinstance(symbol,Expr):
        return None
    if symbol.is_Symbol:
        target=symbol
    else:
        target=Dummy('X0')
        equation=equation.subs(symbol,target)
        if equation.has(symbol):
            return None
    found={}
    def visit(node):
        if not isinstance(node,Expr):
            return
        if isinstance(node,(Derivative,Integral,Subs,Sum,Product)):
            if node.has(target) and node not in found:
                found[node]=Dummy()
            return
        for arg in node.args:
            visit(arg)
    visit(equation)
    if not found:
        return None
    reduced=equation.xreplace(found)
    if reduced.has(target):
        inner={k:v for k,v in flags.items() if k!='dict'}
        values=solve(reduced,target,**inner)
    else:
        values=[]
    restore={v:k for k,v in found.items()}
    values=[v.xreplace(restore) for v in values]
    if target is not symbol:
        values=[v.xreplace({target:symbol}) for v in values]
    return [{symbol:v} for v in values] if flags.get('dict') else values
