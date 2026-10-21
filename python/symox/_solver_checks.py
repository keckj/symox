"""Solution validation and scalar linear isolation through native arithmetic."""
from collections.abc import Mapping
import warnings

from . import _core


def _scalar(value):
    from . import sympify
    value=sympify(value)
    if not isinstance(value,_core.Expr):
        raise TypeError('Solver checks require native scalar expressions')
    return value


def _equation(value, rhs=0):
    from . import Equality
    value=_scalar(value)
    if isinstance(value,Equality):
        if _scalar(rhs)!=0:
            raise ValueError('An Equality requires rhs=0')
        return value.lhs-value.rhs
    return value-_scalar(rhs)


def _illegal(value):
    from . import S, oo
    return value.has(S.NaN,S.ComplexInfinity,oo,-oo)


def _rational_pair(expression):
    """Flatten integer-power rational structure without cancelling poles."""
    from . import S
    cache={}
    def fraction(value,depth=0):
        if depth>256:
            raise NotImplementedError('Rational isolation exceeds 256 nested arithmetic levels')
        if value in cache:return cache[value]
        if value.is_Add:
            pairs=[fraction(arg,depth+1) for arg in value.args]
            denominator=S.One
            for part in dict.fromkeys(d for _,d in pairs):denominator*=part
            result=sum((n*(denominator/d) for n,d in pairs),S.Zero),denominator
        elif value.is_Mul:
            numerator,denominator=S.One,S.One
            for arg in value.args:
                n,d=fraction(arg,depth+1);numerator*=n;denominator*=d
            result=numerator,denominator
        elif value.is_Pow and value.exp.is_Integer:
            n,d=fraction(value.base,depth+1)
            order=int(value.exp)
            result=(n**order,d**order) if order>=0 else (d**(-order),n**(-order))
        else:
            result=value.as_numer_denom()
        cache[value]=result
        return result
    return fraction(expression)


def checksol(f, symbol, sol=None, **flags):
    """Check equations by substitution, returning True, False or None.

    With numerical=True, scalar numerical residuals use an absolute 1e-9
    tolerance. This is approximate validation rather than an identity proof.
    """
    from . import S, Abs, simplify, expand, Matrix
    unknown=set(flags)-{'numerical','minimal','warn','simplify','force'}
    if unknown:
        raise TypeError('Unknown checksol options: '+', '.join(sorted(unknown)))
    for name,value in flags.items():
        if not isinstance(value,bool):
            raise TypeError('%s must be a boolean' % name)
    if sol is None:
        if not isinstance(symbol,Mapping):
            raise ValueError('Supply a symbol and value, or a substitution mapping')
        substitutions={_scalar(key):_scalar(value) for key,value in symbol.items()}
    else:
        substitutions={_scalar(symbol):_scalar(sol)}
    if any(_illegal(value) for value in substitutions.values()):return False
    if isinstance(f,(list,tuple,set,frozenset,Matrix)) or getattr(f,'_symox_sparse',False):
        undecided=False
        for value in f:
            result=checksol(value,substitutions,**flags)
            if result is False:return False
            undecided|=result is None
        return None if undecided else True
    value=_scalar(f)
    if value==S.true:return True
    if value==S.false:return False
    expression=_equation(value)
    value=expression.subs(substitutions)
    if _illegal(value):return False
    minimal=flags.get('minimal',False)
    numerical=flags.get('numerical',True)
    for attempt in range(1 if minimal else 3):
        if value==0 or value.is_zero is True:return True
        if value.is_Rational:return False
        if numerical and value.is_number:
            try:
                # The reference evaluates val.n(18).n(12, chop=True): the
                # chop resolves an exactly cancelling component.
                numeric=value.n(18).n(12,chop=True)
                if _illegal(numeric):return False
                return bool(Abs(numeric)<S('1e-9'))
            except (TypeError,ValueError,NotImplementedError):
                pass
        if value.is_zero is False:return False
        if attempt==0 and not minimal:
            value=expand(value)
        elif flags.get('simplify',True):
            value=simplify(expression.subs(substitutions))
    generators=value.free_symbols
    if (generators and all(isinstance(gen,_core.Expr) and gen.is_Symbol for gen in generators)
            and value.is_rational_function(*generators) is True):
        from . import cancel
        reduced=cancel(value)
        return reduced==0
    if flags.get('warn',False):
        warnings.warn('Could not verify solution %s' % substitutions,RuntimeWarning,stacklevel=2)
    return None


def solve_linear(lhs, rhs=0, symbols=(), exclude=()):
    """Isolate a linear numerator variable while excluding denominator zeros."""
    from . import S, default_sort_key, expand_mul
    from .calculus._analysis import walk
    expression=_equation(lhs,rhs)
    numerator,denominator=_rational_pair(expression)
    if numerator==0:return S.Zero,S.One
    provided=tuple(_scalar(value) for value in symbols)
    if any(not value.is_Symbol for value in provided):
        raise ValueError('solve_linear requires symbols for isolation')
    variables=set(provided) if provided else {s for s in numerator.free_symbols
                                             if isinstance(s,_core.Expr) and s.is_Symbol}
    variables=(variables&numerator.free_symbols)-{_scalar(value) for value in exclude}
    if not variables:return S.Zero,S.One
    derivatives=[]
    denominators=[]
    for node in walk(expression):
        if node.is_Pow and node.exp.is_negative is True:denominators.append(node.base)
        if node.is_Derivative and node.free_symbols&variables:derivatives.append(node)
    evaluated={node:node.doit() for node in derivatives}
    numerator=numerator.xreplace(evaluated)
    all_zero=True
    for symbol in sorted(variables,key=default_sort_key):
        coefficient=numerator.diff(symbol)
        if coefficient==0:continue
        all_zero=False
        if _illegal(coefficient) or coefficient.has(symbol):continue
        candidate=-numerator.subs(symbol,0)/coefficient
        if _illegal(candidate):continue
        invalid=False
        for base in denominators:
            residual=base.subs(symbol,candidate)
            if residual==0 or residual.is_zero is True or _illegal(residual):
                invalid=True
                break
        if not invalid:return symbol,expand_mul(candidate)
    if all_zero:return S.Zero,S.One
    if numerator.is_Symbol:return S.Zero,S.Zero
    return numerator,denominator


def solve_linear_system_LU(matrix, syms):
    from . import Matrix
    matrix=Matrix(matrix)
    symbols=tuple(_scalar(value) for value in syms)
    if matrix.rows!=matrix.cols-1:
        raise ValueError('Rows should be equal to columns - 1')
    if len(symbols)!=matrix.rows or len(set(symbols))!=len(symbols):
        raise ValueError('Provide one distinct symbol per unknown')
    if any(not value.is_Symbol for value in symbols):
        raise ValueError('Linear system unknowns must be symbols')
    if matrix.rows==0:return {}
    solution=matrix[:,:-1].LUsolve(matrix[:,-1:])
    return dict(zip(symbols,solution))


def solve_undetermined_coeffs(equ, coeffs, *syms, **flags):
    """Match independent polynomial or functional coefficients and solve them."""
    from . import cancel, expand, solve, default_sort_key
    coeffs=tuple(_scalar(value) for value in coeffs)
    if not coeffs or any(not value.is_Symbol for value in coeffs):
        raise ValueError('Provide symbols for the unknown coefficients')
    if len(set(coeffs))!=len(coeffs):
        raise ValueError('Unknown coefficients must be distinct')
    expression=_equation(equ)
    normalized=expand(cancel(expression).as_numer_denom()[0])
    active=tuple(value for value in coeffs if value in normalized.free_symbols)
    as_dict=flags.pop('dict',False)
    as_set=flags.pop('set',False)
    if as_dict and as_set:raise ValueError('Choose dictionary or set output')
    if not active:return ([],set()) if as_set else []
    if not syms:
        independent,dependent=normalized.as_independent(*active,as_Add=True)
        variables=(dependent.free_symbols&independent.free_symbols)-set(active)
        if not variables:variables=dependent.free_symbols-set(active)
        syms=tuple(sorted(variables,key=default_sort_key))
    elif len(syms)==1 and isinstance(syms[0],(list,tuple)):
        syms=tuple(syms[0])
    if syms:
        syms=tuple(_scalar(value) for value in syms)
        coefficients=normalized.as_coefficients_dict(*syms)
        if any(generator.has(*active) for generator in coefficients):return None
        equations=list(coefficients.values())
    else:
        equations=[normalized]
    solutions=solve(equations,list(active),dict=True,**flags)
    # Clearing parameter denominators must not introduce undefined values.
    if flags.get('check',True):
        solutions=[solution for solution in solutions if checksol(expression,solution,
                   simplify=flags.get('simplify',True)) is not False]
    if as_set:
        solved=set().union(*(solution.keys() for solution in solutions)) if solutions else set(active)
        ordered=tuple(sorted(solved,key=default_sort_key))
        return list(ordered),{tuple(solution.get(value,value) for value in ordered) for solution in solutions}
    return solutions if as_dict or len(solutions)!=1 else solutions[0]


__all__=['checksol','solve_linear','solve_linear_system_LU','solve_undetermined_coeffs']
