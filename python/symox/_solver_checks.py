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


def _checksol_stock(f, symbol, sol=None, **flags):
    """Check equations by substitution, returning True, False or None.

    With numerical=True, scalar numerical residuals use an absolute 1e-9
    tolerance. This is approximate validation rather than an identity proof.
    """
    from . import S, Abs, simplify, expand, DenseMatrix as Matrix
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


_CHECKSOL_FLAGS=('numerical','minimal','warn','simplify','force')


def _checksol_item(value):
    """The native form of one equation, or None when it is not admitted
    (nested iterables, non-sympifiable or non-native objects)."""
    from .core.sympify import _sympify
    from .utilities.iterables import iterable
    from . import Poly
    if iterable(value):
        return None
    try:
        value=_sympify(value)
    except Exception:
        return None
    if isinstance(value,Poly):
        value=value.as_expr()
    return value if isinstance(value,_core.Expr) else None


def checksol(f, symbol, sol=None, **flags):
    """
    Checks whether sol is a solution of equation f == 0.

    Input can be either a single symbol and corresponding value
    or a dictionary of symbols and values. When given as a dictionary
    and flag ``simplify=True``, the values in the dictionary will be
    simplified. *f* can be a single equation or an iterable of equations.
    A solution must satisfy all equations in *f* to be considered valid;
    if a solution does not satisfy any equation, False is returned; if one or
    more checks are inconclusive (and none are False) then None is returned.

    flags: ``numerical=True``, ``minimal=False``, ``warn=False``,
    ``simplify=True``, ``force=True`` (SymPy 1.14's reading; other flags
    are ignored).

    The computation is the native port of SymPy 1.14's ``checksol``
    (``symox::solve::checksol::checksol_system``); arguments it does not
    admit go to the previous implementation.
    """
    from .core.logic import fuzzy_and
    from .utilities.iterables import iterable
    if sol is not None:
        sol={symbol: sol}
    elif isinstance(symbol,dict):
        sol=symbol
    else:
        msg='Expecting (sym, val) or ({sym: val}, None) but got (%s, %s)'
        raise ValueError(msg % (symbol, sol))
    if iterable(f):
        if not f:
            raise ValueError('no functions to check')
        f=list(f)
        items=[_checksol_item(fi) for fi in f]
        if any(item is None for item in items):
            return fuzzy_and(checksol(fi,sol,**flags) for fi in f)
        result=_checksol_native_call(items,sol,flags)
        if result is _DECLINE:
            return fuzzy_and(checksol(fi,sol,**flags) for fi in f)
        return result
    from .core.sympify import _sympify
    from . import Poly, Dict
    item=_sympify(f)
    if isinstance(item,Poly):
        item=item.as_expr()
    if isinstance(item,Dict):
        return _checksol_dict(item,sol)
    if isinstance(item,_core.Expr):
        result=_checksol_native_call([item],sol,flags)
        if result is not _DECLINE:
            return result
    known={name:bool(flags[name]) for name in _CHECKSOL_FLAGS if name in flags}
    return _checksol_stock(f,sol,**known)


_DECLINE=object()


def _checksol_dict(f, sol):
    """SymPy 1.14's path for a ``Dict`` equation: it is not an expression,
    so after the empty and illegal-value checks the reference fails on the
    missing ``is_constant`` method."""
    from . import S, sympify
    if not f:
        return True
    illegal={S.NaN,S.Infinity,S.NegativeInfinity,S.ComplexInfinity}
    if any(sympify(v).atoms()&illegal for v in sol.values()):
        return False
    if f.subs(sol).atoms()&illegal:
        return False
    raise AttributeError("'Dict' object has no attribute 'is_constant'")


def _checksol_native_call(items, sol, flags):
    from . import sympify
    keys=list(sol)
    try:
        native_keys=[sympify(key) for key in keys]
        native_values=[sympify(sol[key]) for key in keys]
    except Exception:
        return _DECLINE
    if not all(isinstance(v,_core.Expr) for v in native_keys+native_values):
        return _DECLINE
    from . import simplify
    # the reference calls the public simplify (its candidate search and
    # post steps): it is handed to the native check
    code,value,warned,values=_core._checksol_public(items,native_keys,native_values,
        bool(flags.get('numerical',True)),bool(flags.get('force',True)),
        bool(flags.get('simplify',True)),bool(flags.get('minimal',False)),simplify)
    if code=='decline':
        return _DECLINE
    original=[sol[key] for key in keys]
    if values is not None:
        for key,new in zip(keys,values):
            sol[key]=new
    if flags.get('warn',False):
        for snapshot in warned:
            shown=dict(zip(keys,original if snapshot is None else snapshot))
            warnings.warn("\n\tWarning: could not verify solution %s." % shown)
    if code=='attribute':
        raise AttributeError("'%s' object has no attribute 'is_constant'" % type(value).__name__)
    return {'true':True,'false':False,'none':None}[code]


def solve_linear(lhs, rhs=0, symbols=(), exclude=()):
    """SymPy's ``solve_linear``: ``(0, 1)``, ``(0, 0)``, ``(symbol, solution)``
    or ``(n, d)`` for ``lhs - rhs``.

    A port of the reference: a symbol is only isolated when the derivative
    of the numerator with respect to it survives differentiation by its own
    free symbols, so an identity in disguise such as
    ``y*cos(x)**2 + y*sin(x)**2 - y`` is reported as ``(0, 1)``.
    """
    from . import S, Eq, Derivative, Integral, default_sort_key, expand_mul
    from .solvers._algebraic import denoms
    if isinstance(lhs,Eq):
        if rhs:
            raise ValueError('If lhs is an Equality, rhs must be 0 but was %s' % rhs)
        rhs=lhs.rhs
        lhs=lhs.lhs
    dens=None
    eq=lhs-rhs
    n,d=eq.as_numer_denom()
    if not n:
        return S.Zero,S.One
    free=n.free_symbols
    if not symbols:
        symbols=free
    else:
        bad=[s for s in symbols if not s.is_Symbol]
        if bad:
            if len(bad)==1:
                bad=bad[0]
            if len(symbols)==1:
                eg='solve(%s, %s)' % (eq,symbols[0])
            else:
                eg='solve(%s, *%s)' % (eq,list(symbols))
            raise ValueError('solve_linear only handles symbols, not %s. To isolate '
                             'non-symbols use solve, e.g. >>> %s <<<.' % (bad,eg))
        symbols=free.intersection(symbols)
    symbols=set(symbols).difference(exclude)
    if not symbols:
        return S.Zero,S.One
    # derivatives of the symbols of interest are evaluated before analysis
    derivs={}
    for der in n.atoms(Derivative):
        for c in der.free_symbols&symbols:
            derivs.setdefault(c,[]).append(der)
    all_zero=True
    for xi in sorted(symbols,key=default_sort_key):
        newn=n.subs({der:der.doit() for der in derivs.get(xi,())})
        dnewn_dxi=newn.diff(xi)
        # dnewn_dxi can be nonzero if it survives differentation by any
        # of its free symbols
        free=dnewn_dxi.free_symbols
        if dnewn_dxi and (not free or any(dnewn_dxi.diff(s) for s in free) or free==symbols):
            all_zero=False
            if dnewn_dxi is S.NaN:
                break
            if xi not in dnewn_dxi.free_symbols:
                vi=-1/dnewn_dxi*(newn.subs(xi,0))
                if dens is None:
                    dens=set()
                    for den in denoms(eq,symbols):
                        if den.is_Pow and den.exp.is_Number:
                            if den.exp.is_zero:
                                continue
                            den=den.base
                        dens.add(den)
                if not any(checksol(di,{xi:vi},minimal=True) is True for di in dens):
                    irep=[(i,i.doit()) for i in vi.atoms(Integral) if i.function.is_number]
                    return xi,expand_mul(vi.subs(irep))
    if all_zero:
        return S.Zero,S.One
    if n.is_Symbol:
        return S.Zero,S.Zero
    return n,d


def solve_linear_system_LU(matrix, syms):
    from . import DenseMatrix as Matrix
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
