"""Public simplify controls around the native candidate search.

These controls do not assert that the native search implements every SymPy
simplification strategy. Unreviewed strategy/domain gaps remain in the audit.
"""
from . import _core
from ._compat import count_ops


def simplify(expr,ratio=1.7,measure=count_ops,rational=False,inverse=False,doit=True,**kwargs):
    from . import sympify, Float, nfloat, nsimplify, S
    from ._core import Expr
    options=dict(ratio=ratio,measure=measure,rational=rational,inverse=inverse,doit=doit)
    if getattr(expr,'_symox_sparse',False):
        return expr.applyfunc(lambda entry:simplify(entry,**options))
    expression=sympify(expr,rational=rational)
    if not isinstance(expression,Expr):
        if getattr(expression,'_symox_matrixexpr',False):
            if any((ratio!=1.7,measure is not count_ops,rational is not False,inverse,not doit)):
                raise NotImplementedError('custom simplify controls on abstract matrix expressions require scalar mapping')
            return _core.simplify(expression)
        if hasattr(expression,'applyfunc'):
            return expression.applyfunc(lambda entry:simplify(entry,**options))
        method=getattr(expression,'simplify',None)
        if method is not None:return method(**options)
        return _core.simplify(expression)
    # SymPy's simplify returns 0 at once for an expression known to be zero;
    # numeric values decide it exactly (algebraic root objects included).
    if expression.args and expression.is_number and expression.is_zero:
        return S.Zero if not expression.is_Number else expression
    if doit and _core._host_contains(expression)[0]:
        # The reference simplify ends with doit() on unhandled heads; host
        # classes that opt in (held WignerD) are evaluated first here.
        held=expression.find(lambda node:getattr(type(node),'_symox_simplify_doit',False))
        if held:
            expression=expression.xreplace({node:node.doit() for node in held})
            if not isinstance(expression,Expr):return simplify(expression,**options)
    from ._host import simplify_hook
    hook=simplify_hook(expression)
    if hook is not None:
        # SymPy's simplify returns the class's own simplification at once,
        # after the zero test of every Expr.
        if expression.is_zero:return S.Zero if not expression.is_Number else expression
        return hook(**options)
    if not expression.args:return expression
    if getattr(expression.func,'__name__',None) in ('Trace','Determinant'):
        # SymPy simplifies the Expr arguments, then returns done(expr), that
        # is expr.doit() for an unhandled head.
        argument=expression.args[0]
        if getattr(argument,'_symox_matrixexpr',False) and doit and ratio==1.7 and measure is count_ops and not rational and not inverse:
            expression=expression.func(simplify(argument))
        return expression.doit() if doit else expression
    if getattr(expression,'is_Relational',False) and expression.rel_op in ('<','<=','>','>='):
        return _relational_simplify(expression,options)
    from .logic.boolalg import BooleanFunction, Xor
    if isinstance(expression,BooleanFunction) and not isinstance(expression,Xor):
        from . import And as _And, Or as _Or
        from .logic import _pattern_simplify
        lattice='And' if expression.func==_And else ('Or' if expression.func==_Or else None)
        if lattice=='Or':
            # Or._eval_simplify: a <= b becomes Eq(a, b) | (a < b) first
            split=_pattern_simplify.or_split_lege(expression)
            if split is not None:return simplify(split,**options)
        # BooleanFunction._eval_simplify starts with simplify_univariate
        from .logic.boolalg import simplify_univariate
        rv=simplify_univariate(expression)
        if not isinstance(rv,BooleanFunction):
            return rv.simplify(**options) if rv.args else rv
        from . import simplify_logic
        rv=rv.func(*[a.simplify(**options) if a.args else a for a in rv.args])
        expression=simplify_logic(rv)
        if not isinstance(expression,BooleanFunction):return expression
        if lattice=='Or':
            if expression.func!=_Or:return expression
            return _pattern_simplify.or_simplify(expression,measure)
        if lattice=='And':
            if expression.func!=_And:return expression
            return _pattern_simplify.and_simplify(expression,measure,ratio)
    if getattr(expression,'is_Piecewise',False):
        # Piecewise._eval_simplify
        if expression.is_zero:return S.Zero
        from ._public_calculus import piecewise_simplify
        return piecewise_simplify(expression,**options)
    if _core._has_noncommutative_atoms(expression):
        from ._noncommutative import simplify as nc_simplify
        return nc_simplify(expression,**options)
    from . import besselj,bessely,besseli,besselk,hankel1,hankel2,jn,yn
    if expression.func in (besselj,bessely,besseli,besselk,hankel1,hankel2,jn,yn):
        # BesselBase._eval_simplify is besselsimp, without the ratio test.
        if expression.is_zero:return S.Zero
        from . import besselsimp
        return besselsimp(expression)
    original=expression
    if inverse:
        from .simplify.simplify import inversecombine
        expression=inversecombine(expression)
    floats=rational is not False and bool(expression.atoms(Float))
    if floats:expression=nsimplify(expression,rational=True)
    result=_core.simplify(expression)
    if measure is count_ops and result.is_Integer and result==0 and isinstance(ratio,(int,float)) and ratio>=0:
        # Zero has no operations: the default ratio test cannot reject it.
        # Avoid traversing the original expression after the native proof.
        return result
    # A custom objective must select actual candidates, rather than merely be
    # accepted and discarded by a **kwargs facade. Keep the native fast search
    # as the ordinary default; broader strategy equivalence remains unaudited.
    if measure is not count_ops:
        from . import expand,cancel,together,powsimp,factor_terms,expand_log,logcombine
        candidates=[result,expression]
        for seed in (expression,result):
            for transform in (expand,cancel,together,powsimp,factor_terms,expand_log,logcombine):
                try:candidate=transform(seed)
                except (NotImplementedError,ValueError):continue
                if candidate not in candidates:candidates.append(candidate)
        result=min(candidates,key=measure)
    # SymPy's simplify cancels after expand_log(deep=True), which brings
    # log(100) and 2*log(10) to one generator. The native search has no such
    # step; it is used here only to recognize zero, so other forms are kept.
    from . import log
    if result.has(log):
        from . import cancel,expand_log
        try:
            if cancel(expand_log(result,deep=True))==0:return S.Zero
        except (NotImplementedError,ValueError):
            pass
    # SymPy's simplify rationalizes a denominator that is a sum of numeric
    # radicals (max_terms=1, symbolic=False) before the final ratio test.
    numer,denom=result.as_numer_denom()
    if denom.is_Add:
        from . import radsimp,fraction
        n,d=fraction(radsimp(1/denom,symbolic=False,max_terms=1))
        if n!=1:
            result=(numer*n).expand()/d
    from . import Piecewise
    if result.has(Piecewise):
        # the Piecewise block of the reference simplify: fold, simplify the
        # pieces and their conditions, then try factor_terms
        from . import factor_terms
        from ._public_calculus import piecewise_fold,piecewise_simplify
        folded=piecewise_fold(result)
        if folded.is_Piecewise:
            folded=piecewise_simplify(folded,**dict(options,doit=False))
            if folded.has(Piecewise):
                factored=factor_terms(folded)
                if measure(factored)<measure(folded):folded=factored
            return folded
    if not (measure is count_ops and result==original):
        if measure(result)>ratio*measure(original):result=original
    if floats and rational is None:result=nfloat(result,exponent=False)
    if doit:
        from . import Integral,Sum,Product,Derivative,Limit,Subs
        if result.has(Integral,Sum,Product,Derivative,Limit,Subs):result=result.doit()
    return result


def _relational_simplify(self,options):
    """``Relational._eval_simplify`` of the reference for the inequalities:
    evaluate a numeric difference, then write a univariate or linear
    relation with the symbols on the left and keep it when it is smaller."""
    from . import S, Add, Poly, poly, gcd, Expr
    from .core.sorting import ordered
    from .solvers.solveset import linear_coeffs
    from .polys.polyerrors import PolynomialError
    r = self.func(*[i.simplify(**options) for i in self.args])
    if r.is_Relational:
        if not isinstance(r.lhs, Expr) or not isinstance(r.rhs, Expr):
            return r
        dif = r.lhs - r.rhs
        v = None
        if dif.is_comparable:
            v = dif.n(2)
            if any(getattr(i, '_prec', None) == 1 for i in v.as_real_imag()):
                rv, iv = [i.n(2) for i in dif.as_real_imag()]
                v = rv + S.ImaginaryUnit*iv
        elif dif.equals(0):
            v = S.Zero
        if v is not None:
            r = r.func(v, S.Zero)
        if r.is_Relational:
            r = r.canonical
        if r.is_Relational:
            free = list(filter(lambda x: x.is_real is not False, r.free_symbols))
            if len(free) == 1:
                try:
                    x = free.pop()
                    dif = r.lhs - r.rhs
                    m, b = linear_coeffs(dif, x)
                    if m.is_zero is False:
                        if m.is_negative:
                            r = r.func(-b / m, x)
                        else:
                            r = r.func(x, -b / m)
                    else:
                        r = r.func(b, S.Zero)
                except ValueError:
                    try:
                        p = poly(dif, x)
                        c = p.all_coeffs()
                        constant = c[-1]
                        c[-1] = 0
                        scale = gcd(c)
                        c = [ctmp / scale for ctmp in c]
                        r = r.func(Poly.from_list(c, x).as_expr(), -constant / scale)
                    except PolynomialError:
                        pass
            elif len(free) >= 2:
                try:
                    free = list(ordered(free))
                    dif = r.lhs - r.rhs
                    m = linear_coeffs(dif, *free)
                    constant = m[-1]
                    del m[-1]
                    scale = gcd(m)
                    m = [mtmp / scale for mtmp in m]
                    nzm = list(filter(lambda f: f[0] != 0, list(zip(m, free))))
                    if scale.is_zero is False:
                        if constant != 0:
                            newexpr = Add(*[i * j for i, j in nzm])
                            r = r.func(newexpr, -constant / scale)
                        else:
                            lhsterm = nzm[0][0] * nzm[0][1]
                            del nzm[0]
                            newexpr = Add(*[i * j for i, j in nzm])
                            r = r.func(lhsterm, -newexpr)
                    else:
                        r = r.func(constant, S.Zero)
                except ValueError:
                    pass
    if r.is_Relational:
        r = r.canonical
    measure = options['measure']
    if measure(r) < options['ratio'] * measure(self):
        return r
    return self


def expression_simplify(self,**kwargs):
    return simplify(self,**kwargs)


def install(namespace):
    namespace['simplify']=simplify
    _core.Expr.simplify=expression_simplify
