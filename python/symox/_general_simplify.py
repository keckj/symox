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
    from ._host import simplify_hook
    hook=simplify_hook(expression)
    if hook is not None:
        # SymPy's simplify returns the class's own simplification at once,
        # after the zero test of every Expr.
        if expression.is_zero:return S.Zero if not expression.is_Number else expression
        return hook(**options)
    if not expression.args:return expression
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
    if not (measure is count_ops and result==original):
        if measure(result)>ratio*measure(original):result=original
    if floats and rational is None:result=nfloat(result,exponent=False)
    if doit:
        from . import Integral,Sum,Product,Derivative,Limit,Subs
        if result.has(Integral,Sum,Product,Derivative,Limit,Subs):result=result.doit()
    return result


def expression_simplify(self,**kwargs):
    return simplify(self,**kwargs)


def install(namespace):
    namespace['simplify']=simplify
    _core.Expr.simplify=expression_simplify
