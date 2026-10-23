"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'besselsimp hypersimilar hypersimp kroneckersimp logcombine nsimplify posify separatevars signsimp simplify')

import symox as S


def dotprodsimp(expr,withsimp=False):
    from symox import cancel,expand,sympify
    original=sympify(expr);result=cancel(expand(original))
    return (result,result!=original) if withsimp else result


def clear_coefficients(expr,rhs=0):
    expr,rhs=S.sympify(expr),S.sympify(rhs)
    if expr.is_Rational:return S.Integer(0),rhs-expr
    previous=None
    while expr!=previous and expr!=0:
        previous=expr;content,expr=expr.as_content_primitive();rhs/=content
        constant,expr=expr.as_coeff_Add(rational=True);rhs-=constant
    if expr.could_extract_minus_sign():expr,rhs=-expr,-rhs
    return expr,rhs


def factor_sum(self,limits=None,radical=False,clear=False,fraction=False,sign=True):
    expr=S.Sum(self,*limits) if limits is not None else self
    if not expr.is_Sum:return S.factor_terms(expr,radical=radical,clear=clear,fraction=fraction,sign=sign)
    body=S.factor_terms(expr.function,radical=radical,clear=clear,fraction=fraction,sign=sign)
    variables=[limit[0] for limit in expr.limits];coefficient,remainder=body.as_independent(*variables,as_Add=False)
    return coefficient*S.Sum(remainder,*expr.limits)


def _weighted_sum(expr):
    if getattr(expr,'is_Sum',False):return expr
    factors=S.Mul.make_args(expr);sums=[factor for factor in factors if factor.is_Sum]
    if len(sums)!=1:return None
    main=sums[0];coefficient=S.Mul(*(factor for factor in factors if factor!=main))
    if any(coefficient.has(limit[0]) for limit in main.limits):return None
    return S.Sum(coefficient*main.function,*main.limits)


def sum_add(self,other,method=0):
    left,right=_weighted_sum(self),_weighted_sum(other)
    if left is not None and right is not None:
        if method==0 and left.limits==right.limits:return factor_sum(S.Sum(left.function+right.function,*left.limits))
        if method==1 and S.simplify(left.function-right.function)==0 and len(left.limits)==len(right.limits)==1:
            i,a,b=left.limits[0];j,c,d=right.limits[0]
            if i==j:
                if S.simplify(c-b-1)==0:return factor_sum(S.Sum(left.function,(i,a,d)))
                if S.simplify(a-d-1)==0:return factor_sum(S.Sum(left.function,(i,c,b)))
    return self+other


def product_mul(self,other,method=0):
    if getattr(self,'is_Product',False) and getattr(other,'is_Product',False):
        if method==0 and self.limits==other.limits:return S.Product(self.function*other.function,*self.limits)
        if method==1 and self.function==other.function and len(self.limits)==len(other.limits)==1:
            i,a,b=self.limits[0];j,c,d=other.limits[0]
            if i==j:
                if S.simplify(c-b-1)==0:return S.Product(self.function,(i,a,d))
                if S.simplify(a-d-1)==0:return S.Product(self.function,(i,c,b))
    return self*other


def _combine(terms,operation,predicate):
    terms=list(terms)
    for method in (0,1):
        changed=True
        while changed:
            changed=False
            for i in range(len(terms)):
                for j in range(i+1,len(terms)):
                    value=operation(terms[i],terms[j],method)
                    if predicate(value):terms[i]=value;terms.pop(j);changed=True;break
                if changed:break
    return terms

def sum_combine(s_t):return S.Add(*_combine(s_t,sum_add,lambda value:_weighted_sum(value) is not None))
def sum_simplify(s,**kwargs):
    expr=S.sympify(s)
    if not expr.is_Add:return S.bottom_up(expr,lambda node:sum_combine(node.args) if node.is_Add and any(arg.has(S.Sum) for arg in node.args) else node)
    return sum_combine(expr.args)
def product_simplify(s,**kwargs):return S.Mul(*_combine(S.Mul.make_args(S.sympify(s)),product_mul,lambda value:value.is_Product))


def inversecombine(expr):
    """Simplify the composition of a function and its inverse.

    Port of SymPy's ``inversecombine``: no attention is paid to whether the
    inverse is a left inverse or a right inverse, so ``asin(sin(x))`` gives
    ``x`` whatever the branch.
    """
    log,exp=S.log,S.exp
    def f(rv):
        if isinstance(rv,log):
            if isinstance(rv.args[0],exp) or (rv.args[0].is_Pow and rv.args[0].base==S.S.Exp1):
                rv=rv.args[0].exp
        elif rv.is_Function and hasattr(rv,'inverse'):
            if (len(rv.args)==1 and len(rv.args[0].args)==1 and
               isinstance(rv.args[0],rv.inverse(argindex=1))):
                rv=rv.args[0].args[0]
        if rv.is_Pow and rv.base==S.S.Exp1:
            if isinstance(rv.exp,log):
                rv=rv.exp.args[0]
        return rv
    return S.bottom_up(expr,f)


def nc_simplify(expr,deep=True):
    if isinstance(expr,S.Expr) and not getattr(expr,'is_Matrix',False) and not getattr(expr,'_symox_matrixexpr',False):
        from symox._noncommutative import nc_simplify as scalar_nc_simplify
        return scalar_nc_simplify(expr,deep=deep)
    if isinstance(expr,S.MatrixExpr):
        from ._nc_simplify import nc_simplify as matrix_nc_simplify
        return matrix_nc_simplify(expr,deep=deep)
    return S.bottom_up(expr,S.factor_nc) if deep else S.factor_nc(expr)


def _nthroot_solve(p, n, prec):
    """SymPy's helper for ``nthroot``: denest ``p**(1/n)`` through its
    minimal polynomial."""
    from ..polys.numberfields.minpoly import _minimal_polynomial_sq
    from .sqrtdenest import _mexpand
    while n % 2 == 0:
        p = S.sqrtdenest(S.sqrt(p))
        n = n // 2
    if n == 1:
        return p
    pn = p**S.Rational(1, n)
    x = S.Symbol('x')
    f = _minimal_polynomial_sq(p, n, x)
    if f is None:
        return None
    sols = S.solve(f, x)
    for sol in sols:
        if abs(sol - pn).n() < 1./10**prec:
            sol = S.sqrtdenest(sol)
            if _mexpand(sol**n) == p:
                return sol


def nthroot(expr, n, max_len=4, prec=15):
    """SymPy's ``nthroot``: a real nth root of a sum of surds, denested
    with ``nsimplify`` on the surds or through the minimal polynomial."""
    from ..polys.numberfields.minpoly import _is_sum_surds
    from .sqrtdenest import _mexpand
    expr = S.sympify(expr)
    n = S.sympify(n)
    p = expr**S.Rational(1, n)
    if not n.is_integer:
        return p
    if not _is_sum_surds(expr):
        return p
    surds = []
    coeff_muls = [x.as_coeff_Mul() for x in expr.args]
    for x, y in coeff_muls:
        if not x.is_rational:
            return p
        if y is S.S.One:
            continue
        if not (y.is_Pow and y.exp == S.S.Half and y.base.is_integer):
            return p
        surds.append(y)
    surds.sort(key=lambda y: y.base)  # numeric order, as SymPy's sort()
    surds = surds[:max_len]
    if expr < 0 and n % 2 == 1:
        p = (-expr)**S.Rational(1, n)
        a = S.nsimplify(p, constants=surds)
        res = a if _mexpand(a**n) == _mexpand(-expr) else p
        return -res
    a = S.nsimplify(p, constants=surds)
    if _mexpand(a) is not _mexpand(p) and _mexpand(a**n) == _mexpand(expr):
        return _mexpand(a)
    try:
        expr = _nthroot_solve(expr, n, prec)
    except ValueError:
        # the native factorization of the minimal polynomial gave up on
        # its resource limit (SymPy factors it): keep the root
        return p
    if expr is None:
        return p
    return expr

__all__ += ['clear_coefficients','dotprodsimp','factor_sum','inversecombine','nc_simplify','nthroot','product_mul',
 'product_simplify','sum_add','sum_combine','sum_simplify']

# SymPy's simplify module also exposes trigsimp (defined in simplify.trigsimp).
trigsimp = S.trigsimp
