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
    """Cancel formal inverse compositions; branch restrictions are ignored."""
    pairs={'log':'exp','exp':'log','asin':'sin','sin':'asin','acos':'cos','cos':'acos','atan':'tan','tan':'atan',
     'asinh':'sinh','sinh':'asinh','acosh':'cosh','cosh':'acosh','atanh':'tanh','tanh':'atanh'}
    def visit(node):
        name=getattr(node,'class_name','')
        if name in pairs and len(node.args)==1:
            argument=node.args[0]
            if getattr(argument,'class_name','')==pairs[name] and len(argument.args)==1:return argument.args[0]
        return node
    return S.bottom_up(expr,visit)


def nc_simplify(expr,deep=True):
    if isinstance(expr,S.Expr) and not getattr(expr,'is_Matrix',False) and not getattr(expr,'_symox_matrixexpr',False):
        from symox._noncommutative import nc_simplify as scalar_nc_simplify
        return scalar_nc_simplify(expr,deep=deep)
    return S.bottom_up(expr,S.factor_nc) if deep else S.factor_nc(expr)


def nthroot(expr,n,max_len=4,prec=15):
    expr,n=S.sympify(expr),S.sympify(n)
    if n==0:raise ValueError('a root order must be nonzero')
    root=S.real_root(expr,n) if n.is_Integer and int(n)%2 else expr**(1/n)
    constants=sorted([node for node in S.preorder_traversal(expr) if node.is_Pow and node.exp==S.Rational(1,2) and node.is_number],key=S.default_sort_key)[:max_len]
    if constants:
        try:
            candidate=S.nsimplify(root.evalf(prec),constants=constants)
            if S.simplify(S.expand(candidate**n-expr))==0:return candidate
        except (NotImplementedError,ValueError,TypeError):pass
        if not expr.free_symbols:
            import mpmath
            context=mpmath.mp.clone();context.dps=max(60,int(prec)+20)
            values={'r'+str(i):context.mpf(str(value.evalf(context.dps))) for i,value in enumerate(constants)}
            try:
                identified=context.identify(context.mpf(str(root.evalf(context.dps))),constants=values,tol=context.mpf(10)**(-context.dps+15),maxcoeff=10000)
                if identified:
                    candidate=S.parse_expr(identified,local_dict={'r'+str(i):value for i,value in enumerate(constants)})
                    if S.simplify(S.expand(candidate**n-expr))==0:return candidate
            except (ValueError,TypeError,NotImplementedError):pass
    return S.sqrtdenest(root)

__all__ += ['clear_coefficients','dotprodsimp','factor_sum','inversecombine','nc_simplify','nthroot','product_mul',
 'product_simplify','sum_add','sum_combine','sum_simplify']

# SymPy's simplify module also exposes trigsimp (defined in simplify.trigsimp).
trigsimp = S.trigsimp
