from ..._matrixexpr import KroneckerProduct, kronecker_product
from ..._kind import matrix_tag

__all__ = ['KroneckerProduct', 'kronecker_product']


def combine_kronecker(expr):
    from symox import MatMul,MatAdd
    kind=matrix_tag(expr)
    if kind not in ('MatMul','MatAdd'):return expr
    args=list(expr.args)
    if kind=='MatMul':
        i=0
        while i+1<len(args):
            a,b=args[i:i+2]
            if isinstance(a,KroneckerProduct) and isinstance(b,KroneckerProduct) and len(a.args)==len(b.args) and all(x.cols==y.rows for x,y in zip(a.args,b.args)):
                args[i:i+2]=[KroneckerProduct(*(x*y for x,y in zip(a.args,b.args)))]
            else:i+=1
        return MatMul(*args)
    return MatAdd(*args)
__all__.append('combine_kronecker')
