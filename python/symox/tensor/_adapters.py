"""Dispatch native scalar operations when an indexed tensor is encountered."""
from functools import wraps


def install(namespace):
    if namespace.get('_indexed_tensor_adapters_installed',False):return
    namespace['_indexed_tensor_adapters_installed']=True
    from ._expressions import TensExpr,TensAdd,TensMul
    Expr=namespace['Expr']
    def arithmetic(original,name):
        @wraps(original)
        def operation(self,other):
            if isinstance(other,TensExpr):
                if name=='__add__':return TensAdd(self,other)
                if name=='__sub__':return TensAdd(self,-other)
                if name=='__mul__':return TensMul(self,other)
                if name=='__truediv__':raise ValueError('Division by an indexed tensor is undefined')
            return original(self,other)
        return operation
    for name in ('__add__','__sub__','__mul__','__truediv__'):
        setattr(Expr,name,arithmetic(getattr(Expr,name),name))
    def transform(original,name):
        @wraps(original)
        def operation(expression,*args,**kwargs):
            if isinstance(expression,TensExpr):
                if args:raise TypeError('Tensor transformations accept keyword options only')
                return expression.expand(**kwargs) if name=='expand' else expression.simplify(**kwargs)
            return original(expression,*args,**kwargs)
        return operation
    for name in ('simplify','expand'):
        namespace[name]=transform(namespace[name],name)
