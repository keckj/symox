"""Dispatch native scalar operations when an indexed tensor is encountered."""
from functools import wraps
from .._dispatch_arguments import FirstArgument


def install(namespace):
    if namespace.get('_indexed_tensor_adapters_installed',False):return
    namespace['_indexed_tensor_adapters_installed']=True
    from ._expressions import TensExpr,TensAdd,TensMul
    from .. import _core
    Expr=_core.Expr
    def arithmetic(original,name):
        def tensor_operation(self,other):
            if name=='__add__':return TensAdd(self,other)
            if name=='__sub__':return TensAdd(self,-other)
            if name=='__mul__':return TensMul(self,other)
            raise ValueError('Division by an indexed tensor is undefined')
        @wraps(original)
        def operation(self,other):
            if isinstance(other,TensExpr):
                return tensor_operation(self,other)
            return original(self,other)
        # The public priority/context wrapper can perform this same predicate
        # and invoke the selected handler without a second Python frame.
        operation._symox_operand_dispatch=(TensExpr,tensor_operation,original)
        return operation
    for name in ('__add__','__sub__','__mul__','__truediv__'):
        setattr(Expr,name,arithmetic(getattr(Expr,name),name))
    def transform(original,name):
        first=FirstArgument(original)
        @wraps(original)
        def operation(*args,**kwargs):
            expression=first.get(args,kwargs)
            if isinstance(expression,TensExpr):
                expression,rest,options=first.split(args,kwargs)
                if rest:raise TypeError('Tensor transformations accept keyword options only')
                return expression.expand(**options) if name=='expand' else expression.simplify(**options)
            return original(*args,**kwargs)
        return operation
    for name in ('simplify','expand'):
        namespace[name]=transform(namespace[name],name)
