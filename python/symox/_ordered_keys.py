"""Recursive presentation keys for explicit monomial orders."""
from contextvars import ContextVar
from . import _core

_native_scalar_key=_core.Expr.sort_key
_native_matrix_key=_core._MatrixExpr.sort_key
_active=ContextVar('symox_ordered_keys',default=None)


class _Context:
    def __init__(self):
        self.cache={}
        self.depth=0

    def key(self,value,order):
        from . import S, Set, MatrixExpr, Function, default_sort_key
        from ._matrixexpr_methods import _children, _TupleView, _ExprCondPairView, _NameAtom
        if isinstance(value,Set):value=value._e
        identifier=(type(value),getattr(value,'_dag_id',id(value)),
                    order if isinstance(order,str) else id(order))
        if identifier in self.cache:return self.cache[identifier][1]
        if len(self.cache)>=100000:raise NotImplementedError('Ordered keys exceed 100000 distinct nodes')
        if self.depth>=128:raise NotImplementedError('Ordered keys exceed 128 nested levels')
        self.depth+=1
        try:
            one=_native_scalar_key(S.One)
            if isinstance(value,_NameAtom):
                result=(2,0,'Str'),(1,(str(value),)),one,S.One
            elif isinstance(value,_TupleView):
                name='ExprCondPair' if isinstance(value,_ExprCondPairView) else 'Tuple'
                result=(5,0,name),(len(value),tuple(self.key(arg,order) for arg in value)),one,S.One
            elif getattr(value,'_symox_blockgrid',False):
                arguments=(value.rows,value.cols,_TupleView(tuple(value)))
                result=(5,0,'ImmutableDenseMatrix'),(3,tuple(self.key(arg,order) for arg in arguments)),one,S.One
            elif isinstance(value,MatrixExpr):
                arguments=value.as_ordered_terms(order=order) if value.is_MatAdd else _children(value)
                keys=tuple(self.key(arg,order) for arg in arguments)
                if value.kind=='MatrixConjugate':
                    keys=(((5,0,'Transpose'),(1,keys),one,S.One),)
                result=value.class_key(),(len(keys),keys),one,S.One
            elif isinstance(value,_core.Expr):
                if value.is_Atom:
                    result=_native_scalar_key(value)
                else:
                    coefficient,expression=value.as_coeff_Mul()
                    exponent=S.One
                    if expression.is_Pow:
                        expression,exponent=expression.as_base_exp()
                        if expression==S.Exp1:
                            expression=Function('exp')(exponent)
                            exponent=S.One
                    if expression.is_Dummy:
                        arguments=(_native_scalar_key(expression),)
                    elif expression.is_Atom:
                        arguments=(str(expression),)
                    else:
                        if expression.is_Add:children=expression.as_ordered_terms(order=order)
                        elif expression.is_Mul:children=expression.as_ordered_factors(order=order)
                        else:children=_children(expression)
                        arguments=tuple(self.key(child,order) for child in children)
                    result=expression.class_key(),(len(arguments),arguments),self.key(exponent,order),coefficient
            elif isinstance(value,(tuple,list)):
                # Constructor tuple metadata represents a symbolic Tuple.
                result=self.key(_TupleView(value),order)
            else:
                result=default_sort_key(value,order=order)
            self.cache[identifier]=(value,result)
            return result
        finally:self.depth-=1


def _key(value,order):
    context=_active.get()
    if context is not None:return context.key(value,order)
    context=_Context()
    token=_active.set(context)
    try:return context.key(value,order)
    finally:_active.reset(token)


def sort_key(self,order=None):
    return _native_scalar_key(self) if order is None else _key(self,order)


def matrix_sort_key(self,order=None):
    return _native_matrix_key(self) if order is None else _key(self,order)


def install():
    _core.Expr.sort_key=sort_key
    _core._MatrixExpr.sort_key=matrix_sort_key
