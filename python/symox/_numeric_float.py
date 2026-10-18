"""Controlled floating conversion with protected exponent and shape metadata."""
import operator
from collections import defaultdict
from . import _core


def nfloat(expr, n=15, exponent=False, dkeys=False):
    """Float numeric parts, retaining exact exponents unless requested.

    Undefined-function arguments stay exact. Containers retain their type and
    dictionary keys are preserved by default. Native matrices and arrays are
    transformed entrywise; symbolic matrix dimensions and indices stay exact.
    """
    from . import (Expr, MatrixExpr, AlgebraicNumber, Set, Interval, S, Float,
                  Dummy, Pow, sympify, _wrap_set)
    from .tensor.array import NDimArray
    from ._matrixexpr_methods import _map_scalars
    from ._structural_queries import _graph, _NON_FUNCTIONS
    if isinstance(n,bool):raise TypeError('Floating precision must be an integer')
    digits=operator.index(n)
    if digits<1:raise ValueError('Floating precision must be positive')
    if digits>10000:raise NotImplementedError('Floating conversion is limited to 10000 decimal digits')
    cache={}
    active=set()
    def convert(value):
        if isinstance(value,dict):
            identity=id(value)
            if identity in active:raise ValueError('Cyclic floating-conversion container')
            active.add(identity)
            try:
                pairs=[(convert(key) if dkeys else key,convert(item)) for key,item in value.items()]
                if isinstance(value,defaultdict):return type(value)(value.default_factory,pairs)
                return type(value)(pairs)
            finally:active.remove(identity)
        if isinstance(value,(list,tuple,set,frozenset)):
            identity=id(value)
            if identity in active:raise ValueError('Cyclic floating-conversion container')
            active.add(identity)
            try:
                values=[convert(item) for item in value]
                return type(value)(*values) if hasattr(value,'_fields') else type(value)(values)
            finally:active.remove(identity)
        if isinstance(value,NDimArray):return value.applyfunc(convert)
        if isinstance(value,MatrixExpr):return _map_scalars(value,convert,preserve_indices=True)
        if isinstance(value,_core.Matrix) or getattr(value,'is_Matrix',False):return value.applyfunc(convert)
        if isinstance(value,AlgebraicNumber):return value.as_expr().evalf(digits)
        if isinstance(value,Set):
            if not value.args:return value
            if isinstance(value,Interval):
                start=S.Zero if value.start.is_zero is True else convert(value.start)
                end=S.Zero if value.end.is_zero is True else convert(value.end)
                return value.func(start,end,value.left_open,value.right_open)
            return value.func(*(convert(item) for item in value.args))
        if not isinstance(value,Expr):value=sympify(value)
        if not isinstance(value,Expr):raise TypeError('Floating conversion requires native symbolic values')
        wrapped=_wrap_set(value)
        if isinstance(wrapped,Set):return convert(wrapped)
        identifier=value._dag_id
        if identifier in cache:return cache[identifier][1]
        if len(cache)>=100000:raise NotImplementedError('Floating conversion exceeds 100000 distinct nodes')
        result=convert_scalar(value)
        cache[identifier]=(value,result)
        return result

    def convert_scalar(value):
        head=getattr(value.func,'__name__','')
        if value.is_Number:return Float(value,digits)
        if value.is_number is True:
            result=value.evalf(digits)
            return Float(result,digits) if result.is_Number else result
        if value.is_Atom or value.is_UndefinedFunction:return value
        if value.is_Relational:return value.func(*(convert(arg) for arg in value.args))
        if head=='Tuple':return value.func(*(convert(arg) for arg in value.args))
        if head=='Order':return value.func(convert(value.args[0]),*value.args[1:])
        if head=='Lambda':
            body=value.args[1]
            return value.func(value.args[0],S.Zero if body.is_zero is True else convert(body))
        if head=='Subs':
            evaluated=value.doit()
            return value if evaluated==value else convert(evaluated)
        if head=='MatrixElement':return value.func(convert(value.args[0]),*value.args[1:])
        if head in ('Trace','Determinant','DotProduct'):
            return value.func(*(convert(arg) for arg in value.args))
        _,records,order=_graph(value)
        masks={}
        restore={}
        if not exponent:
            for identifier in order:
                node,_=records[identifier]
                if getattr(node,'is_Pow',False):
                    marker=Dummy('float_exponent')
                    base=node.base.xreplace(masks)
                    masks[node]=Pow(base,marker)
                    restore[marker]=node.exp
        transformed=value.xreplace(masks) if masks else value
        result=transformed.evalf(digits)
        if restore:result=result.xreplace(restore)
        if exponent:
            result=result.replace(lambda node:getattr(node,'is_Pow',False) and node.exp.is_Integer,
                                  lambda node:Pow(node.base,Float(node.exp,digits)))
        # evalf deliberately leaves symbolic built-in function arguments
        # exact. nfloat converts these while preserving undefined functions.
        return result.replace(
            lambda node:bool(getattr(node,'is_Function',False)) and not node.is_UndefinedFunction
                and getattr(node.func,'__name__','') not in _NON_FUNCTIONS-{'Min','Max'},
            lambda node:node.func(*(convert(arg) for arg in node.args)))
    return convert(expr)


__all__=['nfloat']
