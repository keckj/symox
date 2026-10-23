"""Optional Aesara computation graphs built from native expressions.

Aesara is imported only when a graph is requested. Its installation must be
compatible with the active NumPy version.
"""
from functools import reduce
import math
from ._legacy import LegacyPrinter as Printer
from .._kind import matrix_tag
_global_cache={}


def dim_handling(inputs,dim=None,dims=None,broadcastables=None):
    if dim is not None:return {arg:(False,)*int(dim) for arg in inputs}
    if dims is not None:
        maximum=max(dims.values(),default=0)
        return {arg:(False,)*rank+(True,)*(maximum-rank) for arg,rank in dims.items()}
    return dict(broadcastables) if broadcastables is not None else {}


class AesaraPrinter(Printer):
    def __init__(self,cache=None,settings=None,**kwargs):super().__init__(settings,**kwargs);self.cache=_global_cache if cache is None else cache
    def doprint(self,expr,dtypes=None,broadcastables=None):
        import aesara.tensor as tensor
        self.tensor=tensor;self.dtypes=dtypes or {};self.broadcastables=broadcastables or {}
        return self._print(expr)
    def _symbol(self,expr):
        dtype=self.dtypes.get(expr,'float64');broadcast=self.broadcastables.get(expr,(False,False) if hasattr(expr,'shape') else ())
        name=getattr(expr,'name',str(expr));key=(name,str(dtype),tuple(broadcast))
        if key not in self.cache:self.cache[key]=self.tensor.TensorType(dtype=dtype,broadcastable=tuple(broadcast))(name)
        return self.cache[key]
    def emptyPrinter(self,expr):
        import symox as S
        t=self.tensor
        if isinstance(expr,(list,tuple)):return [self._print(value) for value in expr]
        if isinstance(expr,S.Matrix):return t.stack([t.stack([self._print(value) for value in row]) for row in expr.tolist()])
        kind=matrix_tag(expr)
        if kind=='MatrixSymbol' or getattr(expr,'is_Symbol',False):return self._symbol(expr)
        if not isinstance(kind,str):kind=getattr(expr,'class_name',type(expr).__name__)
        if kind=='MatMul':return reduce(t.dot,[self._print(value) for value in expr.args])
        if kind=='Inverse':return t.nlinalg.matrix_inverse(self._print(expr.args[0]))
        if kind in ('Transpose','Adjoint'):
            value=self._print(expr.args[0]).T
            return t.conj(value) if kind=='Adjoint' else value
        if kind in ('Trace','Determinant'):return (t.nlinalg.trace if kind=='Trace' else t.nlinalg.det)(self._print(expr.args[0]))
        expr=S.sympify(expr)
        if expr.is_Integer:return t.constant(int(expr))
        if expr.is_Rational or expr.is_Float:return t.constant(float(expr))
        constants={'Pi':math.pi,'Exp1':math.e,'Infinity':math.inf,'NegativeInfinity':-math.inf,'NaN':math.nan,'ImaginaryUnit':1j,'BooleanTrue':True,'BooleanFalse':False}
        if kind in constants:return t.constant(constants[kind])
        if kind=='Piecewise':
            from .mathml import _pairs
            pairs=_pairs(expr);result=t.constant(float('nan'))
            for value,condition in reversed(pairs):result=t.switch(self._print(condition),self._print(value),result)
            return result
        names={'Add':'add','MatAdd':'add','Mul':'mul','HadamardProduct':'mul','Pow':'pow','Abs':'abs','ceiling':'ceil',
         'asin':'arcsin','acos':'arccos','atan':'arctan','atan2':'arctan2','asinh':'arcsinh','acosh':'arccosh','atanh':'arctanh',
         're':'real','im':'imag','arg':'angle','sign':'sgn','conjugate':'conj','loggamma':'gammaln','Max':'maximum','Min':'minimum',
         'Equality':'eq','Unequality':'neq','StrictLessThan':'lt','LessThan':'le','StrictGreaterThan':'gt','GreaterThan':'ge',
         'And':'bitwise_and','Or':'bitwise_or','Xor':'bitwise_xor','Not':'invert'}
        operation=getattr(t,names.get(kind,kind),None)
        if operation is None:raise NotImplementedError('Aesara graph conversion does not support '+str(kind))
        arguments=[self._print(arg) for arg in expr.args]
        if kind in ('Add','MatAdd','Mul','HadamardProduct','Max','Min','And','Or','Xor'):return reduce(operation,arguments)
        return operation(*arguments)


def _require_aesara(what):
    # SymPy imports aesara with import_module and reports its absence (or an
    # import failure) with its own ImportError text.
    try:
        import aesara
    except ImportError:
        raise ImportError(what) from None
    return aesara


def aesara_code(expr,cache=None,**kwargs):
    _require_aesara("aesara is required for aesara_code")
    return AesaraPrinter(cache=cache).doprint(expr,**kwargs)
def aesara_function(inputs,outputs,scalar=False,*,cache=None,dtypes=None,dim=None,dims=None,broadcastables=None,**kwargs):
    aesara=_require_aesara("Aesara is required for aesara_function")
    printer=AesaraPrinter(cache=cache);broadcast=dim_handling(inputs,dim,dims,broadcastables)
    graph_inputs=[printer.doprint(value,dtypes,broadcast) for value in inputs]
    graph_outputs=printer.doprint(outputs,dtypes,broadcast)
    function=aesara.function(graph_inputs,graph_outputs,**kwargs)
    if not scalar:return function
    def wrapped(*args):
        result=function(*args)
        def convert(value):return value.item() if getattr(value,'ndim',None)==0 else value
        return [convert(value) for value in result] if isinstance(result,list) else convert(result)
    return wrapped

__all__=['AesaraPrinter','aesara_code','aesara_function','dim_handling']
