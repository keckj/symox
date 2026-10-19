"""Separate integer bit operations from symbolic Boolean operations."""
import operator
from numbers import Number as PythonNumber
from decimal import Decimal
from . import _core


def _integer(value):
    return isinstance(value,int) or isinstance(value,_core.Expr) and value.is_Integer


def _numeric(value):
    return (isinstance(value,(PythonNumber,Decimal)) and not isinstance(value,bool)
            or isinstance(value,_core.Expr) and value.is_number)


def binary(name,a,b):
    from . import Integer,And,Or,Xor,Implies,Intersection,Union
    if _integer(a) and _integer(b):
        return Integer(getattr(operator,name)(int(a),int(b)))
    if name in ('and_','or_') and (_core._is_set(a) or _core._is_set(b)):
        return (Intersection if name=='and_' else Union)(a,b)
    # Connective constructors accept integer truth values in some positions;
    # Boolean operators themselves follow the stricter Boolean operand rule.
    if _numeric(a) or _numeric(b):return NotImplemented
    a,b=And(a),And(b)
    if name=='and_':return And(a,b)
    if name=='or_':return Or(a,b)
    if name=='xor':return Xor(a,b)
    if name=='rshift':return Implies(a,b)
    if name=='lshift':return Implies(b,a)
    raise AssertionError(name)


def invert(value):
    from . import Integer,Not,And
    if value.is_Integer:return Integer(~int(value))
    if value.is_number:raise TypeError('bad operand type for unary ~')
    return Not(And(value))


def install():
    for name in ('and_','or_','xor','rshift','lshift'):
        special=name.rstrip('_')
        def direct(self,other,_name=name):return binary(_name,self,other)
        def reflected(self,other,_name=name):return binary(_name,other,self)
        setattr(_core.Expr,'__'+special+'__',direct)
        setattr(_core.Expr,'__r'+special+'__',reflected)
    _core.Expr.__invert__=invert
