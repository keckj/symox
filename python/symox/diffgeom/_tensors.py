"""Finite coordinate components for vector fields, forms and mixed tensors."""
from dataclasses import dataclass
from itertools import product,permutations
import math
import operator

from .. import Integer,Rational,simplify,sympify,trigsimp,expand_trig
from ._coordinates import (_native,_scalar_in_chart,_scalar_at,_field_charts,
                           BaseScalarField,CoordSystem,Point)


@dataclass(frozen=True)
class _Slot:
    variance: int
    chart: CoordSystem
    index: int

    def __repr__(self):
        name = self.chart._input_symbols[self.index].name
        return ('e_' if self.variance == 1 else 'd')+name


def _tensor(components,variance,cls=None,**metadata):
    combined = {}
    for slots,coefficient in components.items():
        value = trigsimp(expand_trig(simplify(_native(coefficient))))
        if value != 0:
            combined[tuple(slots)] = value
    if not variance:
        return combined.get((),Integer(0))
    if not combined:
        return Integer(0)
    cls = cls or (_VectorField if tuple(variance) == (1,) else _Tensor)
    obj = object.__new__(cls)
    obj._components,obj.variance = combined,tuple(variance)
    for key,value in metadata.items():
        setattr(obj,key,value)
    return obj


class _Tensor:
    def __repr__(self):
        return ' + '.join('(%s)*%s' % (c,' @ '.join(map(str,slots)))
                          for slots,c in sorted(self._components.items(),key=lambda item:repr(item[0])))

    __str__ = __repr__

    @property
    def args(self):
        return tuple((coefficient,slots) for slots,coefficient in self._components.items())

    @property
    def components(self):
        return dict(self._components)

    @property
    def free_symbols(self):
        return set().union(*(c.free_symbols for c in self._components.values()))

    def __eq__(self,other):
        return isinstance(other,_Tensor) and self.variance == other.variance and self._components == other._components

    __hash__ = None

    def __bool__(self):
        return bool(self._components)

    def __add__(self,other):
        if not isinstance(other,_Tensor):
            if _native(other) == 0:
                return self
            raise ValueError('a scalar cannot be added to a positive-rank tensor')
        if self.variance != other.variance:
            raise ValueError('tensor variance must agree for addition')
        components = self.components
        for slots,c in other._components.items():
            components[slots] = components.get(slots,Integer(0))+c
        return _tensor(components,self.variance)

    __radd__ = __add__

    def __neg__(self):
        return _tensor({slots:-c for slots,c in self._components.items()},self.variance)

    def __sub__(self,other):
        return self+-other

    def __rsub__(self,other):
        return -self+other

    def __mul__(self,scalar):
        if isinstance(scalar,_Tensor):
            raise ValueError('use TensorProduct or WedgeProduct for tensor multiplication')
        scalar = _native(scalar)
        return _tensor({slots:c*scalar for slots,c in self._components.items()},self.variance)

    __rmul__ = __mul__

    def __truediv__(self,scalar):
        scalar = _native(scalar)
        if scalar == 0:
            raise ZeroDivisionError('division by zero')
        return self*(1/scalar)

    def _map(self,method,*args,**kwargs):
        return _tensor({slots:getattr(c,method)(*args,**kwargs) for slots,c in self._components.items()},self.variance)

    def subs(self,*args,**kwargs):
        return self._map('subs',*args,**kwargs)

    def xreplace(self,rule):
        return self._map('xreplace',rule)

    def simplify(self,**kwargs):
        return self._map('simplify',**kwargs)

    def expand(self,**kwargs):
        return self._map('expand',**kwargs)

    def diff(self,*args,**kwargs):
        return self._map('diff',*args,**kwargs)

    def doit(self,**kwargs):
        return self

    def __call__(self,*fields):
        if len(fields) > len(self.variance):
            raise ValueError('too many contraction arguments')
        fields = fields+(None,)*(len(self.variance)-len(fields))
        result = Integer(0)
        remaining_variance = tuple(v for v,f in zip(self.variance,fields) if f is None)
        components = {}
        for slots,c in self._components.items():
            remaining = []
            value = c
            for slot,field in zip(slots,fields):
                if field is None:
                    remaining.append(slot)
                elif slot.variance == 1:
                    if isinstance(field,_Tensor):
                        if field.variance != (-1,):
                            raise ValueError('a vector slot contracts a one-form or acts on a scalar')
                        vector = BaseVectorField(slot.chart,slot.index)
                        value *= field(vector)
                    else:
                        value *= _scalar_in_chart(field,slot.chart).diff(_native(slot.chart.base_scalar(slot.index)))
                else:
                    if not isinstance(field,_Tensor) or field.variance != (1,):
                        if not isinstance(field,_Tensor) and _native(field) == 0:
                            value = Integer(0)
                            continue
                        raise ValueError('one-form slots require vector fields')
                    vector = _express(field,slot.chart,variables=True)
                    value *= vector._components.get((_Slot(1,slot.chart,slot.index),),Integer(0))
            key = tuple(remaining)
            components[key] = components.get(key,Integer(0))+value
        return _tensor(components,remaining_variance)

    rcall = __call__


class _VectorField(_Tensor):
    def __call__(self,scalar):
        if scalar is None:
            return self
        return super().__call__(scalar)

    rcall = __call__


class BaseVectorField(_VectorField):
    def __new__(cls,coord_sys,index,**kwargs):
        index = operator.index(index)
        if kwargs or not 0 <= index < coord_sys.dim:
            raise ValueError('invalid coordinate vector index or options')
        return _tensor({(_Slot(1,coord_sys,index),):Integer(1)},(1,),cls,
                       coord_sys=coord_sys,index=index)

    def __init__(self,*args,**kwargs):
        pass

    patch = property(lambda self:self.coord_sys.patch)
    manifold = property(lambda self:self.coord_sys.manifold)
    dim = property(lambda self:self.coord_sys.dim)


def _chart(tensor):
    if not isinstance(tensor,_Tensor):
        charts = _field_charts(tensor)
    else:
        charts = {slot.chart for slots in tensor._components for slot in slots}
    if len(charts) != 1:
        raise ValueError('an unambiguous single coordinate basis is required')
    return next(iter(charts))


def _express(tensor,chart,variables=False):
    if not isinstance(tensor,_Tensor):
        return _scalar_in_chart(tensor,chart) if variables else _native(tensor)
    components = {}
    for slots,coefficient in tensor._components.items():
        partial = {():_scalar_in_chart(coefficient,chart) if variables else coefficient}
        for slot in slots:
            if slot.chart == chart:
                factors = {slot.index:Integer(1)}
            elif slot.variance == 1:
                jacobian = slot.chart.jacobian(chart,slot.chart.base_scalars())
                factors = {i:jacobian[i,slot.index] for i in range(chart.dim)}
            else:
                jacobian = chart.jacobian(slot.chart,chart.base_scalars())
                factors = {i:jacobian[slot.index,i] for i in range(chart.dim)}
            new = {}
            for prefix,c in partial.items():
                for i,factor in factors.items():
                    if variables:
                        factor = _scalar_in_chart(factor,chart)
                    if factor != 0:
                        key = prefix+(_Slot(slot.variance,chart,i),)
                        new[key] = new.get(key,Integer(0))+c*factor
            partial = new
        for key,c in partial.items():
            components[key] = components.get(key,Integer(0))+c
    return _tensor(components,tensor.variance)


def vectors_in_basis(expr,to_sys):
    if isinstance(expr,_Tensor) and any(v != 1 for v in expr.variance):
        raise ValueError('vectors_in_basis accepts contravariant fields; covariant basis changes are internal')
    return _express(expr,to_sys)


def contravariant_order(expr,_strict=False):
    return expr.variance.count(1) if isinstance(expr,_Tensor) else 0


def covariant_order(expr,_strict=False):
    return expr.variance.count(-1) if isinstance(expr,_Tensor) else 0


class TensorProduct(_Tensor):
    def __new__(cls,*args):
        partial,variance = {():Integer(1)},()
        for factor in args:
            if not isinstance(factor,_Tensor):
                scalar = _native(factor)
                partial = {slots:c*scalar for slots,c in partial.items()}
                continue
            new = {}
            for left,a in partial.items():
                for right,b in factor._components.items():
                    new[left+right] = new.get(left+right,Integer(0))+a*b
            partial = new
            variance += factor.variance
        return _tensor(partial,variance,cls if len(variance) > 1 else None)

    def __init__(self,*args):
        pass


def _parity(permutation):
    return -1 if sum(permutation[i] > permutation[j] for i in range(len(permutation)) for j in range(i+1,len(permutation))) % 2 else 1


class WedgeProduct(TensorProduct):
    def __new__(cls,*args):
        tensor = TensorProduct(*args)
        if not isinstance(tensor,_Tensor):
            return tensor
        if len(set(tensor.variance)) != 1:
            raise ValueError('wedge products require uniformly covariant or contravariant factors')
        n = len(tensor.variance)
        denominator = math.prod(math.factorial(len(factor.variance)) for factor in args if isinstance(factor,_Tensor))
        components = {}
        for permutation in permutations(range(n)):
            multiplier = Rational(_parity(permutation),denominator)
            for slots,c in tensor._components.items():
                key = tuple(slots[i] for i in permutation)
                components[key] = components.get(key,Integer(0))+multiplier*c
        return _tensor(components,tensor.variance,cls)


class Differential(_Tensor):
    def __new__(cls,form_field):
        if isinstance(form_field,Differential):
            return Integer(0)
        if not isinstance(form_field,_Tensor):
            scalar = _native(form_field)
            charts = _field_charts(scalar)
            if not charts:
                return Integer(0)
            chart = sorted(charts,key=lambda c:(c.name,repr(c._key)))[0]
            scalar = _scalar_in_chart(scalar,chart)
            components = {(_Slot(-1,chart,i),):scalar.diff(_native(field)) for i,field in enumerate(chart.base_scalars())}
            return _tensor(components,(-1,),cls,form_field=form_field)
        if any(v != -1 for v in form_field.variance):
            raise ValueError('exterior differentiation requires a scalar or covariant form')
        chart = _chart(form_field)
        form_field = _express(form_field,chart,variables=True)
        degree = len(form_field.variance)
        if degree > 1:
            for slots,c in form_field._components.items():
                for i in range(degree-1):
                    swapped = slots[:i]+(slots[i+1],slots[i])+slots[i+2:]
                    if trigsimp(simplify(c+form_field._components.get(swapped,Integer(0)))) != 0:
                        raise ValueError('exterior differentiation of higher tensors requires an alternating form')
        components = {}
        if degree >= chart.dim:
            return Integer(0)
        for indices in product(range(chart.dim),repeat=degree+1):
            if len(set(indices)) != len(indices):
                continue
            value = Integer(0)
            for i,index in enumerate(indices):
                slots = tuple(_Slot(-1,chart,j) for j in indices[:i]+indices[i+1:])
                coefficient = form_field._components.get(slots,Integer(0))
                value += (-1)**i*coefficient.diff(_native(chart.base_scalar(index)))
            components[tuple(_Slot(-1,chart,j) for j in indices)] = value
        return _tensor(components,(-1,)*(degree+1),cls,form_field=form_field)

    def __init__(self,*args):
        pass


def _vector(value):
    if not isinstance(value,_Tensor) or value.variance != (1,):
        raise ValueError('expected a vector field')
    return value


class Commutator(_VectorField):
    def __new__(cls,v1,v2):
        if not isinstance(v1,_Tensor) and _native(v1) == 0:
            if isinstance(v2,_Tensor):
                _vector(v2)
            elif _native(v2) != 0:
                raise ValueError('expected a vector field')
            return Integer(0)
        if not isinstance(v2,_Tensor) and _native(v2) == 0:
            _vector(v1)
            return Integer(0)
        v1,v2 = _vector(v1),_vector(v2)
        chart = _chart(v1)
        left,right = _express(v1,chart,True),_express(v2,chart,True)
        components = {}
        for i in range(chart.dim):
            slot = (_Slot(1,chart,i),)
            components[slot] = left(right._components.get(slot,Integer(0)))-right(left._components.get(slot,Integer(0)))
        return _tensor(components,(1,),cls,v1=v1,v2=v2)

    def __init__(self,*args):
        pass


class LieDerivative(_Tensor):
    def __new__(cls,v_field,expr):
        vector = _vector(v_field)
        if not isinstance(expr,_Tensor):
            return vector(expr)
        chart = _chart(vector)
        vector,expr = _express(vector,chart,True),_express(expr,chart,True)
        components = {}
        for slots,c in expr._components.items():
            components[slots] = components.get(slots,Integer(0))+vector(c)
            for position,slot in enumerate(slots):
                for index in range(chart.dim):
                    if slot.variance == 1:
                        value = -vector._components.get((_Slot(1,chart,index),),Integer(0)).diff(_native(chart.base_scalar(slot.index)))
                    else:
                        value = vector._components.get((_Slot(1,chart,slot.index),),Integer(0)).diff(_native(chart.base_scalar(index)))
                    key = slots[:position]+(_Slot(slot.variance,chart,index),)+slots[position+1:]
                    components[key] = components.get(key,Integer(0))+c*value
        return _tensor(components,expr.variance,cls if expr.variance != (1,) else _VectorField,v_field=v_field,expr=expr)

    def __init__(self,*args):
        pass
