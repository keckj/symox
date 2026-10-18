"""Levi-Civita connections, curvature and integral curves in local coordinates."""
from itertools import product
import math
import operator

from .. import Integer,Matrix,Function,simplify as _simplify,trigsimp,expand_trig
from ..tensor.array import ImmutableDenseNDimArray
from ._coordinates import _native,_scalar_at,_scalar_in_chart,Point
from ._tensors import (_Tensor,_Slot,_tensor,_chart,_vector,_express,BaseVectorField,
                       covariant_order,contravariant_order)


def simplify(expression):
    return trigsimp(expand_trig(_simplify(expression)))


def twoform_to_matrix(expr):
    if not isinstance(expr,_Tensor) or expr.variance != (-1,-1):
        raise ValueError('expected a covariant tensor of order two')
    chart = _chart(expr)
    # Match SymPy's actual matrix convention: column i, row j evaluates e_i,e_j.
    return Matrix([[expr(chart.base_vector(i),chart.base_vector(j))
                    for i in range(chart.dim)] for j in range(chart.dim)])


def _metric(expr):
    chart = _chart(expr)
    matrix = twoform_to_matrix(expr).applyfunc(lambda c:_scalar_in_chart(c,chart))
    if any(simplify(matrix[i,j]-matrix[j,i]) != 0 for i in range(chart.dim) for j in range(i)):
        raise ValueError('metric tensor must be symmetric')
    return chart,matrix


def metric_to_Christoffel_1st(expr):
    chart,g = _metric(expr)
    x = tuple(_native(f) for f in chart.base_scalars())
    n = chart.dim
    return ImmutableDenseNDimArray([[[simplify((g[i,j].diff(x[k])+g[i,k].diff(x[j])-g[j,k].diff(x[i]))/2)
                                     for k in range(n)] for j in range(n)] for i in range(n)])


def metric_to_Christoffel_2nd(expr):
    chart,g = _metric(expr)
    inverse = g.inv()
    first = metric_to_Christoffel_1st(expr)
    n = chart.dim
    return ImmutableDenseNDimArray([[[simplify(sum(inverse[i,l]*first[l,j,k] for l in range(n)))
                                     for k in range(n)] for j in range(n)] for i in range(n)])


def metric_to_Riemann_components(expr):
    chart = _chart(expr)
    gamma = metric_to_Christoffel_2nd(expr)
    x = tuple(_native(f) for f in chart.base_scalars())
    n = chart.dim
    return ImmutableDenseNDimArray([[[[simplify(gamma[r,s,v].diff(x[u])-gamma[r,s,u].diff(x[v])+
        sum(gamma[r,l,u]*gamma[l,s,v]-gamma[r,l,v]*gamma[l,s,u] for l in range(n)))
        for v in range(n)] for u in range(n)] for s in range(n)] for r in range(n)])


def metric_to_Ricci_components(expr):
    n = _chart(expr).dim
    riemann = metric_to_Riemann_components(expr)
    return ImmutableDenseNDimArray([[simplify(sum(riemann[k,i,k,j] for k in range(n)))
                                    for j in range(n)] for i in range(n)])


class CovarDerivativeOp:
    def __init__(self,wrt,christoffel):
        self.wrt = _vector(wrt)
        self.coord_sys = _chart(wrt)
        self.christoffel = ImmutableDenseNDimArray(christoffel)
        n = self.coord_sys.dim
        if self.christoffel.shape != (n,n,n):
            raise ValueError('Christoffel components must have shape (dim, dim, dim)')

    def __call__(self,field):
        chart = self.coord_sys
        vector = _express(self.wrt,chart,True)
        if not isinstance(field,_Tensor):
            return vector(field)
        field = _express(field,chart,True)
        components = {}
        xvalues = [vector._components.get((_Slot(1,chart,j),),Integer(0)) for j in range(chart.dim)]
        for slots,c in field._components.items():
            components[slots] = components.get(slots,Integer(0))+vector(c)
            for position,slot in enumerate(slots):
                for new in range(chart.dim):
                    if slot.variance == 1:
                        connection = sum(xvalues[j]*self.christoffel[new,j,slot.index] for j in range(chart.dim))
                    else:
                        connection = -sum(xvalues[j]*self.christoffel[slot.index,j,new] for j in range(chart.dim))
                    key = slots[:position]+(_Slot(slot.variance,chart,new),)+slots[position+1:]
                    components[key] = components.get(key,Integer(0))+c*connection
        return _tensor(components,field.variance)


class BaseCovarDerivativeOp(CovarDerivativeOp):
    def __init__(self,coord_sys,index,christoffel):
        self.index = operator.index(index)
        super().__init__(BaseVectorField(coord_sys,self.index),christoffel)


def intcurve_series(vector_field,param,start_point,n=6,coord_sys=None,coeffs=False):
    vector = _vector(vector_field)
    if not isinstance(start_point,Point):
        raise TypeError('initial point must be a differential-geometric Point')
    chart = start_point.coord_sys if coord_sys is None else coord_sys
    vector = _express(vector,chart,True)
    n = operator.index(n)
    if n < 0:
        raise ValueError('series length must be nonnegative')
    param = _native(param)
    current = tuple(_native(f) for f in chart.base_scalars())
    terms = []
    for i in range(n):
        terms.append(Matrix([simplify(_scalar_at(f,start_point)*param**i/math.factorial(i)) for f in current]))
        current = tuple(vector(f) for f in current)
    if coeffs:
        return terms
    return sum(terms,Matrix.zeros(chart.dim,1))


def intcurve_diffequ(vector_field,param,start_point,coord_sys=None):
    vector = _vector(vector_field)
    if not isinstance(start_point,Point):
        raise TypeError('initial point must be a differential-geometric Point')
    chart = start_point.coord_sys if coord_sys is None else coord_sys
    vector = _express(vector,chart,True)
    param = _native(param)
    functions = [Function('f_%s' % i)(param) for i in range(chart.dim)]
    replacements = dict(zip(map(_native,chart.base_scalars()),functions))
    equations = [functions[i].diff(param)-vector(chart.base_scalar(i)).subs(replacements,simultaneous=True)
                 for i in range(chart.dim)]
    initial = [f.subs(param,0)-coordinate for f,coordinate in zip(functions,start_point.coords(chart))]
    return equations,initial
