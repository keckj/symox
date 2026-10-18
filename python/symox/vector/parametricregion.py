"""Parameterized curves, surfaces and volumes, including geometry adapters."""
from .. import Symbol, sympify, pi, sin, cos, sqrt


class ParametricRegion:
    def __init__(self,definition,*bounds):
        if getattr(definition,'is_Tuple',False):
            definition = definition.args
        if not isinstance(definition,(tuple,list)):
            definition = (definition,)
        self.definition = tuple(sympify(v) for v in definition)
        if not 1 <= len(self.definition) <= 3:
            raise ValueError('regions require one to three coordinate expressions')
        self.parameters = ()
        self.limits = {}
        for bound in bounds:
            if getattr(bound,'is_Tuple',False):
                bound = bound.args
            if isinstance(bound,(tuple,list)):
                if len(bound) != 3:
                    raise ValueError('bounds must be (parameter, lower, upper)')
                parameter = sympify(bound[0])
                limits = (sympify(bound[1]),sympify(bound[2]))
                self.limits[parameter] = limits
            else:
                parameter = sympify(bound)
            if not parameter.is_Symbol or parameter in self.parameters:
                raise ValueError('region parameters must be distinct symbols')
            self.parameters += (parameter,)
        if len(self.limits) > len(self.definition):
            raise ValueError('region dimension exceeds its ambient dimension')
        self.dimensions = len(self.limits)
        self.args = (self.definition,)+tuple(bounds)

    def __repr__(self):
        return 'ParametricRegion(%s)' % ', '.join(map(str,self.args))

    def __eq__(self,other):
        return isinstance(other,ParametricRegion) and (self.definition,self.parameters,self.limits) == (other.definition,other.parameters,other.limits)


def parametric_region_list(reg,parameters=None):
    from ..geometry import Point, Segment, Curve, Ellipse, Polygon
    from .implicitregion import ImplicitRegion
    if isinstance(reg,ParametricRegion):
        return [reg]
    t = Symbol('t',real=True) if parameters is None else Symbol(parameters,real=True) if isinstance(parameters,str) else sympify(parameters[0]) if isinstance(parameters,(tuple,list)) else sympify(parameters)
    if isinstance(reg,Point):
        return [ParametricRegion(tuple(reg.args))]
    if isinstance(reg,Segment):
        first,second = reg.points
        return [ParametricRegion(tuple(a+(b-a)*t for a,b in zip(first.args,second.args)),(t,0,1))]
    if isinstance(reg,Curve):
        return [ParametricRegion(tuple(reg.functions),tuple(reg.limits))]
    if isinstance(reg,Ellipse):
        x,y = reg.center.args
        return [ParametricRegion((x+reg.hradius*cos(t),y+reg.vradius*sin(t)),(t,0,2*pi))]
    if isinstance(reg,Polygon):
        return [r for side in reg.sides for r in parametric_region_list(side,t)]
    if isinstance(reg,ImplicitRegion):
        return reg._parametric_regions(parameters)
    raise ValueError('cannot determine a parametrization for this region')
