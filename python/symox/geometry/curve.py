"""Finite-dimensional parametric curves."""
from .entity import GeometrySet, GeometryEntity, construct, scalar
from .point import Point


class Curve(GeometrySet):
    __slots__ = ()

    def __new__(cls, function, limits, **kwargs):
        from .. import Tuple
        kwargs.pop('evaluate', None)
        if kwargs:
            raise TypeError('Unknown curve options')
        if getattr(getattr(function, 'func', None), '__name__', '') == 'Tuple':
            function = function.args
        if getattr(getattr(limits, 'func', None), '__name__', '') == 'Tuple':
            limits = limits.args
        function = tuple(map(scalar, function))
        limits = tuple(map(scalar, limits))
        if len(function) not in (2, 3) or len(limits) != 3 or not limits[0].is_Symbol:
            raise ValueError('A curve needs two or three coordinates and (parameter, start, end)')
        return construct(cls, (Tuple(*function), Tuple(*limits)))

    @property
    def functions(self):
        return self.args[0].args

    @property
    def limits(self):
        return self.args[1].args

    @property
    def parameter(self):
        return self.limits[0]

    @property
    def free_symbols(self):
        result = set().union(*(value.free_symbols for value in self.functions))-{self.parameter}
        return result | self.limits[1].free_symbols | self.limits[2].free_symbols

    @property
    def length(self):
        from .. import sqrt, diff, integrate, S, trigsimp
        speed = sqrt(trigsimp(sum((diff(value, self.parameter)**2 for value in self.functions), S.Zero)))
        return integrate(speed, self.limits)

    def arbitrary_point(self, parameter='t'):
        from .. import Symbol
        if parameter is None:
            return Point(*self.functions)
        parameter = Symbol(parameter, real=True) if isinstance(parameter, str) else scalar(parameter)
        if parameter in self.free_symbols:
            raise ValueError('Parameter conflicts with a free coordinate symbol')
        return Point(*(value.xreplace({self.parameter:parameter}) for value in self.functions))

    def plot_interval(self, parameter='t'):
        from .. import Symbol
        return [Symbol(parameter, real=True) if isinstance(parameter, str) else scalar(parameter), *self.limits[1:]]

    def subs(self,*args,**kwargs):
        from collections.abc import Mapping
        from .. import sympify
        if len(args)==2:
            if args[0]==self:return args[1]
            pairs=[args]
        elif len(args)==1:
            pairs=list(args[0].items()) if isinstance(args[0],Mapping) else list(args[0])
        else:
            return super().subs(*args,**kwargs)
        if any(sympify(old)==self.parameter for old,new in pairs):
            return Point(*(function.subs(*args,**kwargs) for function in self.functions))
        return super().subs(*args,**kwargs)

    def translate(self, x=0, y=0, z=0):
        return Curve(Point(*self.functions).translate(x, y, z).args, self.limits)

    def rotate(self, angle, pt=None):
        return Curve(Point(*self.functions).rotate(angle, pt).args, self.limits)

    def scale(self, x=1, y=1, z=1, pt=None):
        return Curve(Point(*self.functions).scale(x, y, z, pt=pt).args, self.limits)

    def reflect(self, line):
        # SymPy 1.14 returns the curve unchanged (it only maps Point atoms);
        # the mathematically correct image reflects the coordinate functions.
        return Curve(Point(*self.functions, evaluate=False).reflect(line).args, self.limits)
