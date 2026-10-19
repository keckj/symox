"""Focus/directrix parabolas with exact affine constructions."""
from .entity import GeometryEntity, GeometryError, construct, scalar, clean, zero, sign
from .point import Point
from .line import Line, LinearEntity


class Parabola(GeometryEntity):
    __slots__ = ()

    def __new__(cls, focus, directrix, **kwargs):
        kwargs.pop('evaluate', None)
        if kwargs:
            raise TypeError('Unknown parabola options')
        focus = Point(focus, dim=2)
        if not isinstance(directrix, Line) or directrix.ambient_dimension != 2:
            raise ValueError('A parabola directrix must be a planar line')
        if directrix.contains(focus):
            raise GeometryError('The focus cannot lie on the directrix')
        return construct(cls, (focus, directrix))

    @property
    def focus(self):
        return self.args[0]

    @property
    def directrix(self):
        return self.args[1]

    @property
    def vertex(self):
        return self.focus.midpoint(self.directrix.projection(self.focus))

    @property
    def axis(self):
        return self.directrix.perpendicular_line(self.focus)

    @property
    def focal_length(self):
        return self.focus.distance(self.vertex)

    @property
    def eccentricity(self):
        return scalar(1)

    @property
    def p_parameter(self):
        # The signed parameter follows the increasing coordinate direction
        # for the standard horizontal and vertical parabolas.
        direction = self.focus-self.vertex
        if zero(direction.x):
            return direction.y
        if zero(direction.y):
            return direction.x
        return self.focal_length

    @property
    def latus_rectum(self):
        from .line import Segment
        direction = self.directrix.direction/abs(self.directrix.direction)
        return Segment(self.focus-direction*2*self.focal_length, self.focus+direction*2*self.focal_length)

    def scale(self, x=1, y=1, pt=None):
        x, y = scalar(x), scalar(y)
        if not zero(x*x-y*y) or zero(x):
            raise NotImplementedError('Parabola scaling currently requires a nonsingular similarity')
        return Parabola(self.focus.scale(x, y, pt=pt), self.directrix.scale(x, y, pt=pt))

    def equation(self, x='x', y='y'):
        from .. import Symbol
        x = Symbol(x, real=True) if isinstance(x, str) else scalar(x)
        y = Symbol(y, real=True) if isinstance(y, str) else scalar(y)
        point = Point(x, y)
        direction = self.directrix.direction
        normal = Point(-direction.y, direction.x)
        return clean((point-self.directrix.p1).dot(normal)**2/normal.dot(normal)
                     -(point-self.focus).dot(point-self.focus))

    def contains(self, other):
        other = Point(other, dim=2)
        return zero(self.equation(*other))

    def intersection(self, o):
        from .. import Dummy, solve
        if isinstance(o, Point):
            return [o] if self.contains(o) else []
        if not isinstance(o, LinearEntity) or o.ambient_dimension != 2:
            raise NotImplementedError('Parabola intersections currently accept points and planar linear entities')
        parameter = Dummy('parabola_line_parameter', real=True)
        point = o.arbitrary_point(parameter)
        roots = solve(self.equation(*point), parameter)
        result = []
        for root in roots:
            if root.is_real is False:
                continue
            if root.is_real is not True:
                raise NotImplementedError('Parabola intersection reality is undecidable')
            if o._accept(root):
                result.append(o.p1+o.direction*root)
        from .ellipse import _sort_points
        return _sort_points(result)
