"""Axis-aligned ellipses and exact circle intersections."""
from functools import cmp_to_key

from .entity import GeometryEntity, GeometryError, construct, scalar, clean, zero, sign, compare, unique
from .point import Point, Point2D
from .line import LinearEntity, Line, Segment, Ray


def _sort_points(points):
    def order(a, b):
        try:
            for x, y in zip(a, b):
                compared = compare(x, y)
                if compared:
                    return compared
            return 0
        except NotImplementedError:
            from .. import default_sort_key
            akey, bkey = tuple(map(default_sort_key, a)), tuple(map(default_sort_key, b))
            return (akey > bkey)-(akey < bkey)
    return sorted(unique(points), key=cmp_to_key(order))


class Ellipse(GeometryEntity):
    __slots__ = ()

    def __new__(cls, center=None, hradius=None, vradius=None, eccentricity=None, **kwargs):
        from .. import sqrt
        kwargs.pop('evaluate', None)
        if kwargs:
            raise TypeError('Unknown ellipse options')
        center = Point2D(0, 0) if center is None else Point(center, dim=2)
        if eccentricity is not None:
            eccentricity = scalar(eccentricity)
            if hradius is None and vradius is not None:
                hradius = scalar(vradius)/sqrt(1-eccentricity**2)
            elif vradius is None and hradius is not None:
                vradius = scalar(hradius)*sqrt(1-eccentricity**2)
            else:
                raise ValueError('Supply exactly one radius with eccentricity')
        if hradius is None or vradius is None:
            raise ValueError('Both ellipse radii are required')
        hradius, vradius = scalar(hradius), scalar(vradius)
        if hradius.is_real is False or vradius.is_real is False:
            raise ValueError('Ellipse radii must be real')
        if zero(hradius) or zero(vradius):
            raise GeometryError('Degenerate ellipses are not supported')
        if zero(hradius-vradius):
            return Circle(center, hradius)
        return construct(cls, (center, hradius, vradius))

    @property
    def center(self):
        return self.args[0]

    @property
    def hradius(self):
        return self.args[1]

    @property
    def vradius(self):
        return self.args[2]

    @property
    def area(self):
        from .. import pi
        return pi*self.hradius*self.vradius

    @property
    def bounds(self):
        from .. import Abs
        a, b = Abs(self.hradius), Abs(self.vradius)
        return self.center.x-a, self.center.y-b, self.center.x+a, self.center.y+b

    @property
    def major(self):
        from .. import Abs
        a, b = Abs(self.hradius), Abs(self.vradius)
        return a if compare(a, b) >= 0 else b

    @property
    def minor(self):
        from .. import Abs
        a, b = Abs(self.hradius), Abs(self.vradius)
        return a if compare(a, b) <= 0 else b

    @property
    def focus_distance(self):
        from .. import sqrt
        return sqrt(self.major**2-self.minor**2)

    @property
    def eccentricity(self):
        return self.focus_distance/self.major

    @property
    def foci(self):
        from .. import Abs
        offset = Point(self.focus_distance, 0) if compare(Abs(self.hradius), Abs(self.vradius)) >= 0 else Point(0, self.focus_distance)
        return self.center-offset, self.center+offset

    @property
    def circumference(self):
        from .. import elliptic_e
        return 4*self.major*elliptic_e(self.eccentricity**2)

    def equation(self, x='x', y='y'):
        from .. import Symbol
        x = Symbol(x, real=True) if isinstance(x, str) else scalar(x)
        y = Symbol(y, real=True) if isinstance(y, str) else scalar(y)
        return (x-self.center.x)**2/self.hradius**2+(y-self.center.y)**2/self.vradius**2-1

    def contains(self, point):
        if isinstance(point, Ellipse):
            return self.equals(point)
        try:
            point = Point(point, dim=2)
        except (ValueError, TypeError):
            return False
        return zero(self.equation(point.x, point.y))

    def equals(self, other):
        return (isinstance(other, Ellipse) and self.center.equals(other.center)
                and zero(self.hradius**2-other.hradius**2) and zero(self.vradius**2-other.vradius**2))

    def encloses_point(self, point):
        point = Point(point, dim=2)
        return sign(self.equation(point.x, point.y)) < 0

    def arbitrary_point(self, parameter='t'):
        from .. import Symbol, cos, sin
        parameter = Symbol(parameter, real=True) if isinstance(parameter, str) else scalar(parameter)
        if parameter.free_symbols.intersection(self.free_symbols):
            raise ValueError('Parameter conflicts with a coordinate symbol')
        return self.center+Point(self.hradius*cos(parameter), self.vradius*sin(parameter))

    def tangent_lines(self, point):
        from .. import sqrt
        point = Point(point, dim=2)
        level = clean(self.equation(point.x, point.y))
        relation = sign(level)
        if relation < 0:
            return []
        dx, dy = point-self.center
        if relation == 0:
            return [Line(point, point+Point(-dy/self.vradius**2, dx/self.hradius**2))]
        # Normalize to a unit circle and recover its two contact points.
        u, v = dx/self.hradius, dy/self.vradius
        norm = clean(u*u+v*v)
        root = sqrt(norm-1)
        contacts = [self.center+Point(self.hradius*(u-s*v*root)/norm,
                                     self.vradius*(v+s*u*root)/norm) for s in (-1, 1)]
        return [Line(point, contact) for contact in _sort_points(contacts)]

    def is_tangent(self, other):
        result = self.intersection(other)
        return len(result) == 1 and isinstance(result[0], Point)

    def distance(self, other):
        raise NotImplementedError('General point-to-ellipse distance needs algebraic root isolation')

    def is_similar(self, other):
        if not isinstance(other, Ellipse):
            return False
        return zero(self.major*other.minor-self.minor*other.major)

    def intersection(self, other):
        from .. import sqrt, S
        if isinstance(other, Point):
            return [other] if self.contains(other) else []
        if isinstance(other, LinearEntity):
            if other.ambient_dimension != 2:
                raise ValueError('Ellipse intersections require planar lines')
            p, d = other.p1-self.center, other.direction
            a = clean(d.x**2/self.hradius**2+d.y**2/self.vradius**2)
            b = clean(2*(p.x*d.x/self.hradius**2+p.y*d.y/self.vradius**2))
            c = clean(p.x**2/self.hradius**2+p.y**2/self.vradius**2-1)
            disc = clean(b*b-4*a*c)
            if sign(disc) < 0:
                return []
            roots = unique([clean((-b-sqrt(disc))/(2*a)), clean((-b+sqrt(disc))/(2*a))])
            return _sort_points([other.p1+d*t for t in roots if other._accept(t)])
        if self.equals(other):
            return [self]
        if isinstance(other, Ellipse):
            raise NotImplementedError('General ellipse-ellipse intersections require isolated algebraic roots')
        if hasattr(other, 'intersection'):
            return other.intersection(self)
        raise TypeError('Intersection requires a geometric entity')

    def scale(self, x=1, y=1, pt=None):
        return Ellipse(self.center.scale(x, y, pt=pt), self.hradius*scalar(x), self.vradius*scalar(y))

    def rotate(self, angle, pt=None):
        from .. import sin, cos
        angle = scalar(angle)
        if zero(sin(angle)):
            return Ellipse(self.center.rotate(angle, pt), self.hradius, self.vradius)
        if zero(cos(angle)):
            return Ellipse(self.center.rotate(angle, pt), self.vradius, self.hradius)
        raise NotImplementedError('Rotated non-axis-aligned ellipses need an implicit conic representation')


class Circle(Ellipse):
    __slots__ = ()

    def __new__(cls, center, radius=None, third=None, **kwargs):
        from .. import sqrt
        kwargs.pop('evaluate', None)
        if kwargs:
            raise TypeError('Unknown circle options')
        if third is not None:
            a, b, c = Point(center, dim=2), Point(radius, dim=2), Point(third, dim=2)
            determinant = 2*(a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y))
            if zero(determinant):
                raise GeometryError('Three collinear points do not determine a circle')
            aa, bb, cc = a.dot(a), b.dot(b), c.dot(c)
            center = Point(clean((aa*(b.y-c.y)+bb*(c.y-a.y)+cc*(a.y-b.y))/determinant),
                           clean((aa*(c.x-b.x)+bb*(a.x-c.x)+cc*(b.x-a.x))/determinant))
            radius = center.distance(a)
        else:
            center = Point(center, dim=2)
            if radius is None:
                raise ValueError('A circle requires a radius')
            radius = scalar(radius)
        if radius.is_real is False or radius.is_finite is False:
            raise ValueError('Circle radius must be finite and real')
        if zero(radius):
            return center
        return construct(cls, (center, radius))

    @property
    def radius(self):
        return self.args[1]

    @property
    def hradius(self):
        return self.radius

    @property
    def vradius(self):
        return self.radius

    @property
    def circumference(self):
        from .. import pi
        return 2*pi*self.radius

    @property
    def eccentricity(self):
        return scalar(0)

    @property
    def foci(self):
        return self.center, self.center

    def equation(self, x='x', y='y'):
        from .. import Symbol
        x = Symbol(x, real=True) if isinstance(x, str) else scalar(x)
        y = Symbol(y, real=True) if isinstance(y, str) else scalar(y)
        return (x-self.center.x)**2+(y-self.center.y)**2-self.radius**2

    def distance(self, other):
        from .. import Abs
        if isinstance(other, Point) or not isinstance(other, GeometryEntity):
            return Abs(self.center.distance(Point(other, dim=2))-Abs(self.radius))
        raise NotImplementedError('Circle distance currently accepts points')

    def intersection(self, other):
        from .. import sqrt, Abs
        if not isinstance(other, Circle):
            return super().intersection(other)
        offset = other.center-self.center
        dist_squared = offset.dot(offset)
        r1, r2 = Abs(self.radius), Abs(other.radius)
        if zero(dist_squared):
            return [self] if zero(r1-r2) else []
        along = clean((r1*r1-r2*r2+dist_squared)/(2*dist_squared))
        height_squared = clean(r1*r1/dist_squared-along*along)
        if sign(height_squared) < 0:
            return []
        base = self.center+offset*along
        normal = Point(-offset.y, offset.x)*sqrt(height_squared)
        return _sort_points([base-normal, base+normal])

    def rotate(self, angle, pt=None):
        return Circle(self.center.rotate(angle, pt), self.radius)

    def scale(self, x=1, y=1, pt=None):
        from .. import Abs
        x, y = scalar(x), scalar(y)
        center = self.center.scale(x, y, pt=pt)
        if zero(x*x-y*y):
            return Circle(center, self.radius*Abs(x))
        return Ellipse(center, self.radius*x, self.radius*y)
