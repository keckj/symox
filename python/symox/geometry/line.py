"""Exact affine lines, rays and closed segments in two and three dimensions."""
from .entity import GeometrySet, GeometryEntity, GeometryError, construct, scalar, clean, zero, sign, compare
from .point import Point, Point2D, Point3D


def cross2(a, b):
    return clean(a[0]*b[1]-a[1]*b[0])


class LinearEntity(GeometrySet):
    __slots__ = ()
    _kind = 'line'

    def __new__(cls, p1, p2=None, **kwargs):
        from .. import oo
        kwargs.pop('evaluate', None)
        p1 = Point(p1)
        if p2 is None:
            if 'slope' in kwargs:
                if len(p1) != 2:
                    raise ValueError('Slope is only defined in two dimensions')
                slope = scalar(kwargs.pop('slope'))
                p2 = p1 + (Point(0, 1) if slope in (oo, -oo) else Point(1, slope))
            elif 'direction_ratio' in kwargs:
                p2 = p1 + Point(kwargs.pop('direction_ratio'), dim=len(p1))
            else:
                raise ValueError('Supply two points, a slope, or a direction ratio')
        if kwargs:
            raise TypeError('Unknown linear geometry options: '+', '.join(kwargs))
        p2 = Point(p2)
        dim = 2 if cls.__name__.endswith('2D') else 3 if cls.__name__.endswith('3D') else max(len(p1), len(p2))
        p1, p2 = Point(p1, dim=dim), Point(p2, dim=dim)
        if p1.equals(p2):
            if cls._kind == 'segment':
                return p1
            raise ValueError('Two distinct points are required')
        concrete = globals()[{'line':'Line', 'ray':'Ray', 'segment':'Segment'}[cls._kind]+str(dim)+'D']
        return construct(concrete, (p1, p2))

    @property
    def p1(self):
        return self.args[0]

    @property
    def p2(self):
        return self.args[1]

    @property
    def points(self):
        return self.args

    @property
    def ambient_dimension(self):
        return len(self.p1)

    @property
    def direction(self):
        return self.p2-self.p1

    @property
    def direction_ratio(self):
        return list(self.direction)

    @property
    def direction_cosine(self):
        return list(self.direction/abs(self.direction))

    @property
    def length(self):
        from .. import oo
        return self.p1.distance(self.p2) if self._kind == 'segment' else oo

    @property
    def slope(self):
        from .. import oo
        if self.ambient_dimension != 2:
            raise AttributeError('Spatial lines do not have a scalar slope')
        dx, dy = self.direction
        return oo if zero(dx) else clean(dy/dx)

    @property
    def _parameter_bounds(self):
        from .. import oo
        return (scalar(0), scalar(1)) if self._kind == 'segment' else (
            (scalar(0), oo) if self._kind == 'ray' else (-oo, oo))

    def _parameter(self, point):
        point = Point(point, dim=self.ambient_dimension)
        direction = self.direction
        for i, delta in enumerate(direction):
            if not zero(delta):
                parameter = clean((point[i]-self.p1[i])/delta)
                if all(zero(value) for value in point-self.p1-direction*parameter):
                    return parameter
                return None
        raise GeometryError('A line has no nonzero direction')

    def _accept(self, t):
        if self._kind == 'line':
            return True
        if compare(t, scalar(0)) < 0:
            return False
        return self._kind != 'segment' or compare(t, scalar(1)) <= 0

    def contains(self, other):
        if isinstance(other, LinearEntity):
            if other._kind == 'line' and self._kind != 'line':
                return False
            if other._kind == 'ray':
                if self._kind == 'segment' or not self.contains(other.p1):
                    return False
                return self._kind == 'line' or sign(self.direction.dot(other.direction)) > 0
            return self.contains(other.p1) and self.contains(other.p2)
        try:
            parameter = self._parameter(other)
        except (TypeError, ValueError):
            return False
        return parameter is not None and self._accept(parameter)

    def equals(self, other):
        if not isinstance(other, LinearEntity) or self._kind != other._kind:
            return False
        if self._kind == 'segment':
            return (self.p1.equals(other.p1) and self.p2.equals(other.p2)
                    or self.p1.equals(other.p2) and self.p2.equals(other.p1))
        return self.contains(other) and other.contains(self)

    def is_parallel(self, other):
        if not isinstance(other, LinearEntity):
            raise TypeError('Parallelism requires linear entities')
        a, b = self.direction, other.direction
        if len(a) != len(b):
            raise ValueError('Different ambient dimensions')
        return all(zero(a[i]*b[j]-a[j]*b[i]) for i in range(len(a)) for j in range(i+1, len(a)))

    def is_similar(self, other):
        return isinstance(other, LinearEntity) and self._kind == other._kind

    def is_perpendicular(self, other):
        return zero(self.direction.dot(other.direction))

    def angle_between(self, other):
        from .. import acos
        return acos(clean(self.direction.dot(other.direction)/(abs(self.direction)*abs(other.direction))))

    def smallest_angle_between(self, other):
        from .. import acos, Abs
        return acos(Abs(clean(self.direction.dot(other.direction)/(abs(self.direction)*abs(other.direction)))))

    def parallel_line(self, p):
        p = Point(p, dim=self.ambient_dimension)
        return Line(p, p+self.direction)

    def perpendicular_line(self, p):
        p = Point(p, dim=self.ambient_dimension)
        projection = self.projection(p)
        if p.equals(projection):
            if self.ambient_dimension != 2:
                raise NotImplementedError('A spatial perpendicular through an incident point is not unique')
            dx, dy = self.direction
            projection = p+Point(-dy, dx)
        return Line(p, projection)

    def perpendicular_segment(self, p):
        p = Point(p, dim=self.ambient_dimension)
        return Segment(p, self.projection(p))

    def projection(self, other):
        if isinstance(other, LinearEntity):
            a, b = self.projection(other.p1), self.projection(other.p2)
            return a if a.equals(b) else other.func(a, b)
        point = Point(other, dim=self.ambient_dimension)
        direction = self.direction
        parameter = clean((point-self.p1).dot(direction)/direction.dot(direction))
        return self.p1+direction*parameter

    def distance(self, other):
        from .. import Min, Abs
        if isinstance(other, LinearEntity):
            if self.intersection(other):
                return scalar(0)
            if self._kind == other._kind == 'line':
                if self.is_parallel(other) or self.ambient_dimension == 2:
                    return self.distance(other.p1)
                normal = self.direction.cross(other.direction)
                return Abs((other.p1-self.p1).dot(normal))/abs(normal)
            candidates = []
            # Spatial segments and rays can attain their minimum at two
            # interior points on skew supporting lines.
            a, b, w = self.direction, other.direction, self.p1-other.p1
            aa, ab, bb = a.dot(a), a.dot(b), b.dot(b)
            aw, bw = a.dot(w), b.dot(w)
            determinant = clean(aa*bb-ab*ab)
            if not zero(determinant):
                t, u = clean((ab*bw-bb*aw)/determinant), clean((aa*bw-ab*aw)/determinant)
                if self._accept(t) and other._accept(u):
                    candidates.append((self.p1+a*t).distance(other.p1+b*u))
            if self._kind != 'line':
                candidates.append(other.distance(self.p1))
                if self._kind == 'segment':
                    candidates.append(other.distance(self.p2))
            if other._kind != 'line':
                candidates.append(self.distance(other.p1))
                if other._kind == 'segment':
                    candidates.append(self.distance(other.p2))
            return Min(*candidates)
        point = Point(other, dim=self.ambient_dimension)
        projection = self.projection(point)
        t = self._parameter(projection)
        if self._kind != 'line' and compare(t, scalar(0)) < 0:
            return point.distance(self.p1)
        if self._kind == 'segment' and compare(t, scalar(1)) > 0:
            return point.distance(self.p2)
        return clean(point.distance(projection))

    def equation(self, x='x', y='y', z='z'):
        from .. import Symbol
        coordinates = [Symbol(value, real=True) if isinstance(value, str) else scalar(value)
                       for value in (x, y, z)[:self.ambient_dimension]]
        if self.ambient_dimension == 2:
            dx, dy = self.direction
            if zero(dx):
                return coordinates[0]-self.p1.x
            if zero(dy):
                return coordinates[1]-self.p1.y
            return clean(dx*(coordinates[1]-self.p1.y)-dy*(coordinates[0]-self.p1.x))
        base = self.p1
        direction = self.direction
        pivot = next(i for i, value in enumerate(direction) if not zero(value))
        return tuple(clean((coordinates[j]-base[j])*direction[pivot]
                           -(coordinates[pivot]-base[pivot])*direction[j])
                     for j in range(3) if j != pivot)

    def arbitrary_point(self, parameter='t'):
        from .. import Symbol
        parameter = Symbol(parameter, real=True) if isinstance(parameter, str) else scalar(parameter)
        if parameter.free_symbols.intersection(self.free_symbols):
            raise ValueError('Parameter conflicts with a coordinate symbol')
        return self.p1+self.direction*parameter

    def intersection(self, other):
        if isinstance(other, Point):
            return [other] if self.contains(other) else []
        if not isinstance(other, LinearEntity):
            if hasattr(other, 'intersection'):
                return other.intersection(self)
            raise TypeError('Intersection requires a geometric entity')
        if self.ambient_dimension != other.ambient_dimension:
            raise ValueError('Different ambient dimensions')
        a, b, offset = self.direction, other.direction, other.p1-self.p1
        for i in range(len(a)):
            for j in range(i+1, len(a)):
                determinant = clean(a[i]*b[j]-a[j]*b[i])
                if zero(determinant):
                    continue
                t = clean((offset[i]*b[j]-offset[j]*b[i])/determinant)
                u = clean((offset[i]*a[j]-offset[j]*a[i])/determinant)
                point = self.p1+a*t
                if not point.equals(other.p1+b*u):
                    return []
                return [point] if self._accept(t) and other._accept(u) else []
        if self._parameter(other.p1) is None:
            return []
        return self._collinear_intersection(other)

    def _collinear_intersection(self, other):
        from .. import oo
        origin = self._parameter(other.p1)
        pivot = next(i for i, value in enumerate(self.direction) if not zero(value))
        ratio = clean(other.direction[pivot]/self.direction[pivot])
        low, high = self._parameter_bounds
        other_low, other_high = other._parameter_bounds
        if other._kind == 'line':
            other_low, other_high = -oo, oo
        elif other._kind == 'ray':
            other_low, other_high = (origin, oo) if sign(ratio) > 0 else (-oo, origin)
        else:
            end = clean(origin+ratio)
            other_low, other_high = (origin, end) if sign(ratio) > 0 else (end, origin)
        low = low if compare(low, other_low) >= 0 else other_low
        high = high if compare(high, other_high) <= 0 else other_high
        relation = compare(low, high)
        if relation > 0:
            return []
        if relation == 0:
            return [self.p1+self.direction*low]
        if low == -oo and high == oo:
            return [self]
        if low == -oo:
            start = self.p1+self.direction*high
            return [Ray(start, start-self.direction)]
        start = self.p1+self.direction*low
        if high == oo:
            return [Ray(start, start+self.direction)]
        return [Segment(start, self.p1+self.direction*high)]


class LinearEntity2D(LinearEntity):
    __slots__=()


class LinearEntity3D(LinearEntity):
    __slots__=()


class Line(LinearEntity):
    __slots__ = ()


class Ray(LinearEntity):
    __slots__ = ()
    _kind = 'ray'

    @property
    def source(self):
        return self.p1


class Segment(LinearEntity):
    __slots__ = ()
    _kind = 'segment'

    @property
    def midpoint(self):
        return self.p1.midpoint(self.p2)

    def perpendicular_bisector(self, p=None):
        line = Line(self.p1, self.p2).perpendicular_line(self.midpoint)
        if p is not None and line.contains(Point(p)):
            return Segment(self.midpoint, Point(p))
        return line


class Line2D(Line, LinearEntity2D):
    __slots__ = ()


class Line3D(Line, LinearEntity3D):
    __slots__ = ()


class Ray2D(Ray, LinearEntity2D):
    __slots__ = ()

    @property
    def xdirection(self):
        from .. import oo
        return sign(self.direction.x)*oo if not zero(self.direction.x) else scalar(0)

    @property
    def ydirection(self):
        from .. import oo
        return sign(self.direction.y)*oo if not zero(self.direction.y) else scalar(0)


class Ray3D(Ray, LinearEntity3D):
    __slots__ = ()


class Segment2D(Segment, LinearEntity2D):
    __slots__ = ()


class Segment3D(Segment, LinearEntity3D):
    __slots__ = ()
