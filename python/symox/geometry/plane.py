"""Exact affine planes and spatial intersections."""
from .entity import GeometryEntity, GeometryError, construct, scalar, clean, zero
from .point import Point, Point3D
from .line import LinearEntity, Line


class Plane(GeometryEntity):
    __slots__ = ()

    def __new__(cls, p1, a=None, b=None, normal_vector=None, **kwargs):
        from .. import Tuple
        kwargs.pop('evaluate', None)
        if kwargs:
            raise TypeError('Unknown plane options')
        p1 = Point(p1, dim=3)
        if normal_vector is None:
            if b is None:
                normal_vector = a
            else:
                normal_vector = (Point(a, dim=3)-p1).cross(Point(b, dim=3)-p1)
        if normal_vector is None:
            raise ValueError('Supply a normal vector or three noncollinear points')
        normal = Point(normal_vector, dim=3)
        if all(zero(value) for value in normal):
            raise ValueError('A plane normal must be nonzero')
        return construct(cls, (p1, Tuple(*normal)))

    @property
    def p1(self):
        return self.args[0]

    @property
    def normal_vector(self):
        return self.args[1].args

    @property
    def _normal(self):
        return Point(self.normal_vector, dim=3)

    def equation(self, x='x', y='y', z='z'):
        from .. import Symbol
        coordinates = [Symbol(value, real=True) if isinstance(value, str) else scalar(value) for value in (x, y, z)]
        return clean((Point(*coordinates)-self.p1).dot(self._normal))

    def contains(self, other):
        if isinstance(other, LinearEntity):
            return self.contains(other.p1) and self.contains(other.p2)
        if isinstance(other, Plane):
            return self.equals(other)
        try:
            point = Point(other, dim=3)
        except (ValueError, TypeError):
            return False
        return zero((point-self.p1).dot(self._normal))

    def equals(self, o):
        return isinstance(o, Plane) and self.is_parallel(o) and self.contains(o.p1)

    def is_parallel(self, l):
        if isinstance(l, Plane):
            return all(zero(value) for value in self._normal.cross(l._normal))
        if isinstance(l, LinearEntity):
            return zero(self._normal.dot(Point(l.direction, dim=3)))
        raise TypeError('Plane parallelism requires a plane or linear entity')

    def is_perpendicular(self, l):
        if isinstance(l, Plane):
            return zero(self._normal.dot(l._normal))
        if isinstance(l, LinearEntity):
            return all(zero(value) for value in self._normal.cross(Point(l.direction, dim=3)))
        raise TypeError('Plane perpendicularity requires a plane or linear entity')

    def distance(self, o):
        from .. import Abs
        if isinstance(o, Plane):
            return self.distance(o.p1) if self.is_parallel(o) else scalar(0)
        if isinstance(o, LinearEntity):
            if self.intersection(o):
                return scalar(0)
            if self.is_parallel(o):
                return self.distance(o.p1)
            from .. import Min
            candidates = [self.distance(o.p1)]
            if o._kind == 'segment':
                candidates.append(self.distance(o.p2))
            return Min(*candidates)
        point = Point(o, dim=3)
        return clean(Abs((point-self.p1).dot(self._normal))/abs(self._normal))

    def projection(self, pt):
        pt = Point(pt, dim=3)
        normal = self._normal
        return pt-normal*clean((pt-self.p1).dot(normal)/normal.dot(normal))

    def projection_line(self, line):
        if not isinstance(line, LinearEntity):
            raise TypeError('Expected a linear entity')
        a, b = self.projection(line.p1), self.projection(line.p2)
        return a if a.equals(b) else line.func(a, b)

    def perpendicular_line(self, pt):
        pt = Point(pt, dim=3)
        return Line(pt, pt+self._normal)

    def angle_between(self, o):
        # Reference: the signed angle asin(n.d/(|n||d|)) for spatial linear
        # entities, acos of the normals for planes, None for anything else.
        from .. import acos, asin
        if isinstance(o, Plane):
            return acos(clean(self._normal.dot(o._normal)/(abs(self._normal)*abs(o._normal))))
        if isinstance(o, LinearEntity) and o.ambient_dimension == 3:
            return asin(clean(self._normal.dot(o.direction)/(abs(self._normal)*abs(o.direction))))
        return None

    def intersection(self, o):
        if isinstance(o, Point):
            return [o] if self.contains(o) else []
        if isinstance(o, LinearEntity):
            if o.ambient_dimension != 3:
                raise ValueError('Plane intersections require spatial lines')
            denominator = self._normal.dot(o.direction)
            if zero(denominator):
                return [o] if self.contains(o.p1) else []
            t = clean(self._normal.dot(self.p1-o.p1)/denominator)
            return [o.p1+o.direction*t] if o._accept(t) else []
        if not isinstance(o, Plane):
            raise TypeError('Plane intersections accept points, lines and planes')
        normal, second = self._normal, o._normal
        direction = normal.cross(second)
        if all(zero(value) for value in direction):
            return [self] if self.contains(o.p1) else []
        rhs1, rhs2 = self.p1.dot(normal), o.p1.dot(second)
        for i in range(3):
            for j in range(i+1, 3):
                determinant = normal[i]*second[j]-normal[j]*second[i]
                if not zero(determinant):
                    coordinates = [scalar(0)]*3
                    coordinates[i] = clean((rhs1*second[j]-rhs2*normal[j])/determinant)
                    coordinates[j] = clean((normal[i]*rhs2-second[i]*rhs1)/determinant)
                    base = Point(*coordinates)
                    return [Line(base, base+direction)]

    def translate(self, x=0, y=0, z=0):
        return Plane(self.p1.translate(x, y, z), normal_vector=self.normal_vector)

    def scale(self, x=1, y=1, z=1, pt=None):
        factors = list(map(scalar, (x, y, z)))
        if any(zero(factor) for factor in factors):
            raise NotImplementedError('Singular plane scaling is not supported')
        return Plane(self.p1.scale(x, y, z, pt=pt), normal_vector=[v/f for v, f in zip(self.normal_vector, factors)])
