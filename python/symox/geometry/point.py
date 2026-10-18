"""Finite Euclidean points in two and three dimensions."""
from .entity import GeometryEntity, GeometryError, construct, scalar, clean, zero


class Point(GeometryEntity):
    __slots__ = ()

    def __new__(cls, *args, dim=None, evaluate=True, on_morph='ignore'):
        from .. import nsimplify
        if not args and dim is None:
            raise ValueError('A point requires at least two coordinates')
        if len(args) == 1:
            if isinstance(args[0], Point):
                args = args[0].args
            elif getattr(getattr(args[0], 'func', None), '__name__', '') == 'Tuple':
                args = args[0].args
            else:
                try:
                    args = tuple(args[0])
                except TypeError as exc:
                    raise TypeError('Point requires coordinate values') from exc
        if dim is None:
            dim = 2 if cls is Point2D else 3 if cls is Point3D else len(args)
        if dim not in (2, 3):
            raise NotImplementedError('Points currently support dimensions two and three')
        if cls is Point2D and dim != 2 or cls is Point3D and dim != 3:
            raise ValueError('Coordinate dimension disagrees with the point class')
        if len(args) > dim:
            if any(not zero(scalar(arg)) for arg in args[dim:]):
                raise ValueError('Cannot discard nonzero coordinates')
            args = args[:dim]
        args = args + (0,) * (dim - len(args))
        if len(args) < 2:
            raise ValueError('At least two coordinates are required')
        values = tuple(scalar(arg) for arg in args)
        if any(value.is_extended_real is False for value in values):
            raise ValueError('Euclidean coordinates must be real')
        if any(value.is_finite is False for value in values):
            raise ValueError('Euclidean coordinates must be finite')
        if evaluate:
            values = tuple(nsimplify(value, rational=True) for value in values)
        return construct((Point2D if dim == 2 else Point3D) if cls is Point else cls, values)

    def __len__(self):
        return len(self.args)

    def __iter__(self):
        return iter(self.args)

    def __getitem__(self, index):
        return self.args[index]

    @property
    def x(self):
        return self.args[0]

    @property
    def y(self):
        return self.args[1]

    @property
    def ambient_dimension(self):
        return len(self.args)

    @property
    def _norm(self):
        from .. import sqrt
        return sqrt(sum((value**2 for value in self), scalar(0)))

    @property
    def length(self):
        return scalar(0)

    def __abs__(self):
        return self._norm

    def __add__(self, other):
        other = Point(other)
        dim = max(len(self), len(other))
        a, b = Point(self, dim=dim), Point(other, dim=dim)
        return Point(*(clean(x+y) for x, y in zip(a, b)), evaluate=False)

    def __sub__(self, other):
        return self + (-Point(other))

    def __neg__(self):
        return Point(*(-value for value in self), evaluate=False)

    def __mul__(self, factor):
        factor = scalar(factor)
        return Point(*(clean(value*factor) for value in self), evaluate=False)

    __rmul__ = __mul__

    def __truediv__(self, divisor):
        divisor = scalar(divisor)
        if zero(divisor):
            raise ZeroDivisionError('Cannot divide a point by zero')
        return self * (1/divisor)

    def dot(self, other):
        other = Point(other, dim=len(self))
        return clean(sum((a*b for a, b in zip(self, other)), scalar(0)))

    def distance(self, other):
        if isinstance(other, GeometryEntity) and not isinstance(other, Point):
            return other.distance(self)
        return abs(self - Point(other))

    def taxicab_distance(self, other):
        from .. import Abs
        return sum((Abs(v) for v in self-Point(other)), scalar(0))

    def midpoint(self, other):
        return (self + Point(other))/2

    @property
    def unit(self):
        return self/abs(self)

    def contains(self, other):
        try:
            return all(zero(value) for value in self-Point(other))
        except (TypeError, ValueError):
            return False

    def equals(self, other):
        return self.contains(other)

    def intersection(self, other):
        if isinstance(other, GeometryEntity):
            return [self] if other.contains(self) else []
        return [self] if self.contains(other) else []

    @staticmethod
    def is_collinear(*points):
        points = [Point(point) for point in points]
        if len(points) < 3:
            return True
        dimension = max(map(len, points))
        points = [Point(point, dim=dimension) for point in points]
        first = points[0]
        others = [point-first for point in points if point != first]
        if len(others) < 2:
            return True
        a = others[0]
        return all(zero(a[i]*b[j]-a[j]*b[i]) for b in others[1:]
                   for i in range(dimension) for j in range(i+1, dimension))

    @staticmethod
    def is_concyclic(*points):
        from .ellipse import Circle
        from .entity import unique
        points = unique(Point(point) for point in points)
        if len(points) < 4:
            return not Point.is_collinear(*points) if len(points) == 3 else True
        if any(len(point) != 2 for point in points):
            raise NotImplementedError('Concyclic tests require planar points')
        try:
            circle = Circle(*points[:3])
        except GeometryError:
            return False
        return all(circle.contains(point) for point in points[3:])

    def translate(self, x=0, y=0, z=0):
        if len(self) == 2 and not zero(scalar(z)):
            raise ValueError('A planar point cannot translate along z')
        return self + Point(x, y, *([z] if len(self) == 3 else []))

    def rotate(self, angle, pt=None):
        from .. import sin, cos
        if len(self) != 2:
            raise NotImplementedError('Use a three-dimensional rotation matrix for spatial points')
        origin = Point(0, 0) if pt is None else Point(pt, dim=2)
        x, y = self-origin
        angle = scalar(angle)
        return origin+Point(clean(x*cos(angle)-y*sin(angle)), clean(x*sin(angle)+y*cos(angle)))

    def scale(self, x=1, y=1, z=1, pt=None):
        origin = Point(*([0]*len(self))) if pt is None else Point(pt, dim=len(self))
        factors = tuple(map(scalar, (x, y, z)))
        return origin+Point(*(clean(v*factors[i]) for i, v in enumerate(self-origin)), evaluate=False)

    def transform(self, matrix):
        from .. import Matrix
        if matrix.shape != (len(self)+1, len(self)+1):
            raise ValueError('Expected a homogeneous transformation matrix')
        row = Matrix([list(self)+( [1])])*matrix
        return Point(*(clean(row[i]) for i in range(len(self))), evaluate=False)

    def evalf(self, n=15, **options):
        return Point(*(value.evalf(n, **options) for value in self), evaluate=False)

    n = evalf


class Point2D(Point):
    __slots__ = ()


class Point3D(Point):
    __slots__ = ()

    @property
    def z(self):
        return self.args[2]

    def cross(self, other):
        other = Point(other, dim=3)
        a, b, c = self
        d, e, f = other
        return Point3D(clean(b*f-c*e), clean(c*d-a*f), clean(a*e-b*d), evaluate=False)

    @staticmethod
    def are_coplanar(*points):
        from .plane import Plane
        points = [Point(point, dim=3) for point in points]
        if len(points) < 4 or Point.is_collinear(*points):
            return True
        for b in points[1:]:
            for c in points[2:]:
                if not Point.is_collinear(points[0], b, c):
                    plane = Plane(points[0], b, c)
                    return all(plane.contains(point) for point in points)
