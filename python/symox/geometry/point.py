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
                # SymPy keeps an existing point's coordinates unevaluated.
                evaluate = False
                target = dim if dim is not None else 2 if cls is Point2D else 3 if cls is Point3D else len(args[0])
                if len(args[0]) == target:
                    return args[0]
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
        if dim < 2:
            raise ValueError('At least two coordinates are required')
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
        return construct((Point2D if dim == 2 else Point3D if dim == 3 else Point) if cls is Point else cls, values)

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

    def dot(self, p):
        p = Point(p, dim=len(self))
        return clean(sum((a*b for a, b in zip(self, p)), scalar(0)))

    def distance(self, other):
        if isinstance(other, GeometryEntity) and not isinstance(other, Point):
            return other.distance(self)
        return abs(self - Point(other))

    def taxicab_distance(self, p):
        from .. import Abs
        return sum((Abs(v) for v in self-Point(p)), scalar(0))

    def canberra_distance(self, p):
        from .. import Abs
        p = Point(p)
        dim = max(len(self), len(p))
        s, p = Point(self, dim=dim), Point(p, dim=dim)
        if s.is_zero and p.is_zero:
            raise ValueError('Cannot project to the zero vector.')
        return sum((Abs(a-b)/(Abs(a)+Abs(b)) for a, b in zip(s, p)), scalar(0))

    @property
    def is_zero(self):
        values = [value.is_zero for value in self]
        if all(values):
            return True
        if any(value is False for value in values):
            return False
        return None

    @property
    def is_nonzero(self):
        value = self.is_zero
        return None if value is None else not value

    @property
    def origin(self):
        return Point(*([0]*len(self)), evaluate=False)

    @staticmethod
    def project(a, b):
        a, b = Point(a), Point(b)
        dim = max(len(a), len(b))
        a, b = Point(a, dim=dim), Point(b, dim=dim)
        if b.is_zero:
            raise ValueError('Cannot project to the zero vector.')
        return b*(a.dot(b)/b.dot(b))

    @property
    def orthogonal_direction(self):
        dim = len(self)
        if self[0].is_zero:
            return Point(*([1]+(dim-1)*[0]))
        if self[1].is_zero:
            return Point(*([0, 1]+(dim-2)*[0]))
        return Point(*([-self[1], self[0]]+(dim-2)*[0]))

    def reflect(self, line):
        from .line import LinearEntity
        if not isinstance(line, LinearEntity) or line.ambient_dimension != 2:
            raise AttributeError('Reflection requires a planar line')
        a, b = line.p1, line.direction
        normal = Point(-b.y, b.x)
        planar = Point(self.x, self.y)
        image = planar-normal*clean(2*(planar-a).dot(normal)/normal.dot(normal))
        return Point(*(tuple(image)+tuple(self.args[2:])), evaluate=False)

    def midpoint(self, p):
        return (self + Point(p))/2

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
    def affine_rank(*args):
        """Dimension of the smallest affine space containing the points."""
        from .. import Matrix
        if not args:
            return -1
        points = [Point(point) for point in args]
        dimension = max(map(len, points))
        points = [Point(point, dim=dimension) for point in points]
        origin = points[0]
        rows = [list(point-origin) for point in points[1:]]
        if not rows:
            return 0
        return Matrix(rows).rank(iszerofunc=lambda value:
            abs(value.n(2)) < 1e-12 if value.is_number else value.is_zero)

    def is_collinear(self, *args):
        # Reference instance method: ``self`` is one of the tested points
        # (``Point.is_collinear(p, q, r)`` still works unbound).
        points = [Point(point) for point in (self,)+args]
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

    def is_concyclic(self, *args):
        # Reference algorithm: the points are concyclic when they are
        # coplanar and the system |p_i - c|**2 == r**2 is consistent, i.e.
        # the augmented column of squared norms carries no pivot.
        from .. import Matrix
        from .entity import unique
        points = [Point(point) for point in (self,)+args]
        dimension = max(map(len, points))
        points = unique(Point(point, dim=dimension) for point in points)
        if Point.affine_rank(*points) > 2:
            return False
        origin = points[0]
        rows = []
        for point in points:
            delta = point-origin
            rows.append(list(delta)+[delta.dot(delta)])
        return dimension not in Matrix(rows).rref()[1]

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

    def direction_ratio(self, point):
        return [point.x-self.x, point.y-self.y, point.z-self.z]

    def direction_cosine(self, point):
        from .. import sqrt
        a = self.direction_ratio(point)
        b = sqrt(sum((i**2 for i in a), scalar(0)))
        return [(point.x-self.x)/b, (point.y-self.y)/b, (point.z-self.z)/b]

    @staticmethod
    def are_collinear(*points):
        return Point.is_collinear(*points)

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
