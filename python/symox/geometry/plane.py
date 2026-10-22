"""Exact affine planes and spatial intersections."""
from .entity import GeometryEntity, GeometryError, construct, scalar, clean, zero
from .point import Point, Point3D
from .line import LinearEntity, Line


def _circle_parameter(equations, u):
    """Solve the unit-circle parametrisation one component at a time.

    Every candidate from a single trigonometric component is checked exactly
    against all components; an undecidable check raises NotImplementedError.
    """
    from .. import solve, count_ops, simplify
    pending = list(equations)
    for equation in sorted((e for e in pending if u in e.free_symbols), key=count_ops):
        try:
            candidates = solve(equation, u)
        except NotImplementedError:
            continue
        undecided = False
        for candidate in candidates:
            residues = [simplify(e.xreplace({u: candidate})) for e in pending]
            if all(r.is_zero is True for r in residues):
                return [{u: candidate}]
            if any(r.is_zero is None for r in residues):
                undecided = True
        if not undecided:
            return []
    raise NotImplementedError('The circle parameter of this point is undecidable')


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

    @property
    def ambient_dimension(self):
        return 3

    def parallel_plane(self, pt):
        return Plane(pt, normal_vector=self.normal_vector)

    def perpendicular_plane(self, *pts):
        from .line import Line3D
        if len(pts) > 2:
            raise ValueError('No more than 2 pts should be provided.')
        pts = list(pts)
        if not pts:
            pts.append(self.p1)
        nx, ny, nz = self.normal_vector
        if len(pts) == 1:
            direction = (0, 1, 0) if nx == ny == 0 else (0, 0, 1)
            pts.append(Point(pts[0], dim=3)+Point3D(*direction))
        p1, p2 = [Point(i, dim=3) for i in pts]
        l = Line3D(p1, p2)
        n = Line3D(p1, direction_ratio=self.normal_vector)
        if n.contains(l):
            p3 = Point3D(0, 1, 0) if nx == ny == 0 else Point3D(0, 0, 1)
            if l.contains(p3):
                p3 = p3*2
        else:
            p3 = p1+Point3D(*self.normal_vector)
        return Plane(p1, p2, p3)

    @staticmethod
    def are_concurrent(*planes):
        from .entity import unique
        planes = unique(planes)
        for plane in planes:
            if not isinstance(plane, Plane):
                raise ValueError('All objects should be Planes but got %s' % type(plane).__name__)
        if len(planes) < 2:
            return False
        first = planes.pop(0)
        sol = first.intersection(planes[0])
        if sol == []:
            return False
        line = sol[0]
        for plane in planes[1:]:
            l = first.intersection(plane)
            if not l or not line.contains(l[0]):
                return False
        return True

    def is_coplanar(self, o):
        from .. import cancel, Symbol
        from .line import LinearEntity3D
        if isinstance(o, Plane):
            x, y, z = (Symbol(n, real=True) for n in 'xyz')
            return not cancel(self.equation(x, y, z)/o.equation(x, y, z)).has(x, y, z)
        if isinstance(o, Point3D):
            return self.contains(o)
        if isinstance(o, LinearEntity3D):
            # Reference quirk: iterating a plane raises TypeError.
            raise TypeError("'Plane' object is not iterable")
        if isinstance(o, GeometryEntity):
            return all(i == 0 for i in self.normal_vector[:2])
        return None

    def arbitrary_point(self, u=None, v=None):
        from .. import Symbol, S, sqrt, cos, sin
        circle = v is None
        name = lambda value, default: Symbol(value or default, real=True) if value is None or isinstance(value, str) else scalar(value)
        if circle:
            u = name(u, 't')
        else:
            u, v = name(u, 'u'), name(v, 'v')
        x, y, z = self.normal_vector
        a, b, c = self.p1.args
        if x.is_zero and y.is_zero:
            x1, y1, z1 = S.One, S.Zero, S.Zero
        else:
            x1, y1, z1 = -y, x, S.Zero
        x2, y2, z2 = y*z1-z*y1, z*x1-x*z1, x*y1-y*x1
        if circle:
            n1, n2 = sqrt(x1**2+y1**2+z1**2), sqrt(x2**2+y2**2+z2**2)
            x1, y1, z1 = (w/n1 for w in (x1, y1, z1))
            x2, y2, z2 = (w/n2 for w in (x2, y2, z2))
            return Point3D(a+x1*cos(u)+x2*sin(u), b+y1*cos(u)+y2*sin(u), c+z1*cos(u)+z2*sin(u), evaluate=False)
        return Point3D(a+x1*u+x2*v, b+y1*u+y2*v, c+z1*u+z2*v, evaluate=False)

    def parameter_value(self, other, u, v=None):
        from .. import Symbol, solve
        if not isinstance(other, GeometryEntity):
            other = Point(other, dim=3)
        if not isinstance(other, Point):
            raise ValueError('other must be a point')
        if other == self.p1:
            return other
        if isinstance(u, Symbol) and v is None:
            delta = self.arbitrary_point(u)-self.p1
            eq = delta-(other-self.p1).unit
            try:
                sol = solve(list(eq), u, dict=True)
            except NotImplementedError:
                sol = _circle_parameter(list(eq), u)
        elif isinstance(u, Symbol) and isinstance(v, Symbol):
            pt = self.arbitrary_point(u, v)
            sol = solve(list(pt-other), (u, v), dict=True)
        else:
            raise ValueError('expecting 1 or 2 symbols')
        if not sol:
            raise ValueError('Given point is not on Plane')
        return sol[0]

    def random_point(self, seed=None):
        from .. import Rational, Dummy
        from .line import _rng
        rng = _rng(seed)
        u, v = Dummy('u', real=True), Dummy('v', real=True)
        params = {u: 2*Rational(rng.gauss(0, 1))-1, v: 2*Rational(rng.gauss(0, 1))-1}
        return Point3D(*(clean(value.xreplace(params)) for value in self.arbitrary_point(u, v)))

    def translate(self, x=0, y=0, z=0):
        return Plane(self.p1.translate(x, y, z), normal_vector=self.normal_vector)

    def scale(self, x=1, y=1, z=1, pt=None):
        factors = list(map(scalar, (x, y, z)))
        if any(zero(factor) for factor in factors):
            raise NotImplementedError('Singular plane scaling is not supported')
        return Plane(self.p1.scale(x, y, z, pt=pt), normal_vector=[v/f for v, f in zip(self.normal_vector, factors)])
