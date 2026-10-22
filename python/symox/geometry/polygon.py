"""Planar polygons, exact boundary intersections and triangle constructions."""
from .entity import GeometrySet, GeometryEntity, GeometryError, construct, scalar, clean, zero, sign, unique
from .point import Point, Point2D
from .line import Line, Segment, LinearEntity, cross2


class Polygon(GeometrySet):
    __slots__ = ()

    def __new__(cls, *vertices, n=0, **kwargs):
        if n:
            # SymPy returns a virtual regular polygon with n sides.
            args = list(vertices)
            if len(args) == 2:
                args.append(n)
            elif len(args) == 3:
                args.insert(2, n)
            return RegularPolygon(*args, **kwargs)
        kwargs.pop('evaluate', None)
        if kwargs:
            raise TypeError('Unknown polygon options')
        vertices = [Point(vertex, dim=2) for vertex in vertices]
        simplified = []
        for vertex in vertices:
            if not simplified or vertex != simplified[-1]:
                simplified.append(vertex)
        if len(simplified) > 1 and simplified[-1] == simplified[0]:
            simplified.pop()
        vertices = simplified
        if len(vertices) > 4096:
            raise NotImplementedError('Polygon construction is limited to 4096 vertices')
        changed = True
        while changed and len(vertices) > 2:
            changed = False
            for i, vertex in enumerate(vertices):
                before, after = vertices[i-1], vertices[(i+1) % len(vertices)]
                if before != after and Point.is_collinear(before, vertex, after):
                    try:
                        redundant = Segment(before, after).contains(vertex)
                    except NotImplementedError:
                        redundant = False
                    if redundant:
                        vertices.pop(i)
                        changed = True
                        break
        if not vertices:
            raise ValueError('A polygon requires vertices')
        if len(vertices) == 1:
            return vertices[0]
        if len(vertices) == 2:
            return Segment(*vertices)
        if cls is Triangle and len(vertices) != 3:
            raise ValueError('A triangle requires exactly three noncollinear vertices')
        return construct(Triangle if cls is Polygon and len(vertices) == 3 else cls, vertices)

    @property
    def vertices(self):
        return list(self.args)

    @property
    def sides(self):
        vertices = self.vertices
        return [Segment(vertex, vertices[(i+1) % len(vertices)]) for i, vertex in enumerate(vertices)]

    @property
    def area(self):
        vertices = self.vertices
        return clean(sum((cross2(vertex, vertices[(i+1) % len(vertices)])
                          for i, vertex in enumerate(vertices)), scalar(0))/2)

    @property
    def perimeter(self):
        return sum((side.length for side in self.sides), scalar(0))

    @property
    def centroid(self):
        vertices = self.vertices
        denominator = 6*self.area
        if zero(denominator):
            raise GeometryError('Zero signed area has no area centroid')
        coordinates = [sum(((vertex[j]+vertices[(i+1) % len(vertices)][j])
                            *cross2(vertex, vertices[(i+1) % len(vertices)])
                            for i, vertex in enumerate(vertices)), scalar(0))/denominator for j in (0, 1)]
        return Point(*map(clean, coordinates))

    @property
    def angles(self):
        from .. import acos, pi
        vertices = self.vertices
        orientation = sign(self.area)
        result = {}
        # SymPy's dictionary starts at the last vertex.
        for i in range(-1, len(vertices)-1):
            point = vertices[i]
            left, right = vertices[i-1]-point, vertices[(i+1) % len(vertices)]-point
            angle = acos(clean(left.dot(right)/(abs(left)*abs(right))))
            if sign(cross2(left, right))*orientation > 0:
                angle = 2*pi-angle
            result[point] = angle
        return result

    def second_moment_of_area(self, point=None):
        I_xx, I_yy, I_xy = scalar(0), scalar(0), scalar(0)
        args = self.vertices
        for i in range(len(args)):
            x1, y1 = args[i-1].args
            x2, y2 = args[i].args
            v = x1*y2 - x2*y1
            I_xx += (y1**2 + y1*y2 + y2**2)*v
            I_yy += (x1**2 + x1*x2 + x2**2)*v
            I_xy += (x1*y2 + 2*x1*y1 + 2*x2*y2 + x2*y1)*v
        A = self.area
        c_x, c_y = self.centroid.args
        I_xx_c = (I_xx/12) - (A*(c_y**2))
        I_yy_c = (I_yy/12) - (A*(c_x**2))
        I_xy_c = (I_xy/24) - (A*(c_x*c_y))
        if point is None:
            return I_xx_c, I_yy_c, I_xy_c
        point = Point(point, dim=2)
        return (I_xx_c + A*((point[1]-c_y)**2), I_yy_c + A*((point[0]-c_x)**2),
                I_xy_c + A*((point[0]-c_x)*(point[1]-c_y)))

    def first_moment_of_area(self, point=None):
        from .. import oo
        if point is not None:
            point = Point(point, dim=2)
            xc, yc = self.centroid.args
        else:
            point = self.centroid
            xc, yc = point.args
        h_poly = self.cut_section(Line(point, slope=0))
        v_poly = self.cut_section(Line(point, slope=oo))
        poly_1 = h_poly[0] if h_poly[0].area <= h_poly[1].area else h_poly[1]
        poly_2 = v_poly[0] if v_poly[0].area <= v_poly[1].area else v_poly[1]
        return (poly_1.centroid.y - yc)*poly_1.area, (poly_2.centroid.x - xc)*poly_2.area

    def polar_second_moment_of_area(self):
        second_moment = self.second_moment_of_area()
        return second_moment[0] + second_moment[1]

    def section_modulus(self, point=None):
        x_c, y_c = self.centroid.args
        if point is None:
            x_min, y_min, x_max, y_max = self.bounds
            y = max(y_c - y_min, y_max - y_c)
            x = max(x_c - x_min, x_max - x_c)
        else:
            point = Point(point, dim=2)
            y = point.y - y_c
            x = point.x - x_c
        second_moment = self.second_moment_of_area()
        return second_moment[0]/y, second_moment[1]/x

    def cut_section(self, line):
        from .. import Symbol
        intersection_points = self.intersection(line)
        if not intersection_points:
            raise ValueError("This line does not intersect the polygon")
        points = list(self.vertices)
        points.append(points[0])
        x, y = Symbol('x'), Symbol('y')
        eq = line.equation(x, y)
        a, b = eq.coeff(x), eq.coeff(y)
        upper_vertices, lower_vertices = [], []
        prev, prev_point = True, None
        for point in points:
            compare = eq.subs({x: point.x, y: point.y})/b if b else eq.subs(x, point.x)/a
            if compare > 0:
                if not prev:
                    new_point = Line(point, prev_point).intersection(line)
                    upper_vertices.append(new_point[0])
                    lower_vertices.append(new_point[0])
                upper_vertices.append(point)
                prev = True
            else:
                if prev and prev_point:
                    new_point = Line(point, prev_point).intersection(line)
                    upper_vertices.append(new_point[0])
                    lower_vertices.append(new_point[0])
                lower_vertices.append(point)
                prev = False
            prev_point = point
        upper_polygon, lower_polygon = None, None
        if upper_vertices and isinstance(Polygon(*upper_vertices), Polygon):
            upper_polygon = Polygon(*upper_vertices)
        if lower_vertices and isinstance(Polygon(*lower_vertices), Polygon):
            lower_polygon = Polygon(*lower_vertices)
        return upper_polygon, lower_polygon

    @staticmethod
    def _is_clockwise(a, b, c):
        from .. import simplify
        ba, ca = b - a, c - a
        res = simplify(ba.x*ca.y - ca.x*ba.y).is_nonpositive
        if res is None:
            raise ValueError("Can't determine orientation")
        return res

    def bisectors(self, prec=None):
        from .line import Ray
        b = {}
        pts = list(self.args)
        pts.append(pts[0])
        if Polygon._is_clockwise(*pts[:3]):
            pts = list(reversed(pts))
        for v, a in self.angles.items():
            i = pts.index(v)
            ray = Ray(pts[i], pts[i + 1]).rotate(a/2, v)
            direction = ray.direction
            ray = Ray(ray.p1, ray.p1 + direction/direction.distance(Point(0, 0)))
            if prec is not None:
                ray = Ray(ray.p1, ray.p2.n(prec))
            b[v] = ray
        return b

    def arbitrary_point(self, parameter='t'):
        from .line import _symbol
        from .. import Piecewise, And
        t = _symbol(parameter)
        if t.name in (f.name for f in self.free_symbols):
            raise ValueError('Symbol %s already appears in object and cannot be used as a parameter.' % t.name)
        sides = []
        perimeter = self.perimeter
        start = scalar(0)
        for side in self.sides:
            fraction = side.length/perimeter
            end = start + fraction
            pt = side.arbitrary_point(t).subs(t, (t - start)/fraction)
            sides.append((pt, And(t >= start, t < end)))
            start = end
        return _GeometryPiecewise(sides)

    def plot_interval(self, parameter='t'):
        from .. import Symbol
        return [Symbol(parameter, real=True), 0, 1]

    def reflect(self, line):
        return self.func(*(vertex.reflect(line) for vertex in self.vertices))

    def is_convex(self):
        vertices = self.vertices
        signs = {sign(cross2(vertices[i]-vertices[i-1], vertices[(i+1) % len(vertices)]-vertices[i]))
                 for i in range(len(vertices))}
        return not ({-1, 1} <= signs)

    def contains(self, other):
        if isinstance(other, Segment):
            return any(side.contains(other) for side in self.sides)
        if isinstance(other, Polygon):
            return self.equals(other)
        try:
            point = Point(other, dim=2)
        except (TypeError, ValueError):
            return False
        return any(side.contains(point) for side in self.sides)

    def encloses_point(self, p):
        p = Point(p, dim=2)
        if self.contains(p):
            return False
        # Winding number avoids floating rays and retains exact endpoints.
        winding = 0
        for side in self.sides:
            a, b = side.p1-p, side.p2-p
            ay, by = sign(a.y), sign(b.y)
            if ay <= 0 < by and sign(cross2(a, b)) > 0:
                winding += 1
            elif by <= 0 < ay and sign(cross2(a, b)) < 0:
                winding -= 1
        return winding != 0

    def encloses(self, o):
        if isinstance(o, Point):
            return self.encloses_point(o)
        if isinstance(o, (Segment, Polygon)):
            points = o.points if isinstance(o, Segment) else o.vertices
            return all(self.encloses_point(point) for point in points) and not self.intersection(o)
        raise NotImplementedError('Polygon enclosure currently accepts points, segments and polygons')

    def distance(self, o):
        from .. import Min
        if not isinstance(o, GeometryEntity):
            o = Point(o, dim=2)
        if isinstance(o, Point):
            return Min(*(side.distance(o) for side in self.sides))
        if isinstance(o, Polygon):
            return Min(*(side.distance(edge) for side in self.sides for edge in o.sides))
        if isinstance(o, LinearEntity):
            return Min(*(side.distance(o) for side in self.sides))
        raise NotImplementedError('This polygon distance is not implemented')

    def intersection(self, o):
        if isinstance(o, Point):
            return [o] if self.contains(o) else []
        edges = o.sides if isinstance(o, Polygon) else [o]
        result = unique(item for side in self.sides for edge in edges for item in side.intersection(edge))
        segments = [item for item in result if isinstance(item, Segment)]
        return [item for item in result if not isinstance(item, Point) or not any(side.contains(item) for side in segments)]

    def equals(self, o):
        if not isinstance(o, Polygon) or len(self.vertices) != len(o.vertices):
            return False
        a, b = self.vertices, o.vertices
        for oriented in (b, list(reversed(b))):
            for index in range(len(a)):
                if all(left.equals(oriented[(index+i) % len(a)]) for i, left in enumerate(a)):
                    return True
        return False

    def is_similar(self, other):
        if not isinstance(other, Polygon) or len(self.vertices) != len(other.vertices):
            return False
        first, second = self.vertices, other.vertices
        if len(first) > 256:
            raise NotImplementedError('General polygon similarity is limited to 256 vertices')
        original = [point-first[0] for point in first]
        for orientation in (second, list(reversed(second))):
            for start in range(len(first)):
                target = [orientation[(start+i) % len(first)]-orientation[start] for i in range(len(first))]
                for reflection in (1, -1):
                    source = [Point(point.x, reflection*point.y) for point in original]
                    a, b = source[1], target[1]
                    denominator = a.dot(a)
                    if zero(denominator):
                        continue
                    cosine, sine = clean(a.dot(b)/denominator), clean(cross2(a,b)/denominator)
                    if all(Point(clean(cosine*p.x-sine*p.y), clean(sine*p.x+cosine*p.y)).equals(q)
                           for p, q in zip(source,target)):
                        return True
        return False


class _GeometryPiecewise:
    """Piecewise geometric value, mirroring SymPy's Piecewise of points."""
    __slots__ = ('_args',)

    def __init__(self, pairs):
        object.__setattr__(self, '_args', tuple((expr, cond) for expr, cond in pairs))

    @property
    def args(self):
        return self._args

    @property
    def free_symbols(self):
        result = set()
        for expr, cond in self.args:
            result |= expr.free_symbols | getattr(cond, 'free_symbols', set())
        return result

    def _rebuild(self, method, *args, **kwargs):
        from .. import S
        pairs = []
        for expr, cond in self.args:
            cond = getattr(cond, method)(*args, **kwargs) if hasattr(cond, method) else cond
            if cond == S.false or cond is False:
                continue
            expr = getattr(expr, method)(*args, **kwargs)
            if (cond == S.true or cond is True) and not pairs:
                return expr
            pairs.append((expr, cond))
            if cond == S.true or cond is True:
                break
        if not pairs:
            from .. import nan
            return nan
        return _GeometryPiecewise(pairs)

    def subs(self, *args, **kwargs):
        return self._rebuild('subs', *args, **kwargs)

    def xreplace(self, rule):
        return self._rebuild('xreplace', rule)

    def __eq__(self, other):
        return isinstance(other, _GeometryPiecewise) and self.args == other.args

    def __hash__(self):
        return hash(self.args)

    def __repr__(self):
        return 'Piecewise(%s)' % ', '.join('(%s, %s)' % pair for pair in self.args)

    __str__ = __repr__


class Triangle(Polygon):
    __slots__ = ()

    def __new__(cls, *args, **kwargs):
        if len(args) != 3:
            from .. import simplify
            if 'sss' in kwargs:
                return _sss(*[simplify(a) for a in kwargs['sss']])
            if 'asa' in kwargs:
                return _asa(*[simplify(a) for a in kwargs['asa']])
            if 'sas' in kwargs:
                return _sas(*[simplify(a) for a in kwargs['sas']])
            raise GeometryError("Triangle instantiates with three points or a valid keyword.")
        return Polygon.__new__(cls, *args, **kwargs)

    @property
    def medians(self):
        return {vertex:Segment(vertex, Segment(self.vertices[(i+1) % 3], self.vertices[(i+2) % 3]).midpoint)
                for i, vertex in enumerate(self.vertices)}

    @property
    def altitudes(self):
        return {vertex:Line(self.vertices[(i+1) % 3], self.vertices[(i+2) % 3]).perpendicular_segment(vertex)
                for i, vertex in enumerate(self.vertices)}

    @property
    def orthocenter(self):
        vertices = self.vertices
        lines = [Line(vertices[(i+1) % 3], vertices[(i+2) % 3]).perpendicular_line(vertices[i]) for i in (0, 1)]
        return lines[0].intersection(lines[1])[0]

    @property
    def circumcenter(self):
        return self.circumcircle.center

    @property
    def circumradius(self):
        return self.circumcircle.radius

    @property
    def circumcircle(self):
        from .ellipse import Circle
        return Circle(*self.vertices)

    @property
    def incenter(self):
        vertices = self.vertices
        weights = [vertices[(i+1) % 3].distance(vertices[(i+2) % 3]) for i in range(3)]
        return sum((vertex*weight for vertex, weight in zip(vertices, weights)), Point(0, 0))/sum(weights, scalar(0))

    @property
    def inradius(self):
        from .. import Abs
        return clean(2*Abs(self.area)/self.perimeter)

    @property
    def incircle(self):
        from .ellipse import Circle
        return Circle(self.incenter, self.inradius)

    def bisectors(self):
        s = [Line(side) for side in self.sides]
        v = self.vertices
        c = self.incenter
        l1 = Segment(v[0], Line(v[0], c).intersection(s[1])[0])
        l2 = Segment(v[1], Line(v[1], c).intersection(s[2])[0])
        l3 = Segment(v[2], Line(v[2], c).intersection(s[0])[0])
        return {v[0]: l1, v[1]: l2, v[2]: l3}

    def is_right(self):
        return any(zero((self.vertices[(i+1) % 3]-vertex).dot(self.vertices[(i+2) % 3]-vertex))
                   for i, vertex in enumerate(self.vertices))

    def is_equilateral(self):
        lengths = [side.length**2 for side in self.sides]
        return zero(lengths[0]-lengths[1]) and zero(lengths[1]-lengths[2])

    def is_isosceles(self):
        lengths = [side.length**2 for side in self.sides]
        return any(zero(lengths[i]-lengths[j]) for i in range(3) for j in range(i+1, 3))

    def is_scalene(self):
        return not self.is_isosceles()

    def is_similar(self, other):
        from itertools import permutations
        if not isinstance(other, Triangle):
            return False
        left = [side.length**2 for side in self.sides]
        return any(all(zero(left[i]*right[0]-right[i]*left[0]) for i in (1, 2))
                   for right in permutations([side.length**2 for side in other.sides]))


class RegularPolygon(Polygon):
    __slots__ = ()

    def __new__(cls, c, r, n, rot=0, **kwargs):
        from .. import pi
        kwargs.pop('evaluate', None)
        center, radius, n, rot = Point(c, dim=2), scalar(r), scalar(n), scalar(rot)
        if n.is_integer is not True or n.is_number is not True or int(n) < 3:
            raise GeometryError('n must be an integer >= 3, not %s' % n)
        if int(n) > 4096:
            raise NotImplementedError('Regular polygons are limited to 4096 vertices')
        # SymPy accepts any radius; a negative one records a reversed traversal.
        if radius.is_extended_real is False:
            raise GeometryError('r must be a real expression, not %s' % radius)
        if rot.is_number:
            rot = rot % (2*pi/n)
        return construct(cls, (center, radius, n, rot))

    @property
    def center(self):
        return self.args[0]

    @property
    def radius(self):
        return self.args[1]

    @property
    def n(self):
        return self.args[2]

    @property
    def rotation(self):
        return self.args[3]

    @property
    def vertices(self):
        from .. import pi, sin, cos
        from .. import Abs
        c, r, rot, v = self.center, Abs(self.radius), self.rotation, 2*pi/self.n
        return [Point(c.x + r*cos(k*v + rot), c.y + r*sin(k*v + rot)) for k in range(int(self.n))]

    @property
    def inradius(self):
        return self.apothem

    @property
    def circumradius(self):
        return self.radius

    @property
    def length(self):
        from .. import sin, pi
        return self.radius*2*sin(pi/self.n)

    def spin(self, angle):
        """Rotate in place about the center (SymPy mutates the polygon too)."""
        c, r, n, rot = self.args
        object.__setattr__(self, '_args', (c, r, n, rot + scalar(angle)))

    def reflect(self, line):
        from .line import Ray
        c, r, n, rot = self.args
        v = self.vertices[0]
        d = v - c
        cc = c.reflect(line)
        vv = v.reflect(line)
        dd = vv - cc
        ang = Ray((0, 0), dd).closing_angle(Ray((0, 0), d))
        return self.func(cc, -r, n, rot + ang)

    @property
    def centroid(self):
        return self.center

    @property
    def apothem(self):
        from .. import cos, pi
        return self.radius*cos(pi/self.n)

    @property
    def interior_angle(self):
        from .. import pi
        return (self.n-2)*pi/self.n

    @property
    def exterior_angle(self):
        from .. import pi
        return 2*pi/self.n

    @property
    def circumcircle(self):
        from .ellipse import Circle
        return Circle(self.center, self.radius)

    @property
    def incircle(self):
        from .ellipse import Circle
        return Circle(self.center, self.apothem)

    def rotate(self, angle, pt=None):
        return RegularPolygon(self.center.rotate(angle, pt), self.radius, self.n, self.rotation+scalar(angle))

    def scale(self, x=1, y=1, pt=None):
        x, y = scalar(x), scalar(y)
        if zero(x-y) and x.is_positive is True:
            return RegularPolygon(self.center.scale(x, y, pt=pt), self.radius*x, self.n, self.rotation)
        return Polygon(*(vertex.scale(x, y, pt=pt) for vertex in self.vertices))

def _asa(d1, l, d2):
    from .. import tan
    from .util import rad
    xy = Line((0, 0), slope=tan(rad(d1))).intersection(Line((l, 0), slope=tan(rad(180 - d2))))[0]
    return Triangle((0, 0), (l, 0), xy)


def _sss(l1, l2, l3):
    from .ellipse import Circle
    c1 = Circle((0, 0), l3)
    c2 = Circle((l1, 0), l2)
    inter = [a for a in c1.intersection(c2) if a.y.is_nonnegative]
    if not inter:
        return None
    return Triangle((0, 0), (l1, 0), inter[0])


def _sas(l1, d, l2):
    from .. import sin, cos
    from .util import rad
    return Triangle(Point(0, 0), Point(l2, 0), Point(cos(rad(d))*l1, sin(rad(d))*l1))


# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.geometry.util', 'deg rad')
