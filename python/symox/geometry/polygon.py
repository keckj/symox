"""Planar polygons, exact boundary intersections and triangle constructions."""
from .entity import GeometryEntity, GeometryError, construct, scalar, clean, zero, sign, unique
from .point import Point, Point2D
from .line import Line, Segment, LinearEntity, cross2


class Polygon(GeometryEntity):
    __slots__ = ()

    def __new__(cls, *vertices, **kwargs):
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
        for i, point in enumerate(vertices):
            left, right = vertices[i-1]-point, vertices[(i+1) % len(vertices)]-point
            angle = acos(clean(left.dot(right)/(abs(left)*abs(right))))
            if sign(cross2(left, right))*orientation > 0:
                angle = 2*pi-angle
            result[point] = angle
        return result

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

    def encloses_point(self, point):
        point = Point(point, dim=2)
        if self.contains(point):
            return False
        # Winding number avoids floating rays and retains exact endpoints.
        winding = 0
        for side in self.sides:
            a, b = side.p1-point, side.p2-point
            ay, by = sign(a.y), sign(b.y)
            if ay <= 0 < by and sign(cross2(a, b)) > 0:
                winding += 1
            elif by <= 0 < ay and sign(cross2(a, b)) < 0:
                winding -= 1
        return winding != 0

    def encloses(self, other):
        if isinstance(other, Point):
            return self.encloses_point(other)
        if isinstance(other, (Segment, Polygon)):
            points = other.points if isinstance(other, Segment) else other.vertices
            return all(self.encloses_point(point) for point in points) and not self.intersection(other)
        raise NotImplementedError('Polygon enclosure currently accepts points, segments and polygons')

    def distance(self, other):
        from .. import Min
        if not isinstance(other, GeometryEntity):
            other = Point(other, dim=2)
        if isinstance(other, Point):
            return Min(*(side.distance(other) for side in self.sides))
        if isinstance(other, Polygon):
            return Min(*(side.distance(edge) for side in self.sides for edge in other.sides))
        if isinstance(other, LinearEntity):
            return Min(*(side.distance(other) for side in self.sides))
        raise NotImplementedError('This polygon distance is not implemented')

    def intersection(self, other):
        if isinstance(other, Point):
            return [other] if self.contains(other) else []
        edges = other.sides if isinstance(other, Polygon) else [other]
        result = unique(item for side in self.sides for edge in edges for item in side.intersection(edge))
        segments = [item for item in result if isinstance(item, Segment)]
        return [item for item in result if not isinstance(item, Point) or not any(side.contains(item) for side in segments)]

    def equals(self, other):
        if not isinstance(other, Polygon) or len(self.vertices) != len(other.vertices):
            return False
        a, b = self.vertices, other.vertices
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


class Triangle(Polygon):
    __slots__ = ()

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

    @property
    def bisectors(self):
        result = {}
        for i, vertex in enumerate(self.vertices):
            a, b = self.vertices[(i+1) % 3], self.vertices[(i+2) % 3]
            la, lb = vertex.distance(a), vertex.distance(b)
            result[vertex] = Segment(vertex, (a*lb+b*la)/(la+lb))
        return result

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

    def __new__(cls, center, radius, n, rot=0, **kwargs):
        center, radius, n, rot = Point(center, dim=2), scalar(radius), scalar(n), scalar(rot)
        if n.is_integer is not True or n.is_number is not True or int(n) < 3:
            raise ValueError('A regular polygon requires an integer vertex count of at least three')
        if int(n) > 4096:
            raise NotImplementedError('Regular polygons are limited to 4096 vertices')
        if radius.is_positive is False or zero(radius):
            raise ValueError('Regular polygon radius must be positive')
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
        return [self.center+Point(self.radius*cos(self.rotation+2*pi*i/self.n),
                                  self.radius*sin(self.rotation+2*pi*i/self.n)) for i in range(int(self.n))]

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
