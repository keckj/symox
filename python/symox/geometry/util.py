"""Exact Euclidean utilities; ambiguous orderings fail explicitly."""
from functools import cmp_to_key
from .entity import GeometryEntity, GeometryError, scalar, clean, zero, sign, compare, unique
from .point import Point
from .line import Segment, cross2


def rad(degrees):
    from .. import pi
    return scalar(degrees)*pi/180


def deg(radians):
    from .. import pi
    return 180*scalar(radians)/pi


def intersection(*entities, pairwise=False, **kwargs):
    if kwargs:
        raise TypeError('Unknown intersection options')
    if len(entities) < 2:
        return []
    entities = [value if isinstance(value, GeometryEntity) else Point(value) for value in entities]
    if pairwise:
        return unique(item for i, left in enumerate(entities) for right in entities[i+1:]
                      for item in left.intersection(right))
    result = [entities[0]]
    for entity in entities[1:]:
        result = unique(item for current in result for item in (
            [current] if isinstance(current, Point) and entity.contains(current)
            else [] if isinstance(current, Point) else current.intersection(entity)))
    return result


def are_similar(first, second):
    if hasattr(first, 'is_similar'):
        return first.is_similar(second)
    if isinstance(first, Point) and isinstance(second, Point):
        return True
    raise NotImplementedError('Similarity is not implemented for these geometry types')


def centroid(*entities):
    from .polygon import Polygon
    from .. import Abs, S
    if not entities:
        return None
    if all(isinstance(value, Point) for value in entities):
        return sum(entities[1:], entities[0])/len(entities)
    if all(isinstance(value, Segment) for value in entities):
        weighted = [(value.midpoint, value.length) for value in entities]
    elif all(isinstance(value, Polygon) for value in entities):
        weighted = [(value.centroid, Abs(value.area)) for value in entities]
    else:
        raise ValueError('Centroid requires points, segments, or polygons of one common kind')
    total = sum((weight for _, weight in weighted), S.Zero)
    if zero(total):
        raise GeometryError('Cannot weight a centroid by zero total measure')
    start = Point(*([0]*len(weighted[0][0])))
    return sum((point*weight for point, weight in weighted), start)/total


def convex_hull(*entities, polygon=True):
    from .polygon import Polygon
    points = []
    for entity in entities:
        if isinstance(entity, Polygon):
            points.extend(entity.vertices)
        elif isinstance(entity, Segment):
            points.extend(entity.points)
        else:
            points.append(Point(entity, dim=2))
    points = unique(points)
    if not points:
        raise ValueError('A convex hull requires at least one point')
    def order(a, b):
        x = compare(a.x, b.x)
        return x if x else compare(a.y, b.y)
    points.sort(key=cmp_to_key(order))
    if len(points) == 1:
        result = points[0]
    else:
        def half(items):
            result = []
            for item in items:
                while len(result) >= 2 and sign(cross2(result[-1]-result[-2], item-result[-1])) <= 0:
                    result.pop()
                result.append(item)
            return result
        vertices = half(points)[:-1]+half(reversed(points))[:-1]
        result = Polygon(*vertices)
    if polygon:
        return result
    if isinstance(result, Polygon):
        # SymPy returns the upper chain left-to-right and the lower chain
        # left-to-right, with both endpoint vertices included.
        lower, upper = half(points), list(reversed(half(reversed(points))))
        return upper, lower
    return result, None


def _distance_pairs(points, farthest):
    from .ellipse import _sort_points
    points = _sort_points(unique(Point(point) for point in points))
    if len(points) < 2:
        raise ValueError('At least two distinct points are required')
    if len(points) > 1024:
        raise NotImplementedError('Distance-pair search is limited to 1024 points')
    best = None
    result = set()
    for i, left in enumerate(points):
        for right in points[i+1:]:
            delta = left-right
            distance = delta.dot(delta)
            relation = 1 if best is None else compare(distance, best)
            if best is None or relation > 0 and farthest or relation < 0 and not farthest:
                best, result = distance, {(left, right)}
            elif relation == 0:
                result.add((left, right))
    return result


def closest_points(*points):
    return _distance_pairs(points, False)


def farthest_points(*points):
    return _distance_pairs(points, True)


def idiff(eq, y, x, n=1):
    from .. import diff, cancel, Function, Dummy
    from ..solvers._differential import equation
    expression, x = equation(eq), scalar(x)
    if not x.is_Symbol or not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError('Implicit differentiation requires a symbol and positive integer order')
    if n > 256:
        raise NotImplementedError('Implicit differentiation is limited to order 256')
    dependent = tuple(y) if isinstance(y, (tuple, list)) else (y,)
    dependent = tuple(map(scalar, dependent))
    if any(not value.is_Symbol for value in dependent) or x in dependent:
        raise ValueError('Dependent coordinates must be symbols distinct from x')
    first, others = dependent[0], dependent[1:]
    # Other dependent coordinates remain native held derivatives.
    # Dummy display names repeat; their unique identities must distinguish
    # different dependent coordinates before they become function heads.
    from ..solvers._differential import nodes
    occupied = {getattr(node.func,'__name__','') for node in nodes(expression)
                if node.is_UndefinedFunction}
    replacements = {}
    for value in others:
        while True:
            name = 'implicit_'+str(Dummy('coordinate').dummy_index)
            if name not in occupied:
                break
        occupied.add(name)
        replacements[value] = Function(name)(x)
    prepared = expression.xreplace(replacements)
    denominator = diff(prepared, first)
    if zero(denominator):
        raise GeometryError('The equation does not locally determine the requested coordinate')
    slope = cancel(-diff(prepared, x)/denominator)
    value = slope
    for _ in range(n-1):
        value = cancel(diff(value, x)+diff(value, first)*slope)
    return value.xreplace({applied:original for original, applied in replacements.items()})
