"""Exact Euclidean geometry over native symox scalar expressions."""
from .entity import GeometryError, GeometryEntity
from .point import Point, Point2D, Point3D
from .line import Line, Ray, Segment, Line2D, Ray2D, Segment2D, Line3D, Ray3D, Segment3D
from .plane import Plane
from .ellipse import Ellipse, Circle
from .polygon import Polygon, Triangle, RegularPolygon
from .curve import Curve
from .parabola import Parabola
from .util import rad, deg, are_similar, centroid, convex_hull, idiff, intersection, closest_points, farthest_points

__all__ = ['Point', 'Point2D', 'Point3D', 'Line', 'Ray', 'Segment', 'Line2D',
           'Ray2D', 'Segment2D', 'Line3D', 'Ray3D', 'Segment3D', 'Plane', 'Ellipse',
           'Circle', 'Polygon', 'Triangle', 'RegularPolygon', 'Curve', 'Parabola',
           'rad', 'deg', 'are_similar', 'centroid', 'convex_hull', 'idiff',
           'intersection', 'closest_points', 'farthest_points', 'GeometryError']
