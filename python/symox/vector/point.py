"""Points connected by displacement vectors."""
from ._basis import Vector


class Point:
    def __init__(self,name,position=0,parent_point=None):
        if not isinstance(name,str):
            raise TypeError('point name must be a string')
        if position == 0:
            position = Vector.zero
        if not isinstance(position,Vector):
            raise TypeError('position must be a vector')
        if parent_point is not None and not isinstance(parent_point,Point):
            raise TypeError('parent_point must be a Point')
        if parent_point is None and position:
            raise ValueError('a position vector needs a parent point')
        self._name,self._position,self._parent = name,position,parent_point

    def __repr__(self):
        return self._name

    def __hash__(self):
        return hash((self._name,self._parent,frozenset(self._position.components.items())))

    def __eq__(self,other):
        return isinstance(other,Point) and (self._name,self._parent,self._position) == (other._name,other._parent,other._position)

    def locate_new(self,name,position):
        return Point(name,position,self)

    def position_wrt(self,other):
        from .coordsysrect import CoordSys3D
        if isinstance(other,CoordSys3D):
            other = other.origin
        if not isinstance(other,Point):
            raise TypeError('expected a point or coordinate system')
        ancestors = {}
        node,position = self,Vector.zero
        while node is not None:
            ancestors[node] = position
            position += node._position
            node = node._parent
        node,position = other,Vector.zero
        while node is not None:
            if node in ancestors:
                return ancestors[node]-position
            position += node._position
            node = node._parent
        raise ValueError('points have no connecting path')

    def express_coordinates(self,coordinate_system):
        return tuple(self.position_wrt(coordinate_system).to_matrix(coordinate_system))
