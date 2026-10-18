"""Passive direction cosine matrices for axis and angle conventions."""
from .. import Matrix, Integer, sin, cos, sympify, simplify
from ._basis import Vector


def _rot(axis,angle):
    c,s = cos(angle),sin(angle)
    if axis == 1:
        return Matrix([[1,0,0],[0,c,s],[0,-s,c]])
    if axis == 2:
        return Matrix([[c,0,-s],[0,1,0],[s,0,c]])
    return Matrix([[c,s,0],[-s,c,0],[0,0,1]])


class AxisOrienter:
    def __init__(self,angle,axis):
        if not isinstance(axis,Vector) or not axis:
            raise ValueError('rotation axis must be a nonzero vector')
        self.angle,self.axis = sympify(angle),axis

    def rotation_matrix(self,system=None):
        if system is None:
            system = next(iter(self.axis.components)).system
        unit = tuple(self.axis.normalize().to_matrix(system))
        x,y,z = unit
        cross = Matrix([[0,-z,y],[z,0,-x],[-y,x,0]])
        column = Matrix(unit)
        return (cos(self.angle)*Matrix.eye(3)+(1-cos(self.angle))*column*column.T-sin(self.angle)*cross).applyfunc(simplify)


class _ThreeAngle:
    body = True

    def __init__(self,angle1,angle2,angle3,rot_order):
        order = str(rot_order).upper().translate(str.maketrans('XYZ','123'))
        if len(order) != 3 or any(c not in '123' for c in order) or order[0] == order[1] or order[1] == order[2]:
            raise ValueError('rotation order must be a valid three-axis Euler or Tait-Bryan sequence')
        self.angle1,self.angle2,self.angle3 = map(sympify,(angle1,angle2,angle3))
        self.rot_order = order

    def rotation_matrix(self):
        matrices = [_rot(int(i),a) for i,a in zip(self.rot_order,(self.angle1,self.angle2,self.angle3))]
        a,b,c = matrices
        return c*b*a if self.body else a*b*c


class BodyOrienter(_ThreeAngle):
    pass


class SpaceOrienter(_ThreeAngle):
    body = False


class QuaternionOrienter:
    def __init__(self,q0,q1,q2,q3):
        self.q0,self.q1,self.q2,self.q3 = map(sympify,(q0,q1,q2,q3))

    def rotation_matrix(self):
        a,b,c,d = self.q0,self.q1,self.q2,self.q3
        norm = simplify(a*a+b*b+c*c+d*d)
        if norm != 1:
            raise ValueError('quaternion must have unit norm')
        return Matrix([[a*a+b*b-c*c-d*d,2*(b*c+a*d),2*(b*d-a*c)],
                       [2*(b*c-a*d),a*a-b*b+c*c-d*d,2*(c*d+a*b)],
                       [2*(b*d+a*c),2*(c*d-a*b),a*a-b*b-c*c+d*d]])
