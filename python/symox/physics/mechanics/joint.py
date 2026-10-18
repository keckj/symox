"""Kinematic joints connecting bodies through moving attachment frames."""
import symox as S
from ..vector import ReferenceFrame,Point,Vector,dynamicsymbols
from ..vector._core import _TIME
from ._core import Particle,RigidBody,_column


def _symbols(value,count,prefix,name):
    if count==0:return S.zeros(0,1)
    if value is None:value = [dynamicsymbols(prefix+str(i)+'_'+name) for i in range(count)]
    elif isinstance(value,S.Expr):value = [value]
    value = list(value)
    if len(value)!=count:raise ValueError('incorrect number of joint coordinates or speeds')
    return _column(value)


def _attachment(body,value,name):
    if value is None:return body.masscenter
    if isinstance(value,Point):return value
    if isinstance(value,Vector):return body.masscenter.locatenew(name,value)
    raise TypeError('joint attachment must be a Point or displacement Vector')


def _frame(body,value,name):
    base = body.frame if isinstance(body,RigidBody) else None
    if isinstance(value,ReferenceFrame):return value
    if base is None:raise TypeError('particle joint needs an explicit intermediate ReferenceFrame')
    if value is None:return base
    if isinstance(value,Vector):
        axis = value.normalize();rotation_axis = base.x.cross(axis)
        if rotation_axis==0:
            if S.simplify(base.x.dot(axis)-1)==0:return base
            rotation_axis = base.y
        child = ReferenceFrame(name);child.orient_axis(base,rotation_axis,S.acos(base.x.dot(axis)));return child
    raise TypeError('intermediate frame must be a ReferenceFrame or alignment vector')


class Joint:
    def _initialize(self,name,parent,child,count,coordinates,speeds,parent_point,child_point,parent_interframe,child_interframe):
        if not isinstance(parent,(Particle,RigidBody)) or not isinstance(child,(Particle,RigidBody)):raise TypeError('joints require parent and child bodies')
        if parent is child:raise ValueError('parent and child must differ')
        self.name,self.parent,self.child = str(name),parent,child
        self.parent_point = _attachment(parent,parent_point,self.name+'_parent_point');self.child_point = _attachment(child,child_point,self.name+'_child_point')
        self.parent_interframe = _frame(parent,parent_interframe,self.name+'_parent_frame');self.child_interframe = _frame(child,child_interframe,self.name+'_child_frame')
        self.coordinates = _symbols(coordinates,count,'q',self.name);self.speeds = _symbols(speeds,count,'u',self.name)
        self.kdes = self.speeds-self.coordinates.diff(_TIME)
    def __repr__(self):return '%s: %s  parent: %s  child: %s'%(type(self).__name__,self.name,self.parent,self.child)
    def _orient(self,rotation,axis=None,angular_speed=None):
        self.child_interframe.orient_dcm(self.parent_interframe,rotation)
        if axis is not None:self.child_interframe.set_ang_vel(self.parent_interframe,angular_speed*axis)
    def _place(self,displacement):self.child_point.set_pos(self.parent_point,displacement)


class PinJoint(Joint):
    def __init__(self,name,parent,child,coordinates=None,speeds=None,parent_point=None,child_point=None,parent_interframe=None,child_interframe=None,parent_axis=None,child_axis=None,joint_axis=None,parent_joint_pos=None,child_joint_pos=None):
        self._initialize(name,parent,child,1,coordinates,speeds,parent_point if parent_point is not None else parent_joint_pos,child_point if child_point is not None else child_joint_pos,parent_interframe if parent_interframe is not None else parent_axis,child_interframe if child_interframe is not None else child_axis)
        self.joint_axis = self.parent_interframe.x if joint_axis is None else joint_axis.normalize()
        self.child_interframe.orient_axis(self.parent_interframe,self.joint_axis,self.coordinates[0])
        self.child_interframe.set_ang_vel(self.parent_interframe,self.speeds[0]*self.joint_axis);self._place(Vector(0))


class PrismaticJoint(Joint):
    def __init__(self,name,parent,child,coordinates=None,speeds=None,parent_point=None,child_point=None,parent_interframe=None,child_interframe=None,parent_axis=None,child_axis=None,joint_axis=None,parent_joint_pos=None,child_joint_pos=None):
        self._initialize(name,parent,child,1,coordinates,speeds,parent_point if parent_point is not None else parent_joint_pos,child_point if child_point is not None else child_joint_pos,parent_interframe if parent_interframe is not None else parent_axis,child_interframe if child_interframe is not None else child_axis)
        self.joint_axis = self.parent_interframe.x if joint_axis is None else joint_axis.normalize()
        self._orient(S.eye(3));self._place(self.coordinates[0]*self.joint_axis)
        self.child_point.set_vel(self.parent_interframe,self.speeds[0]*self.joint_axis)


class CylindricalJoint(Joint):
    def __init__(self,name,parent,child,rotation_coordinate=None,translation_coordinate=None,rotation_speed=None,translation_speed=None,parent_point=None,child_point=None,parent_interframe=None,child_interframe=None,joint_axis=None):
        coordinates = [rotation_coordinate or dynamicsymbols('q0_'+str(name)),translation_coordinate or dynamicsymbols('q1_'+str(name))]
        speeds = [rotation_speed or dynamicsymbols('u0_'+str(name)),translation_speed or dynamicsymbols('u1_'+str(name))]
        self._initialize(name,parent,child,2,coordinates,speeds,parent_point,child_point,parent_interframe,child_interframe)
        self.rotation_coordinate,self.translation_coordinate = self.coordinates
        self.rotation_speed,self.translation_speed = self.speeds
        self.joint_axis = self.parent_interframe.x if joint_axis is None else joint_axis.normalize()
        self.child_interframe.orient_axis(self.parent_interframe,self.joint_axis,self.coordinates[0]);self.child_interframe.set_ang_vel(self.parent_interframe,self.speeds[0]*self.joint_axis)
        self._place(self.coordinates[1]*self.joint_axis);self.child_point.set_vel(self.parent_interframe,self.speeds[1]*self.joint_axis)


class PlanarJoint(Joint):
    def __init__(self,name,parent,child,rotation_coordinate=None,planar_coordinates=None,rotation_speed=None,planar_speeds=None,parent_point=None,child_point=None,parent_interframe=None,child_interframe=None):
        coordinates = [rotation_coordinate or dynamicsymbols('q0_'+str(name))]+list(planar_coordinates if planar_coordinates is not None else dynamicsymbols('q1_'+str(name)+' q2_'+str(name)))
        speeds = [rotation_speed or dynamicsymbols('u0_'+str(name))]+list(planar_speeds if planar_speeds is not None else dynamicsymbols('u1_'+str(name)+' u2_'+str(name)))
        self._initialize(name,parent,child,3,coordinates,speeds,parent_point,child_point,parent_interframe,child_interframe)
        self.rotation_coordinate,self.planar_coordinates = self.coordinates[0],_column(list(self.coordinates)[1:])
        self.rotation_speed,self.planar_speeds = self.speeds[0],_column(list(self.speeds)[1:])
        self.rotation_axis = self.parent_interframe.x;self.planar_vectors = (self.parent_interframe.y,self.parent_interframe.z)
        self.child_interframe.orient_axis(self.parent_interframe,self.rotation_axis,self.coordinates[0]);self.child_interframe.set_ang_vel(self.parent_interframe,self.speeds[0]*self.rotation_axis)
        self._place(self.coordinates[1]*self.planar_vectors[0]+self.coordinates[2]*self.planar_vectors[1])
        self.child_point.set_vel(self.parent_interframe,self.speeds[1]*self.planar_vectors[0]+self.speeds[2]*self.planar_vectors[1])


class SphericalJoint(Joint):
    def __init__(self,name,parent,child,coordinates=None,speeds=None,parent_point=None,child_point=None,parent_interframe=None,child_interframe=None,rot_type='BODY',amounts=None,rot_order=123):
        self._initialize(name,parent,child,3,coordinates,speeds,parent_point,child_point,parent_interframe,child_interframe)
        self.child_interframe.orient(self.parent_interframe,rot_type,list(self.coordinates) if amounts is None else amounts,rot_order)
        omega = self.child_interframe.ang_vel_in(self.parent_interframe).subs(dict(zip(self.coordinates.diff(_TIME),self.speeds)),simultaneous=True)
        self.child_interframe.set_ang_vel(self.parent_interframe,omega);self._place(Vector(0))


class WeldJoint(Joint):
    def __init__(self,name,parent,child,parent_point=None,child_point=None,parent_interframe=None,child_interframe=None):
        self._initialize(name,parent,child,0,[],[],parent_point,child_point,parent_interframe,child_interframe)
        self._orient(S.eye(3));self._place(Vector(0))


__all__ = ['PinJoint','PrismaticJoint','CylindricalJoint','PlanarJoint','SphericalJoint','WeldJoint']
