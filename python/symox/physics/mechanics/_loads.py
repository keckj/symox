"""Paths, wrapping surfaces and mechanical actuators."""
from abc import ABC,abstractmethod
import symox as S
from ..vector import Point,Vector
from ..vector._core import _TIME,_vector
from ._core import Force,Torque


def _positive(value):
    value = S.sympify(value)
    if value.is_nonpositive is True:raise ValueError('parameter must be positive')
    return value


class PathwayBase(ABC):
    def __init__(self,*attachments):
        if len(attachments)<2 or any(not isinstance(p,Point) for p in attachments):raise TypeError('pathway requires at least two attachment points')
        self.attachments = tuple(attachments)
    @property
    @abstractmethod
    def length(self):raise NotImplementedError
    @property
    def extension_velocity(self):return S.diff(self.length,_TIME)
    @abstractmethod
    def to_loads(self,force):raise NotImplementedError


class LinearPathway(PathwayBase):
    def __init__(self,*attachments):
        if len(attachments)!=2:raise ValueError('linear pathway requires two attachment points')
        super().__init__(*attachments)
    @property
    def length(self):return self.attachments[1].pos_from(self.attachments[0]).magnitude()
    def to_loads(self,force):
        first,second = self.attachments;direction = second.pos_from(first).normalize()
        return [Force(first,-force*direction),Force(second,force*direction)]


class ObstacleSetPathway(PathwayBase):
    def __init__(self,*attachments):
        if len(attachments)<3:raise ValueError('obstacle pathway requires at least three points')
        super().__init__(*attachments)
    @property
    def length(self):return S.Add(*(b.pos_from(a).magnitude() for a,b in zip(self.attachments,self.attachments[1:])))
    def to_loads(self,force):
        loads = {p:Vector(0) for p in self.attachments}
        for first,second in zip(self.attachments,self.attachments[1:]):
            for point,value in LinearPathway(first,second).to_loads(force):loads[point] += value
        return [Force(point,value) for point,value in loads.items()]


class WrappingGeometryBase(ABC):
    def __init__(self,radius,point):
        self.radius = _positive(radius)
        if not isinstance(point,Point):raise TypeError('surface center must be a Point')
        self.point = point
    @abstractmethod
    def geodesic_length(self,point_1,point_2):raise NotImplementedError
    @abstractmethod
    def geodesic_end_vectors(self,point_1,point_2):raise NotImplementedError


class WrappingSphere(WrappingGeometryBase):
    def point_on_surface(self,point):return S.simplify(point.pos_from(self.point).dot(point.pos_from(self.point))-self.radius**2)==0
    def _radii(self,p1,p2):
        a,b = p1.pos_from(self.point),p2.pos_from(self.point)
        for v in (a,b):
            difference = S.simplify(v.dot(v)-self.radius**2)
            if difference.is_zero is False:raise ValueError('attachment point is not on the wrapping sphere')
        return a,b
    def geodesic_length(self,point_1,point_2):
        a,b = self._radii(point_1,point_2);return self.radius*S.acos(S.simplify(a.dot(b)/self.radius**2))
    def geodesic_end_vectors(self,point_1,point_2):
        a,b = self._radii(point_1,point_2);normal = a.cross(b)
        if normal==0:raise ValueError('coincident or antipodal endpoints do not define unique tangents')
        return a.cross(normal).normalize(),normal.cross(b).normalize()


class WrappingCylinder(WrappingGeometryBase):
    def __init__(self,radius,point,axis):super().__init__(radius,point);self.axis = _vector(axis).normalize()
    def _geometry(self,p1,p2):
        a,b = p1.pos_from(self.point),p2.pos_from(self.point)
        ar,br = a-a.dot(self.axis)*self.axis,b-b.dot(self.axis)*self.axis
        for v in (ar,br):
            if S.simplify(v.dot(v)-self.radius**2).is_zero is False:raise ValueError('attachment point is not on the wrapping cylinder')
        angle = S.atan2(self.axis.dot(ar.cross(br)),ar.dot(br));angle = S.Mod(angle,2*S.pi)
        height = (b-a).dot(self.axis)
        return ar,br,angle,height
    def point_on_surface(self,point):
        value = point.pos_from(self.point);radial = value-value.dot(self.axis)*self.axis
        return S.simplify(radial.dot(radial)-self.radius**2)==0
    def geodesic_length(self,point_1,point_2):
        a,b,angle,height = self._geometry(point_1,point_2);return S.sqrt(self.radius**2*angle**2+height**2)
    def geodesic_end_vectors(self,point_1,point_2):
        a,b,angle,height = self._geometry(point_1,point_2)
        return (-angle*self.axis.cross(a)-height*self.axis).normalize(),(angle*self.axis.cross(b)+height*self.axis).normalize()


class WrappingPathway(PathwayBase):
    def __init__(self,attachment1,attachment2,geometry):
        super().__init__(attachment1,attachment2)
        if not isinstance(geometry,WrappingGeometryBase):raise TypeError('a wrapping surface is required')
        self.geometry = geometry
    @property
    def length(self):return self.geometry.geodesic_length(*self.attachments)
    def to_loads(self,force):
        first,second = self.geometry.geodesic_end_vectors(*self.attachments)
        return [Force(self.attachments[0],force*first),Force(self.attachments[1],force*second),Force(self.geometry.point,-force*(first+second))]


class ActuatorBase(ABC):
    @abstractmethod
    def to_loads(self):raise NotImplementedError


class ForceActuator(ActuatorBase):
    def __init__(self,force,pathway):
        if not isinstance(pathway,PathwayBase):raise TypeError('force actuator requires a pathway')
        self.force,self.pathway = S.sympify(force),pathway
    def to_loads(self):return self.pathway.to_loads(self.force)


class LinearSpring(ForceActuator):
    def __init__(self,stiffness,pathway,equilibrium_length=0):
        self.stiffness,self.equilibrium_length = S.sympify(stiffness),S.sympify(equilibrium_length)
        super().__init__(-self.stiffness*(pathway.length-self.equilibrium_length),pathway)


class LinearDamper(ForceActuator):
    def __init__(self,damping,pathway):
        self.damping = S.sympify(damping);super().__init__(-self.damping*pathway.extension_velocity,pathway)


class DuffingSpring(ForceActuator):
    def __init__(self,linear_stiffness,nonlinear_stiffness,pathway,equilibrium_length=0):
        self.linear_stiffness,self.nonlinear_stiffness,self.equilibrium_length = map(S.sympify,(linear_stiffness,nonlinear_stiffness,equilibrium_length))
        extension = pathway.length-self.equilibrium_length
        super().__init__(-self.linear_stiffness*extension-self.nonlinear_stiffness*extension**3,pathway)


class CoulombKineticFriction(ForceActuator):
    def __init__(self,mu_k,f_n,pathway,*,v_s=None,sigma=None,mu_s=None):
        self.mu_k,self.f_n = S.sympify(mu_k),S.sympify(f_n)
        self.v_s,self.sigma,self.mu_s = map(S.sympify,(v_s if v_s is not None else S.Rational(1,100),sigma if sigma is not None else 0,mu_s if mu_s is not None else mu_k))
        velocity = pathway.extension_velocity
        magnitude = (self.mu_k+(self.mu_s-self.mu_k)*S.exp(-(velocity/self.v_s)**2))*self.f_n
        super().__init__(-magnitude*S.sign(velocity)-self.sigma*velocity,pathway)


class TorqueActuator(ActuatorBase):
    def __init__(self,torque,axis,target_frame,reaction_frame=None):
        self.torque,self.axis = S.sympify(torque),_vector(axis).normalize()
        self.target_frame,self.reaction_frame = target_frame,reaction_frame
    def to_loads(self):
        loads = [Torque(self.target_frame,self.torque*self.axis)]
        if self.reaction_frame is not None:loads.append(Torque(self.reaction_frame,-self.torque*self.axis))
        return loads
