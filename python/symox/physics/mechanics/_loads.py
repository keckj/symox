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
    def __repr__(self):return f'{self.__class__.__name__}(radius={self.radius}, point={self.point})'
    def point_on_surface(self,point):
        point_vector = point.pos_from(self.point)
        squared = point_vector.dot(point_vector) if isinstance(point_vector,Vector) else point_vector**2
        return S.Eq(squared,self.radius**2) == True
    def geodesic_length(self,point_1,point_2):
        for point in (point_1,point_2):
            if not self.point_on_surface(point):
                raise ValueError(f'Geodesic length cannot be calculated as point {point} with radius '
                                 f'{point.pos_from(self.point).magnitude()} from the sphere\'s center {self.point} '
                                 f'does not lie on the surface of {self} with radius {self.radius}.')
        point_1_vector = point_1.pos_from(self.point).normalize()
        point_2_vector = point_2.pos_from(self.point).normalize()
        return self.radius*S.acos(point_2_vector.dot(point_1_vector))
    def geodesic_end_vectors(self,point_1,point_2):
        pA,pB,pO = point_1,point_2,self.point
        pA_vec,pB_vec = pA.pos_from(pO),pB.pos_from(pO)
        if pA_vec.cross(pB_vec) == 0:
            raise ValueError(f'Can\'t compute geodesic end vectors for the pair of points {pA} and {pB} on a sphere '
                             f'{self} as they are diametrically opposed, thus the geodesic is not defined.')
        return (pA_vec.cross(pB.pos_from(pA)).cross(pA_vec).normalize(),
                pB_vec.cross(pA.pos_from(pB)).cross(pB_vec).normalize())


def _directional_atan(numerator,denominator):
    """Central angle in [0, 2*pi) along the wrapping direction (reference algorithm)."""
    numerator,denominator = S.sympify(numerator),S.sympify(denominator)
    if numerator.is_number and denominator.is_number:
        angle = S.atan2(numerator,denominator)
        if angle < 0:angle += 2*S.pi
    elif numerator.is_number:
        raise NotImplementedError(f'Cannot compute a directional atan when the numerator {numerator} is numeric and the denominator {denominator} is symbolic.')
    elif denominator.is_number:
        raise NotImplementedError(f'Cannot compute a directional atan when the numerator {numerator} is symbolic and the denominator {denominator} is numeric.')
    else:
        ratio = S.sympify(S.trigsimp(numerator/denominator))
        if isinstance(ratio,S.tan):angle = ratio.args[0]
        elif ratio.is_Mul and ratio.args[0] == S.Integer(-1) and isinstance(ratio.args[1],S.tan):angle = 2*S.pi-ratio.args[1].args[0]
        else:raise NotImplementedError(f'Cannot compute a directional atan for the value {ratio}.')
    return angle


class WrappingCylinder(WrappingGeometryBase):
    def __init__(self,radius,point,axis):super().__init__(radius,point);self.axis = _vector(axis).normalize()
    def __repr__(self):return f'{self.__class__.__name__}(radius={self.radius}, point={self.point}, axis={self.axis})'
    def point_on_surface(self,point):
        relative_position = point.pos_from(self.point)
        point_vector = relative_position-relative_position.dot(self.axis)*self.axis
        squared = point_vector.dot(point_vector) if isinstance(point_vector,Vector) else point_vector**2
        return S.Eq(S.trigsimp(squared),self.radius**2) == True
    def geodesic_length(self,point_1,point_2):
        for point in (point_1,point_2):
            if not self.point_on_surface(point):
                raise ValueError(f'Geodesic length cannot be calculated as point {point} with radius '
                                 f'{point.pos_from(self.point).magnitude()} from the cylinder\'s center {self.point} '
                                 f'does not lie on the surface of {self} with radius {self.radius} and axis {self.axis}.')
        parallel_length = point_2.pos_from(point_1).dot(self.axis)
        p1 = point_1.pos_from(self.point);p1 = (p1-p1.dot(self.axis)*self.axis).normalize()
        p2 = point_2.pos_from(self.point);p2 = (p2-p2.dot(self.axis)*self.axis).normalize()
        central_angle = _directional_atan(S.cancel(p1.cross(p2).dot(self.axis)),S.cancel(p1.dot(p2)))
        return S.sqrt(parallel_length**2+(self.radius*central_angle)**2)
    def geodesic_end_vectors(self,point_1,point_2):
        p1o,p2o = point_1.pos_from(self.point),point_2.pos_from(self.point)
        if p1o == p2o:raise ValueError(f'Cannot compute geodesic end vectors for coincident points {point_1} and {point_2} as no geodesic exists.')
        n1 = p1o-p1o.dot(self.axis)*self.axis;n2 = p2o-p2o.dot(self.axis)*self.axis
        if n1 == n2:perp1 = perp2 = Vector(0)
        else:perp1,perp2 = self.axis.cross(n1).normalize(),-self.axis.cross(n2).normalize()
        geodesic_length = self.geodesic_length(point_1,point_2)
        parallel_length = point_2.pos_from(point_1).dot(self.axis)
        planar_arc_length = S.sqrt(geodesic_length**2-parallel_length**2)
        return ((planar_arc_length*perp1+parallel_length*self.axis).normalize(),
                (planar_arc_length*perp2-parallel_length*self.axis).normalize())


class WrappingPathway(PathwayBase):
    def __init__(self,attachment1,attachment2,geometry):
        super().__init__(attachment1,attachment2)
        if not isinstance(geometry,WrappingGeometryBase):raise TypeError('a wrapping surface is required')
        self.geometry = geometry
    @property
    def length(self):return self.geometry.geodesic_length(*self.attachments)
    def to_loads(self,force):
        pA,pB = self.attachments;pO = self.geometry.point
        pA_force,pB_force = self.geometry.geodesic_end_vectors(pA,pB)
        return [Force(pA,force*pA_force),Force(pB,force*pB_force),Force(pO,force*(-(pA_force+pB_force)))]


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
