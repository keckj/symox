"""Snell refraction, Fresnel amplitudes and paraxial imaging relations."""
import symox as S
from .medium import Medium


def _index(medium):return medium.refractive_index if isinstance(medium,Medium) else S.sympify(medium)


def refraction_angle(incident,medium1,medium2,normal=None,plane=None):
    ratio=_index(medium1)/_index(medium2)
    if isinstance(incident,(int,float,S.Expr)):
        angle=S.sympify(incident)
        if angle.is_number and (angle<0 or angle>S.pi/2):raise ValueError('incidence angle must lie in [0, pi/2]')
        argument=ratio*S.sin(angle)
        if (argument-1).is_positive is True:return S.Integer(0)
        return S.asin(argument)
    from ...geometry import Ray3D
    if isinstance(incident,Ray3D):direction=S.Matrix(incident.direction_ratio)
    elif isinstance(incident,(list,tuple,S.Matrix)):direction=S.Matrix(incident)
    else:raise TypeError('incident ray must be an angle, vector or Ray3D')
    if plane is not None:normal=S.Matrix(plane.normal_vector)
    if normal is None:raise ValueError('vector refraction requires a normal or plane')
    if hasattr(normal,'direction_ratio'):normal=S.Matrix(normal.direction_ratio)
    normal=S.Matrix(normal);direction=direction/S.sqrt(direction.dot(direction));normal=normal/S.sqrt(normal.dot(normal))
    cosine=-normal.dot(direction)
    if cosine.is_negative is True:normal=-normal;cosine=-cosine
    discriminant=1-ratio**2*(1-cosine**2)
    if discriminant.is_negative is True:return S.Integer(0)
    result=ratio*direction+(ratio*cosine-S.sqrt(discriminant))*normal
    return result.applyfunc(S.simplify)


def deviation(incident,medium1,medium2,normal=None,plane=None):
    result=refraction_angle(incident,medium1,medium2,normal,plane)
    if isinstance(incident,(int,float,S.Expr)):return S.sympify(incident)-result
    direction=S.Matrix(incident.direction_ratio if hasattr(incident,'direction_ratio') else incident)
    return S.acos(direction.dot(result)/(S.sqrt(direction.dot(direction))*S.sqrt(result.dot(result))))


def fresnel_coefficients(angle_of_incidence,medium1,medium2):
    angle=S.sympify(angle_of_incidence);n1,n2=_index(medium1),_index(medium2);transmitted=S.asin(n1*S.sin(angle)/n2)
    ci,ct=S.cos(angle),S.cos(transmitted)
    rp=(n2*ci-n1*ct)/(n2*ci+n1*ct);rs=(n1*ci-n2*ct)/(n1*ci+n2*ct)
    tp=2*n1*ci/(n2*ci+n1*ct);ts=2*n1*ci/(n1*ci+n2*ct)
    return [rp,rs,tp,ts]


def brewster_angle(medium1,medium2):return S.atan(_index(medium2)/_index(medium1))
def critical_angle(medium1,medium2):
    n1,n2=_index(medium1),_index(medium2)
    if (n2-n1).is_positive is True:raise ValueError('total internal reflection requires a higher incident refractive index')
    return S.asin(n2/n1)


def lens_makers_formula(n_lens,n_surr,r1,r2,d=0):
    n,ns,r1,r2,d=map(S.sympify,(n_lens,n_surr,r1,r2,d));ratio=n/ns
    return 1/((ratio-1)*(1/r1-1/r2+(ratio-1)*d/(ratio*r1*r2)))


def _imaging(focal_length,u,v,mirror):
    values=[None if value is None else S.sympify(value) for value in (focal_length,u,v)]
    if values.count(None)!=1:raise ValueError('provide exactly two of focal_length, u and v')
    f,u,v=values
    if f is None:return 1/(1/v+(1/u if mirror else -1/u))
    if u is None:return f*v/(v-f) if mirror else f*v/(f-v)
    return f*u/(u-f) if mirror else f*u/(u+f)


def mirror_formula(focal_length=None,u=None,v=None):return _imaging(focal_length,u,v,True)
def lens_formula(focal_length=None,u=None,v=None):return _imaging(focal_length,u,v,False)
def hyperfocal_distance(f,N,c):return S.sympify(f)**2/(S.sympify(N)*S.sympify(c))
def transverse_magnification(si,so):return -S.sympify(si)/S.sympify(so)


__all__=['refraction_angle','deviation','fresnel_coefficients','brewster_angle','critical_angle','lens_makers_formula','mirror_formula','lens_formula','hyperfocal_distance','transverse_magnification']
