"""Predefined real plane and space with common coordinate transformations."""
from .. import symbols,sqrt,sin,cos,acos,atan2
from . import Manifold,Patch,CoordSystem

R2 = Manifold('R^2',2)
R2_origin = Patch('origin',R2)
x,y = symbols('x y',real=True)
r,theta = symbols('rho theta',nonnegative=True)
relations_2d = {
    ('rectangular','polar'): ((x,y),(sqrt(x*x+y*y),atan2(y,x))),
    ('polar','rectangular'): ((r,theta),(r*cos(theta),r*sin(theta))),
}
R2_r = CoordSystem('rectangular',R2_origin,(x,y),relations_2d)
R2_p = CoordSystem('polar',R2_origin,(r,theta),relations_2d)
for target in (R2,R2_origin,R2_r):
    target.x,target.y = R2_r.base_scalars()
    target.e_x,target.e_y = R2_r.base_vectors()
    target.dx,target.dy = R2_r.base_oneforms()
for target in (R2,R2_origin,R2_p):
    target.r,target.theta = R2_p.base_scalars()
    target.e_r,target.e_theta = R2_p.base_vectors()
    target.dr,target.dtheta = R2_p.base_oneforms()

R3 = Manifold('R^3',3)
R3_origin = Patch('origin',R3)
x,y,z = symbols('x y z',real=True)
rho,psi,r,theta,phi = symbols('rho psi r theta phi',nonnegative=True)
relations_3d = {
    ('rectangular','cylindrical'): ((x,y,z),(sqrt(x*x+y*y),atan2(y,x),z)),
    ('cylindrical','rectangular'): ((rho,psi,z),(rho*cos(psi),rho*sin(psi),z)),
    ('rectangular','spherical'): ((x,y,z),(sqrt(x*x+y*y+z*z),acos(z/sqrt(x*x+y*y+z*z)),atan2(y,x))),
    ('spherical','rectangular'): ((r,theta,phi),(r*sin(theta)*cos(phi),r*sin(theta)*sin(phi),r*cos(theta))),
    ('cylindrical','spherical'): ((rho,psi,z),(sqrt(rho*rho+z*z),acos(z/sqrt(rho*rho+z*z)),psi)),
    ('spherical','cylindrical'): ((r,theta,phi),(r*sin(theta),phi,r*cos(theta))),
}
R3_r = CoordSystem('rectangular',R3_origin,(x,y,z),relations_3d)
R3_c = CoordSystem('cylindrical',R3_origin,(rho,psi,z),relations_3d)
R3_s = CoordSystem('spherical',R3_origin,(r,theta,phi),relations_3d)
for system,names in ((R3_r,('x','y','z')),(R3_c,('rho','psi','z')),(R3_s,('r','theta','phi'))):
    for i,name in enumerate(names):
        setattr(system,name,system.base_scalar(i))
        setattr(system,'e_'+name,system.base_vector(i))
        setattr(system,'d'+name,system.base_oneform(i))

__all__ = ['R2','R2_origin','relations_2d','R2_r','R2_p','R3','R3_origin',
           'relations_3d','R3_r','R3_c','R3_s']
