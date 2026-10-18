"""Vector operators, moving-frame derivatives and kinematic equations."""
import symox as S
from ._core import Vector,Dyadic,ReferenceFrame,Point,_TIME,_COORDINATES,_vector,dynamicsymbols


def dot(first,second):return first.dot(second)
def cross(first,second):return first.cross(second)
def outer(first,second):return first.outer(second)


def express(expr,frame,frame2=None,variables=False):
    if isinstance(expr,Dyadic):return expr.express(frame,frame2)
    if isinstance(expr,Vector):return expr.express(frame,variables)
    expr = S.sympify(expr);mapping = {}
    for atom in expr.free_symbols:
        if atom in _COORDINATES:
            source,index = _COORDINATES[atom];mapping.update(source.variable_map(frame))
    return expr.subs(mapping,simultaneous=True)


def time_derivative(expr,frame,order=1):
    if order<0 or int(order)!=order:raise ValueError('derivative order must be a nonnegative integer')
    for i in range(int(order)):expr = expr.dt(frame) if isinstance(expr,(Vector,Dyadic)) else S.diff(expr,_TIME)
    return expr


def partial_velocity(vel_vecs,gen_speeds,frame):return [[velocity.diff(speed,frame) for speed in gen_speeds] for velocity in vel_vecs]


def gradient(scalar_field,frame):
    expression = express(scalar_field,frame)
    return sum((S.diff(expression,x)*axis for x,axis in zip(frame.varlist,frame._axes)),Vector(0))


def divergence(vect,frame):return S.simplify(sum(S.diff(v,x) for v,x in zip(vect.to_matrix(frame),frame.varlist)))


def curl(vect,frame):
    values = vect.to_matrix(frame);variables = frame.varlist
    return sum(((S.diff(values[(i+2)%3],variables[(i+1)%3])-S.diff(values[(i+1)%3],variables[(i+2)%3]))*frame._axes[i] for i in range(3)),Vector(0))


def is_conservative(field,frame):return curl(field,frame).simplify()==0

def is_solenoidal(field,frame):return divergence(field,frame)==0


def scalar_potential(field,frame):
    if not is_conservative(field,frame):raise ValueError('vector field is not conservative')
    result = S.Integer(0)
    for component,variable in zip(field.to_matrix(frame),frame.varlist):result += S.integrate(component-S.diff(result,variable),variable)
    return S.simplify(result)


def scalar_potential_difference(field,frame,point1,point2,origin):
    potential = scalar_potential(field,frame) if isinstance(field,Vector) else S.sympify(field)
    first,second = point1.pos_from(origin).to_matrix(frame),point2.pos_from(origin).to_matrix(frame)
    return S.simplify(potential.subs(dict(zip(frame.varlist,second)),simultaneous=True)-potential.subs(dict(zip(frame.varlist,first)),simultaneous=True))


def get_motion_params(frame,**kwargs):
    acceleration = kwargs.pop('acceleration',None);velocity = kwargs.pop('velocity',None);position = kwargs.pop('position',None)
    t0 = S.sympify(kwargs.pop('timevalue',0))
    if kwargs:raise TypeError('unknown motion parameters: '+', '.join(kwargs))
    if acceleration is not None:
        acceleration = _vector(acceleration);v0 = _vector(velocity or 0);r0 = _vector(position or 0)
        v = acceleration.express(frame).integrate((_TIME,t0,_TIME))+v0
        return acceleration,v,v.express(frame).integrate((_TIME,t0,_TIME))+r0
    if velocity is not None:
        velocity = _vector(velocity);r0 = _vector(position or 0)
        return velocity.dt(frame),velocity,velocity.express(frame).integrate((_TIME,t0,_TIME))+r0
    position = _vector(position or 0)
    return position.dt(frame).dt(frame),position.dt(frame),position


def kinematic_equations(speeds,coords,rot_type,rot_order=''):
    if len(speeds)!=3:raise ValueError('three angular speed components are required')
    if str(rot_type).lower() == 'quaternion':
        if len(coords)!=4:raise ValueError('four quaternion coordinates are required')
        a,b,c,d = map(S.sympify,coords);u,v,w = map(S.sympify,speeds)
        rates = (-(b*u+c*v+d*w)/2,(a*u+c*w-d*v)/2,(a*v+d*u-b*w)/2,(a*w+b*v-c*u)/2)
        return [S.diff(q,_TIME)-r for q,r in zip(coords,rates)]
    if len(coords)!=3:raise ValueError('three rotation coordinates are required')
    temporary = [S.Function('_rotation_coordinate'+str(i))(_TIME) for i in range(3)]
    N = ReferenceFrame('_kinematic_inertial');A = N.orientnew('_kinematic_body',rot_type,temporary,rot_order)
    angular = A.ang_vel_in(N).to_matrix(A);rates = [S.diff(q,_TIME) for q in temporary]
    jacobian = angular.jacobian(rates);rhs = jacobian.inv()*S.Matrix(speeds)
    replacements = dict(zip(temporary,map(S.sympify,coords)))
    return [S.trigsimp(S.diff(q,_TIME)-rhs[i].subs(replacements,simultaneous=True)) for i,q in enumerate(coords)]

# These printing helpers are also public through vector.functions in SymPy.
from .printing import init_vprinting, vlatex, vpprint, vprint, vsprint
