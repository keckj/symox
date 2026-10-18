"""Exact polynomial moments over oriented polygon and polyhedron boundaries.

Edge and triangle integrals use affine polynomial substitution and exact
simplex monomial moments. The boundary divergence formula handles each
homogeneous component separately. Clockwise planar polygons have positive
measure, matching the reference API's orientation convention.
"""
from collections import defaultdict
from functools import cmp_to_key
from itertools import product
from math import factorial,gcd,lcm
from operator import index
import symox as S
x,y,z=S.symbols('x y z')


def norm(point):
    values=point.values() if isinstance(point,dict) else point
    return S.sqrt(sum((S.sympify(c)**2 for c in values),S.Integer(0)))
def cross_product(v1,v2,v3):
    a=[S.sympify(b)-S.sympify(c) for b,c in zip(v2,v1)];b=[S.sympify(c)-S.sympify(d) for c,d in zip(v3,v1)]
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
def _dot(a,b):return sum((u*v for u,v in zip(a,b)),S.Integer(0))
def _points(segment):return segment.points if hasattr(segment,'points') else segment

def _polynomial_components(expr,dims):
    expr=S.sympify(expr);poly=S.Poly(expr,*dims);result=defaultdict(lambda:S.Integer(0))
    for powers,coefficient in poly.terms():result[sum(powers)]+=coefficient*S.Mul(*(symbol**power for symbol,power in zip(dims,powers)))
    return dict(result)

def decompose(expr,separate=False):
    expr=S.expand(S.sympify(expr));variables=sorted(expr.free_symbols,key=str)
    if separate:return set(expr.args) if expr.is_Add else {expr}
    return _polynomial_components(expr,variables) if variables else {0:expr}

def strip(monom):
    monom=S.sympify(monom)
    if monom==0:return S.Integer(0),S.Integer(0)
    variables=sorted(monom.free_symbols,key=str)
    coefficient=S.Poly(monom,*variables).LC() if variables else monom
    return coefficient,S.cancel(monom/coefficient)


def _edge_integral(segment,expr,dims):
    a,b=_points(segment);t=S.Dummy('edge_t',real=True)
    expression=S.sympify(expr).subs({v:S.sympify(p)+(S.sympify(q)-S.sympify(p))*t for v,p,q in zip(dims,a,b)},simultaneous=True)
    polynomial=S.Poly(S.expand(expression),t)
    value=sum((coefficient/S.Integer(power[0]+1) for power,coefficient in polynomial.terms()),S.Integer(0))
    return S.simplify(norm([S.sympify(q)-S.sympify(p) for p,q in zip(a,b)])*value)


def _triangle_integral(points,expr,dims,normal):
    a,b,c=points;u,v=S.Dummy('simplex_u'),S.Dummy('simplex_v')
    expression=S.sympify(expr).subs({symbol:S.sympify(a[i])+u*(S.sympify(b[i])-S.sympify(a[i]))+v*(S.sympify(c[i])-S.sympify(a[i])) for i,symbol in enumerate(dims)},simultaneous=True)
    polynomial=S.Poly(S.expand(expression),u,v);value=S.Integer(0)
    for (i,j),coefficient in polynomial.terms():value+=coefficient*S.Rational(factorial(i)*factorial(j),factorial(i+j+2))
    cross=cross_product(a,b,c);jacobian=_dot(cross,normal)/norm(normal)
    return S.simplify(jacobian*value)


def _rational_gcd(values):
    if not all(getattr(v,'is_Rational',False) for v in values):return S.Integer(1)
    denominator=lcm(*(int(v.q) for v in values));numerator=0
    for value in values:numerator=gcd(numerator,int(value*denominator))
    return S.Rational(numerator,denominator) if numerator else S.Integer(1)


def hyperplane_parameters(poly,vertices=None):
    result=[]
    if isinstance(poly,S.Polygon):
        points=list(poly.vertices)
        for p,q in zip(points,points[1:]+points[:1]):
            normal=(p.y-q.y,q.x-p.x);bound=q.x*p.y-q.y*p.x;factor=_rational_gcd([*normal,bound])
            result.append((tuple(v/factor for v in normal),bound/factor))
    else:
        if vertices is None:raise ValueError('spatial faces require vertex coordinates')
        for facet in poly:
            if len(facet)<3:raise ValueError('a face requires at least three vertices')
            points=[vertices[i] for i in facet[:3]];normal=cross_product(*points);factor=_rational_gcd(normal)
            if norm(normal)==0:raise ValueError('face vertices must span a plane')
            normal=[v/factor for v in normal];result.append((normal,_dot(normal,points[0])))
    return result


def integration_reduction(facets,index,a,b,expr,dims,degree):return _edge_integral(facets[index],expr,tuple(dims))
def lineseg_integrate(polygon,index,line_seg,expr,degree):return _edge_integral(line_seg,expr,(x,y,z))
def left_integral2D(m,index,facets,x0,expr,gens):
    expr=S.sympify(expr)
    return sum((norm([S.sympify(a)-S.sympify(b) for a,b in zip(point,x0)])*expr.subs(dict(zip(gens,point)),simultaneous=True) for point in _points(facets[index])),S.Integer(0))

def distance_to_side(point,line_seg,A):
    p,q=_points(line_seg);direction=[S.sympify(b)-S.sympify(a) for a,b in zip(p,q)]
    unit=[v/norm(direction) for v in direction];reverse=[-S.sympify(v)/norm(A) for v in A]
    sideways=cross_product([0,0,0],reverse,unit)
    return S.simplify(_dot([S.sympify(a)-S.sympify(b) for a,b in zip(p,point)],sideways))

def polygon_integrate(facet,hp_param,index,facets,vertices,expr,degree):
    points=[vertices[i] for i in facet];normal=hp_param[0]
    return S.simplify(sum((_triangle_integral((points[0],points[i],points[i+1]),expr,(x,y,z),normal) for i in range(1,len(points)-1)),S.Integer(0)))

def left_integral3D(facets,index,expr,vertices,hp_param,degree):
    points=[vertices[i] for i in facets[index]];origin=points[0]
    return S.simplify(sum((distance_to_side(origin,(p,q),hp_param[0])*_edge_integral((p,q),expr,(x,y,z)) for p,q in zip(points,points[1:]+points[:1])),S.Integer(0)))

def integration_reduction_dynamic(facets,index,a,b,expr,degree,dims,x_index,y_index,max_index,x0,monomial_values,monom_index,vertices=None,hp_param=None):
    if vertices is None:return integration_reduction(facets,index,a,b,expr,dims,degree)
    return polygon_integrate(facets[index],hp_param or (a,b),index,facets,vertices,expr,degree)


def _monomials(degree,dimension):
    degree=index(degree)
    if degree<0:raise ValueError('moment degree must be nonnegative')
    if degree>64:raise ValueError('explicit moment inventory is limited to degree 64')
    for powers in product(range(degree+1),repeat=dimension):
        if sum(powers)<=degree:yield S.Mul(*(symbol**power for symbol,power in zip((x,y,z),powers)))


def main_integrate(expr,facets,hp_params,max_degree=None):
    if max_degree is not None:return {monomial:main_integrate(monomial,facets,hp_params) for monomial in _monomials(max_degree,2)}
    if isinstance(expr,(list,tuple)):return {S.sympify(e):main_integrate(e,facets,hp_params) for e in expr}
    result=S.Integer(0)
    for degree,homogeneous in _polynomial_components(expr,(x,y)).items():
        result+=sum((S.sympify(b)/norm(a)*_edge_integral(facet,homogeneous,(x,y)) for facet,(a,b) in zip(facets,hp_params)),S.Integer(0))/S.Integer(degree+2)
    return S.simplify(result)


def main_integrate3d(expr,facets,vertices,hp_params,max_degree=None):
    if max_degree is not None:return {monomial:main_integrate3d(monomial,facets,vertices,hp_params) for monomial in _monomials(max_degree,3)}
    if isinstance(expr,(list,tuple)):return {S.sympify(e):main_integrate3d(e,facets,vertices,hp_params) for e in expr}
    result=S.Integer(0)
    for degree,homogeneous in _polynomial_components(expr,(x,y,z)).items():
        result+=sum((S.sympify(hp[1])/norm(hp[0])*polygon_integrate(facet,hp,i,facets,vertices,homogeneous,degree) for i,(facet,hp) in enumerate(zip(facets,hp_params))),S.Integer(0))/S.Integer(degree+3)
    return S.simplify(result)


def polytope_integrate(poly,expr=None,*,clockwise=False,max_degree=None):
    if clockwise:
        if not isinstance(poly,S.Polygon):raise TypeError('clockwise sorting requires a planar polygon')
        poly=S.Polygon(*point_sort(poly.vertices))
    if isinstance(poly,S.Polygon):
        facets=poly.sides;parameters=hyperplane_parameters(poly);integrate=lambda value,degree=None:main_integrate(value,facets,parameters,degree)
    elif poly and len(poly[0])==2 and len(poly[0][0])==2:
        vertices=[intersection(poly[i-1],poly[i],'plane2D') for i in range(len(poly))]
        if any(vertex is None for vertex in vertices):raise ValueError('consecutive hyperplanes must intersect')
        facets=[S.Segment(vertices[i],vertices[(i+1)%len(vertices)]) for i in range(len(vertices))]
        integrate=lambda value,degree=None:main_integrate(value,facets,poly,degree)
    else:
        if not poly or len(poly)<2:raise ValueError('a spatial polytope requires vertices and oriented faces')
        vertices,facets=poly[0],poly[1:];parameters=hyperplane_parameters(facets,vertices)
        integrate=lambda value,degree=None:main_integrate3d(value,facets,vertices,parameters,degree)
    if max_degree is not None:
        if expr is None:return integrate(0,max_degree)
        if not isinstance(expr,(list,tuple)):raise TypeError('a moment request requires a list of polynomials')
        return {S.sympify(value):integrate(value) for value in expr if max(_polynomial_components(value,(x,y,z)),default=0)<=max_degree}
    if expr is None:raise TypeError('a polynomial expression is required')
    return integrate(expr)


def gradient_terms(binomial_power=0,no_of_gens=2):
    n=index(binomial_power)
    if n<0:raise ValueError('gradient order must be nonnegative')
    if no_of_gens==2:return [[x**i*y**j,i,j,0] for i in range(n+1) for j in range(n-i+1)]
    if no_of_gens!=3:raise ValueError('gradient tables support two or three dimensions')
    return [[[[x**i*y**j*z**(degree-i-j),i,j,degree-i-j,degree,i,degree-i-j,0] for j in range(degree-i,-1,-1)] for i in range(degree+1)] for degree in range(n+1)]


def intersection(geom_1,geom_2,intersection_type):
    if intersection_type=='plane2D':
        a,b=geom_1;c,d=geom_2;det=S.sympify(a[0])*c[1]-S.sympify(a[1])*c[0]
        return None if det==0 else (S.cancel((b*c[1]-d*a[1])/det),S.cancel((d*a[0]-b*c[0])/det))
    if intersection_type in ('segment2D','segment3D'):
        first=geom_1 if hasattr(geom_1,'intersection') else S.Segment(*geom_1)
        second=geom_2 if hasattr(geom_2,'intersection') else S.Segment(*geom_2)
        values=first.intersection(second)
        return tuple(values[0]) if len(values)==1 and isinstance(values[0],S.Point) else None
    raise ValueError('unsupported intersection type '+str(intersection_type))

def is_vertex(ent):
    return isinstance(ent,S.Point) or isinstance(ent,(tuple,list)) and len(ent) in (2,3) and all(not hasattr(v,'shape') for v in ent)


def point_sort(poly,normal=None,clockwise=True):
    from ..geometry.entity import sign,compare
    points=list(poly.vertices if hasattr(poly,'vertices') else poly)
    if len(points)<2:return points
    dimension=len(points[0]);center=[sum((S.sympify(p[i]) for p in points),S.Integer(0))/len(points) for i in range(dimension)]
    if dimension==3:
        if normal is None:raise ValueError('spatial point sorting requires a plane normal')
        projection=next((i for i,value in enumerate(normal) if value!=0),None)
        axes=[i for i in range(3) if i!=projection]
        orientation=sign(normal[projection])*(-1 if projection==1 else 1)
        if orientation<0:axes.reverse()
        orientation=1
    else:axes=[0,1];orientation=1
    def order(a,b):
        aa=[S.sympify(a[i])-center[i] for i in axes];bb=[S.sympify(b[i])-center[i] for i in axes]
        half_a=0 if sign(aa[0])>0 or aa[0]==0 and sign(aa[1])>=0 else 1
        half_b=0 if sign(bb[0])>0 or bb[0]==0 and sign(bb[1])>=0 else 1
        if half_a!=half_b:return (half_a-half_b)*(1 if clockwise else -1)
        determinant=aa[0]*bb[1]-aa[1]*bb[0]
        if determinant!=0:return sign(determinant)*(1 if clockwise else -1)*orientation
        return compare(_dot(bb,bb),_dot(aa,aa))
    return sorted(points,key=cmp_to_key(order))


def best_origin(a,b,lineseg,expr):
    p,q=_points(lineseg);polynomial=S.Poly(expr,x,y)
    totals=[sum(powers[i] for powers,coefficient in polynomial.terms()) for i in range(2)]
    for axis in sorted(range(2),key=lambda i:totals[i],reverse=True):
        delta=q[axis]-p[axis]
        if delta!=0:
            parameter=-p[axis]/delta
            if parameter.is_nonnegative is True and (1-parameter).is_nonnegative is True:return tuple(S.simplify(v+parameter*(w-v)) for v,w in zip(p,q))
    return tuple(p)


def plot_polytope(poly):
    from ..plotting.plot import BaseSeries,Plot
    if not isinstance(poly,S.Polygon):raise NotImplementedError('polytope plotting currently requires a planar polygon')
    t=S.Dummy('t');series=[]
    for p,q in [side.points for side in poly.sides]:series.append(BaseSeries((p.x+t*(q.x-p.x),p.y+t*(q.y-p.y)),((t,0,1),),'parametric'))
    return Plot(*series)

def plot_polynomial(expr):
    from ..plotting import plot3d
    return plot3d(expr,(x,-10,10),(y,-10,10),show=False)

__all__=['best_origin','cross_product','decompose','distance_to_side','gradient_terms','hyperplane_parameters',
 'integration_reduction','integration_reduction_dynamic','intersection','is_vertex','left_integral2D','left_integral3D',
 'lineseg_integrate','main_integrate','main_integrate3d','norm','plot_polynomial','plot_polytope','point_sort',
 'polygon_integrate','polytope_integrate','strip']
