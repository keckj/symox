"""Three-hinged parabolic arches, with exact support and section equilibrium."""
import symox as S
from .beam import _solve


class Arch:
    def __init__(self, left_support, right_support, **kwargs):
        self.left_support, self.right_support = tuple(map(S.sympify,left_support)),tuple(map(S.sympify,right_support))
        self.crown_x = S.sympify(kwargs.get('crown_x',(self.left_support[0]+self.right_support[0])/2));self.crown_y = kwargs.get('crown_y')
        if self.crown_y is None:
            # Vertex known in x: y = a (x-crown_x)^2 + crown_y.
            x1,y1=self.left_support;x2,y2=self.right_support
            a=(y2-y1)/((x2-self.crown_x)**2-(x1-self.crown_x)**2);self.crown_y=y1-a*(x1-self.crown_x)**2
        self.crown_y=S.sympify(self.crown_y);x=S.Symbol('x')
        coefficient=(self.left_support[1]-self.crown_y)/(self.left_support[0]-self.crown_x)**2
        self.get_shape_eqn=coefficient*(x-self.crown_x)**2+self.crown_y
        if S.simplify(self.get_shape_eqn.subs(x,self.right_support[0])-self.right_support[1])!=0:raise ValueError('supports and crown do not lie on a vertex-defined parabola')
        self.loads={'concentrated':{},'distributed':{}};self.reaction_force={};self.supports={'left':'hinge','right':'hinge'}
    def apply_load(self, order, label, start, mag, end=None, angle=None):
        if label in ('A','B') or any(label in d for d in self.loads.values()):raise ValueError('reserved or duplicate load label')
        start,mag=S.sympify(start),S.sympify(mag)
        if order == -1:
            if angle is None:raise ValueError('point load requires a direction in degrees')
            angle=S.sympify(angle)*S.pi/180
            self.loads['concentrated'][label]={'x':start,'y':self.get_shape_eqn.subs(S.Symbol('x'),start),'f_x':mag*S.cos(angle),'f_y':mag*S.sin(angle)}
        elif order == 0:
            if end is None or (S.sympify(end)-start).is_nonpositive is True:raise ValueError('distributed load requires end > start')
            self.loads['distributed'][label]={'start':start,'end':S.sympify(end),'f_y':mag}
        else:raise ValueError('supported load orders are -1 and 0')
    def remove_load(self,label):
        for values in self.loads.values():
            if label in values:del values[label];return
        raise ValueError('unknown load')
    def _resultants(self, cut):
        fx=fy=moment=S.Integer(0)
        for load in self.loads['concentrated'].values():
            active=S.Piecewise((1,cut>=load['x']),(0,True));fx+=active*load['f_x'];fy+=active*load['f_y'];moment+=active*(load['x']*load['f_y']-load['y']*load['f_x'])
        for load in self.loads['distributed'].values():
            stop=S.Min(cut,load['end']);span=S.Max(0,stop-load['start']);force=load['f_y']*span
            fy+=force;moment+=force*(load['start']+stop)/2
        return fx,fy,moment
    def solve(self):
        if any(kind!='hinge' for kind in self.supports.values()):raise NotImplementedError('only three-hinged arches without tie members are supported')
        ax,ay,bx,by=S.symbols('R_A_x R_A_y R_B_x R_B_y');xa,ya=self.left_support;xb,yb=self.right_support
        fx,fy,m=self._resultants(xb);lx,ly,lm=self._resultants(self.crown_x)
        equations=[ax+bx+fx,ay+by+fy,xa*ay-ya*ax+xb*by-yb*bx+m,
                   (xa-self.crown_x)*ay-(ya-self.crown_y)*ax+lm-self.crown_x*ly+self.crown_y*lx]
        self.reaction_force=_solve(equations,[ax,ay,bx,by])
        x=S.Symbol('x');y=self.get_shape_eqn;fx,fy,m=self._resultants(x)
        rx,ry=self.reaction_force[ax],self.reaction_force[ay]
        self._bending_moment=-(xa*ry-ya*rx+m-x*(ry+fy)+y*(rx+fx))
        slope=S.diff(y,x);norm=S.sqrt(1+slope**2)
        self._axial_force=-(rx+fx+(ry+fy)*slope)/norm
        self._shear_force=((rx+fx)*slope-(ry+fy))/norm
    def _at(self, name, pos):
        if not hasattr(self,name):raise ValueError('solve the arch first')
        return S.simplify(getattr(self,name).subs(S.Symbol('x'),pos))
    def bending_moment_at(self,pos):return self._at('_bending_moment',pos)
    def shear_force_at(self,pos):return self._at('_shear_force',pos)
    def axial_force_at(self,pos):return self._at('_axial_force',pos)
    def change_support_type(self, **kwargs):self.supports.update(kwargs)

__all__=['Arch']
