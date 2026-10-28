import json,sys
class NoReference:
 def find_spec(self,name,*args):
  if name=="sympy" or name.startswith("sympy."):raise ImportError("reference forbidden in old ground control")
sys.meta_path.insert(0,NoReference())
import symox as S
from symox.polys.polymatrix import PolyMatrix
x,y,z=S.symbols('x y z')
cases=[]
for name,k in [('QQ',S.QQ),('RF',S.QQ.frac_field(y)),('multi',S.QQ.frac_field(y,z)),('ZRF',S.ZZ.frac_field(y))]:
 values=[S.Integer(0),S.Rational(1,2),S.Integer(2),S.Integer(3),S.Integer(4),S.Integer(5)]
 if name!='QQ':values[1]=(y+1)/(y+2);values[4]=4/(y+2)
 if name=='multi':values[1]+=1/(z+3)
 matrix=PolyMatrix(2,3,values,x,ring=k)
 baseline=(matrix.rref(),matrix.nullspace())
 oldstr=S._core.Expr.__dict__['__str__'];oldrepr=S._core.Expr.__dict__['__repr__'];events=[]
 def hostile(self):
  events.append(1);raise RuntimeError('presentation must not be read by ground reduction')
 try:
  S._core.Expr.__str__=hostile;S._core.Expr.__repr__=hostile
  current=(matrix.rref(),matrix.nullspace())
  assert current==baseline and not events
 finally:
  S._core.Expr.__str__=oldstr;S._core.Expr.__repr__=oldrepr
 cases.append({'domain':name,'same_result':True,'old_presentation_callbacks':len(events)})
assert "sympy" not in sys.modules
print(json.dumps(cases))
