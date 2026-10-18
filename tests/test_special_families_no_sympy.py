"""The four completed special-function surfaces operate without SymPy."""
import subprocess
import sys


def test_special_families_without_runtime_sympy():
    subprocess.run([sys.executable,'-c',r'''
import sys
class NoSymPy:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):
   raise RuntimeError('runtime SymPy dependency: '+name)
sys.meta_path.insert(0,NoSymPy())
import symox as s
from symox.functions.special.hyper import HyperRep_power1, TupleArg
from symox.functions.special.mathieu_functions import MathieuBase
from symox.functions.special.bessel import SphericalHankelBase
x,a,q=s.symbols('x a q')
f=s.mathieuc(a,q,x)
assert isinstance(f,MathieuBase)
assert s.expand(f.diff(x,2)+(a-2*q*s.cos(2*x))*f)==0
assert s.mathieus(a,0,x)==s.sin(s.sqrt(a)*x)
assert s.multigamma(x,2)==s.sqrt(s.pi)*s.gamma(x)*s.gamma(x-s.Rational(1,2))
assert s.intlike(2.0) and not s.intlike(2.5)
h=s.hn1(1,x)
assert isinstance(h,SphericalHankelBase)
assert h.diff(x)==s.hn1(0,x)-2*h/x
assert len(s.jn_zeros(1,2))==2
assert s._core._numeric('marcumq',1,0,2).real>0
assert s._core._numeric('appellf1',1,1,1,2,.2,.3).real>1
G=s.meijerg([],[],[0],[],x)
assert G.argument==x and G.bm==TupleArg(0)
assert s.hyperexpand(G)==s.exp(-x)
assert HyperRep_power1(2,x).rewrite('nonrepsmall')==(1-x)**2
assert 'sympy' not in sys.modules
'''],check=True)
