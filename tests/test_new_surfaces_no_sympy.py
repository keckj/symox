"""All seven completed import surfaces must operate without a SymPy runtime."""
import subprocess
import sys


def test_new_surfaces_without_sympy():
    subprocess.run([sys.executable,'-c',r'''
import sys
class NoSymPy:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):
   raise RuntimeError('runtime SymPy dependency: '+name)
sys.meta_path.insert(0,NoSymPy())
import symox as s
from symox import *
from symox.combinatorics import (SymmetricGroup, GrayCode, Prufer, RGS_unrank,
 RGS_rank, Subset, free_group, RubikGroup, cube, Cycle)
from symox._transforms import MellinTransformStripError
R,u,v=ring('u,v',QQ)
assert ((u+v)*(u-v)).exquo(u-v)==u+v
F,a,b=field('a,b',ZZ)
assert (a*a-b*b)/(a-b)==a+b
x=Symbol('x');T=Poly(x*x-5,x)
Z,d=round_two(T)
assert d==5
for P in prime_decomp(11,ZK=Z):assert prime_valuation(P**3,P)==3
assert field_isomorphism(sqrt(2),sqrt(3)) is None
lo,hi=isolate(sqrt(2))
assert (sqrt(2)-lo).is_positive and (hi-sqrt(2)).is_positive
assert RGS_rank(RGS_unrank(12,5))==12
assert Prufer.unrank(5,4).rank==5
assert GrayCode(3,rank=5).current=='111'
assert Subset([1],[0,1,2]).rank_binary==2
G=SymmetricGroup(3);pc=G.polycyclic_group()
for p in G.elements:assert len(pc.collector.exponent_vector(p))==2
W,g,h=free_group('g,h')
assert g*h*h**-1*g**-1==W.identity
assert RubikGroup(3).degree==54 and cube.pgroup.order()==24
assert Cycle(0,1).list()==[1,0]
assert len(PowerSet(FiniteSet(1,2)))==4
assert DisjointUnion(FiniteSet(1),FiniteSet(2)).contains((2,1))==true
assert ComplexRegion(ProductSet(Interval(0,1),Interval(0,1))).contains(I)==true
assert imageset(x,x*x,FiniteSet(-1,1))==FiniteSet(1)
assert epath('/Symbol',x+1)==[x]
assert FU['TR2'](tan(x))==sin(x)/cos(x)
assert kroneckersimp(1+KroneckerDelta(x,0)*KroneckerDelta(x,1))==1
assert 'sympy' not in sys.modules
'''],check=True)
