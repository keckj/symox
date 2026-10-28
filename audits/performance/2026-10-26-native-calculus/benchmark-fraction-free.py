"""Release-only public domain-matrix fraction-free reduction, exact certificates."""
import hashlib,json,statistics,sys,time
from pathlib import Path
class BlockReference:
 def find_spec(self,fullname,*args):
  if fullname=='sympy' or fullname.startswith('sympy.'):
   raise AssertionError('reference import during native calculation')
block=BlockReference();sys.meta_path.insert(0,block)
import symox as S
from symox.polys.domains import ZZ,QQ
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix
from symox.polys.matrices.rref import _dm_rref_choose_method
result={'native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),'method':'9 distinct preconstructed full-row-rank rectangular matrices percase, 3 batches of3; public rref_den output allocation, method selection and conversions inside timers; exact independent checks after alltiming'}
checks=[]
for label,m,K,method in [('ZZ_FF_dense24',24,ZZ,'FF_dense'),('QQ_FF_dense16',16,QQ,'FF_dense'),('ZZ_AUTO32',32,ZZ,'auto'),('QQ_AUTO24',24,QQ,'auto'),('ZZ_AUTO8_GJ_control',8,ZZ,'auto')]:
 n=m+m//2;matrices=[];choices=[]
 for variant in range(9):
  rows=[[K(int(i==j) if j<m else ((3*i+7*j+variant)%17-8),1 if K==ZZ or j<m else j%3+1) for j in range(n)] for i in range(m)]
  for i in range(m-1,0,-1):rows[i]=[a+K((i+variant)%5+1)*b for a,b in zip(rows[i],rows[i-1])]
  rows.reverse();a=DomainMatrix.from_rep(DDM(rows,(m,n),K));matrices.append(a)
  choices.append(_dm_rref_choose_method(a,method,denominator=True))
 samples=[];outputs=[]
 for start in range(0,9,3):
  begin=time.perf_counter();batch=[a.rref_den(method=method) for a in matrices[start:start+3]]
  samples.append((time.perf_counter()-begin)/3);outputs.extend(batch)
 result[label]={'samples_seconds_per_call':samples,'median_seconds_per_call':statistics.median(samples),'selected_methods':choices,'values_digest':hashlib.sha256(str([(r.to_list(),str(d),p) for r,d,p in outputs]).encode()).hexdigest()}
 checks.extend((a,r,d,p,m,K) for a,(r,d,p) in zip(matrices,outputs))
sys.meta_path.remove(block)
import sympy as sp
def matrix(a):return sp.Matrix([[sp.Rational(str(v)) for v in row] for row in a.to_list()])
for a,r,d,p,rank,K in checks:
 expected,pivots=matrix(a).rref();den=sp.Rational(str(d))
 assert den!=0 and matrix(r)==den*expected and tuple(p)==pivots and len(p)==rank
 assert r.domain==K and d.parent()==K
result['independent_exact_checks']=len(checks)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
