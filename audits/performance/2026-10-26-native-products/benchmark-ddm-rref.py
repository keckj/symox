"""Exact release-only DDM rref observations; independent reference checks after timing."""
import hashlib,json,statistics,sys,time
from pathlib import Path
class BlockReference:
 def find_spec(self,fullname,*args):
  if fullname=='sympy' or fullname.startswith('sympy.'):
   raise AssertionError('reference import during native timing')
block=BlockReference();sys.meta_path.insert(0,block)
import symox as S
from symox.polys.domains import QQ
from symox.polys.matrices.ddm import DDM
result={'native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),'method':'12 distinct preconstructed matrices per case; 3 batches of 4; output allocation inside public DDM.rref; independent exact reference RREF outside all timers'}
checks=[]
for label,m,n,rank in [('QQ8',8,8,8),('QQ20',20,20,20),('QQ_rectangular12x24_rank10',12,24,10)]:
 matrices=[]
 for variant in range(12):
  basis=[[QQ(int(i==j) if j<rank else ((3*i+7*j+variant)%17-8),1 if j<rank else j%5+1) for j in range(n)] for i in range(rank)]
  rows=basis+[[sum((basis[k][j]*QQ((i+2*k+variant)%5-2) for k in range(rank)),QQ.zero) for j in range(n)] for i in range(m-rank)]
  for i in range(m-1,0,-1):
   scale=QQ((i+variant)%5+1,3)
   rows[i]=[a+scale*b for a,b in zip(rows[i],rows[i-1])]
  rows.reverse();matrices.append(DDM(rows,(m,n),QQ))
 samples=[];outputs=[]
 for start in range(0,12,4):
  begin=time.perf_counter();batch=[a.rref() for a in matrices[start:start+4]]
  samples.append((time.perf_counter()-begin)/4);outputs.extend(batch)
 result[label]={'samples_seconds_per_call':samples,'median_seconds_per_call':statistics.median(samples),'values_digest':hashlib.sha256(str([(list(map(list,r)),p) for r,p in outputs]).encode()).hexdigest()}
 checks.extend((a,r,p,rank) for a,(r,p) in zip(matrices,outputs))
sys.meta_path.remove(block)
import sympy as sp
def matrix(a):return sp.Matrix([[sp.Rational(str(v)) for v in row] for row in a])
for a,r,p,rank in checks:
 expected,pivots=matrix(a).rref()
 assert matrix(r)==expected and tuple(p)==pivots and len(p)==rank
 assert all(v.parent()==QQ for row in r for v in row)
result['independent_exact_rref_checks']=len(checks)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
