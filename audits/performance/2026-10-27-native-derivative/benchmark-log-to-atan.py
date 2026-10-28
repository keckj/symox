"""Release-only exact atan construction, complete-expression certificates."""
import hashlib,importlib,json,statistics,sys,time
from pathlib import Path
class BlockReference:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):
   raise AssertionError('reference import during native calculation')
block=BlockReference();sys.meta_path.insert(0,block)
import symox as S
rt=importlib.import_module('symox.integrals.rationaltools');x=S.Symbol('x')
cases={'constant_den_degree16':([i%5-2 for i in range(16)]+[1],[3]),'euclid_degree4':([1,-2,0,3,1],[-2,1,0,1]),'euclid_degree8':([1,-3,2,4,-2,1,3,-1,1],[2,1,-2,0,3,1]),'shared_factor':([0,2,-1,3,1],[0,-3,1,1])}
result={'native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),'method':'9 distinct preconstructed QQ pairs percase;3 batches of3; exact full reference-expression and rational derivative checks after all timers; no atan addition identities'};checks=[]
for name,(a,b) in cases.items():
 inputs=[]
 for v in range(9):
  coefficients=a[:];coefficients[1 if name=='shared_factor' else 0]+=v
  f=S.Poly.from_list(list(reversed(coefficients)),x,domain=S.QQ);g=S.Poly.from_list(list(reversed(b)),x,domain=S.QQ);inputs.append((f,g))
 samples=[];outputs=[]
 for start in range(0,9,3):
  begin=time.perf_counter();batch=[rt.log_to_atan(f,g) for f,g in inputs[start:start+3]]
  samples.append((time.perf_counter()-begin)/3);outputs.extend(batch)
 result[name]={'samples_seconds_per_call':samples,'median_seconds_per_call':statistics.median(samples),'values_digest':hashlib.sha256(str(outputs).encode()).hexdigest()}
 checks.extend((f,g,value) for (f,g),value in zip(inputs,outputs))
sys.meta_path.remove(block)
import sympy as R
reference=importlib.import_module('sympy.integrals.rationaltools');rx=R.Symbol('x')
for f,g,value in checks:
 f=R.Poly(str(f.as_expr()),rx,domain=R.QQ);g=R.Poly(str(g.as_expr()),rx,domain=R.QQ)
 actual=R.sympify(str(value));expected=reference.log_to_atan(f,g)
 assert actual==expected
 target=2*(g.as_expr()*f.diff().as_expr()-f.as_expr()*g.diff().as_expr())/(f.as_expr()**2+g.as_expr()**2)
 assert R.cancel(R.diff(actual,rx)-target)==0
result['independent_exact_checks']=len(checks)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
