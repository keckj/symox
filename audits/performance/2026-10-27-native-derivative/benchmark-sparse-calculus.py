"""Sparse EX calculus release comparison; all results have exact certificates."""
import hashlib,json,statistics,sys,time
from pathlib import Path
class BlockReference:
 def find_spec(self,fullname,*args):
  if fullname=='sympy' or fullname.startswith('sympy.'):
   raise AssertionError('reference import during native calculation')
block=BlockReference();sys.meta_path.insert(0,block)
import symox as S
x=S.Symbol('x');result={'native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),'method':'9 distinct preconstructed three-term EX polynomials percase;3 batches of3; exact independent polynomial/value checks after all timers; no rounded/grouped-order corrections used as speedups'};checks=[]
for label,degree,operation in [('diff_degree1000',1000,'diff'),('diff_degree4000',4000,'diff'),('integrate_degree4000',4000,'integrate'),('eval_zero_degree4000',4000,'zero')]:
 inputs=[S.Poly.from_dict({(degree,):S.Rational(v+1,7),(degree//2,):S.Rational(v+2,11),(0,):S.Rational(v+3,13)},x,domain='EX') for v in range(9)]
 def execute(p):
  if operation=='diff':return p.diff(x)
  if operation=='integrate':return p.integrate(x)
  return p.eval(0)
 samples=[];outputs=[]
 for start in range(0,9,3):
  begin=time.perf_counter();batch=[execute(a) for a in inputs[start:start+3]]
  samples.append((time.perf_counter()-begin)/3);outputs.extend(batch)
 result[label]={'samples_seconds_per_call':samples,'median_seconds_per_call':statistics.median(samples),'values_digest':hashlib.sha256(str([v.as_dict() if operation!='zero' else v for v in outputs]).encode()).hexdigest()}
 checks.extend((a,v,operation) for a,v in zip(inputs,outputs))
sys.meta_path.remove(block)
import sympy as R
rx=R.Symbol('x')
for p,actual,operation in checks:
 reference=R.Poly(str(p.as_expr()),rx,domain='EX')
 if operation=='zero':assert R.sympify(str(actual))==reference.eval(0)
 else:
  expected=reference.diff(rx) if operation=='diff' else reference.integrate(rx)
  assert R.Poly(str(actual.as_expr()),rx,domain='EX')==expected and actual.domain==S.EX
result['independent_exact_checks']=len(checks)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
