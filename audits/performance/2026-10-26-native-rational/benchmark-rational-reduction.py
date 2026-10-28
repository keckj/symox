"""Release-only QQ rational reduction with exact independent certificates."""
import hashlib,json,statistics,sys,time
from pathlib import Path
class BlockReference:
 def find_spec(self,fullname,*args):
  if fullname=='sympy' or fullname.startswith('sympy.'):
   raise AssertionError('reference import during native measurement')
block=BlockReference();sys.meta_path.insert(0,block)
import symox as S
from symox.integrals.rationaltools import ratint_ratpart
x=S.Symbol('x')
result={'native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),'method':'9 distinct preconstructed QQ fractions per case; three batches of3; independent exact identity, square-free denominator and reference pair after all timing'}
checks=[]
for label,order in [('repeated_order2',2),('repeated_order6',6),('repeated_order12',12),('square_free_degree12',0)]:
 inputs=[]
 for variant in range(9):
  numerator=(variant+1)*x*x+S.Rational(variant+2,17)*x+1
  denominator=((x+S.Rational(variant+1,17))**order*(x*x+x+2) if order else x**12+x+variant+1)
  inputs.append((S.Poly(numerator,x,domain='QQ'),S.Poly(denominator,x,domain='QQ')))
 samples=[];outputs=[]
 for start in range(0,9,3):
  begin=time.perf_counter();batch=[ratint_ratpart(f,g,x) for f,g in inputs[start:start+3]]
  samples.append((time.perf_counter()-begin)/3);outputs.extend(batch)
 result[label]={'samples_seconds_per_call':samples,'median_seconds_per_call':statistics.median(samples),'values_digest':hashlib.sha256(str(outputs).encode()).hexdigest()}
 checks.extend((f,g,a,b) for (f,g),(a,b) in zip(inputs,outputs))
sys.meta_path.remove(block)
import sympy as sp
from sympy.integrals.rationaltools import ratint_ratpart as reference
rx=sp.Symbol('x')
for f,g,a,b in checks:
 f=sp.Poly(str(f.as_expr()),rx);g=sp.Poly(str(g.as_expr()),rx)
 a=sp.sympify(str(a));b=sp.sympify(str(b))
 expected=reference(f,g,rx)
 assert all(sp.cancel(v-w)==0 for v,w in zip((a,b),expected))
 assert sp.cancel(sp.diff(a,rx)+b-f.as_expr()/g.as_expr())==0
 denominator=sp.Poly(sp.fraction(sp.cancel(b))[1],rx)
 assert sp.gcd(denominator,denominator.diff()).degree()==0
result['independent_exact_checks']=len(checks)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
