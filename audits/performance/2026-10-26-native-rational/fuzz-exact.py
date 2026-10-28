import hashlib,json,random,sys
from pathlib import Path
class NoReference:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):
   raise AssertionError('reference import in native calculation')
block=NoReference();sys.meta_path.insert(0,block)
import symox as S
x=S.Symbol('x');random.seed(20261026);results=[]
for case in range(96):
 a,b,c=[random.randrange(-5,6) for _ in range(3)]
 g=S.Rational(random.randrange(1,10),random.randrange(1,10))*(x+a)**(case%4+1)*(x*x+b*x+c*c+1)**(case%2+1)
 degree=S.degree(g,x);f=sum(S.Rational(random.randrange(-20,21),random.randrange(1,20))*x**i for i in range(int(degree)))
 p,q=S.Poly(f,x,domain='QQ'),S.Poly(g,x,domain='QQ')
 result=S._core._ratint_ratpart(p,q,x);assert result is not None
 results.append((str(f),str(g),[str(v) for v in result]))
sys.meta_path.remove(block)
import sympy as R
from sympy.integrals.rationaltools import ratint_ratpart
rx=R.Symbol('x')
for f,g,actual in results:
 f,g=R.Poly(f,rx),R.Poly(g,rx);a,b=map(R.sympify,actual)
 expected=ratint_ratpart(f,g,rx)
 assert all(R.cancel(v-w)==0 for v,w in zip((a,b),expected))
 assert R.cancel(R.diff(a,rx)+b-f.as_expr()/g.as_expr())==0
 d=R.Poly(R.fraction(R.cancel(b))[1],rx);assert R.gcd(d,d.diff()).degree()==0
result={'phase':'private debug candidate before release freeze','native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),'seed':20261026,'independent_exact_checks':len(results),'reference':'SymPy1.14.0','outputs_sha256':hashlib.sha256(json.dumps(results).encode()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
