import pathlib,json,sympy as S
from sympy.integrals.rationaltools import ratint_logpart
x,t=S.symbols('x t');g=(x*x+1)*(x*x+x+1)*(x-2);f=S.cancel(2*(2*x)*g/(x*x+1)+3*(2*x+1)*g/(x*x+x+1)-g/(x-2))
cases={'quadratic':[(1,x*x+x+1,'ZZ')],'quartic':[(1,x**4+x+1,'ZZ')],'nonconstant_numerator':[(x**3+1,x**4+x+1,'ZZ')],'multiplicity_groups':[(f,g,'QQ')],'distinct_arguments':[(x+k,x**4+(k+1)*x+k*k+1,'ZZ') for k in range(1,9)]}
def terms(p):return [[list(m),str(c.p),str(c.q)] for m,c in p.terms()]
def payload(value):return [[terms(S.Poly(h.as_expr(),x,t,domain='QQ')),str(h.domain),terms(S.Poly(q.as_expr(),t,domain='QQ')),str(q.domain)] for h,q in value]
out={}
for name,items in cases.items():
 expected=[]
 for f,g,domain in items:
  a,b=S.Poly(f,x,domain=domain),S.Poly(g,x,domain=domain)
  assert S.gcd(a,b).degree()==0 and S.gcd(b,b.diff()).degree()==0
  expected.append(payload(ratint_logpart(a,b,x,t)))
 out[name]={'inputs':[[str(f),str(g),domain] for f,g,domain in items],'expected':expected}
pathlib.Path(__file__).with_name('fixtures.json').write_text(json.dumps(out,indent=2)+'\n')
