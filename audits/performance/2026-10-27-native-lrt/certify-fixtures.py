"""Independent exact rational trace certificates, outside all timed workers."""
import pathlib,json,sympy as S
root=pathlib.Path(__file__).resolve().parent;fixture=json.loads((root/'fixtures.json').read_text());x,t=S.symbols('x t');K=S.QQ.frac_field(x);checks=[]
def expression(terms,variables):return sum(S.Rational(n,d)*S.prod(v**power for v,power in zip(variables,monomial)) for monomial,n,d in terms)
for name,case in fixture.items():
 for (f,g,domain),pairs in zip(case['inputs'],case['expected']):
  total=0
  for hterms,hdomain,qterms,qdomain in pairs:
   h=expression(hterms,(x,t));q=S.Poly(expression(qterms,(t,)),t,domain=K)
   inverse=S.invert(S.Poly(h,t,domain=K),q)
   r=(S.Poly(t*S.diff(h,x),t,domain=K)*inverse).rem(q)
   total+=sum((r*S.Poly(t**j,t,domain=K)).rem(q).nth(j) for j in range(q.degree()))
  assert S.cancel(total-S.sympify(f)/S.sympify(g))==0
  checks.append({'case':name,'f':f,'g':g,'trace_residual':'0'})
(root/'exact-trace-certificates.json').write_text(json.dumps({'reference_version':S.__version__,'certificate':'Trace of t*h_x/h in QQ(x)[t]/q sums exactly to f/g; full pair coefficients/domains preserve branch representatives','checks':checks},indent=2)+'\n')
print(len(checks),'exact traces passed')
