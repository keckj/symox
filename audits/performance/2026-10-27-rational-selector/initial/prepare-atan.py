import json,pathlib,sympy as S
from sympy.integrals.rationaltools import log_to_atan
p=pathlib.Path(__file__).resolve().parent;x=S.Symbol('x')
cases={'linear_branch':[('x','x+1','ZZ')], 'quadratic':[('x**2+1','x-1','QQ')], 'quartic':[('x**4+x+1','x**3+2','QQ')], 'distinct_arguments':[(f'x**3+{i}*x+1',f'x**2+{i+1}','QQ') for i in range(1,5)]}
fixtures={};certificates=[]
for name,inputs in cases.items():
 expected=[]
 for a,b,domain in inputs:
  f,g=S.Poly(a,x,domain=domain),S.Poly(b,x,domain=domain)
  result=log_to_atan(f,g);arguments=[];trace=[];F,G=f,g
  while True:
   if F.degree()<G.degree():F,G=-G,F
   F,G=F.to_field(),G.to_field();q,r=F.div(G)
   if r.is_zero:arguments.append(q.as_expr());break
   s,t,h=G.gcdex(-F);u=(F*s+G*t).quo(h);arguments.append(u.as_expr())
   assert G*s-F*t==h
   trace.append({'f':str(F.as_expr()),'g':str(G.as_expr()),'s':str(s.as_expr()),'t':str(t.as_expr()),'h':str(h.as_expr()),'u':str(u.as_expr())})
   F,G=s,t
  assert sum(2*S.atan(u) for u in arguments)==result
  residual=S.cancel(S.diff(result,x)-2*(f.diff()*g-f*g.diff()).as_expr()/(f*f+g*g).as_expr())
  assert residual==0
  expected.append(str(result));certificates.append({'case':name,'f':a,'g':b,'arguments_in_recurrence_order':list(map(str,arguments)),'euclid_trace':trace,'rational_derivative_residual':'0','principal_branch_representative':str(result)})
 fixtures[name]={'inputs':inputs,'expected':expected}
(p/'atan-fixtures.json').write_text(json.dumps(fixtures,indent=2)+'\n')
(p/'atan-exact-certificates.json').write_text(json.dumps({'reference':S.__version__,'checks':certificates},indent=2)+'\n')
print(len(certificates),'Atan exact recurrence/derivative certificates prepared')
