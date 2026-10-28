"""Only exact rational identities/affine integration-constant changes; no branch rewrites."""
import gzip,json
from pathlib import Path
import sympy as S
root=Path(__file__).parent
rows=[row for p in sorted((root/'ode-570f').glob('[0-9][0-9].json.gz')) for row in json.loads(gzip.decompress(p.read_bytes())) if row['verdict']=='different_output_review_required']
proof=[]
x=S.Symbol('x');C=S.Symbol('C1')
for row in rows:
 record={'id':row['id'],'method':row['method'],'classification':'unresolved_output_difference'}
 try:
  r=S.sympify(row['sympy']['text']);n=S.sympify(row['symox']['text'])
  if isinstance(r,S.Equality) and isinstance(n,S.Equality) and r.lhs==n.lhs:
   symbols=sorted((r.rhs.free_symbols|n.rhs.free_symbols),key=str)
   if r.rhs.is_rational_function(*symbols) and n.rhs.is_rational_function(*symbols):
    if S.cancel(r.rhs-n.rhs)==0:
     record.update(classification='exact_rational_identity_on_common_domain',proof='Exact QQ(I) rational-function cancellation is zero; no transcendental or branch simplification.')
    elif C in symbols and S.diff(r.rhs,C,2)==0 and S.diff(n.rhs,C,2)==0:
     rn=S.cancel(S.diff(r.rhs,C)/S.diff(n.rhs,C));shift=S.cancel((r.rhs.subs(C,0)-n.rhs.subs(C,0))/S.diff(n.rhs,C))
     if not (rn.free_symbols|shift.free_symbols) and rn.is_zero is False and S.cancel(n.rhs.subs(C,rn*C+shift)-r.rhs)==0:
      record.update(classification='bijective_affine_constant_reparametrization_on_common_domain',native_C1=str(rn*C+shift),proof='Nonzero constant slope, constant intercept and exact rational identity after substitution. Local/common domain only; no global completeness inferred.')
 except Exception as error:record['review_note']=type(error).__name__+': '+str(error)[:200]
 proof.append(record)
from collections import Counter
result={'counts':dict(Counter(row['classification'] for row in proof)),'cases':proof,'scope':'53 differing-output observations only. Certificates concern identities on the common expression domain; singular points and complete families require separate review.'}
(root/'ode-rational-triage.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['counts']))
