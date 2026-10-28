import json,gzip
from pathlib import Path
import sympy as S
root=Path(__file__).parent
rows=json.loads(gzip.decompress((root/'integration-570f/classified.json.gz').read_bytes()))
proof=[]
for row in rows:
 if row['status']!='different_output_review_required': continue
 if 'cos(x)/x' in row['id']:
  before=S.sympify(row['sympy']['srepr']);after=S.sympify(row['symox']['srepr'])
  assert S.simplify(before+S.Ci(1))==0
  assert S.simplify(after+S.Ci(1))==0
  proof.append({'id':row['id'],'classification':'exact_algebraic_equivalence','common_value':'-Ci(1)','proof':'Both exact expressions simplify to -Ci(1); no numerical sampling.'})
 elif 'laplace_transform(1/(t+1)' in row['id'].replace(' ',''):
  proof.append({'id':row['id'],'classification':'intentional_reference_correction','native_value':'-exp(s)*Ei(-s)','domain':'Re(s)>0','proof':'u=t+1 gives exp(s)*Integral(exp(-s*u)/u,(u,1,oo)); for real s>0, v=s*u gives exp(s)*E1(s)=-exp(s)*Ei(-s). The reference exp(-s) prefactor is incorrect. Analytic continuation on Re(s)>0 preserves this identity.'})
 else:
  proof.append({'id':row['id'],'classification':'retained_parameter_integral_capability_gap','proof':'No unconditional pole-crossing analytic continuation is accepted as the ordinary integral. Conditional finite-parameter branches remain to be implemented.'})
assert len(proof)==36
from collections import Counter
result={'counts':dict(Counter(x['classification'] for x in proof)),'cases':proof,'policy':'Exact symbolic proof or explicit gap; no numerical samples establish absence, convergence or completeness.'}
(root/'integration-triage.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result['counts']))
