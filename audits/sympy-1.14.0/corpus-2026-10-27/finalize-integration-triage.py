"""Exact coverage and scoped manual-proof classification, after frozen replay."""
from pathlib import Path
from collections import Counter
import hashlib,json,gzip,importlib.util
import sympy as S
root=Path(__file__).resolve().parent;out=root/'integration-r1-cb5c0c98'
rows=json.loads((out/'classified.json').read_text());cases=json.loads((out/'cases.json').read_text());summary=json.loads((out/'summary.json').read_text())
expected=Counter(c['id'] for c in cases);raw=[]
for index in range(16):raw+=json.loads((out/f'{index:02d}.json').read_text())
assert len(expected)==1210 and all(n==1 for n in expected.values())
assert Counter(r['id'] for r in rows)==expected==Counter(r['id'] for r in raw)
assert all(r[b]['ok'] for r in raw for b in ('sympy','symox'))
identities=[]
for index in range(16):
 ids=json.loads((out/f'{index:02d}-identities.json').read_text())
 assert ids['sympy'] and ids['symox']
 assert all(i['version']=='1.14.0' and '/site-packages/sympy/' in i['module'] for i in ids['sympy'])
 assert all(i['native_sha256']==summary['native_sha256'] and i['module']=='/home/keckj/symox-ci/wt/native-derivative-r1b-20261027/python/symox/__init__.py' for i in ids['symox'])
 identities.append(ids)
coverage={'cases':1210,'backend_calls':2420,'every_case_once_per_backend':True,'unique_case_ids':1210,'all_outputs_successful':True,'native_identities':sum(len(i['symox']) for i in identities),'reference_identities':sum(len(i['sympy']) for i in identities),'reference_version':'1.14.0','native_sha256':summary['native_sha256'],'limits':summary['limits'],'jobs':16,'native_sympy_import':'blocked by unchanged versioned worker meta_path finder'}
(root/'coverage-check.json').write_text(json.dumps(coverage,indent=2)+'\n')
controls=json.loads((root/'parameter-branch-controls.json').read_text());assert controls['passed_cases']==21 and controls['failed_cases']==0 and not controls['reference_imported'];checked={r['id']:r for r in controls['checks']};plan=json.loads((root/'parameter-gap-classification-plan.json').read_text());planned={r['id']:r for r in plan['cases']}
proof=[]
ns=dict(vars(S));from sympy.functions.elementary.piecewise import ExprCondPair
ns['ExprCondPair']=ExprCondPair
for row in rows:
 if row['status']!='different_output_review_required':continue
 ident=row['id'];item={'id':ident,'native_output':row['symox'],'reference_output':row['sympy']}
 if ident in checked:
  control=checked[ident];item.update(classification=control['classification'],proof_file='parameter-branch-proof.md',family=planned[ident]['family'],substitution_controls=len(control['controls']))
  if 'formula' in control:item.update(formula=control['formula'],condition=control['condition'])
  else:item['remaining_gap']='Free assumed-nonzero affine slope retains the complete Integral because its zero branch would be removed before contrary substitutions.'
 elif 'cos(x)/x' in ident:
  a,b=(eval(row[name]['srepr'],ns) for name in ('sympy','symox'))
  assert S.simplify(a+S.Ci(1))==0 and S.simplify(b+S.Ci(1))==0
  item.update(classification='exact_algebraic_equivalence',common_value='-Ci(1)',proof='Exact symbolic cancellation to -Ci(1); no numerical sample.')
 elif 'laplace_transform(1/(t+1)' in ident.replace(' ',''):
  facts=ident.split('[',1)[1].split(']',1)[0]
  if facts in ('none','cplx'):
   item.update(classification='native_complex_laplace_branch_gap_with_reference_prefactor_error',correct_value='exp(s)*expint(1,s)',domain='Re(s)>0',native_error='native minus correct = +I*pi*exp(s) for Im(s)>0 and -I*pi*exp(s) for Im(s)<0',proof_file='laplace-branch-review.md')
  else:item.update(classification='intentional_reference_prefactor_correction_on_real_domain',correct_value='-exp(s)*Ei(-s)',domain='s real, s>0',proof_file='laplace-branch-review.md')
 else:raise AssertionError('unclassified difference '+ident)
 proof.append(item)
assert len(proof)==36
report={'status':'final scoped exact triage','native_sha256':summary['native_sha256'],'counts':dict(Counter(r['classification'] for r in proof)),'former_parameter_gaps':{'total':21,'certified_conditional_results':20,'retained_affine_nonzero_gap':1},'cases':proof,'interpretation':'36 different raw outputs, with explicit domain/branch proofs or retained gaps.20newconditional capabilities are not global output equality. Four inherited complex Laplace gaps narrow an overbroad old classification; old570f archive unchanged. No numerical sampling establishes equivalence, convergence or completeness.'}
(root/'integration-triage-r1.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report['counts'],indent=2))
