"""Exact E1 formula/metadata certificates on recorded R2 outputs, no reference."""
import hashlib,json,pathlib,sys,resource
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
class NoReference:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('reference forbidden in native E1 controls')
sys.meta_path.insert(0,NoReference())
import symox as S
root=pathlib.Path(__file__).resolve().parent;out=root/'integration-r3-8a8a8e12'
rows=json.loads((out/'classified.json').read_text());cases={c['id']:c for c in json.loads((out/'cases.json').read_text())}
ns0={};exec('from symox import *',ns0)
checks=[]
for row in rows:
 if 'laplace_transform(1/(t+1)' not in row['id'].replace(' ',''):continue
 case=cases[row['id']];ns=dict(ns0);exec(case['base'],ns);exec(case['setup'],ns)
 value=eval(row['symox']['srepr'],ns);s=ns['s'];expected=S.exp(s)*S.expint(1,s)
 if row['id'].startswith('tr-nc['):formula=value;assert value==expected
 else:assert value==S.Tuple(expected,0,S.true);formula=value[0]
 assert not formula.has(S.Ei)
 assert S.cancel(formula.diff(s)-formula+1/s)==0
 controls=[]
 for point in (1+S.I,1-S.I):
  actual=formula.subs(s,point);target=S.exp(point)*S.expint(1,point)
  assert actual==target
  controls.append({'frequency':S.srepr(point),'value':S.srepr(actual)})
 facts=row['id'].split('[',1)[1].split(']',1)[0]
 checks.append({'id':row['id'],'formula':S.srepr(formula),'domain':'Re(s)>0','classification':'principal_E1_reference_correction','repairs_preexisting_complex_gap':facts in ('none','cplx'),'exact_ode':'F_prime-F+1/s=0','specialization_controls':controls,'analytic_proof':'laplace-e1-proof.md: integral substitution, half-plane analyticity and identity theorem; no Ei-cut shortcut'})
assert len(checks)==10 and sum(c['repairs_preexisting_complex_gap'] for c in checks)==4
assert 'sympy' not in sys.modules
report={'native_module':S.__file__,'native_sha256':hashlib.sha256((pathlib.Path(S.__file__).parent/'_core.abi3.so').read_bytes()).hexdigest(),'native_reference_imported':False,'checked_records':10,'repaired_complex_records':4,'checks':checks,'interpretation':'Exact reconstruction and analytic E1 certificate on recorded outputs. No original corpus call is re-executed; point substitutions are regression controls, not a proof by samples.'}
(root/'laplace-e1-record-controls.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ('native_sha256','native_reference_imported','checked_records','repaired_complex_records')},indent=2))
