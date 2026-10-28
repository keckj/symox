from pathlib import Path
import hashlib,json
import sympy as sp
from sympy.solvers.ode.tests import test_single
assert sp.__version__=='1.14.0'
rows=test_single._get_all_examples();assert len(rows)==315
cases=[]
for row in rows:
 cases.append({'id':row['example_name'],'equation':sp.srepr(row['eq'],order='none'),'function':sp.srepr(row['func'],order='none'),'hint':row['hint'],'simplify_flag':row['simplify_flag'],'reference_slow_metadata':{k:row[k] for k in ['slow','dsolve_too_slow','checkodesol_too_slow','checkodesol_XFAIL']}})
assert len({r['id'] for r in cases})==315
record={'reference_version':sp.__version__,'reference_dataset_file':test_single.__file__,'reference_dataset_sha256':hashlib.sha256(Path(test_single.__file__).read_bytes()).hexdigest(),'cases':cases,'methods':{'classify':'classify_ode(eq,func)','hint':'dsolve(eq,func,hint=case_hint,simplify=case_simplify_flag)','default':'dsolve(eq,func)'},'interpretation':'Presentation and error observations, not solution-completeness certification. All 315 cases included regardless of reference slow metadata.'}
Path('/home/keckj/.cache/symox-validation/corpus-replay-20261026/ode-cases.json').write_text(json.dumps(record,indent=2))
print('Exported315 reference cases')
