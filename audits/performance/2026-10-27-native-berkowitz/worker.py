from pathlib import Path
import argparse,hashlib,json,os,statistics,sys,time
p=argparse.ArgumentParser();p.add_argument('--sha',required=True);p.add_argument('--out',required=True);p.add_argument('--cpu',type=int,default=24);a=p.parse_args();os.sched_setaffinity(0,{a.cpu})
class NoRef:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('reference forbidden during native calculation')
sys.meta_path.insert(0,NoRef())
import symox as S
from symox.polys.matrices.ddm import DDM
import symox._core as core
assert hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest()==a.sha
fixture=json.loads((Path(__file__).parent/'fixtures.json').read_text());report={'native_sha256':a.sha,'symox_file':S.__file__,'affinity':sorted(os.sched_getaffinity(0)),'iterations':3,'rounds':5,'cases':{}}
def payload(results):
 return [[[str(v.as_expr().as_numer_denom()[0]),str(v.as_expr().as_numer_denom()[1])] for v in row] for row in results]
for name,case in fixture['cases'].items():
 K=getattr(S,case['domain']);n=case['n'];matrices=[DDM([[K(int(a),int(b)) if case['domain']=='QQ' else K(int(a)) for a,b in row] for row in v['matrix']],(n,n),K) for v in case['inputs']]
 def operation():return [m.charpoly() for m in matrices]
 expected=[v['coefficients'] for v in case['inputs']];initial=operation();assert payload(initial)==expected,name
 original_rows=[[id(row) for row in m] for m in matrices];original_values=[[[id(v) for v in row] for row in m] for m in matrices]
 selections=None
 if hasattr(core,'_ddm_berk'):
  old=core._ddm_berk;calls=[]
  def spy(*args):
   result=old(*args);calls.append(result is not None);return result
  core._ddm_berk=spy
  try:operation()
  finally:core._ddm_berk=old
  assert calls==[True]*len(matrices);selections=calls
 rounds=[]
 for _ in range(5):
  start=time.perf_counter_ns()
  for _ in range(3):last=operation()
  rounds.append((time.perf_counter_ns()-start)/3/len(matrices))
 assert payload(last)==expected
 assert original_rows==[[id(row) for row in m] for m in matrices]
 assert original_values==[[[id(v) for v in row] for row in m] for m in matrices]
 report['cases'][name]={'median_ns':statistics.median(rounds),'rounds_ns':rounds,'inputs':len(matrices),'digest':hashlib.sha256(json.dumps(expected).encode()).hexdigest(),'native_selections_outside_timers':selections}
assert 'sympy' not in sys.modules;report['reference_imported']=False
Path(a.out).write_text(json.dumps(report,indent=2)+'\n')
