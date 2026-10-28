"""Complete integration corpus replay with isolated no-reference native workers.

Only exact reconstructed outputs/held forms are classified. Different forms
need mathematical review; no numeric sample is treated as an equivalence proof.
"""
import argparse,concurrent.futures,hashlib,importlib.util,json,os,sys,time
from collections import Counter
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('ode_worker',HERE/'replay-ode.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
WORKER=r'''
import hashlib,importlib,json,resource,sys,time
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
backend=sys.argv[1]
if backend=='symox':
 class NoReference:
  def find_spec(self,name,*args):
   if name=='sympy' or name.startswith('sympy.'):raise ImportError('reference forbidden in native corpus worker')
 sys.meta_path.insert(0,NoReference())
m=importlib.import_module(backend);ns0={};exec('from '+backend+' import *',ns0)
print(json.dumps({'ready':True,'module':m.__file__,'version':m.__version__,'native_sha256':hashlib.sha256((Path(m.__file__).parent/'_core.abi3.so').read_bytes()).hexdigest() if backend=='symox' else None}),flush=True)
for line in sys.stdin:
 case=json.loads(line);ns=dict(ns0);start=time.monotonic()
 try:
  exec(case['base'],ns);exec(case['setup'],ns);value=eval(case['call'],ns)
  if type(value) is tuple:value=m.Tuple(*value)
  result={'ok':True,'text':str(value),'srepr':m.srepr(value,order='none')}
 except Exception as error:result={'ok':False,'error':type(error).__name__,'text':str(error)}
 result['elapsed_seconds']=time.monotonic()-start
 print(json.dumps(result),flush=True)
'''
helper.WORKER=WORKER

def main():
 p=argparse.ArgumentParser();p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--jobs',type=int,default=16);p.add_argument('--timeout',type=float,default=60);p.add_argument('--native-sha',required=True);args=p.parse_args()
 source=args.snapshot/'tools/integration_corpus.py';spec=importlib.util.spec_from_file_location('integration_cases',source);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 cases=[dict(c,base=module.BASE) for c in module.build_cases()];assert len({c['id'] for c in cases})==len(cases)
 args.out.mkdir(parents=True,exist_ok=False);(args.out/'cases.json').write_text(json.dumps(cases,indent=1));started=time.time()
 def shard(index):
  selected=cases[index::args.jobs];rows=[]
  with (args.out/f'{index:02d}-sympy.stderr').open('w') as er,(args.out/f'{index:02d}-symox.stderr').open('w') as en:
   workers={'sympy':helper.Worker('sympy',args.timeout,er),'symox':helper.Worker('symox',args.timeout,en)}
   try:
    for serial,case in enumerate(selected):
     pair={}
     for backend in (('sympy','symox') if serial%2==0 else ('symox','sympy')):pair[backend]=workers[backend].run(case)
     rows.append(dict(id=case['id'],**pair));(args.out/f'{index:02d}.json').write_text(json.dumps(rows,indent=1));print(index,serial+1,len(selected),flush=True)
    assert all(v['native_sha256']==args.native_sha for v in workers['symox'].identities)
    assert all(v['version']=='1.14.0' for v in workers['sympy'].identities)
    (args.out/f'{index:02d}-identities.json').write_text(json.dumps({b:w.identities for b,w in workers.items()},indent=2))
   finally:
    for w in workers.values():w.stop()
  return rows
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:results=sum(pool.map(shard,range(args.jobs)),[])
 assert Counter(r['id'] for r in results)==Counter(c['id'] for c in cases)
 # Reference interpretation is kept out of every native calculation worker.
 import sympy as sp
 namespace=module.srepr_namespace()
 for row in results:
  ref,nat=row['sympy'],row['symox']
  if any(v.get('error') in ('DEADLINE','MemoryError') for v in (ref,nat)):status='inconclusive_resource_limit'
  elif any(v.get('error') in ('WORKER_EXIT','PROTOCOL') for v in (ref,nat)):status='worker_failure_review_required'
  elif ref['ok'] and nat['ok']:
   try:
    a,b=(sp.sympify(eval(v['srepr'],namespace)) for v in (ref,nat))
    if a==b:status='exact_output_match'
    elif not module.evaluated(a) and not module.evaluated(b):status='both_held_observation'
    else:status='different_output_review_required'
   except Exception as error:status='comparison_error';row['comparison_error']=type(error).__name__+':'+str(error)
  elif not ref['ok'] and not nat['ok'] and (ref['error'],ref['text'])==(nat['error'],nat['text']):status='equal_error_observation'
  else:status='different_error_behavior_review_required'
  row['status']=status
 summary={'reference_version':sp.__version__,'native_sha256':args.native_sha,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'coverage':{'expected':len(cases),'finished':len(results),'exactly_once':True},'counts':dict(Counter(r['status'] for r in results)),'elapsed_seconds':time.time()-started,'limits':{'seconds_per_call':args.timeout,'memory_gib_per_worker':8},'interpretation':'Exact output and held/error observations, not completeness certification. No numerical sampling for equivalence or absence.'}
 (args.out/'classified.json').write_text(json.dumps(results,indent=1));(args.out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
