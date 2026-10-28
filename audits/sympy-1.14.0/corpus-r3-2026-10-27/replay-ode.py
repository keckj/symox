"""315 reference ODE inputs, isolated backends, exact presentation observations.

Matches are observations of output/error behavior, not completeness proofs.
Every input is included; deadlines/memory limits are explicitly inconclusive.
"""
import argparse,concurrent.futures,hashlib,json,os,select,subprocess,sys,time
from collections import Counter
from pathlib import Path

WORKER=r'''
import hashlib,importlib,json,sys,time
from pathlib import Path
import resource
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
backend=sys.argv[1]
if backend=='symox':
 class NoReference:
  def find_spec(self,fullname,path=None,target=None):
   if fullname=='sympy' or fullname.startswith('sympy.'):
    raise ImportError('reference forbidden in native corpus worker')
 sys.meta_path.insert(0,NoReference())
m=importlib.import_module(backend)
ns=dict(vars(m));roots=importlib.import_module(backend+'.polys.rootoftools')
ns.setdefault('ComplexRootOf',roots.ComplexRootOf)
print(json.dumps({'ready':True,'module':m.__file__,'version':m.__version__,'native_sha256':hashlib.sha256((Path(m.__file__).parent/'_core.abi3.so').read_bytes()).hexdigest() if backend=='symox' else None}),flush=True)
for line in sys.stdin:
 item=json.loads(line);start=time.monotonic()
 try:
  eq=eval(item['equation'],ns);func=eval(item['function'],ns)
  method=item['method']
  if method=='classify':value=m.classify_ode(eq,func)
  elif method=='hint':value=m.dsolve(eq,func,hint=item['hint'],simplify=item['simplify_flag'])
  else:value=m.dsolve(eq,func)
  result={'ok':True,'text':str(value)}
 except Exception as error:
  result={'ok':False,'error':type(error).__name__,'text':str(error)}
 result['elapsed_seconds']=time.monotonic()-start
 print(json.dumps(result),flush=True)
'''

class Worker:
 def __init__(self,backend,timeout,log):
  self.backend,self.timeout,self.log=backend,timeout,log;self.process=None;self.identities=[]
 def start(self):
  self.process=subprocess.Popen([sys.executable,'-u','-c',WORKER,self.backend],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=self.log,text=True,env=os.environ.copy())
  if not select.select([self.process.stdout],[],[],45)[0]:
   self.stop();raise RuntimeError('worker startup timeout: '+self.backend)
  ready=json.loads(self.process.stdout.readline());assert ready['ready'];self.identities.append(ready)
 def stop(self):
  if self.process is not None:
   self.process.kill();self.process.wait();self.process=None
 def run(self,item):
  if self.process is None or self.process.poll() is not None:self.start()
  self.process.stdin.write(json.dumps(item)+'\n');self.process.stdin.flush()
  if not select.select([self.process.stdout],[],[],self.timeout)[0]:
   self.stop();return {'ok':False,'error':'DEADLINE','text':'inconclusive resource deadline'}
  line=self.process.stdout.readline()
  if not line:
   code=self.process.poll();self.stop();return {'ok':False,'error':'WORKER_EXIT','returncode':code,'text':'worker terminated; requires diagnosis'}
  try:return json.loads(line)
  except json.JSONDecodeError:
   self.stop();return {'ok':False,'error':'PROTOCOL','text':line[:1000]}

def verdict(a,b):
 if any(v.get('error') in ('WORKER_EXIT','PROTOCOL') for v in (a,b)):return 'worker_failure_review_required'
 if any(v.get('error') in ('DEADLINE','MemoryError') for v in (a,b)):return 'inconclusive_resource_limit'
 if a['ok'] and b['ok']:return 'equal_output' if a['text']==b['text'] else 'different_output_review_required'
 if not a['ok'] and not b['ok'] and (a['error'],a['text'])==(b['error'],b['text']):return 'equal_error_observation'
 return 'different_error_behavior_review_required'

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--cases',required=True,type=Path);parser.add_argument('--out',required=True,type=Path);parser.add_argument('--jobs',type=int,default=16);parser.add_argument('--timeout',type=float,default=60);parser.add_argument('--native-sha',required=True)
 args=parser.parse_args();source=json.loads(args.cases.read_text());cases=source['cases'];assert len(cases)==315
 args.out.mkdir(parents=True,exist_ok=False);started=time.time()
 def shard(index):
  selected=cases[index::args.jobs];rows=[]
  with (args.out/f'{index:02d}-sympy.stderr').open('w') as referr,(args.out/f'{index:02d}-symox.stderr').open('w') as naterr:
   workers={'sympy':Worker('sympy',args.timeout,referr),'symox':Worker('symox',args.timeout,naterr)}
   try:
    for serial,case in enumerate(selected):
     for method in ('classify','hint','default'):
      item=dict(case,method=method);results={}
      for backend in (('sympy','symox') if serial%2==0 else ('symox','sympy')):results[backend]=workers[backend].run(item)
      rows.append({'id':case['id'],'method':method,'verdict':verdict(results['sympy'],results['symox']),**results})
     (args.out/f'{index:02d}.json').write_text(json.dumps(rows,indent=1))
     print(index,serial+1,len(selected),flush=True)
    assert all(i['native_sha256']==args.native_sha for i in workers['symox'].identities)
    assert all(i['version']=='1.14.0' for i in workers['sympy'].identities)
    (args.out/f'{index:02d}-identities.json').write_text(json.dumps({b:w.identities for b,w in workers.items()},indent=2))
   finally:
    for worker in workers.values():worker.stop()
  return rows
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
  results=sum(executor.map(shard,range(args.jobs)),[])
 expected=Counter((c['id'],m) for c in cases for m in ('classify','hint','default'));actual=Counter((r['id'],r['method']) for r in results);assert actual==expected
 summary={'reference_version':'1.14.0','native_sha256':args.native_sha,'input_sha256':hashlib.sha256(args.cases.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'dataset_sha256':source['reference_dataset_sha256'],'methods':source['methods'],'coverage':{'expected':945,'finished':len(results),'exactly_once':actual==expected},'counts':{m:dict(Counter(r['verdict'] for r in results if r['method']==m)) for m in ('classify','hint','default')},'elapsed_seconds':time.time()-started,'limits':{'timeout_per_backend_call':args.timeout,'memory_gib_per_worker':8},'interpretation':'Exact output/error observations only. Different forms require mathematical review; equal forms do not establish complete solution families. Resource results are inconclusive.'}
 (args.out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
