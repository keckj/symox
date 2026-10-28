import argparse,hashlib,json,os,pathlib,statistics,sys,time
p=argparse.ArgumentParser();p.add_argument('--sha',required=True);p.add_argument('--out',required=True);p.add_argument('--variant',choices=['plain','vector'],required=True);p.add_argument('--cpu',type=int,default=24);a=p.parse_args();os.sched_setaffinity(0,{a.cpu})
class NoRef:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('reference forbidden during native benchmark')
sys.meta_path.insert(0,NoRef())
import symox as S
from symox.integrals import rationaltools as rt
if a.variant=='vector':import symox.vector
assert hashlib.sha256(pathlib.Path(S._core.__file__).read_bytes()).hexdigest()==a.sha
root=pathlib.Path(__file__).resolve().parent;x,t=S.symbols('x t')
def terms(poly):
 output=[]
 for monomial,c in sorted(poly.as_dict().items(),reverse=True):
  n,d=c.as_numer_denom();output.append([list(monomial),str(n),str(d)])
 return output
def lrt_payload(result):return [[terms(S.Poly(h.as_expr(),x,t,domain='QQ')),str(h.domain),terms(S.Poly(q.as_expr(),t,domain='QQ')),str(q.domain)] for h,q in result]
def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':'),sort_keys=True).encode()).hexdigest()
report={'sha256':a.sha,'module':S.__file__,'variant':a.variant,'cpu':a.cpu,'affinity':sorted(os.sched_getaffinity(0)),'reference_imported':False,'iterations':20,'rounds':5,'cases':{}}
for kind in ('lrt','atan'):
 fixture=json.loads((root/(kind+'-fixtures.json')).read_text())
 for name,case in fixture.items():
  pairs=[(S.Poly(S.sympify(f,locals={'x':x}),x,domain=d),S.Poly(S.sympify(g,locals={'x':x}),x,domain=d)) for f,g,d in case['inputs']]
  function=rt.ratint_logpart if kind=='lrt' else rt.log_to_atan
  def operation():return [function(f,g,x,t) if kind=='lrt' else function(f,g) for f,g in pairs]
  initial=operation()
  if kind=='lrt':certificate=[lrt_payload(v) for v in initial];assert certificate==case['expected']
  else:
   assert initial==[S.sympify(v,locals={'x':x}) for v in case['expected']]
   certificate=[S.srepr(v) for v in initial]
  native_name='_ratint_logpart' if kind=='lrt' else '_log_to_atan';old=getattr(S._core,native_name);calls=[]
  def spy(*args):
   value=old(*args);calls.append(value is not None);return value
  setattr(S._core,native_name,spy)
  try:assert operation()==initial
  finally:setattr(S._core,native_name,old)
  assert not calls or len(calls)==len(pairs) and all(calls)
  elapsed=[]
  for _ in range(5):
   start=time.perf_counter_ns()
   for _ in range(20):last=operation()
   elapsed.append((time.perf_counter_ns()-start)/20/len(pairs))
  final=[lrt_payload(v) for v in last] if kind=='lrt' else [S.srepr(v) for v in last]
  assert final==certificate
  report['cases'][kind+'/'+name]={'median_ns_per_pair':statistics.median(elapsed),'rounds_ns_per_pair':elapsed,'digest':digest(certificate),'pairs_per_batch':len(pairs),'native_selected_outside_timers':calls}
assert 'sympy' not in sys.modules
pathlib.Path(a.out).write_text(json.dumps(report,indent=2)+'\n')
