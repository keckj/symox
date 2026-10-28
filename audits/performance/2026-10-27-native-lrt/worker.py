import argparse,hashlib,json,os,pathlib,statistics,sys,time
p=argparse.ArgumentParser();p.add_argument('--sha',required=True);p.add_argument('--out',required=True);p.add_argument('--cpu',type=int,default=24);a=p.parse_args()
os.sched_setaffinity(0,{a.cpu})
class NoRef:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('reference forbidden during native benchmark')
sys.meta_path.insert(0,NoRef())
import symox as S
from symox.integrals.rationaltools import ratint_logpart
import symox._core as core
assert hashlib.sha256(pathlib.Path(core.__file__).read_bytes()).hexdigest()==a.sha
root=pathlib.Path(__file__).resolve().parent;fixture=json.loads((root/'fixtures.json').read_text());x,t=S.symbols('x t')
def terms(poly):
 output=[]
 for monomial,c in sorted(poly.as_dict().items(),reverse=True):
  n,d=c.as_numer_denom();output.append([list(monomial),str(n),str(d)])
 return output

def payload(result):return [[terms(S.Poly(h.as_expr(),x,t,domain='QQ')),str(h.domain),terms(S.Poly(q.as_expr(),t,domain='QQ')),str(q.domain)] for h,q in result]
def digest(value):return hashlib.sha256(json.dumps(value,separators=(',',':'),sort_keys=True).encode()).hexdigest()
report={'sha256':a.sha,'module':S.__file__,'cpu':a.cpu,'affinity':sorted(os.sched_getaffinity(0)),'reference_imported':False,'iterations':50,'rounds':5,'cases':{}}
for name,case in fixture.items():
 pairs=[(S.Poly(S.sympify(f,locals={'x':x}),x,domain=domain),S.Poly(S.sympify(g,locals={'x':x}),x,domain=domain)) for f,g,domain in case['inputs']]
 def operation():return [ratint_logpart(f,g,x,t) for f,g in pairs]
 initial=operation();certificate=[payload(r) for r in initial];assert certificate==case['expected'],name
 selected=None
 if hasattr(core,'_ratint_logpart'):
  old=core._ratint_logpart;calls=[]
  def spy(*args):
   value=old(*args);calls.append(value is not None);return value
  core._ratint_logpart=spy
  try:operation()
  finally:core._ratint_logpart=old
  assert len(calls)==len(pairs) and all(calls),calls
  selected=calls
 elapsed=[]
 for _ in range(5):
  start=time.perf_counter_ns()
  for _ in range(50):last=operation()
  elapsed.append((time.perf_counter_ns()-start)/50/len(pairs))
 assert [payload(r) for r in last]==certificate
 report['cases'][name]={'median_ns_per_pair':statistics.median(elapsed),'rounds_ns_per_pair':elapsed,'digest':digest(certificate),'pairs_per_batch':len(pairs),'native_selected_outside_timers':selected}
assert 'sympy' not in sys.modules
pathlib.Path(a.out).write_text(json.dumps(report,indent=2)+'\n')
