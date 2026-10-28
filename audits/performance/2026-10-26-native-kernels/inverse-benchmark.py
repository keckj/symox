import argparse, hashlib, json, os, pathlib, statistics, subprocess, sys, time
parser=argparse.ArgumentParser()
parser.add_argument('--root', required=True)
parser.add_argument('--output')
parser.add_argument('--case')
args=parser.parse_args()
if args.case:
 sys.path.insert(0,args.root)
 import symox as s
 n,p,kind=map(str,args.case.split(':'));n=int(n);p=int(p)
 def matrix(shift):
  def entry(i,j):
   q=s.Rational(((i+2)*(j+3)+shift)%17-8,32)
   real=s.Float(n+2+s.Rational(shift,100),precision=p) if i==j else s.Float(q,precision=p)
   return real+(s.I*s.Float(s.Rational((i-j+shift)%7-3,64),precision=p) if kind=='CC' else 0)
  return s.Matrix(n,n,lambda i,j:entry(i,j))
 matrices=[matrix(j) for j in range(8)]
 def run(m):
  t=time.perf_counter_ns();v=m.inv();elapsed=(time.perf_counter_ns()-t)/1e9
  assert v.shape==(n,n)
  return elapsed,v
 cold,first=run(matrices[0])
 warm=[run(matrices[0])[0] for _ in range(7)]
 changing=[run(m) for m in matrices[1:]]
 changed=[t for t,v in changing]
 # All independent certificates are outside every timed call.
 import sympy as sp
 def exact_matrix(M):
  def convert(v):
   v=s.to_sympy(v)
   return v.xreplace({f:sp.Rational(f) for f in v.atoms(sp.Float)})
  return sp.Matrix(M.rows,M.cols,lambda i,j:convert(M[i,j]))
 for M,V in zip(matrices,[first]+[v for t,v in changing]):
  product=(exact_matrix(M).to_DM()*exact_matrix(V).to_DM()).to_Matrix()-sp.eye(n)
  for value in product:
   value=sp.expand_complex(value)
   assert sp.re(value)**2+sp.im(value)**2 < sp.Rational(2)**(-2*(p-10))

 print(json.dumps({'cold':cold,'warm':statistics.median(warm),'changed':statistics.median(changed),'exact_binary_certificates':len(matrices)}))
else:
 results={}
 env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
 for n in (3,10,20):
  for p in (53,140):
   for kind in ('RR','CC'):
    case=f'{n}:{p}:{kind}'
    rounds=[]
    for _ in range(3):
     output=subprocess.check_output([sys.executable,'-B',__file__,'--root',args.root,'--case',case],env=env,text=True)
     rounds.append(json.loads(output))
    results[case]={'rounds':rounds,**{mode:statistics.median(r[mode] for r in rounds) for mode in ('cold','warm','changed')}}
    print(case,json.dumps({k:round(results[case][k],6) for k in ('cold','warm','changed')}),flush=True)
 package=pathlib.Path(args.root)/'symox'
 native=next(package.glob('_core*.so'))
 payload={'root':args.root,'native_sha256':hashlib.sha256(native.read_bytes()).hexdigest(),'python':sys.version,'rounds':3,'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'warm_repeats':7,'changed_repeats':7,'results':results}
 pathlib.Path(args.output).write_text(json.dumps(payload,indent=2)+'\n')
