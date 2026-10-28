import argparse,cProfile,hashlib,json,pathlib,pstats,sys
p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--sha',required=True);a=p.parse_args()
class NoRef:
 def find_spec(self,name,*args):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('reference forbidden')
sys.meta_path.insert(0,NoRef())
import symox as S
import symox.vector
from symox.integrals import rationaltools as rt
assert hashlib.sha256(pathlib.Path(S._core.__file__).read_bytes()).hexdigest()==a.sha
x=S.Symbol('x');f=S.Poly(x,x,domain='ZZ');g=S.Poly(x+1,x,domain='ZZ')
expected=2*S.atan(2*x+1)+S.pi/2
assert rt.log_to_atan(f,g)==expected
profile=cProfile.Profile();profile.enable()
for i in range(500):value=rt.log_to_atan(f,g)
profile.disable();assert value==expected and 'sympy' not in sys.modules
stats=pstats.Stats(profile)
rows=[dict(file=key[0],line=key[1],function=key[2],primitive_calls=v[0],calls=v[1],own_seconds=v[2],cumulative_seconds=v[3]) for key,v in stats.stats.items()]
rows.sort(key=lambda r:-r['cumulative_seconds'])
out=pathlib.Path(a.out);out.write_text(json.dumps({'native_sha256':a.sha,'package':S.__file__,'iterations':500,'reference_imported':False,'interpretation':'cProfile attribution for one vector linear atan fixture. Instrumented values are not benchmark times or speedup estimates.','functions':rows},indent=2)+'\n')
