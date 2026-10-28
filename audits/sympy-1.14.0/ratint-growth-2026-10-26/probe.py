import hashlib,importlib,json,resource,sys,time
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
backend=sys.argv[1]
if backend=='symox':
 class NoReference:
  def find_spec(self,name,*args):
   if name=='sympy' or name.startswith('sympy.'):raise AssertionError('reference import in native computation')
 sys.meta_path.insert(0,NoReference())
S=importlib.import_module(backend);rt=importlib.import_module(backend+'.integrals.rationaltools');x=S.Symbol('x')
start=time.perf_counter();value=rt.ratint(1/(x**3+x+1),x);duration=time.perf_counter()-start
start=time.perf_counter();printed=str(value);printing=time.perf_counter()-start
result={'backend':backend,'version':S.__version__,'package':S.__file__,'call':'ratint(1/(x**3+x+1),x)','compute_seconds':duration,'print_seconds':printing,'characters':len(printed),'text_sha256':hashlib.sha256(printed.encode()).hexdigest(),'head':value.func.__name__,'interpretation':'Expression-size observation only; no equivalence or completeness claim from string length.'}
if backend=='symox':result['native_sha256']=hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest()
Path(sys.argv[2]).write_text(json.dumps(result,indent=2));Path(sys.argv[2]).with_suffix('.expression.txt').write_text(printed);print(json.dumps(result))
