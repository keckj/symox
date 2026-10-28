"""Fresh-process metadata audits; construct/import/AST preparation precede hook."""
import ast,importlib,json,sys,types
from pathlib import Path
import symox as s
nf=importlib.import_module('symox.polys.matrices.normalforms')
dense=importlib.import_module('symox.polys.matrices.dense')
from symox.polys.matrices.domainmatrix import DomainMatrix
route,arm=sys.argv[1:]
name={'hnf':'_hermite_normal_form','lu':'ddm_ilu','solve':'ddm_ilu_solve'}[route];module=nf if route=='hnf' else dense
node=next(n for n in ast.parse(Path(module.__file__).read_text()).body if isinstance(n,ast.FunctionDef) and n.name==name)
def native_node(n):
 if isinstance(n,ast.Assign):return any(isinstance(t,ast.Name) and t.id=='native' for t in n.targets)
 if isinstance(n,ast.If):
  return (isinstance(n.test,ast.Compare) and isinstance(n.test.left,ast.Name) and n.test.left.id=='native') or (isinstance(n.test,ast.Call) and isinstance(n.test.func,ast.Name) and n.test.func.id=='_native_ilu_solve')
 return False
node.body=[n for n in node.body if not native_node(n)]
code=compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),module.__file__,'exec')
old=types.FunctionType(next(c for c in code.co_consts if type(c)is types.CodeType),module.__dict__)
function=old if arm=='old' else getattr(module,name)
q=lambda rows:[[s.QQ(v) for v in row] for row in rows]
if route=='hnf':
 a=DomainMatrix.from_list([[1,2],[3,4]],s.ZZ);invoke=lambda:function(a);values=lambda out:[[int(v) for v in row] for row in out.to_list()]
elif route=='lu':
 a=q([[2,1],[4,3]]);invoke=lambda:function(a);values=lambda out:{'matrix':[[str(v) for v in row] for row in a],'result':out}
else:
 x=[[None],[None]];l=q([[1,0],[2,1]]);u=q([[2,1],[0,1]]);b=q([[4],[10]])
 invoke=lambda:function(x,l,u,[],b);values=lambda out:{'matrix':[[str(v) for v in row] for row in x],'result':out}
events=[]
def audit(event,args):
 if event=='object.__getattr__' and len(args)>1 and args[1] in ('__code__','__defaults__','__kwdefaults__','__closure__','__class__'):
  events.append(args[1]);raise RuntimeError('metadata audit refusal: '+args[1])
sys.addaudithook(audit)
try:out=invoke()
except Exception as error:result={'status':'raised','type':type(error).__name__,'message':str(error)}
else:result={'status':'returned','value':values(out)}
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
print(json.dumps({'route':route,'arm':arm,'events':events,'result':result,'reference_imported':False}))
