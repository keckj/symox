import marshal,hashlib,json
real=marshal.dumps;events=[]
def changed(*args,**kwargs):
    events.append('dumps')
    return real(*args,**kwargs)
marshal.dumps=changed
import importlib,symox as s
from pathlib import Path
r=importlib.import_module('symox.integrals.risch');rt=importlib.import_module('symox.integrals.rationaltools');helper=importlib.import_module('symox.integrals._risch_native')
t=s.Symbol('probe_marshal_t');p=s.Poly(t**4+t,t);de=r.DifferentialExtension(extension={'D':[s.Poly(t*t+1,t)]})
importlib.reload(rt)
events.clear();actual=r.polynomial_reduce(p,de);observed=events[:]
r._native_polynomial_reduce=lambda *args:None
events.clear();expected=r.polynomial_reduce(p,de);old=events[:]
marshal.dumps=real
print(json.dumps({'same_values':actual==expected,'actual_callbacks':observed,'old_callbacks':old,'captured_foreign_marshal':helper._CODE is changed,'native_sha256':hashlib.sha256(Path(s._core.__file__).read_bytes()).hexdigest()}))
