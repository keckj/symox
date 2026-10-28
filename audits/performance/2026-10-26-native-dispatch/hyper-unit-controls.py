import hashlib,json,statistics,sys,time
from pathlib import Path
import symox as s

def family(name,count):
    for i in range(count):
        a=s.Rational(1,3)+s.Rational(i,1001)
        if name=='gamma': yield s.gamma(a)
        elif name=='gauss': yield s.hyper([a,s.Rational(2,3)],[2],1)
        else: yield s.hyper([-1-i%9,a],[s.Rational(3,2)],1)
result={'native_sha256':hashlib.sha256(Path(s._core.__file__).read_bytes()).hexdigest(),
        'method':'5 batches x 12 distinct preconstructed inputs, evalf(30), construction excluded'}
for name in ['gamma','gauss','terminating']:
    exprs=list(family(name,60)); samples=[];values=[]
    for start in range(0,60,12):
        t=time.perf_counter()
        batch=[e.evalf(30) for e in exprs[start:start+12]]
        samples.append((time.perf_counter()-t)/12);values.extend(batch)
    result[name]={'median_seconds_per_eval':statistics.median(samples), 'samples':samples,
                  'finite':sum(v.is_finite is True for v in values),
                  'values_digest':hashlib.sha256('\n'.join(map(str,values)).encode()).hexdigest()}
Path(sys.argv[1]).write_text(json.dumps(result,indent=2))
