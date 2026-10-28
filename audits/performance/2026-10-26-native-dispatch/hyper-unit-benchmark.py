import hashlib,json,statistics,sys,time
from pathlib import Path
import symox as s

def family(name, count):
    for i in range(count):
        if name == 'gauss':
            a=s.Rational(1,3)+s.Rational(i,1001)
            yield s.hyper([a,s.Rational(2,3)],[2],1)
        elif name == 'balanced':
            a=s.Rational(1,3)+s.Rational(i,1001)
            b=s.Rational(1,4)
            yield s.hyper([a,b],[a+b],1)
        elif name == 'one_f_zero':
            a=s.Rational(1,3)+s.Rational(i,1001)
            yield s.hyper([a],[],1)
        else:
            a=-s.Rational(1,3)+s.Rational(i,1001)
            b=-s.Rational(1,4)
            yield s.hyper([a,b],[a+b],1)

result={'native_sha256':hashlib.sha256(Path(s._core.__file__).read_bytes()).hexdigest(),
        'method':'5 batches, 12 distinct preconstructed inputs each, evalf(30); no input repeated, construction excluded'}
for name in ['one_f_zero','balanced','signed_balanced','gauss']:
    exprs=list(family(name,60))
    samples=[]; values=[]; errors=[]
    for start in range(0,60,12):
        t=time.perf_counter()
        batch=[]
        for e in exprs[start:start+12]:
            try: batch.append(e.evalf(30))
            except Exception as error:
                errors.append(type(error).__name__+': '+str(error))
        samples.append((time.perf_counter()-t)/12)
        values.extend(batch)
    result[name]={'seconds_per_eval_samples':samples,'median_seconds_per_eval':statistics.median(samples),
                  'output_examples':[str(v) for v in values[:2]],
                  'counts':{'positive_infinity':sum(v==s.oo for v in values),'negative_infinity':sum(v==-s.oo for v in values),
                            'held_hyper':sum(v.has(s.hyper) for v in values),'finite':sum(v.is_finite is True for v in values)},'errors':errors}
Path(sys.argv[1]).write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
