import hashlib,json,statistics,sys,time
from pathlib import Path
import mpmath as mp
import symox as s

result={'native_sha256':hashlib.sha256(Path(s._core.__file__).read_bytes()).hexdigest(),
        'method':'24 independently allocated preconstructed expressions per family/precision; 4 batches x 6 evalf calls; construction and independent high-precision validation excluded'}
for digits in [30,80]:
    for reflected in [False,True]:
        for distinct in [False,True]:
            key=('%sd_%s_%s'%(digits,'reflected' if reflected else 'mid','distinct' if distinct else 'repeated'))
            numer,denom=(1,3) if reflected else (11,7)
            args=[s.Rational(numer,denom)+(s.Rational(i,1001) if distinct else 0) for i in range(24)]
            exprs=[s.gamma(x,evaluate=False) for x in args]
            assert len({e._dag_id for e in exprs})==len(exprs)
            samples=[];values=[];cold=[]
            for start in range(0,len(exprs),6):
                t=time.perf_counter();batch=[]
                for e in exprs[start:start+6]:
                    q=time.perf_counter();batch.append(e.evalf(digits));cold.append(time.perf_counter()-q)
                samples.append((time.perf_counter()-t)/6);values.extend(batch)
            with mp.workdps(digits+80):
                errors=[]
                for i,v in enumerate(values):
                    argument=mp.mpf(numer)/denom+(mp.mpf(i)/1001 if distinct else 0)
                    reference=mp.gamma(argument);actual=mp.mpf(str(v))
                    errors.append(abs((actual-reference)/reference))
                assert max(errors)<mp.mpf(10)**(1-digits)
                error=str(max(errors))
            result[key]={'seconds_per_eval_samples':samples,'median_seconds_per_eval':statistics.median(samples),
                         'first_eval_seconds':cold[0],'max_relative_error':error,
                         'values_digest':hashlib.sha256('\n'.join(map(str,values)).encode()).hexdigest()}
Path(sys.argv[1]).write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
