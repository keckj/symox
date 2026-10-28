"""Release-only public/raw Poisson measurements, validated after all timings."""
import hashlib,json,statistics,sys,time
from pathlib import Path
class BlockReference:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import during native measurement')
block=BlockReference();sys.meta_path.insert(0,block)
import symox as S
from symox import stats
legacy='--legacy' in sys.argv
result={'native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),
        'legacy_control':legacy,
        'method':'18 distinct preconstructed inputs per case; 3 batches of 6; constructor and validation costs excluded; no completed value reused'}
validations=[]
for label,order,query in [('raw12',12,'raw'),('raw36',36,'raw'),('raw72',72,'raw'),
                          ('public12',12,'public'),('independent_mix4',4,'mixed'),
                          ('normal_control12',12,'normal')]:
    operations=[]
    for i in range(18):
        rate=S.Symbol('benchmark_rate_'+label+'_'+str(i),positive=True)
        X=stats.Poisson('benchmark_P_'+label+'_'+str(i),rate)
        model=stats.density(X)
        if legacy:
            model._moment=lambda n,rate=rate:S.Add(*(S.stirling(n,j,kind=2)*rate**j for j in range(int(n)+1))) if n.is_Integer else None
        n=S.Integer(order)
        if query=='raw':operation=lambda model=model,n=n:model.moment(n)
        elif query=='public':
            expression=X**n
            operation=lambda expression=expression:stats.E(expression)
        else:
            N=stats.Normal('benchmark_N_'+label+'_'+str(i),0,1)
            expression=(X+N)**n if query=='mixed' else N**n
            operation=lambda expression=expression:stats.E(expression)
        operations.append(operation)
        validations.append((label,rate,order,query))
    samples=[];values=[]
    for start in range(0,18,6):
        begin=time.perf_counter();batch=[operation() for operation in operations[start:start+6]]
        samples.append((time.perf_counter()-begin)/6);values.extend(batch)
    result[label]={'samples_seconds_per_call':samples,'median_seconds_per_call':statistics.median(samples),
                   'values_digest':hashlib.sha256('\n'.join(map(str,values)).encode()).hexdigest()}
    for row,value in zip(validations[-18:],values):
        validations[validations.index(row)]=(*row,value)
sys.meta_path.remove(block)
import sympy as p
for label,rate,order,query,value in validations:
    lam=S.to_sympy(rate)
    if query in ('raw','public'):expected=p.bell(order,lam)
    elif query=='normal':expected=p.factorial2(order-1)
    else:expected=sum(p.binomial(order,j)*p.bell(order-j,lam)*p.factorial2(j-1) for j in range(0,order+1,2))
    assert p.expand(S.to_sympy(value)-expected)==0,(label,value,expected)
result['independent_exact_checks']=len(validations)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
