"""Release-only correlated-normal recurrence timings; correctness outside timing."""
import hashlib,json,statistics,sys,time
from pathlib import Path
class BlockReference:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import during native measurement')
block=BlockReference();sys.meta_path.insert(0,block)
import symox as S
from symox import stats
from symox.stats._joint import _normal_raw
result={'native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),
        'method':'18 distinct preconstructed inputs per case, 3 batches of 6; construction and correctness excluded'}
checks=[]
for label,powers,query in [('raw_pair44',(4,4),'exact'),('raw_triple642',(6,4,2),'exact'),
                          ('raw_float64',(6,4),'float'),('function_fallback42',(4,2),'function'),
                          ('public_pair42',(4,2),'public')]:
    operations=[];metadata=[]
    for i in range(18):
        dimension=len(powers)
        mean=[S.Rational(i+1,7),S.Rational(-i-2,5),S.Rational(i+3,11)][:dimension]
        cov=[[S.Integer(2),S.Rational(1,3),S.Integer(0)],
             [S.Rational(1,3),S.Integer(3),S.Rational(1,5)],
             [S.Integer(0),S.Rational(1,5),S.Integer(4)]]
        cov=[row[:dimension] for row in cov[:dimension]]
        if query=='float':
            mean=[S.N(value,80) for value in mean]
            cov=[[S.N(value,80) for value in row] for row in cov]
        if query=='function':mean[0]=S.sin(S.Symbol('normal_fallback_'+str(i),real=True))
        mean_matrix=S.ImmutableMatrix(mean);cov_matrix=S.ImmutableMatrix(cov)
        if query=='public':
            X=stats.MultivariateNormal('normal_bench_'+str(i),mean,cov)
            expression=S.Mul(*(X[j]**n for j,n in enumerate(powers)))
            operation=lambda expression=expression:stats.E(expression)
        else:operation=lambda powers=powers,mean=mean_matrix,cov=cov_matrix:_normal_raw(powers,mean,cov)
        operations.append(operation);metadata.append((powers,mean,cov,query))
    samples=[];values=[]
    for start in range(0,18,6):
        begin=time.perf_counter();batch=[operation() for operation in operations[start:start+6]]
        samples.append((time.perf_counter()-begin)/6);values.extend(batch)
    result[label]={'samples_seconds_per_call':samples,'median_seconds_per_call':statistics.median(samples),
                   'values_digest':hashlib.sha256('\n'.join(map(str,values)).encode()).hexdigest()}
    checks.extend((*row,value) for row,value in zip(metadata,values))
sys.meta_path.remove(block)
import sympy as p
import mpmath as mp

def coefficient(powers,mean,cov):
    dimension=len(powers);zero=(0,)*dimension;poly={zero:p.S.One};factors=[]
    for i,mu in enumerate(mean):
        factors.append({tuple(n if j==i else 0 for j in range(dimension)):mu**n/p.factorial(n) for n in range(powers[i]+1)})
    for i in range(dimension):
        for j in range(i,dimension):
            factor={};a=cov[i][j]/(2 if i==j else 1)
            for n in range((powers[i]//2 if i==j else min(powers[i],powers[j]))+1):
                exponent=list(zero);exponent[i]+=n;exponent[j]+=n
                factor[tuple(exponent)]=a**n/p.factorial(n)
            factors.append(factor)
    for factor in factors:
        result={}
        for left,a in poly.items():
            for right,b in factor.items():
                exponent=tuple(x+y for x,y in zip(left,right))
                if all(n<=maximum for n,maximum in zip(exponent,powers)):
                    result[exponent]=result.get(exponent,0)+a*b
        poly=result
    return p.expand(poly.get(tuple(powers),0)*p.prod(p.factorial(n) for n in powers))

for powers,mean,cov,query,value in checks:
    mean=[S.to_sympy(v) for v in mean];cov=[[S.to_sympy(v) for v in row] for row in cov]
    if query=='float':
        mean=list(map(p.Rational,mean));cov=[list(map(p.Rational,row)) for row in cov]
    expected=coefficient(powers,mean,cov)
    if query=='float':
        with mp.workdps(130):
            reference=mp.mpf(str(expected.evalf(130)));actual=mp.mpf(str(value))
            assert abs(actual-reference)<=abs(reference)*mp.mpf('1e-77')
    else:assert p.expand(S.to_sympy(value)-expected)==0
result['independent_checks']=len(checks)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
