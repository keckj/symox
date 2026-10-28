"""Release-only DDM timings; exact product certificates outside every timer."""
import hashlib,json,statistics,sys,time
from pathlib import Path
class BlockReference:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import during native measurement')
block=BlockReference();sys.meta_path.insert(0,block)
import symox as S
from symox.polys.domains import ZZ,QQ
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.dense import ddm_imatmul

result={'native_sha256':hashlib.sha256(Path(S._core.__file__).read_bytes()).hexdigest(),
        'method':'18 distinct preconstructed inputs per case, 3 batches of 6; public DDM allocates output inside call; raw accumulations preallocate nonzero targets; exact verification after all timing'}
checks=[]
for label,n,domain,raw in [('int4',4,int,True),('int32',32,int,True),
                           ('ZZ8',8,ZZ,False),('ZZ32',32,ZZ,False),
                           ('QQ8',8,QQ,False),('QQ24',24,QQ,False),('QQ_accumulate16',16,QQ,True)]:
    operations=[];metadata=[]
    for variant in range(18):
        def element(p,q=1):return p if domain is int else domain(p,q)
        b=[[element((3*i+2*j+variant)%13-6,1 if domain is not QQ else j%5+1) for j in range(n)] for i in range(n)]
        c=[[element((2*i-3*j+variant)%11-5,1 if domain is not QQ else i%7+1) for j in range(n)] for i in range(n)]
        if raw:
            a=[[element(i-j+variant,1 if domain is not QQ else 3) for j in range(n)] for i in range(n)]
            initial=[row[:] for row in a];identities=[id(row) for row in a]
            def operation(a=a,b=b,c=c):ddm_imatmul(a,b,c);return a
        else:
            b=DDM(b,(n,n),domain);c=DDM(c,(n,n),domain)
            a=None;initial=None;identities=None
            operation=lambda b=b,c=c:b.matmul(c)
        operations.append(operation)
        metadata.append((b,c,initial,identities,domain))
    samples=[];outputs=[]
    for start in range(0,18,6):
        begin=time.perf_counter();batch=[operation() for operation in operations[start:start+6]]
        samples.append((time.perf_counter()-begin)/6);outputs.extend(batch)
    result[label]={'samples_seconds_per_call':samples,'median_seconds_per_call':statistics.median(samples),
                   'values_digest':hashlib.sha256('\n'.join(str([list(row) for row in value]) for value in outputs).encode()).hexdigest()}
    checks.extend((*row,value) for row,value in zip(metadata,outputs))
sys.meta_path.remove(block)
import sympy as sp

def matrix(rows):return sp.Matrix([[sp.Rational(str(value)) for value in row] for row in rows])
for b,c,initial,identities,domain,actual in checks:
    expected=matrix(b)*matrix(c)
    if initial is not None:
        expected+=matrix(initial)
        assert [id(row) for row in actual]==identities
    assert matrix(actual)==expected
    if domain is int:assert all(type(value) is int for row in actual for value in row)
    else:assert all(value.parent()==domain for row in actual for value in row)
result['independent_exact_product_checks']=len(checks)
Path(sys.argv[1]).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
