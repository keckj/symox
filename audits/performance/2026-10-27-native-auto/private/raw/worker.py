import hashlib,json,os,statistics,sys,time
from fractions import Fraction
from pathlib import Path
import symox as s
from symox.polys.matrices import DomainMatrix
import importlib
r=importlib.import_module('symox.polys.matrices.rref')
arm=sys.argv[1]
if arm=='private-disabled':r._native_choose_method=lambda *args:None
specs=[('QQ-CD',n,'dense',operation) for n in [3,10,20] for operation in ['rref','rref_den']]
specs +=[(kind,10,fmt,operation) for kind,fmt in [('QQ-CD','sparse'),('QQ-FF','dense'),('ZZ-FF','dense'),('QQ-GJ','dense')] for operation in ['rref','rref_den']]
def fraction(e):
    if type(e) is int:return Fraction(e)
    return Fraction(int(e.numerator),int(e.denominator))
results=[]
for kind,n,fmt,operation in specs:
    matrices=[];inputs=[]
    for sample in range(8):
        k=s.ZZ if kind=='ZZ-FF' else s.QQ
        factor=(1<<211)+1 if kind=='ZZ-FF' else ((1<<80)+1 if kind=='QQ-FF' else 1)
        rows=[];exact=[]
        for i in range(n):
            values=[((i*13+j*17+sample*19)%23-11 or 1)*factor for j in range(n)]
            values[i]=sum(abs(v) for j,v in enumerate(values) if j!=i)+factor
            den=1 if kind=='ZZ-FF' else ((1<<55) if kind=='QQ-FF' else ((1<<400) if kind=='QQ-GJ' else [2,3,5,7][(i+sample)%4]))
            assert abs(Fraction(values[i],den))>sum(abs(Fraction(v,den)) for j,v in enumerate(values) if j!=i)
            exact.append([Fraction(v,den) for v in values])
            rows.append([k(v) if kind=='ZZ-FF' else k(v,den) for v in values])
        data=rows if fmt=='dense' else {i:{j:v for j,v in enumerate(row)} for i,row in enumerate(rows)}
        matrices.append(DomainMatrix(data,(n,n),k));inputs.append(exact)
    expected=kind.split('-')[1]
    choices=[r._dm_rref_choose_method(a,'auto')[0] for a in matrices]
    assert choices==[expected]*8,(kind,choices)
    digests=[];samples=[]
    for phase in ['cold','warm']:
        elapsed=[]
        for a,original in zip(matrices,inputs):
            start=time.perf_counter_ns();out=getattr(a,operation)();elapsed.append(time.perf_counter_ns()-start)
            reduced=out[0];den=Fraction(1) if operation=='rref' else fraction(out[1]);pivots=out[-1]
            assert den and pivots==tuple(range(n)) and reduced.shape==(n,n) and reduced.rep.fmt==fmt
            values=[[fraction(e) for e in row] for row in reduced.to_list()]
            assert all(e==(den if i==j else 0) for i,row in enumerate(values) for j,e in enumerate(row))
            assert [[fraction(e) for e in row] for row in a.to_list()]==original
            # Strict row dominance proves invertibility. Thus den*I is the
            # complete exact RREF certificate, not only a residual witness.
            payload=[str(reduced.domain),fmt,list(pivots),str(den),[[str(e) for e in row] for row in values]]
            digests.append(hashlib.sha256(json.dumps(payload).encode()).hexdigest())
        samples.append({'phase':phase,'nanoseconds':elapsed,'median_ns':statistics.median(elapsed)})
    assert digests[:8]==digests[8:]
    results.append({'case':f'{kind}-{n}-{fmt}-{operation}','choices':choices,'digests':digests,'timings':samples})
core=Path(s._core.__file__)
print(json.dumps({'arm':arm,'pid':os.getpid(),'affinity':sorted(os.sched_getaffinity(0)),'python':sys.version,'package':s.__file__,'native_sha256':hashlib.sha256(core.read_bytes()).hexdigest(),'certificates':len(results)*16,'cases':results},sort_keys=True))
