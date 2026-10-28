"""Identical scalar operation counts, repeated and distinct trees, exact certificates."""
import hashlib,json,os,statistics,sys,time
from pathlib import Path
import symox as s
import sympy as sp

def expression(backend,kind,variant):
    x,y=backend.symbols('x y')
    if kind=='scalar':
        return sum((x+backend.Rational(i+1,3))**(i%5+2)/(y+i+1)+backend.sin(x+i) for i in range(variant+12))
    f=backend.Function('count_dag')
    e=x+variant
    for _ in range(int(kind[3:])):
        e=f(e,e)
    return e
rows=[]
for kind in ['scalar','dag8','dag12']:
    inputs=[expression(s,kind,i) for i in range(8)]
    expected=[int(sp.count_ops(expression(sp,kind,i))) for i in range(8)]
    start=time.perf_counter_ns(); first=s.count_ops(inputs[0]); cold=time.perf_counter_ns()-start
    assert int(first)==expected[0]
    results=[]; times=[]
    for _ in range(31):
        start=time.perf_counter_ns(); value=s.count_ops(inputs[0]); times.append(time.perf_counter_ns()-start); results.append(value)
    changed=[]; changed_times=[]
    for item in inputs[1:]:
        start=time.perf_counter_ns(); value=s.count_ops(item); changed_times.append(time.perf_counter_ns()-start); changed.append(int(value))
    assert all(int(v)==expected[0] for v in results)
    assert changed==expected[1:]
    rows.append(dict(case=kind,first_ns=cold,repeated_ns=times,changed_ns=changed_times,expected=expected,repeated_median_ns=statistics.median(times),changed_median_ns=statistics.median(changed_times)))
core=Path(s.__file__).parent/'_core.abi3.so'
print(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),native_sha256=hashlib.sha256(core.read_bytes()).hexdigest(),symox_file=s.__file__,python=sys.version,cpu_affinity=sorted(os.sched_getaffinity(0)),cases=rows),indent=2))
