"""Rational box/knapsack programs with exact primal/dual certificates."""
import hashlib,json,os,statistics,sys,time
from pathlib import Path
import symox as s
from symox.solvers.simplex import linprog

def problem(n,variant):
    capacities=[s.Rational(i%5+2+variant,3) for i in range(n)]
    budget=sum(capacities)*s.Rational(2,3)
    rows=[[s.Integer(i==j) for j in range(n)] for i in range(n)]+[[s.Integer(1)]*n]
    return [s.Integer(-1)]*n,rows,capacities+[budget]

def certify(problem,result):
    cost,rows,bounds=problem;optimum,x=result
    assert len(x)==len(cost) and all(value>=0 for value in x)
    assert all(sum(a*v for a,v in zip(row,x))<=b for row,b in zip(rows,bounds))
    assert sum(c*v for c,v in zip(cost,x))==optimum
    # Dual multiplier -1 for the final sum constraint, zero for every
    # individual bound. Feasibility and equality certify the exact optimum.
    assert optimum==-bounds[-1]

records=[]
for n in (10,30,60):
    inputs=[problem(n,v) for v in range(4)]
    start=time.perf_counter_ns();first=linprog(*inputs[0]);cold=time.perf_counter_ns()-start
    certify(inputs[0],first)
    repeated=[]
    for _ in range(7):
        start=time.perf_counter_ns();result=linprog(*inputs[0]);repeated.append(time.perf_counter_ns()-start)
        certify(inputs[0],result)
    changed=[]
    for item in inputs[1:]:
        start=time.perf_counter_ns();result=linprog(*item);changed.append(time.perf_counter_ns()-start)
        certify(item,result)
    records.append(dict(n=n,first_ns=cold,repeated_ns=repeated,changed_ns=changed,repeated_median_ns=statistics.median(repeated),changed_median_ns=statistics.median(changed),distinct_certified_inputs=len(inputs)))
print(json.dumps(dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),native_sha256=hashlib.sha256((Path(s.__file__).parent/'_core.abi3.so').read_bytes()).hexdigest(),symox_file=s.__file__,python=sys.version,cpu_affinity=sorted(os.sched_getaffinity(0)),cases=records),indent=2))
