"""Median over runs of first-call and repeated timings: before vs after."""
import json, statistics, sys
def load(p):
    rows=[json.loads(l) for l in open(p) if l.strip()]
    keys=set().union(*rows)
    return {k:(statistics.median(r[k][0] for r in rows), statistics.median(r[k][1] for r in rows)) for k in keys if all(k in r for r in rows)}
b,a=load(sys.argv[1]),load(sys.argv[2])
print('| case | first before (ms) | first after (ms) | repeated before (ms) | repeated after (ms) |')
print('|---|---|---|---|---|')
fmt=lambda d,k,i,p: ('%.'+str(p)+'f')%(d[k][i]*1e3) if k in d else 'n/a'
for k in a:
    print('| %s | %s | %s | %s | %s |'%(k,fmt(b,k,0,2),fmt(a,k,0,2),fmt(b,k,1,3),fmt(a,k,1,3)))
