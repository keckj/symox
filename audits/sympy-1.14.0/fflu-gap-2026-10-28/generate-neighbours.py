from pathlib import Path
from itertools import permutations
import json
folder=Path(__file__).resolve().parent
cases=json.loads((folder/'fixtures-witness.json').read_text())['cases'];original=cases[:]
columns=[(0,1,2,3),(3,2,1,0),(1,2,3,0),(2,3,0,1),(1,0,2,3),(0,1,3,2)]
for case in original[:2]:
 rows=case['rows'];m,n=case['shape']
 for rp in permutations(range(m)):
  for cp in columns:
   cases.append({'id':case['id']+'-r'+''.join(map(str,rp))+'-c'+''.join(map(str,cp)),'group':'permuted-witness','shape':[m,n],'rows':[[rows[i][j] for j in cp] for i in rp]})
for m in range(5):
 for n in range(6):
  for r in range(min(m,n)+1):
   b=[[int(i==k) if i<r else (i+1)*(k+1)+(-1)**(i+k) for k in range(r)] for i in range(m)]
   c=[[int(k==j) if j<r else (k+2)*(j+1)+(-1)**(j+k) for j in range(n)] for k in range(r)]
   a=[[sum(b[i][k]*c[k][j] for k in range(r)) for j in range(n)] for i in range(m)]
   for variant in range(2):
    rp=list(range(m));cp=list(range(n))
    if variant:
     rp=rp[1:]+rp[:1];cp=cp[1:]+cp[:1]
    cases.append({'id':f'rank-m{m}-n{n}-r{r}-v{variant}','group':'all-small-shapes-ranks','shape':[m,n],'expected_rank':r,'rows':[[a[i][j] for j in cp] for i in rp]})
# Additional zero gaps are required independently of random dependent columns.
for m,n in [(1,4),(2,5),(3,5),(4,5)]:
 for r in range(1,min(m,n)+1):
  a=[[int(j==n-r+i)*(i+2) if i<r else (i+1)*(j+1) if j>=n-r else 0 for j in range(n)] for i in range(m)]
  cases.append({'id':f'delayed-m{m}-n{n}-r{r}','group':'late-pivots','shape':[m,n],'expected_rank':r,'rows':a})
(folder/'fixtures-neighbours.json').write_text(json.dumps({'scope':'ZZ/QQ exact public/private shape-rank and pivot permutation investigation; no timings','cases':cases},indent=2)+'\n')
print(len(cases))
