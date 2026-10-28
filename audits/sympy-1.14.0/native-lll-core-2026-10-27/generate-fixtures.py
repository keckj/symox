from pathlib import Path
import random,json,hashlib
import sympy.polys.matrices.lll as reference_lll
from fractions import Fraction
import sympy as S
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.exceptions import DMRankError,DMShapeError,DMValueError
assert S.__version__=='1.14.0'
rng=random.Random(20261027);records=[]
def add(rows,n,delta):
 m=len(rows);reference_error=None;unpatched_changed=False;unpatched=None
 dm=DomainMatrix.from_list(rows,S.ZZ) if m else DomainMatrix.zeros((0,n),S.ZZ)
 try:
  unpatched=None
  try:
   b,t=dm.lll_transform(delta)
   unpatched={'basis':[[int(c) for c in row] for row in b.to_list()],'transform':[[int(c) for c in row] for row in t.to_list()]}
  except (AssertionError,OverflowError) as error:
   reference_error=type(error).__name__
  floor=reference_lll.mfloor
  reference_lll.mfloor=lambda q:q.numerator//q.denominator
  try:b,t=dm.lll_transform(delta)
  finally:reference_lll.mfloor=floor
  unpatched_changed=unpatched is not None and (unpatched['basis']!=b.to_list() or unpatched['transform']!=t.to_list())
  assert t.matmul(dm)==b
  if m:assert abs(t.det())==1
  stars=[];norms=[];mus=[]
  for row in b.to_list():
   original=list(map(Fraction,row));star=original[:];mu=[]
   for previous,norm in zip(stars,norms):
    coefficient=sum(a*c for a,c in zip(original,previous))/norm
    mu.append(coefficient);star=[a-coefficient*c for a,c in zip(star,previous)]
   stars.append(star);norms.append(sum(c*c for c in star));mus.append(mu)
  assert all(abs(c)<=Fraction(1,2) for mu in mus for c in mu)
  assert all(norms[i]>=(Fraction(delta.numerator,delta.denominator)-mus[i][i-1]**2)*norms[i-1] for i in range(1,m))
  status='ok';basis=[int(c) for row in b.to_list() for c in row];transform=[int(c) for row in t.to_list() for c in row]
 except (DMRankError,DMShapeError,DMValueError) as e:
  status={DMRankError:'rank',DMShapeError:'shape',DMValueError:'delta'}[type(e)];basis=[];transform=[]
 records.append(dict(rows=rows,columns=n,delta=[int(delta.numerator),int(delta.denominator)],status=status,basis=basis,transform=transform,unpatched_reference_error=reference_error,unpatched_reference_output_changed=unpatched_changed,unpatched_reference_result=unpatched))
for m in range(7):
 for extra in range(3):
  n=m+extra
  for seed in range(10):
   rows=[[rng.randrange(-99,100) for _ in range(n)] for _ in range(m)]
   for delta in [S.QQ(1,3),S.QQ(3,4),S.QQ(99,100)]:add(rows,n,delta)
for rows,n in [([[0,0]],2),([[0,0],[0,0]],2),([[1,2],[2,4]],2),([[1],[2]],1)]:
 for d in [S.QQ(1,4),S.QQ(3,4),S.QQ(1)]:add(rows,n,d)
for bits in [65,127,1024]:
 h=1<<bits
 for d in [S.QQ(1,3),S.QQ(3,4),S.QQ(99,100)]:
  add([[h+1,h-1,3],[h,2*h+1,-7],[1,0,h+3]],3,d)
folder=Path(__file__).parent
(folder/'reference-fixtures.json').write_text(json.dumps({'reference':'SymPy1.14.0','seed':20261027,'cases':records},indent=2)+'\n')
lines=[]
for r in records:
 fields=[str(len(r['rows'])),str(r['columns']),*map(str,r['delta']),r['status'],','.join(str(x) for row in r['rows'] for x in row),','.join(map(str,r['basis'])),','.join(map(str,r['transform']))]
 lines.append('\t'.join(fields))
(folder/'reference-fixtures.tsv').write_text('\n'.join(lines)+'\n')
print({'cases':len(records),'statuses':{s:sum(r['status']==s for r in records) for s in ['ok','rank','shape','delta']}})
