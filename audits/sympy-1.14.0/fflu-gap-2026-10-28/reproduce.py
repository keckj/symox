"""Exact FFLU inspection; no timers, builds or engine modifications."""
import argparse,hashlib,importlib,json
from fractions import Fraction
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--library',choices=['sympy','symox'],required=True);p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--witness-only',action='store_true');args=p.parse_args()
s=importlib.import_module(args.library)
DDM=importlib.import_module(args.library+'.polys.matrices.ddm').DDM
DomainMatrix=importlib.import_module(args.library+'.polys.matrices.domainmatrix').DomainMatrix
if args.library=='sympy':assert s.__version__=='1.14.0'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
root=Path(s.__file__).resolve().parent
sources={str(p.relative_to(root)):sha(p) for p in [root/'polys/matrices/ddm.py',root/'polys/matrices/domainmatrix.py',root/'polys/matrices/sdm.py']}
native=sha(Path(s._core.__file__)) if args.library=='symox' else None

def values(a):return [[Fraction(str(v)) for v in row] for row in (a.to_list() if isinstance(a,DomainMatrix) else a)]
def encoded(a):return [[str(v) for v in row] for row in a]
def rank(a,m,n):
 a=[row[:] for row in a];r=0
 for j in range(n):
  pivot=next((i for i in range(r,m) if a[i][j]),None)
  if pivot is None:continue
  a[r],a[pivot]=a[pivot],a[r];d=a[r][j]
  for k in range(j,n):a[r][k]/=d
  for i in range(r+1,m):
   q=a[i][j]
   for k in range(j,n):a[i][k]-=q*a[r][k]
  r+=1
  if r==m:break
 return r

def certify(a,m,n,p,l,d,u):
 required={'P':[m,m],'L':[m,m],'D':[m,m],'U':[m,n]}
 shapes=all(len(v)==required[k][0] and all(len(row)==required[k][1] for row in v) for k,v in [('P',p),('L',l),('D',d),('U',u)])
 if not shapes:return {'shapes':False,'complete_certificate':False}
 permutation=all(v in (0,1) for row in p for v in row) and all(sum(row)==1 for row in p) and all(sum(p[i][j] for i in range(m))==1 for j in range(m))
 lower=all(l[i][j]==0 for i in range(m) for j in range(i+1,m));diagonal=all(d[i][j]==0 for i in range(m) for j in range(m) if i!=j)
 invertible_d=all(d[i][i]!=0 for i in range(m));invertible_l=all(l[i][i]!=0 for i in range(m))
 upper=all(u[i][j]==0 for i in range(m) for j in range(min(i,n)))
 pivots=[];echelon=True
 for i,row in enumerate(u):
  col=next((j for j,v in enumerate(row) if v),None)
  if col is None:continue
  if i!=len(pivots):echelon=False
  if any(u[k][col] for k in range(i+1,m)):echelon=False
  pivots.append(col)
 echelon=echelon and pivots==sorted(set(pivots))
 pa=[[sum(p[i][k]*a[k][j] for k in range(m)) for j in range(n)] for i in range(m)]
 lhs=None;identity=False;wrong=[]
 if diagonal and invertible_d:
  lhs=[[sum(l[i][k]*u[k][j]/d[k][k] for k in range(m)) for j in range(n)] for i in range(m)]
  identity=pa==lhs
  wrong=[{'index':[i,j],'P_A':str(pa[i][j]),'L_Dinv_U':str(lhs[i][j])} for i in range(m) for j in range(n) if pa[i][j]!=lhs[i][j]]
 ra,ru=rank(a,m,n),rank(u,m,n);ranks=ra==ru==len(pivots)
 return {'shapes':shapes,'permutation':permutation,'L_lower':lower,'D_diagonal':diagonal,'D_invertible':invertible_d,'L_invertible':invertible_l,'U_upper':upper,'U_echelon':echelon,'pivot_columns':pivots,'rank_A':ra,'rank_U':ru,'rank_complete':ranks,'P_A_equals_L_Dinv_U':identity,'mismatches':wrong,'complete_certificate':all([shapes,permutation,lower,diagonal,invertible_d,invertible_l,upper,echelon,ranks,identity])}

def packed_factors(lu,perm,m,n):
 u=[row[:] for row in lu];p=[[Fraction(int(perm[i]==j)) for j in range(m)] for i in range(m)];l=[[Fraction(0) for j in range(m)] for i in range(m)];i=0
 for j in range(n):
  if i==m:break
  if u[i][j]:
   l[i][i]=u[i][j]
   for k in range(i+1,m):l[k][i]=u[k][j];u[k][j]=Fraction(0)
   i+=1
 for k in range(i,m):l[k][k]=Fraction(1)
 d=[[Fraction(0) for j in range(m)] for i in range(m)]
 if m:d[0][0]=l[0][0]
 for k in range(1,m):d[k][k]=l[k-1][k-1]*l[k][k]
 return p,l,d,u

fixtures=json.loads(args.fixtures.read_text());records=[]
for case in fixtures['cases']:
 if args.witness_only and case['group']!='witness':continue
 m,n=case['shape'];a=[[Fraction(v) for v in row] for row in case['rows']]
 for domain in case.get('domains',['ZZ','QQ']):
  k=getattr(s,domain)
  def cells():return [[k(v.numerator,v.denominator) if domain=='QQ' else k(v.numerator) for v in row] for row in a]
  for entry in ['private','ddm','dense','sparse']:
   raw=DDM(cells(),(m,n),k)
   source=raw if entry in ['private','ddm'] else DomainMatrix.from_rep(raw)
   if entry=='sparse':source=source.to_sparse()
   expected=values(source);ids=[]
   if entry in ['private','ddm','dense']:
    rep=source.rep if entry=='dense' else source;ids=[id(rep),[(id(row),[id(v) for v in row]) for row in rep]]
   record={'case':case['id'],'group':case['group'],'entry':entry,'domain':domain,'shape':[m,n]}
   try:out=source._fflu() if entry=='private' else source.fflu()
   except Exception as error:record.update(status='raised',error={'type':type(error).__name__,'message':str(error)})
   else:
    record['status']='returned'
    if entry=='private':
     lu,perm=out;packed=values(lu);record.update(permutation=list(perm),packed_LU=encoded(packed));factors=packed_factors(packed,perm,m,n);record['certificate']=certify(a,m,n,*factors)
     record['private_shapes_domain']=list(lu.shape)==[m,n] and lu.domain==k and sorted(perm)==list(range(m))
    else:
     assert len(out)==4
     factors=[values(v) for v in out];record['factors']={name:encoded(value) for name,value in zip('PLDU',factors)};record['factor_shapes']=[list(v.shape) for v in out];record['domain_preserved']=all(v.domain==k for v in out);record['output_format']=[v.rep.fmt if isinstance(v,DomainMatrix) else v.fmt for v in out];record['certificate']=certify(a,m,n,*factors)
   record['input_unchanged']=values(source)==expected
   if ids:
    rep=source.rep if entry=='dense' else source
    record['input_rows_cells_unchanged']=ids==[id(rep),[(id(row),[id(v) for v in row]) for row in rep]]
   assert record['input_unchanged'] and record.get('input_rows_cells_unchanged',True),record
   records.append(record)
assert sources=={p:sha(root/p) for p in sources}
if args.library=='symox':assert native==sha(Path(s._core.__file__))
summary={'library':args.library,'version':getattr(s,'__version__',None),'package':str(root),'native_sha256':native,'sources_before_after':sources,'sources_unchanged':True,'no_timing_fields':True,'calls':len(records),'returned':sum(r['status']=='returned' for r in records),'raised':sum(r['status']=='raised' for r in records),'certified':sum(r.get('certificate',{}).get('complete_certificate',False) for r in records),'returned_but_uncertified':sum(r['status']=='returned' and not r['certificate']['complete_certificate'] for r in records),'records':records}
args.out.write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='records'}))
