"""Executable exact mathematics only; no engine source changes or timings."""
from pathlib import Path
from fractions import Fraction
import hashlib,json
from certificates import certify,encoded,packed_factors

def determinant(a):
 n=len(a);a=[[Fraction(v) for v in row] for row in a];result=Fraction(1)
 for j in range(n):
  pivot=next((i for i in range(j,n) if a[i][j]),None)
  if pivot is None:return Fraction(0)
  if pivot!=j:a[j],a[pivot]=a[pivot],a[j];result=-result
  p=a[j][j];result*=p
  for i in range(j+1,n):
   q=a[i][j]/p
   for k in range(j+1,n):a[i][k]-=q*a[j][k]
 return result

def corrected(a,m,n,integer=False):
 lu=[[int(v) if integer else Fraction(v) for v in row] for row in a];perm=list(range(m));r=0;previous=1;columns=[];proof=[];divisions=0
 for j in range(n):
  if r==m:break
  pivot=next((i for i in range(r,m) if lu[i][j]),None)
  if pivot is None:continue
  lu[r],lu[pivot]=lu[pivot],lu[r];perm[r],perm[pivot]=perm[pivot],perm[r]
  p=lu[r][j]
  assert p==determinant([[a[i][k] for k in columns+[j]] for i in perm[:r+1]])
  detail={'rank':r,'column':j,'pivot':str(p),'previous':str(previous),'permutation':perm[:],'division_count':0}
  for i in range(r+1,m):
   multiplier=lu[i][j]
   for k in range(j+1,n):
    numerator=p*lu[i][k]-lu[r][k]*multiplier
    assert previous!=0
    if integer:
     q,rem=divmod(numerator,previous);assert rem==0
    else:q=numerator/previous
    assert previous*q==numerator
    selected_rows=perm[:r+1]+[perm[i]];selected_columns=columns+[j,k]
    assert q==determinant([[a[ii][kk] for kk in selected_columns] for ii in selected_rows])
    lu[i][k]=q;divisions+=1;detail['division_count']+=1
   # The original multiplier object/value stays in the packed pivot column.
   assert lu[i][j]==multiplier
  columns.append(j);previous=p;r+=1;proof.append(detail)
 factors=packed_factors([[Fraction(v) for v in row] for row in lu],perm,m,n)
 certificate=certify([[Fraction(v) for v in row] for row in a],m,n,*factors)
 assert certificate['complete_certificate'] and certificate['pivot_columns']==columns
 return lu,perm,factors,certificate,proof,divisions

if __name__=='__main__':
 folder=Path(__file__).resolve().parent;paths=[folder/'fixtures-neighbours.json',folder/'fixtures-rational.json'];cases=[case for p in paths for case in json.loads(p.read_text())['cases']];records=[]
 for case in cases:
  m,n=case['shape'];a=[[Fraction(v) for v in row] for row in case['rows']]
  for domain in case.get('domains',['ZZ','QQ']):
   lu,perm,factors,certificate,trace,divisions=corrected(a,m,n,domain=='ZZ')
   if 'expected_rank' in case:assert certificate['rank_A']==case['expected_rank']
   records.append({'case':case['id'],'domain':domain,'shape':[m,n],'factors':{name:encoded(value) for name,value in zip('PLDU',factors)},'packed_LU':encoded(lu),'permutation':perm,'certificate':certificate,'pivot_minor_trace':trace,'exact_divisions_certified':divisions})
 result={'status':'All exact mathematical prototype certificates passed; not a production port','case_count':len(cases),'domain_cases':len(records),'integer_or_rational_divisions_certified':sum(r['exact_divisions_certified'] for r in records),'every_pivot_and_update_matches_independent_minor_determinant':True,'all_factor_identities_ranks_shapes_certified':True,'fixture_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'records':records}
 (folder/'correction-math.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'}))
