"""Independent Fraction-only factor/rank certificates; no engine imports."""
from fractions import Fraction

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
