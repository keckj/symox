import sys
mod = sys.argv[1]
S = __import__(mod)
from importlib import import_module
kr = import_module(mod+'.matrices.expressions.kronecker')
KP = kr.KroneckerProduct; ck = kr.combine_kronecker
n,m,x,k = S.symbols('n m x k')
Z=S.MatrixSymbol('Z',n,n); W=S.MatrixSymbol('W',m,m); A=S.MatrixSymbol('A',n,m); B=S.MatrixSymbol('B',n,m)
e = KP(W,Z)*KP(W.I,Z.I)
print(1, type(e).__name__, e, e.is_Identity)
print(2, repr(KP(S.eye(2),S.eye(3))), type(KP(S.eye(2),S.eye(3))).__name__)
print(3, KP(S.I*A,B).conjugate())
X=S.MatrixSymbol('X',n,n);Y=S.MatrixSymbol('Y',n,n)
print(4, ck(KP(X,Y)**x), ck(x*KP(X,Y)**2))
print(5, KP(X+Y,Y+Z).expand(kroneckerproduct=True))
i,j,o,p=S.symbols('i j o p')
A2=S.MatrixSymbol('A',n,m);B2=S.MatrixSymbol('B',o,p)
print(6, KP(A2,B2)._entry(i,j) if hasattr(KP(A2,B2),'_entry') else 'noentry')
print(7, KP(A2,B2)[i,j])
print(8, (KP(W,Z)*KP(W.I,Z.I)).doit())
