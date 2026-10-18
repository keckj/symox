"""Exact linear polynomial systems over native sparse coefficient rings."""
import symox as S
from .rings import sring
from .matrices.domainmatrix import DomainMatrix


class PolyNonlinearError(ValueError):pass
class RawMatrix:
    """Rectangular mutable matrix preserving raw coefficient objects."""
    def __init__(self,*args):
        if len(args)==1:
            rows=args[0].tolist() if hasattr(args[0],'tolist') else [list(row) for row in args[0]]
            self.rows=len(rows);self.cols=len(rows[0]) if rows else 0
        elif len(args)==3:
            self.rows,self.cols=map(int,args[:2]);data=args[2]
            rows=[[data(i,j) if callable(data) else data[i*self.cols+j] for j in range(self.cols)] for i in range(self.rows)]
        else:raise TypeError('RawMatrix requires rows or (nrows,ncols,elements)')
        if any(len(row)!=self.cols for row in rows):raise ValueError('ragged matrix')
        self._rows=rows
    shape=property(lambda self:(self.rows,self.cols))
    def __getitem__(self,index):
        i,j=index if isinstance(index,tuple) else divmod(index,self.cols)
        return self._rows[i][j]
    def __setitem__(self,index,value):
        i,j=index if isinstance(index,tuple) else divmod(index,self.cols);self._rows[i][j]=value
    def __iter__(self):return (v for row in self._rows for v in row)
    def tolist(self):return [list(row) for row in self._rows]
    def as_mutable(self):return S.Matrix([[v.as_expr() if hasattr(v,'as_expr') else v for v in row] for row in self._rows])
    def rref(self):
        m,p=self.as_mutable().rref();return RawMatrix(m.tolist()),p


def eqs_to_matrix(eqs_coeffs,eqs_rhs,gens,domain):
    gens=list(gens);rows=[]
    for coefficients,rhs in zip(eqs_coeffs,eqs_rhs):rows.append([domain.convert(coefficients.get(g,domain.zero)) for g in gens]+[-domain.convert(rhs)])
    return DomainMatrix(rows,(len(rows),len(gens)+1),domain)

def sympy_eqs_to_ring(eqs,symbols):
    ring,expressions=sring(list(eqs),list(symbols),field=True,extension=True)
    return expressions,ring

def solve_lin_sys(eqs,ring,_raw=True):
    symbols=list(ring.symbols);expressions=[e.as_expr() if hasattr(e,'as_expr') else S.sympify(e) for e in eqs]
    try:matrix,rhs=S.linear_eq_to_matrix(expressions,symbols)
    except ValueError as exc:raise PolyNonlinearError(str(exc)) from exc
    result=S.linsolve((matrix,rhs),symbols)
    if result==S.EmptySet:return None
    solutions=list(result)
    if not solutions:return {}
    values=list(solutions[0].args if getattr(solutions[0],'is_Tuple',False) else solutions[0]);answer={}
    for i,(variable,value) in enumerate(zip(symbols,values)):
        if value==variable:continue
        if _raw:
            converted=ring(value) if value.free_symbols&set(symbols) else ring.domain.convert(value)
            answer[ring.gens[i]]=converted
        else:answer[variable]=value
    return answer
__all__=['PolyNonlinearError','RawMatrix','eqs_to_matrix','sympy_eqs_to_ring','solve_lin_sys']
