"""Exact linear polynomial systems over native sparse coefficient rings."""
import symox as S
from .rings import sring
from .matrices.domainmatrix import DomainMatrix


class PolyNonlinearError(ValueError):pass
class RawMatrix:
    """Rectangular mutable matrix preserving raw coefficient objects."""
    def __init__(self,*args):
        from ..utilities.exceptions import sympy_deprecation_warning
        sympy_deprecation_warning(
            """
            The RawMatrix class is deprecated. Use either DomainMatrix or
            Matrix instead.
            """,
            deprecated_since_version="1.9",
            active_deprecations_target="deprecated-rawmatrix",
        )
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
    def __eq__(self,other):
        if not hasattr(other,'tolist') or getattr(other,'shape',None)!=self.shape:return NotImplemented
        return all(a==b for r,q in zip(self._rows,other.tolist()) for a,b in zip(r,q))
    __hash__=None
    def as_mutable(self):return S.Matrix([[v.as_expr() if hasattr(v,'as_expr') else v for v in row] for row in self._rows])
    def rref(self):
        m,p=self.as_mutable().rref()
        # Internal result: no deprecation warning for the user.
        r=object.__new__(RawMatrix);r.rows,r.cols=m.shape;r._rows=m.tolist()
        return r,p


def eqs_to_matrix(eqs_coeffs,eqs_rhs,gens,domain):
    gens=list(gens);rows=[]
    for coefficients,rhs in zip(eqs_coeffs,eqs_rhs):rows.append([domain.convert(coefficients.get(g,domain.zero)) for g in gens]+[-domain.convert(rhs)])
    return DomainMatrix(rows,(len(rows),len(gens)+1),domain)

def sympy_eqs_to_ring(eqs,symbols):
    ring,expressions=sring(list(eqs),list(symbols),field=True,extension=True)
    return expressions,ring.to_domain()

def _inexact(domain):
    """An inexact ground (``RR``, ``CC`` or a fraction field over them)."""
    while domain is not None:
        if getattr(domain,'is_Exact',True) is False:return True
        try:inner=domain.domain
        except Exception:return False  # a ground domain has no coefficient domain
        if inner is domain:return False
        domain=inner
    return False


def _solve_inexact(matrix,rhs,symbols,ring,_raw):
    """The reference ``_solve_lin_sys_component`` over an inexact ground:
    the reduced row echelon form of the augmented matrix computed with the
    arithmetic of ``ring.domain`` (``RR(k)`` cancels its fractions), not
    over expressions, where float rounding keeps spurious common factors
    (``(k - c)/(k**2 - 2*c*k + c**2)`` instead of ``1/(k - c)``)."""
    K=ring.domain;n=len(symbols)
    rows=[[K.convert(matrix[i,j]) for j in range(n)]+[K.convert(rhs[i,0])] for i in range(matrix.rows)]
    prepared=_native_inexact_rref(rows,K,n+1)
    if prepared is not None:
        rows,pivots=prepared
    else:
        pivots=[];r=0
        for col in range(n+1):
            piv=next((i for i in range(r,len(rows)) if rows[i][col]),None)
            if piv is None:continue
            rows[r],rows[piv]=rows[piv],rows[r]
            inv=K.one/rows[r][col];rows[r]=[v*inv for v in rows[r]]
            for i in range(len(rows)):
                if i!=r and rows[i][col]:
                    f=rows[i][col];rows[i]=[a-f*b for a,b in zip(rows[i],rows[r])]
            pivots.append(col);r+=1
            if r==len(rows):break
    if pivots and pivots[-1]==n:return None
    answer={}
    for i,p in enumerate(pivots):
        value=K.to_sympy(rows[i][n])-S.Add(*[K.to_sympy(rows[i][j])*symbols[j] for j in range(p+1,n) if rows[i][j]])
        if _raw:
            answer[ring.gens[p]]=ring(value) if value.free_symbols&set(symbols) else ring.domain.convert(value)
        else:answer[symbols[p]]=value
    return answer


def solve_lin_sys(eqs,ring,_raw=True):
    symbols=list(ring.symbols);expressions=[e.as_expr() if hasattr(e,'as_expr') else S.sympify(e) for e in eqs]
    try:matrix,rhs=S.linear_eq_to_matrix(expressions,symbols)
    except ValueError as exc:raise PolyNonlinearError(str(exc)) from exc
    if _inexact(getattr(ring,'domain',None)):return _solve_inexact(matrix,rhs,symbols,ring,_raw)
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


from ._inexact_linear import register as _register_inexact_rref
_native_inexact_rref = _register_inexact_rref(globals())
del _register_inexact_rref
