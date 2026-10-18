"""Linear-system expression preserving its solve intent until code generation."""
import symox as S

class MatrixSolve:
    _symox_matrixexpr=True
    def __init__(self,A,b):
        if A.rows!=A.cols or A.rows!=b.rows:raise S.ShapeError('linear-system shapes do not match')
        self.A,self.b=A,b;self.args=(A,b)
    @property
    def shape(self):return self.b.shape
    rows=property(lambda self:self.b.rows)
    cols=property(lambda self:self.b.cols)
    free_symbols=property(lambda self:self.A.free_symbols|self.b.free_symbols)
    @property
    def func(self):return type(self)
    def _symox_(self):return S.MatrixExpr(self.A).inv()*S.MatrixExpr(self.b)
    def doit(self,**hints):
        if isinstance(self.A,S.Matrix) and isinstance(self.b,S.Matrix):return self.A.LUsolve(self.b)
        return self
    def as_explicit(self):return self.doit() if isinstance(self.A,S.Matrix) else self._symox_().as_explicit()
    def rewrite(self,*args,**kwargs):return self._symox_()
    def __getitem__(self,key):return self._symox_()[key]
    def __mul__(self,value):return self._symox_()*value
    def __rmul__(self,value):return value*self._symox_()
    def __add__(self,value):return self._symox_()+value
    __radd__=__add__
    def __sub__(self,value):return self._symox_()-value
    def subs(self,*args,**kwargs):return type(self)(self.A.subs(*args,**kwargs),self.b.subs(*args,**kwargs))
    def __repr__(self):return 'MatrixSolve(%s, %s)'%(self.A,self.b)
    def __eq__(self,other):return isinstance(other,MatrixSolve) and self.args==other.args
    def __hash__(self):return hash((type(self),self.args))
__all__=['MatrixSolve']
