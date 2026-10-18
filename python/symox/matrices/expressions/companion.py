"""Companion matrices of monic univariate polynomials."""
from ... import Poly, Matrix, MatrixExpr

class CompanionMatrix:
    _symox_matrixexpr=True
    def __init__(self,poly):
        self.poly=poly if isinstance(poly,Poly) else Poly(poly)
        if len(self.poly.gens)!=1 or self.poly.degree()<1 or self.poly.LC()!=1:
            raise ValueError('a companion matrix requires a monic univariate polynomial of positive degree')
    args=property(lambda self:(self.poly,))
    shape=property(lambda self:(self.poly.degree(),)*2)
    rows=cols=property(lambda self:self.poly.degree())
    free_symbols=property(lambda self:set().union(*(c.free_symbols for c in self.poly.all_coeffs())))
    func=property(lambda self:type(self))
    def as_explicit(self):
        coefficients=list(reversed(self.poly.all_coeffs()));n=self.rows
        return Matrix(n,n,lambda i,j:-coefficients[i] if j==n-1 else 1 if i==j+1 else 0)
    def _symox_(self):return MatrixExpr(self.as_explicit())
    def __getitem__(self,index):return self.as_explicit()[index]
    def __mul__(self,other):return self._symox_()*other
    def __rmul__(self,other):return other*self._symox_()
    def __add__(self,other):return self._symox_()+other
    __radd__=__add__
    def __pow__(self,n):return self._symox_()**n
    def det(self):return self.as_explicit().det()
    def charpoly(self,*args):return self.as_explicit().charpoly(*args)
    def __repr__(self):return 'CompanionMatrix(%s)'%self.poly
    def __eq__(self,other):return isinstance(other,CompanionMatrix) and self.poly==other.poly
    def __hash__(self):return hash((type(self),self.poly))

__all__=['CompanionMatrix']
