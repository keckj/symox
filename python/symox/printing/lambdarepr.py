"""Expression strings for generated Python numerical callables."""
from .pycode import PythonCodePrinter,MpmathPrinter,ArrayPrinter
from .precedence import precedence


class LambdaPrinter(PythonCodePrinter):
    printmethod='_lambdacode'
    # SymPy's PythonCodePrinter prints a Rational as a bare ``p/q`` (Mul
    # precedence decides the parentheses).
    def _print_Rational(self,expr):return '%s/%s'%(expr.p,expr.q)
    def _print_Half(self,expr):return self._print_Rational(expr)
    def _print_Pow(self,expr,rational=False):
        # SymPy's LambdaPrinter falls back to StrPrinter._print_Pow.
        import symox as _S
        PREC=precedence(expr)
        if expr.exp is _S.S.Half and not rational:
            return 'sqrt(%s)'%self._print(expr.base)
        if expr.is_commutative:
            if -expr.exp is _S.S.Half and not rational:
                return '%s/sqrt(%s)'%(self._print(_S.S.One),self._print(expr.base))
            if expr.exp is _S.S.NegativeOne:
                return '%s/%s'%(self._print(_S.S.One),self.parenthesize(expr.base,PREC,strict=False))
        return '%s**%s'%(self.parenthesize(expr.base,PREC,strict=False),
                         self.parenthesize(expr.exp,PREC,strict=False))


    def emptyPrinter(self,expr):
        # SymPy's PythonCodePrinter._print_MatrixBase on the sympified matrix.
        if getattr(expr,'is_Matrix',False) and hasattr(expr,'tolist'):
            # symox has a single dense Matrix class; doprint sympifies it to
            # the immutable class in SymPy.
            name='ImmutableSparseMatrix' if 'Sparse' in type(expr).__name__ else 'ImmutableDenseMatrix'
            return '%s(%s)'%(name,self._print(expr.tolist()))
        return super().emptyPrinter(expr)


class NumPyPrinter(ArrayPrinter):pass
class IntervalPrinter(MpmathPrinter):_module='mpmath.iv'
def lambdarepr(expr,**settings):return LambdaPrinter(settings).doprint(expr)
__all__=['LambdaPrinter','NumPyPrinter','MpmathPrinter','PythonCodePrinter','IntervalPrinter','lambdarepr']
