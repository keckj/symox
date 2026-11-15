"""Expression strings for generated Python numerical callables."""
from .pycode import PythonCodePrinter,MpmathPrinter,ArrayPrinter
from .numpy import NumPyPrinter
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


_NUMEXPR_PYTHON_KNOWN=('Abs','Max','Min','Sqrt','acos','acosh','asin','asinh','atan','atan2','atanh','ceiling',
    'cos','cosh','erf','erfc','exp','expm1','factorial','floor','gamma','hypot','isinf','isnan','ln','log',
    'log10','log1p','log2','loggamma','sin','sinh','tan','tanh')


class NumExprPrinter(LambdaPrinter):
    """Port of SymPy's ``NumExprPrinter``: the expression is printed inside
    ``numexpr.evaluate('...', truediv=True)``."""
    # key, value pairs correspond to SymPy name and numexpr name
    # functions not appearing in this dict will raise a TypeError
    printmethod='_numexprcode'
    _numexpr_functions={
        'sin':'sin','cos':'cos','tan':'tan','asin':'arcsin','acos':'arccos','atan':'arctan',
        'atan2':'arctan2','sinh':'sinh','cosh':'cosh','tanh':'tanh','asinh':'arcsinh',
        'acosh':'arccosh','atanh':'arctanh','ln':'log','log':'log','exp':'exp','sqrt':'sqrt',
        'Abs':'abs','conjugate':'conj','im':'imag','re':'real','where':'where',
        'complex':'complex','contains':'contains',
    }
    module='numexpr'
    def _print_ImaginaryUnit(self,expr):return '1j'
    def _print_seq(self,seq,delimiter=', '):
        s=[self._print(item) for item in seq]
        return delimiter.join(s) if s else ""
    def _print_Function(self,e):
        func_name=getattr(e,'class_name',None) or e.func.__name__
        nstr=self._numexpr_functions.get(func_name,None)
        if nstr is None:
            # check for implemented_function
            if hasattr(e,'_imp_'):
                return "(%s)"%self._print(e._imp_(*e.args))
            if func_name in _NUMEXPR_PYTHON_KNOWN:
                # Printed by the PythonCodePrinter known functions.
                return PythonCodePrinter.emptyPrinter(self,e)
            raise TypeError("numexpr does not support function '%s'"%func_name)
        return "%s(%s)"%(nstr,self._print_seq(e.args))
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Matrix',False):return self.blacklisted(expr)
        if getattr(expr,'is_Function',False):return self._print_Function(expr)
        return super().emptyPrinter(expr)
    def _print_Piecewise(self,expr):
        "Piecewise function printer"
        exprs=[self._print(arg.expr) for arg in expr.args]
        conds=[self._print(arg.cond) for arg in expr.args]
        ans=[];parenthesis_count=0;is_last_cond_True=False
        for cond,value in zip(conds,exprs):
            if cond=='True':
                ans.append(value);is_last_cond_True=True
                break
            ans.append('where(%s, %s, '%(cond,value));parenthesis_count+=1
        if not is_last_cond_True:
            # See https://github.com/pydata/numexpr/issues/298
            ans.append('log(-1)')
        return ''.join(ans)+')'*parenthesis_count
    def _print_ITE(self,expr):
        import symox as _S
        return self._print(expr.rewrite(_S.Piecewise))
    def blacklisted(self,expr):
        from ._classes import class_names
        name=class_names(expr)[0] if hasattr(expr,'args') else expr.__class__.__name__
        raise TypeError("numexpr cannot be used with %s"%name)
    # blacklist some Python expressions
    _print_list=_print_tuple=_print_dict=_print_Dict=blacklisted
    def _print_Tuple(self,expr):return self.blacklisted(expr)
    def _print_NumExprEvaluate(self,expr):
        evaluate=self.module+".evaluate"
        return "%s('%s', truediv=True)"%(evaluate,self._print(expr.expr))
    def doprint(self,expr,assign_to=None):
        from ..codegen.ast import CodegenAST
        from ..codegen.pynodes import NumExprEvaluate
        if not isinstance(expr,CodegenAST):
            expr=NumExprEvaluate(expr)
        return super().doprint(expr,assign_to)
    def _print_Return(self,expr):
        from ..codegen.pynodes import NumExprEvaluate
        r,=expr.args
        if not isinstance(r,NumExprEvaluate):r=NumExprEvaluate(r)
        return 'return %s'%self._print(r)
    def _print_Assignment(self,expr):
        from ..codegen.pynodes import NumExprEvaluate
        lhs,rhs=expr.lhs,expr.rhs
        if not isinstance(rhs,NumExprEvaluate):rhs=NumExprEvaluate(rhs)
        return '%s = %s'%(self._print(lhs),self._print(rhs))
    def _print_CodeBlock(self,expr):
        # SymPy 1.14 passes the printer twice to the parent method, which
        # raises TypeError.
        raise TypeError("CodePrinter._print_CodeBlock() takes 2 positional arguments but 3 were given")
    def _print_FunctionDefinition(self,fd):
        # SymPy's PythonCodePrinter prints the body statements one by one.
        from ._nodes import indent
        body='\n'.join(self._print(arg) for arg in fd.body.args)
        return "def {name}({parameters}):\n{body}".format(
            name=self._print(fd.name),parameters=', '.join(self._print(var.symbol) for var in fd.parameters),
            body=indent(body))


class IntervalPrinter(MpmathPrinter):_module='mpmath.iv'
def lambdarepr(expr,**settings):return LambdaPrinter(settings).doprint(expr)
__all__=['LambdaPrinter','NumPyPrinter','NumExprPrinter','MpmathPrinter','PythonCodePrinter','IntervalPrinter','lambdarepr']
