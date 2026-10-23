"""Extensible Python, NumPy and arbitrary-precision expression printers."""
import ast
import keyword
from ._legacy import LegacyCXX17CodePrinter as CXX17CodePrinter
from .precedence import PRECEDENCE
from ._nodes import NodePrinter


def _indexed_symbol_name(name):
    """Recognize the array references used by downstream symbolic printers."""
    try:
        expression=ast.parse(name, mode='eval').body
    except (SyntaxError, ValueError):
        return False
    if not isinstance(expression, ast.Subscript):
        return False
    allowed=(ast.Subscript,ast.Name,ast.Load,ast.Tuple,ast.Constant,ast.Slice,
             ast.UnaryOp,ast.USub,ast.UAdd)
    return all(isinstance(node,allowed) for node in ast.walk(expression))


class AbstractPythonCodePrinter(NodePrinter,CXX17CodePrinter):
    _style='python'
    language='Python';printmethod='_pythoncode';_module='math'
    def _print_Symbol(self,expr):
        mapping=self._settings.get('symbol_names',{})
        if expr in mapping:return mapping[expr]
        name=expr.name
        if (not name.isidentifier() or keyword.iskeyword(name)) and not _indexed_symbol_name(name):
            raise ValueError('invalid Python identifier: '+name)
        return name
    def _print_Dummy(self,expr):
        mapping=self._settings.get('symbol_names',{})
        if expr in mapping:return mapping[expr]
        # SymPy's CodePrinter._print_Dummy.
        if expr.name.startswith('Dummy_'):return '_'+expr.name
        return expr.name+'_'+str(expr.dummy_index)
    def _print_Rational(self,expr):return '('+str(expr.p)+'/'+str(expr.q)+')'
    def _print_Float(self,expr):return str(expr)
    def _print_Pow(self,expr):
        return '('+self.parenthesize(expr.base,PRECEDENCE['Pow'])+'**'+self.parenthesize(expr.exp,PRECEDENCE['Pow'])+')'
    def _print_BooleanTrue(self,expr):return 'True'
    def _print_BooleanFalse(self,expr):return 'False'
    def _print_ImaginaryUnit(self,expr):return '1j'
    def _print_Piecewise(self,expr):
        from .mathml import _pairs
        pairs=_pairs(expr);tail="float('nan')"
        for value,condition in reversed(pairs):
            tail=self._print(value) if str(condition) in ('True','true') else '('+self._print(value)+' if '+self._print(condition)+' else '+tail+')'
        return tail
    def _print_Assignment(self,expr):return self._print(expr.lhs)+' = '+self._print(expr.rhs)
    def _print_CodeBlock(self,expr):return '\n'.join(self._print(v) for v in expr.args)
    def _print_Comment(self,expr):return '# '+expr.text
    def _print_And(self,expr):return '('+' and '.join(self._print(a) for a in expr.args)+')'
    def _print_Or(self,expr):return '('+' or '.join(self._print(a) for a in expr.args)+')'
    def _print_Not(self,expr):return '(not '+self._print(expr.args[0])+')'
    def _relation(self,expr,op):return '('+self._print(expr.lhs)+' '+op+' '+self._print(expr.rhs)+')'
    def _print_Equality(self,expr):return self._relation(expr,'==')
    def _print_Unequality(self,expr):return self._relation(expr,'!=')
    def _print_LessThan(self,expr):return self._relation(expr,'<=')
    def _print_StrictLessThan(self,expr):return self._relation(expr,'<')
    def _print_GreaterThan(self,expr):return self._relation(expr,'>=')
    def _print_StrictGreaterThan(self,expr):return self._relation(expr,'>')
    def _print_Tuple(self,expr):return '('+', '.join(self._print(v) for v in expr.args)+(',' if len(expr.args)==1 else '')+')'
    def _print_list(self,expr):return '['+', '.join(self._print(v) for v in expr)+']'
    def _print_tuple(self,expr):return '('+', '.join(self._print(v) for v in expr)+(',' if len(expr)==1 else '')+')'
    def _function(self,name):
        aliases={'Abs':'abs','ceiling':'ceil','Min':'min','Max':'max','ln':'log'}
        name=aliases.get(name,name)
        return name if name in ('abs','min','max') or not self._module else self._module+'.'+name
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Matrix',False) and hasattr(expr,'tolist'):
            return '['+', '.join('['+', '.join(self._print(v) for v in row)+']' for row in expr.tolist())+']'
        name=getattr(expr,'class_name',type(expr).__name__)
        constants={'Pi':self._module+'.pi','Exp1':self._module+'.e','Infinity':"float('inf')",'NegativeInfinity':"-float('inf')",'NaN':"float('nan')",'ComplexInfinity':"complex(float('inf'),float('inf'))"}
        if name in constants:return constants[name]
        if getattr(expr,'is_Function',False):
            custom=self._settings.get('user_functions',{}).get(name)
            if custom is None and name in _CFUNCTION_NAMES and type(expr).__module__=='symox.codegen.cfunctions':
                custom=self._cfunction(expr,name)
            if custom is None:custom=name if getattr(expr,'is_UndefinedFunction',False) else self._function(name)
            return custom+'('+', '.join(self._print(a) for a in expr.args)+')'
        return super().emptyPrinter(expr)


    def _cfunction(self,expr,name):
        # The ``sympy.codegen.cfunctions`` names the reference Python printers
        # know: math has no cbrt, exp2 or fma, numpy has exp2 but no cbrt or fma.
        known={'Sqrt':'sqrt','log2':'log2','log10':'log10','expm1':'expm1','log1p':'log1p',
               'hypot':'hypot','isnan':'isnan','isinf':'isinf'}
        if self._module=='numpy':known['exp2']='exp2'
        if self._module in ('math','numpy') and name in known:return self._module+'.'+known[name]
        if self._module in ('math','numpy') and self._settings.get('strict',True) is not False:
            from .codeprinter import PrintMethodNotImplementedError
            raise PrintMethodNotImplementedError(
                f"Unsupported by {type(self)}: {name}" +
                "\nSet the printer option 'strict' to False in order to generate partially printed code.")
        return None


_CFUNCTION_NAMES=('Sqrt','Cbrt','exp2','log2','log10','expm1','log1p','fma','hypot','isnan','isinf')


class PythonCodePrinter(AbstractPythonCodePrinter):pass
class CmathPrinter(PythonCodePrinter):_module='cmath'
class MpmathPrinter(PythonCodePrinter):
    """Python with mpmath, keeping the precision of Float literals."""
    _module='mpmath'
    def _print_Float(self,expr):
        # SymPy prints the exact binary mantissa of the Float.
        from .._host import _float_mpf
        return 'mpmath.mpf('+str(tuple(map(int,_float_mpf(expr))))+')'
    def _print_Rational(self,expr):return 'mpmath.mpf('+str(expr.p)+')/mpmath.mpf('+str(expr.q)+')'
    def _print_Half(self,expr):return self._print_Rational(expr)
    def _print_Pow(self,expr):
        import symox as _S
        base,exp=expr.args
        if exp==_S.Rational(1,2):return 'mpmath.sqrt('+self._print(base)+')'
        if exp==_S.Rational(-1,2) and expr.is_commutative:return '1/mpmath.sqrt('+self._print(base)+')'
        return self.parenthesize(base,PRECEDENCE['Pow'])+'**'+self.parenthesize(exp,PRECEDENCE['Pow'])
    def emptyPrinter(self,expr):
        name=getattr(expr,'class_name',type(expr).__name__)
        constants={'GoldenRatio':'phi','EulerGamma':'euler','Catalan':'catalan','NaN':'nan','Infinity':'inf','NegativeInfinity':'ninf'}
        if name in constants:return 'mpmath.'+constants[name]
        return super().emptyPrinter(expr)
    def _function(self,name):
        # mpmath spellings of the functions SymPy's lambdify translates.
        from ..utilities.lambdify import _MPMATH_TRANSLATIONS
        if name!='Abs' and name in _MPMATH_TRANSLATIONS:return 'mpmath.'+_MPMATH_TRANSLATIONS[name]
        return super()._function(name)
    _print_Infinity=lambda self,expr:'mpmath.inf'
    _print_NegativeInfinity=lambda self,expr:'mpmath.ninf'
    _print_NaN=lambda self,expr:'mpmath.nan'
    _print_EulerGamma=lambda self,expr:'mpmath.euler'
    _print_Catalan=lambda self,expr:'mpmath.catalan'
    _print_GoldenRatio=lambda self,expr:'mpmath.phi'
class SymPyPrinter(PythonCodePrinter):_module='symox'
class ArrayPrinter(PythonCodePrinter):
    _module='numpy'
    def _function(self,name):
        return 'numpy.'+{'Abs':'abs','ceiling':'ceil','Min':'minimum','Max':'maximum','asin':'arcsin','acos':'arccos','atan':'arctan','asinh':'arcsinh','acosh':'arccosh','atanh':'arctanh'}.get(name,name)
    def _print_Piecewise(self,expr):
        from .mathml import _pairs
        pairs=_pairs(expr)
        return 'numpy.select(['+', '.join(self._print(c) for _,c in pairs)+'], ['+', '.join(self._print(v) for v,_ in pairs)+"], default=numpy.nan)"
    def _print_And(self,expr):return self._boolean(expr,'logical_and')
    def _print_Or(self,expr):return self._boolean(expr,'logical_or')
    def _print_Not(self,expr):return 'numpy.logical_not('+self._print(expr.args[0])+')'
    def _print_MatrixSolve(self,expr):return '%s(%s, %s)'%(self._module+'.linalg.solve',self._print(expr.matrix),self._print(expr.vector))
    def _boolean(self,expr,op):
        value=self._print(expr.args[0])
        for arg in expr.args[1:]:value='numpy.'+op+'('+value+', '+self._print(arg)+')'
        return value
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Matrix',False) and hasattr(expr,'tolist'):return 'numpy.array('+super().emptyPrinter(expr)+')'
        return super().emptyPrinter(expr)


def pycode(expr,assign_to=None,**settings):return PythonCodePrinter(settings).doprint(expr,assign_to)
__all__=['pycode','AbstractPythonCodePrinter','PythonCodePrinter','CmathPrinter','MpmathPrinter','SymPyPrinter','ArrayPrinter']
