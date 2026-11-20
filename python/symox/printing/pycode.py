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
    def _print_Tuple(self,expr):
        from ..codegen.abstract_nodes import List
        # SymPy's AbstractPythonCodePrinter._print_List (codegen List nodes).
        if isinstance(expr,List):return '['+', '.join(self._print(v) for v in expr.args)+']'
        return '('+', '.join(self._print(v) for v in expr.args)+(',' if len(expr.args)==1 else '')+')'
    def doprint(self,expr,assign_to=None):
        # SymPy's CodePrinter.doprint: the nodes printed by
        # ``_print_not_supported`` are listed in front of the code.
        self._not_supported=set()
        result=super().doprint(expr,assign_to)
        not_supported,self._not_supported=self._not_supported,set()
        if not_supported and self._settings.get('human',True) and isinstance(result,str):
            front=['  # Not supported in {}:'.format(self.language)]
            front+=['  # {}'.format(type(e).__name__) for e in sorted(not_supported,key=str)]
            result='\n'.join(front+result.splitlines())
        return result
    def _print_Derivative(self,expr):
        # SymPy's CodePrinter._print_Derivative.
        obj,*wrt_order_pairs=expr.args
        for func_arg in obj.args:
            if not func_arg.is_Symbol:
                raise ValueError("%s._print_Derivative(...) only supports functions with symbols as arguments."%
                                 self.__class__.__name__)
        meth_name='_print_Derivative_%s'%obj.func.__name__
        pmeth=getattr(self,meth_name,None)
        if pmeth is None:
            strict=self._settings.get('strict',True)
            if strict is None:strict=self._settings.get('human',True)==True
            if strict:
                from .codeprinter import PrintMethodNotImplementedError
                raise PrintMethodNotImplementedError(
                    f"Unsupported by {type(self)}: {type(expr)}"+
                    f"\nPrinter has no method: {meth_name}"+
                    "\nSet the printer option 'strict' to False in order to generate partially printed code.")
            return self._print_not_supported(expr)
        orders=dict(wrt_order_pairs)
        seq_orders=[orders[arg] for arg in obj.args]
        return pmeth(obj.args,seq_orders)
    def _print_not_supported(self,expr):
        # SymPy's CodePrinter._print_not_supported (strict defaults to
        # ``human``); the non-strict text is the repr of the node.
        strict=self._settings.get('strict',True)
        if strict is None:strict=self._settings.get('human',True)==True
        if strict:
            from .codeprinter import PrintMethodNotImplementedError
            kind=expr.func if getattr(expr,'is_Function',False) else type(expr)
            raise PrintMethodNotImplementedError(
                f"Unsupported by {type(self)}: {kind}" +
                "\nSet the printer option 'strict' to False in order to generate partially printed code.")
        try:self.__dict__.setdefault('_not_supported',set()).add(expr)
        except TypeError:pass
        return repr(expr)
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Matrix',False) and hasattr(expr,'tolist'):
            return '['+', '.join('['+', '.join(self._print(v) for v in row)+']' for row in expr.tolist())+']'
        name=getattr(expr,'class_name',type(expr).__name__)
        if type(expr).__module__ in _UNSUPPORTED_NODE_MODULES and not (
                getattr(expr,'is_Function',False) and name in self._settings.get('user_functions',{})):
            # The numpy/scipy/python codegen nodes the printer has no method for.
            return self._print_not_supported(expr)
        constants={'Pi':self._module+'.pi','Exp1':self._module+'.e','Infinity':"float('inf')",'NegativeInfinity':"-float('inf')",'NaN':"float('nan')",'ComplexInfinity':"complex(float('inf'),float('inf'))"}
        if name in constants:return constants[name]
        if getattr(expr,'is_Function',False):
            custom=self._settings.get('user_functions',{}).get(name)
            if custom is None and name in _CFUNCTION_NAMES and type(expr).__module__=='symox.codegen.cfunctions':
                custom=self._cfunction(expr,name)
            if custom is None:custom=name if getattr(expr,'is_UndefinedFunction',False) else self._function(name)
            return custom+'('+', '.join(self._print(a) for a in expr.args)+')'
        try:
            return super().emptyPrinter(expr)
        except NotImplementedError as exc:
            # The native printer refuses this very node (e.g. Derivative):
            # SymPy's printers go through _print_not_supported.
            if type(exc) is NotImplementedError and str(exc)=='Unsupported by the code printer: '+name:
                return self._print_not_supported(expr)
            raise


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


_UNSUPPORTED_NODE_MODULES=('symox.codegen.numpy_nodes','symox.codegen.scipy_nodes','symox.codegen.pynodes')
_CFUNCTION_NAMES=('Sqrt','Cbrt','exp2','log2','log10','expm1','log1p','fma','hypot','isnan','isinf')


class PythonCodePrinter(AbstractPythonCodePrinter):pass
class CmathPrinter(PythonCodePrinter):_module='cmath'
class MpmathPrinter(PythonCodePrinter):
    """Python with mpmath, keeping the precision of Float literals."""
    _module='mpmath'
    language='Python with mpmath'
    def _print_Float(self,expr):
        # SymPy prints the exact binary mantissa of the Float.
        from .._host import _float_mpf
        return 'mpmath.mpf('+str(tuple(map(int,_float_mpf(expr))))+')'
    def _print_Rational(self,expr):return 'mpmath.mpf('+str(expr.p)+')/mpmath.mpf('+str(expr.q)+')'
    def _print_Half(self,expr):return self._print_Rational(expr)
    def _print_Derivative_zeta(self,args,seq_orders):
        arg,=args
        deriv_order,=seq_orders
        return '{}({}, derivative={})'.format('mpmath.zeta',self._print(arg),deriv_order)
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
    language='Python with numpy'
    def _function(self,name):
        return 'numpy.'+{'Abs':'abs','ceiling':'ceil','Min':'minimum','Max':'maximum','asin':'arcsin','acos':'arccos','atan':'arctan','asinh':'arcsinh','acosh':'arccosh','atanh':'arctanh'}.get(name,name)
    def _print_Piecewise(self,expr):
        from .mathml import _pairs
        pairs=_pairs(expr)
        return 'numpy.select(['+', '.join(self._print(c) for _,c in pairs)+'], ['+', '.join(self._print(v) for v,_ in pairs)+"], default=numpy.nan)"
    def _print_And(self,expr):return self._boolean(expr,'logical_and')
    def _print_Or(self,expr):return self._boolean(expr,'logical_or')
    def _print_Not(self,expr):return 'numpy.logical_not('+self._print(expr.args[0])+')'
    def _print_sinc(self,expr):
        # NumPyPrinter._print_sinc: numpy.sinc is the normalized sinc.
        import symox as _S
        return self._module+'.sinc('+self._print(expr.args[0]/_S.pi)+')'
    def _print_MatrixSolve(self,expr):return '%s(%s, %s)'%(self._module+'.linalg.solve',self._print(expr.matrix),self._print(expr.vector))
    # The ``sympy.codegen.numpy_nodes`` functions (NumPyPrinter).
    def _print_logaddexp(self,expr):return self._module+'.logaddexp('+', '.join(self._print(a) for a in expr.args)+')'
    def _print_logaddexp2(self,expr):return self._module+'.logaddexp2('+', '.join(self._print(a) for a in expr.args)+')'
    def _helper_minimum_maximum(self,op,*args):
        if len(args)==0:raise NotImplementedError(f"Need at least one argument for {op}")
        if len(args)==1:return self._print(args[0])
        return 'functools.reduce(%s, [%s])'%(op,', '.join(self._print(a) for a in args))
    # SymPy's NumPyPrinter._print_Min/_print_Max: a binary ufunc folded over
    # the sorted arguments (a third positional argument would be ``out``).
    def _print_Min(self,expr):return self._helper_minimum_maximum(self._module+'.minimum',*expr.args)
    def _print_Max(self,expr):return self._helper_minimum_maximum(self._module+'.maximum',*expr.args)
    def _print_minimum(self,expr):return self._helper_minimum_maximum(self._module+'.minimum',*expr.args)
    def _print_maximum(self,expr):return self._helper_minimum_maximum(self._module+'.maximum',*expr.args)
    def _print_amin(self,expr):return '{}({}, axis={})'.format(self._module+'.amin',self._print(expr.array),self._print(expr.axis))
    def _print_amax(self,expr):return '{}({}, axis={})'.format(self._module+'.amax',self._print(expr.array),self._print(expr.axis))
    def _boolean(self,expr,op):
        value=self._print(expr.args[0])
        for arg in expr.args[1:]:value='numpy.'+op+'('+value+', '+self._print(arg)+')'
        return value
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Matrix',False) and hasattr(expr,'tolist'):return 'numpy.array('+super().emptyPrinter(expr)+')'
        return super().emptyPrinter(expr)


def pycode(expr,assign_to=None,**settings):return PythonCodePrinter(settings).doprint(expr,assign_to)
__all__=['pycode','AbstractPythonCodePrinter','PythonCodePrinter','CmathPrinter','MpmathPrinter','SymPyPrinter','ArrayPrinter']
