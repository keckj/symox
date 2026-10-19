"""Extensible Python, NumPy and arbitrary-precision expression printers."""
import ast
import keyword
from .cxx import CXX17CodePrinter
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
    def _print_Dummy(self,expr):return self._settings.get('symbol_names',{}).get(expr,expr.name+'_'+str(expr.dummy_index))
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
            if custom is None:custom=name if getattr(expr,'is_UndefinedFunction',False) else self._function(name)
            return custom+'('+', '.join(self._print(a) for a in expr.args)+')'
        return super().emptyPrinter(expr)


class PythonCodePrinter(AbstractPythonCodePrinter):pass
class CmathPrinter(PythonCodePrinter):_module='cmath'
class MpmathPrinter(PythonCodePrinter):
    _module='mpmath'
    def _print_Float(self,expr):return "mpmath.mpf('"+str(expr)+"')"
    def _print_Rational(self,expr):return '(mpmath.mpf('+str(expr.p)+')/mpmath.mpf('+str(expr.q)+'))'
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
    def _boolean(self,expr,op):
        value=self._print(expr.args[0])
        for arg in expr.args[1:]:value='numpy.'+op+'('+value+', '+self._print(arg)+')'
        return value
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Matrix',False) and hasattr(expr,'tolist'):return 'numpy.array('+super().emptyPrinter(expr)+')'
        return super().emptyPrinter(expr)


def pycode(expr,assign_to=None,**settings):return PythonCodePrinter(settings).doprint(expr,assign_to)
__all__=['pycode','AbstractPythonCodePrinter','PythonCodePrinter','CmathPrinter','MpmathPrinter','SymPyPrinter','ArrayPrinter']
