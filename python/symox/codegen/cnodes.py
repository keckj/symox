"""C language statements, aggregates and side-effecting expressions."""
from .ast import (CodegenAST,Token,none,_Expression,restrict,volatile,Attribute,Declaration,Node,
 String,Type,FunctionCall,CodeBlock)
import symox as S

void=Type('void')
static=Attribute('static')

class _Unary(CodegenAST,_Expression):
    _fields=('arg',);operator='';prefix=True
    def _symox_(self):return S.Symbol(self.operator+str(self.arg) if self.prefix else str(self.arg)+self.operator)
class PreIncrement(_Unary):operator='++'
class PreDecrement(_Unary):operator='--'
class PostIncrement(_Unary):operator='++';prefix=False
class PostDecrement(_Unary):operator='--';prefix=False
class CommaOperator(CodegenAST,_Expression):
    def __init__(self,*args):object.__setattr__(self,'_values',tuple(args))
    args=property(lambda self:self._values)
    def _symox_(self):return S._core._apply('CommaOperator',*(S.sympify(v) for v in self.args),evaluate=False)
class Label(Token):
    _fields=('name','body');defaults={'body':none}
    # SymPy's ``_construct_name = String`` and ``_construct_body``.
    def _coerce(self,field,value):
        if field=='name':return value if isinstance(value,String) else String(value)
        if field=='body' and value != none and not isinstance(value,CodeBlock):return CodeBlock(*value)
        return super()._coerce(field,value)
class goto(Token):
    _fields=('label',)
    # SymPy's ``_construct_label = Label``.
    def _coerce(self,field,value):
        if field=='label':return value if isinstance(value,Label) else Label(value)
        return super()._coerce(field,value)
def alignof(arg):
    """Generate a FunctionCall instance for calling 'alignof'."""
    from .ast import FunctionCall,String
    return FunctionCall('alignof',[String(arg) if isinstance(arg,str) else arg])
def sizeof(arg):
    """Generate a FunctionCall instance for calling 'sizeof'."""
    from .ast import FunctionCall,String
    return FunctionCall('sizeof',[String(arg) if isinstance(arg,str) else arg])
class struct(Token):_fields=('name','declarations');defaults={'name':none}
class union(struct):pass
__all__=['Attribute','Declaration','Node','String','Token','Type','none','FunctionCall','CodeBlock','void','static','restrict','volatile','PreIncrement','PreDecrement','PostIncrement','PostDecrement','CommaOperator','Label','goto','sizeof','alignof','struct','union']
