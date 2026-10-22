"""C language statements, aggregates and side-effecting expressions."""
from .ast import CodegenAST,Token,none,_Expression,restrict,volatile
import symox as S

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
class Label(Token):_fields=('name','body');defaults={'body':none}
class goto(Token):_fields=('label',)
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
__all__=['restrict','volatile','PreIncrement','PreDecrement','PostIncrement','PostDecrement','CommaOperator','Label','goto','sizeof','alignof','struct','union']
