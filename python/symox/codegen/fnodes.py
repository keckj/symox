"""Fortran declarations, procedures, array expressions and intrinsics."""
from .ast import Token,Node,Variable,Attribute,CodegenAST,none,real,untyped
import symox as S

class ArrayConstructor(Token):_fields=('elements',)
class Do(Token):
    _fields=('body','counter','first','last','step','concurrent');defaults={'step':S.Integer(1),'concurrent':False}
class FortranReturn(Token):_fields=('return_value',);defaults={'return_value':none}
class GoTo(Token):_fields=('labels','expr');defaults={'expr':none}
class ImpliedDoLoop(Token):_fields=('expr','counter','first','last','step');defaults={'step':S.Integer(1)}
class Module(Token):_fields=('name','declarations','definitions');defaults={'declarations':()}
class Program(Token):_fields=('name','body')
class Subroutine(Node):_fields=('name','parameters','body','attrs');defaults={'attrs':()}
class SubroutineCall(Token):_fields=('name','subroutine_args')
class use(Token):_fields=('namespace','rename','only');defaults={'rename':none,'only':none}
class use_rename(Token):_fields=('local','original')
class sum_(Token):_fields=('array','dim','mask');defaults={'dim':none,'mask':none}
class product_(sum_):pass

class FFunction:
    def __new__(cls,*args,**kwargs):return S._core._apply(cls.__name__,*(S.sympify(a) for a in args),evaluate=False)
class F95Function(FFunction):pass
class Extent(FFunction):pass
class kind(FFunction):pass
class cmplx(FFunction):pass
class dsign(FFunction):pass
class isign(FFunction):pass
class merge(FFunction):pass

def dimension(*dims):return Attribute('dimension',dims)
def bind_C(name=None):return Attribute('bind_C',() if name is None else (str(name),))
def array(symbol,dim,intent=None,*,attrs=(),value=None,type=None):
    if not isinstance(dim,(tuple,list)):dim=(dim,)
    extra=[dimension(*dim)]
    if intent is not None:extra.append(Attribute('intent',(intent,)))
    return Variable(symbol,type or real,none if value is None else value,tuple(attrs)+tuple(extra))

def _intrinsic(name,*args):
    from .ast import FunctionCall
    return FunctionCall(name,tuple(a for a in args if a is not None))
def allocated(array):return _intrinsic('allocated',array)
def lbound(array,dim=None,kind=None):return _intrinsic('lbound',array,dim,kind)
def ubound(array,dim=None,kind=None):return _intrinsic('ubound',array,dim,kind)
def shape(source,kind=None):return _intrinsic('shape',source,kind)
def size(array,dim=None,kind=None):return _intrinsic('size',array,dim,kind)
def reshape(source,shape,pad=None,order=None):return _intrinsic('reshape',source,shape,pad,order)
def literal_dp(num,dps=None,precision=None):return S.Float(str(num),dps=dps or 17) if precision is None else S.Float(str(num),precision=precision)
def literal_sp(num,dps=None,precision=None):return S.Float(str(num),dps=dps or 8) if precision is None else S.Float(str(num),precision=precision)

__all__=['ArrayConstructor','Do','FortranReturn','GoTo','ImpliedDoLoop','Module','Program','Subroutine','SubroutineCall',
 'use','use_rename','sum_','product_','FFunction','F95Function','Extent','kind','cmplx','dsign','isign','merge','dimension',
 'bind_C','array','allocated','lbound','ubound','shape','size','reshape','literal_dp','literal_sp']
