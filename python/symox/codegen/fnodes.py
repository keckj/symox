"""Fortran declarations, procedures, array expressions and intrinsics."""
from .ast import _Expression,Token,Node,Variable,Attribute,CodegenAST,CodeBlock,String,none,NoneToken,real,untyped
import symox as S

intent_in=Attribute('intent_in');intent_out=Attribute('intent_out');intent_inout=Attribute('intent_inout')
allocatable=Attribute('allocatable')
pure=Attribute('pure');elemental=Attribute('elemental')

class ArrayConstructor(Token):_fields=('elements',)
class Do(Token):
    _fields=('body','counter','first','last','step','concurrent');defaults={'step':S.Integer(1),'concurrent':False}
class FortranReturn(Token):_fields=('return_value',);defaults={'return_value':none}
class GoTo(Token):_fields=('labels','expr');defaults={'expr':none}
class ImpliedDoLoop(Token):_fields=('expr','counter','first','last','step');defaults={'step':S.Integer(1)}
def _name(arg):
    return str(arg.name) if hasattr(arg,'name') else str(arg)
class Module(Token):
    _fields=('name','declarations','definitions')
    def _coerce(self,field,value):
        if field=='declarations':return CodeBlock(*[String(a) if isinstance(a,str) else a for a in value])
        if field=='definitions':return value if isinstance(value,CodeBlock) else CodeBlock(*value)
        return super()._coerce(field,value)
class Program(Token):_fields=('name','body')
class Subroutine(Node):
    _fields=('name','parameters','body','attrs');defaults={'attrs':()}
    def _coerce(self,field,value):
        if field=='parameters':return tuple(Variable.deduced(p) for p in value)
        if field=='body' and not isinstance(value,CodeBlock):return CodeBlock(*value)
        return super()._coerce(field,value)
class SubroutineCall(Token):
    _fields=('name','subroutine_args')
    def _coerce(self,field,value):
        if field=='name':return _name(value)
        if field=='subroutine_args':return tuple(a if isinstance(a,CodegenAST) else S.sympify(a) for a in value)
        return super()._coerce(field,value)
class use_rename(Token):
    _fields=('local','original')
    def _coerce(self,field,value):return String(value) if isinstance(value,str) else value
class use(Token):
    _fields=('namespace','rename','only');defaults={'rename':none,'only':none}
    def _coerce(self,field,value):
        if field=='namespace':return _name(value)
        if field=='rename' and not isinstance(value,NoneToken):return tuple(a if isinstance(a,use_rename) else use_rename(*a) for a in value)
        if field=='only' and not isinstance(value,NoneToken):return tuple(a if isinstance(a,use_rename) else _name(a) for a in value)
        return super()._coerce(field,value)
class sum_(Token,_Expression):
    _fields=('array','dim','mask');defaults={'dim':none,'mask':none}
    def _symox_(self):
        if self.mask is not none:raise NotImplementedError("arithmetic on a masked %s is not supported"%type(self).__name__)
        args=(self.array,) if self.dim is none else (self.array,self.dim)
        return S.Function(type(self).__name__.rstrip('_'))(*args)
class product_(sum_):pass

class FFunction(S.Function):
    """A Fortran intrinsic, printed by ``fcode`` once the standard allows it."""
    _required_standard=77
    _symox_user_hooks=True
    def _fcode(self,printer):
        name=self.__class__.__name__
        if printer._settings['standard']<self._required_standard:
            raise NotImplementedError("%s requires Fortran %d or newer"%(name,self._required_standard))
        return '{}({})'.format(name,', '.join(map(printer._print,self.args)))
class F95Function(FFunction):
    _required_standard=95
class Extent:
    """A Fortran array extent ``low:high``, or ``:`` for an assumed shape."""
    def __init__(self,*args):
        if len(args)==2:self.args=tuple(S.sympify(a) for a in args)
        elif len(args)==0 or (len(args)==1 and args[0] in (':',None)):self.args=()
        else:raise ValueError("Expected 0 or 2 args (or one argument == None or ':')")
    def __eq__(self,other):return type(other) is Extent and other.args==self.args
    def __hash__(self):return hash((Extent,self.args))
    def __str__(self):return ':' if not self.args else ':'.join(str(a) for a in self.args)
    def __repr__(self):return 'Extent(%s)'%', '.join(str(a) for a in self.args)
assumed_extent=Extent()
class kind(FFunction):
    """Fortran kind function."""
    nargs=1
class cmplx(FFunction):
    """Fortran complex conversion function."""
    nargs=2
class dsign(FFunction):
    """Fortran sign intrinsic for double precision arguments."""
    nargs=2
class isign(FFunction):
    """Fortran sign intrinsic for integer arguments."""
    nargs=2
class merge(F95Function):
    """Fortran merge function."""
    nargs=3
class _fortran_literal(S.Function):
    """A double precision literal printed as SymPy's ``literal_dp`` (``1d0``).

    Symox's Float cannot be subclassed, so the Fortran printer uses this
    node where SymPy builds ``literal_dp`` values internally.
    """
    nargs=1
    _symox_user_hooks=True
    def _fcode(self,printer,*args,**kwargs):
        mantissa,sgnd_ex=('%.17e'%float(self.args[0])).split('e')
        mantissa=mantissa.strip('0').rstrip('.')
        ex_sgn,ex_num=sgnd_ex[0],sgnd_ex[1:].lstrip('0')
        ex_sgn='' if ex_sgn=='+' else ex_sgn
        return (mantissa or '0')+'d'+ex_sgn+(ex_num or '0')

def dimension(*args):
    if len(args)>7:raise ValueError("Fortran only supports up to 7 dimensional arrays")
    parameters=[]
    for arg in args:
        if isinstance(arg,Extent):parameters.append(arg)
        elif isinstance(arg,str):parameters.append(Extent() if arg==':' else String(arg))
        elif isinstance(arg,(tuple,list)):parameters.append(Extent(*arg))
        else:parameters.append(S.sympify(arg))
    if len(args)==0:raise ValueError("Need at least one dimension")
    return Attribute('dimension',tuple(parameters))
assumed_size=dimension('*')
def bind_C(name=None):return Attribute('bind_C',() if name is None else (String(name),))
def array(symbol,dim,intent=None,*,attrs=(),value=None,type=None):
    if isinstance(dim,Attribute):
        if str(dim.name)!='dimension':raise ValueError("Got an unexpected Attribute argument as dim: %s"%str(dim))
    else:dim=dimension(*dim)
    attrs=list(attrs)+[dim]
    if intent is not None:
        if intent not in (intent_in,intent_out,intent_inout):
            intent={'in':intent_in,'out':intent_out,'inout':intent_inout}[intent]
        attrs.append(intent)
    if type is None:return Variable.deduced(symbol,value=value,attrs=attrs)
    return Variable(symbol,type,none if value is None else value,attrs)

def _intrinsic(name,*args):
    from .ast import FunctionCall
    return FunctionCall(name,tuple(a for a in args if a is not None))
def allocated(array):return _intrinsic('allocated',array)
def lbound(array,dim=None,kind=None):return _intrinsic('lbound',array,dim,kind)
def ubound(array,dim=None,kind=None):return _intrinsic('ubound',array,dim,kind)
def shape(source,kind=None):return _intrinsic('shape',source,kind)
def size(array,dim=None,kind=None):return _intrinsic('size',array,dim,kind)
def reshape(source,shape,pad=None,order=None):return _intrinsic('reshape',source,shape,pad,order)
def literal_dp(num,dps=None,precision=None):return S.Float(num,dps=dps,precision=precision)
def literal_sp(num,dps=None,precision=None):return S.Float(num,dps=dps,precision=precision)

__all__=['pure','elemental','intent_in','intent_out','intent_inout','allocatable','assumed_extent','assumed_size','ArrayConstructor','Do','FortranReturn','GoTo','ImpliedDoLoop','Module','Program','Subroutine','SubroutineCall',
 'use','use_rename','sum_','product_','FFunction','F95Function','Extent','kind','cmplx','dsign','isign','merge','dimension',
 'bind_C','array','allocated','lbound','ubound','shape','size','reshape','literal_dp','literal_sp']
