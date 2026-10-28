"""Protocol evidence for the exact LLL kernel; no matrix arithmetic here."""
import abc as _abc
import _abc as _abc_native
import builtins as _builtins
import sys as _sys
from types import FunctionType,ModuleType,MethodDescriptorType,WrapperDescriptorType,GetSetDescriptorType,BuiltinFunctionType,CodeType
from symox import _core,_basic
from symox.polys.domains import _methods

# Retain evidence across reloads. A removed/reimported helper cannot bless
# callbacks installed after package initialization.
if '_CAPTURED' not in globals():
    _CAPTURED=False
    _TYPE,_ALL,_ANY,_LEN,_ZIP,_DICT,_STR,_OBJECT=type,all,any,len,zip,dict,str,object
    _TUPLE,_LIST,_GETATTR=tuple,list,getattr
    _SCALARS=(bool,int,float,complex,str,bytes)


def register(ddm,lll):
    global _CAPTURED,_DDM,_DDM_STATE,_LLL,_MODULES,_FUNCTIONS,_LOOKUPS,_CLASS_FUNCTIONS,prefix,ready
    global _DOMAIN,_ELEMENT,_DOMAIN_DESCRIPTORS,_ELEMENT_DESCRIPTORS,_IMPORT,_FLOOR,_NATIVE_CONVERT
    global _BASIC,_BASIC_META,_BASIC_CHECK,_BASIC_CODE,_BASIC_CELLS,_ABC_CHECK,_ABC_CODE,_ABC_DELEGATE,_DUMP,_TOKEN
    global _RANK_ERROR,_RANK_MRO,_RANK_CLASSES,_RANK_FUNCTIONS,_MATRIX,_ABC_META,_ZZ,_QQ,_ELEMENT_NAME,_DDM_MODULE,_DENSE_MODULE,_PY_DIV,_PY_DIV_CODE,_PY_DIV_CLOSURE,_PY_DIV_CELLS,_NATIVE_CALL,_COLD,_BOOTSTRAP_STOCK
    if _CAPTURED:
        _register_reload(ddm,lll)
        return
    root=_sys.modules.get('symox')
    if _TYPE(root) is not ModuleType or not getattr(root.__dict__.get('__spec__'),'_initializing',False):return
    ddm_module=_sys.modules.get(ddm.__module__)
    dense_module=_sys.modules.get('symox.polys.matrices._dense_poly')
    if _TYPE(ddm_module) is not ModuleType or _TYPE(dense_module) is not ModuleType:return
    _DDM,_LLL,_DDM_MODULE,_DENSE_MODULE=ddm,lll,ddm_module,dense_module
    _DDM_STATE=ddm.__dict__['__dict__']
    _DOMAIN,_ELEMENT=_core.Domain,_core.DomainElement
    _MATRIX,_ABC_META=_core.Matrix,_abc.ABCMeta
    _RANK_ERROR=lll.DMRankError
    _RANK_MRO=_RANK_ERROR.__mro__
    _RANK_CLASSES=tuple((cls,tuple(cls.__dict__.items())) for cls in _RANK_MRO
                        if not _ANY(cls is immutable for immutable in (Exception,BaseException,object)))
    _RANK_FUNCTIONS=[]
    for cls,attributes in _RANK_CLASSES:
        for name,value in attributes:
            function=value.__func__ if _TYPE(value) in (staticmethod,classmethod) else value
            if _TYPE(function) is FunctionType:
                _RANK_FUNCTIONS.append((function,function.__code__,function.__defaults__,function.__kwdefaults__,
                    function.__closure__,tuple(c.cell_contents for c in function.__closure__ or ())))
    _ZZ,_QQ=lll.ZZ,lll.QQ
    _ELEMENT_NAME=_ELEMENT.__name__
    _DOMAIN_DESCRIPTORS=tuple((n,_DOMAIN.__dict__.get(n)) for n in ('of_type','zero','one','is_EXRAW','is_EX','is_ZZ_I','is_QQ_I'))
    _ELEMENT_DESCRIPTORS=tuple((n,_ELEMENT.__dict__.get(n)) for n in (
        '__add__','__radd__','__sub__','__rsub__','__mul__','__rmul__',
        '__pow__','__abs__','__float__','__floor__','__lt__','__le__','__gt__','__ge__','__bool__','parent',
        '__iadd__','__isub__','__imul__','__itruediv__','__ipow__'))
    _IMPORT,_FLOOR=_builtins.__import__,lll.mfloor
    _NATIVE_CONVERT=_methods._native_convert
    _NATIVE_CALL=_methods._native_call
    _PY_DIV=_ELEMENT.__dict__['__truediv__'];_PY_DIV_CODE=_PY_DIV.__code__
    _PY_DIV_CLOSURE=_PY_DIV.__closure__;_PY_DIV_CELLS=tuple(c.cell_contents for c in _PY_DIV_CLOSURE)

    _MODULES=(_core,_basic,_methods,_abc,_abc_native,_builtins,_sys,lll,ddm_module,dense_module)
    _FUNCTIONS=tuple((module,name,module.__dict__[name],module.__dict__[name].__code__,
                      module.__dict__[name].__closure__,module.__dict__[name].__defaults__,module.__dict__[name].__kwdefaults__)
        for module,names in ((_methods,('_call','convert_from','convert','_not_a_coeff','_inexact_value','_float_from_sympy')),
                             (dense_module,('_native_elements',))) for name in names)
    _CLASS_FUNCTIONS=[]
    for n in ('__init__','copy','zeros','eye'):
        descriptor=ddm.__dict__[n];f=getattr(descriptor,'__func__',descriptor)
        _CLASS_FUNCTIONS.append((n,descriptor,f.__code__,f.__defaults__,f.__kwdefaults__,f.__closure__,
                                 tuple(c.cell_contents for c in f.__closure__ or ())))
    _LOOKUPS=tuple((module.__dict__,n,module.__dict__.get(n,_builtins.__dict__.get(n)))
        for module in (lll,_methods,_basic,ddm_module,dense_module)
        for n in ('isinstance','getattr','hasattr','type','all','any','dict','set','frozenset','bool','float','complex',
                  'str','bytes','bytearray','list','tuple','int','len','range','min','max','enumerate','sum','abs','super'))
    _BASIC,_BASIC_META=_basic.Basic,_basic._BasicMeta
    _BASIC_CHECK=_BASIC_META.__instancecheck__;_BASIC_CODE=_BASIC_CHECK.__code__
    _BASIC_CELLS=tuple(cell.cell_contents for cell in _BASIC_CHECK.__closure__ or ())
    _ABC_CHECK=_abc.ABCMeta.__instancecheck__;_ABC_CODE=_ABC_CHECK.__code__
    _ABC_DELEGATE=_ABC_CHECK.__globals__['_abc_instancecheck']
    _DUMP,_TOKEN=_abc_native._get_dump,_abc_native.get_cache_token
    _COLD=_core.__dict__.get('_lll_cache_is_cold')
    auth=_core.__dict__.get('_lll_bootstrap_stock')
    # C names and module ownership cannot be forged by a Python wrapper.
    # This is immutable bootstrap evidence, not a memoized live proof: the
    # original lookup identity checks below still run on every admission.
    if not (_TYPE is BuiltinFunctionType.__class__ and _TYPE(_COLD) is BuiltinFunctionType
            and _COLD.__name__=='_lll_cache_is_cold' and _COLD.__self__ is _core):_COLD=None
    aliases=((None,'type',_TYPE),(None,'all',_ALL),(None,'any',_ANY),(None,'len',_LEN),
             (None,'zip',_ZIP),(None,'dict',_DICT),(None,'str',_STR),(None,'object',_OBJECT),
             (None,'tuple',_TUPLE),(None,'list',_LIST),(None,'getattr',_GETATTR))
    _BOOTSTRAP_STOCK=(_TYPE is BuiltinFunctionType.__class__ and _TYPE(auth) is BuiltinFunctionType
            and auth.__name__=='_lll_bootstrap_stock' and auth.__self__ is _core and auth(_LOOKUPS+aliases))
    _BOOTSTRAP_STOCK=(_BOOTSTRAP_STOCK and _TYPE(_abc_native) is ModuleType
            and _TYPE(_DUMP) is BuiltinFunctionType and _DUMP.__self__ is _abc_native
            and _DUMP.__name__=='_get_dump'
            and _TYPE(_TOKEN) is BuiltinFunctionType and _TOKEN.__self__ is _abc_native
            and _TOKEN.__name__=='get_cache_token')
    _CAPTURED=True
    # Bind actual admission to native code only after authentic bootstrap.
    # Validate every key before looking up either newly consumed C entry.
    if (_BOOTSTRAP_STOCK and _TYPE(_core) is ModuleType
            and _ALL(_TYPE(k) is _STR for k in _core.__dict__)):
        supported=_core.__dict__.get('_lll_first_supported')
        first=_core.__dict__.get('_lll_prefix_admit')
        after=_core.__dict__.get('_lll_ready_admit')
        if (_TYPE(supported) is BuiltinFunctionType and supported.__self__ is _core
                and supported.__name__=='_lll_first_supported'
                and _TYPE(first) is BuiltinFunctionType and first.__self__ is _core
                and first.__name__=='_lll_prefix_admit'
                and _TYPE(after) is BuiltinFunctionType and after.__self__ is _core
                and after.__name__=='_lll_ready_admit' and supported()):
            prefix,ready=first,after


def _protocol():
    if not _CAPTURED:return False
    if not _BOOTSTRAP_STOCK:return False
    if not _ALL(_TYPE(m) is ModuleType for m in _MODULES):return False
    modules=_sys.__dict__.get('modules')
    if (_TYPE(modules) is not _DICT
            or modules.get('symox.polys.matrices.lll') is not _LLL
            or modules.get('symox.polys.matrices.ddm') is not _DDM_MODULE):return False
    # A custom exception constructor may inspect the partially reduced local
    # matrices. These two tiny heap exception classes must remain stock in
    # every attribute; their inherited native exception bases are immutable.
    if _LLL.__dict__.get('DMRankError') is not _RANK_ERROR or _TYPE(_RANK_ERROR) is not _TYPE:return False
    current=_RANK_ERROR.__mro__
    if _LEN(current)!=_LEN(_RANK_MRO) or not _ALL(a is b for a,b in _ZIP(current,_RANK_MRO)):return False
    for cls,attributes in _RANK_CLASSES:
        if _TYPE(cls) is not _TYPE:return False
        state=cls.__dict__
        if (_LEN(state)!=_LEN(attributes) or not _ALL(_TYPE(n) is _STR for n in state)
                or not _ALL(state.get(n) is value for n,value in attributes)):return False
    for function,code,defaults,kwdefaults,closure,contents in _RANK_FUNCTIONS:
        if (function.__code__ is not code or function.__defaults__ is not defaults
                or function.__kwdefaults__ is not kwdefaults or function.__closure__ is not closure):return False
        try:
            if not _ALL(c.cell_contents is v for c,v in _ZIP(closure or (),contents)):return False
        except ValueError:return False
    if not _ALL(ns.get(n,_builtins.__dict__.get(n)) is value for ns,n,value in _LOOKUPS):return False
    if (_TYPE(_IMPORT) is not BuiltinFunctionType or _IMPORT.__module__!='builtins' or _IMPORT.__name__!='__import__'
            or _builtins.__dict__.get('__import__') is not _IMPORT):return False
    if (_LLL.__dict__.get('mfloor') is not _FLOOR or _TYPE(_FLOOR) is not BuiltinFunctionType
            or _FLOOR.__module__!='math' or _FLOOR.__name__!='floor'
            or _LLL.__dict__.get('ZZ') is not _ZZ or _LLL.__dict__.get('QQ') is not _QQ):return False
    if (_core.__dict__.get('Domain') is not _DOMAIN or _core.__dict__.get('DomainElement') is not _ELEMENT
            or _core.__dict__.get('Matrix') is not _MATRIX or _TYPE(_MATRIX) is not _TYPE
            or _basic.__dict__.get('_core') is not _core
            or _DDM_MODULE.__dict__.get('DDM') is not _DDM
            or _DDM_MODULE.__dict__.get('_native_elements') is not _DENSE_MODULE.__dict__.get('_native_elements')):return False
    if (_TYPE(_DDM) is not _TYPE or (_LEN(_DDM.__bases__)!=1 or _DDM.__bases__[0] is not _builtins.list)
            or '__getattribute__' in _DDM.__dict__ or '__getattr__' in _DDM.__dict__
            or _DDM.__dict__.get('__dict__') is not _DDM_STATE):return False
    if not _ALL(_DDM.__dict__.get(n,_builtins.list.__dict__[n]) is _builtins.list.__dict__[n]
                for n in ('__iter__','__getitem__','__setitem__','__len__','__new__')):return False
    if _ANY(name in _DDM.__dict__ for name in ('shape','rows','cols','domain','__del__')):return False
    for name,descriptor,code,defaults,kwdefaults,closure,contents in _CLASS_FUNCTIONS:
        if _DDM.__dict__.get(name) is not descriptor:return False
        f=descriptor.__func__ if _TYPE(descriptor) is classmethod else descriptor
        if (_TYPE(f) is not FunctionType or f.__code__ is not code or f.__globals__ is not _DDM_MODULE.__dict__
                or f.__defaults__ is not defaults or f.__kwdefaults__ is not kwdefaults or f.__closure__ is not closure):return False
        try:
            if not _ALL(c.cell_contents is v for c,v in _ZIP(closure or (),contents)):return False
        except ValueError:return False
    for cls in (_DOMAIN,_ELEMENT):
        if ((_LEN(cls.__bases__)!=1 or cls.__bases__[0] is not _OBJECT) or '__getattribute__' in cls.__dict__ or '__getattr__' in cls.__dict__
                or '__class__' in cls.__dict__):return False
    if '_symox_monogenic_extension' in _DOMAIN.__dict__ or _ELEMENT.__name__ is not _ELEMENT_NAME:return False
    if not _ALL(_DOMAIN.__dict__.get(n) is value for n,value in _DOMAIN_DESCRIPTORS):return False
    if not _ALL(_ELEMENT.__dict__.get(n) is value for n,value in _ELEMENT_DESCRIPTORS):return False
    # Authenticate the original native descriptors, not an arbitrary capture.
    if not _ALL(_TYPE(v) in (MethodDescriptorType,WrapperDescriptorType,GetSetDescriptorType)
                and v.__objclass__ is _DOMAIN and v.__name__==n for n,v in _DOMAIN_DESCRIPTORS):return False
    if not _ALL(v is None or (_TYPE(v) in (WrapperDescriptorType,MethodDescriptorType) and v.__objclass__ is _ELEMENT and v.__name__==n)
                for n,v in _ELEMENT_DESCRIPTORS):return False
    if (_DOMAIN.__dict__.get('__call__') is not _methods.__dict__.get('_call')
            or _methods.__dict__.get('_native_call') is not _NATIVE_CALL
            or _TYPE(_NATIVE_CALL) is not WrapperDescriptorType or _NATIVE_CALL.__objclass__ is not _DOMAIN
            or _NATIVE_CALL.__name__!='__call__'):return False
    if (_ELEMENT.__dict__.get('__truediv__') is not _PY_DIV or _PY_DIV.__code__ is not _PY_DIV_CODE
            or _PY_DIV.__closure__ is not _PY_DIV_CLOSURE):return False
    try:
        if not _ALL(c.cell_contents is v for c,v in _ZIP(_PY_DIV_CLOSURE,_PY_DIV_CELLS)):return False
    except ValueError:return False
    native,op,reflected=_PY_DIV_CELLS
    if (_TYPE(native) is not WrapperDescriptorType or native.__objclass__ is not _ELEMENT
            or native.__name__!='__truediv__' or _TYPE(op) is not _STR or op!='truediv' or reflected is not False):return False
    if (_DOMAIN.__dict__.get('convert') is not _methods.__dict__.get('convert')
            or _DOMAIN.__dict__.get('convert_from') is not _methods.__dict__.get('convert_from')
            or _methods.__dict__.get('_native_convert') is not _NATIVE_CONVERT
            or _TYPE(_NATIVE_CONVERT) is not MethodDescriptorType
            or _NATIVE_CONVERT.__objclass__ is not _DOMAIN or _NATIVE_CONVERT.__name__!='convert'):return False
    for module,name,function,code,closure,defaults,kwdefaults in _FUNCTIONS:
        f=module.__dict__.get(name)
        if (f is not function or _TYPE(f) is not FunctionType or f.__code__ is not code or f.__closure__ is not closure
                or f.__defaults__ is not defaults or f.__kwdefaults__ is not kwdefaults):return False
    return True


def _matrix(value):
    if _TYPE(value) is not _DDM:return False
    state=_DDM_STATE.__get__(value,_DDM)
    if not (_TYPE(state) is _DICT and _LEN(state)==4 and _ALL(_TYPE(k) is _STR for k in state)
            and _ALL(n in state for n in ('shape','rows','cols','domain'))):return False
    shape=state['shape']
    return (_TYPE(shape) is _builtins.tuple and _LEN(shape)==2
            and _ALL(_TYPE(v) is _builtins.int for v in shape)
            and _TYPE(state['rows']) is _builtins.int and _TYPE(state['cols']) is _builtins.int
            and state['rows']==shape[0] and state['cols']==shape[1]
            and _builtins.list.__len__(value)==shape[0] and _TYPE(state['domain']) is _DOMAIN)


def prefix(x,delta,transform):
    # A first conversion can invoke a classifier that captures scratch matrices
    # and restores itself. Certify it BEFORE construction, not after the hook.
    # Cold/invalid caches retain the complete old math for this public call.
    if not _CAPTURED or not _BOOTSTRAP_STOCK:return False
    if _TYPE(_core) is not ModuleType:return False
    if _COLD is not None:
        if _core.__dict__.get('_lll_cache_is_cold') is not _COLD:return False
        if _COLD(_BASIC,_ELEMENT,_DUMP,_abc_native) is True:return False
    return (_TYPE(transform) is _builtins.bool and _TYPE(delta) is _ELEMENT and ready(x,x,None))


def ready(x,y,transform):
    if not _protocol() or not _matrix(x) or not _matrix(y) or (transform is not None and not _matrix(transform)):return False
    if (_basic.__dict__.get('Basic') is not _BASIC or _basic.__dict__.get('_BasicMeta') is not _BASIC_META
            or _BASIC_META.__dict__.get('__instancecheck__') is not _BASIC_CHECK or _BASIC_CHECK.__code__ is not _BASIC_CODE):return False
    try:
        if not _ALL(c.cell_contents is v for c,v in _ZIP(_BASIC_CHECK.__closure__ or (),_BASIC_CELLS)):return False
    except ValueError:return False
    if (_abc.__dict__.get('ABCMeta') is not _ABC_META
            or _LEN(_BASIC_META.__bases__)!=1 or _BASIC_META.__bases__[0] is not _ABC_META
            or _ABC_META.__dict__.get('__instancecheck__') is not _ABC_CHECK or _ABC_CHECK.__code__ is not _ABC_CODE
            or _ABC_CHECK.__globals__.get('_abc_instancecheck') is not _ABC_DELEGATE
            or _abc_native.__dict__.get('_get_dump') is not _DUMP or _abc_native.__dict__.get('get_cache_token') is not _TOKEN):return False
    _,_,negative,version=_DUMP(_BASIC)
    if version!=_TOKEN() or not _ANY(ref() is _ELEMENT for ref in negative):return False
    root=_sys.modules.get('symox');mp=_sys.modules.get('mpmath')
    return (_TYPE(root) is ModuleType and root.__dict__.get('Basic') is _BASIC
            and _ALL(n in root.__dict__ for n in ('Expr','Set','Float','I','CoercionFailed','S'))
            and _TYPE(mp) is ModuleType and 'libmp' in mp.__dict__)


def apply(y,columns,delta,transform,first,rank_error):
    native=_core.__dict__.get('_ddm_lll_native')
    if native is None:return False
    result=native(y,columns,delta,transform,first,_DDM,rank_error)
    if result is None:return False
    basis,rows=result
    y[:]=basis
    if transform is not None:transform[:]=rows
    return True


# A real reload recompiles methods. Compare immutable bytecode/constants and
# the actual global/closure/delegate environment, never __wrapped__ or names.
def _same_constant(a,b):
    if _TYPE(a) is not _TYPE(b):return False
    if _TYPE(a) is CodeType:return _same_code(a,b)
    if _TYPE(a) is _TUPLE:return _LEN(a)==_LEN(b) and _ALL(_same_constant(x,y) for x,y in _ZIP(a,b))
    if a is None:return True
    if _ANY(_TYPE(a) is t for t in _SCALARS):return a==b
    return False


def _same_code(a,b):
    return _ALL(_same_constant(_GETATTR(a,n),_GETATTR(b,n)) for n in (
        'co_code','co_consts','co_names','co_varnames','co_freevars','co_cellvars','co_argcount',
        'co_posonlyargcount','co_kwonlyargcount','co_nlocals','co_stacksize','co_flags'))


def _register_reload(ddm,lll):
    global _DDM,_DDM_STATE,_CLASS_FUNCTIONS
    if ddm is _DDM:return
    if (_TYPE(ddm) is not _TYPE or lll is not _LLL or _DDM_MODULE.__dict__.get('DDM') is not ddm
            or _LEN(ddm.__bases__)!=1 or ddm.__bases__[0] is not _LIST):return
    state=ddm.__dict__.get('__dict__')
    if (_TYPE(state) is not GetSetDescriptorType or state.__objclass__ is not ddm or state.__name__!='__dict__'):return
    methods=[]
    for name,old,code,defaults,kwdefaults,closure,contents in _CLASS_FUNCTIONS:
        current=ddm.__dict__.get(name)
        if _TYPE(current) is not _TYPE(old):return
        f=current.__func__ if _TYPE(current) is classmethod else current
        if (_TYPE(f) is not FunctionType or f.__globals__ is not _DDM_MODULE.__dict__
                or f.__builtins__ is not _builtins.__dict__ or not _same_code(f.__code__,code)
                or f.__defaults__ is not defaults or f.__kwdefaults__ is not kwdefaults):return
        actual=f.__closure__ or ()
        if _LEN(actual)!=_LEN(contents):return
        try:
            if not _ALL(c.cell_contents is (ddm if v is _DDM else v) for c,v in _ZIP(actual,contents)):return
        except ValueError:return
        methods.append((name,current,f.__code__,f.__defaults__,f.__kwdefaults__,f.__closure__,
                        tuple(c.cell_contents for c in actual)))
    _DDM,_DDM_STATE,_CLASS_FUNCTIONS=ddm,state,methods


# Preserve the established guard functions for explicit validation controls.
# Normal public selection uses the native bindings above.
_LEGACY_PREFIX,_LEGACY_READY=prefix,ready
