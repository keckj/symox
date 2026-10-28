"""Admission for the ordinary exact ZZ HNF loop; no matrix arithmetic here."""
import builtins as _builtins
from types import ModuleType,WrapperDescriptorType,MethodDescriptorType,GetSetDescriptorType
from symox import _core
from ..domains import _algebraic_methods as _am,_element_print as _ep
from .ddm import DDM

def _expected_gcdex(a, b):
    r"""
    This supports the functions that compute Hermite Normal Form.

    Explanation
    ===========

    Let x, y be the coefficients returned by the extended Euclidean
    Algorithm, so that x*a + y*b = g. In the algorithms for computing HNF,
    it is critical that x, y not only satisfy the condition of being small
    in magnitude -- namely that |x| <= |b|/g, |y| <- |a|/g -- but also that
    y == 0 when a | b.

    """
    x, y, g = ZZ.gcdex(a, b)
    if a != 0 and b % a == 0:
        y = 0
        x = -1 if a < 0 else 1
    return x, y, g

def _expected_add_columns(m, i, j, a, b, c, d):
    # replace m[:, i] by a*m[:, i] + b*m[:, j]
    # and m[:, j] by c*m[:, i] + d*m[:, j]
    for k in range(len(m)):
        e = m[k][i]
        m[k][i] = a*e + b*m[k][j]
        m[k][j] = c*e + d*m[k][j]

if '_CAPTURED' not in globals():
    _CAPTURED=False
    _TYPE=False.__class__.__class__
    _STR=''.__class__
    _BOOL=False.__class__
    _INT=(0).__class__
    _BYTES=b''.__class__
    _LIST=[].__class__
    _TUPLE=().__class__
    _DICT={}.__class__
    _INT=(0).__class__
    _FUNCTION=(lambda:None).__class__
    _ANCHORS_PLAIN=True
    for _module in (_core,_am,_ep,_builtins):
        if _TYPE(_module) is not ModuleType:
            _ANCHORS_PLAIN=False;break
        for _key in _module.__dict__:
            if _TYPE(_key) is not _STR:
                _ANCHORS_PLAIN=False;break
        if not _ANCHORS_PLAIN:break
    _READ_FUNCTION=(_core.__dict__.get('_ddm_function_metadata') if _ANCHORS_PLAIN else None)
    _BUILTIN_FUNCTION=[].append.__class__
    _METADATA_OK=(_ANCHORS_PLAIN and _TYPE(_am) is ModuleType and _TYPE(_ep) is ModuleType
                  and _TYPE(_READ_FUNCTION) is _BUILTIN_FUNCTION
                  and _READ_FUNCTION.__self__ is _core and _READ_FUNCTION.__name__=='_ddm_function_metadata')
    _BUILTIN=[].append.__class__
    _OBJECT=_LIST.__bases__[0]
    _DOMAIN,_ELEMENT=(_core.Domain,_core.DomainElement) if _ANCHORS_PLAIN else (None,None)
    _ZZ=_core.ZZ if _ANCHORS_PLAIN else None
    _RANGE,_LEN=range,len
    _DDM=DDM
    _HOOKS=[(name,_ELEMENT.__dict__.get(name)) for name in (('__lt__','__neg__','__mul__','__rmul__','__add__','parent') if _TYPE(_ELEMENT) is _TYPE else ())]
    _DOMAIN_HOOKS=[(name,_DOMAIN.__dict__.get(name)) for name in (('gcdex','is_FiniteField','is_ZZ_I','is_QQ_I') if _TYPE(_DOMAIN) is _TYPE else ())]
    _WRAPPER_CODE=None
    _factory_metadata=_READ_FUNCTION(_am._gaussian_binary) if _METADATA_OK else None
    _factory_code=_factory_metadata[0] if _factory_metadata is not None else None
    _CODE=_factory_code.__class__
    for _constant in _factory_code.co_consts if _factory_code is not None else ():
        if _TYPE(_constant) is _CODE and _constant.co_name=='method':_WRAPPER_CODE=_constant
    _congruent=_ep._int_congruent if _ANCHORS_PLAIN else None
    _PRINT_FUNCTIONS=[]
    for _name in ('_int_congruent','_eq','_ne') if _ANCHORS_PLAIN else ():
        _function=_ep.__dict__[_name]
        _metadata=_READ_FUNCTION(_function) if _METADATA_OK else None
        if _metadata is not None:_PRINT_FUNCTIONS.append((_name,_function,_metadata[0]))
    _PRINT_BUILTINS=_ep.__dict__.get('__builtins__') if _ANCHORS_PLAIN else None
    # The two independent literal templates above are never executed. A
    # forged/reimported peer module cannot supply ordinary-HNF's authority.
    _NORMALFORM_CODES=()
    _CODE_READERS=()
    if _METADATA_OK and _factory_code is not None:
        _gcd_code=_READ_FUNCTION(_expected_gcdex)
        _columns_code=_READ_FUNCTION(_expected_add_columns)
        if _gcd_code is not None and _columns_code is not None:
            _NORMALFORM_CODES=(('_gcdex',_gcd_code[0]),('add_columns',_columns_code[0]))
            _CODE_READERS=_TUPLE(_CODE.__dict__[_name].__get__ for _name in
                ('co_argcount','co_posonlyargcount','co_kwonlyargcount','co_nlocals','co_stacksize','co_flags',
                 'co_code','co_consts','co_names','co_varnames','co_freevars','co_cellvars'))


def _core_plain():
    if _TYPE(_core) is not ModuleType:return False
    for key in _core.__dict__:
        if _TYPE(key) is not _STR:return False
    return True


def _function_metadata(function):
    if (not _METADATA_OK or not _core_plain()
            or _core.__dict__.get('_ddm_function_metadata') is not _READ_FUNCTION
            or _TYPE(function) is not _FUNCTION or _TYPE(_builtins) is not ModuleType
            or function.__builtins__ is not _builtins.__dict__):return None
    return _READ_FUNCTION(function)


def register(namespace):
    global _CAPTURED,_NAMESPACE,_NS_BUILTINS,_FUNCTIONS
    if _CAPTURED:return
    if _WRAPPER_CODE is None or _PRINT_FUNCTIONS.__len__()!=3:return
    if _TYPE(namespace) is not _DICT:return
    for key in namespace:
        if _TYPE(key) is not _STR:return
    funcs=[]
    for name in ('_gcdex','add_columns'):
        function=namespace.get(name)
        if _TYPE(function) is not _FUNCTION:return
        metadata=_function_metadata(function)
        if metadata is None or not _normalform_code(metadata[0],name):return
        funcs.append((name,function,metadata[0]))
    _FUNCTIONS=funcs;_NAMESPACE=namespace;_NS_BUILTINS=namespace.get('__builtins__');_CAPTURED=True


def _method(value,cls,name,kinds):
    return _TYPE(value) in kinds and value.__objclass__ is cls and value.__name__==name


def _builtins_entry(entry):
    return entry is _builtins or entry is _builtins.__dict__


def _plain(function,code,namespace):
    metadata=_function_metadata(function)
    return (metadata is not None and metadata[0] is code and metadata[1] is namespace
            and metadata[2] is None and metadata[3] is None and metadata[4] is None)


def _wrapper(name):
    function=_ELEMENT.__dict__.get(name)
    metadata=_function_metadata(function)
    if (metadata is None or metadata[0] is not _WRAPPER_CODE
            or metadata[1] is not _am.__dict__ or metadata[2] is not None
            or metadata[3] is not None):return False
    closure=metadata[4]
    if _TYPE(closure) is not _TUPLE or closure.__len__()!=3:return False
    try:native,op,reflected=(cell.cell_contents for cell in closure)
    except ValueError:return False
    return (_method(native,_ELEMENT,name,(WrapperDescriptorType,)) and _TYPE(op) is ''.__class__
            and op==name.strip('_') and reflected is False)


def try_hnf(matrix,dimensions,dependencies):
    if not _CAPTURED or not _core_plain() or _TYPE(dependencies) is not _TUPLE or dependencies.__len__()!=10:return None
    for module in (_core,_am,_ep,_builtins):
        if _TYPE(module) is not ModuleType:return None
    if (_core.__dict__.get('Domain') is not _DOMAIN or _core.__dict__.get('DomainElement') is not _ELEMENT
            or _TYPE(_ZZ) is not _DOMAIN):return None
    if (_TYPE(_RANGE) is not _TYPE or _RANGE.__flags__ & (1<<9)
            or _RANGE.__name__!='range' or _RANGE.__module__!='builtins'
            or _TYPE(_LEN) is not _BUILTIN or _LEN.__self__ is not _builtins or _LEN.__name__!='len'):return None
    expected=(_RANGE,_LEN,_ZZ,_FUNCTIONS[0][1],_FUNCTIONS[1][1],_ZZ,_RANGE,_LEN,_TYPE,_INT)
    index=0
    for value in expected:
        if dependencies[index] is not value:return None
        index+=1
    if (_NAMESPACE.get('__builtins__') is not _NS_BUILTINS or not _builtins_entry(_NS_BUILTINS)
            or _ep.__dict__.get('__builtins__') is not _PRINT_BUILTINS or not _builtins_entry(_PRINT_BUILTINS)):return None
    for name,function,code in _FUNCTIONS:
        if not _plain(function,code,_NAMESPACE):return None
    for name,function,code in _PRINT_FUNCTIONS:
        if _ep.__dict__.get(name) is not function or not _plain(function,code,_ep.__dict__):return None
    if (_ELEMENT.__dict__.get('__eq__') is not _PRINT_FUNCTIONS[1][1]
            or _ELEMENT.__dict__.get('__ne__') is not _PRINT_FUNCTIONS[2][1]
            or not _method(_ep._native_eq,_ELEMENT,'__eq__',(WrapperDescriptorType,))
            or not _method(_ep._native_ne,_ELEMENT,'__ne__',(WrapperDescriptorType,))):return None
    for cls in (_DOMAIN,_ELEMENT):
        if _TYPE(cls) is not _TYPE or cls.__bases__!=(_OBJECT,):return None
        for name in ('__getattribute__','__getattr__','__class__'):
            if name in cls.__dict__:return None
    for name,value in _HOOKS:
        kinds=(MethodDescriptorType,) if name=='parent' else (WrapperDescriptorType,)
        if _ELEMENT.__dict__.get(name) is not value or not _method(value,_ELEMENT,name,kinds):return None
    for name,value in _DOMAIN_HOOKS:
        kinds=(MethodDescriptorType,) if name=='gcdex' else (GetSetDescriptorType,)
        if _DOMAIN.__dict__.get(name) is not value or not _method(value,_DOMAIN,name,kinds):return None
    if not _wrapper('__floordiv__') or not _wrapper('__mod__'):return None
    if _TYPE(_DDM) is not _TYPE or _DDM.__bases__!=(_LIST,):return None
    for name in ('__iter__','__len__','__getitem__','__setitem__','__new__'):
        if _DDM.__dict__.get(name,_LIST.__dict__[name]) is not _LIST.__dict__[name]:return None
    return _core._ddm_hnf(matrix,dimensions,_DDM)


def _literal_equal(left,right):
    kind=_TYPE(left)
    if kind is not _TYPE(right):return False
    if left is None:return True
    if kind is _BOOL or kind is _INT or kind is _STR or kind is _BYTES:return left==right
    if kind is _TUPLE:
        if left.__len__()!=right.__len__():return False
        i=0
        for value in left:
            if not _literal_equal(value,right[i]):return False
            i+=1
        return True
    return False

def _normalform_code(code,name):
    if _TYPE(code) is not _CODE or not _CODE_READERS:return False
    for key,expected in _NORMALFORM_CODES:
        if key==name:
            for reader in _CODE_READERS:
                if not _literal_equal(reader(code),reader(expected)):return False
            return True
    return False
