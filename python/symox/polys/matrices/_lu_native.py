"""Protocol admission for QQ LU; all elimination arithmetic lives in Rust."""
import builtins as _builtins
from types import ModuleType, WrapperDescriptorType, MethodDescriptorType, GetSetDescriptorType
from symox import _core
from ..domains import _algebraic_methods

# Reloads cannot bless changes made after the package's first initialization.
if '_CAPTURED' not in globals():
    _CAPTURED=True
    _TYPE=False.__class__.__class__
    _STR=''.__class__
    _LIST=[].__class__
    _TUPLE=().__class__
    _OBJECT=_LIST.__bases__[0]
    _FUNCTION=(lambda:None).__class__
    _ANCHORS_PLAIN=True
    for _module in (_core,_algebraic_methods):
        if _TYPE(_module) is not ModuleType:
            _ANCHORS_PLAIN=False;break
        for _key in _module.__dict__:
            if _TYPE(_key) is not _STR:
                _ANCHORS_PLAIN=False;break
        if not _ANCHORS_PLAIN:break
    _READ_FUNCTION=(_core.__dict__.get('_ddm_function_metadata') if _ANCHORS_PLAIN else None)
    _BUILTIN_FUNCTION=[].append.__class__
    _METADATA_OK=(_ANCHORS_PLAIN and _TYPE(_algebraic_methods) is ModuleType
                  and _TYPE(_READ_FUNCTION) is _BUILTIN_FUNCTION
                  and _READ_FUNCTION.__self__ is _core and _READ_FUNCTION.__name__=='_ddm_function_metadata')
    _BUILTIN=[].append.__class__
    _DOMAIN,_ELEMENT=(_core.Domain,_core.DomainElement) if _ANCHORS_PLAIN else (None,None)
    _ELEMENT_HOOKS=[(name,_ELEMENT.__dict__.get(name)) for name in
                        (('__bool__','__mul__','__sub__','parent') if _TYPE(_ELEMENT) is _TYPE else ())]
    _DOMAIN_HOOKS=[(name,_DOMAIN.__dict__.get(name)) for name in (('is_ZZ_I','is_QQ_I') if _TYPE(_DOMAIN) is _TYPE else ())]
    _LEN,_RANGE,_MIN=len,range,min
    _DIV_CODE=None
    _factory_metadata=_READ_FUNCTION(_algebraic_methods._gaussian_binary) if _METADATA_OK else None
    _factory_code=_factory_metadata[0] if _factory_metadata is not None else None
    _CODE=_factory_code.__class__
    for _constant in _factory_code.co_consts if _factory_code is not None else ():
        if _TYPE(_constant) is _CODE and _constant.co_name=='method':
            _DIV_CODE=_constant


def _core_plain():
    if _TYPE(_core) is not ModuleType:return False
    for key in _core.__dict__:
        if _TYPE(key) is not _STR:return False
    return True


def _function_metadata(function):
    if (not _METADATA_OK or not _core_plain()
            or _core.__dict__.get('_ddm_function_metadata') is not _READ_FUNCTION
            or _TYPE(function) is not _FUNCTION):return None
    return _READ_FUNCTION(function)


def _builtin(value,name):
    return (_TYPE(value) is _BUILTIN and value.__self__ is _builtins
            and value.__name__==name)


def _native_method(value,cls,name,types):
    return (_TYPE(value) in types and value.__objclass__ is cls and value.__name__==name)


def _protocol(dependencies):
    if not _core_plain():return False
    if (_TYPE(_core) is not ModuleType or _TYPE(_algebraic_methods) is not ModuleType
            or _TYPE(_builtins) is not ModuleType):return False
    if not _builtin(_LEN,'len') or not _builtin(_MIN,'min'):return False
    if (_TYPE(_RANGE) is not _TYPE or _RANGE.__name__!='range'
            or _RANGE.__module__!='builtins' or _RANGE.__flags__ & (1<<9)):return False
    if (_TYPE(dependencies) is not _TUPLE or dependencies.__len__()!=2
            or dependencies[0] is not _RANGE or dependencies[1] is not _MIN):return False
    if (_core.__dict__.get('Domain') is not _DOMAIN
            or _core.__dict__.get('DomainElement') is not _ELEMENT):return False
    for cls in (_DOMAIN,_ELEMENT):
        if _TYPE(cls) is not _TYPE or cls.__bases__!=(_OBJECT,):return False
        for name in ('__getattribute__','__getattr__','__class__'):
            if name in cls.__dict__:return False
    if '__isub__' in _ELEMENT.__dict__:return False
    for name,original in _ELEMENT_HOOKS:
        if _ELEMENT.__dict__.get(name) is not original:return False
        kinds=(MethodDescriptorType,) if name=='parent' else (WrapperDescriptorType,)
        if not _native_method(original,_ELEMENT,name,kinds):return False
    for name,original in _DOMAIN_HOOKS:
        if (_DOMAIN.__dict__.get(name) is not original
                or not _native_method(original,_DOMAIN,name,(GetSetDescriptorType,))):return False
    divide=_ELEMENT.__dict__.get('__truediv__')
    metadata=_function_metadata(divide)
    if (_DIV_CODE is None or metadata is None or metadata[0] is not _DIV_CODE
            or metadata[1] is not _algebraic_methods.__dict__
            or metadata[2] is not None or metadata[3] is not None):return False
    closure=metadata[4]
    if _TYPE(closure) is not _TUPLE or closure.__len__()!=3:return False
    try:
        native,operation,reflected=(cell.cell_contents for cell in closure)
    except ValueError:return False
    return (_native_method(native,_ELEMENT,'__truediv__',(WrapperDescriptorType,))
            and _TYPE(operation) is ''.__class__ and operation=='truediv' and reflected is False)


def try_ilu(matrix,dependencies,stock,dimensions=None):
    if not _protocol(dependencies):return None
    if stock is not None:
        if _TYPE(stock) is not _TYPE or stock.__bases__!=(_LIST,):return None
        for name in ('__iter__','__len__','__getitem__','__setitem__','__new__'):
            if stock.__dict__.get(name,_LIST.__dict__[name]) is not _LIST.__dict__[name]:return None
    return _core._ddm_ilu(matrix,stock,dimensions)
