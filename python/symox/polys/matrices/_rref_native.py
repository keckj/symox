"""Protocol-only selection of stock ZZ/QQ AUTO statistics in Rust."""
from types import FunctionType, MethodDescriptorType, GetSetDescriptorType, WrapperDescriptorType
from marshal import dumps as _CODE
import builtins as _builtins
from symox import _core

_TYPE=type
_ALL=all
_LEN=len
_TUPLE=tuple
_BUILTINS=_builtins.__dict__
_DICT=type.__dict__['__dict__']
_BASES=type.__dict__['__bases__']
_DOMAIN=_core.Domain
_ELEMENT=_core.DomainElement
_PROTOCOL=globals().get('_PROTOCOL')
_NATIVE_BASES=globals().get('_NATIVE_BASES',((_DOMAIN,_BASES.__get__(_DOMAIN)),(_ELEMENT,_BASES.__get__(_ELEMENT))))


def _keywords(current,entries):
    if entries is None:return current is None
    return (_TYPE(current) is dict and _LEN(current)==_LEN(entries)
            and _ALL(_TYPE(key) is str for key in current)
            and _ALL(current.get(key) is value for key,value in entries))


def _function(function):
    if _TYPE(function) is not FunctionType:return None
    keywords=function.__kwdefaults__
    return (function,function.__code__,function.__defaults__,
            None if keywords is None else _TUPLE(keywords.items()),
            _TUPLE(c.cell_contents for c in function.__closure__ or ()),
            _CODE(function.__code__,0))


def _same_function(current,state,reloaded=False):
    function,code,defaults,keywords,closure,serialized=state
    if (_TYPE(current) is not FunctionType or current.__globals__ is not function.__globals__
            or current.__defaults__ is not defaults
            or not _keywords(current.__kwdefaults__,keywords)):return False
    try:
        return ((current.__code__ is code or (reloaded and _CODE(current.__code__,0)==serialized))
                and _LEN(current.__closure__ or ())==_LEN(closure)
                and _ALL(cell.cell_contents is old for cell,old in zip(current.__closure__ or (),closure)))
    except (ValueError,TypeError):return False


# Capture only native descriptors of methods actually used by the heuristic.
# An optional package can add unrelated methods without invalidating selection.
_DOMAIN_NAMES=('is_ZZ','is_QQ','one','lcm','of_type')
_ELEMENT_NAMES=('__bool__','__int__','numerator','denominator')
_DOMAIN_HOOKS=_TUPLE((name,_DICT.__get__(_DOMAIN).get(name)) for name in _DOMAIN_NAMES)
_ELEMENT_HOOKS=_TUPLE((name,_DICT.__get__(_ELEMENT).get(name)) for name in _ELEMENT_NAMES)
_NATIVE=(_ALL(_TYPE(value) in (MethodDescriptorType,GetSetDescriptorType)
             and value.__objclass__ is _DOMAIN and value.__name__==name for name,value in _DOMAIN_HOOKS)
         and _ALL(_TYPE(value) in (WrapperDescriptorType,GetSetDescriptorType)
                  and value.__objclass__ is _ELEMENT and value.__name__==name for name,value in _ELEMENT_HOOKS))


def register_auto_protocol(dm,ddm,sdm,namespace):
    """Capture completed stock matrix protocols, before user callbacks."""
    global _PROTOCOL
    if _PROTOCOL is not None:return
    modules=(namespace,__import__(ddm.__module__,fromlist=['']).__dict__,
             __import__(sdm.__module__,fromlist=['']).__dict__)
    element_function=modules[2]['_native_elements']
    modules+=(element_function.__globals__,)
    records=[]
    for cls,base,names in (
            (dm,object,('to_flat_nz',)),
            (ddm,list,('to_sdm','to_flat_nz','__iter__','__len__','__getitem__')),
            (sdm,dict,('to_sdm','from_list','to_flat_nz','to_dok','__init__','values','items'))):
        names+=('__getattribute__','__getattr__')
        entries=_TUPLE((name,_DICT.__get__(cls).get(name)) for name in names)
        functions=[]
        for _,value in entries:
            function=value.__func__ if _TYPE(value) is classmethod else value
            state=_function(function)
            if state is not None:functions.append(state)
        records.append((cls,base,entries,_TUPLE(functions),_DICT.__get__(cls)['__dict__']))
    names=('_dm_rref_choose_method_ZZ','_dm_rref_choose_method_QQ','_dm_row_density','_dm_elements','_dm_QQ_numers_denoms')
    functions=_TUPLE((name,_function(namespace[name])) for name in names)
    aliases=((namespace,'ZZ',namespace['ZZ']),)
    aliases+=_TUPLE((module,name,module.get(name)) for module in modules[1:3] for name in ('SDM','DDM','_native_elements'))
    builtin_names=('sum','map','len','max','min','int','list','dict','tuple','range','enumerate','all','any','zip','type','isinstance','next','super','TypeError','AttributeError')
    lookups=_TUPLE((module,name,module.get(name,_BUILTINS[name])) for module in modules for name in builtin_names)
    lookups+=_TUPLE((_BUILTINS,name,_BUILTINS[name]) for name in builtin_names+('__import__',))
    _PROTOCOL=(_TUPLE(records),namespace,functions,aliases+lookups,_function(element_function))


def native_choose_method(matrix):
    protocol=_PROTOCOL
    if protocol is None or not _NATIVE:return None
    records,namespace,functions,lookups,element_function=protocol
    dm=records[0][0]
    if _TYPE(matrix) is not dm:return None
    # Check builtin aliases without importing or invoking a replaced builtin.
    if not _ALL(module.get(name,_BUILTINS.get(name)) is value for module,name,value in lookups):return None
    if not _ALL(_same_function(namespace.get(name),state,True) for name,state in functions):return None
    if not _same_function(element_function[0],element_function):return None
    for cls,base,entries,states,_ in records:
        current=_DICT.__get__(cls);parents=_BASES.__get__(cls)
        if (_LEN(parents)!=1 or parents[0] is not base
                or any(name in current for name in ('rep','shape','domain'))
                or not _ALL(current.get(name) is value for name,value in entries)
                or not _ALL(_same_function(state[0],state) for state in states)):return None
    for cls,expected in _NATIVE_BASES:
        parents=_BASES.__get__(cls)
        if _LEN(parents)!=_LEN(expected) or not _ALL(a is b for a,b in zip(parents,expected)):return None
    for cls,entries in ((_DOMAIN,_DOMAIN_HOOKS),(_ELEMENT,_ELEMENT_HOOKS)):
        current=_DICT.__get__(cls)
        if ('__getattribute__' in current or '__getattr__' in current
                or not _ALL(current.get(name) is value for name,value in entries)):return None
    # Dense conversion fetches convert but native of_type succeeds for every
    # admitted, already typed scalar. Fetching a function/native descriptor is
    # inert; a custom descriptor must retain its observable __get__ instead.
    convert=_DICT.__get__(_DOMAIN).get('convert')
    if (_TYPE(convert) is not FunctionType
            and not (_TYPE(convert) is MethodDescriptorType
                     and convert.__objclass__ is _DOMAIN and convert.__name__=='convert')):return None
    attrs=records[0][4].__get__(matrix,dm)
    if (_TYPE(attrs) is not dict or not _ALL(_TYPE(key) is str for key in attrs)
            or not _ALL(name in attrs for name in ('rep','shape','domain'))
            or not _ALL(name not in attrs for name,_ in records[0][2])):return None
    rep,shape,domain=attrs['rep'],attrs['shape'],attrs['domain']
    if _TYPE(domain) is not _DOMAIN:return None
    info=next((record for record in records[1:] if _TYPE(rep) is record[0]),None)
    if info is None:return None
    rep_attrs=info[4].__get__(rep,info[0])
    if (_TYPE(rep_attrs) is not dict or not _ALL(_TYPE(key) is str for key in rep_attrs)
            or rep_attrs.get('shape') is not shape or rep_attrs.get('domain') is not domain
            or not _ALL(name not in rep_attrs for name,_ in info[2])):return None
    return _core._dm_rref_auto(rep,shape,domain,info[0])
