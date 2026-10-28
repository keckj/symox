"""Protocol-only selection of exact QQ sparse field Gauss--Jordan."""
from types import WrapperDescriptorType, BuiltinFunctionType
import builtins as _builtins
from symox import _core

_TYPE=type
_ALL=all
_LEN=len
_ANY=any
_DICTTYPE=dict
_OBJECT=object
_DICT=type.__dict__['__dict__']
_BASES=type.__dict__['__bases__']
_BUILTINS=_builtins.__dict__
_ELEMENT=_core.DomainElement
_OPERATIONS=('__pow__','__mul__','__sub__','__neg__','__bool__')
_DESCRIPTORS=tuple((name,_DICT.__get__(_ELEMENT).get(name)) for name in _OPERATIONS)
_NATIVE=all(type(value) is WrapperDescriptorType and value.__objclass__ is _ELEMENT
            and value.__name__==name for name,value in _DESCRIPTORS)
_PARENT_BASES=_BASES.__get__(_ELEMENT)
_MARKERS=('_symox_domainmatrix','_symox_quotient_element','is_Matrix',
          '_symox_sparse','_symox_ring_element')
_PROTOCOL=globals().get('_PROTOCOL')


def _stock_aliases(aliases):
    """Authenticate the initial aliases, before trusting identity snapshots."""
    for _,name,value in aliases:
        if name in ('min','sorted','len','all','any'):
            if (_TYPE(value) is not BuiltinFunctionType
                    or value.__self__ is not _builtins
                    or value.__name__!=name):return False
        else:
            # These factories are static CPython types. Heap types, including
            # subclasses carrying copied names/modules, are custom protocols.
            owner='collections' if name=='defaultdict' else 'builtins'
            if (_TYPE(value) is not _TYPE or value.__flags__ & 512
                    or value.__name__!=name or value.__module__!=owner):return False
    return True


def register_irref(stock,namespace):
    """Keep first bootstrap evidence; a reload cannot bless new callbacks."""
    global _PROTOCOL
    if _PROTOCOL is not None and _PROTOCOL[2][0][0] is not namespace:
        # A fresh module owns different global lookups. Never authorize it
        # using snapshots of a module removed from sys.modules. Ordinary
        # reload reuses its dictionary and keeps the existing proof.
        return _decline_irref
    if _PROTOCOL is None:
        names=('sorted','min','set','dict','enumerate','len')
        aliases=tuple((namespace,name,namespace.get(name,_BUILTINS[name])) for name in names)
        aliases+=tuple((_BUILTINS,name,_BUILTINS[name]) for name in names)
        # default factory is read even on the empty legacy branch.
        aliases+=((namespace,'defaultdict',namespace['defaultdict']),)
        # Selection itself must not invoke a guard helper captured while it
        # was replaced: the general algorithm does not consume all or any.
        bootstrap=_stock_aliases(aliases+((namespace,'all',_ALL),(namespace,'any',_ANY)))
        _PROTOCOL=(stock,_DICT.__get__(stock)['__dict__'],aliases,bootstrap)
    return native_irref


def _decline_irref(matrix):
    return None


def native_irref(matrix):
    protocol=_PROTOCOL
    if protocol is None or not _NATIVE:return None
    stock,state,aliases,bootstrap=protocol
    if not bootstrap:return None
    if not _ALL(module.get(name,_BUILTINS.get(name)) is value for module,name,value in aliases):return None
    current=_DICT.__get__(_ELEMENT);parents=_BASES.__get__(_ELEMENT)
    if (_LEN(parents)!=_LEN(_PARENT_BASES)
            or _LEN(parents)!=1 or parents[0] is not _OBJECT
            or not _ALL(current.get(name) is value for name,value in _DESCRIPTORS)
            or _ANY(name in current for name in ('__getattribute__','__getattr__','__imul__')+_MARKERS)):return None
    if _TYPE(matrix) is not _DICTTYPE:
        if _TYPE(matrix) is not stock:return None
        parent=_BASES.__get__(stock);methods=_DICT.__get__(stock)
        if (_LEN(parent)!=1 or parent[0] is not _DICTTYPE
                or '__getattribute__' in methods or '__getattr__' in methods
                or methods.get('values',_DICTTYPE.__dict__['values']) is not _DICTTYPE.__dict__['values']):return None
        # An instance values override is observed by A.values on the old path.
        attrs=state.__get__(matrix,stock)
        if (_TYPE(attrs) is not _DICTTYPE or not _ALL(_TYPE(k) is str for k in attrs)
                or 'values' in attrs):return None
    # The binding checks actual containers/keys, exact typed QQ payloads and
    # explicit stored zeros without invoking any foreign scalar protocol.
    # Gaussian division wrappers, Domain attributes, Basic/ABC caches and
    # imports are not consumed by QQ pow(-1)/mul/sub/neg/bool; no such probe.
    return _core._sdm_irref(matrix,stock)
