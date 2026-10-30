"""Protocol-only selection of native exact PolyMatrix ground reduction."""
from types import FunctionType,CodeType,ModuleType
from sys import modules as _MODULES
from marshal import dumps as _CODE
from .._native_dispatch import strip_slot_cache as _strip_slot_cache
from .. import _core

_TYPE=type
_ALL=all
_LEN=len
_TUPLE=tuple
_ROOT=None
_FORWARDING=None
_FORWARDING_METADATA=None
_STATS=None
_HOST_FREE_SYMBOLS=None
_METADATA=None
_DIFFGEOM_AUX=None
_RATIONAL_BUILTINS=None
_RELOADABLE={}
_TYPE_DICT=type.__dict__['__dict__']
_TYPE_BASES=type.__dict__['__bases__']
_TYPE_MRO=type.__dict__['__mro__']
_STR=''.__class__
_DICT={}.__class__
_LIST=[].__class__
_BUILTIN=[].append.__class__


_STR_KEYED=_core.__dict__.get('_str_keyed')


def _plain_namespace(namespace):
    # Exact dict whose keys are all exact str (native key scan when built).
    if _STR_KEYED is not None:return _STR_KEYED(namespace)
    if _TYPE(namespace) is not _DICT:return False
    for key in namespace:
        if _TYPE(key) is not _STR:return False
    return True


def _metadata_reader():
    if _METADATA is None:return None
    owner,reader=_METADATA
    if (_TYPE(owner) is not ModuleType or _core is not owner
            or not _plain_namespace(owner.__dict__)
            or owner.__dict__.get('_ddm_function_metadata') is not reader):return None
    return reader



_FUNCTION_ATTRIBUTE=(lambda:None).__class__.__getattribute__
_FUNCTION_FIELDS=('__code__','__globals__','__defaults__','__kwdefaults__','__closure__')


def _metadata_field(function,index,reader):
    if reader is None:return _FUNCTION_ATTRIBUTE(function,_FUNCTION_FIELDS[index])
    values=reader(function)
    return None if values is None else values[index]


def _function(function):
    if _TYPE(function) is not FunctionType:return None
    return (function,function.__code__,function.__defaults__,function.__kwdefaults__,
            _TUPLE(cell.cell_contents for cell in (function.__closure__ or ())))


def _unchanged(record,_metadata=None):
    function,code,defaults,keywords,closure=record
    if _metadata is not None:
        current=_metadata(function)
        if (current is None or current[0] is not code or current[2] is not defaults
                or current[3] is not keywords):return False
        if (defaults is not None and _TYPE(defaults) is not _TUPLE
                or keywords is not None and _TYPE(keywords) is not dict):return False
        cells=current[4]
        if cells is None:cells=()
        elif _TYPE(cells) is not _TUPLE:return False
        return _LEN(cells)==_LEN(closure) and _ALL(c.cell_contents is value for c,value in zip(cells,closure))
    if (function.__code__ is not code or function.__defaults__ is not defaults
            or function.__kwdefaults__ is not keywords):return False
    current=function.__closure__ or ()
    return _LEN(current)==_LEN(closure) and _ALL(c.cell_contents is value for c,value in zip(current,closure))


def _lookup_unchanged(module,name,value,_metadata=None):
    current=module.get(name)
    if current is value:return True
    saved=_RELOADABLE.get(value) if _TYPE(value) is FunctionType else None
    if saved is None or _TYPE(current) is not FunctionType:return False
    code,defaults,keywords,closure=saved
    if _metadata is not None:
        values=_metadata(current)
        if (values is None or values[4] is not None and _TYPE(values[4]) is not _TUPLE
                or values[2] is not None and _TYPE(values[2]) is not _TUPLE
                or values[3] is not None and _TYPE(values[3]) is not _DICT):return False
    if (_metadata_field(current,1,_metadata) is not module or _metadata_field(current,2,_metadata) is not defaults
            or _metadata_field(current,3,_metadata) is not keywords):return False
    # A stock module reload creates new function identities. Authenticate its
    # executable code, without calling equality/reduce on user constants.
    try:
        return (_CODE(_metadata_field(current,0,_metadata),0)==code
                and _LEN(_metadata_field(current,4,_metadata) or ())==_LEN(closure)
                and _ALL(cell.cell_contents is old for cell,old in zip(_metadata_field(current,4,_metadata) or (),closure)))
    except (ValueError,TypeError):return False


def _class(cls):
    namespace=_TYPE_DICT.__get__(cls)
    items=_TUPLE(namespace.items())
    functions=[]
    for name,value in items:
        # Exact rational ground arithmetic never presents an expression.
        # Do not retain accidental executable checks for late print wrappers.
        if cls is _core.Expr and name in ('__str__','__repr__'):continue
        function=value.__func__ if _TYPE(value) in (classmethod,staticmethod) else value
        record=_function(function)
        if record is not None:functions.append(record)
    return cls,_TYPE_BASES.__get__(cls),items,_TUPLE(functions)


def _stats_installed(cls,old,current,_metadata=None):
    if _STATS is None or cls is not _STATS[0] or _LEN(_STATS[1])!=2:return False
    names=('doit','rewrite','_stats_doit_installed')
    if (any(name in old for name in names) or _LEN(current)!=_LEN(old)+3
            or current.get('_stats_doit_installed') is not True
            or not _ALL(current.get(name) is value for name,value in old.items())):return False
    module=_MODULES.get('symox.stats._symbolic')
    if _TYPE(module) is not ModuleType:return False
    namespace=module.__dict__;installer=namespace.get('install')
    if (_TYPE(installer) is not FunctionType
            or _metadata_field(installer,0,_metadata) is not namespace.get('_STATS_INSTALL_CODE')):return False
    # These methods are not used by commutative exact ground arithmetic.
    # Nevertheless, only the two official wrappers over the bootstrap
    # delegates are admitted; arbitrary Expr overrides retain the full path.
    for name,freevar,delegate in _STATS[1]:
        code=namespace.get('_STATS_'+name.upper()+'_CODE');function=current.get(name)
        if (_TYPE(code) is not CodeType or not any(value is code for value in _metadata_field(installer,0,_metadata).co_consts)
                or _TYPE(function) is not FunctionType or _metadata_field(function,0,_metadata) is not code
                or _metadata_field(function,1,_metadata) is not namespace or _metadata_field(function,2,_metadata) is not None
                or _metadata_field(function,3,_metadata) is not None or code.co_freevars!=(freevar,)):return False
        closure=_metadata_field(function,4,_metadata)
        if _metadata is not None and closure is not None and _TYPE(closure) is not _TUPLE:return False
        closure=closure or ()
        try:
            if _LEN(closure)!=1 or closure[0].cell_contents is not delegate:return False
        except ValueError:return False
    return True


def _host_free_symbols_installed(cls,old,current,_metadata=None):
    saved=_HOST_FREE_SYMBOLS
    if saved is None or cls is not saved[0]:return False
    _,module,namespace,getter,installer,property_type,functions,lookups,hooks=saved
    if ('free_symbols' in old or _LEN(current)!=_LEN(old)+1
            or not _ALL(current.get(name) is value for name,value in old.items())):return False
    if (_TYPE(module) is not ModuleType or _MODULES.get('symox._host') is not module
            or module.__dict__ is not namespace
            or not _ALL(target.get(name) is value for target,name,value in lookups)
            or not _ALL(_unchanged(function,_metadata) for function in functions)
            or _metadata_field(getter,1,_metadata) is not namespace or _metadata_field(installer,1,_metadata) is not namespace
            or any(name in namespace for name in ('property','set'))
            or hooks.get('free_symbols') is not True):return False
    # A shadowed builtin at package import must not bless a property subtype.
    # Its descriptor lookups can themselves execute Python callbacks.
    if (_TYPE(property_type) is not _TYPE or property_type.__flags__ & (1<<9)
            or property_type.__name__!='property' or property_type.__module__!='builtins'):return False
    value=current.get('free_symbols')
    # This single official late addition is below the canonical scalar classes'
    # own free_symbols descriptor in their MRO. Ground arithmetic cannot consume
    # it. Keep rejecting foreign properties and every other class mutation.
    return (_TYPE(value) is property_type and value.fget is getter
            and value.fset is None and value.fdel is None)


def _class_unchanged(record, classes=(),_metadata=None):
    cls,bases,items,functions=record
    current=_TYPE_DICT.__get__(cls);parents=_TYPE_BASES.__get__(cls)
    relevant=items
    if classes and cls is _core.Expr:
        # Ignoring values cannot let a removed print key hide a new hook.
        if not _ALL(name in current for name,_ in items):return False
        relevant=_TUPLE((name,value) for name,value in items if name not in ('__str__','__repr__'))
    if (_LEN(parents)!=_LEN(bases) or not _ALL(a is b for a,b in zip(parents,bases))
            or not _ALL(_unchanged(function,_metadata) for function in functions)):return False
    if (_LEN(current)==_LEN(items)
            and _ALL(current.get(name) is value for name,value in relevant)):return True
    # copy.copy/pickle of an instance caches __slotnames__ in the class
    # dictionary (copyreg); only pickling reads it. Compare without it.
    stripped=_strip_slot_cache(cls,current,items)
    if stripped is None:return False
    if stripped[0] is not current or stripped[1] is not items:
        current,items=stripped
        relevant=(_TUPLE((name,value) for name,value in items if name not in ('__str__','__repr__'))
                  if classes and cls is _core.Expr else items)
        if (_LEN(current)==_LEN(items)
                and _ALL(current.get(name) is value for name,value in relevant)):return True
    if classes and _host_free_symbols_installed(cls,dict(items),current,_metadata):return True
    if classes and _STATS is not None and cls is _STATS[0] and _stats_installed(cls,dict(items),current,_metadata):return True
    if cls is not _core.Expr or not classes or _FORWARDING is None:return False
    old=dict(items)
    markers=('_coordinate_vector_dispatch_installed','_ore_dispatch_installed','_diffgeom_dispatch_installed')
    added=_TUPLE(name for name in markers if name not in old and name in current)
    if (_LEN(current)!=_LEN(items)+_LEN(added)
            or not _ALL(current[name] is True for name in added)):return False
    arithmetic=('__add__','__sub__','__mul__','__truediv__')
    exceptions=arithmetic
    if '_diffgeom_dispatch_installed' in added:
        if _DIFFGEOM_AUX is None or not _DIFFGEOM_AUX(cls,old,current):return False
        exceptions+=('atoms','rcall')
    if not _ALL(current.get(name) is value or name in exceptions for name,value in relevant):return False
    # Optional vector/Ore imports install certified scalar forwarders. Reuse
    # their existing code/closure/delegate verifier, with the exact bootstrap
    # descriptors and function states. No optional module is imported here.
    verify,operators,states=_FORWARDING
    try:
        if _metadata is None:return verify(None,None,operators,states)
        return _FORWARDING_METADATA(None,None,operators,states,_metadata)
    except (AttributeError,TypeError,ValueError):return False


def register_root(namespace):
    """Called once at the end of package installation, before user hooks."""
    global _ROOT,_FORWARDING,_FORWARDING_METADATA,_STATS,_HOST_FREE_SYMBOLS,_METADATA,_DIFFGEOM_AUX,_RATIONAL_BUILTINS
    import builtins
    from .. import _noncommutative as nc, _node_types as nodes, _host as host
    from ..tensor.array import _adapters as arrays
    from . import _finite_field_options as finite, _inexact_functions as inexact
    from importlib import import_module
    sympification=import_module('symox._sympify')
    core_sympify=import_module('symox.core.sympify')
    decimal=import_module('decimal')
    rationaltools=import_module('symox.integrals.rationaltools')
    if _TYPE(_core) is ModuleType and _plain_namespace(_core.__dict__):
        reader=_core.__dict__.get('_ddm_function_metadata')
        if (_TYPE(reader) is _BUILTIN and reader.__self__ is _core
                and reader.__name__=='_ddm_function_metadata'):_METADATA=(_core,reader)
    from .._basic import Basic
    from .._public_core import Dict
    host_namespace=host.__dict__
    host_getter=host._tree_free_symbols;host_installer=host._install_free_symbols
    host_lookups=((host_namespace,'_tree_free_symbols',host_getter),
                  (host_namespace,'_install_free_symbols',host_installer),
                  (host_namespace,'HostTree',host.HostTree),
                  (host_namespace,'_hooks',host._hooks),
                  (host_namespace,'_core',_core),
                  (_core.__dict__,'_host_contains',_core._host_contains),
                  (builtins.__dict__,'property',property),
                  (builtins.__dict__,'set',set))
    # A module override of either builtin also changes the real function lookup.
    host_lookups+=_TUPLE((host_namespace,name,host_namespace.get(name))
                         for name in ('property','set'))
    _HOST_FREE_SYMBOLS=(host.HostTree,host,host_namespace,host_getter,host_installer,
                       property,(_function(host_getter),_function(host_installer)),
                       host_lookups,host._hooks)
    cancel=namespace['cancel']
    functions=[];lookups=[];seen=set()
    modules=(nc.__dict__,arrays.__dict__,finite.__dict__,inexact.__dict__,nodes.__dict__,sympification.__dict__,core_sympify.__dict__,rationaltools.__dict__)
    def visit(function):
        if _TYPE(function) is not FunctionType or id(function) in seen:return
        seen.add(id(function));functions.append(_function(function))
        if function.__globals__ is rationaltools.__dict__:
            state=functions[-1]
            _RELOADABLE[function]=(_CODE(function.__code__,0),state[2],state[3],state[4])
        for cell in function.__closure__ or ():
            value=cell.cell_contents
            if _TYPE(value) is FunctionType:visit(value)
        for name in function.__code__.co_names:
            if name in function.__globals__:
                value=function.__globals__[name]
                lookups.append((function.__globals__,name,value))
                if _TYPE(value) is FunctionType and any(value.__globals__ is module for module in modules):visit(value)
    visit(cancel)
    visit(rationaltools._atan_stock_operators)
    visit(rationaltools._atan_diffgeom_auxiliary)
    _DIFFGEOM_AUX=rationaltools._atan_diffgeom_auxiliary
    operators=_TUPLE((_core.Expr,name,_TYPE_DICT.__get__(_core.Expr)[name])
                     for name in ('__add__','__sub__','__mul__','__truediv__'))
    _FORWARDING=(rationaltools._atan_stock_operators,operators,
                 rationaltools._atan_function_states(value for _,_,value in operators))
    _FORWARDING_METADATA=rationaltools._atan_stock_operators_with_metadata
    for function in (nc.is_nc,nc._expand_float_logs,finite._rewrite,finite._relabel_roots,inexact.selected,nodes.is_native_node,core_sympify._sympify):visit(function)
    # Every intermediate expression on this exact rational route has one of
    # these native classes. No user subtype/host callback is admitted.
    names=('Add','Mul','Pow','Symbol','Dummy','Integer','Rational','Half','Zero','One','NegativeOne')
    classes=(_core.Expr,)+_TUPLE(nodes.node_class(name) for name in names)
    classes=_TUPLE(dict.fromkeys(classes))
    ancestors=dict.fromkeys(base for cls in classes for base in _TYPE_MRO.__get__(cls))
    ancestors[nodes.NodeType]=None
    for cls in (Basic,Dict,core_sympify.CantSympify):
        for base in _TYPE_MRO.__get__(cls):ancestors[base]=None
        ancestors[_TYPE(cls)]=None
    class_records=_TUPLE(_class(cls) for cls in ancestors)
    stats_class=nodes.node_class('Expr');delegates=[]
    for name,freevar in (('doit','original'),('rewrite','original_rewrite')):
        for base in _TYPE_MRO.__get__(stats_class):
            raw=_TYPE_DICT.__get__(base)
            if name in raw:
                delegates.append((name,freevar,raw[name]));break
    _STATS=(stats_class,_TUPLE(delegates))
    aliases=_TUPLE((namespace,name,namespace[name]) for name in ('cancel','log','Float','Expr','sympify'))
    aliases+=_TUPLE((_core.__dict__,name,_core.__dict__[name]) for name in ('cancel','_has_noncommutative_atoms','Expr','Matrix'))
    aliases+=((decimal.__dict__,'Decimal',decimal.Decimal),)
    builtin_names=('__import__','isinstance','type','all','any','len','tuple','list','set','sorted','enumerate','range','zip','getattr','hasattr','sum','dict','object','id','next','str','TypeError','ValueError','IndexError')
    builtins_snapshot=_TUPLE((builtins.__dict__,name,builtins.__dict__[name]) for name in builtin_names)
    rational_namespace=rationaltools.__dict__
    _RATIONAL_BUILTINS=(rational_namespace,_TUPLE((name,name in rational_namespace,rational_namespace.get(name))for name in builtin_names))
    for name,present,value in _RATIONAL_BUILTINS[1]:
        if present and value is not builtins.__dict__[name]:_METADATA=None
    _ROOT=(cancel,_TUPLE(functions),_TUPLE(lookups)+aliases+builtins_snapshot,class_records,classes,
           builtin_names,builtins.__dict__,sympification.__dict__,sympification.converter,
           _TUPLE(sympification.converter.items()))


def register(cls,namespace):
    record=_class(cls)
    attributes=_TYPE_DICT.__get__(cls)['__dict__']
    functions=_TUPLE(value for name in ('_cancel','_sympify') if (value:=_function(namespace[name])) is not None)
    aliases=_TUPLE((name,namespace[name]) for name in ('_cancel','_sympify','Expr','Poly'))
    original_builtins=namespace.get('__builtins__')
    def reduce(matrix,basis=False):
        root=_ROOT
        if root is None or _TYPE(matrix) is not cls:return None
        cancel,root_functions,lookups,class_records,classes,builtin_names,builtins_dictionary,conversion_namespace,converters,converter_entries=root
        metadata=_metadata_reader()
        # An unavailable/foreign reader requests the entire legacy body.
        if metadata is None:return None
        scopes=[]
        for scope,_,_ in lookups:
            found=False
            for previous in scopes:
                if scope is previous:found=True;break
            if not found:scopes.append(scope)
        scopes.append(namespace)
        for scope in scopes:
            if not _plain_namespace(scope):return None
        rational_namespace,rational_builtins=_RATIONAL_BUILTINS
        if not _plain_namespace(rational_namespace):return None
        for name,present,value in rational_builtins:
            if (name in rational_namespace) is not present or rational_namespace.get(name) is not value:return None
        for scope,name,value in lookups:
            if scope is builtins_dictionary and scope.get(name) is not value:return None
        if not _ALL(_lookup_unchanged(module,name,value,metadata) for module,name,value in lookups):return None
        if (conversion_namespace.get('converter') is not converters
                or _LEN(converters)!=_LEN(converter_entries)
                or not _ALL(key is old_key and value is old_value for (key,value),(old_key,old_value) in zip(converters.items(),converter_entries))
                or not _ALL(_lookup_unchanged(module,name,value,metadata) for module,name,value in lookups)
                or not _ALL(namespace.get(name,builtins_dictionary[name]) is builtins_dictionary[name] for name in builtin_names)
                or not _class_unchanged(record,_metadata=metadata)
                or namespace.get('__builtins__') is not original_builtins
                or not _ALL(namespace.get(name) is value for name,value in aliases)
                or not _ALL(_unchanged(function,metadata) for function in functions)
                or not _ALL(_unchanged(function,metadata) for function in root_functions)
                or not _ALL(_class_unchanged(record,classes,metadata) for record in class_records)):
            return None
        attrs=attributes.__get__(matrix,cls)
        if (_TYPE(attrs) is not dict or not _ALL(_TYPE(key) is str for key in attrs)
                or not _ALL(name in attrs for name in ('_rows','_cols','_items','gens','domain'))
                or _ALL(name not in attrs for name,_ in record[2]) is False):return None
        rows,cols,items=attrs['_rows'],attrs['_cols'],attrs['_items']
        if (_TYPE(rows) is not int or _TYPE(cols) is not int or rows<0 or cols<0
                or _TYPE(attrs['gens']) is not tuple or _TYPE(items) is not _LIST):return None
        if _TYPE_DICT.__get__(classes[0]).get('_diffgeom_dispatch_installed') is True:
            # The bare extension type can report a user-replaced public facade
            # through __class__. Nominal stock nodes have the proved MRO.
            if any(_TYPE(value) is classes[0] for value in items):return None
        return _core._polymatrix_ground(items,rows,cols,attrs['domain'],basis,classes)
    return reduce
