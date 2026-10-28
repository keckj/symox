"""Protocol admission only for the exact native polynomial_reduce loop."""
import builtins as _builtins
from types import FunctionType,ModuleType,BuiltinFunctionType
import marshal as _marshal
from marshal import dumps as _CODE
from .. import _core
from ..polys import _polymatrix_native as _ground

_TYPE=type
_ALL=all
_ANY=any
_LEN=len
_ZIP=zip
_TUPLE=tuple
_DICT=dict
_ROOT=None


def register_root(namespace):
    """Capture final stock constructors/functions before external hooks."""
    global _ROOT
    if _ROOT is not None:return
    # Unlike a name-only bootstrap capture, foreign wrappers cannot bless
    # themselves as canonical helpers when installed before package import.
    if _TYPE is not _core.Poly.__class__:return
    if (_TYPE(_marshal) is not ModuleType or _TYPE(_CODE) is not BuiltinFunctionType
            or _CODE.__self__ is not _marshal or _CODE.__name__!='dumps'
            or _marshal.__dict__.get('__name__')!='marshal'
            or _marshal.__dict__.get('dumps') is not _CODE):return
    serialization=(_marshal,_CODE,globals())
    for name in ('all','any','len','getattr','isinstance','sum','__import__'):
        value=_builtins.__dict__[name]
        if (_TYPE(value) is not BuiltinFunctionType or value.__self__ is not _builtins
                or value.__name__!=name):return
    for name in ('type','list','tuple','dict','zip'):
        value=_builtins.__dict__[name]
        if (_TYPE(value) is not _TYPE or value.__flags__&512
                or value.__module__!='builtins' or value.__name__!=name):return
    from importlib import import_module
    risch=import_module('symox.integrals.risch')
    from .. import _poly_constructor as constructor,_node_types as nodes
    from ..polys import _poly_methods as methods,_evaluation as evaluation
    from . import rationaltools
    cls=risch.DifferentialExtension
    slots=tuple((name,cls.__dict__[name])for name in ('t','d','T','D','level'))
    preflight=_core.__dict__.get('_risch_decline_no_step')
    if (_TYPE(preflight) is not BuiltinFunctionType or preflight.__self__ is not _core
            or preflight.__name__!='_risch_decline_no_step'):return
    preflight=(_core,preflight,slots[1][1])
    # Slot reads below use the original native member descriptors directly.
    classes=(_ground._class(_core.Poly),_ground._class(cls),
             _ground._class(nodes.node_class('NegativeInfinity')))
    functions=[];lookups=[];mutable=[];seen=set()
    modules=(risch.__dict__,constructor.__dict__,methods.__dict__,evaluation.__dict__)
    def visit(function):
        if _TYPE(function) is not FunctionType or id(function) in seen:return
        seen.add(id(function));functions.append(_ground._function(function))
        for cell in function.__closure__ or ():
            value=cell.cell_contents
            if _TYPE(value) is FunctionType:visit(value)
        for name in function.__code__.co_names:
            ns=function.__globals__
            if name not in ns:continue
            value=ns[name];lookups.append((ns,name,value))
            if _TYPE(value) is dict:
                mutable.append((value,tuple(value.items())))
            if _TYPE(value) is FunctionType and _ANY(value.__globals__ is m for m in modules):visit(value)
    visit(risch.derivation)
    for name in ('__new__','as_poly','as_expr','mul','degree','diff','LC','__add__','__sub__','__iadd__'):
        visit(getattr(_core.Poly,name,None))
    # Constructor paths and scalar wrappers have an existing separately
    # authenticated bootstrap proof. Also authenticate that verifier itself.
    verify=rationaltools._stock_ratint_protocol
    verify_state=_ground._function(verify)
    lookup_classes=tuple((_TYPE(m),m,m.__dict__)for m in (risch,constructor,methods,evaluation,rationaltools))
    aliases=tuple((_builtins.__dict__,name,_builtins.__dict__[name])for name in
                  ('type','all','any','len','getattr','isinstance','sum','zip','tuple','dict','list','__import__'))
    aliases+=tuple((risch.__dict__,name,risch.__dict__.get(name,_builtins.__dict__.get(name)))for name in ('Poly','len','zip','derivation'))
    _ROOT=(cls,slots,classes,tuple(functions),tuple(lookups),tuple(mutable),
           verify,verify_state,rationaltools,lookup_classes,aliases,_ground._ROOT,serialization,preflight)


def reduce(p,de,initial):
    root=_ROOT
    if root is None:return None
    cls,slots,classes,functions,lookups,mutable,verify,verify_state,rationaltools,modules,aliases,ground,serialization,preflight=root
    core_owner,decline,d_slot=preflight
    if (_core is not core_owner or _TYPE(core_owner) is not ModuleType
            or core_owner.__dict__.get('_risch_decline_no_step') is not decline):return None
    if decline(p,de,d_slot):return None
    owner,serializer,selector_namespace=serialization
    if (_TYPE(_ground) is not ModuleType
            or _ground.__dict__.get('_CODE') is not serializer
            or selector_namespace.get('_CODE') is not serializer
            or selector_namespace.get('_marshal') is not owner
            or _TYPE(owner) is not ModuleType
            or owner.__dict__.get('dumps') is not serializer):return None
    if (_TYPE(_core) is not ModuleType or _TYPE(_ground) is not ModuleType
            or _TYPE(de) is not cls or _TYPE(p) is not _core.Poly):return None
    if _ANY(_TYPE(m) is not kind or m.__dict__ is not ns for kind,m,ns in modules):return None
    if _ANY(ns.get(name,_builtins.__dict__.get(name)) is not value for ns,name,value in aliases):return None
    if _ANY(ns.get(name) is not value for ns,name,value in lookups):return None
    for value,items in mutable:
        if (_LEN(value)!=_LEN(items) or any(k is not old or v is not prior
                for (k,v),(old,prior) in _ZIP(value.items(),items))):return None
    try:
        if (not _ALL(_ground._class_unchanged(c)for c in classes)
                or not _ALL(_ground._unchanged(f)for f in functions)
                or not _ground._unchanged(verify_state)):return None
        current=rationaltools.__dict__.get('_stock_ratint_protocol')
        if current is not verify:
            # A stock reload creates a new function identity. Compare its full
            # code graph (including nested constants), real globals, defaults
            # and closure, then keep calling the authenticated bootstrap
            # delegate. Never capture the replacement as a new trust anchor.
            if (_TYPE(current) is not FunctionType
                    or current.__globals__ is not rationaltools.__dict__
                    or current.__defaults__ is not verify.__defaults__
                    or current.__kwdefaults__ is not verify.__kwdefaults__
                    or current.__closure__ is not verify.__closure__):return None
            codes=[]
            for code in (current.__code__,verify_state[1]):
                # Check before each use, then call the kept C delegate rather
                # than a restorable module/global alias.
                if (selector_namespace.get('_CODE') is not serializer
                        or selector_namespace.get('_marshal') is not owner
                        or _TYPE(owner) is not ModuleType
                        or owner.__dict__.get('dumps') is not serializer):return None
                codes.append(serializer(code,0))
            if codes[0]!=codes[1]:return None
        if _ground._ROOT is not ground or ground is None:return None
        _,root_functions,root_lookups,records,scalar_classes,_,_,_,_,_=ground
        if (not _ALL(_ground._lookup_unchanged(ns,n,v)for ns,n,v in root_lookups)
                or not _ALL(_ground._unchanged(f)for f in root_functions)
                or not _ALL(_ground._class_unchanged(c,scalar_classes)for c in records)
                or not verify(None)):return None
    except (AttributeError,TypeError,ValueError):return None
    # The exact class and all stored member descriptors were proved above.
    # No mutable attribute lookup or foreign instance hook is invoked here.
    try:fields={name:slot.__get__(de,cls)for name,slot in slots}
    except AttributeError:return None
    return _core._risch_polynomial_reduce(p,fields['d'],fields['t'],fields['T'],
                                        fields['D'],fields['level'],initial)


def _prove_fields(p,de,metadata,verify_with_metadata,lookup):
    """Shared stock-field proof for splitfactor, using a required C reader.

    The original polynomial_reduce function above is deliberately unchanged.
    Its bootstrap state remains usable by existing reload/error contracts.
    The caller authenticates this function, all namespaces and the reader
    before entry; a changed/reloaded verifier requests the full old body.
    """
    root=_ROOT
    if root is None:return None
    cls,slots,classes,functions,lookups,mutable,verify,verify_state,rationaltools,modules,aliases,ground,serialization,preflight=root
    owner,serializer,selector_namespace=serialization
    if (_TYPE(owner) is not ModuleType
            or owner.__dict__.get('dumps') is not serializer
            or selector_namespace.get('_CODE') is not serializer
            or selector_namespace.get('_marshal') is not owner
            or _ground.__dict__.get('_CODE') is not serializer):return None
    if (_TYPE(de) is not cls or _TYPE(p) is not _core.Poly
            or _ANY(_TYPE(m) is not kind or m.__dict__ is not ns for kind,m,ns in modules)
            or _ANY(ns.get(name,_builtins.__dict__.get(name)) is not value for ns,name,value in aliases)
            or _ANY(ns.get(name) is not value for ns,name,value in lookups)):return None
    for value,items in mutable:
        if (_TYPE(value) is not _DICT or _LEN(value)!=_LEN(items)
                or _ANY(k is not old or v is not prior
                    for (k,v),(old,prior)in _ZIP(value.items(),items))):return None
    try:
        if (not _ALL(_ground._class_unchanged(c,_metadata=metadata)for c in classes)
                or not _ALL(_ground._unchanged(f,metadata)for f in functions)
                or not _ground._unchanged(verify_state,metadata)
                or not lookup(rationaltools.__dict__,'_stock_ratint_protocol',verify)
                or _ground.__dict__.get('_ROOT') is not ground or ground is None):return None
        _,root_functions,root_lookups,records,scalar_classes,_,_,_,_,_=ground
        if (not _ALL(lookup(ns,n,v)for ns,n,v in root_lookups)
                or not _ALL(_ground._unchanged(f,metadata)for f in root_functions)
                or not _ALL(_ground._class_unchanged(c,scalar_classes,metadata)for c in records)
                or not verify_with_metadata(None,metadata)):return None
        fields={name:slot.__get__(de,cls)for name,slot in slots}
    except (AttributeError,TypeError,ValueError):return None
    return fields
