"""Admit the closed-QQ splitfactor kernel through stock, read-only protocols.

Refusal resumes the complete Python body, after its real kinv/One prefix.
The older polynomial_reduce selector and its executable body remain intact.
"""
from types import ModuleType,CodeType,CellType
from .. import _core
from . import _risch_native as _shared
from ..polys import _polymatrix_native as _ground

_TYPE=(0).__class__.__class__
_FUNCTION=(lambda:None).__class__
_BUILTIN=[].append.__class__
_TUPLE=().__class__
_DICT={}.__class__
_STR=''.__class__
_ROOT=None


_STR_KEYED=_core.__dict__.get('_str_keyed')


def _plain(namespace):
    # Exact dict whose keys are all exact str; the native check reads the keys
    # without hashing, equality or other Python callbacks.
    if _STR_KEYED is not None:return _STR_KEYED(namespace)
    if _TYPE(namespace) is not _DICT:return False
    for key in namespace:
        if _TYPE(key) is not _STR:return False
    return True


def _metadata(function):
    root=_ROOT
    if root is None or _TYPE(function) is not _FUNCTION:return None
    owner,reader=root['reader']
    if (_TYPE(owner) is not ModuleType or _core is not owner
            or not _plain(owner.__dict__)
            or owner.__dict__.get('_ddm_function_metadata') is not reader):return None
    values=reader(function)
    if (values is None
            or values[2] is not None and _TYPE(values[2]) is not _TUPLE
            or values[3] is not None and _TYPE(values[3]) is not _DICT
            or values[4] is not None and _TYPE(values[4]) is not _TUPLE):return None
    return values


def _states_match(states,reader,builtins):
    for function,code,namespace,defaults,keywords,contents,attributes in states:
        if _TYPE(function) is not _FUNCTION:return False
        values=reader(function)
        if (values is None or values[0] is not code or values[1] is not namespace
                or values[2] is not defaults or values[3] is not keywords
                or values[2] is not None and _TYPE(values[2]) is not _TUPLE
                or values[3] is not None and _TYPE(values[3]) is not _DICT
                or function.__builtins__ is not builtins
                or not _plain(namespace)
                or _TYPE(function.__dict__) is not _DICT):return False
        current=function.__dict__
        if len(current)!=len(attributes):return False
        # Identity comparison avoids equality callbacks on arbitrary values.
        for (key,value),(old,prior) in zip(current.items(),attributes):
            if key is not old or value is not prior:return False
        closure=values[4]
        if closure is None:closure=()
        elif _TYPE(closure) is not _TUPLE:return False
        if len(closure)!=len(contents):return False
        try:
            for cell,value in zip(closure,contents):
                if cell.cell_contents is not value:return False
        except ValueError:return False
    return True


def _same_defaults(current,original):
    if current is original:return True
    if _TYPE(current) is not _TUPLE or _TYPE(original) is not _TUPLE:return False
    return len(current)==len(original) and all(a is b for a,b in zip(current,original))


def _same_keywords(current,original):
    if current is original:return True
    if not _plain(current) or not _plain(original) or len(current)!=len(original):return False
    return all(name in current and current[name] is value for name,value in original.items())


def _equivalent(current,original,reader,builtins,namespace,depth=0):
    """Authenticate reload data against the kept bootstrap, without marshal.

    Exact builtin/code containers only; arbitrary constants never get equality,
    hashing or reduction calls. A function must use the real original globals
    and actual stock f_builtins, including after a globals entry is restored.
    """
    if current is original:return True
    if depth>64 or _TYPE(current) is not _TYPE(original):return False
    kind=_TYPE(original)
    if kind is _FUNCTION:
        a,b=reader(current),reader(original)
        if (a is None or b is None or a[1] is not namespace or b[1] is not namespace
                or current.__builtins__ is not builtins or original.__builtins__ is not builtins
                or not _same_defaults(a[2],b[2])
                or not _same_keywords(a[3],b[3])
                or _TYPE(current.__dict__) is not _DICT or _TYPE(original.__dict__) is not _DICT):return False
        return (_equivalent(a[0],b[0],reader,builtins,namespace,depth+1)
                and _equivalent(a[4],b[4],reader,builtins,namespace,depth+1)
                and _equivalent(current.__dict__,original.__dict__,reader,builtins,namespace,depth+1))
    if kind is CodeType:
        a=(current.co_argcount,current.co_posonlyargcount,current.co_kwonlyargcount,
           current.co_nlocals,current.co_stacksize,current.co_flags,current.co_code,
           current.co_consts,current.co_names,current.co_varnames,current.co_filename,
           current.co_name,current.co_firstlineno,current.co_lnotab,
           current.co_linetable,current.co_freevars,current.co_cellvars)
        b=(original.co_argcount,original.co_posonlyargcount,original.co_kwonlyargcount,
           original.co_nlocals,original.co_stacksize,original.co_flags,original.co_code,
           original.co_consts,original.co_names,original.co_varnames,original.co_filename,
           original.co_name,original.co_firstlineno,original.co_lnotab,
           original.co_linetable,original.co_freevars,original.co_cellvars)
        return _equivalent(a,b,reader,builtins,namespace,depth+1)
    if kind is _TUPLE:
        return len(current)==len(original) and all(_equivalent(a,b,reader,builtins,namespace,depth+1)for a,b in zip(current,original))
    if kind is _DICT:
        if not _plain(current) or not _plain(original) or len(current)!=len(original):return False
        return all(name in current and _equivalent(current[name],value,reader,builtins,namespace,depth+1)for name,value in original.items())
    if kind is CellType:
        try:return _equivalent(current.cell_contents,original.cell_contents,reader,builtins,namespace,depth+1)
        except ValueError:return False
    if kind is frozenset:
        # The compiled membership constant used by the vector proof contains
        # only exact Unicode names. Validate both sides before set equality.
        if not all(_TYPE(name) is _STR for name in current):return False
        if not all(_TYPE(name) is _STR for name in original):return False
        return current==original
    if any(kind is scalar for scalar in (str,bytes,int,bool,float,complex)):return current==original
    # Other values (classes, modules, native delegates and sentinel objects)
    # retain their original identity. Foreign constants simply decline.
    return False


def _lookup(namespace,name,original):
    current=namespace.get(name)
    if current is original:return True
    root=_ROOT
    return (namespace is root['rational_namespace']
            and _equivalent(current,original,_metadata,root['builtins'],namespace))


def register_root(namespace):
    """Keep bootstrap delegates; never bless a later replacement snapshot."""
    global _ROOT
    if _ROOT is not None:return
    if (_TYPE(_core) is not ModuleType or _TYPE(_ground) is not ModuleType
            or _TYPE(_shared) is not ModuleType or not _plain(_core.__dict__)
            or _shared._ROOT is None or _ground._ROOT is None):return
    reader=_core.__dict__.get('_ddm_function_metadata')
    kernel=_core.__dict__.get('_risch_splitfactor')
    if (_TYPE(reader) is not _BUILTIN or reader.__self__ is not _core
            or reader.__name__!='_ddm_function_metadata'
            or _TYPE(kernel) is not _BUILTIN or kernel.__self__ is not _core
            or kernel.__name__!='_risch_splitfactor'):return
    precheck=_core.__dict__.get('_risch_splitfactor_precheck')
    if (_TYPE(precheck) is not _BUILTIN or precheck.__self__ is not _core
            or precheck.__name__!='_risch_splitfactor_precheck'):precheck=None
    from importlib import import_module
    risch=import_module('symox.integrals.risch')
    from .. import _poly_constructor as constructor
    from ..polys import _poly_methods as methods,_evaluation as evaluation
    rationaltools=import_module('symox.integrals.rationaltools')
    builtins=_shared._builtins.__dict__
    globals_builtin=builtins.get('globals')
    if (_TYPE(globals_builtin) is not _BUILTIN
            or globals_builtin.__self__ is not _shared._builtins
            or globals_builtin.__name__!='globals'):return
    for name in ('bytes','int','bool','float','complex','frozenset'):
        value=builtins.get(name)
        if (_TYPE(value) is not _TYPE or value.__flags__&512
                or value.__module__!='builtins' or value.__name__!=name):return
    self_namespace=_plain.__globals__
    data_functions=tuple(rationaltools.__dict__[name]for name in
        ('cancel','ratint_ratpart','ratint_logpart','log_to_real','log_to_atan'))
    modules=[_core,_ground,_shared,risch,constructor,methods,evaluation,rationaltools,_shared._builtins,_shared._marshal]
    allowed=tuple(m.__dict__ for m in modules if _TYPE(m) is ModuleType)
    states=[];lookups=[];seen=set()
    def visit(function):
        if _TYPE(function) is not _FUNCTION or id(function) in seen:return True
        seen.add(id(function));values=reader(function)
        if values is None:return False
        code,scope,defaults,keywords,closure=values
        if (not _plain(scope) or function.__builtins__ is not builtins
                or defaults is not None and _TYPE(defaults) is not _TUPLE
                or keywords is not None and _TYPE(keywords) is not _DICT
                or closure is not None and _TYPE(closure) is not _TUPLE
                or _TYPE(function.__dict__) is not _DICT):return False
        if closure is None:closure=()
        try:contents=tuple(cell.cell_contents for cell in closure)
        except ValueError:return False
        states.append((function,code,scope,defaults,keywords,contents,tuple(function.__dict__.items())))
        for value in contents:
            if _TYPE(value) is _FUNCTION and not visit(value):return False
        # These are code-state data for the stock verifier, never executed by
        # it. Their mathematical globals are outside splitfactor's consumed
        # path; the actual callables are still authenticated as code states.
        if scope is rationaltools.__dict__ and any(function is value for value in data_functions):return True
        for name in code.co_names:
            if scope is self_namespace and name=='_ROOT':continue
            present=name in scope;value=scope.get(name,builtins.get(name))
            lookups.append((scope,name,present,value))
            if (_TYPE(value) is _FUNCTION
                    and any(value.__globals__ is m for m in allowed)
                    and not visit(value)):return False
        return True
    # All new/transitive proof functions are authenticated via the C reader.
    verify=rationaltools._stock_ratint_protocol_with_metadata
    prove=_shared._prove_fields
    for function in (prove,verify,_plain,_metadata,_states_match,_same_defaults,_same_keywords,_equivalent,_lookup,reduce):
        if not visit(function):return
    lookups.append((_shared.__dict__,'_prove_fields',True,prove))
    lookups.append((rationaltools.__dict__,'_stock_ratint_protocol_with_metadata',True,verify))
    for name in ('_unchanged','_lookup_unchanged','_class_unchanged'):
        if not visit(_ground.__dict__[name]):return
    for record in _shared._ROOT[3]+_ground._ROOT[1]:
        if not visit(record[0]):return
    for record in _shared._ROOT[2]+_ground._ROOT[3]:
        for state in record[3]:
            if not visit(state[0]):return
    if not visit(risch.splitfactor):return
    for name in ('get_domain','as_poly','exquo','gcd','to_field','diff'):
        if not visit(getattr(_core.Poly,name,None)):return
    for cls in _ground._ROOT[4]:
        if not visit(getattr(cls,'has',None)):return
    # Check every consumed global namespace before performing keyed lookups.
    scopes=[]
    for scope,_,_,_ in lookups:
        if not any(scope is old for old in scopes):scopes.append(scope)
    for scope,_,_ in _shared._ROOT[4]+_shared._ROOT[10]+_ground._ROOT[2]:
        if not any(scope is old for old in scopes):scopes.append(scope)
    for scope in scopes:
        if not _plain(scope):return
    aliases=tuple((m,m.__dict__)for m in modules)
    builtin_names=_ground._ROOT[5]+('globals','bytes','int','bool','float','complex','frozenset')
    builtin_aliases=tuple((name,builtins[name])for name in builtin_names)
    _ROOT={'reader':(_core,reader),'kernel':kernel,'precheck':precheck,'shared':_shared,
           'shared_root':_shared._ROOT,'ground_root':_ground._ROOT,
           'prove':prove,'verify':verify,'states':tuple(states),
           'lookups':tuple(lookups),'scopes':tuple(scopes),'modules':aliases,
           'builtins':builtins,'builtin_aliases':builtin_aliases,
           'rational_namespace':rationaltools.__dict__,
           'local_aliases':tuple((name,self_namespace[name])for name in
                ('_plain','_metadata','_states_match','_same_defaults','_same_keywords','_equivalent','_lookup'))}


def reduce(p,de,one,coefficientD,z):
    root=_ROOT
    if root is None or coefficientD is not False or z is not None:return None
    owner,reader=root['reader']
    # Exact owner/key checks precede consumed lookups, including hash collisions.
    if _TYPE(owner) is not ModuleType or _core is not owner:return None
    if not _plain(owner.__dict__):return None
    if (owner.__dict__.get('_ddm_function_metadata') is not reader
            or owner.__dict__.get('_risch_splitfactor') is not root['kernel']):return None
    # The kernel's own input conditions on p and One, read natively before the
    # authentication below: when they fail the kernel returns None anyway.
    precheck=root['precheck']
    if precheck is not None and (owner.__dict__.get('_risch_splitfactor_precheck') is not precheck
            or not precheck(p,one)):return None
    for module,scope in root['modules']:
        if _TYPE(module) is not ModuleType or module.__dict__ is not scope:return None
    for scope in root['scopes']:
        if not _plain(scope):return None
    builtins=root['builtins']
    for name,value in root['builtin_aliases']:
        if builtins.get(name) is not value:return None
    local_namespace=_plain.__globals__
    for name,value in root['local_aliases']:
        if local_namespace.get(name) is not value:return None
    if (_shared is not root['shared'] or _shared.__dict__.get('_ROOT') is not root['shared_root']
            or _ground.__dict__.get('_ROOT') is not root['ground_root']):return None
    # Authenticate the helpers before invoking them, without Python metadata reads.
    for function,code,namespace,defaults,keywords,contents,attributes in root['states']:
        values=reader(function)
        if (values is None or values[0] is not code or values[1] is not namespace
                or values[2] is not defaults or values[3] is not keywords
                or values[2] is not None and _TYPE(values[2]) is not _TUPLE
                or values[3] is not None and _TYPE(values[3]) is not _DICT
                or values[4] is not None and _TYPE(values[4]) is not _TUPLE):return None
    if not _states_match(root['states'],reader,builtins):return None
    for scope,name,present,value in root['lookups']:
        if (name in scope) is not present:return None
        current=scope.get(name,builtins.get(name))
        if current is not value and not (scope is root['rational_namespace']
                and _equivalent(current,value,_metadata,builtins,scope)):return None
    fields=root['prove'](p,de,_metadata,root['verify'],_lookup)
    if fields is None:return None
    return root['kernel'](p,fields['d'],fields['t'],fields['T'],fields['D'],fields['level'],one)
