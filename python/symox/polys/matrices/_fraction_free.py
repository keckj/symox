"""Stock domain/operator protocol guards for the native FF reductions."""
from types import CodeType,FunctionType,WrapperDescriptorType,MethodDescriptorType,GetSetDescriptorType
from symox import _core
from ..domains import _methods,_algebraic_methods

_TYPE=type
_ALL=all
_DOMAIN=_core.Domain
_ELEMENT=_core.DomainElement
_DOMAIN_HOOKS=tuple((name,_DOMAIN.__dict__.get(name)) for name in
                    ('one','zero','is_one','quo','is_Composite'))
_ELEMENT_HOOKS=tuple((name,_ELEMENT.__dict__.get(name)) for name in
                     ('__bool__','__mul__','__sub__','__add__','__neg__'))
_EXQUO_CODE=globals().get('_EXQUO_CODE',_algebraic_methods._exquo.__code__)
_EXACT_CODE=globals().get('_EXACT_CODE',next((c for c in _methods.install.__code__.co_consts
                       if _TYPE(c) is CodeType and c.co_freevars==('native_exact',)),None))
_NATIVE=(_ALL(_TYPE(method) in (MethodDescriptorType,GetSetDescriptorType)
              and method.__objclass__ is _DOMAIN and method.__name__==name
              for name,method in _DOMAIN_HOOKS)
         and _ALL(_TYPE(method) is WrapperDescriptorType
                  and method.__objclass__ is _ELEMENT and method.__name__==name
                  for name,method in _ELEMENT_HOOKS))


def _native_domain_method(method,name):
    return (_TYPE(method) in (MethodDescriptorType,GetSetDescriptorType)
            and method.__objclass__ is _DOMAIN and method.__name__==name)


def _stock_exact(descriptor):
    if _native_domain_method(descriptor,'is_Exact'):return True
    if _TYPE(descriptor) is not property:return False
    function=descriptor.fget
    if _TYPE(function) is not FunctionType or function.__code__ is not _EXACT_CODE:return False
    try:native,=function.__closure__
    except (TypeError,ValueError):return False
    try:native=native.cell_contents
    except ValueError:return False
    return _native_domain_method(native,'is_Exact')


def _stock_exquo(method):
    if _native_domain_method(method,'exquo'):return True
    return (_TYPE(method) is FunctionType and method.__code__ is _EXQUO_CODE
            and method.__globals__ is _algebraic_methods.__dict__
            and _native_domain_method(_algebraic_methods._native_exquo,'exquo'))


def stock_domain(domain):
    # The installed exactness property and exquo formatter delegate unchanged
    # on QQ/ZZ. Their real code/delegates must remain stock; names or copied
    # functools metadata alone cannot authorize a user callback.
    return (_TYPE(domain) is _DOMAIN and _NATIVE
            and _ALL(_DOMAIN.__dict__.get(name) is method for name,method in _DOMAIN_HOOKS)
            and _stock_exact(_DOMAIN.__dict__.get('is_Exact'))
            and _stock_exquo(_DOMAIN.__dict__.get('exquo'))
            and _ALL(_ELEMENT.__dict__.get(name) is method for name,method in _ELEMENT_HOOKS)
            and '__getattribute__' not in _DOMAIN.__dict__
            and '__getattr__' not in _DOMAIN.__dict__
            and '__imul__' not in _ELEMENT.__dict__)

# QQ denominator preparation has a separate, conservative protocol boundary.
# Register only after DomainMatrix is complete; domain wrappers are installed
# later, so their canonical functions/delegates are checked when called.
_CLEAR_FUNCTIONS = tuple((name, getattr(_methods, name), getattr(_methods, name).__code__)
                         for name in ('convert', 'convert_from', '_has_assoc', '_not_a_coeff'))
_CLEAR_ASSOC_CODE = next(c for c in _methods.install.__code__.co_consts
                         if _TYPE(c) is CodeType and '_has_assoc' in c.co_names
                         and 'get_ring' in c.co_consts)
_CLEAR_DOMAIN_NAMES = ('denom', 'numer', 'lcm', 'get_ring', 'is_EX', 'is_EXRAW', 'of_type')
_CLEAR_ELEMENT_NAMES = ('__getattribute__', '__getattr__', '__eq__', 'parent')
_CLEAR_DOMAIN_OTHER_HOOKS = tuple((name, _DOMAIN.__dict__.get(name)) for name in ('__eq__', '__ne__'))
_CLEAR_ELEMENT_HOOKS = tuple((name, _ELEMENT.__dict__.get(name)) for name in _CLEAR_ELEMENT_NAMES)
_CLEAR_PROTOCOL = None


def register_clear_protocol(dm, ddm, sdm, namespace):
    global _CLEAR_PROTOCOL
    import builtins
    from . import _dense_poly
    classes = []
    for cls, base, names in (
            (dm, object, ('to_dod', 'to_flat_nz', 'from_dod_like', 'from_flat_nz', 'from_rep', 'diag')),
            (ddm, list, ('to_dod', 'to_flat_nz', 'to_sdm', 'to_dok', 'from_dod', 'from_flat_nz', '__init__')),
            (sdm, dict, ('to_dod', 'to_flat_nz', 'to_dok', 'from_dod', 'from_flat_nz', 'to_ddm', 'diag', 'from_list', '__init__'))):
        names += ('__getattribute__', '__getattr__', '__bool__', '__iter__', '__len__', '__getitem__', 'items', 'values')
        entries = tuple((name, cls.__dict__.get(name)) for name in names)
        codes = []
        for _, value in entries:
            function = value.__func__ if _TYPE(value) is classmethod else value
            if _TYPE(function) is FunctionType:
                codes.append((function, function.__code__))
        classes.append((cls, base, entries, tuple(codes), cls.__dict__['__dict__']))
    funcs = tuple((name, getattr(_dense_poly, name), getattr(_dense_poly, name).__code__)
                  for name in ('dup_clear_denoms', 'dup_convert', 'dup_strip', '_native_elements'))
    modules = (namespace, _dense_poly.__dict__, _methods.__dict__,
               __import__(ddm.__module__, fromlist=['']).__dict__,
               __import__(sdm.__module__, fromlist=['']).__dict__)
    aliases = tuple((module, name, module.get(name)) for module in modules[3:]
                    for name in ('_native_elements', 'SDM', 'DDM'))
    aliases += ((namespace, 'ZZ', namespace['ZZ']),)
    builtin_names = ('zip', 'dict', 'list', 'enumerate', 'range', 'len', 'bool', 'type', 'getattr', 'isinstance', 'all', 'any')
    lookups = tuple((module, name, module.get(name, getattr(builtins, name)))
                    for module in modules for name in builtin_names)
    _CLEAR_PROTOCOL = (tuple(classes), funcs, namespace, namespace['dup_clear_denoms'], lookups, builtins, aliases)


def _stock_clear_domain(domain):
    if not stock_domain(domain):
        return False
    if not _ALL(_native_domain_method(_DOMAIN.__dict__.get(name), name) for name in _CLEAR_DOMAIN_NAMES):
        return False
    if not _ALL(_DOMAIN.__dict__.get(name) is value for name, value in _CLEAR_DOMAIN_OTHER_HOOKS):
        return False
    if not _ALL(_ELEMENT.__dict__.get(name) is value for name, value in _CLEAR_ELEMENT_HOOKS):
        return False
    if '_symox_monogenic_extension' in _DOMAIN.__dict__:
        return False
    if not _ALL(getattr(_methods, name) is function and function.__code__ is code
                for name, function, code in _CLEAR_FUNCTIONS):
        return False
    if (_DOMAIN.__dict__.get('convert') is not _methods.convert
            or _DOMAIN.__dict__.get('convert_from') is not _methods.convert_from
            or not _native_domain_method(_methods._native_convert, 'convert')):
        return False
    assoc = _DOMAIN.__dict__.get('has_assoc_Ring')
    return (_TYPE(assoc) is property and _TYPE(assoc.fget) is FunctionType
            and assoc.fget.__code__ is _CLEAR_ASSOC_CODE
            and assoc.fget.__globals__ is _methods.__dict__)


def native_clear_denoms(matrix, convert, rowwise):
    protocol = _CLEAR_PROTOCOL
    if protocol is None or _TYPE(convert) is not bool:
        return None
    classes, funcs, namespace, clear, lookups, builtins, aliases = protocol
    dm = classes[0][0]
    if _TYPE(matrix) is not dm or namespace.get('dup_clear_denoms') is not clear:
        return None
    from . import _dense_poly
    if not _ALL(getattr(_dense_poly, name) is function and function.__code__ is code
                for name, function, code in funcs):
        return None
    if not _ALL(module.get(name, getattr(builtins, name)) is value for module, name, value in lookups):
        return None
    if not _ALL(module.get(name) is value for module, name, value in aliases):
        return None
    for cls, base, entries, codes, _ in classes:
        if any(name in cls.__dict__ for name in ('rep', 'shape', 'domain')):
            return None
        if (cls.__bases__ != (base,) or not _ALL(cls.__dict__.get(name) is value for name, value in entries)
                or not _ALL(function.__code__ is code for function, code in codes)):
            return None
    attrs = classes[0][4].__get__(matrix, dm)
    if (_TYPE(attrs) is not dict or not _ALL(_TYPE(key) is str for key in attrs)
            or not _ALL(name in attrs for name in ('rep', 'shape', 'domain'))):
        return None
    if _ALL(name not in attrs for name, _ in classes[0][2]) is False:
        return None
    rep, shape, domain = attrs['rep'], attrs['shape'], attrs['domain']
    info = next((info for info in classes[1:] if _TYPE(rep) is info[0]), None)
    if info is None or not _stock_clear_domain(domain):
        return None
    rep_attrs = info[4].__get__(rep, info[0])
    if (_TYPE(rep_attrs) is not dict or not _ALL(_TYPE(key) is str for key in rep_attrs)
            or rep_attrs.get('shape') is not shape or rep_attrs.get('domain') is not domain
            or not _ALL(name not in rep_attrs for name, _ in info[2])):
        return None
    # CAPI container reads and exact QQ scalar checks occur only after every
    # skipped public protocol has been proved stock. Output factories stay in
    # DomainMatrix and retain their established type/zero/format conventions.
    return _core._dm_clear_denoms(rep, shape, domain, convert, rowwise, info[0])
