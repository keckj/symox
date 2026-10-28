"""Stock protocol boundary for the native inexact linear-system reduction."""
from types import WrapperDescriptorType, MethodDescriptorType, FunctionType, CodeType
from .. import _core
from .matrices._fraction_free import _stock_clear_domain
from .domains import _algebraic_methods

_TYPE=type
_ALL=all
_DOMAIN=_core.Domain
_ELEMENT=_core.DomainElement
_ELEMENT_DICT=_ELEMENT.__dict__
_PROBES=('_symox_domainmatrix','_symox_quotient_element','is_Matrix','_symox_sparse','_symox_ring_element')
_OPERATORS=tuple((name,_ELEMENT_DICT.get(name)) for name in ('__bool__','__mul__','__sub__'))
_NATIVE=_ALL(_TYPE(method) is WrapperDescriptorType and method.__objclass__ is _ELEMENT and method.__name__==name
             for name,method in _OPERATORS)

_DIV_CODE=next(code for code in _algebraic_methods._gaussian_binary.__code__.co_consts
               if _TYPE(code) is CodeType and code.co_name=='method')
_DIV_DOMAIN=tuple((name,_DOMAIN.__dict__.get(name)) for name in ('is_ZZ_I','is_QQ_I'))


def _stock_division():
    method=_ELEMENT_DICT.get('__truediv__')
    if _TYPE(method) is WrapperDescriptorType:
        return method.__objclass__ is _ELEMENT and method.__name__=='__truediv__'
    if (_TYPE(method) is not FunctionType or method.__code__ is not _DIV_CODE
            or method.__globals__ is not _algebraic_methods.__dict__
            or not _ALL(_DOMAIN.__dict__.get(name) is value for name,value in _DIV_DOMAIN)):
        return False
    try:
        values=dict(zip(method.__code__.co_freevars,(c.cell_contents for c in method.__closure__)))
    except (ValueError,TypeError):
        return False
    native=values.get('native')
    return (values.get('op')=='truediv' and values.get('reflected') is False
            and _TYPE(native) is WrapperDescriptorType
            and native.__objclass__ is _ELEMENT and native.__name__=='__truediv__')


def register(namespace):
    import builtins
    dictionary=builtins.__dict__
    original_builtins=namespace.get('__builtins__')
    names=('len','range','next','enumerate','zip','any','all','getattr','type','property')
    methods=tuple((name,dictionary[name]) for name in names)
    guard_globals=_stock_clear_domain.__globals__
    domain_to_sympy=_DOMAIN.__dict__.get('to_sympy')
    native_to_sympy=(_TYPE(domain_to_sympy) is MethodDescriptorType
                     and domain_to_sympy.__objclass__ is _DOMAIN and domain_to_sympy.__name__=='to_sympy')

    def reduce(rows,domain,width):
        # All conversion calls and the public answer assembly remain in their
        # original locations. Only the loop whose scalar protocols are proved
        # stock is skipped. No user getter is evaluated speculatively.
        if (not _NATIVE or not _stock_division() or not native_to_sympy or namespace.get('__builtins__') is not original_builtins
                or not _ALL(namespace.get(name,dictionary[name]) is method for name,method in methods)
                or not _ALL(guard_globals.get(name,dictionary[name]) is method for name,method in methods)
                or not _stock_clear_domain(domain)
                or _DOMAIN.__dict__.get('to_sympy') is not domain_to_sympy
                or not _ALL(_ELEMENT_DICT.get(name) is method for name,method in _OPERATORS)
                or '__getattribute__' in _ELEMENT_DICT or '__getattr__' in _ELEMENT_DICT
                or any(name in _ELEMENT_DICT for name in _PROBES)):
            return None
        return _core._inexact_linear_rref(rows,domain,width)
    return reduce
