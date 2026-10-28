"""This module implements tools for integrating rational functions. """

from symox.core.function import Lambda
from symox.core.numbers import I
from symox.core.singleton import S
from symox.core.symbol import (Dummy, Symbol, symbols)
from symox.functions.elementary.exponential import log
from symox.functions.elementary.trigonometric import atan
from symox.polys.polyerrors import DomainError
from symox.polys.polyroots import roots
from symox.polys.polytools import cancel
from symox.polys.rootoftools import RootSum
from symox.polys import Poly, resultant, ZZ
from symox.polys import _wrt_option as _lrt_wrt
import importlib as _real_importlib
_real_radsimp_module = _real_importlib.import_module('symox.simplify.radsimp')
_real_accum_module = _real_importlib.import_module('symox.calculus._accum_eval')
_real_subs_module = _real_importlib.import_module('symox._function_subs')
_real_compat_module = _real_importlib.import_module('symox._compat')
import builtins as _real_builtins


def _atan_function_states(functions):
    from types import FunctionType
    pending=list(functions);seen=set();states=[]
    while pending:
        function=pending.pop()
        if type(function) is not FunctionType or id(function) in seen:continue
        seen.add(id(function));closure=function.__closure__
        contents=tuple(cell.cell_contents for cell in closure or ())
        states.append((function,function.__code__,closure,contents,function.__defaults__,function.__kwdefaults__))
        pending.extend(value for value in contents if type(value) is FunctionType)
    return tuple(states)


# These guards only select a native route. The full general algorithm below
# remains responsible for custom constructors, Poly methods and held evaluation.
from symox import _core as _rational_core
from symox import _poly_constructor as _lrt_constructor
import symox as _lrt_public
import types as _lrt_types
_ATAN_MODULE_TYPE=_lrt_types.ModuleType
from symox.core.sympify import sympify as _lrt_stock_sympify
from symox import _node_types as _atan_nodes
from symox._public_core import _evaluation as _atan_evaluation
_ATAN_HEAD = globals().get('_ATAN_HEAD', atan)
_ATAN_POLY = globals().get('_ATAN_POLY', Poly)
_ATAN_TYPE = globals().get('_ATAN_TYPE', type)
_ATAN_GETATTR = globals().get('_ATAN_GETATTR', getattr)
_ATAN_FACTORY = globals().get('_ATAN_FACTORY', atan.__dict__.get('_symox_factory'))
_ATAN_CALL = globals().get('_ATAN_CALL', type(atan).__dict__.get('__call__'))
_ATAN_AS_EXPR = globals().get('_ATAN_AS_EXPR', Poly.as_expr)
_ATAN_APPLY = globals().get('_ATAN_APPLY', _rational_core._apply)
_ATAN_OPERATORS = globals().get('_ATAN_OPERATORS', tuple(
    (cls, name, getattr(cls, name))
    for cls in (_rational_core.Expr, atan, _atan_nodes.node_class('atan'), _atan_nodes.node_class('Add'), _atan_nodes.node_class('Mul'), type(S.Zero), type(S.Pi))
    for name in ('__add__', '__radd__', '__mul__', '__rmul__')))
_LRT_RESULTANT = globals().get('_LRT_RESULTANT', _rational_core.resultant)
if '_LRT_WRT_CODE' not in globals():
    _lrt_wrt_factory = _lrt_wrt.__dict__.get('_wrap')
    _LRT_WRT_CODE = next((code for code in (
        _lrt_wrt_factory.__code__.co_consts if type(_lrt_wrt_factory) is _lrt_types.FunctionType else ())
        if type(code) is _lrt_types.CodeType and code.co_name == 'wrapper'), None)
    del _lrt_wrt_factory
_LRT_WRT_GLOBALS = globals().get('_LRT_WRT_GLOBALS', _lrt_wrt.__dict__)
_LRT_CONSTRUCT = _lrt_constructor._native_construct
_LRT_SYMPIFY = globals().get('_LRT_SYMPIFY', _lrt_stock_sympify)
_LRT_AS_POLY_CODE = next(code for code in _lrt_public._install_poly_methods.__code__.co_consts
                         if isinstance(code, _lrt_types.CodeType) and code.co_name == 'as_poly')
_ATAN_ALL = globals().get('_ATAN_ALL', all)
_ATAN_BOOTSTRAP_STATES = globals().get('_ATAN_BOOTSTRAP_STATES',
    _atan_function_states(original for _,_,original in _ATAN_OPERATORS))
# Symox imports this adapter during bootstrap, before later user overrides.
# A removed/reimported module cannot bless arbitrary newly captured operators.
import sys as _atan_sys
_ATAN_BOOTSTRAP_CAPTURE = globals().get('_ATAN_BOOTSTRAP_CAPTURE',
    bool(getattr(getattr(_atan_sys.modules.get('symox'),'__spec__',None),'_initializing',False)))
_ATAN_GENS = globals().get('_ATAN_GENS', Poly.gens)
_ATAN_SCALARS = globals().get('_ATAN_SCALARS', tuple(_atan_nodes.node_class(name) for name in
    ('Symbol','Dummy','Integer','Rational','Zero','One','NegativeOne','Half','Pi','Add','Mul','Pow','atan')))
_ATAN_OBJECT_GETATTRIBUTE = globals().get('_ATAN_OBJECT_GETATTRIBUTE', object.__getattribute__)
_ATAN_MRO_DESCRIPTOR = globals().get('_ATAN_MRO_DESCRIPTOR', type.__dict__['__mro__'])
_ATAN_DICT_DESCRIPTOR = globals().get('_ATAN_DICT_DESCRIPTOR', type.__dict__['__dict__'])
_ATAN_MISSING = globals().get('_ATAN_MISSING', object())


def _atan_raw_attribute(cls,name):
    # Inspect the stored definition; neither instance nor metaclass
    # descriptors may run during route selection.
    for base in _ATAN_MRO_DESCRIPTOR.__get__(cls):
        namespace=_ATAN_DICT_DESCRIPTOR.__get__(base)
        if name in namespace:return namespace[name]
    return _ATAN_MISSING


_ATAN_PRIORITIES = globals().get('_ATAN_PRIORITIES', tuple(
    (cls,_atan_raw_attribute(cls,'_op_priority')) for cls in _ATAN_SCALARS))
_LRT_OPERATORS = globals().get('_LRT_OPERATORS', _ATAN_OPERATORS + tuple(
    (cls,name,getattr(cls,name)) for cls in (_rational_core.Expr,)+_ATAN_SCALARS
    for name in ('__add__','__radd__','__sub__','__rsub__','__mul__','__rmul__',
                 '__truediv__','__rtruediv__','__pow__','__rpow__','__neg__')
    if all(cls is not previous or name!=operation for previous,operation,_ in _ATAN_OPERATORS)))
_LRT_BOOTSTRAP_STATES = globals().get('_LRT_BOOTSTRAP_STATES',
    _atan_function_states(original for _,_,original in _LRT_OPERATORS))




_ATAN_NATIVE_STATES_MATCH = globals().get(
    '_ATAN_NATIVE_STATES_MATCH', _rational_core._atan_states_match)


def _atan_states_match_with_metadata(states,_metadata):
    from types import FunctionType
    if _metadata is None and _rational_core._atan_states_match is _ATAN_NATIVE_STATES_MATCH:
        matched = _ATAN_NATIVE_STATES_MATCH(
            states, FunctionType, (_ATAN_TYPE, _ATAN_ALL, len, zip, tuple))
        if matched is not None:return matched
    if _ATAN_TYPE(states) is not tuple:return False
    try:
        for state in states:
            if _ATAN_TYPE(state) is not tuple or len(state)!=6:return False
            function,code,closure,contents,defaults,kwdefaults=state
            _function_values=_metadata(function) if _metadata is not None else None
            if _metadata is not None and _function_values is None:return False
            if _metadata is not None and (closure is not None and _DIFFGEOM_TYPE(closure) is not _DIFFGEOM_TUPLE
                    or _DIFFGEOM_TYPE(contents) is not _DIFFGEOM_TUPLE
                    or _function_values[4] is not None and _DIFFGEOM_TYPE(_function_values[4]) is not _DIFFGEOM_TUPLE):return False
            if (_ATAN_TYPE(function) is not FunctionType or _ATAN_TYPE(contents) is not tuple
                    or (_function_values[0] if _metadata is not None else function.__code__) is not code or (_function_values[4] if _metadata is not None else function.__closure__) is not closure
                    or (_function_values[2] if _metadata is not None else function.__defaults__) is not defaults or (_function_values[3] if _metadata is not None else function.__kwdefaults__) is not kwdefaults
                    or len(closure or ())!=len(contents)
                    or not _ATAN_ALL(cell.cell_contents is value for cell,value in zip(closure or (),contents))):return False
    except (AttributeError,ValueError,TypeError):return False
    return True


def _atan_forwarding_proof_with_metadata(proof,operation_code,name,_metadata):
    import builtins
    from types import FunctionType
    operation,delegate,states,namespace,ore,mro,isinstance_,getattr_=proof
    _operation_values=_metadata(operation) if _metadata is not None else None
    if _metadata is not None and _operation_values is None:return False
    if _metadata is not None and (_operation_values[4] is not None and _DIFFGEOM_TYPE(_operation_values[4]) is not _DIFFGEOM_TUPLE
            or not _atan_diffgeom_namespace(namespace)):return False
    if (_ATAN_TYPE(operation) is not FunctionType or (_operation_values[0] if _metadata is not None else operation.__code__) is not operation_code
            or (_operation_values[1] if _metadata is not None else operation.__globals__) is not namespace or _ATAN_TYPE(namespace) is not dict
            or namespace.get('_OreOperator') is not ore
            or _ATAN_TYPE(ore) is not _ATAN_TYPE or ore.__mro__ is not mro
            or namespace.get('isinstance',builtins.isinstance) is not isinstance_
            or namespace.get('getattr',builtins.getattr) is not getattr_):return False
    try:
        cells=dict(zip(operation_code.co_freevars,(_operation_values[4] if _metadata is not None else operation.__closure__) or ()))
        if (set(cells)!= {'original','name'} or cells['original'].cell_contents is not delegate
                or _ATAN_TYPE(cells['name'].cell_contents) is not str
                or cells['name'].cell_contents!=name):return False
    except (ValueError,TypeError):return False
    return (_atan_states_match(states) if _metadata is None else _atan_states_match_with_metadata(states,_metadata))


def _atan_vector_proof_with_metadata(proof,operation_code,name,vector_types,_metadata):
    import builtins
    from types import FunctionType
    operation,delegate,states,namespace,basis,mro,isinstance_,operation_name=proof
    _operation_values=_metadata(operation) if _metadata is not None else None
    if _metadata is not None and _operation_values is None:return False
    if _metadata is not None and (_operation_values[4] is not None and _DIFFGEOM_TYPE(_operation_values[4]) is not _DIFFGEOM_TUPLE
            or _DIFFGEOM_TYPE(operation_name) is not _DIFFGEOM_STR
            or not _atan_diffgeom_namespace(namespace)):return False
    if (_ATAN_TYPE(operation) is not FunctionType or (_operation_values[0] if _metadata is not None else operation.__code__) is not operation_code
            or (_operation_values[1] if _metadata is not None else operation.__globals__) is not namespace or _ATAN_TYPE(namespace) is not dict
            or _ATAN_TYPE(basis) is not _ATAN_TYPE or basis.__mro__ is not mro
            or namespace.get('isinstance',builtins.isinstance) is not isinstance_
            or operation_name!=name or not (_atan_states_match(states) if _metadata is None else _atan_states_match_with_metadata(states,_metadata))):return False
    try:
        cells=dict(zip(operation_code.co_freevars,(_operation_values[4] if _metadata is not None else operation.__closure__) or ()))
        if (set(cells)!={'original','name','_BasisDependent'}
                or cells['original'].cell_contents is not delegate
                or cells['_BasisDependent'].cell_contents is not basis
                or _ATAN_TYPE(cells['name'].cell_contents) is not str
                or cells['name'].cell_contents!=name):return False
    except (ValueError,TypeError):return False
    # Exclude these basis classes once for every scalar MRO below, after
    # checking all operator layers. No state survives this selector call.
    vector_types[id(basis)]=basis
    return True


_DIFFGEOM_TYPE=False.__class__.__class__
_DIFFGEOM_STR=''.__class__
_DIFFGEOM_TUPLE=().__class__
_DIFFGEOM_DICT={}.__class__
_DIFFGEOM_FUNCTION=(lambda:None).__class__
_DIFFGEOM_BUILTIN=[].append.__class__
_DIFFGEOM_MODULE_TYPE=_DIFFGEOM_TYPE(_real_builtins)
_DIFFGEOM_SET={0}.__class__


def _atan_diffgeom_namespace(namespace):
    if _DIFFGEOM_TYPE(namespace) is not _DIFFGEOM_DICT:return False
    return _ATAN_ALL(_DIFFGEOM_TYPE(key) is _DIFFGEOM_STR for key in namespace)


def _atan_diffgeom_data():
    # No optional import or new lookup when geometry has not installed a layer.
    raw=_ATAN_DICT_DESCRIPTOR.__get__(_ATAN_OPERATORS[0][0])
    if raw.get('_diffgeom_dispatch_installed') is not True:return None
    import sys
    if not _atan_diffgeom_namespace(sys.modules):return False
    module=sys.modules.get('symox.diffgeom._adapters')
    if (_DIFFGEOM_TYPE(_DIFFGEOM_MODULE_TYPE) is not _DIFFGEOM_TYPE
            or _DIFFGEOM_MODULE_TYPE.__flags__ & (1<<9)
            or _DIFFGEOM_MODULE_TYPE.__name__!='module'
            or _DIFFGEOM_MODULE_TYPE.__module__!='builtins'
            or _DIFFGEOM_TYPE(module) is not _DIFFGEOM_MODULE_TYPE):return False
    namespace=module.__dict__
    if not _atan_diffgeom_namespace(namespace):return False
    anchor=namespace.get('_SCALAR_METADATA')
    if _DIFFGEOM_TYPE(anchor) is not _DIFFGEOM_TUPLE or len(anchor)!=2:return False
    owner,reader=anchor
    if (owner is not _rational_core or _DIFFGEOM_TYPE(owner) is not _DIFFGEOM_MODULE_TYPE
            or not _atan_diffgeom_namespace(owner.__dict__)
            or owner.__dict__.get('_ddm_function_metadata') is not reader
            or _DIFFGEOM_TYPE(reader) is not _DIFFGEOM_BUILTIN
            or reader.__self__ is not owner or reader.__name__!='_ddm_function_metadata'):return False
    installer=namespace.get('install')
    if _DIFFGEOM_TYPE(installer) is not _DIFFGEOM_FUNCTION:return False
    metadata=reader(installer)
    if (metadata is None or metadata[0] is not namespace.get('_SCALAR_INSTALL_CODE')
            or metadata[1] is not namespace or metadata[2] is not None
            or metadata[3] is not None or metadata[4] is not None
            or not _atan_diffgeom_namespace(installer.__dict__)):return False
    if (_DIFFGEOM_TYPE(_real_builtins) is not _DIFFGEOM_MODULE_TYPE
            or not _atan_diffgeom_namespace(_real_builtins.__dict__)
            or installer.__builtins__ is not _real_builtins.__dict__):return False
    proofs=namespace.get('_SCALAR_DISPATCH_PROOFS');aux=namespace.get('_SCALAR_AUX_PROOFS')
    if (_DIFFGEOM_TYPE(proofs) is not _DIFFGEOM_TUPLE or len(proofs)!=4
            or _DIFFGEOM_TYPE(aux) is not _DIFFGEOM_TUPLE or len(aux)!=2
            or installer.__dict__.get('_stock_proofs') is not proofs
            or installer.__dict__.get('_stock_aux_proofs') is not aux):return False
    return namespace,reader,metadata[0],proofs,aux


def _atan_diffgeom_states(states,reader):
    if _DIFFGEOM_TYPE(states) is not _DIFFGEOM_TUPLE:return False
    try:
        for state in states:
            if _DIFFGEOM_TYPE(state) is not _DIFFGEOM_TUPLE or len(state)!=6:return False
            function,code,closure,contents,defaults,keywords=state
            if (_DIFFGEOM_TYPE(function) is not _DIFFGEOM_FUNCTION
                    or _DIFFGEOM_TYPE(contents) is not _DIFFGEOM_TUPLE):return False
            current=reader(function)
            if (current is None or current[0] is not code or current[4] is not closure
                    or current[2] is not defaults or current[3] is not keywords
                    or closure is not None and _DIFFGEOM_TYPE(closure) is not _DIFFGEOM_TUPLE
                    or defaults is not None and _DIFFGEOM_TYPE(defaults) is not _DIFFGEOM_TUPLE
                    or keywords is not None and _DIFFGEOM_TYPE(keywords) is not _DIFFGEOM_DICT
                    or len(closure or ())!=len(contents)
                    or not _ATAN_ALL(cell.cell_contents is value for cell,value in zip(closure or (),contents))):return False
    except ValueError:return False
    return True


def _atan_diffgeom_template(data,factory,name):
    namespace,reader,installer,proofs,aux=data
    code_type=_DIFFGEOM_TYPE(installer)
    for constant in installer.co_consts:
        if _DIFFGEOM_TYPE(constant) is not code_type:continue
        if factory is None and constant.co_name==name:return constant
        if constant.co_name==factory:
            for nested in constant.co_consts:
                if _DIFFGEOM_TYPE(nested) is code_type and nested.co_name==name:return nested
    return None


def _atan_diffgeom_operation(function,template,states,data,attribute,proof):
    namespace,reader,installer,proofs,aux=data
    if _DIFFGEOM_TYPE(function) is not _DIFFGEOM_FUNCTION:return False
    metadata=reader(function)
    if (template is None or metadata is None or metadata[0] is not template
            or metadata[1] is not namespace or metadata[2] is not None or metadata[3] is not None
            or function.__builtins__ is not _real_builtins.__dict__
            or not _atan_diffgeom_namespace(function.__dict__)
            or function.__dict__.get(attribute) is not proof
            or not _atan_diffgeom_states(states,reader)):return False
    return metadata


def _atan_diffgeom_proof(proof,data,name,scalar_exclusions):
    if _DIFFGEOM_TYPE(proof) is not _DIFFGEOM_TUPLE or len(proof)!=8:return False
    operation,delegate,states,namespace,tensor,mro,isinstance_,operation_name=proof
    if (namespace is not data[0] or _DIFFGEOM_TYPE(operation_name) is not _DIFFGEOM_STR
            or operation_name!=name):return False
    template=_atan_diffgeom_template(data,'arithmetic','operation')
    metadata=_atan_diffgeom_operation(operation,template,states,data,'_symox_diffgeom_proof',proof)
    if metadata is False:return False
    if (_DIFFGEOM_TYPE(tensor) is not _DIFFGEOM_TYPE or _DIFFGEOM_TYPE(mro) is not _DIFFGEOM_TUPLE
            or len(tensor.__mro__)!=len(mro)or not _ATAN_ALL(a is b for a,b in zip(tensor.__mro__,mro))
            or _DIFFGEOM_TYPE(isinstance_) is not _DIFFGEOM_BUILTIN
            or isinstance_.__self__ is not _real_builtins or isinstance_.__name__!='isinstance'
            or namespace.get('isinstance',_real_builtins.__dict__.get('isinstance')) is not isinstance_
            or template.co_freevars!=('_Tensor','name','original')):return False
    closure=metadata[4]or ()
    try:
        if (len(closure)!=3 or closure[0].cell_contents is not tensor
                or _DIFFGEOM_TYPE(closure[1].cell_contents) is not _DIFFGEOM_STR
                or closure[1].cell_contents!=name or closure[2].cell_contents is not delegate):return False
    except ValueError:return False
    scalar_exclusions[id(tensor)]=tensor
    return True


def _atan_diffgeom_auxiliary(cls,old,current):
    if cls is not _ATAN_OPERATORS[0][0]:return False
    data=_atan_diffgeom_data()
    if data is None or data is False:return False
    namespace,reader,installer,proofs,aux=data
    seen=[]
    for proof in aux:
        if _DIFFGEOM_TYPE(proof) is not _DIFFGEOM_TUPLE or len(proof)!=5:return False
        name,operation,delegate,states,scope=proof
        if (_DIFFGEOM_TYPE(name) is not _DIFFGEOM_STR or name not in ('atoms','rcall')
                or name in seen or scope is not namespace or delegate is not old.get(name)
                or current.get(name) is not operation):return False
        seen.append(name)
        template=_atan_diffgeom_template(data,None,'atoms')if name=='atoms'else _atan_diffgeom_template(data,'evaluation','operation')
        metadata=_atan_diffgeom_operation(operation,template,states,data,'_symox_diffgeom_aux_proof',proof)
        if metadata is False:return False
        closure=metadata[4]or ()
        freevar='original_atoms'if name=='atoms'else'original'
        try:
            index=template.co_freevars.index(freevar)
            if closure[index].cell_contents is not delegate:return False
        except (ValueError,IndexError):return False
    return len(seen)==2


def _atan_diffgeom_atoms(method,original):
    # The rational driver's real branch calls atoms() without selectors.
    # The certified geometry wrapper forwards that exact call to its delegate.
    data=_atan_diffgeom_data()
    if data is None or data is False:return False
    for proof in data[4]:
        if _DIFFGEOM_TYPE(proof) is not _DIFFGEOM_TUPLE or len(proof)!=5:return False
        name,operation,delegate,states,scope=proof
        if _DIFFGEOM_TYPE(name) is not _DIFFGEOM_STR:return False
        if name!='atoms':continue
        if operation is not method or delegate is not original or scope is not data[0]:return False
        code=_atan_diffgeom_template(data,None,'atoms')
        metadata=_atan_diffgeom_operation(operation,code,states,data,'_symox_diffgeom_aux_proof',proof)
        if metadata is False:return False
        try:
            index=code.co_freevars.index('original_atoms')
            return metadata[4][index].cell_contents is original
        except (ValueError,IndexError):return False
    return False


def _atan_diffgeom_class(cls,bases,items):
    if (_DIFFGEOM_TYPE(cls) is not _DIFFGEOM_TYPE
            or _DIFFGEOM_TYPE(bases) is not _DIFFGEOM_TUPLE
            or _DIFFGEOM_TYPE(items) is not _DIFFGEOM_TUPLE):return False
    current=cls.__dict__;parents=cls.__bases__
    if len(parents)!=len(bases)or not _ATAN_ALL(a is b for a,b in zip(parents,bases)):return False
    if len(current)!=len(items):return False
    for item in items:
        if (_DIFFGEOM_TYPE(item) is not _DIFFGEOM_TUPLE or len(item)!=2
                or _DIFFGEOM_TYPE(item[0]) is not _DIFFGEOM_STR
                or current.get(item[0]) is not item[1]):return False
    return True


def _atan_diffgeom_holder(data):
    # RF intermediates are polynomials in the binding's non-host parameters.
    # Every canonical field atom below is host: Symbol equality rejects
    # Some(host)/None before any HostClass.same callback. No field wrapper
    # can therefore be selected by wrap for these non-host intermediates.
    import sys
    namespace,reader,installer,proofs,aux=data
    if not _atan_diffgeom_namespace(sys.modules):return False
    coordinates=sys.modules.get('symox.diffgeom._coordinates')
    if (_DIFFGEOM_TYPE(coordinates) is not _DIFFGEOM_MODULE_TYPE
            or _DIFFGEOM_TYPE(_lrt_public) is not _DIFFGEOM_MODULE_TYPE):return False
    scope=coordinates.__dict__;root=_lrt_public.__dict__
    if not _atan_diffgeom_namespace(scope)or not _atan_diffgeom_namespace(root):return False
    proof=scope.get('_FIELD_ATOM_PROOF')
    if _DIFFGEOM_TYPE(proof) is not _DIFFGEOM_TUPLE or len(proof)!=4:return False
    cls,bases,items,saved=proof
    holder=scope.get('_field_atoms')
    if (scope.get('_FieldAtomSet') is not cls or _DIFFGEOM_TYPE(holder) is not cls
            or root.get('_diffgeom_field_atoms') is not holder
            or not _atan_diffgeom_class(cls,bases,items)
            or len(bases)!=1 or bases[0] is not _DIFFGEOM_SET):return False
    contains=cls.__dict__.get('__contains__')
    if (_DIFFGEOM_TYPE(saved) is not _DIFFGEOM_TUPLE or len(saved)!=5
            or _DIFFGEOM_TYPE(contains) is not _DIFFGEOM_FUNCTION):return False
    current=reader(contains)
    if (current is None or current[0] is not saved[0] or current[1] is not scope
            or current[2] is not None or current[3] is not None or current[4] is not None
            or contains.__builtins__ is not _real_builtins.__dict__
            or scope.get('set',_real_builtins.__dict__.get('set')) is not _DIFFGEOM_SET
            or _real_builtins.__dict__.get('set') is not _DIFFGEOM_SET
            or scope.get('TypeError',_real_builtins.__dict__.get('TypeError')) is not _real_builtins.__dict__.get('TypeError')):return False
    field=namespace.get('_SCALAR_FIELD_CLASS')
    if _DIFFGEOM_TYPE(field) is not _DIFFGEOM_TUPLE or len(field)!=4:return False
    field_cls,field_bases,field_items,field_states=field
    if (root.get('_diffgeom_callable_type') is not field_cls
            or not _atan_diffgeom_class(field_cls,field_bases,field_items)
            or len(field_bases)!=1 or field_bases[0] is not _ATAN_OPERATORS[0][0]
            or not _atan_diffgeom_states(field_states,reader)):return False
    host_reader=namespace.get('_SCALAR_HOST_READER')
    if (_DIFFGEOM_TYPE(host_reader) is not _DIFFGEOM_BUILTIN
            or host_reader.__self__ is not _rational_core or host_reader.__name__!='_host_contains'
            or not _atan_diffgeom_namespace(_rational_core.__dict__)
            or _rational_core.__dict__.get('_host_contains') is not host_reader):return False
    for atom in _DIFFGEOM_SET.__iter__(holder):
        if _DIFFGEOM_TYPE(atom) is not field_cls:return False
        flags=host_reader(atom)
        if (_DIFFGEOM_TYPE(flags) is not _DIFFGEOM_TUPLE or len(flags)!=2
                or flags[0] is not True or flags[1] is not True):return False
    return True


def _atan_forwarding_chain_with_metadata(current,original,name,ore_proofs,ore_code,vector_proofs,vector_code,vector_types,diffgeom,_metadata):
    seen=set();forwarded=0
    while id(current) not in seen:
        seen.add(id(current))
        ore=next((proof for proof in ore_proofs if current is proof[0]),None)
        vector=next((proof for proof in vector_proofs if current is proof[0]),None)
        geometric=next((proof for proof in diffgeom[3] if _DIFFGEOM_TYPE(proof) is _DIFFGEOM_TUPLE and len(proof)==8 and current is proof[0]),None)if diffgeom is not None else None
        if ore is not None:
            if (_ATAN_GETATTR(current,'_symox_ore_proof',None) is not ore
                    or not (_atan_forwarding_proof(ore,ore_code,name) if _metadata is None else _atan_forwarding_proof_with_metadata(ore,ore_code,name,_metadata))):return False,forwarded
            delegate=ore[1];forwarded|=1
        elif vector is not None:
            if (_ATAN_GETATTR(current,'_symox_vector_proof',None) is not vector
                    or not (_atan_vector_proof(vector,vector_code,name,vector_types) if _metadata is None else _atan_vector_proof_with_metadata(vector,vector_code,name,vector_types,_metadata))):return False,forwarded
            delegate=vector[1];forwarded|=2
        elif geometric is not None:
            if not _atan_diffgeom_proof(geometric,diffgeom,name,vector_types):return False,forwarded
            delegate=geometric[1];forwarded|=4
        else:return current is original,forwarded
        if current is original:return True,forwarded
        current=delegate
    return False,forwarded


def _atan_stock_operators_with_metadata(f,g,operators,states,_metadata):
    import sys
    from types import CodeType,FunctionType,GetSetDescriptorType
    if operators is None:operators=_ATAN_OPERATORS
    if states is None:states=_ATAN_BOOTSTRAP_STATES
    if _ATAN_TYPE(operators) is not tuple or not (_atan_states_match(states) if _metadata is None else _atan_states_match_with_metadata(states,_metadata)):return False
    diffgeom=_atan_diffgeom_data()
    if diffgeom is False:return False
    if diffgeom is not None and not _atan_diffgeom_holder(diffgeom):return False
    recurrence=sys.modules.get('symox.holonomic.recurrence')
    proofs=()
    if recurrence is not None:
        if _metadata is not None and _ATAN_TYPE(recurrence) is not _ATAN_MODULE_TYPE:return False
        if _metadata is not None and not _atan_diffgeom_namespace(recurrence.__dict__):return False
        proofs=recurrence.__dict__.get('_SCALAR_DISPATCH_PROOFS')
        installer=recurrence.__dict__.get('_install_scalar_dispatch')
        _installer_values=_metadata(installer) if _metadata is not None else None
        if _metadata is not None and _installer_values is None:return False
        if (_ATAN_TYPE(proofs) is not tuple or _ATAN_TYPE(installer) is not FunctionType
                or proofs is not _ATAN_GETATTR(installer,'_stock_proofs',None)):return False
        factory=next((code for code in (_installer_values[0] if _metadata is not None else installer.__code__).co_consts
                      if _ATAN_TYPE(code) is CodeType and code.co_name=='factory'),None)
        if factory is None:return False
        operation_code=next((code for code in factory.co_consts
                            if _ATAN_TYPE(code) is CodeType and code.co_name=='operation'),None)
        if operation_code is None:return False
    vector=sys.modules.get('symox.vector._adapters')
    vector_proofs=();vector_code=None
    if vector is not None:
        if _metadata is not None and _ATAN_TYPE(vector) is not _ATAN_MODULE_TYPE:return False
        if _metadata is not None and not _atan_diffgeom_namespace(vector.__dict__):return False
        namespace=vector.__dict__;vector_proofs=namespace.get('_SCALAR_DISPATCH_PROOFS')
        installer=namespace.get('install');vector_code=namespace.get('_SCALAR_OPERATION_CODE')
        _installer_values=_metadata(installer) if _metadata is not None else None
        if _metadata is not None and _installer_values is None:return False
        if (_ATAN_TYPE(vector_proofs) is not tuple or _ATAN_TYPE(installer) is not FunctionType
                or (_installer_values[0] if _metadata is not None else installer.__code__) is not namespace.get('_SCALAR_INSTALL_CODE')
                or vector_proofs is not _ATAN_GETATTR(installer,'_stock_proofs',None)
                or _ATAN_TYPE(vector_code) is not CodeType):return False
    if (not _ATAN_ALL(_ATAN_TYPE(p) is tuple and len(p)==8 for p in proofs)
            or not _ATAN_ALL(_ATAN_TYPE(p) is tuple and len(p)==8 for p in vector_proofs)):return False
    verified={};forwarded=0;vector_types={}
    for cls,name,original in operators:
        current=_ATAN_GETATTR(cls,name)
        key=(id(current),id(original),name)
        if key in verified:
            if not verified[key]:return False
            continue
        valid,has_forwarder=(_atan_forwarding_chain_with_metadata(current,original,name,proofs,
                                                  operation_code if recurrence is not None else None,
                                                  vector_proofs,vector_code,vector_types,diffgeom,None) if _metadata is None else _atan_forwarding_chain_with_metadata(current,original,name,proofs,
                                                  operation_code if recurrence is not None else None,
                                                  vector_proofs,vector_code,vector_types,diffgeom,_metadata))
        verified[key]=valid
        if not valid:return False
        forwarded|=has_forwarder
    # Inspect each scalar MRO and each base dictionary once. Priority and
    # attribute hooks apply to all delegates; holonomic markers apply only
    # when an actual Ore layer was authenticated. A vector-only chain never
    # reads that marker in its normal scalar path.
    has_ore=forwarded&1
    bases=tuple(vector_types.values())
    for cls,priority in _ATAN_PRIORITIES:
        mro=_ATAN_MRO_DESCRIPTOR.__get__(cls)
        value=_ATAN_MISSING;attribute=_ATAN_MISSING
        found_value=False;found_attribute=False
        for base in mro:
            namespace=_ATAN_DICT_DESCRIPTOR.__get__(base)
            if not found_value and '_op_priority' in namespace:
                value=namespace['_op_priority'];found_value=True
            if not found_attribute and '__getattribute__' in namespace:
                attribute=namespace['__getattribute__'];found_attribute=True
            if '__getattr__' in namespace:return False
            if has_ore and '_symox_holonomic_function' in namespace:return False
            if bases:
                for basis in bases:
                    if base is basis:return False
        if value is not priority or attribute is not _ATAN_OBJECT_GETATTRIBUTE:return False
        if has_ore and _ATAN_GETATTR(cls,'__getattribute__') is not _ATAN_OBJECT_GETATTRIBUTE:return False
    if not forwarded:return True
    if f is None and g is None:return True
    if (_ATAN_TYPE(f) is not _ATAN_POLY or _ATAN_TYPE(g) is not _ATAN_POLY
            or _ATAN_TYPE(_ATAN_GENS) is not GetSetDescriptorType or Poly.gens is not _ATAN_GENS):return False
    # Read the stock native descriptor for selection only. Host generators
    # retain the full algorithm, including their arithmetic/attribute hooks.
    for polynomial in (f,g):
        for gen in _ATAN_GENS.__get__(polynomial,Poly):
            cls=_ATAN_TYPE(gen)
            if cls is not _ATAN_SCALARS[0] and cls is not _ATAN_SCALARS[1]:return False
    return True


def _stock_lrt_as_poly():
    # The public Expr facade can override a scalar method without changing
    # the storage base's descriptor. Check the methods actually resolved by
    # the native scalar classes produced by LRT, without invoking descriptors.
    for cls in (_rational_core.Expr,)+_ATAN_SCALARS:
        method=_atan_raw_attribute(cls,'as_poly')
        if (type(method) is not _lrt_types.FunctionType
                or method.__code__ is not _LRT_AS_POLY_CODE
                or method.__globals__ is not _lrt_public.__dict__):return False
    return True


def _stock_lrt_resultant():
    if resultant is _LRT_RESULTANT:return True
    # The ordering-options adapter is installed after this module's first
    # bootstrap import. Reload then sees that official wrapper instead of
    # its original builtin. LRT supplies only includePRS: the verified
    # wrapper's first branch forwards directly, without invoking any hook.
    if (type(resultant) is not _lrt_types.FunctionType
            or resultant.__code__ is not _LRT_WRT_CODE
            or resultant.__globals__ is not _LRT_WRT_GLOBALS
            or resultant.__defaults__ is not None
            or resultant.__kwdefaults__ is not None):return False
    try:
        cells=dict(zip(_LRT_WRT_CODE.co_freevars,resultant.__closure__ or ()))
        return (set(cells)=={'name','native'}
                and type(cells['name'].cell_contents) is str
                and cells['name'].cell_contents=='resultant'
                and cells['native'].cell_contents is _LRT_RESULTANT)
    except (AttributeError,ValueError,TypeError):return False


def _stock_atan_protocol(f=None,g=None):
    return (_ATAN_BOOTSTRAP_CAPTURE and atan is _ATAN_HEAD and Poly is _ATAN_POLY and _atan_evaluation.get()
            and _ATAN_TYPE(atan).__dict__.get('__call__') is _ATAN_CALL
            and atan.__dict__.get('_symox_factory') is _ATAN_FACTORY
            and Poly.as_expr is _ATAN_AS_EXPR
            and _rational_core._apply is _ATAN_APPLY
            and _atan_stock_operators(f,g)
            and log_to_atan is _ATAN_RECURSION)


def ratint(f, x, **flags):
    """
    Performs indefinite integration of rational functions.

    Explanation
    ===========

    Given a field :math:`K` and a rational function :math:`f = p/q`,
    where :math:`p` and :math:`q` are polynomials in :math:`K[x]`,
    returns a function :math:`g` such that :math:`f = g'`.

    Examples
    ========

    >>> from symox.integrals.rationaltools import ratint
    >>> from symox.abc import x

    >>> ratint(36/(x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2), x)
    (12*x + 6)/(x**2 - 1) + 4*log(x - 2) - 4*log(x + 1)

    References
    ==========

    .. [1] M. Bronstein, Symbolic Integration I: Transcendental
       Functions, Second Edition, Springer-Verlag, 2005, pp. 35-70

    See Also
    ========

    symox.integrals.integrals.Integral.doit
    symox.integrals.rationaltools.ratint_logpart
    symox.integrals.rationaltools.ratint_ratpart

    """
    if _stock_ratint_protocol(f):
        native = _rational_core._ratint_driver(f,x,flags)
        if native is not None:
            value,logarithmic = native
            if not logarithmic or _stock_ratint_log_protocol():
                return value

    if isinstance(f, tuple):
        p, q = f
    else:
        p, q = f.as_numer_denom()

    p, q = Poly(p, x, composite=False, field=True), Poly(q, x, composite=False, field=True)

    coeff, p, q = p.cancel(q)
    poly, p = p.div(q)

    result = poly.integrate(x).as_expr()

    if p.is_zero:
        return coeff*result

    g, h = ratint_ratpart(p, q, x)

    P, Q = h.as_numer_denom()

    P = Poly(P, x)
    Q = Poly(Q, x)

    q, r = P.div(Q)

    result += g + q.integrate(x).as_expr()

    if not r.is_zero:
        symbol = flags.get('symbol', 't')

        if not isinstance(symbol, Symbol):
            t = Dummy(symbol)
        else:
            t = symbol.as_dummy()

        L = ratint_logpart(r, Q, x, t)

        real = flags.get('real')

        if real is None:
            if isinstance(f, tuple):
                p, q = f
                atoms = p.atoms() | q.atoms()
            else:
                atoms = f.atoms()

            for elt in atoms - {x}:
                if not elt.is_extended_real:
                    real = False
                    break
            else:
                real = True

        eps = S.Zero

        if not real:
            for h, q in L:
                _, h = h.primitive()
                eps += RootSum(
                    q, Lambda(t, t*log(h.as_expr())), quadratic=True)
        else:
            for h, q in L:
                _, h = h.primitive()
                R = log_to_real(h, q, x, t)

                if R is not None:
                    eps += R
                else:
                    eps += RootSum(
                        q, Lambda(t, t*log(h.as_expr())), quadratic=True)

        result += eps

    return coeff*result


def ratint_ratpart(f, g, x):
    """
    Horowitz-Ostrogradsky algorithm.

    Explanation
    ===========

    Given a field K and polynomials f and g in K[x], such that f and g
    are coprime and deg(f) < deg(g), returns fractions A and B in K(x),
    such that f/g = A' + B and B has square-free denominator.

    Examples
    ========

        >>> from symox.integrals.rationaltools import ratint_ratpart
        >>> from symox.abc import x, y
        >>> from sympy import Poly
        >>> ratint_ratpart(Poly(1, x, domain='ZZ'),
        ... Poly(x + 1, x, domain='ZZ'), x)
        (0, 1/(x + 1))
        >>> ratint_ratpart(Poly(1, x, domain='EX'),
        ... Poly(x**2 + y**2, x, domain='EX'), x)
        (0, 1/(x**2 + y**2))
        >>> ratint_ratpart(Poly(36, x, domain='ZZ'),
        ... Poly(x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2, x, domain='ZZ'), x)
        ((12*x + 6)/(x**2 - 1), 12/(x**2 - x - 2))

    See Also
    ========

    ratint, ratint_logpart
    """
    f = Poly(f, x)
    g = Poly(g, x)

    from symox._core import _ratint_ratpart
    native = _ratint_ratpart(f, g, x)
    if native is not None:
        return native

    from symox.solvers.solvers import solve

    u, v, _ = g.cofactors(g.diff())

    n = u.degree()
    m = v.degree()

    A_coeffs = [ Dummy('a' + str(n - i)) for i in range(0, n) ]
    B_coeffs = [ Dummy('b' + str(m - i)) for i in range(0, m) ]

    C_coeffs = A_coeffs + B_coeffs

    A = Poly(A_coeffs, x, domain=ZZ[C_coeffs])
    B = Poly(B_coeffs, x, domain=ZZ[C_coeffs])

    H = f - A.diff()*v + A*(u.diff()*v).quo(u) - B*u

    result = solve(H.coeffs(), C_coeffs)

    A = A.as_expr().subs(result)
    B = B.as_expr().subs(result)

    rat_part = cancel(A/u.as_expr(), x)
    log_part = cancel(B/v.as_expr(), x)

    return rat_part, log_part


def ratint_logpart(f, g, x, t=None):
    r"""
    Lazard-Rioboo-Trager algorithm.

    Explanation
    ===========

    Given a field K and polynomials f and g in K[x], such that f and g
    are coprime, deg(f) < deg(g) and g is square-free, returns a list
    of tuples (s_i, q_i) of polynomials, for i = 1..n, such that s_i
    in K[t, x] and q_i in K[t], and::

                           ___    ___
                 d  f   d  \  `   \  `
                 -- - = --  )      )   a log(s_i(a, x))
                 dx g   dx /__,   /__,
                          i=1..n a | q_i(a) = 0

    Examples
    ========

    >>> from symox.integrals.rationaltools import ratint_logpart
    >>> from symox.abc import x
    >>> from sympy import Poly
    >>> ratint_logpart(Poly(1, x, domain='ZZ'),
    ... Poly(x**2 + x + 1, x, domain='ZZ'), x)
    [(Poly(x + 3*_t/2 + 1/2, x, domain='QQ[_t]'),
    ...Poly(3*_t**2 + 1, _t, domain='ZZ'))]
    >>> ratint_logpart(Poly(12, x, domain='ZZ'),
    ... Poly(x**2 - x - 2, x, domain='ZZ'), x)
    [(Poly(x - 3*_t/8 - 1/2, x, domain='QQ[_t]'),
    ...Poly(-_t**2 + 16, _t, domain='ZZ'))]

    See Also
    ========

    ratint, ratint_ratpart
    """
    f, g = Poly(f, x), Poly(g, x)

    t = t or Dummy('t')
    if (_ATAN_BOOTSTRAP_CAPTURE and Poly is _rational_core.Poly and Poly is _ATAN_POLY
            and _stock_lrt_resultant()
            and Poly.__init__ is object.__init__
            and _lrt_public.Poly is _rational_core.Poly
            and _lrt_constructor.Poly is _rational_core.Poly
            and _lrt_constructor.construct is _LRT_CONSTRUCT
            and _lrt_public.sympify is _LRT_SYMPIFY
            and type(_LRT_RESULTANT) is _lrt_types.BuiltinFunctionType
            and _LRT_RESULTANT.__name__ == 'resultant'
            and _LRT_RESULTANT.__module__ == 'symox._core'
            and _stock_lrt_as_poly()
            and _atan_evaluation.get()
            and (_ATAN_TYPE(t) is _ATAN_SCALARS[0] or _ATAN_TYPE(t) is _ATAN_SCALARS[1])
            and _atan_stock_operators(f,g,_LRT_OPERATORS,_LRT_BOOTSTRAP_STATES)):
        native = _rational_core._ratint_logpart(f, g, x, t)
        if native is not None:
            return native
    a, b = g, f - g.diff()*Poly(t, x)

    res, R = resultant(a, b, includePRS=True)
    res = Poly(res, t, composite=False)

    assert res, "BUG: resultant(%s, %s) cannot be zero" % (a, b)

    R_map, H = {}, []

    for r in R:
        R_map[r.degree()] = r

    def _include_sign(c, sqf):
        if c.is_extended_real and (c < 0) == True:
            h, k = sqf[0]
            c_poly = c.as_poly(h.gens)
            sqf[0] = h*c_poly, k

    C, res_sqf = res.sqf_list()
    _include_sign(C, res_sqf)

    for q, i in res_sqf:
        _, q = q.primitive()

        if g.degree() == i:
            H.append((g, q))
        else:
            h = R_map[i]
            h_lc = Poly(h.LC(), t, field=True)

            c, h_lc_sqf = h_lc.sqf_list(all=True)
            _include_sign(c, h_lc_sqf)

            for a, j in h_lc_sqf:
                h = h.quo(Poly(a.gcd(q)**j, x))

            inv, coeffs = h_lc.invert(q), [S.One]

            for coeff in h.coeffs()[1:]:
                coeff = coeff.as_poly(inv.gens)
                T = (inv*coeff).rem(q)
                coeffs.append(T.as_expr())

            h = Poly(dict(list(zip(h.monoms(), coeffs))), x)

            H.append((h, q))

    return H


def log_to_atan(f, g):
    """
    Convert complex logarithms to real arctangents.

    Explanation
    ===========

    Given a real field K and polynomials f and g in K[x], with g != 0,
    returns a sum h of arctangents of polynomials in K[x], such that:

                   dh   d         f + I g
                   -- = -- I log( ------- )
                   dx   dx        f - I g

    Examples
    ========

        >>> from symox.integrals.rationaltools import log_to_atan
        >>> from symox.abc import x
        >>> from sympy import Poly, sqrt, S
        >>> log_to_atan(Poly(x, x, domain='ZZ'), Poly(1, x, domain='ZZ'))
        2*atan(x)
        >>> log_to_atan(Poly(x + S(1)/2, x, domain='QQ'),
        ... Poly(sqrt(3)/2, x, domain='EX'))
        2*atan(2*sqrt(3)*x/3 + sqrt(3)/3)

    See Also
    ========

    log_to_real
    """
    if _stock_atan_protocol(f,g):
        native = _rational_core._log_to_atan(f, g)
        if native is not None:
            return native

    if f.degree() < g.degree():
        f, g = -g, f

    f = f.to_field()
    g = g.to_field()

    if g.degree() == 0 and f.domain.is_EX:
        # Division by a nonzero constant is exact. SymPy's EX elements are
        # canonicalised by cancel, so its f.div(g) leaves a zero remainder;
        # symox's EX arithmetic does not always recognise (a/c)*c - a as 0
        # for nested radicals, so divide by the constant directly.
        c = g.LC()
        p = Poly([coeff/c for coeff in f.all_coeffs()], *f.gens, domain=f.domain)
        return 2*atan(p.as_expr())

    p, q = f.div(g)

    if q.is_zero:
        return 2*atan(p.as_expr())
    else:
        s, t, h = g.gcdex(-f)
        u = (f*s + g*t).quo(h)
        A = 2*atan(u.as_expr())

        return A + log_to_atan(s, t)


_ATAN_RECURSION = log_to_atan


def _complete_real_roots(f, rs):
    """Whether ``rs`` (``roots(f, filter='R')``) holds every real root.

    Parametric coefficients keep the reference behaviour: realness of their
    roots depends on the parameters and is not decided here.
    """
    if f.free_symbols - set(f.gens):
        return True
    allr = roots(f)
    if sum(allr.values()) != f.degree():
        return False
    for r in allr:
        if r in rs or r.is_real is False:
            continue
        # The reference real-part tests (`D.evalf(chop=True)`) are numeric.
        im = r.evalf(50).as_real_imag()[1]
        if not im.is_Number or abs(im) < 1e-40:
            return False
    return True


def _get_real_roots(f, x):
    """get real roots of f if possible"""
    rs = roots(f, filter='R')

    try:
        num_roots = f.count_roots()
    except DomainError:
        # The real-root count is not available over this domain. SymPy 1.14
        # returns the filtered dict anyway, so a missed real root silently
        # drops terms (ratint(1/(x**8 + 1)) is 0). The dict is only used
        # when every root was found and each dropped one is non-real.
        return rs if _complete_real_roots(f, rs) else None
    else:
        if len(rs) == num_roots:
            return rs
        else:
            return None


def log_to_real(h, q, x, t):
    r"""
    Convert complex logarithms to real functions.

    Explanation
    ===========

    Given real field K and polynomials h in K[t,x] and q in K[t],
    returns real function f such that:
                          ___
                  df   d  \  `
                  -- = --  )  a log(h(a, x))
                  dx   dx /__,
                         a | q(a) = 0

    Examples
    ========

        >>> from symox.integrals.rationaltools import log_to_real
        >>> from symox.abc import x, y
        >>> from sympy import Poly, S
        >>> log_to_real(Poly(x + 3*y/2 + S(1)/2, x, domain='QQ[y]'),
        ... Poly(3*y**2 + 1, y, domain='ZZ'), x, y)
        2*sqrt(3)*atan(2*sqrt(3)*x/3 + sqrt(3)/3)/3
        >>> log_to_real(Poly(x**2 - 1, x, domain='ZZ'),
        ... Poly(-2*y + 1, y, domain='ZZ'), x, y)
        log(x**2 - 1)/2

    See Also
    ========

    log_to_atan
    """
    from symox.simplify.radsimp import collect
    if _stock_real_protocol(h,q,collect):
        native = _rational_core._log_to_real(h,q,x,t)
        if native is not None:
            return native
    u, v = symbols('u,v', cls=Dummy)

    H = h.as_expr().xreplace({t: u + I*v}).expand()
    Q = q.as_expr().xreplace({t: u + I*v}).expand()

    H_map = collect(H, I, evaluate=False)
    Q_map = collect(Q, I, evaluate=False)

    a, b = H_map.get(S.One, S.Zero), H_map.get(I, S.Zero)
    c, d = Q_map.get(S.One, S.Zero), Q_map.get(I, S.Zero)

    R = Poly(resultant(c, d, v), u)

    R_u = _get_real_roots(R, u)

    if R_u is None:
        return None

    result = S.Zero

    for r_u in R_u.keys():
        C = Poly(c.xreplace({u: r_u}), v)
        if not C:
            # t was split into real and imaginary parts
            # and denom Q(u, v) = c + I*d. We just found
            # that c(r_u) is 0 so the roots are in d
            C = Poly(d.xreplace({u: r_u}), v)
            # we were going to reject roots from C that
            # did not set d to zero, but since we are now
            # using C = d and c is already 0, there is
            # nothing to check
            d = S.Zero

        R_v = _get_real_roots(C, v)

        if R_v is None:
            return None

        R_v_paired = [] # take one from each pair of conjugate roots
        for r_v in R_v:
            if r_v not in R_v_paired and -r_v not in R_v_paired:
                if r_v.is_negative or r_v.could_extract_minus_sign():
                    R_v_paired.append(-r_v)
                elif not r_v.is_zero:
                    R_v_paired.append(r_v)

        for r_v in R_v_paired:

            D = d.xreplace({u: r_u, v: r_v})

            if D.evalf(chop=True) != 0:
                continue

            A = Poly(a.xreplace({u: r_u, v: r_v}), x)
            B = Poly(b.xreplace({u: r_u, v: r_v}), x)

            AB = (A**2 + B**2).as_expr()

            result += r_u*log(AB) + r_v*log_to_atan(A, B)

    R_q = _get_real_roots(q, t)

    if R_q is None:
        return None

    for r in R_q.keys():
        result += r*log(h.as_expr().subs(t, r))

    return result


# Dispatch identities are captured during bootstrap and retained on reload.
# These checks cover protocols exercised by the unchanged general path,
# including the temporary roots/collect step skipped by this bounded kernel.
_REAL_CAPTURE = globals().get('_REAL_CAPTURE', (
    symbols,Dummy,roots,log,
    _real_radsimp_module.collect,log.__dict__.get('_symox_factory'),
    _ATAN_TYPE(log).__dict__.get('__call__')))
_REAL_SCALAR_METHODS = globals().get('_REAL_SCALAR_METHODS', tuple(
    (cls,name,_atan_raw_attribute(cls,name))
    for cls in _ATAN_SCALARS
    for name in ('subs','xreplace','expand','evalf','could_extract_minus_sign','is_negative',
                 'args','func','is_Pow','is_Add','is_Mul','is_Function','__hash__','__eq__','__ne__')))
_REAL_DUMMY_METHODS = globals().get('_REAL_DUMMY_METHODS', tuple(
    (cls,name,_atan_raw_attribute(cls,name))
    for cls,names in ((Dummy,('__new__','__init__')),(_ATAN_TYPE(Dummy),('__call__','__instancecheck__','__subclasscheck__')))
    for name in names))
_REAL_PROTOCOL_GLOBALS = globals().get('_REAL_PROTOCOL_GLOBALS', tuple(
    (namespace,name,namespace.get(name,_real_builtins.__dict__.get(name)))
    for namespace,names in (
        (symbols.__globals__,('_symbols','Wild')),
        (_real_compat_module.__dict__,('_split_extra_assumptions','_symbol_instance','_core')),
        (_rational_core.__dict__,('Dummy','Expr')),
        (_real_subs_module.__dict__,('xreplace_rule','function_class','Mapping','_core','isinstance','any')),
        (_real_accum_module.__dict__,('getattr','isinstance','any','len','type','set')))
    for name in names))
_REAL_FUNCTION_STATES = globals().get('_REAL_FUNCTION_STATES',
    _atan_function_states((*_REAL_CAPTURE,_ATAN_AS_EXPR,*(value for _,_,value in _REAL_PROTOCOL_GLOBALS),
        *(value.__func__ if _ATAN_TYPE(value) is staticmethod else value
          for _,_,value in _REAL_DUMMY_METHODS))))
_REAL_SCALAR_STATES = globals().get('_REAL_SCALAR_STATES',
    _atan_function_states(method for _,_,method in _REAL_SCALAR_METHODS))
_REAL_LOCAL_ROOTS = _get_real_roots
_REAL_LOCAL_ATAN = log_to_atan
_REAL_LOCAL_STATES = _atan_function_states((_get_real_roots,log_to_atan))
_REAL_ACCUM_CODE = globals().get('_REAL_ACCUM_CODE', next(
    code for code in _real_accum_module.install_subs.__code__.co_consts
    if _ATAN_TYPE(code) is _lrt_types.CodeType and code.co_name=='subs'))
_REAL_ACCUM_HELPERS = globals().get('_REAL_ACCUM_HELPERS', tuple(
    (name,_real_accum_module.__dict__.get(name))
    for name in ('_needs','_substitutes_illegal','_illegal_values','reevaluate','HEADS','_core')))
_REAL_ACCUM_STATES = globals().get('_REAL_ACCUM_STATES',
    _atan_function_states(value for _,value in _REAL_ACCUM_HELPERS))


def _stock_real_subs(current,original):
    if current is original:return True
    if (_ATAN_TYPE(current) is not _lrt_types.FunctionType
            or current.__code__ is not _REAL_ACCUM_CODE
            or current.__globals__ is not _real_accum_module.__dict__
            or current.__defaults__ is not None or current.__kwdefaults__ is not None
            or not _atan_states_match(_REAL_ACCUM_STATES)
            or not _ATAN_ALL(_real_accum_module.__dict__.get(name) is value
                             for name,value in _REAL_ACCUM_HELPERS)):return False
    try:
        cells=dict(zip(_REAL_ACCUM_CODE.co_freevars,current.__closure__ or ()))
        created=cells['_CREATED'].cell_contents
        return (set(cells)=={'original','_CREATED'}
                and cells['original'].cell_contents is original
                and _ATAN_TYPE(created) is list and len(created)==1
                and _ATAN_TYPE(created[0]) is bool)
    except (KeyError,ValueError,TypeError):return False


def _stock_real_protocol(h,q,collect):
    return (_ATAN_BOOTSTRAP_CAPTURE and _atan_evaluation.get()
            and (symbols is _REAL_CAPTURE[0] and Dummy is _REAL_CAPTURE[1]
                 and roots is _REAL_CAPTURE[2] and log is _REAL_CAPTURE[3]
                 and collect is _REAL_CAPTURE[4]
                 and _get_real_roots is _REAL_LOCAL_ROOTS and log_to_atan is _REAL_LOCAL_ATAN)
            and log.__dict__.get('_symox_factory') is _REAL_CAPTURE[5]
            and _ATAN_TYPE(log).__dict__.get('__call__') is _REAL_CAPTURE[6]
            and _atan_states_match(_REAL_FUNCTION_STATES)
            and _ATAN_ALL(namespace.get(name,_real_builtins.__dict__.get(name)) is value
                          for namespace,name,value in _REAL_PROTOCOL_GLOBALS)
            and _ATAN_ALL(_atan_raw_attribute(cls,name) is value
                          for cls,name,value in _REAL_DUMMY_METHODS)
            and _atan_states_match(_REAL_LOCAL_STATES)
            and _atan_states_match(_REAL_SCALAR_STATES)
            and _stock_atan_protocol(h,q)
            and _atan_stock_operators(h,q,_LRT_OPERATORS,_LRT_BOOTSTRAP_STATES)
            and _ATAN_ALL((_stock_real_subs(_atan_raw_attribute(cls,name),original)
                           if name=='subs' else _atan_raw_attribute(cls,name) is original)
                          for cls,name,original in _REAL_SCALAR_METHODS))


# The driver publishes only a fully certified result. Basic construction and
# final scalar arithmetic are guarded for every route; log-specific protocols
# are checked only when Rust reports that the old logarithmic branch was used.
_ratint_modular = _real_importlib.import_module('symox._modular_poly')
_ratint_iterables = _real_importlib.import_module('symox.utilities.iterables')
_ratint_sympify = _real_importlib.import_module('symox._sympify')
_ratint_dispatch = _real_importlib.import_module('symox._native_dispatch')
_RATINT_INTEGER_CONVERTERS = globals().get('_RATINT_INTEGER_CONVERTERS',_ratint_dispatch.stock_integer_converters)
_RATINT_BASE_MODULES = globals().get('_RATINT_BASE_MODULES',tuple(
    (namespace,name,namespace.get(name,_real_builtins.__dict__.get(name)))
    for namespace,names in (
        (_lrt_constructor.__dict__,('Poly','construct','_options','_native_construct')),
        (_lrt_public.__dict__,('Poly','sympify')),
        (_ratint_modular.__dict__,('ring_domain',)),
        (_ratint_iterables.__dict__,('iterable',)),
        (_rational_core.__dict__,('_apply','_ratint_ratpart','_ratint_logpart','_log_to_atan','_log_to_real')),
        (globals(),('isinstance','tuple')),
        (_real_builtins.__dict__,('__import__',)))
    for name in names))
_RATINT_BASE_METHODS = globals().get('_RATINT_BASE_METHODS',tuple(
    (cls,name,_atan_raw_attribute(cls,name)) for cls in _ATAN_SCALARS
    for name in ('as_numer_denom','is_Matrix','is_Equality','_iterable','_symox_matrixexpr','__iter__')))
_RATINT_LOG_METHODS = globals().get('_RATINT_LOG_METHODS',tuple(
    (cls,name,_atan_raw_attribute(cls,name)) for cls in _ATAN_SCALARS
    for name in ('atoms','is_extended_real','as_dummy','name','is_Symbol','is_commutative')))
_RATINT_BASE_STATES = globals().get('_RATINT_BASE_STATES',_atan_function_states((
    *(value for _,_,value in _RATINT_BASE_MODULES),
    *(value for _,_,value in _RATINT_BASE_METHODS),
    _ATAN_AS_EXPR,_RATINT_INTEGER_CONVERTERS)))
_RATINT_LOG_STATES = globals().get('_RATINT_LOG_STATES',
    _atan_function_states(value for _,_,value in _RATINT_LOG_METHODS))
# These function definitions are intentionally renewed by an ordinary reload;
# user constructor/method captures above remain the original bootstrap values.

def _atan_states_match(states):
    return _atan_states_match_with_metadata(states,None)


def _atan_forwarding_proof(proof,operation_code,name):
    return _atan_forwarding_proof_with_metadata(proof,operation_code,name,None)


def _atan_vector_proof(proof,operation_code,name,vector_types):
    return _atan_vector_proof_with_metadata(proof,operation_code,name,vector_types,None)


def _atan_forwarding_chain(current,original,name,ore_proofs,ore_code,vector_proofs,vector_code,vector_types):
    return _atan_forwarding_chain_with_metadata(current,original,name,ore_proofs,ore_code,vector_proofs,vector_code,vector_types,None,None)


def _atan_stock_operators(f,g,operators=None,states=None):
    return _atan_stock_operators_with_metadata(f,g,operators,states,None)


def _stock_ratint_protocol(f):
    # Check the validators themselves before invoking any of them. This keeps
    # a real replacement/modified-code callback in the old public path once.
    for name,value,code,closure,defaults,kwdefaults in _RATINT_GUARD_FIELDS:
        if (globals().get(name) is not value or value.__code__ is not code
                or value.__closure__ is not closure or value.__defaults__ is not defaults
                or value.__kwdefaults__ is not kwdefaults):return False
    if (not _ATAN_BOOTSTRAP_CAPTURE or not _atan_evaluation.get()
            or _rational_core.__dict__.get('_ratint_driver') is None
            or Poly is not _ATAN_POLY or Poly is not _rational_core.Poly
            or Poly.__init__ is not object.__init__
            or not _ATAN_ALL(namespace.get(name,_real_builtins.__dict__.get(name)) is value
                            for namespace,name,value in _RATINT_BASE_MODULES)
            or (ratint_ratpart is not _RATINT_LOCAL_HELPERS[0]
                or ratint_logpart is not _RATINT_LOCAL_HELPERS[1]
                or log_to_real is not _RATINT_LOCAL_HELPERS[2]
                or log_to_atan is not _RATINT_LOCAL_HELPERS[3])
            or not _atan_states_match(_RATINT_BASE_STATES)
            or not _atan_states_match(_RATINT_LOCAL_STATES)
            or not _atan_states_match(_RATINT_GUARD_STATES)
            or not _atan_stock_operators(None,None,_LRT_OPERATORS,_LRT_BOOTSTRAP_STATES)
            or not _ATAN_ALL(_atan_raw_attribute(cls,name) is value
                            for cls,name,value in _RATINT_BASE_METHODS)):
        return False
    # Exact tuple metadata only: never iterate a generator or tuple subclass.
    # Native expressions return before converter lookup in public sympify.
    if _ATAN_TYPE(f) is tuple:
        if len(f)!=2:return False
        if (_ATAN_TYPE(f[0]) is int or _ATAN_TYPE(f[1]) is int):
            if not _RATINT_INTEGER_CONVERTERS(_ratint_sympify.converter):return False
    return True


_RATINT_LOCAL_HELPERS = (ratint_ratpart,ratint_logpart,log_to_real,log_to_atan)
_RATINT_GUARD_HELPERS = tuple((name,globals()[name]) for name in (
    '_atan_states_match','_atan_raw_attribute','_atan_stock_operators',
    '_atan_forwarding_chain','_atan_forwarding_proof','_atan_vector_proof',
    '_stock_lrt_resultant','_stock_lrt_as_poly','_stock_real_protocol','_stock_atan_protocol'))
_RATINT_GUARD_FIELDS = tuple((name,value,value.__code__,value.__closure__,value.__defaults__,value.__kwdefaults__)
    for name,value in _RATINT_GUARD_HELPERS)
_RATINT_GUARD_STATES = _atan_function_states(value for _,value in _RATINT_GUARD_HELPERS)
_RATINT_REGISTRY = globals().get('_RATINT_REGISTRY',S)
_RATINT_ZERO = globals().get('_RATINT_ZERO',S.Zero)
_RATINT_REGISTRY_METHODS = globals().get('_RATINT_REGISTRY_METHODS',tuple(
    (name,_atan_raw_attribute(_ATAN_TYPE(S),name)) for name in ('Zero','__dict__','__getattribute__')))
_RATINT_LOCAL_STATES = _atan_function_states(_RATINT_LOCAL_HELPERS)


def _stock_ratint_protocol_with_metadata(f,_metadata):
    # Check the validators themselves before invoking any of them. This keeps
    # a real replacement/modified-code callback in the old public path once.
    for name,value,code,closure,defaults,kwdefaults in _RATINT_GUARD_FIELDS:
        _guard_values=_metadata(value) if _metadata is not None else None
        if _metadata is not None and _guard_values is None:return False
        if (globals().get(name) is not value or (_guard_values[0] if _metadata is not None else value.__code__) is not code
                or (_guard_values[4] if _metadata is not None else value.__closure__) is not closure or (_guard_values[2] if _metadata is not None else value.__defaults__) is not defaults
                or (_guard_values[3] if _metadata is not None else value.__kwdefaults__) is not kwdefaults):return False
    if (not _ATAN_BOOTSTRAP_CAPTURE or not _atan_evaluation.get()
            or _rational_core.__dict__.get('_ratint_driver') is None
            or Poly is not _ATAN_POLY or Poly is not _rational_core.Poly
            or Poly.__init__ is not object.__init__
            or not _ATAN_ALL(namespace.get(name,_real_builtins.__dict__.get(name)) is value
                            for namespace,name,value in _RATINT_BASE_MODULES)
            or (ratint_ratpart is not _RATINT_LOCAL_HELPERS[0]
                or ratint_logpart is not _RATINT_LOCAL_HELPERS[1]
                or log_to_real is not _RATINT_LOCAL_HELPERS[2]
                or log_to_atan is not _RATINT_LOCAL_HELPERS[3])
            or not (_atan_states_match(_RATINT_BASE_STATES) if _metadata is None else _atan_states_match_with_metadata(_RATINT_BASE_STATES,_metadata))
            or not (_atan_states_match(_RATINT_LOCAL_STATES) if _metadata is None else _atan_states_match_with_metadata(_RATINT_LOCAL_STATES,_metadata))
            or not (_atan_states_match(_RATINT_GUARD_STATES) if _metadata is None else _atan_states_match_with_metadata(_RATINT_GUARD_STATES,_metadata))
            or not (_atan_stock_operators(None,None,_LRT_OPERATORS,_LRT_BOOTSTRAP_STATES) if _metadata is None else _atan_stock_operators_with_metadata(None,None,_LRT_OPERATORS,_LRT_BOOTSTRAP_STATES,_metadata))
            or not _ATAN_ALL(_atan_raw_attribute(cls,name) is value
                            for cls,name,value in _RATINT_BASE_METHODS)):
        return False
    # Exact tuple metadata only: never iterate a generator or tuple subclass.
    # Native expressions return before converter lookup in public sympify.
    if _ATAN_TYPE(f) is tuple:
        if len(f)!=2:return False
        if (_ATAN_TYPE(f[0]) is int or _ATAN_TYPE(f[1]) is int):
            if not _RATINT_INTEGER_CONVERTERS(_ratint_sympify.converter):return False
    return True


def _stock_ratint_log_protocol():
    if S is not _RATINT_REGISTRY:return False
    if not _ATAN_ALL(_atan_raw_attribute(_ATAN_TYPE(S),name) is value
                    for name,value in _RATINT_REGISTRY_METHODS):return False
    namespace=_ATAN_OBJECT_GETATTRIBUTE(S,'__dict__')
    if _ATAN_TYPE(namespace) is not dict or not _ATAN_ALL(_ATAN_TYPE(key) is str for key in namespace):return False
    if namespace.get('Zero',_RATINT_ZERO) is not _RATINT_ZERO:return False
    return (_stock_lrt_resultant() and _stock_lrt_as_poly()
            and _stock_real_protocol(None,None,_real_radsimp_module.collect)
            and _lrt_public.Dummy is Dummy
            and _atan_states_match(_RATINT_LOG_STATES)
            and _ATAN_ALL((_atan_raw_attribute(cls,name) is value
                           or name=='atoms'and _atan_diffgeom_atoms(_atan_raw_attribute(cls,name),value))
                          for cls,name,value in _RATINT_LOG_METHODS))
