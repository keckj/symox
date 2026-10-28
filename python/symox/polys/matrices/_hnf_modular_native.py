"""Admission for exact ZZ modulo-D Hermite math in Rust; no arithmetic driver."""
from types import ModuleType,WrapperDescriptorType,MethodDescriptorType,GetSetDescriptorType
import builtins as _builtins
from symox import _core
from ..domains import _element_print as _ep,_algebraic_methods as _am

if '_INITIALISED' not in globals():
    _INITIALISED=True
    _TYPE=False.__class__.__class__
    _STR=''.__class__
    _INT=(0).__class__
    _BOOL=False.__class__
    _BYTES=b''.__class__
    _TUPLE=().__class__
    _LIST=[].__class__
    _DICT={}.__class__
    _OBJECT=_LIST.__bases__[0]
    _FUNCTION=(lambda:None).__class__
    _BUILTIN=[].append.__class__
    _RANGE,_LEN=range,len


def _module_dict(module):
    if _TYPE(module) is not ModuleType:return None
    namespace=module.__dict__
    for key in namespace:
        if _TYPE(key) is not _STR:return None
    return namespace

def _class_dict(cls):
    if _TYPE(cls) is not _TYPE:return None
    bases=cls.__bases__
    if _TYPE(bases) is not _TUPLE or bases!=(_OBJECT,):return None
    namespace=cls.__dict__
    for key in namespace:
        if _TYPE(key) is not _STR:return None
    for key in ('__getattribute__','__getattr__','__class__'):
        if key in namespace:return None
    return namespace

def _descriptor(value,cls,name,kind):
    return _TYPE(value) is kind and value.__objclass__ is cls and value.__name__==name

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

def _same_code(left,right):
    if _TYPE(left) is not _TYPE(right):return False
    for name,reader in _CODE_READERS:
        # Actual code objects are immutable, non-subclassable builtins. The
        # fields used here are not Python function.__code__ audit operations.
        if not _literal_equal(reader(left),reader(right)):return False
    return True

def _expected_int_congruent(self, other):
    """The reference ``ModularInteger == int``: an ``int`` is compared
    modulo the characteristic; ``None`` when that rule does not apply."""
    if type(other) is int or (type(other) is not type(self)
                              and getattr(other, 'is_Integer', False) is True):
        domain = self.parent()
        if domain.is_FiniteField:
            return (int(self) - int(other)) % int(domain.characteristic()) == 0
    return None

def _expected_eq(self, other):
    result = _int_congruent(self, other)
    return _native_eq(self, other) if result is None else result

def _expected_ne(self, other):
    result = _int_congruent(self, other)
    return _native_ne(self, other) if result is None else not result

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

def _expected_symmetric_residue(a, m):
    """The residue of a mod m in ``(-m/2, m/2]``."""
    if a <= m // 2:
        return a
    return a - m

def _expected_gaussian_binary(name):
    native = _native_ops[name]
    reflected = name.startswith('__r')
    op = name.strip('_').lstrip('r') if reflected else name.strip('_')

    def method(self, other):
        K = self.parent()
        if K.is_ZZ_I or K.is_QQ_I:
            if op == 'truediv' and K.is_ZZ_I:
                from ... import QQ_I
                K = QQ_I
            o = _gaussian_operand(K, other)
            if o is not None:
                a = K.convert(self)
                a, b = (o, a) if reflected else (a, o)
                _zero_check(a, b, op)
                if op == 'truediv':
                    return K.exquo(a, b)
                q, r = K.div(a, b)
                return {'floordiv': q, 'mod': r, 'divmod': (q, r)}[op]
            return NotImplemented
        if native is None:
            return NotImplemented
        return native(self, other)
    method.__name__ = name
    return method

if '_CAPTURED' not in globals():
    _CAPTURED=True
    _READY=False
    _REGISTERED=False
    _CORE_NAMESPACE=_module_dict(_core)
    _PRINT_NAMESPACE=_module_dict(_ep)
    _ALGEBRAIC_NAMESPACE=_module_dict(_am)
    _BUILTINS_NAMESPACE=_module_dict(_builtins)
    _READ_FUNCTION=_CORE_NAMESPACE.get('_ddm_function_metadata') if _CORE_NAMESPACE is not None else None
    _METADATA_OK=(_TYPE(_READ_FUNCTION) is _BUILTIN and _READ_FUNCTION.__self__ is _core
                  and _READ_FUNCTION.__name__=='_ddm_function_metadata')
    _DOMAIN=_CORE_NAMESPACE.get('Domain') if _CORE_NAMESPACE is not None else None
    _ELEMENT=_CORE_NAMESPACE.get('DomainElement') if _CORE_NAMESPACE is not None else None
    _ZZ=_CORE_NAMESPACE.get('ZZ') if _CORE_NAMESPACE is not None else None
    _CODE_READERS=[]
    _CODES=[]
    _WRAPPER_CODE=None
    if (_METADATA_OK and _PRINT_NAMESPACE is not None and _ALGEBRAIC_NAMESPACE is not None
            and _BUILTINS_NAMESPACE is not None and _class_dict(_DOMAIN) is not None and _class_dict(_ELEMENT) is not None):
        for _name,_template in (('_int_congruent',_expected_int_congruent),('_eq',_expected_eq),('_ne',_expected_ne),
                                ('_gcdex',_expected_gcdex),('add_columns',_expected_add_columns),('symmetric_residue',_expected_symmetric_residue)):
            _data=_READ_FUNCTION(_template)
            if _data is None:break
            _CODES.append((_name,_data[0]))
        if _CODES.__len__()==6:
            _CODE_TYPE=_TYPE(_CODES[0][1])
            for _name in ('co_argcount','co_posonlyargcount','co_kwonlyargcount','co_nlocals','co_stacksize','co_flags',
                          'co_code','co_consts','co_names','co_varnames','co_freevars','co_cellvars'):
                _CODE_READERS.append((_name,_CODE_TYPE.__dict__[_name].__get__))
            _factory=_READ_FUNCTION(_expected_gaussian_binary)
            if _factory is not None:
                for _constant in _factory[0].co_consts:
                    if _TYPE(_constant) is _CODE_TYPE and _constant.co_name=='method':_WRAPPER_CODE=_constant
            _READY=_WRAPPER_CODE is not None
    _CODES=_TUPLE(_CODES)
    _CODE_READERS=_TUPLE(_CODE_READERS)
    _congruent=_PRINT_NAMESPACE.get('_int_congruent') if _PRINT_NAMESPACE is not None else None


def _metadata(function,namespace,code):
    if _TYPE(function) is not _FUNCTION or function.__builtins__ is not _BUILTINS_NAMESPACE:return False
    data=_READ_FUNCTION(function)
    return (data is not None and _same_code(data[0],code) and data[1] is namespace
            and data[2] is None and data[3] is None and data[4] is None)


def register(namespace):
    global _REGISTERED,_NAMESPACE,_DEFAULTDICT,_SYMMETRIC_NAMESPACE
    if _REGISTERED or not _READY or _TYPE(namespace) is not _DICT:return
    for key in namespace:
        if _TYPE(key) is not _STR:return
    kind=namespace.get('defaultdict')
    if (_TYPE(kind) is not _TYPE or kind.__flags__ & (1<<9) or kind.__name__!='defaultdict'
            or kind.__module__!='collections' or kind.__bases__!=(_DICT,)):return
    symmetric=namespace.get('symmetric_residue')
    if _TYPE(symmetric) is not _FUNCTION:return
    data=_READ_FUNCTION(symmetric)
    if data is None or _TYPE(data[1]) is not _DICT:return
    _NAMESPACE=namespace;_DEFAULTDICT=kind;_SYMMETRIC_NAMESPACE=data[1];_REGISTERED=True


def _wrapper(element,algebraic,name):
    function=element.get(name)
    if _TYPE(function) is not _FUNCTION or function.__builtins__ is not _BUILTINS_NAMESPACE:return False
    data=_READ_FUNCTION(function)
    if (data is None or not _same_code(data[0],_WRAPPER_CODE) or data[1] is not algebraic
            or data[2] is not None or data[3] is not None):return False
    closure=data[4]
    if _TYPE(closure) is not _TUPLE or closure.__len__()!=3:return False
    try:native,operation,reflected=(cell.cell_contents for cell in closure)
    except ValueError:return False
    return (_descriptor(native,_ELEMENT,name,WrapperDescriptorType)
            and _TYPE(operation) is _STR and operation==name.strip('_') and reflected is False)


def try_hnf_modular(matrix,output,modulus,shape,dependencies):
    if not _READY or not _REGISTERED:return False
    core,printing,algebraic,builtins=(_module_dict(module) for module in (_core,_ep,_am,_builtins))
    if core is None or printing is None or algebraic is None or builtins is None:return False
    for key in _NAMESPACE:
        if _TYPE(key) is not _STR:return False
    if (core.get('_ddm_function_metadata') is not _READ_FUNCTION or core.get('Domain') is not _DOMAIN
            or core.get('DomainElement') is not _ELEMENT or core.get('ZZ') is not _ZZ or _TYPE(_ZZ) is not _DOMAIN):return False
    element=_class_dict(_ELEMENT);domain=_class_dict(_DOMAIN)
    if element is None or domain is None or '__del__' in element or '__ifloordiv__' in element:return False
    for name in ('__lt__','__le__','__neg__','__mul__','__rmul__','__add__','__sub__'):
        if not _descriptor(element.get(name),_ELEMENT,name,WrapperDescriptorType):return False
    if not _descriptor(element.get('parent'),_ELEMENT,'parent',MethodDescriptorType):return False
    if not _descriptor(domain.get('gcdex'),_DOMAIN,'gcdex',MethodDescriptorType):return False
    for name in ('is_FiniteField','is_ZZ_I','is_QQ_I'):
        if not _descriptor(domain.get(name),_DOMAIN,name,GetSetDescriptorType):return False
    if not _wrapper(element,algebraic,'__floordiv__') or not _wrapper(element,algebraic,'__mod__'):return False
    if _TYPE(dependencies) is not _TUPLE or dependencies.__len__()!=11:return False
    expected=(_RANGE,_LEN,_ZZ,_NAMESPACE.get('_gcdex'),_NAMESPACE.get('add_columns'),_NAMESPACE.get('symmetric_residue'),
              _ZZ,_RANGE,_LEN,_TYPE,_INT)
    i=0
    for value in expected:
        if dependencies[i] is not value:return False
        i+=1
    if (_TYPE(_RANGE) is not _TYPE or _RANGE.__flags__ & (1<<9) or _RANGE.__name__!='range' or _RANGE.__module__!='builtins'
            or _TYPE(_LEN) is not _BUILTIN or _LEN.__self__ is not _builtins or _LEN.__name__!='len'):return False
    for name,code in _CODES:
        namespace=printing if name in ('_int_congruent','_eq','_ne') else _NAMESPACE
        owner=_SYMMETRIC_NAMESPACE if name=='symmetric_residue' else namespace
        function=namespace.get(name)
        if not _metadata(function,owner,code):return False
        if name in ('_eq','_ne') and element.get('__'+name.strip('_')+'__') is not function:return False
    for name,slot in (('_native_eq','__eq__'),('_native_ne','__ne__')):
        if not _descriptor(printing.get(name),_ELEMENT,slot,WrapperDescriptorType):return False
    return _core._ddm_hnf_modular(matrix,output,modulus,shape,_DEFAULTDICT)


