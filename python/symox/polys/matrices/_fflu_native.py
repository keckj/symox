"""Admission/publication only; exact FFLU arithmetic lives in Rust."""
from types import ModuleType,WrapperDescriptorType,MethodDescriptorType,GetSetDescriptorType
import builtins as _builtins
from symox import _core
from ..domains import _element_print as _ep, _algebraic_methods as _am

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
    _RANGE,_ALL,_ENUMERATE=range,all,enumerate
    _DDM=None


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


# Immutable semantic templates prevent reimport from blessing changed code.
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


def _expected_exquo(self, a, b):
    try:
        return _native_exquo(self, a, b)
    except Exception as exc:
        from ... import ExactQuotientFailed
        if not isinstance(exc, ExactQuotientFailed):
            raise
        if not getattr(self, 'is_ZZ_I', False):
            # SymPy's ExactQuotientFailed(a, b, dom) text.
            try:
                text = '%s does not divide %s in %s' % (self.convert(b), self.convert(a), self)
            except Exception:
                raise exc from None
            raise ExactQuotientFailed(text) from None
        a, b = self.convert(a), self.convert(b)
        text = '(%d + %d*I) does not divide (%d + %d*I) in ZZ_I' % (
            int(b.x), int(b.y), int(a.x), int(a.y))
        raise ExactQuotientFailed(text) from None

if '_CAPTURED' not in globals():
    _CAPTURED=True
    _READY=False
    _CORE_NAMESPACE=_module_dict(_core)
    _PRINT_NAMESPACE=_module_dict(_ep)
    _ALGEBRAIC_NAMESPACE=_module_dict(_am)
    _BUILTINS_NAMESPACE=_module_dict(_builtins)
    _READ_FUNCTION=_CORE_NAMESPACE.get('_ddm_function_metadata') if _CORE_NAMESPACE is not None else None
    _METADATA_OK=(_TYPE(_READ_FUNCTION) is _BUILTIN and _READ_FUNCTION.__self__ is _core
                  and _READ_FUNCTION.__name__=='_ddm_function_metadata')
    _DOMAIN=_CORE_NAMESPACE.get('Domain') if _CORE_NAMESPACE is not None else None
    _ELEMENT=_CORE_NAMESPACE.get('DomainElement') if _CORE_NAMESPACE is not None else None
    _DOMAIN_DICT=_class_dict(_DOMAIN)
    _ELEMENT_DICT=_class_dict(_ELEMENT)
    _EXPECTED=[]
    _CODE_ANCHORS=[]
    _EXQUO_CODE=None
    _CODE_READERS=[]
    if (_METADATA_OK and _PRINT_NAMESPACE is not None and _ALGEBRAIC_NAMESPACE is not None and _BUILTINS_NAMESPACE is not None
            and _DOMAIN_DICT is not None and _ELEMENT_DICT is not None):
        for _name,_template in (('_int_congruent',_expected_int_congruent),('_eq',_expected_eq),('_ne',_expected_ne)):
            _data=_READ_FUNCTION(_template)
            if _data is None:break
            _EXPECTED.append((_name,_data[0]))
        if _EXPECTED.__len__()==3:
            # Fixed code descriptors avoid global getattr interposition.
            _CODE_TYPE=_TYPE(_EXPECTED[0][1])
            for _name in ('co_argcount','co_posonlyargcount','co_kwonlyargcount','co_nlocals','co_stacksize','co_flags',
                          'co_code','co_consts','co_names','co_varnames','co_freevars','co_cellvars'):
                _CODE_READERS.append((_name,_CODE_TYPE.__dict__[_name].__get__))
            # Code objects are immutable anchors, not cached function data.
            # Verify original semantics first so fresh imports cannot bless a
            # changed __code__; still reread all function fields on every call.
            for _name,_code in _EXPECTED:
                _actual=_PRINT_NAMESPACE.get(_name)
                if _TYPE(_actual) is _FUNCTION:
                    _data=_READ_FUNCTION(_actual)
                    if _data is not None and _same_code(_data[0],_code):_CODE_ANCHORS.append((_name,_data[0]))
            _data=_READ_FUNCTION(_expected_exquo)
            if _data is not None:
                _EXQUO_CODE=_data[0]
                _actual=_ALGEBRAIC_NAMESPACE.get('_exquo')
                if _TYPE(_actual) is _FUNCTION:
                    _data=_READ_FUNCTION(_actual)
                    if _data is not None and _same_code(_data[0],_EXQUO_CODE):_CODE_ANCHORS.append(('_exquo',_data[0]))
                _READY=True
    _EXPECTED=_TUPLE(_EXPECTED)
    _CODE_READERS=_TUPLE(_CODE_READERS)
    _CODE_ANCHORS=_TUPLE(_CODE_ANCHORS)
    _congruent=_PRINT_NAMESPACE.get('_int_congruent') if _PRINT_NAMESPACE is not None else None


def register(stock):
    global _DDM
    if _DDM is None:_DDM=stock


def _anchor(name):
    for key,code in _CODE_ANCHORS:
        if key==name:return code
    return None


def _metadata(function,namespace,code,canonical):
    if _TYPE(function) is not _FUNCTION:return False
    # A function retains its own f_builtins even when its globals dictionary's
    # __builtins__ entry is restored later. This public readonly member is
    # non-audited on the verified CPython3.10 interpreters.
    if function.__builtins__ is not _BUILTINS_NAMESPACE:return False
    data=_READ_FUNCTION(function)
    return (data is not None and (data[0] is canonical or _same_code(data[0],code)) and data[1] is namespace
            and data[2] is None and data[3] is None and data[4] is None)


def _protocol(k,stock,dependencies,phase):
    if not _READY or stock is not _DDM:return False
    core,printing,algebraic,builtins=(_module_dict(module) for module in (_core,_ep,_am,_builtins))
    if core is None or printing is None or algebraic is None or builtins is None:return False
    if (core.get('_ddm_function_metadata') is not _READ_FUNCTION or core.get('Domain') is not _DOMAIN
            or core.get('DomainElement') is not _ELEMENT or _TYPE(k) is not _DOMAIN):return False
    element=_class_dict(_ELEMENT);domain=_class_dict(_DOMAIN)
    if element is None or domain is None or '__del__' in element:return False
    for name in ('__mul__','__sub__'):
        if not _descriptor(element.get(name),_ELEMENT,name,WrapperDescriptorType):return False
    for name in ('zero','one'):
        if not _descriptor(domain.get(name),_DOMAIN,name,GetSetDescriptorType):return False
    exquo=domain.get('exquo')
    if not _descriptor(exquo,_DOMAIN,'exquo',MethodDescriptorType):
        if (exquo is not algebraic.get('_exquo') or not _metadata(exquo,algebraic,_EXQUO_CODE,_anchor('_exquo'))
                or not _descriptor(algebraic.get('_native_exquo'),_DOMAIN,'exquo',MethodDescriptorType)):return False
    if _TYPE(stock) is not _TYPE:return False
    bases=stock.__bases__
    if _TYPE(bases) is not _TUPLE or bases!=(_LIST,):return False
    slots=stock.__dict__
    for key in slots:
        if _TYPE(key) is not _STR:return False
    for name in ('__iter__','__len__','__getitem__','__setitem__','__new__'):
        if slots.get(name,_LIST.__dict__[name]) is not _LIST.__dict__[name]:return False
    for name in ('__getattribute__','__getattr__','__class__'):
        if name in slots:return False
    if _TYPE(dependencies) is not _TUPLE:return False
    expected=(_RANGE,_ALL,_TYPE,_INT) if phase==0 else (_ENUMERATE,_TYPE,_INT) if phase==1 else (_RANGE,_TYPE,_INT)
    if dependencies.__len__()!=expected.__len__():return False
    i=0
    for value in expected:
        if dependencies[i] is not value:return False
        i+=1
    if (_TYPE(_ALL) is not _BUILTIN or _ALL.__self__ is not _builtins or _ALL.__name__!='all'
            or _TYPE(_RANGE) is not _TYPE or _RANGE.__name__!='range' or _RANGE.__module__!='builtins'
            or _RANGE.__flags__ & (1<<9) or _TYPE(_ENUMERATE) is not _TYPE
            or _ENUMERATE.__name__!='enumerate' or _ENUMERATE.__module__!='builtins' or _ENUMERATE.__flags__ & (1<<9)):return False
    for name,code in _EXPECTED:
        function=printing.get(name)
        if not _metadata(function,printing,code,_anchor(name)):return False
        if name!='_int_congruent' and element.get('__'+name.strip('_')+'__') is not function:return False
    for name in ('_native_eq','_native_ne'):
        if not _descriptor(printing.get(name),_ELEMENT,'__'+name[8:]+'__',WrapperDescriptorType):return False
    owned_builtins=printing.get('__builtins__')
    if owned_builtins is not _builtins and owned_builtins is not builtins:return False
    return True


def try_eliminate(a,perm,k,shape,stock,dependencies):
    if not _protocol(k,stock,dependencies,0):return False
    return _core._ddm_fflu(a,perm,k,shape,stock)


def try_permutation(p,perm,k,shape,stock,dependencies):
    if not _protocol(k,stock,dependencies,1):return False
    return _core._ddm_fflu_p(p,perm,k,shape,stock)


def try_lower(u,l,k,shape,stock,dependencies):
    if not _protocol(k,stock,dependencies,2):return False
    return _core._ddm_fflu_l(u,l,k,shape,stock)


def try_diagonal(l,d,k,shape,stock,dependencies):
    if not _protocol(k,stock,dependencies,3):return False
    return _core._ddm_fflu_d(l,d,k,shape,stock)
