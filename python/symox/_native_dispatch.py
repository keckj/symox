"""Recognize stock scalar dispatchers without authorizing user callbacks."""


def stock_scalar_dispatch(current,standard,name):
    """Prove every optional Tensor/Vector dispatcher reaches the canonical API.

    These factories delegate ordinary native scalar expressions unchanged.
    Compare actual code and closure at each link: public names and copied
    functools.wraps metadata do not establish a safe scalar delegate.
    """
    if current is standard:return True
    import sys
    from types import CodeType,FunctionType
    codes=[]
    for module_name in ('symox.diffgeom._adapters','symox.vector._adapters'):
        adapter=sys.modules.get(module_name)
        if adapter is None or type(adapter.install) is not FunctionType:continue
        transform=next((c for c in adapter.install.__code__.co_consts
                        if type(c) is CodeType and c.co_name=='transform'),None)
        if transform is None:continue
        operation=next((c for c in transform.co_consts
                        if type(c) is CodeType and c.co_name=='operation'),None)
        if operation is not None:codes.append(operation)
    seen=set()
    while current is not standard:
        if (type(current) is not FunctionType or id(current) in seen
                or not any(current.__code__ is code for code in codes)):return False
        seen.add(id(current))
        try:
            closure=dict(zip(current.__code__.co_freevars,
                             (cell.cell_contents for cell in current.__closure__)))
        except ValueError:return False
        delegate=closure.get('original')
        if (type(closure.get('name')) is not str or closure['name']!=name
                or current.__dict__.get('__wrapped__') is not delegate):return False
        current=delegate
    return True


def stock_integer_converters(converters):
    """No applicable int/MRO converter and no key with custom lookup behavior.

    Inspect registry keys by identity only; looking up int in an arbitrary
    converter dictionary could itself invoke an unrelated key's equality hook.
    """
    if type(converters) is not dict:return False
    for cls,conversion in converters.items():
        if type(cls) is not type:return False
        if (cls is int or cls is object) and conversion is not None:return False
    return True


_TYPE=type
_STR=str
_LIST=list
_TUPLE=tuple
_LEN=len
_TYPE_DICT=type.__dict__['__dict__']
_TYPE_MRO=type.__dict__['__mro__']
_TYPE_NAME=type.__dict__['__name__']
_SLOTNAMES='__slotnames__'
_ISINSTANCE=isinstance
_ZIP=zip


def _stock_slotnames(cls):
    """``copyreg._slotnames(cls)`` computed from the class dictionaries only.

    None when a ``__slots__`` entry is not a plain str/tuple/list of str (the
    cache of such a class is not recognized). No attribute lookup, equality
    or iteration hook of a user object runs here.
    """
    names=[]
    for base in _TYPE_MRO.__get__(cls):
        namespace=_TYPE_DICT.__get__(base)
        if '__slots__' not in namespace:continue
        slots=namespace['__slots__']
        if _TYPE(slots) is _STR:slots=(slots,)
        elif _TYPE(slots) is not _TUPLE and _TYPE(slots) is not _LIST:return None
        for name in slots:
            if _TYPE(name) is not _STR:return None
            if name=='__dict__' or name=='__weakref__':continue
            if name.startswith('__') and not name.endswith('__'):
                stripped=_TYPE_NAME.__get__(base).lstrip('_')
                if stripped:name='_'+stripped+name
            names.append(name)
    return names


def copyreg_slot_cache(cls,value):
    """Whether ``value`` is the ``__slotnames__`` entry ``copyreg`` caches in
    the dictionary of ``cls`` on the first ``copy.copy``/``pickle`` of an
    instance: exactly the list ``copyreg._slotnames`` computes.

    The interpreter reads that entry only when copying/pickling, never when a
    native kernel or the Python algorithm runs, so a class whose dictionary
    differs from its snapshot by this cache alone keeps its protocol.
    """
    if _TYPE(value) is not _LIST or not _ISINSTANCE(cls,_TYPE):return False
    expected=_stock_slotnames(cls)
    if expected is None or _LEN(value)!=_LEN(expected):return False
    for found,name in _ZIP(value,expected):
        if _TYPE(found) is not _STR or found!=name:return False
    return True


def strip_slot_cache(cls,current,items):
    """``(current, items)`` without copyreg's ``__slotnames__`` cache.

    ``current`` is the present dictionary of ``cls`` and ``items`` its
    captured ``(name, value)`` snapshot. Both are returned unchanged when
    neither holds the entry; None when ``current`` holds a ``__slotnames__``
    entry that is not the cache (a real change of the class).
    """
    found=_SLOTNAMES in current
    if found:
        if not copyreg_slot_cache(cls,current[_SLOTNAMES]):return None
        current={name:value for name,value in current.items()
                 if not (_TYPE(name) is _STR and name==_SLOTNAMES)}
    kept=_TUPLE((name,value) for name,value in items if not (_TYPE(name) is _STR and name==_SLOTNAMES))
    if _LEN(kept)!=_LEN(items):items=kept
    return current,items
