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
