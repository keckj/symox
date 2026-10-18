from .. import vectorize
__all__ = ['vectorize']


def _iterable(value):return hasattr(value,'__iter__') and not isinstance(value,(str,bytes))
def iter_copy(structure):return [iter_copy(value) if _iterable(value) else value for value in structure]
def structure_copy(structure):return structure.copy() if hasattr(structure,'copy') else iter_copy(structure)
def apply_on_element(function,args,kwargs,n):
    args=list(args);kwargs=dict(kwargs);structure=args[n] if isinstance(n,int) else kwargs[n]
    def recurse(value):
        if _iterable(value):return [recurse(item) for item in value]
        if isinstance(n,int):args[n]=value
        else:kwargs[n]=value
        return function(*args,**kwargs)
    return [recurse(value) for value in structure]
__all__ += ['apply_on_element','iter_copy','structure_copy']
