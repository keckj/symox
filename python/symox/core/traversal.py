"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox.utilities.traversal', 'bottom_up postorder_traversal preorder_traversal use')


def iterargs(expr):
    from collections import deque
    queue=deque([expr])
    while queue:
        value=queue.popleft();yield value;queue.extend(getattr(value,'args',()))

def iterfreeargs(expr,_first=True):
    # sympy/core/traversal.py: the arguments of a binding node are taken
    # from its canonical dummy form, without those involving the dummies.
    args=[expr]
    for i in args:
        yield i
        if _first and hasattr(i,'bound_symbols'):
            void=i.canonical_variables.values()
            for i in iterfreeargs(i.as_dummy(),_first=False):
                if not i.has(*void):
                    yield i
        args.extend(i.args)

def walk(e,*target):
    if isinstance(e,target):
        yield e
        for arg in e.args:yield from walk(arg,*target)
__all__ += ['iterargs','iterfreeargs','walk']
