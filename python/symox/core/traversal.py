"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox.utilities.traversal', 'bottom_up postorder_traversal preorder_traversal use')


def iterargs(expr):
    from collections import deque
    queue=deque([expr])
    while queue:
        value=queue.popleft();yield value;queue.extend(getattr(value,'args',()))

def iterfreeargs(expr,_first=True):
    from collections import deque
    queue=deque([(expr,frozenset())])
    while queue:
        value,bound=queue.popleft()
        if not (getattr(value,'free_symbols',set())&bound):yield value
        new_bound=bound|frozenset(getattr(value,'bound_symbols',()))
        for arg in getattr(value,'args',()):queue.append((arg,new_bound))

def walk(e,*target):
    if isinstance(e,target):
        yield e
        for arg in e.args:yield from walk(arg,*target)
__all__ += ['iterargs','iterfreeargs','walk']
