"""Expression traversal and transformations on the native public tree."""
import operator


def _children(node, keys):
    if hasattr(node,'args'):
        children=node.args
    elif not isinstance(node,(str,bytes)):
        try:
            children=tuple(node)
        except TypeError:
            return ()
    else:
        return ()
    if keys:
        from .sorting import ordered
        return tuple(ordered(children,keys=None if keys is True else keys))
    return children


class preorder_traversal:
    """Iterate parent before children; ``skip()`` omits the last node's subtree."""
    def __init__(self, node, keys=None):
        self._stack=[iter((node,))]
        self._keys=keys
        self._previous=None
        self._expand=False

    def __iter__(self):
        return self

    def __next__(self):
        if self._expand:
            self._stack.append(iter(_children(self._previous,self._keys)))
        self._expand=False
        while self._stack:
            try:
                node=next(self._stack[-1])
            except StopIteration:
                self._stack.pop()
            else:
                self._previous=node
                self._expand=True
                return node
        raise StopIteration

    def skip(self):
        self._expand=False


def postorder_traversal(node, keys=None):
    stack=[(node,False)]
    while stack:
        value,visited=stack.pop()
        if visited:
            yield value
        else:
            stack.append((value,True))
            stack.extend((child,False) for child in reversed(tuple(_children(value,keys))))


def bottom_up(rv, F, atoms=False, nonbasic=False):
    args=getattr(rv,'args',None)
    if args is None:
        if nonbasic:
            try:
                return F(rv)
            except TypeError:
                pass
        return rv
    if not args:
        return F(rv) if atoms else rv
    changed=tuple(bottom_up(arg,F,atoms,nonbasic) for arg in args)
    if changed!=args:
        rv=rv.func(*changed)
    return F(rv)


def use(expr, func, level=0, args=(), kwargs=None):
    from .. import sympify
    level=operator.index(level)
    if level<0:
        raise ValueError('Transformation level must be nonnegative')
    kwargs={} if kwargs is None else kwargs
    def apply(value,remaining):
        if remaining==0:
            return func(value,*args,**kwargs)
        children=getattr(value,'args',())
        if not children:
            return value
        transformed=tuple(apply(child,remaining-1) for child in children)
        return value if transformed==children else value.func(*transformed)
    return apply(sympify(expr),level)


__all__=['preorder_traversal','postorder_traversal','bottom_up','use']
