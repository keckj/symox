"""SymPy's per-object fact cache ``_assumptions``.

SymPy stores the facts known about an object in ``obj._assumptions``, a
``StdFactKB``: for a Symbol the facts deduced from its assumptions, for other
objects the facts computed so far. Native nodes compute their facts in the
kernel, so the cache is read as a snapshot: a Symbol's deduced assumptions,
and no cached facts for other nodes. Instances of user classes may assign
their own cache (Devito's symbols do), which then shadows the snapshot.
"""
from . import _core


class _FactCache:
    """Non-data descriptor: an instance attribute takes precedence."""

    def __get__(self, obj, cls=None):
        if obj is None:
            return self
        from .core.assumptions import StdFactKB
        if getattr(obj, 'is_Symbol', False):
            return StdFactKB(obj.assumptions0)
        return StdFactKB({})


def install():
    if '_assumptions' not in _core.Expr.__dict__:
        _core.Expr._assumptions = _FactCache()
