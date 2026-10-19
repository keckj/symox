"""Container protocols for Tuple nodes sharing the native scalar DAG type."""
import operator


def _iter(self):
    if not self.is_Tuple:
        raise TypeError(f"'{self.class_name}' object is not iterable")
    return iter(self.args)


def _len(self):
    if not self.is_Tuple:
        raise TypeError(f"object of type '{self.class_name}' has no len()")
    return len(self.args)


def _getitem(self, index):
    if not self.is_Tuple:
        raise TypeError(f"'{self.class_name}' object is not subscriptable")
    if isinstance(index, slice):
        from . import Tuple
        return Tuple(*self.args[index])
    return self.args[operator.index(index)]


def install(expr_type):
    expr_type.__iter__ = _iter
    expr_type.__len__ = _len
    expr_type.__getitem__ = _getitem
    expr_type._iterable = property(lambda self: self.is_Tuple)
