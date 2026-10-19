"""Container protocols for Tuple nodes sharing the native scalar DAG type."""
import operator


def _iter(self):
    if not self.is_Tuple:
        raise TypeError(f"'{self.class_name}' object is not iterable")
    return iter(self.args)


def _len(self):
    if not self.is_Tuple:
        raise TypeError(f"object of type '{self.class_name}' has no len()")
    return self._tuple_len()


def _getitem(self, index):
    if not self.is_Tuple:
        raise TypeError(f"'{self.class_name}' object is not subscriptable")
    if isinstance(index, slice):
        from . import Tuple
        return Tuple(*self.args[index])
    return self._tuple_item(operator.index(index))


def install(expr_type):
    from . import Tuple
    # Sequence slots belong to the actual Tuple subtype. Scalar Expr objects
    # must neither satisfy Iterable nor create a sequence iterator.
    Tuple.__iter__ = _iter
    Tuple.__len__ = _len
    Tuple.__getitem__ = _getitem
    expr_type._iterable = property(lambda self: self.is_Tuple)
