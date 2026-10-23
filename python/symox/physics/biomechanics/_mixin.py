"""Mixin classes for the biomechanics package, ported from SymPy 1.14."""
__all__ = ['_NamedMixin']

class _NamedMixin:

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        if hasattr(self, '_name'):
            msg = f"Can't set attribute `name` to {repr(name)} as it is immutable."
            raise AttributeError(msg)
        if not isinstance(name, str):
            msg = f'Name {repr(name)} passed to `name` was of type {type(name)}, must be {str}.'
            raise TypeError(msg)
        if name in {''}:
            msg = f'Name {repr(name)} is invalid, must be a nonzero length {type(str)}.'
            raise ValueError(msg)
        self._name = name
