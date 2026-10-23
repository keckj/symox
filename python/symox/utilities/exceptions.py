"""Visible deprecation warnings and scoped warning filters (SymPy 1.14 port)."""
import contextlib
from textwrap import dedent
import warnings


class SymPyDeprecationWarning(DeprecationWarning):
    """A warning for deprecated features of SymPy."""

    def __init__(self, message, *, deprecated_since_version, active_deprecations_target):
        super().__init__(message, deprecated_since_version, active_deprecations_target)
        self.message = message
        if not isinstance(deprecated_since_version, str):
            raise TypeError(f"'deprecated_since_version' should be a string, got {deprecated_since_version!r}")
        self.deprecated_since_version = deprecated_since_version
        self.active_deprecations_target = active_deprecations_target
        if any(i in active_deprecations_target for i in '()='):
            raise ValueError("active_deprecations_target be the part inside of the '(...)='")
        self.full_message = f"""

{dedent(message).strip()}

See https://docs.sympy.org/latest/explanation/active-deprecations.html#{active_deprecations_target}
for details.

This has been deprecated since SymPy version {deprecated_since_version}. It
will be removed in a future version of SymPy.
"""

    def __str__(self):
        return self.full_message

    def __repr__(self):
        return (f"{self.__class__.__name__}({self.message!r}, "
                f"deprecated_since_version={self.deprecated_since_version!r}, "
                f"active_deprecations_target={self.active_deprecations_target!r})")

    def __eq__(self, other):
        return isinstance(other, SymPyDeprecationWarning) and self.args == other.args

    @classmethod
    def _new(cls, message, deprecated_since_version, active_deprecations_target):
        return cls(message, deprecated_since_version=deprecated_since_version,
                   active_deprecations_target=active_deprecations_target)

    def __reduce__(self):
        return (self._new, (self.message, self.deprecated_since_version,
                            self.active_deprecations_target))


# As in SymPy, every deprecation is shown once, wherever it is triggered.
warnings.simplefilter("once", SymPyDeprecationWarning)


def sympy_deprecation_warning(message, *, deprecated_since_version,
                              active_deprecations_target, stacklevel=3):
    """Warn that a feature is deprecated (see SymPyDeprecationWarning)."""
    w = SymPyDeprecationWarning(message, deprecated_since_version=deprecated_since_version,
                                active_deprecations_target=active_deprecations_target)
    warnings.warn(w, stacklevel=stacklevel)


@contextlib.contextmanager
def ignore_warnings(warningcls):
    """Context manager to suppress warnings of class *warningcls*; other
    warnings are re-emitted."""
    with warnings.catch_warnings(record=True) as warnrec:
        warnings.simplefilter("always", warningcls)
        yield
    for w in warnrec:
        if not issubclass(w.category, warningcls):
            warnings.warn_explicit(w.message, w.category, w.filename, w.lineno)


__all__ = ['SymPyDeprecationWarning', 'sympy_deprecation_warning', 'ignore_warnings']
