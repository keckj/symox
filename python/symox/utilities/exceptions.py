"""Visible deprecation warnings and scoped warning filters."""
from contextlib import contextmanager
import warnings
class SymPyDeprecationWarning(DeprecationWarning):
    def __init__(self,message,*,deprecated_since_version,active_deprecations_target):
        self.message=str(message);self.deprecated_since_version=deprecated_since_version;self.active_deprecations_target=active_deprecations_target
        super().__init__(self.message+'\nDeprecated since '+str(deprecated_since_version)+'. See '+str(active_deprecations_target)+'.')
def sympy_deprecation_warning(message,*,deprecated_since_version,active_deprecations_target,stacklevel=3):
    warnings.warn(SymPyDeprecationWarning(message,deprecated_since_version=deprecated_since_version,active_deprecations_target=active_deprecations_target),stacklevel=stacklevel)
@contextmanager
def ignore_warnings(warningcls):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',warningcls);yield
__all__=['SymPyDeprecationWarning','sympy_deprecation_warning','ignore_warnings']
