"""Optional dependency loading without import-time hard dependencies."""
import importlib
import warnings


def import_module(module,min_module_version=None,min_python_version=None,warn_not_installed=None,
                  warn_old_version=None,module_version_attr='__version__',
                  module_version_attr_call_args=None,import_kwargs=None,catch=()):
    import sys
    if min_python_version is not None and sys.version_info[:len(min_python_version)] < tuple(min_python_version):return None
    try:
        loaded = __import__(module,**import_kwargs) if import_kwargs else importlib.import_module(module)
    except (ImportError,)+tuple(catch):
        if warn_not_installed:warnings.warn('optional dependency is not installed: '+module,RuntimeWarning,stacklevel=2)
        return None
    if min_module_version is not None:
        version = getattr(loaded,module_version_attr,None)
        if module_version_attr_call_args is not None and callable(version):version = version(*module_version_attr_call_args)
        import re
        def version_tuple(text):return tuple(int(n) for n in re.findall(r'\d+',str(text)))
        if version is None or version_tuple(version) < version_tuple(min_module_version):
            if warn_old_version:warnings.warn('optional dependency version is too old: '+module,RuntimeWarning,stacklevel=2)
            return None
    return loaded


__all__ = ['import_module']
