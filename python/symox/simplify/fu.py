"""Public Fu rules backed by the existing native transformation registry."""
from .._simptools import fu_rules as _rules
__all__ = ['L', 'TR0', 'TR1', 'TR10', 'TR10i', 'TR11', 'TR111', 'TR12', 'TR12i', 'TR13', 'TR14', 'TR15', 'TR16', 'TR2', 'TR22', 'TR2i', 'TR3', 'TR4', 'TR5', 'TR6', 'TR7', 'TR8', 'TR9', 'TRmorrie', 'TRpower', 'as_f_sign_1', 'fu', 'hyper_as_trig', 'process_common_addends', 'sincos_to_sum', 'trig_split']

def __getattr__(name):
    if name in __all__:return getattr(_rules, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():return sorted(set(globals()) | set(__all__))
