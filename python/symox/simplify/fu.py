"""Public Fu rules backed by the existing native transformation registry."""
from .._simptools import fu_rules as _rules
__all__ = ['L', 'TR0', 'TR1', 'TR10', 'TR10i', 'TR11', 'TR111', 'TR12', 'TR12i', 'TR13', 'TR14', 'TR15', 'TR16', 'TR2', 'TR22', 'TR2i', 'TR3', 'TR4', 'TR5', 'TR6', 'TR7', 'TR8', 'TR9', 'TRmorrie', 'TRpower', 'as_f_sign_1', 'fu', 'hyper_as_trig', 'process_common_addends', 'sincos_to_sum', 'trig_split']

# Names the reference module imports into its namespace (``from
# sympy.simplify.fu import TR1, sec, csc`` in the TR1 docstring).
_REFERENCE_NAMES = frozenset(('Add', 'Dummy', 'Expr', 'Factors', 'HyperbolicFunction',
    'I', 'Mul', 'Pow', 'S', 'TrigonometricFunction', 'binomial', 'bottom_up', 'cacheit',
    'cos', 'cosh', 'cot', 'coth', 'csc', 'csch', 'expand_mul', 'factor', 'factor_terms',
    'gcd_terms', 'identity', 'ordered', 'perfect_power', 'pi', 'sec', 'sech',
    'sin', 'sinh', 'sqrt', 'sympify', 'tan', 'tanh'))

# Private helpers that the reference tests import from the module.
_PRIVATE_RULES = frozenset(('_TR11', '_TR56', '_osborne', '_osbornei'))

def __getattr__(name):
    if name in __all__ or name in _PRIVATE_RULES:return getattr(_rules, name)
    if name in _REFERENCE_NAMES:
        import symox
        for source in ('symox', 'symox.simplify.simplify', 'symox.core.exprtools',
                       'symox.strategies', 'symox.strategies.core', 'symox.utilities.iterables',
                       'symox.ntheory', 'symox.core.cache', 'symox.core.function',
                       'symox.functions.elementary.trigonometric',
                       'symox.functions.elementary.hyperbolic', 'symox.strategies.tree'):
            try:
                module = __import__(source, fromlist=['_'])
            except ImportError:
                continue
            if hasattr(module, name):
                return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():return sorted(set(globals()) | set(__all__))
