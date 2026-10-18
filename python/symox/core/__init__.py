"""Core symbolic types and evaluation controls."""
from .. import (Basic,Expr,Symbol,Dummy,Wild,WildFunction,Number,Integer,Rational,
    Float,RealNumber,NumberSymbol,Add,Mul,Pow,Tuple,Dict,Function,FunctionClass,
    Derivative,Lambda,Subs,sympify,SympifyError,evaluate,expand,expand_multinomial,
    arity,comp,seterr)
__all__ = [name for name in dir() if not name.startswith('_')]

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'AlgebraicNumber Atom AtomicExpr BooleanKind Catalan E Eq Equality EulerGamma Ge GoldenRatio GreaterThan Gt I Le LessThan Lt Mod N Ne NumberKind PoleError PrecisionExhausted Rel S StrictGreaterThan StrictLessThan TribonacciConstant UndefinedKind Unequality UnevaluatedExpr bottom_up cacheit check_assumptions count_ops default_sort_key diff evalf expand_complex expand_func expand_log expand_mul expand_power_base expand_power_exp expand_trig factor_nc factor_terms failing_assumptions gcd_terms igcd ilcm integer_log integer_nthroot mod_inverse nan nfloat num_digits oo ordered pi postorder_traversal preorder_traversal prod symbols trailing use var vectorize zoo')

from importlib import import_module as _import_module
# Complete homonymous child imports before restoring callable package exports.
_import_module('.evalf', __name__)
_import_module('.sympify', __name__)
from symox import sympify
from .assumptions import assumptions, common_assumptions
__all__ += ['assumptions', 'common_assumptions']
