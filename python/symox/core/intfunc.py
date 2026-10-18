"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox.ntheory._intfunc', 'igcd igcdex ilcm integer_log integer_nthroot isqrt mod_inverse num_digits trailing')


def igcd_lehmer(a,b):
    """Exact integer gcd using Python's multiple-precision implementation."""
    import math,operator
    return math.gcd(operator.index(a),operator.index(b))
__all__.append('igcd_lehmer')
