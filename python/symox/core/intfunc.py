"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox.ntheory._intfunc', 'igcd igcdex ilcm integer_log integer_nthroot isqrt mod_inverse num_digits trailing')
