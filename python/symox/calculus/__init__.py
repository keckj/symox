"""Native symbolic calculus adapters."""
from .finite_diff import finite_diff_weights, apply_finite_diff, differentiate_finite
from .singularities import (singularities, is_increasing, is_strictly_increasing,
    is_decreasing, is_strictly_decreasing, is_monotonic)
from .util import continuous_domain, stationary_points, periodicity, is_convex, minimum, maximum, function_range, not_empty_in

__all__=['finite_diff_weights','apply_finite_diff','differentiate_finite',
    'singularities','continuous_domain','stationary_points','periodicity',
    'is_increasing','is_strictly_increasing','is_decreasing',
    'is_strictly_decreasing','is_monotonic','is_convex','minimum','maximum','function_range','not_empty_in']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'AccumBounds euler_equations')
