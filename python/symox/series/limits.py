"""Compatible imports for the native exact limit engine."""
from .. import Limit,limit


def heuristics(e, z, z0, dir):
    """``sympy.series.limits.heuristics``: computed natively (``None`` when undecided)."""
    from .. import _core, sympify
    return _core._limit_heuristics(sympify(e), sympify(z), sympify(z0), str(dir))


__all__=['Limit','limit','heuristics']
