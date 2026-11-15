"""Exact limits using the Rust most-rapidly-varying comparison engine."""


def gruntz(e,z,z0,dir='+'):
    from .. import _core,PoleError,sympify
    z=sympify(z)
    if not z.is_symbol:
        raise NotImplementedError("Second argument must be a Symbol")
    result=_core._limit(e,z,z0,str(dir))
    if result is None:
        # SymPy's mrv_leadterm gives the accumulation bounds of an
        # oscillating leading coefficient (`gruntz(sin(x), x, oo)` is
        # `AccumBounds(-1, 1)`), as `limit` computes them.
        from ..calculus._accum_limit import oscillating_limit
        result=oscillating_limit(e,z,z0,str(dir))
    if result is None:
        raise PoleError('The native limit engine could not determine this limit')
    return result


__all__=['gruntz']
