"""Exact limits using the Rust most-rapidly-varying comparison engine."""


def gruntz(e,z,z0,dir='+'):
    from .. import _core,PoleError
    result=_core._limit(e,z,z0,str(dir))
    if result is None:
        raise PoleError('The native limit engine could not determine this limit')
    return result


__all__=['gruntz']
