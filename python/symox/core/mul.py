"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'Mul prod')

class NC_Marker:
    """Private sentinel used to mark a noncommutative factor group."""

def expand_2arg(e):
    from symox import Mul,Add
    if not getattr(e,'is_Mul',False) or len(e.args)!=2:return e
    a,b=e.args
    if a.is_Add:return Add(*(term*b for term in a.args))
    if b.is_Add:return Add(*(a*term for term in b.args))
    return e
__all__ += ['NC_Marker','expand_2arg']


def _mulsort(args):
    # in-place sorting of args
    from .basic import _args_sortkey
    args.sort(key=_args_sortkey)


def _unevaluated_Mul(*args):
    """Return a well-formed unevaluated Mul: Numbers are collected and
    put in slot 0, any arguments that are Muls will be flattened, and args
    are sorted. Use this when args have changed but you still want to return
    an unevaluated Mul.
    """
    from symox import Mul, S
    cargs = []
    ncargs = []
    args = list(args)
    co = S.One
    for a in args:
        if a.is_Mul:
            a_c, a_nc = a.args_cnc()
            args.extend(a_c)  # grow args
            ncargs.extend(a_nc)
        elif a.is_Number:
            co *= a
        elif a.is_commutative:
            cargs.append(a)
        else:
            ncargs.append(a)
    _mulsort(cargs)
    if co is not S.One:
        cargs.insert(0, co)
    return Mul._from_args(cargs+ncargs)


def _keep_coeff(coeff, factors, clear=True, sign=False):
    """Return ``coeff*factors`` unevaluated if necessary.

    If ``clear`` is False, do not keep the coefficient as a factor
    if it can be distributed on a single factor such that one or
    more terms will still have integer coefficients.

    If ``sign`` is True, allow a coefficient of -1 to remain factored out.
    """
    from symox import Add, Mul, S
    if not coeff.is_Number:
        if factors.is_Number:
            factors, coeff = coeff, factors
        else:
            return coeff*factors
    if factors is S.One:
        return coeff
    if coeff is S.One:
        return factors
    elif coeff is S.NegativeOne and not sign:
        return -factors
    elif factors.is_Add:
        if not clear and coeff.is_Rational and coeff.q != 1:
            args = [i.as_coeff_Mul() for i in factors.args]
            args = [(_keep_coeff(c, coeff), m) for c, m in args]
            if any(c.is_Integer for c, _ in args):
                return Add._from_args([Mul._from_args(
                    i[1:] if i[0] == 1 else i) for i in args])
        return Mul(coeff, factors, evaluate=False)
    elif factors.is_Mul:
        margs = list(factors.args)
        if margs[0].is_Number:
            margs[0] *= coeff
            if margs[0] == 1:
                margs.pop(0)
        else:
            margs.insert(0, coeff)
        return Mul._from_args(margs)
    else:
        m = coeff*factors
        if m.is_Number and not factors.is_Number:
            m = Mul._from_args((coeff, factors))
        return m
__all__ += ['_mulsort', '_unevaluated_Mul', '_keep_coeff']
