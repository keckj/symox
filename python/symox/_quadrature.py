"""Numerical evaluation of an unevaluated ``Integral``, as SymPy 1.14's
``evalf_integral``.

The native evaluator calls :func:`evalf_integral` for an ``Integral`` with a
single ``(x, a, b)`` limit and no free symbol that ``doit`` leaves
unevaluated (always, for ``quad='osc'``, as SymPy never calls ``doit``). The
quadrature is mpmath's, as in SymPy: ``quadts`` (tanh-sinh, with its error
estimate), or ``quadosc`` for ``quad='osc'``. The integrand is evaluated
natively at each node. The working precision is raised with SymPy's schedule
until the accuracy estimate reaches the target or ``maxprec`` is reached.
"""

import mpmath
from mpmath import mp, workprec
from mpmath.libmp import fzero, from_man_exp

_MINUS_INF = float('-inf')
_INF = float('inf')
# SymPy's default maxn = 100 digits: int(100*LG10) bits.
_DEFAULT_MAXPREC = 332


class _Unavailable(Exception):
    """The integrand or a limit has no numerical value."""


def _fastlog(x):
    if not x or x == fzero:
        return _MINUS_INF
    return x[2] + x[3]


def _signed(value):
    sign, man, exp, _ = value
    return (-man if sign else man), exp


def _point(t):
    if isinstance(t, mpmath.mpc):
        re_man, re_exp = _signed(t.real._mpf_)
        im_man, im_exp = _signed(t.imag._mpf_)
    else:
        re_man, re_exp = _signed(mpmath.mpf(t)._mpf_)
        im_man, im_exp = 0, 0
    return re_man, re_exp, im_man, im_exp


def _parts(value):
    if value is None:
        raise _Unavailable
    re = from_man_exp(value[0], value[1])
    im = from_man_exp(value[2], value[3])
    return (None if re == fzero else re), (None if im == fzero else im)


def _as_mpmath(value, prec):
    from symox._core import _quadrature_point
    from symox import S
    if value == S.Infinity:
        return mpmath.inf
    if value == S.NegativeInfinity:
        return mpmath.ninf
    re, im = _parts(_quadrature_point(value, None, None, prec))
    if im is not None:
        return mp.make_mpc((re or fzero, im))
    return mp.make_mpf(re or fzero)


def _complex_accuracy(result):
    re, im, re_acc, im_acc = result
    if not im:
        if not re:
            return _INF
        return re_acc
    if not re:
        return im_acc
    re_size = _fastlog(re)
    im_size = _fastlog(im)
    absolute_error = max(re_size - re_acc, im_size - im_acc)
    relative_error = absolute_error - max(re_size, im_size)
    return -relative_error


def _do_integral(expr, prec, osc):
    from symox._core import _quadrature_point
    from symox import S, Wild, cos, sin
    func = expr.function
    x, xlow, xhigh = expr.limits[0]
    if xlow == xhigh:
        xlow = xhigh = S.Zero

    with workprec(prec + 5):
        xlow = _as_mpmath(xlow, prec + 15)
        xhigh = _as_mpmath(xhigh, prec + 15)

        have_part = [False, False]
        max_real_term = _MINUS_INF
        max_imag_term = _MINUS_INF

        def f(t):
            nonlocal max_real_term, max_imag_term
            re, im = _parts(_quadrature_point(func, x, _point(t), mp.prec))
            have_part[0] = re or have_part[0]
            have_part[1] = im or have_part[1]
            max_real_term = max(max_real_term, _fastlog(re))
            max_imag_term = max(max_imag_term, _fastlog(im))
            if im:
                return mp.make_mpc((re or fzero, im))
            return mp.make_mpf(re or fzero)

        if osc:
            A = Wild('A', exclude=[x])
            B = Wild('B', exclude=[x])
            D = Wild('D')
            m = func.match(cos(A*x + B)*D)
            if not m:
                m = func.match(sin(A*x + B)*D)
            if not m:
                raise ValueError("An integrand of the form sin(A*x+B)*f(x) "
                  "or cos(A*x+B)*f(x) is required for oscillatory quadrature")
            period = _as_mpmath(2*S.Pi/m[A], prec + 15)
            result = mpmath.quadosc(f, [xlow, xhigh], period=period)
            # mpmath's quadosc has no error estimate.
            quadrature_error = _MINUS_INF
        else:
            result, quadrature_err = mpmath.quadts(f, [xlow, xhigh], error=1)
            quadrature_error = _fastlog(quadrature_err._mpf_)

        # the parts are read at the working precision (a conversion after
        # the block would round them to the default 53 bits)
        result_re = result.real._mpf_
        result_im = result.imag._mpf_
    if have_part[0]:
        re = result_re
        if re == fzero:
            # SymPy returns a scaled zero without a correct bit.
            re, re_acc = fzero, -1
        else:
            re_acc = int(-max(max_real_term - _fastlog(re) - prec, quadrature_error))
    else:
        re, re_acc = None, None
    if have_part[1]:
        im = result_im
        if im == fzero:
            im, im_acc = fzero, -1
        else:
            im_acc = int(-max(max_imag_term - _fastlog(im) - prec, quadrature_error))
    else:
        im, im_acc = None, None
    return re, im, re_acc, im_acc


def evalf_integral(expr, prec, osc=False, maxprec=None):
    """SymPy's ``evalf_integral`` at ``prec`` bits.

    Returns ``(re_man, re_exp, im_man, im_exp)`` when the accuracy estimate
    reaches ``prec`` bits, None otherwise (the integral then stays
    unevaluated, where SymPy returns a Float of the lower accuracy). The
    errors of mpmath's quadrature propagate, as in SymPy.
    """
    limits = expr.limits
    if len(limits) != 1 or len(limits[0]) != 3:
        return None
    maxprec = max(prec, _DEFAULT_MAXPREC if maxprec is None else maxprec)
    workprec_ = prec
    i = 0
    try:
        while 1:
            result = _do_integral(expr, workprec_, osc)
            accuracy = _complex_accuracy(result)
            if accuracy >= prec:
                break
            if workprec_ >= maxprec:
                return None
            if accuracy == -1:
                workprec_ *= 2
            else:
                workprec_ += max(prec, 2**i)
            workprec_ = min(workprec_, maxprec)
            i += 1
    except _Unavailable:
        return None
    re, im, _, _ = result
    if re is None and im is None:
        return 0, 0, 0, 0
    re_man, re_exp = _signed(re) if re is not None else (0, 0)
    im_man, im_exp = _signed(im) if im is not None else (0, 0)
    return re_man, re_exp, im_man, im_exp
