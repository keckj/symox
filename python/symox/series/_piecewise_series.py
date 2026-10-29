"""Series of expressions holding a ``Piecewise``.

Ports of the pieces of SymPy 1.14 that ``Expr.series`` reaches when an
expression contains a ``Piecewise``: ``Piecewise._eval_nseries`` (each
branch with ``_eval_nseries`` semantics), ``Add._eval_nseries`` (term by
term), ``Mul._eval_nseries`` (leading terms of the factors, where a
``Piecewise`` contributes the leading term of its first branch whose
condition holds at 0, ``Piecewise._eval_as_leading_term``) and the order
term handling of ``Expr.series``. Expressions without a ``Piecewise`` keep
the native series.
"""


def _has_pw(e, x):
    from symox import Piecewise
    return any(x in p.free_symbols for p in e.atoms(Piecewise))


def _native(e, x, n):
    from symox import _series_native
    return _series_native(e, x, 0, n, "+")


def _top_orders(s):
    from symox import Add
    return [t for t in Add.make_args(s) if t.is_Order]


def _remove_top_o(s):
    from symox import Add
    return Add(*[t for t in Add.make_args(s) if not t.is_Order])


def _getn(s, x):
    """``Expr.getn``: the exponent of a top-level ``O(x**k)``, else None."""
    from symox import Mul
    os_ = _top_orders(s)
    if not os_:
        return None
    o = os_[0].expr
    if o.is_Symbol:
        return 1
    if o.is_Pow and o.base == x and o.exp.is_Rational:
        return o.exp
    if o.is_Mul:
        for t in Mul.make_args(o):
            if t.is_Symbol:
                return 1
            if t.is_Pow and t.base.is_Symbol and t.exp.is_Rational:
                return t.exp
    return None


def nseries(e, x, n):
    """``e.nseries(x, n=n)`` at 0 from the right (``_eval_nseries``)."""
    if x not in e.free_symbols:
        return e
    return eval_nseries(e, x, n)


def eval_nseries(e, x, n):
    from symox import Add, Order
    if e == x:
        return e
    if e.is_Piecewise:
        return piecewise_nseries(e, x, n)
    if e.is_Add:
        return Add(*[nseries(t, x, n) for t in e.args])
    if e.is_Mul and _has_pw(e, x):
        return mul_nseries(e, x, n)
    if e.is_Pow and e.base == x and e.exp.is_number and e.exp.is_real:
        # Pow._eval_nseries: maxpow = n - m*e < 0 gives O(x**(m*e)),
        # otherwise x**e is its own expansion
        if (n - e.exp).is_negative:
            return Order(e, x)
        return e
    return _native(e, x, n)


def piecewise_nseries(e, x, n):
    """``Piecewise._eval_nseries``: expand each branch."""
    return e.func(*[(nseries(b, x, n), c) for b, c in e.args])


def _as_leading_term(t, x):
    from symox import true
    if t.is_Piecewise:
        # Piecewise._eval_as_leading_term
        for b, c in t.args:
            if c == true or c.subs(x, 0) == true:
                return _as_leading_term(b, x)
        raise NotImplementedError("as_leading_term(%s, %s)" % (t, x))
    if _has_pw(t, x):
        raise NotImplementedError("as_leading_term(%s, %s)" % (t, x))
    if t.is_Add and _top_orders(t):
        # Expr.as_leading_term works on removeO() when it is nonzero
        rest = _remove_top_o(t)
        if rest != 0:
            t = rest
    return t.as_leading_term(x)


def _leadterm(t, x):
    from symox import Dummy, log
    lt = _as_leading_term(t, x)
    d = Dummy("logx")
    if lt.has(log(x)):
        lt = lt.subs(log(x), d)
    c, ex = lt.as_coeff_exponent(x)
    if x in c.free_symbols:
        raise ValueError("cannot compute leadterm(%s, %s)" % (t, x))
    return c.subs(d, log(x)), ex


def _coeff_exp(term, x):
    from symox import S
    if term.is_Piecewise:
        lt = (term, S.Zero)
    else:
        lt = term.as_coeff_exponent(x)
    if lt[0].has(x):
        try:
            lt = _leadterm(term, x)
        except (ValueError, NotImplementedError):
            return term, S.Zero
    return lt


def mul_nseries(e, x, n):
    """``Mul._eval_nseries`` for a product holding a ``Piecewise``."""
    from itertools import product
    from symox import Add, Mul, Order, S, ceiling, powsimp, expand_mul, sympify
    from symox import PoleError
    ords = []
    try:
        for t in e.args:
            coeff, ex = _leadterm(t, x)
            if coeff.has(x):
                raise ValueError
            ords.append((t, ex))
        n0 = sum(m for _, m in ords if m.is_number)
        facs = []
        for t, m in ords:
            n1 = int(ceiling(n - n0 + (m if m.is_number else 0)))
            s = nseries(t, x, n1)
            ns = _getn(s, x)
            if ns is not None and ns < n1:
                n -= n1 - ns
            facs.append(s)
    except (ValueError, NotImplementedError, TypeError, PoleError):
        n0 = sympify(sum(m for _, m in ords if m.is_number))
        if n0.is_nonnegative:
            n0 = S.Zero
        facs = [nseries(t, x, int(ceiling(n - n0))) for t in e.args]
        res = powsimp(Mul(*facs).expand(), combine="exp", deep=True)
        if res.has(Order):
            res += Order(x**n, x)
        return res
    res = S.Zero
    if any(t.is_Mul for t in facs):
        facs = [t.expand() for t in facs]
    for fac in product(*[Add.make_args(f) for f in facs]):
        if any(t.is_Order for t in fac):
            continue
        coeffs, powers = zip(*[_coeff_exp(t, x) for t in fac])
        power = Add(*powers)
        if (power - n).is_negative:
            res += Mul(*coeffs)*x**power
    res = expand_mul(res)
    if res != e:
        try:
            zero = (e - res).subs(x, 0) == S.Zero
        except (ValueError, NotImplementedError, TypeError):
            zero = False
        if zero and n > 0:
            try:
                lt = _as_leading_term(e, x)
            except NotImplementedError:
                lt = None
            if lt == S.Zero:
                return res
        res += Order(x**n, x)
    return res


def series(expr, x, x0, n, dir):
    """``Expr.series`` for an expression holding a ``Piecewise`` in ``x``."""
    from symox import Dummy, Order, Rational, S, collect, sympify
    x0 = sympify(x0)
    if x0 in (S.Infinity, S.NegativeInfinity):
        # Expr.series at an infinity: the expansion of e(sgn/x) at 0+,
        # substituted back (the conditions stay as they are when they do
        # not decide a branch, Piecewise((x, y < 1), (x**5, True)) at oo).
        from symox import series as _series
        sgn = 1 if x0 is S.Infinity else -1
        return _series(expr.subs(x, sgn/x), x, 0, n, "+").subs(x, sgn/x)
    if x0.is_infinite:
        return None
    if x0 or dir == "-":
        if dir == "-":
            rep, back = -x + x0, -x + x0
        else:
            rep, back = x + x0, x - x0
        from symox import series as _series
        s = _series(expr.subs(x, rep), x, 0, n, "+")
        return s.subs(x, back)
    if (x.is_positive is None and x.is_negative is None) or not x.is_Symbol:
        xpos = Dummy("x", positive=True)
        from symox import series as _series
        return _series(expr.subs(x, xpos), xpos, 0, n, "+").subs(xpos, x)
    s1 = eval_nseries(expr, x, n)
    o = _top_orders(s1)
    if o:
        o = o[0]
        ngot = _getn(s1, x)
        if ngot is not None and ngot > n:
            if n != 0:
                s1 += o.subs(x, x**Rational(n, ngot))
            else:
                s1 += Order(1, x)
        elif ngot is not None and ngot < n:
            from symox import ceiling
            for more in range(1, 9):
                s1 = eval_nseries(expr, x, n + more)
                newn = _getn(s1, x)
                if newn != ngot:
                    ndo = n + int(ceiling((n - ngot)*more/(newn - ngot)))
                    s1 = eval_nseries(expr, x, ndo)
                    while _getn(s1, x) < n:
                        s1 = eval_nseries(expr, x, ndo)
                        ndo += 1
                    break
            else:
                raise ValueError("Could not calculate %s terms for %s"
                                 % (str(n), expr))
            s1 += Order(x**n, x)
        o = _top_orders(s1)
        o = o[0] if o else 0
        s1 = _remove_top_o(s1)
    elif s1.has(Order):
        # asymptotic expansion
        return s1
    else:
        o = Order(x**n, x)
        s1done = s1.doit()
        try:
            if _remove_top_o(s1done + o) == s1done:
                o = 0
        except NotImplementedError:
            return s1
    try:
        return collect(s1, x) + o
    except NotImplementedError:
        return s1 + o
