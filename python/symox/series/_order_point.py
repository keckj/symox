"""Moving point of ``Order._eval_subs`` for a non-affine replacement."""


def order_point(old, new, var, point):
    """SymPy 1.14 ``Order._eval_subs``: solve ``old = new(d)`` for ``d``
    with ``solveset``, take the first solution and its limit at
    ``old -> point``. ``None`` keeps the substitution in the arguments."""
    from symox import Dummy, solveset, Complement, ImageSet, S, limit
    d = Dummy()
    try:
        sol = solveset(old - new.subs(var, d), d)
        if isinstance(sol, Complement):
            sol = set(sol.args[0]) - set(sol.args[1])
        if isinstance(sol, ImageSet) and sol.base_set == S.Integers:
            # ImageSet.__iter__ starts with the image of 0 for Integers
            sol = [sol.lamda(S.Zero)]
        res = dict(zip((d,), sol))
        return limit(d.subs(res), old, point)
    except (TypeError, ValueError, NotImplementedError):
        return None
