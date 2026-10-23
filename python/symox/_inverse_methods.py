"""SymPy's ``Function.inverse(self, argindex=1)`` methods: the head of the
inverse function, as defined by each SymPy 1.14 class (``exp`` -> ``log``,
``asin`` -> ``sin``, ``erf`` -> ``erfinv``, ...). Heads that SymPy does not
give an ``inverse`` (``sin``, ``cos``, ...) keep raising AttributeError."""

_INVERSES = {
    'acos': 'cos', 'acosh': 'cosh', 'acot': 'cot', 'acoth': 'coth', 'acsc': 'csc',
    'acsch': 'csch', 'asec': 'sec', 'asech': 'sech', 'asin': 'sin', 'asinh': 'sinh',
    'atan': 'tan', 'atanh': 'tanh', 'conjugate': 'conjugate', 'cot': 'acot',
    'coth': 'acoth', 'erf': 'erfinv', 'erfc': 'erfcinv', 'erfcinv': 'erfc',
    'erfinv': 'erf', 'exp': 'log', 'exp_polar': 'log', 'log': 'exp',
    'sinh': 'asinh', 'tan': 'atan', 'tanh': 'atanh',
}


def _make_inverse(target):
    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        import symox
        return getattr(symox, target)
    return inverse


def install(namespace):
    for name, target in _INVERSES.items():
        head = namespace.get(name)
        if head is None or 'inverse' in getattr(head, '__dict__', {}):
            continue
        try:
            head.inverse = _make_inverse(target)
        except (AttributeError, TypeError):
            pass


# ``f._eval_rewrite_as_exp(*args)``: the private rewrite hook SymPy's tests
# call directly (``cos(x)._eval_rewrite_as_exp(x)``); it is the exp rewrite
# of the function applied to ``args``.
_EXP_REWRITES = ('sin', 'cos', 'tan', 'cot', 'sec', 'csc',
                 'sinh', 'cosh', 'tanh', 'coth', 'sech', 'csch')


def _eval_rewrite_as_exp(self, *args, **kwargs):
    import symox
    return self.func(*args).rewrite(symox.exp)


def install_rewrites(namespace):
    for name in _EXP_REWRITES:
        head = namespace.get(name)
        if head is None or '_eval_rewrite_as_exp' in getattr(head, '__dict__', {}):
            continue
        try:
            head._eval_rewrite_as_exp = _eval_rewrite_as_exp
        except (AttributeError, TypeError):
            pass
