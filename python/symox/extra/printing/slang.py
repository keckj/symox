"""Standalone scalar Slang expressions, with optional custom printer hooks.

These expressions use only Slang builtins. ``float32`` is the default;
``float64`` emits double literals, whose availability depends on the target.
All inputs are real scalars. Batched kernels apply this scalar code per item.
"""

from .pytorch import PyTorchPrinter


def _is_slang_identifier(name):
    from symox import _core
    return isinstance(name, str) and _core._is_slang_identifier(name)


def _normalize_slang_dtype(dtype):
    aliases = {"float": "float32", "float32": "float32", "double": "float64", "float64": "float64"}
    try:
        return aliases[dtype]
    except (KeyError, TypeError):
        raise ValueError("Slang dtype must be float32 or float64") from None


class SlangCodePrinter(PyTorchPrinter):
    """Native real scalar Slang printer, supporting nested ``_print_*`` hooks.

    ``symbol_names`` maps symbols (including their assumptions) or textual
    names to safe Slang identifiers. ``user_functions`` maps symbolic function
    names to caller-provided Slang functions, optionally namespace-qualified.
    Complex and unsupported special functions raise ``NotImplementedError``.
    """

    language = "Slang"
    printmethod = "_slangcode"
    _native_method = "_slangcode"
    _default_settings = {"order": None, "dtype": "float32", "symbol_names": {}, "user_functions": {}}

    _normalize_dtype = staticmethod(_normalize_slang_dtype)

    def __init__(self, settings=None, **kwargs):
        super().__init__(settings, **kwargs)
        self._settings["dtype"] = type(self)._normalize_dtype(self._settings["dtype"])
        self._custom_methods = {
            name: getattr(type(self), name)
            for name in dir(type(self))
            if name.startswith("_print_")
            and getattr(type(self), name) is not getattr(SlangCodePrinter, name, None)
        }

    def doprint(self, expr, assign_to=None):
        from symox import sympify
        result = self._print(sympify(expr))
        if assign_to is not None:
            name = assign_to if isinstance(assign_to, str) else getattr(assign_to, "name", "")
            if not self._is_identifier(name):
                raise ValueError("assign_to must be a safe Slang identifier")
            result = f"{name} = {result};"
        return result

    def _is_identifier(self, name):
        return _is_slang_identifier(name)

    def emptyPrinter(self, expr):
        from symox import _core
        overrides = {}
        if self._custom_methods:
            stack, seen = [expr], set()
            while stack:
                item = stack.pop()
                if item in seen:
                    continue
                seen.add(item)
                if item not in self._active_overrides and self._custom_method(item) is not None:
                    overrides[item] = self._print(item)
                else:
                    stack.extend(getattr(item, "args", ()))
        return getattr(_core, self._native_method)(
            expr,
            dtype=self._settings["dtype"],
            symbol_names=self._settings["symbol_names"],
            user_functions=self._settings["user_functions"],
            expression_overrides=overrides or None,
        )

    _print_Symbol = _print_Integer = _print_Rational = _print_Float = emptyPrinter
    _print_Add = _print_Mul = _print_Pow = _print_Piecewise = emptyPrinter
    _print_UnevaluatedExpr = emptyPrinter


SlangPrinter = SlangCodePrinter


def slangcode(expr, assign_to=None, **settings):
    """Print a real scalar expression using portable Slang builtin functions."""
    return SlangCodePrinter(settings).doprint(expr, assign_to)


__all__ = ["SlangCodePrinter", "SlangPrinter", "slangcode"]
