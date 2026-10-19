"""Native PyTorch expression printing.

``torchcode`` emits tensor expressions which use ``torch`` and the helper
``_sx_tensor(value, reference)``. Use :func:`symox.extra.torch_lambdify` or
:func:`symox.extra.torchcode_function` for complete, executable functions with that
helper, device/dtype inference, broadcasting, and autograd support.

The helper preserves existing tensors and tensorizes scalar values using the
optional reference tensor's device and dtype (booleans retain their dtype).
An expression printed without ``reference=`` uses ``None`` as its reference.
"""

from ...printing.printer import Printer


class PyTorchPrinter(Printer):
    """PyTorch printer backed by the native expression tree.

    Settings: ``symbol_names``, ``user_functions``, ``reference``, ``mask``.
    ``loggamma`` requires a provably positive real argument, since the
    real ``lgamma`` primitive computes log(abs(gamma)) on negative inputs.
    Explicit function mappings are dotted Python names; unknown functions
    raise ``NotImplementedError``. ``_print_<Class>`` subclass hooks and the
    expression ``_torchcode`` hook work through the shared printer protocol.
    """

    language = "Python with PyTorch"
    printmethod = "_torchcode"
    _native_method = "_torchcode"
    _default_settings = {
        "order": None,
        "symbol_names": {},
        "user_functions": {},
        "reference": None,
        "mask": None,
    }

    def __init__(self, settings=None, **kwargs):
        super().__init__(settings, **kwargs)
        allowed = set(self._default_settings) | {"fully_qualified_modules", "strict"}
        unknown = set(self._settings) - allowed
        if unknown:
            raise TypeError(f"Unknown {self.language} printer setting: {sorted(unknown)[0]}")
        self._active_overrides = set()
        self._custom_methods = {
            name: getattr(type(self), name)
            for name in dir(type(self))
            if name.startswith("_print_")
            and getattr(type(self), name) is not getattr(PyTorchPrinter, name, None)
        }

    def _print(self, expr, **kwargs):
        # Mark the dispatched node while hooks call super(), so a native
        # fallback does not redispatch the same custom hook recursively.
        was_active = expr in self._active_overrides
        self._active_overrides.add(expr)
        try:
            return super()._print(expr, **kwargs)
        finally:
            if not was_active:
                self._active_overrides.remove(expr)

    def doprint(self, expr, assign_to=None):
        from symox import sympify
        result = self._print(sympify(expr))
        if assign_to is not None:
            import keyword
            name = assign_to if isinstance(assign_to, str) else getattr(assign_to, "name", "")
            if not name.isidentifier() or keyword.iskeyword(name):
                raise ValueError("assign_to must be a valid Python identifier")
            result = f"{name} = {result}"
        return result

    def _custom_method(self, expr):
        name = getattr(expr, "class_name", type(expr).__name__)
        method = self._custom_methods.get("_print_" + name)
        if method is None and name in ("Zero", "One", "NegativeOne", "Half"):
            method = self._custom_methods.get("_print_" + ("Rational" if name == "Half" else "Integer"))
        if method is None and name in ("Dummy", "Wild"):
            method = self._custom_methods.get("_print_Symbol")
        if method is None:
            for cls in type(expr).__mro__:
                method = self._custom_methods.get("_print_" + cls.__name__)
                if method is not None:
                    break
        return method

    def emptyPrinter(self, expr):
        from symox import _core
        overrides = {}
        if self._native_method == "_torchcode":
            from symox import I
            reference = self._settings["reference"] or "None"
            # Keep constants on the reference device/dtype and avoid Python
            # complex operands in Inductor's fallback kernels during backward.
            overrides[I] = f"_sx_tensor(1j, {reference})"
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
        return getattr(_core, self._native_method)(expr, expression_overrides=overrides or None, **{
            key: self._settings[key]
            for key in ("symbol_names", "user_functions", "reference", "mask")
        })

    # Common protocol methods permit subclass implementations to call super().
    _print_Symbol = _print_Integer = _print_Rational = _print_Float = emptyPrinter
    _print_Add = _print_Mul = _print_Pow = _print_Piecewise = emptyPrinter
    _print_UnevaluatedExpr = emptyPrinter

    def parenthesize(self, item, level, strict=False):
        from ...printing.precedence import precedence
        text = self._print(item)
        current = precedence(item)
        return f"({text})" if current < level or (not strict and current == level) else text


TorchPrinter = PyTorchPrinter


def torchcode(expr, assign_to=None, **settings):
    """Return a native PyTorch expression (see the module helper contract)."""
    return PyTorchPrinter(settings).doprint(expr, assign_to)


__all__ = ["PyTorchPrinter", "TorchPrinter", "torchcode"]
