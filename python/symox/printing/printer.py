"""Small extensible printer protocol, shared by the compatibility printers."""


class Printer:
    _default_settings = {"order": None}
    printmethod = None

    def __init__(self, settings=None, **kwargs):
        self._settings = dict(self._default_settings)
        if settings:
            self._settings.update(settings)
        self._settings.update(kwargs)

    @property
    def order(self):
        return self._settings.get("order")

    def _print(self, expr, **kwargs):
        if self.printmethod and hasattr(expr, self.printmethod):
            return getattr(expr, self.printmethod)(self, **kwargs)
        name = getattr(expr, "class_name", type(expr).__name__)
        method = getattr(self, "_print_" + name, None)
        if method is not None:
            return method(expr, **kwargs)
        # Rust nodes share the Expr Python type, while SymPy's singleton
        # numbers inherit the usual Integer/Rational printer methods.
        if name in ("Zero", "One", "NegativeOne"):
            method = getattr(self, "_print_Integer", None)
        elif name == "Half":
            method = getattr(self, "_print_Rational", None)
        if method is not None:
            return method(expr, **kwargs)
        for cls in type(expr).__mro__:
            method = getattr(self, "_print_" + cls.__name__, None)
            if method is not None:
                return method(expr, **kwargs)
        return self.emptyPrinter(expr)

    def doprint(self, expr):
        return self._print(expr)

    def emptyPrinter(self, expr):
        return str(expr)
