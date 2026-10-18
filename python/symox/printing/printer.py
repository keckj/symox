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
        undefined = (getattr(expr, "is_UndefinedFunction", False)
                     or any(cls.__name__ == "AppliedUndef" for cls in type(expr).__mro__))
        if undefined:
            # Native nodes share Expr's Python class.  A coincidentally named
            # undefined function must not inherit a builtin printer handler.
            # Explicit handlers supplied by a user subclass remain available.
            if method is not None and not getattr(method, "__module__", "").startswith("symox.printing."):
                return method(expr, **kwargs)
            for head in ("UndefinedFunction", "AppliedUndef", "Function"):
                function_method = getattr(self, "_print_" + head, None)
                if function_method is not None:
                    return function_method(expr, **kwargs)
            return self.emptyPrinter(expr)
        if method is not None:
            return method(expr, **kwargs)
        # Rust nodes share the Expr Python type, while SymPy's singleton
        # numbers inherit the usual Integer/Rational printer methods.
        if name in ("Zero", "One", "NegativeOne"):
            method = getattr(self, "_print_Integer", None)
        elif name == "Half":
            method = getattr(self, "_print_Rational", None)
        elif name in ("Dummy", "Wild"):
            method = getattr(self, "_print_Symbol", None)
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

from contextlib import contextmanager

@contextmanager
def printer_context(printer,**kwargs):
    original=dict(getattr(printer,'_context',{}))
    printer._context={**original,**kwargs}
    try:yield
    finally:printer._context=original


def print_function(print_cls):
    def decorate(function):
        import inspect
        signature=inspect.signature(function)
        parameters=[p for p in signature.parameters.values() if p.kind!=p.VAR_KEYWORD]
        existing={p.name for p in parameters}
        for name,default in print_cls._default_settings.items():
            if name not in existing and name.isidentifier():parameters.append(inspect.Parameter(name,inspect.Parameter.KEYWORD_ONLY,default=default))
        function.__signature__=signature.replace(parameters=parameters)
        return function
    return decorate
