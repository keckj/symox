"""Extensible C++ printers using the SymPy ``_print_Class`` protocol.

The native printers remain available as ``symox.cxxcode``.  This layer
walks expressions in Python so subclasses can customize every child node.
"""

from .precedence import PRECEDENCE, precedence
from .printer import Printer


class CXX17CodePrinter(Printer):
    standard = "C++17"
    language = "C++"
    printmethod = "_cxxcode"
    _default_settings = {"order": None, "precision": 17, "user_functions": {}, "strict": True}

    def __init__(self, settings=None, **kwargs):
        super().__init__(settings, **kwargs)
        self.known_functions = dict(self._settings.get("user_functions", {}))

    def doprint(self, expr, assign_to=None):
        from symox import sympify
        from symox.codegen.ast import Assignment
        if assign_to is not None:
            expr = Assignment(assign_to, expr)
        if not hasattr(expr, "args"):
            expr = sympify(expr)
        return self._print(expr)

    def parenthesize(self, item, level, strict=False):
        value = self._print(item)
        p = precedence(item)
        return f"({value})" if p < level or (not strict and p == level) else value

    def _print_Symbol(self, expr):
        return expr.name

    def _print_Dummy(self, expr):
        from symox import _core
        return _core.cxxcode(expr)

    def _print_MatrixElement(self, expr):
        parent, row, column = expr.args
        if getattr(parent, "kind", None) != "MatrixSymbol":
            raise NotImplementedError("Materialize the matrix operation before printing a scalar entry")
        return f"{parent.name}[{self._print(row * parent.cols + column)}]"

    def _print_Integer(self, expr):
        return str(int(expr))

    def _print_Rational(self, expr):
        return f"{expr.p}.0/{expr.q}.0"

    def _print_Float(self, expr):
        from symox import _core
        return _core.cxxcode(expr)

    def _print_Add(self, expr):
        args = expr.args
        if self.order not in ("none", "old"):
            args = expr.as_ordered_terms()
        parts = []
        for term in args:
            printed = self._print(term)
            if printed.startswith("-"):
                parts.append((" - " if parts else "-") + printed[1:])
            else:
                parts.append((" + " if parts else "") + printed)
        return "".join(parts)

    def _print_Mul(self, expr):
        from symox import Integer
        args = expr.args
        if self.order not in ("none", "old"):
            args = expr.as_ordered_factors()
        numerator, denominator = [], []
        sign = ""
        for item in args:
            if getattr(item, "is_Number", False) and item.is_negative:
                sign = "" if sign else "-"
                item = -item
                if item == 1:
                    continue
            if (item.is_commutative and item.is_Pow and item.exp.is_Rational
                    and item.exp.is_negative):
                denominator.append(item.args[0] ** (-item.exp))
            else:
                numerator.append(item)
        if not numerator:
            numerator = [Integer(1)]
        level = precedence(expr)
        if len(numerator) == 1 and sign:
            top = self.parenthesize(numerator[0], (PRECEDENCE["Mul"] + PRECEDENCE["Pow"]) / 2)
        else:
            top = "*".join(self.parenthesize(item, level) for item in numerator)
        if not denominator:
            return sign + top
        bottom = "*".join(self.parenthesize(item, level) for item in denominator)
        if len(denominator) > 1:
            bottom = f"({bottom})"
        return sign + top + "/" + bottom

    def _print_Pow(self, expr):
        base, exponent = expr.args
        if exponent == -1:
            return "1.0/" + self.parenthesize(base, PRECEDENCE["Pow"])
        from symox import Rational
        if exponent == Rational(1, 2):
            return f"std::sqrt({self._print(base)})"
        if exponent == Rational(1, 3):
            return f"std::cbrt({self._print(base)})"
        return f"std::pow({self._print(base)}, {self._print(exponent)})"

    def _print_UnevaluatedExpr(self, expr):
        return self._print(expr.args[0])

    def _print_Assignment(self, expr):
        return f"{self._print(expr.lhs)} = {self._print(expr.rhs)};"

    def _print_Comment(self, expr):
        return "/* " + expr.text + " */"

    def _print_CodeBlock(self, expr):
        # SymPy's block printer terminates comments as statements, too.
        return "\n".join(self._print(arg).rstrip(";") + ";" for arg in expr.args)

    def _print_IndexedBase(self, expr):
        return self._print(expr.label)

    def _print_Indexed(self, expr):
        from symox import Integer
        if len(expr.indices) == 1:
            index = expr.indices[0]
        elif expr.base.shape is not None:
            index = Integer(0)
            stride = Integer(1)
            for idx, size in zip(reversed(expr.indices), reversed(expr.base.shape)):
                index += idx * stride
                stride *= size
        else:
            raise ValueError("Multidimensional indexed references need a shape")
        return f"{self._print(expr.base)}[{self._print(index)}]"

    def emptyPrinter(self, expr):
        name = getattr(expr, "class_name", type(expr).__name__)
        if getattr(expr, "is_Function", False):
            function = self.known_functions.get(name)
            if function is None:
                if (getattr(expr, "is_UndefinedFunction", False)
                        or any(cls.__name__ == "AppliedUndef" for cls in type(expr).__mro__)):
                    raise TypeError(f"Undefined function {name} needs an explicit user_functions mapping")
                function = {"Abs": "std::fabs", "ceiling": "std::ceil"}.get(name, "std::" + name)
            if not isinstance(function, str):
                function = next((f for test, f in function if test(*expr.args)), name)
            return f"{function}({', '.join(self._print(a) for a in expr.args)})"
        # Constants and less common nodes use the existing Rust code printer.
        from symox import _core
        try:
            return _core.cxxcode(expr)
        except (TypeError, ValueError):
            if self._settings.get("strict", True):
                raise TypeError(f"Unsupported expression in C++ printer: {name}") from None
            return str(expr)


# The supported expression subset shares its formatting across standards.
CXX11CodePrinter = CXX17CodePrinter
CXX98CodePrinter = CXX17CodePrinter


def cxxcode(expr, assign_to=None, standard="c++11", **settings):
    return CXX17CodePrinter(settings).doprint(expr, assign_to)
