"""The lenient printer protocol of Symox's own language printers.

SymPy's Printer, StrPrinter and CodePrinter are ported in printer.py,
str.py, codeprinter.py, c.py and cxx.py. The LaTeX, MathML, pretty, repr,
Fortran, Rust, Python and scalar language printers of Symox were written
against this smaller protocol: it accepts settings as keyword arguments,
dispatches on the printed class name and falls back to the native printers.
"""
from .precedence import PRECEDENCE, precedence
from .printer import Printer
from ._nodes import NodePrinter


def _code_error():
    from .codeprinter import PrintMethodNotImplementedError
    return PrintMethodNotImplementedError


class LegacyPrinter(Printer):
    _default_settings = {"order": None}
    printmethod = None

    def __init__(self, settings=None, **kwargs):
        self._str = str
        self._context = {}
        self._print_level = 0
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
        if name=='_MatrixExpr':
            name=expr._kind_tag
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
        from .. import _core
        for cls in type(expr).__mro__:
            method = getattr(self, "_print_" + cls.__name__, None)
            if method is None and cls is _core.Matrix:
                # The native dense base (DenseMatrix) keeps the printers'
                # historical ``_print_Matrix`` handlers.
                method = getattr(self, "_print_Matrix", None)
            if method is not None:
                return method(expr, **kwargs)
        return self.emptyPrinter(expr)

    def doprint(self, expr):
        return self._print(expr)

    def _as_ordered_terms(self, expr, order=None):
        """A compatibility function for ordering terms in Add. """
        order = order or self.order
        if order == 'old':
            from functools import cmp_to_key
            from symox import Add
            return sorted(Add.make_args(expr), key=cmp_to_key(_compare_pretty))
        elif order == 'none':
            return list(expr.args)
        else:
            return expr.as_ordered_terms(order=order)

    def emptyPrinter(self, expr):
        return str(expr)


def _compare_pretty(a, b):
    """SymPy's ``Basic._compare_pretty``: the ``order='old'`` comparison."""
    from symox import Order, Wild
    if isinstance(a, Order) and not isinstance(b, Order):
        return 1
    if not isinstance(a, Order) and isinstance(b, Order):
        return -1
    if a.is_Rational and b.is_Rational:
        l = a.p * b.q
        r = b.p * a.q
        return (l > r) - (l < r)
    else:
        p1, p2, p3 = Wild("p1"), Wild("p2"), Wild("p3")
        r_a = a.match(p1 * p2**p3)
        if r_a and p3 in r_a:
            a3 = r_a[p3]
            r_b = b.match(p1 * p2**p3)
            if r_b and p3 in r_b:
                b3 = r_b[p3]
                c = a3.compare(b3)
                if c != 0:
                    return c
    return a.compare(b)


class LegacyCXX17CodePrinter(LegacyPrinter):
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

    def _print_MatrixSymbol(self, expr):
        return self._print_Symbol(expr)

    def _print_Dummy(self, expr):
        from symox import _core
        return _core.cxxcode(expr)

    def _print_MatrixElement(self, expr):
        parent, row, column = expr.args
        if getattr(parent, "kind", None) in ("FunctionMatrix", "ElementwiseApplyFunction"):
            return self._print(parent[row, column])
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
        from symox import MatrixSymbol
        from symox.codegen.ast import Assignment
        if isinstance(expr.lhs,MatrixSymbol):
            if expr.rhs.shape!=expr.lhs.shape:
                raise ValueError('Dimensions of lhs and rhs do not align.')
            return '\n'.join(self._print(Assignment(expr.lhs[i,j],expr.rhs[i,j]))
                             for i,j in self._traverse_matrix_indices(expr.lhs))
        return f"{self._print(expr.lhs)} = {self._print(expr.rhs)};"

    def _traverse_matrix_indices(self, matrix):
        rows,columns=matrix.shape
        return ((i,j) for i in range(rows) for j in range(columns))

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


class LegacyCodePrinter(NodePrinter,LegacyCXX17CodePrinter):
    language='C';printmethod='_ccode'
    def __init__(self,settings=None,**kwargs):
        super().__init__(settings,**kwargs);self.headers=set();self.libraries=set();self.macros=set()
    def _print_Pow(self,expr):
        from symox import Rational
        base,exponent=expr.args
        if exponent==Rational(1,2):return 'sqrt('+self._print(base)+')'
        if exponent==Rational(1,3):return 'cbrt('+self._print(base)+')'
        return 'pow('+self._print(base)+', '+self._print(exponent)+')'
    def emptyPrinter(self,expr):
        if getattr(expr,'is_Function',False):
            name=expr.class_name;function=self.known_functions.get(name)
            if function is None:
                if getattr(expr,'is_UndefinedFunction',False):raise _code_error()('unknown function '+name+' requires user_functions')
                function={'Abs':'fabs','ceiling':'ceil','Max':'fmax','Min':'fmin'}.get(name,name)
            if not isinstance(function,str):function=next((value for condition,value in function if condition(*expr.args)),name)
            return function+'('+', '.join(self._print(a) for a in expr.args)+')'
        from symox import _core
        return _core.ccode(expr)



def _c_relation(self,expr,op):return '('+self._print(expr.lhs)+' '+op+' '+self._print(expr.rhs)+')'
LegacyCodePrinter._print_And=lambda self,e:'('+' && '.join(self._print(a) for a in e.args)+')'
LegacyCodePrinter._print_Or=lambda self,e:'('+' || '.join(self._print(a) for a in e.args)+')'
LegacyCodePrinter._print_Not=lambda self,e:'!('+self._print(e.args[0])+')'
LegacyCodePrinter._print_Xor=lambda self,e:'('+' ^ '.join(self._print(a) for a in e.args)+')'
LegacyCodePrinter._print_ITE=lambda self,e:'('+self._print(e.args[0])+' ? '+self._print(e.args[1])+' : '+self._print(e.args[2])+')'
for _name,_op in [('Equality','=='),('Unequality','!='),('LessThan','<='),('StrictLessThan','<'),('GreaterThan','>='),('StrictGreaterThan','>')]:
    setattr(LegacyCodePrinter,'_print_'+_name,lambda self,e,op=_op:_c_relation(self,e,op))
def _c_minmax(self,expr):
    function='fmin' if expr.class_name=='Min' else 'fmax';value=self._print(expr.args[0])
    for arg in expr.args[1:]:value=function+'('+value+', '+self._print(arg)+')'
    return value
LegacyCodePrinter._print_Min=LegacyCodePrinter._print_Max=_c_minmax
