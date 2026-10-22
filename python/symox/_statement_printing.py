"""Statement-level printing of codegen nodes for the native ``pycode``.

The native ``pycode`` prints expressions only. This wrapper prints the
``sympy.codegen.ast`` statement nodes with ports of the corresponding methods
of SymPy's ``PythonCodePrinter`` and delegates every expression they contain
to the native printer with the same settings.
"""

# Function names the reference printers know for the C99 functions of
# ``sympy.codegen.cfunctions`` (``_known_functions_math`` of pycode).
_PYTHON_CFUNCTIONS = {'Sqrt': 'sqrt', 'log2': 'log2', 'log10': 'log10', 'expm1': 'expm1',
                      'log1p': 'log1p', 'hypot': 'hypot', 'isnan': 'isnan', 'isinf': 'isinf'}


def _is_node(value):
    from .codegen.ast import CodegenAST
    from .codegen.fnodes import Extent
    if isinstance(value, (list, tuple)):
        return any(_is_node(v) for v in value)
    return isinstance(value, (CodegenAST, Extent))


def _none(value):
    return value is None or type(value).__name__ == 'NoneToken'


def _with_functions(settings, functions):
    user = dict(settings.get('user_functions') or {})
    for name, target in functions.items():
        user.setdefault(name, target)
    settings = dict(settings)
    settings['user_functions'] = user
    return settings


class _StatementPrinter:
    def __init__(self, native, settings):
        self._native = native
        self._settings = settings

    def _expr(self, expr):
        return self._native(expr, **self._settings)

    def _print(self, expr):
        from .codegen.ast import CodegenAST
        from .codegen.fnodes import Extent
        if isinstance(expr, str):
            return expr
        if isinstance(expr, Extent):
            return str(expr)
        if isinstance(expr, CodegenAST):
            for cls in type(expr).__mro__:
                method = getattr(self, '_print_' + cls.__name__, None)
                if method is not None:
                    return method(expr)
            raise NotImplementedError('Unsupported by %s: %s' % (type(self).__name__, type(expr).__name__))
        return self._expr(expr)

    def _print_String(self, expr):
        return str(expr.text)

    def _print_QuotedString(self, expr):
        return '"%s"' % expr.text

    def _print_Variable(self, expr):
        return self._print(expr.symbol)

    def _print_CodeBlock(self, expr):
        return '\n'.join(self._print(i) for i in expr.args)

    def _print_FunctionCall(self, expr):
        return '%s(%s)' % (expr.name, ', '.join(self._print(arg) for arg in expr.function_args))


class _PythonStatements(_StatementPrinter):
    """Port of the statement methods of SymPy's ``PythonCodePrinter``."""

    def _expr(self, expr):
        return self._native(expr, **_with_functions(self._settings, {k: 'math.' + v for k, v in _PYTHON_CFUNCTIONS.items()}))

    @staticmethod
    def _indent_codestring(codestring):
        return '\n'.join('    ' + line for line in codestring.split('\n'))

    def _print_Comment(self, expr):
        return '  # ' + str(expr)

    def _print_NoneToken(self, expr):
        return 'None'

    def _print_Assignment(self, expr):
        return '%s = %s' % (self._print(expr.lhs), self._print(expr.rhs))

    def _print_AugmentedAssignment(self, expr):
        return '{} {} {}'.format(self._print(expr.lhs), expr.op, self._print(expr.rhs))

    def _print_FunctionDefinition(self, fd):
        body = '\n'.join(self._print(arg) for arg in fd.body.args)
        return 'def {name}({parameters}):\n{body}'.format(
            name=fd.name, parameters=', '.join(self._print(var.symbol) for var in fd.parameters),
            body=self._indent_codestring(body))

    def _print_While(self, whl):
        body = '\n'.join(self._print(arg) for arg in whl.body.args)
        return 'while {cond}:\n{body}'.format(cond=self._print(whl.condition), body=self._indent_codestring(body))

    def _print_Declaration(self, decl):
        return '%s = %s' % (self._print(decl.variable.symbol), self._print(decl.variable.value))

    def _print_BreakToken(self, expr):
        return 'break'

    def _print_ContinueToken(self, expr):
        return 'continue'

    def _print_Return(self, ret):
        return 'return %s' % self._print(ret.return_)

    def _print_Raise(self, rs):
        return 'raise %s' % self._print(rs.exception)

    def _print_RuntimeError_(self, re):
        return 'RuntimeError(%s)' % self._print(re.message)

    def _print_Print(self, prnt):
        print_args = ', '.join(self._print(arg) for arg in prnt.print_args)
        if not _none(prnt.format_string):
            print_args = '{} % ({}), end=""'.format(self._print(prnt.format_string), print_args)
        if not _none(prnt.file):
            print_args += ', file=%s' % self._print(prnt.file)
        return 'print(%s)' % print_args

    def _print_Stream(self, strm):
        if str(strm.name) == 'stdout':
            return 'sys.stdout'
        elif str(strm.name) == 'stderr':
            return 'sys.stderr'
        return self._print(strm.name)


def _python_wrapper(native):
    def pycode(expr, assign_to=None, **settings):
        if assign_to is None and _is_node(expr):
            from .codegen.ast import CodeBlock
            if isinstance(expr, (list, tuple)):
                expr = CodeBlock(*expr)
            return _PythonStatements(native, settings)._print(expr)
        settings = _with_functions(settings, {k: 'math.' + v for k, v in _PYTHON_CFUNCTIONS.items()})
        try:
            return native(expr, assign_to, **settings)
        except NotImplementedError as exc:
            # SymPy's PythonCodePrinter raises PrintMethodNotImplementedError.
            prefix = 'Unsupported by the code printer: '
            message = str(exc)
            if type(exc) is not NotImplementedError or not message.startswith(prefix):
                raise
            from .printing.codeprinter import PrintMethodNotImplementedError
            from .printing.pycode import PythonCodePrinter
            raise PrintMethodNotImplementedError(
                f"Unsupported by {PythonCodePrinter}: {message[len(prefix):]}"
                "\nSet the printer option 'strict' to False in order to generate partially printed code."
            ) from None
    pycode.__doc__ = native.__doc__
    pycode.__wrapped__ = native
    return pycode
