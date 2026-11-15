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
            return self._print_not_supported(expr)
        return self._expr(expr)

    # The reference printer class named by unsupported-node errors.
    _reference_printer = None

    def _print_not_supported(self, expr):
        # SymPy's CodePrinter._print_not_supported: strict (the default when
        # ``human`` is true) raises, otherwise the node is printed by the
        # empty printer and listed in front of the code.
        strict = self._settings.get('strict', None)
        if strict is None:
            strict = self._settings.get('human', True) == True
        if strict:
            from .printing.codeprinter import PrintMethodNotImplementedError
            reference = self._reference_printer() if self._reference_printer else type(self)
            raise PrintMethodNotImplementedError(
                f"Unsupported by {reference}: {type(expr)}"
                "\nSet the printer option 'strict' to False in order to generate partially printed code.")
        not_supported = self.__dict__.setdefault('_not_supported', set())
        try:
            not_supported.add(expr)
        except TypeError:
            pass
        return repr(expr)

    def doprint(self, expr):
        """SymPy's ``CodePrinter.doprint`` front lines for unsupported nodes."""
        self._not_supported = set()
        lines = self._print(expr).splitlines()
        if self._not_supported and self._settings.get('human', True):
            front = [self._get_comment('Not supported in Python:')]
            front += [self._get_comment(type(e).__name__) for e in sorted(self._not_supported, key=str)]
            lines = front + lines
        return '\n'.join(lines)

    @staticmethod
    def _get_comment(text):
        return '  # {}'.format(text)

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


def _python_code_printer():
    from .printing.pycode import PythonCodePrinter
    return PythonCodePrinter


class _PythonStatements(_StatementPrinter):
    """Port of the statement methods of SymPy's ``PythonCodePrinter``."""
    _reference_printer = staticmethod(_python_code_printer)

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


_UNSUPPORTED_PREFIX = 'Unsupported by the code printer: '


def _unsupported_names(native, expr, assign_to, settings):
    """The function names the native printer refuses in strict mode."""
    names = []
    probe = dict(settings)
    probe['strict'] = True
    for _ in range(64):
        try:
            native(expr, assign_to, **probe)
            break
        except NotImplementedError as exc:
            message = str(exc)
            if type(exc) is not NotImplementedError or not message.startswith(_UNSUPPORTED_PREFIX):
                break
            name = message[len(_UNSUPPORTED_PREFIX):]
            if name in names:
                break
            names.append(name)
            probe['user_functions'] = dict(probe.get('user_functions') or {}, **{name: name})
        except Exception:
            break
    return names


def _function_name(node):
    func = getattr(node, 'func', None)
    return getattr(func, '__name__', None) or getattr(node, 'class_name', None)


def _not_supported_front(expr, names):
    """SymPy's ``CodePrinter.doprint`` front lines: one per distinct
    unsupported node (the printer does not descend into it), sorted by str."""
    found = {}
    stack = [expr]
    while stack:
        node = stack.pop()
        if not hasattr(node, 'args') or isinstance(node, type):
            continue
        name = _function_name(node)
        if name in names:
            found.setdefault(node, name)
            continue
        stack.extend(node.args)
    return ['  # Not supported in Python:'] + ['  # ' + found[node] for node in sorted(found, key=str)]


def _python_wrapper(native):
    def pycode(expr, assign_to=None, **settings):
        # SymPy's AbstractPythonCodePrinter.__init__ check.
        if settings.get('standard') not in (None, 'python3'):
            raise ValueError('Only Python 3 is supported.')
        if assign_to is None and _is_node(expr):
            from .codegen.ast import CodeBlock
            if isinstance(expr, (list, tuple)):
                expr = CodeBlock(*expr)
            return _PythonStatements(native, settings).doprint(expr)
        from .codegen.abstract_nodes import List
        if assign_to is None and isinstance(expr, List):
            # SymPy's AbstractPythonCodePrinter._print_List.
            return '[%s]' % ', '.join(pycode(arg, **settings) for arg in expr.args)
        settings = _with_functions(settings, {k: 'math.' + v for k, v in _PYTHON_CFUNCTIONS.items()})
        strict = settings.get('strict', True)
        if strict is None:
            strict = settings.get('human', True) == True
        if not strict and settings.get('human', True):
            # Non-strict output lists the unsupported functions in front.
            code = native(expr, assign_to, **settings)
            names = _unsupported_names(native, expr, assign_to, settings)
            if names and isinstance(code, str):
                code = '\n'.join(_not_supported_front(expr, names) + code.splitlines())
            return code
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
