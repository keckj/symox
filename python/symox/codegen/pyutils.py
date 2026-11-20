""" This module collects utilities for rendering Python code. """
import io
import tokenize


def _imported_modules(source):
    # SymPy's PythonCodePrinter records ``module_imports`` while printing;
    # symox's native printer emits the same fully qualified names, so the
    # modules are read back from the ``module.name`` references in the code.
    modules = []
    try:
        tokens = list(tokenize.generate_tokens(io.StringIO(source).readline))
    except (tokenize.TokenError, IndentationError, SyntaxError):
        tokens = []
    for prev, tok, nxt in zip(tokens, tokens[1:], tokens[2:]):
        if (tok.type == tokenize.NAME and tok.string in ('math', 'cmath')
                and nxt.type == tokenize.OP and nxt.string == '.'
                and not (prev.type == tokenize.OP and prev.string == '.')
                and tok.string not in modules):
            modules.append(tok.string)
    return modules


def render_as_module(content, standard='python3'):
    """Renders Python code as a module (with the required imports).

    Parameters
    ==========

    standard :
        See the parameter ``standard`` in
        :meth:`symox.printing.pycode.pycode`
    """
    from symox import pycode
    pystr = pycode(content, standard=standard)
    module_imports_str = '\n'.join('import %s' % k for k in _imported_modules(pystr))
    return module_imports_str + '\n\n' + pystr
