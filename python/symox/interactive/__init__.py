"""Explicit interactive display and session initialization."""
import sys


def init_printing(pretty_print=True,order=None,use_unicode=None,use_latex=None,wrap_line=None,
                  num_columns=None,no_global=False,ip=None,euler=False,forecolor=None,backcolor=None,
                  fontsize=None,latex_mode='plain',print_builtin=True,str_printer=None,
                  pretty_printer=None,latex_printer=None,scale=1.0,**settings):
    from .. import Expr,pretty,sstr,latex
    plain = str_printer or sstr
    prettify = pretty_printer or pretty
    latexify = latex_printer or latex
    def render(expr):
        options = dict(settings)
        if use_unicode is not None:options['use_unicode'] = use_unicode
        if num_columns is not None:options['num_columns'] = num_columns
        if wrap_line is not None:options['wrap_line'] = wrap_line
        if order is not None:options['order'] = order
        return prettify(expr,**options) if pretty_print else plain(expr)
    if ip is None:
        try:
            from IPython import get_ipython
            ip = get_ipython()
        except ImportError:pass
    if ip is not None:
        ip.display_formatter.formatters['text/plain'].for_type(Expr,lambda expr,p,cycle:p.text(render(expr)))
        if use_latex:
            ip.display_formatter.formatters['text/latex'].for_type(Expr,lambda expr:'$'+latexify(expr)+'$')
    if not no_global:
        def displayhook(expr):
            if expr is None:return
            import builtins
            builtins._ = None
            print(render(expr))
            builtins._ = expr
        sys.displayhook = displayhook


def init_session(ipython=None,pretty_print=True,order=None,use_unicode=None,use_latex=None,
                 quiet=False,auto_symbols=False,auto_int_to_Integer=False,str_printer=None,
                 pretty_printer=None,latex_printer=None,argv=None):
    if auto_symbols or auto_int_to_Integer:
        raise NotImplementedError('automatic source transformations are not supported by the session console')
    import symox
    namespace = {name:getattr(symox,name) for name in symox.__all__}
    exec('x, y, z, t = symbols("x y z t")\nk, m, n = symbols("k m n", integer=True)\nf, g, h = symbols("f g h", cls=Function)',namespace)
    init_printing(pretty_print,order,use_unicode,use_latex,str_printer=str_printer,
                  pretty_printer=pretty_printer,latex_printer=latex_printer)
    banner = '' if quiet else 'Symox %s interactive session' % symox.__version__
    if ipython:
        from IPython import start_ipython
        start_ipython(argv=argv or [],user_ns=namespace)
    else:
        import code
        code.InteractiveConsole(namespace).interact(banner=banner)


__all__ = ['init_printing','init_session']

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'interactive_traversal')
