"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox.interactive', 'init_session')


def int_to_Integer(s):
    import io,tokenize
    tokens=[]
    for token in tokenize.generate_tokens(io.StringIO(s).readline):
        kind,text=token.type,token.string
        if kind==tokenize.NUMBER:
            try:int(text.replace('_',''),0)
            except ValueError:
                try:is_integer=text.replace('_','').isdigit()
                except AttributeError:is_integer=False
            else:is_integer=True
            if is_integer:
                tokens.extend([(tokenize.NAME,'Integer'),(tokenize.OP,'('),(kind,text),(tokenize.OP,')')]);continue
        tokens.append((kind,text))
    return tokenize.untokenize(tokens)


def enable_automatic_int_sympification(shell):
    from symox import Integer
    shell.user_ns.setdefault('Integer',Integer)
    def transform(lines):return int_to_Integer(''.join(lines)).splitlines(keepends=True)
    shell.input_transformers_cleanup.append(transform)
    return transform


def enable_automatic_symbols(shell):
    """Create missing loaded identifiers through an IPython AST transformer."""
    import ast,builtins
    from symox import Symbol
    class AutomaticSymbols(ast.NodeTransformer):
        def visit_Module(self,node):
            bound={n.id for n in ast.walk(node) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Store)}
            known=set(shell.user_ns)|set(vars(builtins))|bound
            for part in ast.walk(node):
                if isinstance(part,ast.Name) and isinstance(part.ctx,ast.Load) and part.id not in known:
                    shell.user_ns[part.id]=Symbol(part.id);known.add(part.id)
            return node
    transform=AutomaticSymbols();shell.ast_transformers.append(transform);return transform


def init_ipython_session(shell=None,argv=None,auto_symbols=False,auto_int_to_Integer=False):
    if shell is None:
        from IPython.terminal.interactiveshell import TerminalInteractiveShell
        shell=TerminalInteractiveShell.instance()
    if auto_symbols:enable_automatic_symbols(shell)
    if auto_int_to_Integer:enable_automatic_int_sympification(shell)
    return shell


def init_python_session():
    import code,symox
    return code.InteractiveConsole({name:getattr(symox,name) for name in symox.__all__})

__all__ += ['int_to_Integer','enable_automatic_int_sympification','enable_automatic_symbols','init_ipython_session','init_python_session']
