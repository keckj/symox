"""Native printing helpers adapted from SymPy 1.14; see LICENSE-SYMPY."""
from symox import Basic
from symox import Expr
from symox import Symbol
from symox import Integer, Rational, Float
from symox.printing.repr import srepr
__all__ = ['dotprint']
default_styles = ((Basic, {'color': 'blue', 'shape': 'ellipse'}), (Expr, {'color': 'black'}))
slotClasses = (Symbol, Integer, Rational, Float)

def purestr(x, with_args=False):
    """A string that follows ```obj = type(obj)(*obj.args)``` exactly."""
    sargs = ()
    if not isinstance(x, Basic):
        rv = str(x)
    elif not x.args:
        rv = srepr(x)
    else:
        args = x.args
        sargs = tuple(map(purestr, args))
        rv = '%s(%s)' % (getattr(x, 'class_name', type(x).__name__), ', '.join(sargs))
    if with_args:
        rv = (rv, sargs)
    return rv

def styleof(expr, styles=default_styles):
    """Merge style dictionaries in order"""
    style = {}
    for (typ, sty) in styles:
        if isinstance(expr, typ):
            style.update(sty)
    return style

def attrprint(d, delimiter=', '):
    """Print a dictionary of attributes"""
    import json
    return delimiter.join(json.dumps(str(k),ensure_ascii=False)+"="+json.dumps(str(v),ensure_ascii=False) for k,v in sorted(d.items()))

def dotnode(expr, styles=default_styles, labelfunc=str, pos=(), repeat=True):
    """String defining a node"""
    style = styleof(expr, styles)
    if isinstance(expr, Basic) and (not expr.is_Atom):
        label = getattr(expr, 'class_name', type(expr).__name__)
    else:
        label = labelfunc(expr)
    style['label'] = label
    expr_str = purestr(expr)
    if repeat:
        expr_str += '_%s' % str(pos)
    import json
    return '%s [%s];' % (json.dumps(expr_str,ensure_ascii=False), attrprint(style))

def dotedges(expr, atom=lambda x: not isinstance(x, Basic), pos=(), repeat=True):
    """List of strings for all expr->expr.arg pairs"""
    if atom(expr):
        return []
    else:
        (expr_str, arg_strs) = purestr(expr, with_args=True)
        if repeat:
            expr_str += '_%s' % str(pos)
            arg_strs = ['%s_%s' % (a, str(pos + (i,))) for (i, a) in enumerate(arg_strs)]
        import json
        return [json.dumps(expr_str,ensure_ascii=False)+" -> "+json.dumps(a,ensure_ascii=False)+";" for a in arg_strs]
template = 'digraph{\n\n# Graph style\n%(graphstyle)s\n\n#########\n# Nodes #\n#########\n\n%(nodes)s\n\n#########\n# Edges #\n#########\n\n%(edges)s\n}'
_graphstyle = {'rankdir': 'TD', 'ordering': 'out'}

def dotprint(expr, styles=default_styles, atom=lambda x: not isinstance(x, Basic), maxdepth=None, repeat=True, labelfunc=str, **kwargs):
    """DOT description of a Symox expression tree"""
    graphstyle = _graphstyle.copy()
    graphstyle.update(kwargs)
    nodes = []
    edges = []

    def traverse(e, depth, pos=()):
        nodes.append(dotnode(e, styles, labelfunc=labelfunc, pos=pos, repeat=repeat))
        if maxdepth is not None and depth >= maxdepth:
            return
        edges.extend(dotedges(e, atom=atom, pos=pos, repeat=repeat))
        [traverse(arg, depth + 1, pos + (i,)) for (i, arg) in enumerate(e.args) if not atom(arg)]
    traverse(expr, 0)
    return template % {'graphstyle': attrprint(graphstyle, delimiter='\n'), 'nodes': '\n'.join(nodes), 'edges': '\n'.join(edges)}
