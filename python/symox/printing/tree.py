"""Native printing helpers adapted from SymPy 1.14; see LICENSE-SYMPY."""
def pprint_nodes(subtrees):
    """Prettyprints systems of nodes."""

    def indent(s, type=1):
        x = s.split('\n')
        r = '+-%s\n' % x[0]
        for a in x[1:]:
            if a == '':
                continue
            if type == 1:
                r += '| %s\n' % a
            else:
                r += '  %s\n' % a
        return r
    if not subtrees:
        return ''
    f = ''
    for a in subtrees[:-1]:
        f += indent(a)
    f += indent(subtrees[-1], 2)
    return f

def print_node(node, assumptions=True):
    """Returns information about the "node"."""
    s = '%s: %s\n' % (getattr(node, 'class_name', type(node).__name__), str(node))
    if assumptions:
        d = getattr(node, 'assumptions0', {})
    else:
        d = None
    if d:
        for a in sorted(d):
            v = d[a]
            if v is None:
                continue
            s += '%s: %s\n' % (a, v)
    return s

def tree(node, assumptions=True):
    """Returns a tree representation of "node" as a string."""
    subtrees = []
    for arg in node.args:
        subtrees.append(tree(arg, assumptions=assumptions))
    s = print_node(node, assumptions=assumptions) + pprint_nodes(subtrees)
    return s

def print_tree(node, assumptions=True):
    """Prints a tree representation of "node"."""
    print(tree(node, assumptions=assumptions))
