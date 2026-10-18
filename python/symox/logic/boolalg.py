"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'And Or Not Xor Nand Nor Implies Equivalent ITE true false to_cnf to_dnf to_nnf POSform SOPform simplify_logic bool_map')


def gateinputcount(expr):
    """Count wires entering Boolean gates, including nested gate occurrences."""
    from symox import sympify
    expr=sympify(expr)
    gates={'And','Or','Not','Xor','Nand','Nor','ITE','Implies','Equivalent'}
    if expr.class_name not in gates:return 0
    return len(expr.args)+sum(gateinputcount(a) for a in expr.args)
__all__ += ['gateinputcount']
