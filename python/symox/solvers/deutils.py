"""Structural differential-equation analysis on native expression graphs."""


def ode_order(expr, func):
    """Return the highest total held derivative order involving ``func``.

    Derivative powers do not multiply the order. Derivatives of composite
    expressions add their own order to the highest order in their body.
    Native derivative counts are read without enumerating differentiation
    variables, and repeated expression subgraphs are visited once.
    """
    from .. import _core, sympify, Set, _wrap_set
    from .._structural_queries import _graph
    expr=sympify(expr)
    func=sympify(func)
    if not isinstance(expr,_core.Expr) or isinstance(_wrap_set(expr),Set):
        raise TypeError('Differential order requires a scalar expression or equation')
    if not isinstance(func,_core.Expr) or isinstance(_wrap_set(func),Set):
        raise TypeError('The dependent function must be a native scalar expression')
    root,records,postorder=_graph(expr)
    contains={}
    orders={}
    for identifier in postorder:
        value,children=records[identifier]
        present=value==func or any(contains[child] for child in children)
        contains[identifier]=present
        if not present:orders[identifier]=0
        elif getattr(value,'is_Derivative',False):
            count=sum(int(count) for _,count in value.variable_count)
            if value.expr==func:orders[identifier]=count
            else:
                _,body_children=records[children[0]]
                orders[identifier]=count+max((orders[child] for child in body_children),default=0)
        else:orders[identifier]=max((orders[child] for child in children),default=0)
    return orders[root]


def homogeneous_order(eq, *symbols):
    """Prove a common degree under positive simultaneous scaling.

    Arithmetic combinations, constant powers and degree-zero function
    arguments are handled structurally. A missing proof returns ``None``.
    Applied functions can be independent coordinates when all their arguments
    are also supplied. Principal-power branches are preserved.
    """
    from .. import S, Dummy, re, conjugate
    from ..calculus._analysis import scalar, walk, name
    from .._structural_queries import _graph
    if not symbols:raise ValueError('Supply at least one homogeneous coordinate')
    expression=scalar(eq)
    selected={scalar(symbol) for symbol in symbols}
    nodes=list(walk(expression))
    if any(node.is_Order or node.is_Derivative for node in nodes):return None
    if expression.is_number is True:return S.Zero
    replacements={}
    coordinates=set()
    for symbol in selected:
        if symbol.is_Symbol:coordinates.add(symbol)
        elif symbol.is_Function or symbol.is_UndefinedFunction:
            if not set(symbol.args).issubset(selected):return None
            replacements[symbol]=Dummy('homogeneous_coordinate')
            coordinates.add(replacements[symbol])
        else:return None
    expression=expression.xreplace(replacements)
    if not expression.free_symbols.intersection(coordinates):return None
    root,records,postorder=_graph(expression)
    degrees={}
    dependent={}
    for identifier in postorder:
        value,children=records[identifier]
        varying=value in coordinates or any(dependent[child] for child in children)
        dependent[identifier]=varying
        degree=None
        if not varying:degree=S.Zero
        elif value in coordinates:degree=S.One
        elif getattr(value,'is_Add',False):
            values=[degrees[child] for child in children]
            if values and all(item is not None and item==values[0] for item in values):degree=values[0]
        elif getattr(value,'is_Mul',False):
            values=[degrees[child] for child in children]
            if all(item is not None for item in values):degree=sum(values,S.Zero)
        elif getattr(value,'is_Pow',False):
            base,exponent=children
            base_degree=degrees[base]
            if not dependent[exponent] and base_degree is not None:
                # A real degree scales the base by a positive real number;
                # an integer outer power also permits complex scaling.
                if base_degree.is_real is True or value.exp.is_integer is True:
                    degree=base_degree*value.exp
        elif getattr(value,'is_Function',False) or getattr(value,'is_UndefinedFunction',False):
            values=[degrees[child] for child in children]
            head=name(value)
            if values and all(item==0 for item in values):degree=S.Zero
            elif len(values)==1 and values[0] is not None:
                inner=values[0]
                if head=='Abs':degree=re(inner)
                elif head=='conjugate':degree=conjugate(inner)
                elif head in ('re','im') and inner.is_real is True:degree=inner
                elif head in ('sign','arg') and inner.is_real is True:degree=S.Zero
        degrees[identifier]=degree
    return degrees[root]


__all__=['ode_order','homogeneous_order']
