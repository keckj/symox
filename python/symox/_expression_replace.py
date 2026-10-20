"""Bottom-up replacement with native matching and constructor validation."""


def replace(self, query, value, map=False, simultaneous=True, exact=None):
    from . import Expr, MatrixExpr, AlgebraicNumber, sympify, _is_function_head
    from ._core import Head
    from ._matrixexpr_methods import _TupleView, _ExprCondPairView
    from ._structural_queries import _graph, matcher, _SPECIAL_HEADS, _RELATIONS
    from . import __dict__ as namespace

    expressions=(Expr, MatrixExpr, AlgebraicNumber)
    if isinstance(query,str) or not callable(query) and not isinstance(query,expressions):
        try:query=sympify(query)
        except (TypeError,ValueError):
            raise TypeError('Replacement query must be a type, expression or callable') from None
    if not callable(value) and not isinstance(value,expressions):
        value=sympify(value)

    constructor=(isinstance(query,(type,Head)) or _is_function_head(query)
                 or bool(getattr(query,'_symox_head_name',None))
                 or (getattr(query,'__name__',None) in _SPECIAL_HEADS|_RELATIONS.keys()
                     and namespace.get(query.__name__) is query))
    pattern=isinstance(query,expressions)
    if constructor or not pattern:
        if not callable(value) or isinstance(value,expressions):
            raise TypeError('A type or predicate query requires a callable replacement')
        test=matcher(query,matching=False)
    if pattern and exact is None:
        _,records,_=_graph(query)
        wilds={node for node,_ in records.values() if getattr(node,'is_Symbol',False)
               and str(node)==node.name+'_'}
        exact=len(wilds)>1

    root,records,order=_graph(self)
    results={}
    changes={}
    for identifier in order:
        original,children=records[identifier]
        oldargs=tuple(records[child][0] for child in children)
        args=tuple(results[child] for child in children)
        node=original
        if any(new is not old for new,old in zip(args,oldargs)):
            if isinstance(original,_ExprCondPairView):node=type(original)(*args)
            elif isinstance(original,_TupleView):node=type(original)(args)
            elif getattr(original,'_symox_blockgrid',False):
                from . import _core
                node=_core._BlockGrid(original.rows,original.cols,args)
            elif getattr(original,'is_Piecewise',False):
                from . import Piecewise
                node=Piecewise(*args)
            elif getattr(original,'is_Order',False):
                node=original.func(args[0],*(item for pair in args[1:] for item in pair))
            elif getattr(original,'func',None) is not None:
                node=original.func(*(tuple(arg) if isinstance(arg,_TupleView) else arg for arg in args))
        # A reconstruction can collapse to one of the original children.
        # In simultaneous mode that child must not be transformed twice.
        skip=simultaneous and any(node==old and new!=old for old,new in zip(oldargs,args))
        replacement=None
        if not skip:
            if pattern:
                method=getattr(node,'match',None)
                if isinstance(query,AlgebraicNumber) and not isinstance(node,AlgebraicNumber):
                    bindings=None
                else:bindings=method(query) if method else ({} if node==query else None)
                if bindings is not None and (not exact or all(item!=0 for item in bindings.values())):
                    if callable(value) and not isinstance(value,expressions):
                        replacement=value(**{key.name:item for key,item in bindings.items()})
                    else:replacement=value.subs(bindings)
            elif test(node):
                replacement=value(*getattr(node,'args',())) if constructor else value(node)
        if replacement is not None:
            if not isinstance(replacement,expressions+(_TupleView,)):
                replacement=sympify(replacement)
            if replacement!=node:
                changes[node]=replacement
                node=replacement
        results[identifier]=node
    result=results[root]
    return (result,changes) if map else result
