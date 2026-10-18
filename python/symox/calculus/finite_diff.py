"""Finite differences on arbitrary grids using native scalar arithmetic."""
import itertools

from .. import _core

MAX_POINTS = 256
MAX_TABLE = 100_000


def _scalar(value):
    from .. import sympify
    value=sympify(value)
    if not isinstance(value,_core.Expr):
        raise TypeError('Finite differences require scalar expressions')
    return value


def _order(value):
    value=_scalar(value)
    if not value.is_number or value.is_finite is False:
        raise ValueError('Derivative order must be a nonnegative integer')
    if value.is_real is False or value<0 or value>MAX_POINTS:
        raise ValueError('Finite difference order must lie between 0 and 256')
    try:
        integer=int(value)
    except (ValueError,TypeError,OverflowError):
        raise ValueError('Derivative order must be a nonnegative integer') from None
    if (value-integer).is_zero is not True:
        raise ValueError('Derivative order must be a nonnegative integer')
    if integer>MAX_POINTS:
        raise ValueError('Finite difference order exceeds 256')
    return integer


def _points(values):
    try:
        values=list(itertools.islice(iter(values),MAX_POINTS+1))
    except TypeError:
        raise TypeError('Grid points must be an iterable') from None
    if not values:
        raise ValueError('At least one grid point is required')
    if len(values)>MAX_POINTS:
        raise ValueError('Finite difference grids are limited to 256 points')
    values=[_scalar(value) for value in values]
    for i,value in enumerate(values):
        if value.is_finite is False:
            raise ValueError('Grid points must be finite values or symbolic scalars')
        if any((value-previous).is_zero is True for previous in values[:i]):
            raise ValueError('Grid points must be distinct')
    return values


def finite_diff_weights(order, x_list, x0=1):
    """Weights for every derivative through order and every grid prefix.

    The incremental Fornberg recurrence retains exact symbolic coefficients.
    Unknown symbolic point differences are assumed nonzero.
    """
    from .. import S
    order=_order(order)
    points=_points(x_list)
    x0=_scalar(x0)
    size=len(points)
    if (order+1)*size*size>MAX_TABLE:
        raise ValueError('Finite difference weight table exceeds 100000 coefficients')
    table=[[[S.Zero]*size for _ in points] for _ in range(order+1)]
    table[0][0][0]=S.One
    previous_product=S.One
    for n in range(1,size):
        product=S.One
        for j in range(n):
            difference=points[n]-points[j]
            product *= difference
            for derivative in range(min(n,order)+1):
                lower=table[derivative-1][n-1][j] if derivative else S.Zero
                table[derivative][n][j]=((points[n]-x0)*table[derivative][n-1][j]
                                        -derivative*lower)/difference
        for derivative in range(min(n,order)+1):
            lower=table[derivative-1][n-1][n-1] if derivative else S.Zero
            table[derivative][n][n]=previous_product/product*(
                derivative*lower-(points[n-1]-x0)*table[derivative][n-1][n-1])
        previous_product=product
    return table


def apply_finite_diff(order, x_list, y_list, x0=0):
    """Apply arbitrary-grid finite difference weights to scalar samples."""
    from .. import S
    points=_points(x_list)
    samples=list(itertools.islice(iter(y_list),MAX_POINTS+1))
    if len(samples)!=len(points):
        raise ValueError('Grid points and sample values must have equal lengths')
    order=_order(order)
    weights=finite_diff_weights(order,points,x0)[order][-1]
    return sum((weight*_scalar(value) for weight,value in zip(weights,samples)),S.Zero)


def _point_sequence(points):
    if isinstance(points,_core.Expr):
        return points.args if getattr(points.func,'__name__','')=='Tuple' else None
    if isinstance(points,(str,bytes)):
        return None
    try:
        return iter(points)
    except TypeError:
        return None


def _map_derivatives(expression, points, x0, wrt):
    if expression.is_Derivative:
        return _as_finite_diff(expression,points,x0,wrt)
    replacements={}
    pending=list(expression.args)
    while pending:
        node=pending.pop()
        if not isinstance(node,_core.Expr):
            continue
        if node.is_Derivative:
            replacements[node]=_as_finite_diff(node,points,x0,wrt)
        else:
            pending.extend(node.args)
    return expression.xreplace(replacements) if replacements else expression


def _as_finite_diff(derivative, points=1, x0=None, wrt=None):
    """Convert a derivative, or derivative nodes in an expression, to a stencil."""
    from .. import S, Derivative
    derivative=_scalar(derivative)
    sequence=_point_sequence(points)
    if sequence is not None:
        points=_points(sequence)
    if not derivative.is_Derivative:
        return _map_derivatives(derivative,points,x0,wrt)
    counts=derivative.variable_count
    if wrt is None:
        result=derivative
        for variable,_ in counts:
            result=_map_derivatives(result,points,x0,variable)
        return result
    wrt=_scalar(wrt)
    order=sum(count for variable,count in counts if variable==wrt)
    if not order:
        return derivative
    order=_order(order)
    center=wrt if x0 is None else _scalar(x0)
    if sequence is None:
        step=_scalar(points)
        if step.is_zero is True:
            raise ValueError('Finite difference step must be nonzero')
        if step.is_Function and step.has(wrt):
            step=step.subs(wrt,center)
        grid=[center+step*S(i)/2 for i in range(-order,order+1,2)]
    else:
        grid=points
    if len(grid)<order+1:
        raise ValueError('Too few points for derivative order %d' % order)
    other_variables=[(variable,count) for variable,count in counts if variable!=wrt]
    samples=[]
    for point in grid:
        value=derivative.expr.subs(wrt,point)
        if other_variables:
            value=Derivative(value,*other_variables,evaluate=False)
        samples.append(value)
    return apply_finite_diff(order,grid,samples,center)


def _convert_substituted_derivatives(expression, points):
    """Discretize derivatives evaluated at shifted stencil positions."""
    from .. import Dummy
    pending=[expression]
    replacements={}
    while pending:
        node=pending.pop()
        if not isinstance(node,_core.Expr):
            continue
        if node.is_Subs:
            variables,values=node.args[1:]
            variables=variables.args if isinstance(variables,_core.Expr) else variables
            values=values.args if isinstance(values,_core.Expr) else values
            if len(variables)!=1:
                raise NotImplementedError('Finite differences of multivariate Subs require explicit nesting')
            variable,point=variables[0],values[0]
            inner,_=_convert_substituted_derivatives(node.expr,points)
            # Protect stencil centers from a second substitution. In
            # Subs(Derivative(f(x),x)**2,x,x+h), replacing the derivative
            # first and then x would otherwise shift each stencil twice.
            protected,stencils={},{}
            children=[inner]
            while children:
                child=children.pop()
                if not isinstance(child,_core.Expr):
                    continue
                if child.is_Derivative:
                    token=Dummy('stencil')
                    protected[child]=token
                    if any(v==variable for v,_ in child.variable_count):
                        stencils[token]=_as_finite_diff(child,points,point,variable)
                    else:
                        stencils[token]=child.subs(variable,point)
                else:
                    children.extend(child.args)
            converted=inner.xreplace(protected).subs(variable,point).xreplace(stencils)
            if converted!=node:
                replacements[node]=converted
        else:
            pending.extend(node.args)
    result=expression.xreplace(replacements) if replacements else expression
    # A shifted derivative can itself contain a shifted higher derivative.
    return result,replacements


def differentiate_finite(expr, *symbols, points=1, x0=None, wrt=None, evaluate=False):
    """Differentiate an expression using finite difference stencils."""
    from .. import Derivative
    expr=_scalar(expr)
    if not symbols:
        variables=expr.free_symbols
        if len(variables)!=1:
            raise ValueError('Supply differentiation variables explicitly')
        symbols=(next(iter(variables)),)
    if evaluate:
        pending=[expr]
        while pending:
            node=pending.pop()
            if node.is_Derivative:
                evaluate=False
                break
            pending.extend(child for child in node.args if isinstance(child,_core.Expr))
    derivative=Derivative(expr,*symbols,evaluate=bool(evaluate))
    result=_as_finite_diff(derivative,points,x0,wrt)
    for _ in range(MAX_POINTS):
        result,replacements=_convert_substituted_derivatives(result,points)
        if not replacements:
            return result
    raise NotImplementedError('Nested stencil substitution exceeds 256 levels')


def install():
    _core.Expr.as_finite_difference=_as_finite_diff


__all__=['finite_diff_weights','apply_finite_diff','differentiate_finite']
