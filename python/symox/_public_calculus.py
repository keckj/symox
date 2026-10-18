"""Piecewise algebra, exact B-splines and variational/line integration."""
from functools import lru_cache
from itertools import product
import operator


def _pairs(expr):
    return tuple(zip(expr.args[::2],expr.args[1::2]))


def piecewise_fold(expr,evaluate=True):
    from . import Expr,Piecewise,And,Or,sympify,true
    expr = sympify(expr)
    def visit(node):
        if not isinstance(node,Expr) or not node.args:
            return node
        if node.is_Piecewise:
            pieces = []
            for value,condition in _pairs(node):
                value,condition = visit(value),visit(condition)
                if condition.is_Piecewise:
                    condition = Or(*(And(c,v) for v,c in _pairs(condition)))
                if value.is_Piecewise:
                    pieces.extend((v,And(condition,c)) for v,c in _pairs(value))
                else:
                    pieces.append((value,condition))
            return Piecewise(*pieces,evaluate=evaluate)
        args = tuple(visit(a) for a in node.args)
        choices = [_pairs(a) if isinstance(a,Expr) and a.is_Piecewise else ((a,true),) for a in args]
        if any(isinstance(a,Expr) and a.is_Piecewise for a in args):
            return Piecewise(*((node.func(*(v for v,c in combination)),And(*(c for v,c in combination)))
                               for combination in product(*choices)),evaluate=evaluate)
        return node if args == node.args else node.func(*args)
    return visit(expr)


def piecewise_exclusive(expr,*,skip_nan=False,deep=True):
    from . import Expr,Piecewise,And,Or,Not,sympify,false,nan,simplify_logic
    expr = sympify(expr)
    def visit(node):
        if not isinstance(node,Expr) or not node.args:
            return node
        if node.is_Piecewise:
            pieces,covered = [],false
            for value,condition in _pairs(node):
                exclusive = simplify_logic(And(condition,Not(covered)))
                if exclusive != false:
                    pieces.append((visit(value) if deep else value,exclusive))
                covered = Or(covered,condition)
            remainder = simplify_logic(Not(covered))
            if not skip_nan and remainder != false:
                pieces.append((nan,remainder))
            return Piecewise(*pieces,evaluate=False)
        if not deep:
            return node
        args = tuple(visit(a) for a in node.args)
        return node if args == node.args else node.func(*args)
    return visit(expr)


def _knots(values):
    from . import sympify
    values = values.args if getattr(values,'is_Tuple',False) else values
    result = tuple(map(sympify,values))
    if any(v.is_real is False for v in result):
        raise ValueError('spline knots must be real')
    for a,b in zip(result,result[1:]):
        if not bool(a <= b):
            raise ValueError('spline knots must be nondecreasing')
    return result


def _basis_pieces(degree,knots,index,x):
    from . import Integer,expand
    @lru_cache(maxsize=None)
    def basis(d,i,span):
        if d == 0:
            return Integer(int(i == span))
        left = knots[i+d]-knots[i]
        right = knots[i+d+1]-knots[i+1]
        result = Integer(0)
        if left != 0:
            result += (x-knots[i])*basis(d-1,i,span)/left
        if right != 0:
            result += (knots[i+d+1]-x)*basis(d-1,i+1,span)/right
        return expand(result)
    return [(basis(degree,index,span),knots[span],knots[span+1])
            for span in range(index,index+degree+1) if knots[span] < knots[span+1]]


def bspline_basis(d,knots,n,x):
    from . import sympify,Piecewise,And
    d,n = operator.index(d),operator.index(n)
    knots,x = _knots(knots),sympify(x)
    if d < 0 or n < 0 or n+d+1 >= len(knots):
        raise ValueError('invalid spline degree or basis index')
    pieces = [(value,And(x>=a,x<=b)) for value,a,b in _basis_pieces(d,knots,n,x) if value != 0]
    return Piecewise(*pieces,(0,True))


def bspline_basis_set(d,knots,x):
    d = operator.index(d)
    knots = _knots(knots)
    if d < 0 or len(knots) < d+2:
        raise ValueError('not enough knots for the spline degree')
    return [bspline_basis(d,knots,i,x) for i in range(len(knots)-d-1)]


def interpolating_spline(d,x,X,Y):
    from . import Matrix,Piecewise,And,sympify,expand
    d = operator.index(d)
    x,X,Y = sympify(x),tuple(map(sympify,X)),tuple(map(sympify,Y))
    if d < 1 or len(X) != len(Y) or len(X) < d+1 or any(not bool(a < b) for a,b in zip(X,X[1:])):
        raise ValueError('spline degree and strictly increasing interpolation data are incompatible')
    if d % 2:
        j = (d+1)//2
        interior = X[j:-j]
    else:
        j = d//2
        interior = tuple((a+b)/2 for a,b in zip(X[j:-j-1],X[j+1:-j]))
    knots = (X[0],)*(d+1)+tuple(interior)+(X[-1],)*(d+1)
    basis = bspline_basis_set(d,knots,x)
    matrix = Matrix([[b.subs(x,v) for b in basis] for v in X])
    coefficients = matrix.inv()*Matrix(Y)
    polynomials = [_basis_pieces(d,knots,i,x) for i in range(len(basis))]
    intervals = sorted({(a,b) for pieces in polynomials for value,a,b in pieces},key=lambda pair:float(pair[0]))
    result = []
    for a,b in intervals:
        value = sum((coefficients[i]*next((v for v,l,r in pieces if (l,r)==(a,b)),sympify(0)) for i,pieces in enumerate(polynomials)),sympify(0))
        result.append((expand(value),And(x>=a,x<=b)))
    return Piecewise(*result)


def line_integrate(field,curve,vars):
    from . import sympify,integrate,sqrt,simplify,trigsimp
    from .geometry import Curve
    if not isinstance(curve,Curve):
        raise ValueError('line integration requires a geometry Curve')
    variables = tuple(map(sympify,vars))
    coordinates = curve.functions
    coordinates = coordinates.args if getattr(coordinates,'is_Tuple',False) else tuple(coordinates)
    limits = curve.limits
    limits = limits.args if getattr(limits,'is_Tuple',False) else tuple(limits)
    parameter = limits[0]
    if len(variables) != len(coordinates) or parameter in variables:
        raise ValueError('field coordinates and curve parameter are incompatible')
    value = sympify(field).subs(dict(zip(variables,coordinates)),simultaneous=True)
    metric = sqrt(trigsimp(sum(c.diff(parameter)**2 for c in coordinates)))
    return integrate(simplify(value*metric),limits)


def euler_equations(L,funcs=(),vars=()):
    from . import sympify,Expr,Dummy,Symbol,Eq,simplify
    L = sympify(L)
    def sequence(values):
        if isinstance(values,(tuple,list)):
            return tuple(values)
        if getattr(values,'is_Tuple',False):
            return values.args
        return (values,)
    functions = tuple(map(sympify,sequence(funcs)))
    nodes,seen,stack = [],set(),[L]
    while stack:
        node = stack.pop()
        if not isinstance(node,Expr) or node in seen:
            continue
        seen.add(node)
        nodes.append(node)
        stack.extend(node.args)
    if not functions:
        functions = tuple(sorted((node for node in nodes if node.is_UndefinedFunction),key=str))
    if not functions or any(not f.is_Function for f in functions):
        raise TypeError('supply dependent applied functions for the Lagrangian')
    variables = tuple(map(sympify,sequence(vars))) or functions[0].args
    if any(not isinstance(v,Symbol) for v in variables) or any(f.args != variables for f in functions):
        raise ValueError('independent variables must match every dependent function')
    jets = {f:() for f in functions}
    for node in nodes:
        if node.is_Derivative and node.expr in functions:
            jets[node] = tuple(node.variables)
    replacements = {jet:Dummy('jet') for jet in jets}
    inverse = {value:jet for jet,value in replacements.items()}
    algebraic = L.xreplace(replacements)
    equations = []
    for function in functions:
        residual = sympify(0)
        for jet,derivatives in jets.items():
            if jet != function and (not jet.is_Derivative or jet.expr != function):
                continue
            term = algebraic.diff(replacements[jet]).xreplace(inverse)
            if derivatives:
                term = term.diff(*derivatives)
            residual += (-1)**len(derivatives)*term
        equation = Eq(simplify(residual),0)
        if equation.is_Relational:
            equations.append(equation)
    return equations
