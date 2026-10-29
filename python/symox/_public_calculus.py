"""Piecewise algebra, exact B-splines and variational/line integration."""
from functools import lru_cache
from itertools import product
import operator


def _pairs(expr):
    arguments=expr.args
    if arguments and isinstance(arguments[0],(tuple,list)):return arguments
    return tuple(zip(arguments[::2],arguments[1::2]))


def _bound_piecewise_variables(expr):
    """None unless ``expr`` is a Sum, Product or Integral; otherwise whether a
    Piecewise condition in its function depends on one of its bound variables."""
    from . import Sum, Product, Integral, Piecewise
    if not isinstance(expr, (Sum, Product, Integral)):
        return None
    variables = set()
    for limit in expr.args[1:]:
        items = limit.args if hasattr(limit, 'args') and limit.args else (limit,)
        variables.add(items[0])
    for piece in expr.args[0].atoms(Piecewise):
        for _, condition in _pairs(piece):
            if getattr(condition, 'free_symbols', set()) & variables:
                return True
    return False


def _has_piecewise(expr):
    """``expr.has(Piecewise)``; atoms(Function) is much cheaper than the
    structural has() query when there is no Piecewise."""
    from . import Piecewise, Function
    try:
        return any(isinstance(a, Piecewise) for a in expr.atoms(Function))
    except NotImplementedError:
        return expr.has(Piecewise)


def piecewise_fold(expr, evaluate=True):
    """SymPy's ``piecewise_fold``: move Piecewise outward through the
    surrounding function, collecting branches with identical conditions."""
    from . import Piecewise, And, ITE, sympify, true, simplify_logic
    from ._core import Expr
    from .core.sorting import ordered
    from .utilities.iterables import sift, common_prefix
    expr = sympify(expr)
    if not isinstance(expr, Expr) or not _has_piecewise(expr):
        return expr
    # an unevaluated Piecewise keeps duplicate conditions
    held = evaluate if evaluate is not None else False
    new_args = []
    if expr.is_Piecewise or type(expr).__name__ == 'ExprCondPair':
        for e, c in _pairs(expr):
            if not e.is_Piecewise:
                e = piecewise_fold(e)
            if isinstance(c, ITE):
                c = c.to_nnf()
                c = simplify_logic(c, form='cnf')
            if e.is_Piecewise:
                new_args.extend([(piecewise_fold(ei), And(ci, c))
                    for ei, ci in _pairs(e)])
            else:
                new_args.append((e, c))
    else:
        bound = _bound_piecewise_variables(expr)
        if bound is not None:
            # A condition on a summation or integration variable cannot
            # leave its binder: Sum(Piecewise((f, x > 1), (0, True)), (x, ..))
            # is not Piecewise((Sum(f, (x, ..)), x > 1), ...). The reference
            # folds it anyway (a wrong value); keep the binder intact.
            if bound:
                return expr.func(piecewise_fold(expr.args[0], evaluate), *expr.args[1:])
        # f(P1, P2) with branch conditions in common keeps them as the
        # leading conditions when f is a commutative Add or Mul
        if expr.is_Add or expr.is_Mul and expr.is_commutative:
            p, args = sift(expr.args, lambda x: x.is_Piecewise, binary=True)
            pc = sift(p, lambda x: tuple([c for e, c in _pairs(x)]))
            for c in list(ordered(pc)):
                if len(pc[c]) > 1:
                    pargs = [list(_pairs(i)) for i in pc[c]]
                    com = common_prefix(*[
                        [cond for _, cond in j] for j in pargs])
                    n = len(com)
                    collected = []
                    for i in range(n):
                        collected.append((
                            expr.func(*[ai[i][0] for ai in pargs]),
                            com[i]))
                    remains = []
                    for a in pargs:
                        if n == len(a):
                            continue
                        if a[n][1] == True:  # noqa: E712
                            remains.append(a[n][0])
                        else:
                            remains.append(
                                Piecewise(*a[n:], evaluate=False))
                    if remains:
                        collected.append((expr.func(*remains), True))
                    args.append(Piecewise(*collected, evaluate=False))
                    continue
                args.extend(pc[c])
        else:
            args = expr.args
        folded = list(map(piecewise_fold, args))
        for ec in product(*[
                (_pairs(i) if isinstance(i, Expr) and i.is_Piecewise else
                 [(i, true)]) for i in folded]):
            e, c = zip(*ec)
            new_args.append((expr.func(*e), And(*c)))
    if evaluate is None:
        # don't return duplicate conditions, otherwise don't evaluate
        new_args = list(reversed([(e, c) for c, e in {
            c: e for e, c in reversed(new_args)}.items()]))
    rv = Piecewise(*new_args, evaluate=held)
    if evaluate is None and len(rv.args) == 1 and _pairs(rv)[0][1] == True:  # noqa: E712
        return _pairs(rv)[0][0]
    if not (isinstance(rv, Expr) and rv.is_Piecewise):
        return rv
    if any(e.has(Piecewise) for q in rv.atoms(Piecewise) for e, _ in _pairs(q)):
        return piecewise_fold(rv)
    return rv


def piecewise_exclusive(expr, *, skip_nan=False, deep=True):
    """SymPy's ``piecewise_exclusive``: each condition excludes the earlier
    ones, simplified with ``Boolean.simplify``."""
    from . import Piecewise, And, Or, Not, sympify, false, true, nan
    expr = sympify(expr)

    def make_exclusive(*pwargs):
        pwargs = [tuple(a) for a in pwargs]
        cumcond = false
        newargs = []
        for expr_i, cond_i in pwargs[:-1]:
            cancond = And(cond_i, Not(cumcond)).simplify()
            cumcond = Or(cond_i, cumcond).simplify()
            newargs.append((expr_i, cancond))
        expr_n, cond_n = pwargs[-1]
        cancond_n = And(cond_n, Not(cumcond)).simplify()
        newargs.append((expr_n, cancond_n))
        if not skip_nan:
            cumcond = Or(cond_n, cumcond).simplify()
            if cumcond is not true and cumcond != true:
                newargs.append((nan, Not(cumcond).simplify()))
        return Piecewise(*newargs, evaluate=False)

    if deep:
        return expr.replace(Piecewise, make_exclusive)
    if getattr(expr, 'is_Piecewise', False):
        return make_exclusive(*_pairs(expr))
    return expr


def _diff_wrt(e):
    """``Basic._diff_wrt`` of the reference: Symbols and applied undefined
    functions."""
    flag = getattr(e, '_diff_wrt', None)
    if flag is not None:
        return flag
    if getattr(e, 'is_Symbol', False):
        return True
    from .core.function import AppliedUndef
    return isinstance(e, AppliedUndef)


def _blessed(e):
    from . import Rational, NumberSymbol
    return _diff_wrt(e.lhs) and (
        _diff_wrt(e.rhs) or
        isinstance(e.rhs, (Rational, NumberSymbol)))


def _piecewise_simplify_eq_and(args):
    from . import And, Eq
    from .core.sorting import ordered
    from .utilities.iterables import sift
    for i, (expr, cond) in enumerate(args):
        if isinstance(cond, And):
            eqs, other = sift(cond.args, lambda i: isinstance(i, Eq), binary=True)
        elif isinstance(cond, Eq):
            eqs, other = [cond], []
        else:
            eqs = other = []
        if eqs:
            eqs = list(ordered(eqs))
            for j, e in enumerate(eqs):
                if _blessed(e):
                    expr = expr.subs(*e.args)
                    eqs[j + 1:] = [ei.subs(*e.args) for ei in eqs[j + 1:]]
                    other = [ei.subs(*e.args) for ei in other]
            cond = And(*(eqs + other))
            args[i] = (expr, cond)
    return args


def _piecewise_simplify_equal_to_next_segment(args):
    from . import And, Eq
    from .core.sorting import ordered
    from .utilities.iterables import sift
    prevexpr = None
    for i, (expr, cond) in reversed(list(enumerate(args))):
        if prevexpr is not None:
            if isinstance(cond, And):
                eqs, other = sift(cond.args, lambda i: isinstance(i, Eq), binary=True)
            elif isinstance(cond, Eq):
                eqs, other = [cond], []
            else:
                eqs = other = []
            _prevexpr = prevexpr
            _expr = expr
            if eqs and not other:
                eqs = list(ordered(eqs))
                for e in eqs:
                    if len(args) == 2 or _blessed(e):
                        _prevexpr = _prevexpr.subs(*e.args)
                        _expr = _expr.subs(*e.args)
            if _prevexpr == _expr:
                args[i] = (args[i + 1][0], cond)
            else:
                prevexpr = expr
        else:
            prevexpr = expr
    return args


def piecewise_simplify_arguments(expr, **kwargs):
    """SymPy's ``piecewise_simplify_arguments``.  The conditions come from
    ``Piecewise._intervals`` on the whole real line, so a periodic
    condition is rewritten only when it is never or always true (the
    reference keeps the solutions of one period, which is wrong)."""
    from . import S, Piecewise, Eq, And, Interval, Basic, nan, simplify
    from . import _core
    pairs = list(_pairs(expr))
    f1 = pairs[0][1].free_symbols
    args = None
    if len(f1) == 1 and not expr.atoms(Eq):
        x = f1.pop()
        ok, abe_ = _core._piecewise_intervals(expr, x, True)

        def include(c, x, a):
            try:
                return c.subs(x, a) == True  # noqa: E712
            except TypeError:
                return False
        if ok:
            args = []
            covered = S.EmptySet
            for a, b, e, i in abe_:
                c = pairs[i][1]
                incl_a = include(c, x, a)
                incl_b = include(c, x, b)
                iv = Interval(a, b, not incl_a, not incl_b)
                cset = iv - covered
                if not cset:
                    continue
                try:
                    a = cset.inf
                except NotImplementedError:
                    pass
                else:
                    incl_a = include(c, x, a)
                if incl_a and incl_b:
                    if a.is_infinite and b.is_infinite:
                        c = S.true
                    elif b.is_infinite:
                        c = (x > a) if a in covered else (x >= a)
                    elif a.is_infinite:
                        c = (x <= b)
                    elif a in covered:
                        c = And(a < x, x <= b)
                    else:
                        c = And(a <= x, x <= b)
                elif incl_a:
                    if a.is_infinite:
                        c = (x < b)
                    elif a in covered:
                        c = And(a < x, x < b)
                    else:
                        c = And(a <= x, x < b)
                elif incl_b:
                    if b.is_infinite:
                        c = (x > a)
                    else:
                        c = And(a < x, x <= b)
                else:
                    if a in covered:
                        c = (x < b)
                    else:
                        c = And(a < x, x < b)
                covered |= iv
                if a is S.NegativeInfinity and incl_a:
                    covered |= {S.NegativeInfinity}
                if b is S.Infinity and incl_b:
                    covered |= {S.Infinity}
                args.append((e, c))
            if not S.Reals.is_subset(covered):
                args.append((nan, True))
    if args is None:
        args = list(pairs)
        for i in range(len(args)):
            e, c = args[i]
            if isinstance(c, Basic):
                c = simplify(c, **kwargs)
            args[i] = (e, c)
    doit = kwargs.pop('doit', None)
    for i in range(len(args)):
        e, c = args[i]
        if isinstance(e, Basic):
            newe = simplify(e, doit=False, **kwargs)
            if newe != e:
                e = newe
        args[i] = (e, c)
    if doit is not None:
        kwargs['doit'] = doit
    return Piecewise(*args)


def piecewise_simplify(expr, **kwargs):
    """SymPy's ``piecewise_simplify``."""
    from . import Piecewise
    expr = piecewise_simplify_arguments(expr, **kwargs)
    if not getattr(expr, 'is_Piecewise', False):
        return expr
    args = [tuple(a) for a in _pairs(expr)]
    args = _piecewise_simplify_eq_and(args)
    args = _piecewise_simplify_equal_to_next_segment(args)
    return Piecewise(*args)


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
    from . import sympify, Dummy, Symbol, Eq, simplify
    from ._core import Expr
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


def _piecewise_integrate(self, x, **kwargs):
    """``Piecewise.piecewise_integrate``: each piece integrated on its own."""
    from . import Piecewise, integrate
    return Piecewise(*[(integrate(e, x, **kwargs), c) for e, c in _pairs(self)])


def _as_expr_set_pairs(self, domain=None):
    """``Piecewise.as_expr_set_pairs``: each expression with the part of
    ``domain`` (the reals by default) where its condition first holds."""
    from . import S, Eq, Ne
    from .core.relational import Relational
    from ._public_tools import filldedent
    if domain is None:
        domain = S.Reals
    exp_sets = []
    U = domain
    complex = not domain.is_subset(S.Reals)
    cond_free = set()
    for expr, cond in _pairs(self):
        cond_free |= cond.free_symbols
        if len(cond_free) > 1:
            raise NotImplementedError(filldedent('''
                multivariate conditions are not handled.'''))
        if complex:
            for i in cond.atoms(Relational):
                if not isinstance(i, (Eq, Ne)):
                    raise ValueError(filldedent('''
                        Inequalities in the complex domain are
                        not supported. Try the real domain by
                        setting domain=S.Reals'''))
        cond_int = U.intersect(cond.as_set())
        U = U - cond_int
        if cond_int != S.EmptySet:
            exp_sets.append((expr, cond_int))
    return exp_sets


def _install_piecewise_methods():
    from ._node_types import node_class
    cls = node_class('Piecewise')
    if cls is not None:
        type.__setattr__(cls, 'piecewise_integrate', _piecewise_integrate)
        type.__setattr__(cls, 'as_expr_set_pairs', _as_expr_set_pairs)


_install_piecewise_methods()
