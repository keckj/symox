"""First-order Lie point symmetries, proved by the determining equation.

Heuristics follow the determining-equation ansatzes used by SymPy. Every
returned generator is independently substituted into that equation. All
arithmetic, integration, polynomial extraction and nullspaces are native.
"""
from .._differential import equation, dependent, nodes, linear_data
from ..deutils import ode_order

lie_heuristics = ('abaco1_simple','abaco1_product','abaco2_similar',
                  'abaco2_unique_unknown','abaco2_unique_general','linear',
                  'function_sum','bivariate','chi')


def _simplify(value):
    from ... import simplify, cancel, powsimp, trigsimp, exp
    value = simplify(trigsimp(powsimp(cancel(value))))
    # exp(r)**z = exp(r*z) for real r, since exp(r) is positive. Keep
    # genuinely complex exponential bases on their principal-power branch.
    replacements = {node:exp(node.base.args[0]*node.exp) for node in nodes(value)
                    if node.is_Pow and getattr(node.base.func,'__name__','')=='exp'
                    and node.base.args[0].is_real is True}
    return value.xreplace(replacements) if replacements else value


def _zero(value):
    value = _simplify(value)
    return value == 0 or value.is_zero is True


def _illegal(value):
    from ... import S, oo
    return value.has(S.NaN,S.ComplexInfinity,oo,-oo)


def _D(value,variable,count=1,_depth=0):
    """Native differentiation with the chain rule for held substitutions.

    The scalar engine intentionally leaves some derivatives of Subs held.
    Symmetry tests need their explicit chain rule, especially for unknown
    forcing functions evaluated at a characteristic expression.
    """
    from ... import diff, Dummy, Subs, S
    if _depth>64 or count>256:
        raise NotImplementedError('Symmetry differentiation exceeds its structural bounds')
    for _ in range(count):
        raw = diff(value,variable)
        replacements = {}
        for original in nodes(raw):
            if not original.is_Derivative:
                continue
            node = original.xreplace(replacements)
            if not node.expr.is_Subs:
                continue
            result = node.expr
            for coordinate,power in node.variable_count:
                if int(power)>256:
                    raise NotImplementedError('Symmetry derivative order exceeds 256')
                for _ in range(int(power)):
                    if not result.is_Subs:
                        result = _D(result,coordinate,_depth=_depth+1)
                        continue
                    body = result.expr
                    variables,points = tuple(result.variables),tuple(result.point)
                    if any(not bound.is_Symbol for bound in variables):
                        raise NotImplementedError('Symmetry substitution variables must be symbols')
                    if coordinate in variables:
                        renaming = {bound:Dummy('lie_bound_coordinate') for bound in variables}
                        body = body.xreplace(renaming)
                        variables = tuple(renaming[bound] for bound in variables)
                    derivative = _D(body,coordinate,_depth=_depth+1)
                    total = S.Zero if derivative==0 else Subs(derivative,variables,points)
                    for bound,point in zip(variables,points):
                        scale = _D(point,coordinate,_depth=_depth+1)
                        if scale==0:
                            continue
                        partial = _D(body,bound,_depth=_depth+1)
                        if partial!=0:
                            total += scale*Subs(partial,variables,points)
                    result = total
            replacements[original] = result
        value = raw.xreplace(replacements)
    return value


def _problem(eq, func=None, order=None, match=None):
    from ... import Dummy, Derivative, solve, sympify
    expression = equation(eq)
    func = dependent(expression,func)
    actual_order = ode_order(expression,func)
    if actual_order != 1 or order not in (None,0,1):
        raise NotImplementedError("Infinitesimals for only first order ODE's have been implemented")
    x = func.args[0]
    if match is not None:
        try:
            y,h = sympify(match['y']),sympify(match['h'])
        except (TypeError,KeyError) as error:
            raise ValueError('A Lie match must supply the dependent coordinate and slope') from error
        if not y.is_Symbol or y == x:
            raise ValueError('The dependent coordinate must be a distinct symbol')
    else:
        y = Dummy('lie_dependent_coordinate')
        derivative = Derivative(func,x)
        result = linear_data(expression,[derivative])
        if result is not None and not _zero(result[0][0]):
            h = -result[1]/result[0][0]
        else:
            placeholder = Dummy('lie_derivative')
            roots = solve(expression.xreplace({derivative:placeholder}),placeholder)
            if not roots:
                raise NotImplementedError('Cannot isolate a first derivative for symmetry analysis')
            h = roots[0]
        h = h.xreplace({func:y})
    h = _simplify(h)
    if _illegal(h) or h.has(func):
        raise NotImplementedError('Cannot isolate a scalar first-order slope')
    if any(node.is_Derivative and node.expr.has(y) for node in nodes(h)):
        raise NotImplementedError('The isolated slope still contains a dependent derivative')
    return expression,func,x,y,h


def _residual(h, x, y, xi, eta):
    from ... import powsimp
    value = _simplify(_D(eta,x)+(_D(eta,y)-_D(xi,x))*h
                      -_D(xi,y)*h*h-xi*_D(h,x)-eta*_D(h,y))
    return value if value==0 else _simplify(powsimp(value.expand()))


def _symmetry_checkinfsol(eq, infinitesimals, func=None, order=None):
    """Check each supplied generator by its exact determining residual."""
    from ... import Function, sympify
    expression,func,x,y,h = _problem(eq,func,order)
    xi,eta = Function('xi')(x,func),Function('eta')(x,func)
    results = []
    for candidate in infinitesimals:
        if xi not in candidate or eta not in candidate:
            raise ValueError('Each generator requires xi(x, func) and eta(x, func) entries')
        left = sympify(candidate[xi]).xreplace({func:y})
        right = sympify(candidate[eta]).xreplace({func:y})
        residual = _residual(h,x,y,left,right)
        valid = not _illegal(residual) and _zero(residual)
        results.append((valid,residual.xreplace({y:func})))
    return results


def _exp_integral(value, variable):
    from ... import exp, integrate
    return _simplify(exp(integrate(_simplify(value),variable)))


def _simple(h,x,y):
    hx,hy = _D(h,x),_D(h,y)
    if not hy.has(y):
        yield 0,_exp_integral(hy,x)
    if _zero(h):
        yield 1,0
        return
    factor = _simplify(hy/h)
    if not factor.has(x):
        yield 0,_exp_integral(factor,y)
    factor = _simplify(-hx/h)
    if not factor.has(y):
        yield _exp_integral(factor,x),0
    factor = _simplify(-hx/(h*h))
    if not factor.has(x):
        yield _exp_integral(factor,y),0


def _inverse(h,x,y):
    return _simplify((1/h).xreplace({x:y,y:x}))


def _product(h,x,y):
    from ... import separatevars
    if _zero(h):
        return
    for inverted in (False,True):
        slope = _inverse(h,x,y) if inverted else h
        # The mixed logarithmic derivative is evaluated rationally, avoiding
        # needless logarithm branch rewrites.
        factor = _simplify(_D(_D(slope,y)/slope,x)/slope**2)
        separated = separatevars(factor,symbols=[x,y],dict=True)
        if separated is None or _zero(separated['coeff']):
            continue
        first = separated[x]
        integrand = _simplify(first*_D(1/(first*slope),x))
        if integrand.has(x):
            continue
        candidate = _simplify(first*_exp_integral(integrand,y))
        yield (0,candidate.xreplace({x:y,y:x})) if inverted else (candidate,0)


def _similar(h,x,y):
    from ... import exp, S
    if _zero(h):
        return
    # Proportional coordinate speeds are a frequent special case of this
    # ansatz. Recover their ratio from invariant function arguments before
    # forming the much larger second-derivative invariants.
    hx,hy = _D(h,x),_D(h,y)
    ratios = []
    for node in nodes(h):
        if not node.is_Function and not node.is_UndefinedFunction:
            continue
        for argument in node.args:
            if not argument.has(x,y) or _zero(_D(argument,y)):
                continue
            ratio = _simplify(-_D(argument,x)/_D(argument,y))
            if not ratio.has(x,y) and ratio not in ratios:
                ratios.append(ratio)
    for ratio in ratios+[S.Zero,S.One,-S.One]:
        if _zero(ratio-h):
            continue
        numerator = _simplify(hx+ratio*hy)
        rate_x = _simplify(numerator/(ratio-h))
        if not rate_x.has(y):
            speed = _exp_integral(rate_x,x)
            yield speed,ratio*speed
        rate_y = _simplify(numerator/(h*(ratio-h)))
        if not rate_y.has(x):
            speed = _exp_integral(rate_y,y)
            yield speed,ratio*speed
    for inverted in (False,True):
        slope = _inverse(h,x,y) if inverted else h
        hy,hyy = _D(slope,y),_D(slope,y,2)
        if _zero(hyy):
            continue
        ratio = _simplify(hy/hyy)
        dx,dy = _D(ratio,x),_D(ratio,y)
        if not ratio.has(x,y) and not _zero(ratio):
            # h=A(x)+B(x)*exp(y/ratio). Recover A and B by differentiation
            # instead of pattern matching arbitrary-function Wild objects.
            A = _simplify(slope-ratio*hy)
            B = _simplify(ratio*hy*exp(-y/ratio))
            if A.has(y) or B.has(y) or _zero(B):
                continue
            xi = _simplify(_exp_integral(-A/ratio,x)/B)
            eta = _simplify(A*xi)
        elif not _zero(dy):
            gamma = _simplify(dx/dy)
            if gamma.has(y) or _zero(slope+gamma):
                continue
            factor = _simplify((gamma*hy-_D(gamma,x)-_D(slope,x))/(slope+gamma))
            if factor.has(y):
                continue
            xi = _exp_integral(factor,x)
            eta = _simplify(-gamma*xi)
        else:
            continue
        if inverted:
            yield eta.xreplace({x:y,y:x}),xi.xreplace({x:y,y:x})
        else:
            yield xi,eta


def _unique_unknown(h,x,y):
    from ... import separatevars, S
    functions = [node for node in nodes(h) if (node.is_UndefinedFunction
                 or node.is_Pow and node.exp.is_integer is not True)
                 and node.has(x) and node.has(y)]
    for value in functions:
        if _zero(_D(value,x)):
            continue
        ratio = _simplify(_D(value,y)/_D(value,x))
        if _zero(ratio):
            continue
        separated = separatevars(ratio,symbols=[x,y],dict=True)
        if separated is not None and not _zero(separated['coeff']):
            xi = separated[x]
            eta = -1/(separated[y]*separated['coeff'])
            yield xi,eta
            if not _zero(xi) and not _zero(eta):
                yield 1/eta,1/xi
        else:
            yield S.One,-1/ratio
            yield -ratio,S.One


def _unique_general(h,x,y):
    hx,hy = _D(h,x),_D(h,y)
    A,B,C = _D(hx,y),_D(hy,y)+hy**2,_D(hx,x)-hx**2
    if any(_zero(value) for value in (A,B,C)):
        return
    Ax,Ay = _D(A,x),_D(A,y)
    D = _simplify((2*_D(Ax,y)+hx*Ay-Ax*hy+(hx*hy+2*A)*A)*A-3*Ax*Ay)
    if _zero(D):
        E1 = _simplify(3*Ax**2+((hx**2+2*C)*A-2*_D(Ax,x))*A)
        if _zero(E1):
            return
        E2 = _simplify((2*_D(Ay,y)+(2*B-hy**2)*A)*A-3*Ay**2)
        E3 = _simplify(E1*((28*Ax+4*hx*A)*A**3-E1*(hy*A+Ay))-_D(E1,x)*8*A**4)
        if not _zero(E2) or not _zero(E3):
            return
        factor = _simplify((4*A**3*(Ax-hx*A)+E1*(hy*A-Ay))/(2*A*E1))
        if factor.has(x):
            return
        eta = _exp_integral(factor,y)
        xi = _simplify(-4*A**3*eta/E1)
    else:
        E1 = _simplify((2*_D(Ay,y)+(2*B-hy**2)*A)*A-3*Ay**2)
        if _zero(E1):
            return
        E2 = _simplify(4*A**3*D-D**2+E1*((2*_D(Ax,x)-(hx**2+2*C)*A)*A-3*Ax**2))
        E3 = _simplify(-(A*D)*_D(E1,y)+((_D(E1,x)-hy*D)*A+3*Ay*D+(A*hx-3*Ax)*E1)*E1)
        if not _zero(E2) or not _zero(E3):
            return
        factor = _simplify(((A*hx-Ax)*E1-(Ay+A*hy)*D)/(2*A*D))
        if factor.has(x):
            return
        eta = _exp_integral(factor,y)
        xi = _simplify(-E1*eta/D)
    if not xi.has(y):
        yield xi,eta


def _function_sum(h,x,y):
    from ... import separatevars, Dummy, integrate, solve
    if _zero(h):
        return
    for inverted in (False,True):
        slope = _inverse(h,x,y) if inverted else h
        factor = _simplify(slope*_D(1/slope,x,2))
        if _zero(factor):
            continue
        separated = separatevars(_simplify(_D(1/factor,y)),symbols=[x,y],dict=True)
        if (separated is None or _zero(separated['coeff'])
                or not separated[x].has(x) or not separated[y].has(y)):
            continue
        k = Dummy('lie_separation_constant')
        gy = k*integrate(separated[y],y)
        second = 1/(k*separated[x]*separated['coeff'])
        fx = _simplify(second/factor-gy)
        check = _simplify(_D(fx,x,2)-second)
        choices = [1] if _zero(check) else solve(check,k)
        for choice in choices:
            if getattr(choice,'has',lambda *args:False)(x,y):
                continue
            candidate = _simplify((fx+gy).subs(k,choice))
            if inverted:
                yield 0,candidate.xreplace({x:y,y:x})
            else:
                yield candidate,0


def _basis_solutions(residual, coefficients, x, y, polynomial):
    from ... import Poly, Add, S, linear_eq_to_matrix
    numerator = _simplify(residual).as_numer_denom()[0].expand()
    if polynomial:
        conditions = list(Poly(numerator,x,y).as_dict().values())
    else:
        groups = {}
        for term in Add.make_args(numerator):
            independent,part = term.as_independent(x,y,as_Add=False)
            groups[part] = groups.get(part,S.Zero)+independent
        conditions = list(groups.values())
    if not conditions:
        from ... import eye
        return [list(column) for column in eye(len(coefficients)).columnspace()]
    matrix,rhs = linear_eq_to_matrix(conditions,coefficients)
    if any(value != 0 for value in rhs):
        return []
    return [list(column) for column in matrix.nullspace()]


def _separated_generators(h,x,y):
    """Search xi=f(x), eta=g(y) directly as separate univariate ansatzes."""
    from ... import Dummy, S
    for degree in range(4):
        left = [Dummy('separated_xi') for _ in range(degree+1)]
        right = [Dummy('separated_eta') for _ in range(degree+1)]
        xi = sum((coefficient*x**power for power,coefficient in enumerate(left)),S.Zero)
        eta = sum((coefficient*y**power for power,coefficient in enumerate(right)),S.Zero)
        residual = _D(eta,y)*h-_D(xi,x)*h-xi*_D(h,x)-eta*_D(h,y)
        for vector in _basis_solutions(residual,left+right,x,y,False):
            replace = dict(zip(left+right,vector))
            yield xi.xreplace(replace),eta.xreplace(replace)


def _general(h,x,y):
    yield from _unique_general(h,x,y)
    # Exact coefficient solving also covers regular degenerate limits of the
    # closed-form differential invariants, without division by a vanishing E.
    yield from _separated_generators(h,x,y)


def _polynomial(h,x,y,mode):
    from ... import Poly, Dummy, S
    if mode == 'linear':
        maximum = 1
        polynomial = False
    else:
        if not h.is_rational_function(x,y):
            return
        maximum = min(6,max(1,Poly(_simplify(h).as_numer_denom()[0],x,y).total_degree()+2))
        polynomial = True
    degrees = [1] if mode == 'linear' else range(maximum+1)
    for degree in degrees:
        monomials = [x**i*y**(total-i) for total in range(degree+1) for i in range(total+1)]
        count = len(monomials)
        coefficients = [Dummy('lie_coefficient') for _ in range(count*(1 if mode=='chi' else 2))]
        if mode == 'chi':
            xi = S.Zero
            eta = sum((c*m for c,m in zip(coefficients,monomials)),S.Zero)
        else:
            xi = sum((c*m for c,m in zip(coefficients[:count],monomials)),S.Zero)
            eta = sum((c*m for c,m in zip(coefficients[count:],monomials)),S.Zero)
        residual = _D(eta,x)+(_D(eta,y)-_D(xi,x))*h-_D(xi,y)*h*h-xi*_D(h,x)-eta*_D(h,y)
        basis = _basis_solutions(residual,coefficients,x,y,polynomial)
        found = False
        for vector in basis:
            replacements = dict(zip(coefficients,vector))
            left,right = _simplify(xi.xreplace(replacements)),_simplify(eta.xreplace(replacements))
            if not _zero(left) or not _zero(right):
                yield left,right
                found = True
        if found:
            return


_HEURISTICS = {'abaco1_simple':_simple,'abaco1_product':_product,
               'abaco2_similar':_similar,'abaco2_unique_unknown':_unique_unknown,
               'abaco2_unique_general':_general,'function_sum':_function_sum}


def _symmetry_infinitesimals(eq, func=None, order=None, hint='default', match=None):
    """Find and verify first-order point symmetry generators.

    Named heuristics search their respective ansatz; ``all`` combines the
    search results. Polynomial ansatzes use exact native linear nullspaces.
    """
    from ... import Function, sympify
    expression,func,x,y,h = _problem(eq,func,order,match)
    if hint not in ('default','all')+lie_heuristics:
        raise ValueError('Heuristic not recognized: '+str(hint))
    xi,eta = Function('xi')(x,func),Function('eta')(x,func)
    selected = lie_heuristics if hint in ('default','all') else (hint,)
    output = []
    for name in selected:
        candidates = _polynomial(h,x,y,name) if name in ('linear','bivariate','chi') else _HEURISTICS[name](h,x,y)
        try:
            for left,right in candidates:
                left,right = _simplify(sympify(left)),_simplify(sympify(right))
                if _illegal(left) or _illegal(right) or _zero(left) and _zero(right):
                    continue
                if not _zero(_residual(h,x,y,left,right)):
                    continue
                from .._native_pde import _at
                result = {xi:_at(left,y,func,{}),eta:_at(right,y,func,{})}
                if result not in output:
                    output.append(result)
                if hint == 'default':
                    return output
                if name not in ('abaco1_simple','abaco1_product'):
                    # These ansatzes return their first verified generator,
                    # matching the reference search's early-exit semantics.
                    break
        except (NotImplementedError,ZeroDivisionError):
            if hint not in ('default','all'):
                raise
    if output:
        return output
    if hint in lie_heuristics:
        raise ValueError('Infinitesimals could not be found using '+hint)
    raise NotImplementedError('No verified Lie point symmetry was found')


def solve_lie(eq,func=None,xi=None,eta=None,ics=None):
    """Integrate the exact differential form supplied by a point symmetry.

    For Q=eta-xi*h, the determining equation implies that
    (dy-h*dx)/Q is closed. Its native potential is a first integral of the ODE.
    """
    from ... import Function, Eq, S, sympify, integrate, solve, Symbol
    from .._native_pde import _at
    expression,func,x,y,h = _problem(eq,func)
    if (xi is None)!=(eta is None):
        raise ValueError('Supply both xi and eta for a custom Lie generator')
    left_key,right_key = Function('xi')(x,func),Function('eta')(x,func)
    if xi is not None:
        left = sympify(xi).xreplace({func:y})
        right = sympify(eta).xreplace({func:y})
        if not _zero(_residual(h,x,y,left,right)):
            raise ValueError('The supplied infinitesimals do not satisfy the determining equation')
        searches = [(left,right)]
    else:
        def generated():
            for heuristic in lie_heuristics:
                try:
                    found = _symmetry_infinitesimals(expression,func,hint=heuristic)
                except (ValueError,NotImplementedError):
                    continue
                for candidate in found:
                    yield candidate[left_key].xreplace({func:y}),candidate[right_key].xreplace({func:y})
        searches = generated()
    occupied = {symbol.name for symbol in expression.free_symbols}
    conditions = {} if ics is None else {sympify(key):sympify(value) for key,value in ics.items()}
    occupied.update(symbol.name for pair in conditions.items() for value in pair for symbol in value.free_symbols)
    index = 1
    while 'C'+str(index) in occupied:
        index += 1
    constant = Symbol('C'+str(index))
    for left,right in searches:
        characteristic = _simplify(right-left*h)
        if _zero(characteristic):
            continue
        mu = _simplify(1/characteristic)
        potential = None
        # Both orders integrate the same closed form. Trying x first also
        # avoids unnecessary radical integration for horizontal generators.
        for variable,first,other,second in ((x,-h*mu,y,mu),(y,mu,x,-h*mu)):
            try:
                primitive = integrate(first,variable)
                remainder = _simplify(second-_D(primitive,other))
                candidate = _simplify(primitive+integrate(remainder,other))
            except (NotImplementedError,ValueError):
                continue
            if _zero(_D(candidate,y)-mu) and _zero(_D(candidate,x)+h*mu):
                potential = candidate
                break
        if potential is None:
            continue
        chosen = constant
        if conditions:
            constraints = []
            for condition,value in conditions.items():
                if condition.func!=func.func or len(condition.args)!=1:
                    raise NotImplementedError('Lie initial conditions require function values')
                evaluated = _at(_at(potential,y,value,{}),x,condition.args[0],{})
                constraints.append(evaluated-constant)
            values = solve(constraints,[constant],dict=True)
            if isinstance(values,list) and len(values)==1:
                values = values[0]
            if not isinstance(values,dict) or constant not in values:
                raise ValueError('The initial conditions are inconsistent or cannot isolate the constant')
            chosen = values[constant]
        try:
            roots = solve(potential-chosen,y)
        except NotImplementedError:
            roots = []
        if roots:
            from ._constants import constantsimp
            solutions = [Eq(func,constantsimp(_simplify(value),{constant}) if chosen==constant else _simplify(value))
                         for value in roots]
            return solutions[0] if len(solutions)==1 else solutions
        return Eq(_at(potential,y,func,{}),chosen)
    raise NotImplementedError('No verified symmetry produced an integrable first integral')


# Reference API of ``sympy/solvers/ode/lie_group.py`` (SymPy 1.14). The
# heuristics, ``infinitesimals``, ``_lie_group_remove`` and the solver run in
# ``symox::ode::lie_group``; this module admits the arguments, hands over the
# public callables and keeps the reference body (``_lie_heuristics``) for a
# native decline.

def _native(name):
    from ... import _core
    return getattr(_core, name, None)


def _expr_args(*values):
    from ... import _core
    return all(isinstance(v, _core.Expr) for v in values)


def _hooks():
    from . import _lie_heuristics as ref
    return {'simplify': ref.simplify, 'expand': ref.expand, 'cancel': ref.cancel,
            'integrate': ref.integrate, 'solve': ref.solve, 'div': ref.div,
            'pdsolve': ref.pdsolve}


def infinitesimals(eq, func=None, order=None, hint='default', match=None):
    from . import _lie_heuristics as ref
    from ...core.relational import Equality
    from ._ode_helpers import _preprocess
    if isinstance(eq, Equality):
        eq = eq.lhs - eq.rhs
    if not func:
        eq, func = _preprocess(eq)
    variables = func.args
    if len(variables) != 1:
        raise ValueError("ODE's have only one independent variable")
    x = variables[0]
    if not order:
        order = ode_order(eq, func)
    if order != 1:
        raise NotImplementedError("Infinitesimals for only "
            "first order ODE's have been implemented")
    native = _native('ode_lie_infinitesimals')
    if native is not None and isinstance(hint, str) and _expr_args(eq, func, x):
        given = (None, None)
        if match:
            given = (match.get('h'), match.get('y'))
        if (given[0] is None) == (given[1] is None) and (not match or _expr_args(*given)):
            value = native(eq, func, x, hint, given[0], given[1], _hooks())
            if value is not None:
                return value
    return ref.infinitesimals(eq, func, order, hint, match)


def _heuristic(name, match, comp):
    from . import _lie_heuristics as ref
    native = _native('ode_lie_heuristic')
    if native is not None and isinstance(match, dict):
        keys = ('h', 'hx', 'hy', 'y', 'hinv', 'func')
        if all(k in match for k in keys) and _expr_args(*(match[k] for k in keys)) \
                and len(match['func'].args) == 1:
            value = native(name, *(match[k] for k in keys), match['func'].args[0],
                           bool(comp), _hooks())
            if value is not None:
                return value[0]
    return getattr(ref, 'lie_heuristic_' + name)(match, comp)


def lie_heuristic_abaco1_simple(match, comp=False):
    return _heuristic('abaco1_simple', match, comp)


def lie_heuristic_abaco1_product(match, comp=False):
    return _heuristic('abaco1_product', match, comp)


def lie_heuristic_bivariate(match, comp=False):
    return _heuristic('bivariate', match, comp)


def lie_heuristic_chi(match, comp=False):
    return _heuristic('chi', match, comp)


def lie_heuristic_function_sum(match, comp=False):
    return _heuristic('function_sum', match, comp)


def lie_heuristic_abaco2_similar(match, comp=False):
    return _heuristic('abaco2_similar', match, comp)


def lie_heuristic_abaco2_unique_unknown(match, comp=False):
    return _heuristic('abaco2_unique_unknown', match, comp)


def lie_heuristic_abaco2_unique_general(match, comp=False):
    return _heuristic('abaco2_unique_general', match, comp)


def lie_heuristic_linear(match, comp=False):
    return _heuristic('linear', match, comp)


def _lie_group_remove(coords):
    from . import _lie_heuristics as ref
    native = _native('ode_lie_group_remove')
    if native is not None and _expr_args(coords):
        return native(coords)
    return ref._lie_group_remove(coords)


def __getattr__(name):
    # checkinfsol is the reference one of ``ode.py``
    if name == 'checkinfsol':
        from .ode import checkinfsol
        return checkinfsol
    raise AttributeError(name)


__all__ = ['infinitesimals','checkinfsol','lie_heuristics']
