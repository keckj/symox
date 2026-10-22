"""D-finite functions represented by native polynomial differential equations.

Closure operations use cyclic vectors in finite dimensional differential
systems. No symbolic reference engine is used to discover annihilators.
"""
import math
import operator
from .recurrence import _OreAlgebra, _OreOperator, _expr, _base, HolonomicSequence, RecurrenceOperators


class DifferentialOperator(_OreOperator):
    def is_singular(self, point):
        return self.coefficients[-1].subs(self.parent.x, point) == 0


class DifferentialOperatorAlgebra(_OreAlgebra):
    differential = True
    operator_class = DifferentialOperator


def DifferentialOperators(base, generator=None):
    parent = DifferentialOperatorAlgebra(base, generator)
    return parent, parent.derivative_operator


def _normal(coefficients, parent):
    from .. import cancel, expand, lcm, gcd
    denominator = _expr(1)
    for c in coefficients:
        denominator = lcm(denominator, cancel(c).as_numer_denom()[1])
    values = [expand(cancel(c*denominator)) for c in coefficients]
    common = _expr(0)
    for value in values:
        common = gcd(common,value)
    if common != 0:
        values = [cancel(value/common) for value in values]
    return DifferentialOperator(values, parent)


def _step(row, matrix, x):
    from .. import cancel
    return [cancel(row[j].diff(x) + sum(row[i]*matrix[i][j] for i in range(len(row))))
            for j in range(len(row))]


def _cyclic(matrix, row, x, parent):
    from .. import Matrix
    if all(c == 0 for c in row):
        return parent.derivative_operator
    rows = [row]
    for _ in range(len(row)):
        rows.append(_step(rows[-1], matrix, x))
        null = Matrix(rows).T.nullspace()
        if null:
            return _normal(list(null[0]), parent)
    raise ValueError('unable to find a differential relation in the supplied coefficient field')


class HolonomicFunction:
    _symox_holonomic_function = True
    def __init__(self, annihilator, x, x0=0, y0=None):
        from .. import sympify
        if not isinstance(annihilator, DifferentialOperator):
            raise TypeError('annihilator must be a DifferentialOperator')
        if annihilator.order == 0:
            raise ValueError('a holonomic function requires a nonzero positive-order differential equation')
        self.annihilator = annihilator
        self.x, self.x0 = sympify(x), sympify(x0)
        if self.x != annihilator.parent.x:
            raise ValueError('independent variable must match the operator base ring')
        if isinstance(y0, dict):
            raise NotImplementedError('Frobenius initial conditions are not supported; supply ordinary derivatives')
        self.y0 = None if y0 is None else [sympify(v) for v in y0]
        self._expr = None

    def __repr__(self):
        return 'HolonomicFunction(%r, %s, %s, %r)' % (self.annihilator, self.x, self.x0, self.y0)

    def __eq__(self, other):
        return isinstance(other, HolonomicFunction) and (self.annihilator, self.x, self.x0, self.y0) == (other.annihilator, other.x, other.x0, other.y0)

    def _have_init_cond(self):
        return bool(self.y0)

    def is_singularics(self):
        return self.annihilator.is_singular(self.x0)

    def _initial_count(self):
        """Include every nonnegative integral indicial exponent in the jet.

        At an apparent singularity, the differential order alone need not
        determine an analytic solution: x*sin(x) has two zero initial values.
        """
        from .. import Poly, Dummy, expand
        r = self.annihilator.order
        if not self.is_singularics():
            return r
        polynomials = [Poly(c.subs(self.x,self.x+self.x0).expand(),self.x).as_dict()
                       for c in self.annihilator.coefficients]
        valuation = min(j-i for i,p in enumerate(polynomials) for (j,),c in p.items() if c != 0)
        k = Dummy('exponent')
        indicial = _expr(0)
        for i,p in enumerate(polynomials):
            for (j,),c in p.items():
                if j-i == valuation:
                    for v in range(i):
                        c *= k-v
                    indicial += c
        roots = Poly(expand(indicial),k).ground_roots()
        return max([r]+[int(root)+1 for root in roots if root.is_Integer and int(root) >= 0])

    def degree(self):
        from .. import Poly
        return max(Poly(c,self.x).degree() for c in self.annihilator.coefficients if c != 0)

    def _system(self):
        from .. import cancel, Integer
        r = self.annihilator.order
        if r == 0:
            raise ValueError('order-zero annihilator has no companion differential system')
        matrix = [[Integer(0) for j in range(r)] for i in range(r)]
        for i in range(r-1):
            matrix[i][i+1] = Integer(1)
        c = self.annihilator.coefficients
        matrix[-1] = [cancel(-a/c[-1]) for a in c[:-1]]
        return matrix

    def _from_system(self, matrix, row, expression=None, x0=None, y0=None):
        result = HolonomicFunction(_cyclic(matrix, row, self.x, self.annihilator.parent), self.x,
                                   self.x0 if x0 is None else x0, y0)
        result._expr = expression
        if expression is not None and y0 is None:
            result.y0 = _conditions(expression, self.x, result.x0, result._initial_count())
        elif expression is not None and y0 is not None:
            if list(y0) != _conditions(expression,self.x,result.x0,len(y0)):
                result._expr = None
        return result

    def _coerce(self, other):
        if not isinstance(other, HolonomicFunction):
            other = expr_to_holonomic(other, self.x, self.x0, domain=self.annihilator.parent.base.domain)
        if self.x != other.x:
            raise ValueError('holonomic functions must use the same independent variable')
        return other

    def unify(self, other):
        other = self._coerce(other)
        domain = self.annihilator.parent.base.domain.unify(other.annihilator.parent.base.domain)
        parent, _ = DifferentialOperators(_base(domain,self.x), str(self.annihilator.parent.gen_symbol))
        result = []
        for value in (self, other):
            item = HolonomicFunction(DifferentialOperator(value.annihilator.coefficients, parent), self.x, value.x0, value.y0)
            item._expr = value._expr
            result.append(item)
        return tuple(result)

    def _binary(self, other, product=False):
        from .. import Integer
        left, right = self.unify(other)
        a, b = left._system(), right._system()
        m, n = len(a), len(b)
        size = m*n if product else m+n
        matrix = [[Integer(0) for j in range(size)] for i in range(size)]
        if product:
            for i in range(m):
                for j in range(n):
                    for k in range(m):
                        matrix[i*n+j][k*n+j] += a[i][k]
                    for k in range(n):
                        matrix[i*n+j][i*n+k] += b[j][k]
        else:
            for i in range(m):
                for j in range(m):
                    matrix[i][j] = a[i][j]
            for i in range(n):
                for j in range(n):
                    matrix[m+i][m+j] = b[i][j]
        row = [Integer(0)]*size
        row[0] = Integer(1)
        if not product:
            row[m] = Integer(1)
        expression = None
        if left._expr is not None and right._expr is not None:
            expression = left._expr*right._expr if product else left._expr+right._expr
        result = left._from_system(matrix, row, expression)
        if expression is None and left.y0 is not None and right.y0 is not None:
            if left.x0 != right.x0:
                raise ValueError('initial conditions must be at the same point; use change_ics first')
            count = result._initial_count()
            lv, rv = left._derivatives(count), right._derivatives(count)
            result.y0 = [sum(math.comb(i,j)*lv[j]*rv[i-j] for j in range(i+1)) if product else lv[i]+rv[i] for i in range(count)]
        return result

    def __add__(self, other):
        return self._binary(other)

    __radd__ = __add__

    def __mul__(self, other):
        return self._binary(other, True)

    __rmul__ = __mul__

    def __neg__(self):
        result = HolonomicFunction(self.annihilator, self.x, self.x0, None if self.y0 is None else [-v for v in self.y0])
        result._expr = None if self._expr is None else -self._expr
        return result

    def __sub__(self, other):
        return self + -self._coerce(other)

    def __rsub__(self, other):
        return self._coerce(other) + -self

    def __truediv__(self, other):
        from .. import sympify
        if isinstance(other, HolonomicFunction):
            raise NotImplementedError('D-finite functions are not closed under division')
        other = sympify(other)
        if other == 0:
            raise ZeroDivisionError('division by zero')
        return self * expr_to_holonomic(1/other, self.x, self.x0)

    def __pow__(self, power):
        power = operator.index(power)
        if power < 0:
            raise NotImplementedError('negative powers of arbitrary D-finite functions need not be D-finite')
        result = expr_to_holonomic(1, self.x, self.x0)
        for _ in range(power):
            result = result*self
        return result

    def diff(self, *args, **kwargs):
        from .. import Integer
        if kwargs:
            raise TypeError('unexpected differentiation options')
        if not args:
            count = 1
        elif args[0] == self.x:
            count = operator.index(args[1]) if len(args) == 2 else 1
        elif len(args) == 1:
            count = operator.index(args[0])
        else:
            raise ValueError('differentiate with respect to the independent variable')
        if count < 0:
            raise ValueError('derivative order must be nonnegative')
        if count == 0:
            return self
        matrix = self._system()
        row = [Integer(1)] + [Integer(0)]*(len(matrix)-1)
        for _ in range(count):
            row = _step(row, matrix, self.x)
        result = self._from_system(matrix, row, None if self._expr is None else self._expr.diff(self.x, count))
        if self._expr is None and self.y0 is not None:
            result.y0 = self._derivatives(count+result._initial_count())[count:]
        return result

    def integrate(self, limits):
        from .. import integrate
        definite = isinstance(limits, (tuple, list)) or getattr(limits, 'is_Tuple', False)
        limits = tuple(limits.args) if getattr(limits, 'is_Tuple', False) else limits
        if definite:
            if len(limits) != 3 or limits[0] != self.x:
                raise ValueError('expected (independent variable, lower bound, upper bound)')
            if limits[2] != self.x:
                return self.integrate(self.x).to_expr().subs(self.x, limits[2]) - self.integrate(self.x).to_expr().subs(self.x, limits[1])
            base = _expr(limits[1])
        else:
            if limits != self.x:
                raise ValueError('integration variable must be the independent variable')
            base = self.x0
        parent = self.annihilator.parent
        result = HolonomicFunction(self.annihilator*parent.derivative_operator, self.x, base)
        if self._expr is not None:
            primitive = integrate(self._expr, self.x)
            result._expr = primitive - primitive.subs(self.x, base)
            result.y0 = _conditions(result._expr, self.x, base, result._initial_count())
        elif self.y0 is not None and base == self.x0:
            result.y0 = [_expr(0)] + self.y0
        return result

    def composition(self, expression, *args, **kwargs):
        from .. import Integer, cancel
        expression = _expr(expression)
        numerator, denominator = expression.as_numer_denom()
        if not numerator.is_polynomial(self.x) or not denominator.is_polynomial(self.x):
            raise NotImplementedError('composition currently supports rational substitutions')
        x0 = args[0] if args else kwargs.pop('x0', self.x0)
        y0 = args[1] if len(args) > 1 else kwargs.pop('y0', None)
        if len(args) > 2 or kwargs:
            raise TypeError('expected composition(expression, x0=None, y0=None)')
        matrix = [[cancel(c.subs(self.x, expression)*expression.diff(self.x)) for c in row] for row in self._system()]
        result = self._from_system(matrix, [Integer(1)]+[Integer(0)]*(len(matrix)-1),
                                   None if self._expr is None else self._expr.subs(self.x, expression), x0, y0)
        if result._expr is None and y0 is None and self.y0 is not None and expression.subs(self.x, x0) == self.x0:
            # A sufficiently long Taylor jet determines the chain rule exactly.
            count = result._initial_count()
            jet = self.series(n=count+1, order=False)
            result.y0 = _conditions(jet.subs(self.x, expression), self.x, _expr(x0), count)
        return result

    def shift_x(self, amount):
        return self.composition(self.x+amount, self.x0-_expr(amount))

    def change_x(self, x):
        parent, _ = DifferentialOperators(_base(self.annihilator.parent.base.domain,x), str(self.annihilator.parent.gen_symbol))
        result = HolonomicFunction(DifferentialOperator([c.subs(self.x, x) for c in self.annihilator.coefficients], parent), x, self.x0, self.y0)
        result._expr = None if self._expr is None else self._expr.subs(self.x, x)
        return result

    def _coefficients(self, count):
        from .. import Poly, cancel, Integer, simplify
        count = operator.index(count)
        if count < 0:
            raise ValueError('series length must be nonnegative')
        if self.y0 is None:
            raise ValueError('series requires initial derivative values')
        r = self.annihilator.order
        if len(self.y0) < r:
            raise ValueError('insufficient initial derivatives')
        values = [v/math.factorial(i) for i,v in enumerate(self.y0[:count])]
        polynomials = [Poly(c.subs(self.x, self.x+self.x0).expand(), self.x).as_dict() for c in self.annihilator.coefficients]
        relations = {}
        degree = max((j for poly in polynomials for (j,) in poly), default=0)
        for m in range(count+degree+1):
            relation = {}
            for i, poly in enumerate(polynomials):
                for (j,), coefficient in poly.items():
                    if j <= m:
                        index = m-j+i
                        factor = coefficient*math.factorial(index)/math.factorial(m-j)
                        relation[index] = relation.get(index,Integer(0))+factor
            relation = {k:cancel(v) for k,v in relation.items() if cancel(v) != 0}
            if relation:
                relations.setdefault(max(relation),[]).append(relation)
        for target in range(len(values), count):
            for known,candidates in relations.items():
                if known < len(values):
                    for relation in candidates:
                        residual = simplify(sum((factor*values[index] for index,factor in relation.items()),Integer(0)))
                        if residual.is_number and residual != 0:
                            raise ValueError('initial derivatives are inconsistent with the differential equation')
            candidates = relations.get(target)
            if not candidates:
                if self._expr is not None:
                    extra = _conditions(self._expr,self.x,self.x0,target+1)
                    if extra is not None:
                        values.append(extra[target]/math.factorial(target))
                        continue
                raise NotImplementedError('Taylor coefficient %s is undetermined; supply more initial derivatives (Frobenius branches are not implemented)' % target)
            relation = candidates[0]
            total = sum((factor*values[index] for index,factor in relation.items() if index != target),Integer(0))
            values.append(cancel(-total/relation[target]))
        for target,candidates in relations.items():
            if target < len(values):
                for relation in candidates:
                    residual = simplify(sum((factor*values[index] for index,factor in relation.items()),Integer(0)))
                    if residual.is_number and residual != 0:
                        raise ValueError('initial derivatives are inconsistent with the differential equation')
        return values

    def _derivatives(self, count):
        return [c*math.factorial(i) for i,c in enumerate(self._coefficients(count))]

    def series(self, n=6, coefficient=False, order=True, _recur=None):
        from .. import Order
        if _recur is not None:
            raise NotImplementedError('custom Frobenius recurrence input is not implemented')
        coefficients = self._coefficients(n)
        if coefficient:
            return coefficients
        result = sum((a*(self.x-self.x0)**i for i,a in enumerate(coefficients)), _expr(0))
        return result + Order((self.x-self.x0)**n, (self.x,self.x0)) if order else result

    def to_sequence(self, lb=True):
        from .. import Symbol, QQ, Poly, expand
        n = Symbol('n', integer=True)
        polys = [Poly(c.subs(self.x, self.x+self.x0).expand(), self.x).as_dict() for c in self.annihilator.coefficients]
        shifts = [i-j for i,p in enumerate(polys) for (j,),c in p.items() if c != 0]
        lower, upper = min(shifts), max(shifts)
        coefficients = [_expr(0)]*(upper-lower+1)
        for i,p in enumerate(polys):
            for (j,),c in p.items():
                factor = c
                for k in range(1,i+1):
                    factor *= n-lower-j+k
                coefficients[i-j-lower] += expand(factor)
        parent, _ = RecurrenceOperators(_base(self.annihilator.parent.base.domain,n), 'Sn')
        if self.is_singularics():
            count, smallest = max(upper-lower, upper, self._initial_count()), max(0, lower)
        else:
            # Reference: the recurrence holds for n >= lower + degree and
            # beyond the integer roots of its leading coefficient; it needs
            # order + smallest_n initial terms.
            from .. import roots
            degree = max(j for p in polys for (j,),c in p.items() if c != 0)
            smallest = lower + degree
            integer_roots = list(roots(coefficients[-1], n, filter='Z'))
            if integer_roots:
                smallest = max(max(integer_roots) + 1, smallest)
            count = max(upper - lower + smallest, len(self.y0 or ()))
        sequence = HolonomicSequence(parent.operator_class(coefficients, parent), self._coefficients(count))
        return [(sequence, smallest)] if lb else [sequence]

    def to_expr(self):
        from .. import Function, Eq, dsolve, solve
        if self._expr is not None:
            return self._expr
        f = Function('y')(self.x)
        equation = Eq(sum(c*f.diff(self.x,i) for i,c in enumerate(self.annihilator.coefficients)), 0)
        solution = dsolve(equation, f)
        expression = solution.rhs
        if self.y0 is not None:
            parameters = set().union(*(c.free_symbols for c in self.annihilator.coefficients))
            constants = sorted(expression.free_symbols-parameters-{self.x}, key=str)
            conditions = [expression.diff(self.x,i).subs(self.x,self.x0)-v for i,v in enumerate(self.y0)]
            values = solve(conditions, constants, dict=True)
            if not values:
                raise ValueError('unable to resolve the initial conditions')
            expression = expression.subs(values[0])
        self._expr = expression
        return expression

    def to_hyper(self, as_list=False):
        from .. import hyper, Rational
        def convert(value):
            name = str(value.func)
            if name in ('hyper','meijerg') or not value.has(self.x):
                return value
            if name == 'exp':
                return hyper([],[],value.args[0])
            if name in ('sin','cos','sinh','cosh'):
                arg = value.args[0]
                odd = name in ('sin','sinh')
                return (arg if odd else 1)*hyper([],[Rational(3,2) if odd else Rational(1,2)],(1 if name in ('sinh','cosh') else -1)*arg**2/4)
            if value.is_Add or value.is_Mul:
                items = [convert(a) for a in value.args]
                from .. import Add,Mul
                return Add(*items) if value.is_Add else Mul(*items)
            if value.is_polynomial(self.x):
                return value
            raise NotImplementedError('hypergeometric reconstruction currently supports exponentials and trigonometric/hyperbolic functions with polynomial prefactors')
        expression = convert(self.to_expr())
        if as_list:
            return [expression]
        return expression

    def to_meijerg(self):
        if self._expr is not None and str(self._expr.func) == 'meijerg':
            return self._expr
        raise NotImplementedError('general Meijer G reconstruction from an annihilator is not implemented')

    def change_ics(self, b, lenics=None):
        b = _expr(b)
        count = self.annihilator.order if lenics is None else operator.index(lenics)
        result = HolonomicFunction(self.annihilator, self.x, b)
        result._expr = self._expr
        if self._expr is not None:
            result.y0 = _conditions(self._expr, self.x, b, count)
        else:
            if count != self.annihilator.order:
                raise NotImplementedError('numeric initial-condition transfer returns exactly the equation order')
            result.y0 = self.evalf(b, derivatives=True)
        return result

    def evalf(self, points, method='RK4', h=0.05, derivatives=False):
        from .. import Float
        if method not in ('RK4', 'Euler'):
            raise ValueError('method must be RK4 or Euler')
        if h <= 0:
            raise ValueError('step size must be positive')
        if self.y0 is None or len(self.y0) < self.annihilator.order:
            raise ValueError('numeric continuation requires complete initial derivatives')
        if self.is_singularics():
            raise ValueError('cannot start numeric continuation at a singular point')
        matrix = self._system()
        def numeric(value):
            real,imag = value.as_real_imag()
            return complex(float(real.evalf()),float(imag.evalf()))
        def rhs(t, y):
            native = Float(t.real)+Float(t.imag)*_expr('I')
            if self.annihilator.is_singular(native):
                raise ValueError('numeric path crosses a singularity')
            return [sum(numeric(c.subs(self.x,native))*v for c,v in zip(row,y)) for row in matrix]
        def advance(start, target, state, steps):
            delta = (target-start)/steps
            for k in range(steps):
                t = start+k*delta
                k1 = rhs(t,state)
                if method == 'Euler':
                    state = [v+delta*d for v,d in zip(state,k1)]
                else:
                    k2 = rhs(t+delta/2,[v+delta*d/2 for v,d in zip(state,k1)])
                    k3 = rhs(t+delta/2,[v+delta*d/2 for v,d in zip(state,k2)])
                    k4 = rhs(t+delta,[v+delta*d for v,d in zip(state,k3)])
                    state = [v+delta*(a+2*b+2*c+d)/6 for v,a,b,c,d in zip(state,k1,k2,k3,k4)]
            return state
        def output(state):
            result = [Float(v.real)+Float(v.imag)*_expr('I') for v in state]
            return result if derivatives else result[0]
        start = numeric(self.x0)
        state = [numeric(v) for v in self.y0[:len(matrix)]]
        if isinstance(points,(list,tuple)):
            # Reference: the given points are the mesh, one Euler/RK4 step
            # from each point to the next (starting at x0); h is unused.
            results = []
            for p in points:
                target = numeric(_expr(p))
                state = advance(start, target, state, 1)
                start = target
                results.append(output(state))
            return results
        # A single endpoint is reached exactly with steps of length <= h.
        target = numeric(_expr(points))
        steps = max(1, math.ceil(abs(target-start)/h))
        if steps > 100000:
            raise ValueError('numeric path exceeds 100000 steps; increase h')
        return output(advance(start, target, state, steps))


def _conditions(expression, x, x0, count):
    from .. import oo, zoo, nan
    if str(expression.func) == 'meijerg':
        from .. import hyperexpand
        expression = hyperexpand(expression)
    def special_zero(value):
        replacements = {}
        def visit(node):
            if str(node.func) == 'hyper' and node.args[-1] == 0:
                if any(b.is_Integer and int(b) <= 0 for b in node.bq):
                    return
                replacements[node] = _expr(1)
            for arg in node.args:
                visit(arg)
        visit(value)
        return value.xreplace(replacements)
    values = [special_zero(expression.diff(x,i).subs(x,x0)) for i in range(count)]
    if any(v.has(oo,zoo,nan) for v in values):
        return None
    return values


def expr_to_holonomic(func, x=None, x0=0, y0=None, lenics=None, domain=None, initcond=True):
    from .. import QQ, Symbol, cancel
    func = _expr(func)
    if x is None:
        symbols = sorted(func.free_symbols, key=str)
        if len(symbols) > 1:
            raise ValueError('specify the independent variable for multivariate expressions')
        x = symbols[0] if symbols else Symbol('x')
    if domain is None:
        parameters = sorted(func.free_symbols-{x},key=str)
        domain = QQ.frac_field(*parameters) if parameters else QQ
    parent, D = DifferentialOperators(_base(domain,x), 'Dx')
    name = str(func.func)
    # Prefer a first-order logarithmic derivative for rational/algebraic powers
    # and exponentials whose logarithmic derivative is rational.
    ratio = cancel(func.diff(x)/func) if func != 0 else _expr(0)
    num, den = ratio.as_numer_denom()
    if num.is_polynomial(x) and den.is_polynomial(x):
        result = HolonomicFunction(_normal([-ratio, _expr(1)], parent), x, x0)
    elif func.is_Add or func.is_Mul:
        items = [expr_to_holonomic(a,x,x0,domain=domain,initcond=False) for a in func.args]
        result = items[0]
        for item in items[1:]:
            result = result+item if func.is_Add else result*item
    elif name in ('sin','cos','sinh','cosh','log','erf','erfc','besselj','bessely','besseli','besselk','airyai','airybi'):
        arg = func.args[-1]
        if name in ('sin','cos'):
            base = D**2+1
        elif name in ('sinh','cosh'):
            base = D**2-1
        elif name == 'log':
            base = x*D**2+D
        elif name in ('erf','erfc'):
            base = D**2+2*x*D
        elif name in ('airyai','airybi'):
            base = D**2-x
        else:
            nu = func.args[0]
            base = x**2*D**2+x*D+((x**2 if name in ('besselj','bessely') else -x**2)-nu**2)
        result = HolonomicFunction(base,x,x0)
        if arg != x:
            result = result.composition(arg)
    elif name == 'hyper':
        return _special(func, x, x0, y0, lenics, domain, initcond, False)
    elif name == 'meijerg':
        return _special(func, x, x0, y0, lenics, domain, initcond, True)
    elif func.is_Pow and func.args[1].is_Integer and int(func.args[1]) >= 0:
        result = expr_to_holonomic(func.args[0],x,x0,domain=domain,initcond=False)**int(func.args[1])
    else:
        raise NotImplementedError('no D-finite conversion implemented for %s' % func)
    result._expr = func
    count = result._initial_count() if lenics is None else operator.index(lenics)
    result.y0 = [_expr(v) for v in y0] if y0 is not None else _conditions(func,x,_expr(x0),count) if initcond else None
    if y0 is not None and result.y0 != _conditions(func,x,_expr(x0),len(result.y0)):
        result._expr = None
    return result


def _special(func, x, x0, y0, lenics, domain, initcond, meijer):
    parent, D = DifferentialOperators(_base(domain,x), 'Dx')
    theta = x*D
    if meijer:
        a, b = tuple(func.ap), tuple(func.bq)
        left = parent(x*(-1)**(len(func.bm)+len(func.an)-len(a)))
        right = parent(1)
        for value in a:
            left = left*(theta+1-value)
        for value in b:
            right = right*(theta-value)
    else:
        a, b = tuple(func.ap), tuple(func.bq)
        left, right = parent(x), theta
        for value in a:
            left = left*(theta+value)
        for value in b:
            right = right*(theta+value-1)
    result = HolonomicFunction(_normal((left-right).coefficients,parent),x,x0)
    argument = func.args[-1]
    if argument != x:
        result = result.composition(argument)
    result._expr = func
    count = result._initial_count() if lenics is None else operator.index(lenics)
    result.y0 = [_expr(v) for v in y0] if y0 is not None else _conditions(func,x,_expr(x0),count) if initcond else None
    if y0 is not None and result.y0 != _conditions(func,x,_expr(x0),len(result.y0)):
        result._expr = None
    return result


def from_hyper(func, x0=0, evalf=False):
    if str(func.func) != 'hyper':
        raise TypeError('expected a hypergeometric function')
    symbols = tuple(func.args[-1].free_symbols)
    if len(symbols) != 1:
        raise ValueError('hypergeometric argument must have exactly one independent variable')
    result = expr_to_holonomic(func,x=symbols[0],x0=x0)
    if evalf and result.y0 is not None:
        result.y0 = [v.evalf() for v in result.y0]
    return result


def from_meijerg(func, x0=0, evalf=False, initcond=True, domain=None):
    if str(func.func) != 'meijerg':
        raise TypeError('expected a Meijer G-function')
    symbols = tuple(func.args[-1].free_symbols)
    if len(symbols) != 1:
        raise ValueError('Meijer G argument must have exactly one independent variable')
    result = expr_to_holonomic(func,x=symbols[0],x0=x0,domain=domain,initcond=initcond)
    if evalf and result.y0 is not None:
        result.y0 = [v.evalf() for v in result.y0]
    return result
