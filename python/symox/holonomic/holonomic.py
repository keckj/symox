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
            # Reference singular (Frobenius) format {p: [c0, c1, ...]}:
            # the function is sum_p x**p*(c0 + c1*x + ...) about x0.
            self.y0 = {sympify(k): [sympify(v) for v in vals] for k, vals in y0.items()}
        else:
            self.y0 = None if y0 is None else [sympify(v) for v in y0]
        self._expr = None

    def __repr__(self):
        if self._have_init_cond():
            return 'HolonomicFunction(%r, %s, %s, %r)' % (self.annihilator, self.x, self.x0, self.y0)
        return 'HolonomicFunction(%r, %s)' % (self.annihilator, self.x)

    def __eq__(self, other):
        return isinstance(other, HolonomicFunction) and (self.annihilator, self.x, self.x0, self.y0) == (other.annihilator, other.x, other.x0, other.y0)

    def _have_init_cond(self):
        return bool(self.y0)

    def is_singularics(self):
        # Reference: True for initial conditions in the singular dictionary
        # format, False for the ordinary list format, None without any.
        if isinstance(self.y0, dict):
            return True
        if isinstance(self.y0, list):
            return False
        return None

    def _singular_point(self):
        return self.annihilator.is_singular(self.x0)

    def _singularics_to_ord(self):
        from .. import factorial, Integer
        a = list(self.y0)[0]
        b = self.y0[a]
        if len(self.y0) == 1 and a == int(a) and a > 0:
            a = int(a)
            y0 = [Integer(0)]*a + [j*factorial(a + i) for i, j in enumerate(b)]
            return HolonomicFunction(self.annihilator, self.x, self.x0, y0)

    def _initial_count(self):
        """Include every nonnegative integral indicial exponent in the jet.

        At an apparent singularity, the differential order alone need not
        determine an analytic solution: x*sin(x) has two zero initial values.
        """
        from .. import Poly, Dummy, expand
        r = self.annihilator.order
        if not self._singular_point():
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
        if isinstance(left.y0, dict) or isinstance(right.y0, dict):
            result.y0 = None
            if left.y0 and right.y0 and left.x0 == right.x0:
                result.y0 = _singular_combination(left.y0, right.y0, product)
            return result
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
        if isinstance(self.y0, dict):
            return HolonomicFunction(self.annihilator, self.x, self.x0, {k: [-v for v in c] for k, c in self.y0.items()})
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
        if isinstance(self.y0, dict):
            # Reference: singular initial conditions are not differentiated.
            result.y0 = None
        elif self._expr is None and self.y0 is not None:
            result.y0 = self._derivatives(count+result._initial_count())[count:]
        return result

    def integrate(self, limits):
        from .. import integrate
        if isinstance(self.y0, dict):
            # Reference: singular conditions integrate termwise.
            ordinary = self._singularics_to_ord()
            if ordinary:
                return ordinary.integrate(limits)
            y0 = {}
            for i, c in self.y0.items():
                c2 = []
                for j, cj in enumerate(c):
                    if cj == 0:
                        c2.append(_expr(0))
                    elif i + j + 1 == 0:
                        raise NotImplementedError("logarithmic terms in the series are not supported")
                    else:
                        c2.append(cj/_expr(i + j + 1))
                y0[i + 1] = c2
            if hasattr(limits, "__iter__"):
                raise NotImplementedError("Definite integration for singular initial conditions")
            return HolonomicFunction(self.annihilator*self.annihilator.parent.derivative_operator, self.x, self.x0, y0)
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
        if result._expr is None and y0 is None and isinstance(self.y0, list) and expression.subs(self.x, x0) == self.x0:
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
        if _recur is None and isinstance(self.y0, list):
            try:
                coefficients = self._coefficients(n)
            except NotImplementedError:
                # Undetermined Taylor coefficients: the reference keeps them
                # as C_j symbols through the recurrence of to_sequence.
                pass
            else:
                if coefficient:
                    return coefficients
                result = sum((a*(self.x-self.x0)**i for i,a in enumerate(coefficients)), _expr(0))
                return result + Order((self.x-self.x0)**n, (self.x,self.x0)) if order else result
        if _recur is None and self.y0 is None and not self._singular_point():
            raise ValueError('series requires initial derivative values')
        return self._reference_series(n, coefficient, order, _recur)

    def _reference_series(self, n=6, coefficient=False, order=True, _recur=None):
        # Reference HolonomicFunction.series, from the recurrence of
        # to_sequence (Frobenius exponents included).
        from .. import Order, cancel, Integer
        recurrence = self.to_sequence() if _recur is None else _recur
        if isinstance(recurrence, tuple) and len(recurrence) == 2:
            recurrence, constantpower = recurrence[0], 0
        elif isinstance(recurrence, tuple) and len(recurrence) == 3:
            constantpower, recurrence = recurrence[1], recurrence[0]
        elif len(recurrence) == 1 and len(recurrence[0]) == 2:
            recurrence, constantpower = recurrence[0][0], 0
        elif len(recurrence) == 1 and len(recurrence[0]) == 3:
            constantpower, recurrence = recurrence[0][1], recurrence[0][0]
        else:
            return [self._reference_series(_recur=i) for i in recurrence]
        n = n - int(constantpower)
        l = len(recurrence.u0) - 1
        k = recurrence.recurrence.order
        x, x0 = self.x, self.x0
        seq = recurrence.recurrence.coefficients
        sub = [cancel(-seq[i]/seq[k]) for i in range(k)]
        sol = list(recurrence.u0)
        if l + 1 < n:
            for i in range(l + 1 - k, n - k):
                sol.append(sum((sub[j].subs(recurrence.n, i)*sol[i + j]
                                for j in range(k) if i + j >= 0), Integer(0)))
        if coefficient:
            return sol
        ser = sum((x**(i + constantpower)*j for i, j in enumerate(sol)), Integer(0))
        if order:
            ser += Order(x**(n + int(constantpower)), x)
        if x0 != 0:
            return ser.subs(x, x - x0)
        return ser

    def _shift_to_zero(self):
        # Reference shift_x(x0): the annihilator moved to 0, y0 kept.
        coefficients = [c.subs(self.x, self.x + self.x0).expand() for c in self.annihilator.coefficients]
        return HolonomicFunction(DifferentialOperator(coefficients, self.annihilator.parent), self.x, 0, self.y0)

    def to_sequence(self, lb=True):
        if isinstance(self.y0, list) and not self._singular_point():
            try:
                return self._ordinary_sequence(lb)
            except NotImplementedError:
                pass
        if self.x0 != 0:
            return self._shift_to_zero().to_sequence(lb)
        if self.annihilator.is_singular(self.x0):
            return self._frobenius(lb=lb)
        if self.y0 is None:
            raise ValueError('series requires initial derivative values')
        return self._reference_ordinary_sequence(lb)

    def _ordinary_sequence(self, lb=True):
        from .. import Symbol, Poly, expand, roots
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
        # Reference: the recurrence holds for n >= lower + degree and
        # beyond the integer roots of its leading coefficient; it needs
        # order + smallest_n initial terms.
        degree = max(j for p in polys for (j,),c in p.items() if c != 0)
        smallest = lower + degree
        integer_roots = list(roots(coefficients[-1], n, filter='Z'))
        if integer_roots:
            smallest = max(max(integer_roots) + 1, smallest)
        count = max(upper - lower + smallest, len(self.y0 or ()))
        sequence = HolonomicSequence(parent.operator_class(coefficients, parent), self._coefficients(count))
        return [(sequence, smallest)] if lb else [sequence]

    def _recurrence_terms(self, p, singular):
        # Reference: x**k*Dx**i contributes rf(n - k + 1 + p, i) to the
        # coefficient of index (i - k) (and degree k, or k - i at a
        # regular singular point).
        from .. import Poly, rf, expand
        n = self._n
        dict1 = {}
        for i, c in enumerate(self.annihilator.coefficients):
            for (k,), coeff in Poly(c, self.x).as_dict().items():
                if coeff == 0:
                    continue
                key = (i - k, k - i) if singular else (i - k, k)
                dict1[key] = dict1.get(key, _expr(0)) + coeff*expand(rf(n - k + 1 + p, i))
        return dict1

    @property
    def _n(self):
        from .. import Symbol
        return Symbol('n', integer=True)

    def _recurrence_operator(self, dict1):
        from .. import roots, expand, Integer
        n = self._n
        keylist = [i[0] for i in dict1]
        lower, upper = min(keylist), max(keylist)
        sol = []
        for j in range(lower, upper + 1):
            if j in keylist:
                sol.append(expand(sum((v.subs(n, n - lower) for key, v in dict1.items() if key[0] == j), Integer(0))))
            else:
                sol.append(Integer(0))
        parent, _ = RecurrenceOperators(_base(self.annihilator.parent.base.domain, n), 'Sn')
        sol = parent.operator_class(sol, parent)
        all_roots = list(roots(sol.coefficients[-1], n, filter='Z'))
        return sol, lower, all_roots

    def _solve_missing(self, dict1, u0, order, lower_index, upper_index, letter):
        # Reference: equate the coefficients of the low powers of x to zero
        # to find the missing initial terms; the rest stay symbols.
        from .. import Symbol, solve, Integer
        n = self._n
        dummys, eqs, unknowns = {}, [], []
        for i in range(lower_index, upper_index):
            eq = Integer(0)
            for j in dict1:
                if i + j[0] < 0:
                    dummys[i + j[0]] = Integer(0)
                elif i + j[0] < len(u0):
                    dummys[i + j[0]] = u0[i + j[0]]
                elif not i + j[0] in dummys:
                    dummys[i + j[0]] = Symbol(letter + '_%s' % (i + j[0]))
                    unknowns.append(dummys[i + j[0]])
                if j[1] <= i:
                    eq += dict1[j].subs(n, i)*dummys[i + j[0]]
            eqs.append(eq)
        soleqs = solve(eqs, *unknowns)
        if isinstance(soleqs, dict):
            for i in range(len(u0), order):
                if i not in dummys:
                    dummys[i] = Symbol(letter + '_%s' % i)
                u0.append(soleqs[dummys[i]] if dummys[i] in soleqs else dummys[i])
            return True
        for i in range(len(u0), order):
            if i not in dummys:
                dummys[i] = Symbol(letter + '_%s' % i)
            found = False
            for j in soleqs:
                if dummys[i] in j:
                    u0.append(j[dummys[i]])
                    found = True
            if not found:
                u0.append(dummys[i])
        return False

    def _reference_ordinary_sequence(self, lb=True):
        from .. import factorial
        dict1 = self._recurrence_terms(0, False)
        sol, lower, all_roots = self._recurrence_operator(dict1)
        smallest_n = lower + self.degree()
        if all_roots:
            smallest_n = max(max(all_roots) + 1, smallest_n)
        order = sol.order + smallest_n
        y0 = _extend_y0(self, order)
        u0 = [j/factorial(i) for i, j in enumerate(y0)]
        if len(u0) < order:
            self._solve_missing(dict1, u0, order, 0, self.degree(), 'C')
        sequence = HolonomicSequence(sol, u0)
        return [(sequence, smallest_n)] if lb else [sequence]

    def _indicial(self):
        # Reference HolonomicFunction._indicial.
        from .. import Poly, roots, Integer
        if self.x0 != 0:
            return self._shift_to_zero()._indicial()
        x = self.x
        coefficients = [Poly(c, x) for c in self.annihilator.coefficients]
        degree = max(max(q.degree(), 0) for q in coefficients)
        inf = 10*(max(1, degree) + max(1, self.annihilator.order))
        def pole_degree(q):
            return inf if q.is_zero else roots(q.as_expr(), x, filter='Z').get(0, 0)
        b = min(pole_degree(q) - j for j, q in enumerate(coefficients))
        s, y = Integer(0), Integer(1)
        for i, q in enumerate(coefficients):
            if 0 <= i + b <= max(q.degree(), 0):
                s = s + q.nth(i + b)*y
            y *= x - i
        return roots(s.expand(), x)

    def _frobenius(self, lb=True):
        # Reference HolonomicFunction._frobenius: a recurrence for each
        # selected root p of the indicial equation, for the series
        # x**p*sum(u(n)*x**n).
        from .. import ordered, factorial
        indicialroots = self._indicial()
        reals, compl = [], []
        for i in ordered(indicialroots.keys()):
            if i.is_real:
                reals.extend([i]*indicialroots[i])
            else:
                a, b = i.as_real_imag()
                compl.extend([(i, a, b)]*indicialroots[i])
        compl.sort(key=lambda v: v[1])
        compl.sort(key=lambda v: v[2])
        reals.sort()
        grp = []
        for i in reals:
            for j in grp:
                if (j[0] - i).is_integer:
                    j.append(i)
                    break
            else:
                grp.append([i])
        independent = all(len(i) == 1 for i in grp)
        allpos = all(i >= 0 for i in reals)
        allint = all(i.is_integer for i in reals)
        if self.is_singularics() == True:
            rootstoconsider = []
            for i in ordered(self.y0.keys()):
                for j in ordered(indicialroots.keys()):
                    if (j - i) == 0:
                        rootstoconsider.append(i)
        elif allpos and allint:
            rootstoconsider = [min(reals)]
        elif independent:
            rootstoconsider = [i[0] for i in grp] + [j[0] for j in compl]
        elif not allint:
            rootstoconsider = [i for i in reals if not int(i) == i]
        elif not allpos:
            if not self._have_init_cond() or _expr(self.y0[0]).is_finite == False:
                rootstoconsider = [min(reals)]
            else:
                rootstoconsider = [min([i for i in reals if i >= 0])]
        finalsol = []
        char = ord('C')
        for p in rootstoconsider:
            dict1 = self._recurrence_terms(p, True)
            sol, lower, all_roots = self._recurrence_operator(dict1)
            degree = max(i[1] for i in dict1)
            degree2 = min(i[1] for i in dict1)
            smallest_n = lower + degree
            if all_roots:
                smallest_n = max(max(all_roots) + 1, smallest_n)
            order = sol.order + smallest_n
            u0 = []
            if self.is_singularics() == True:
                u0 = list(self.y0[p])
            elif self.is_singularics() == False and p >= 0 and int(p) == p and len(rootstoconsider) == 1:
                y0 = _extend_y0(self, order + int(p))
                if len(y0) > int(p):
                    u0 = [y0[i]/factorial(i) for i in range(int(p), len(y0))]
            if len(u0) < order:
                if self._solve_missing(dict1, u0, order, degree2, degree, chr(char)):
                    finalsol.append((HolonomicSequence(sol, u0), p, smallest_n) if lb else (HolonomicSequence(sol, u0), p))
                    continue
            finalsol.append((HolonomicSequence(sol, u0), p, smallest_n) if lb else (HolonomicSequence(sol, u0), p))
            char += 1
        return finalsol

    def to_expr(self):
        from .. import Function, Eq, dsolve, solve, hyperexpand
        if self._expr is not None:
            return self._expr
        # Reference: hyperexpand(self.to_hyper()).simplify().
        try:
            expression = hyperexpand(self._reference_hyper()).simplify()
        except (NotHyperSeriesError, NotImplementedError, ValueError):
            if isinstance(self.y0, dict):
                raise
            expression = None
        if expression is not None:
            self._expr = expression
            return expression
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

    def _reference_hyper(self, as_list=False, _recur=None):
        # Reference HolonomicFunction.to_hyper.
        from .. import hyper, roots, Poly, Symbol, nsimplify, ordered, Integer
        recurrence = self.to_sequence() if _recur is None else _recur
        if isinstance(recurrence, tuple) and len(recurrence) == 2:
            smallest_n, recurrence, constantpower = recurrence[1], recurrence[0], 0
        elif isinstance(recurrence, tuple) and len(recurrence) == 3:
            smallest_n, constantpower, recurrence = recurrence[2], recurrence[1], recurrence[0]
        elif len(recurrence) == 1 and len(recurrence[0]) == 2:
            smallest_n, recurrence, constantpower = recurrence[0][1], recurrence[0][0], 0
        elif len(recurrence) == 1 and len(recurrence[0]) == 3:
            smallest_n, constantpower, recurrence = recurrence[0][2], recurrence[0][1], recurrence[0][0]
        else:
            sol = self._reference_hyper(as_list=as_list, _recur=recurrence[0])
            for i in recurrence[1:]:
                sol += self._reference_hyper(as_list=as_list, _recur=i)
            return sol
        u0 = list(recurrence.u0)
        r = recurrence.recurrence
        coefficients = r.coefficients
        nn = recurrence.n
        x, x0 = self.x, self.x0
        m = r.order
        if m == 0:
            nonzeroterms = roots(coefficients[0], nn, filter='R')
            sol = Integer(0)
            for j, i in enumerate(nonzeroterms):
                if i < 0 or not i.is_integer:
                    continue
                i = int(i)
                if i < len(u0):
                    sol += u0[i]*x**i
                else:
                    sol += Symbol('C_%s' % j)*x**i
            sol = sol*x**constantpower
            if as_list:
                return [(sol.subs(x, x - x0),)] if x0 != 0 else [(sol,)]
            return sol.subs(x, x - x0) if x0 != 0 else sol
        if smallest_n + m > len(u0):
            raise NotImplementedError("Can't compute sufficient Initial Conditions")
        if any(c != 0 for c in coefficients[1:-1]):
            raise NotHyperSeriesError(self, self.x0)
        a, b = Poly(coefficients[0], nn), Poly(coefficients[-1], nn)
        c = -(a.LC()*m**(a.degree()))/(b.LC()*m**(b.degree()))
        sol = Integer(0)
        arg1 = roots(a.as_expr(), nn)
        arg2 = roots(b.as_expr(), nn)
        listofsol = []
        for i in range(smallest_n + m):
            if i < smallest_n:
                if as_list:
                    listofsol.append(((_expr(u0[i])*x**(i + constantpower)).subs(x, x - x0),))
                else:
                    sol += _expr(u0[i])*x**i
                continue
            if _expr(u0[i]) == 0:
                continue
            ap, bq = [], []
            for k in ordered(arg1.keys()):
                ap.extend([nsimplify((i - k)/m)]*arg1[k])
            for k in ordered(arg2.keys()):
                bq.extend([nsimplify((i - k)/m)]*arg2[k])
            if 1 in bq:
                bq.remove(1)
            else:
                ap.append(1)
            if as_list:
                listofsol.append(((_expr(u0[i])*x**(i + constantpower)).subs(x, x - x0),
                                  hyper(ap, bq, c*x**m).subs(x, x - x0)))
            else:
                sol += _expr(u0[i])*hyper(ap, bq, c*x**m)*x**i
        if as_list:
            return listofsol
        sol = sol*x**constantpower
        return sol.subs(x, x - x0) if x0 != 0 else sol

    def to_hyper(self, as_list=False, _recur=None):
        from .. import hyper, Rational
        try:
            return self._reference_hyper(as_list=as_list, _recur=_recur)
        except NotHyperSeriesError:
            # A known expression that is a sum of hypergeometric terms is
            # still converted (the reference raises).
            if self._expr is None or _recur is not None:
                raise
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
            raise NotHyperSeriesError(self, self.x0)
        expression = convert(self._expr)
        if as_list:
            return [(expression,)]
        return expression

    def to_meijerg(self):
        if self._expr is not None and str(self._expr.func) == 'meijerg':
            return self._expr
        # Reference: each hypergeometric term of to_hyper(as_list=True).
        from .. import Integer
        sol = Integer(0)
        for i in self._reference_hyper(as_list=True):
            if len(i) == 1:
                sol += i[0]
            elif len(i) == 2:
                sol += i[0]*_hyper_to_meijerg(i[1])
        return sol

    def change_ics(self, b, lenics=None):
        b = _expr(b)
        count = self.annihilator.order if lenics is None else operator.index(lenics)
        if isinstance(self.y0, dict) and self._expr is None:
            # Reference: through to_expr and expr_to_holonomic at b.
            return expr_to_holonomic(self.to_expr(), x=self.x, x0=b, lenics=lenics,
                                     domain=self.annihilator.parent.base.domain)
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
        if self.is_singularics() or self._singular_point():
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
    if y0 is None and initcond and _expr(x0) == 0 and result.annihilator.is_singular(0) \
            and (lenics is None or lenics <= 1):
        singular = _algebraic_singular_ics(func, x)
        if singular is not None:
            result.y0 = singular
            return result
    count = result._initial_count() if lenics is None else operator.index(lenics)
    if isinstance(y0, dict):
        result.y0 = {_expr(k): [_expr(v) for v in c] for k, c in y0.items()}
        result._expr = None
        return result
    result.y0 = [_expr(v) for v in y0] if y0 is not None else _conditions(func,x,_expr(x0),count) if initcond else None
    if y0 is not None and result.y0 != _conditions(func,x,_expr(x0),len(result.y0)):
        result._expr = None
    if y0 is None and initcond and (result.y0 is None or _fractional_exponent(result)):
        _undefined_conditions(result, func, x, lenics)
    return result


def _fractional_exponent(result):
    # A single simple non-integer indicial root: ordinary values at x0 do
    # not determine the function (x**(3/2)*exp(x) has y(0) = 0 for every
    # multiple), so the reference singular conditions are used.
    if not result.annihilator.is_singular(result.x0):
        return False
    r = result._indicial()
    keys = list(r)
    return len(r) == 1 and r[keys[0]] == 1 and not keys[0].is_integer


def _algebraic_singular_ics(func, x):
    # Reference _convert_poly_rat_alg for basepoly**(m/n) singular at 0:
    # {i*m/n: [c_i**(m/n)]} for the lowest nonzero coefficient c_i.
    from .. import Poly
    if not func.is_Pow:
        return None
    base, exponent = func.args
    if not (exponent.is_Rational and not exponent.is_Integer) or not base.is_polynomial(x):
        return None
    for i, (c,) in enumerate(reversed([[v] for v in Poly(base, x).all_coeffs()])):
        if c != 0:
            return {_expr(i)*exponent: [_expr(c)**exponent]}
    return None


def _undefined_conditions(result, func, x, lenics):
    # Reference expr_to_holonomic when the ordinary conditions are undefined
    # at x0: a single simple indicial root r gives the singular conditions
    # of func/(x - x0)**r, otherwise x0 moves to x0 + 1, x0 + 2, ...
    from .. import factorial
    x0 = result.x0
    count = result.annihilator.order if not lenics else lenics
    if result.annihilator.is_singular(x0):
        r = result._indicial()
        keys = list(r)
        if len(r) == 1 and r[keys[0]] == 1:
            g = func/(x - x0)**keys[0]
            ics = _find_conditions(g, x, x0, count)
            if ics:
                result.y0 = {keys[0]: [j/factorial(i) for i, j in enumerate(ics)]}
                return
    for _ in range(64):
        x0 += 1
        ics = _find_conditions(func, x, x0, count)
        if ics:
            result.x0, result.y0 = x0, ics
            return


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
    if isinstance(y0, dict):
        result.y0 = {_expr(k): [_expr(v) for v in c] for k, c in y0.items()}
        result._expr = None
        return result
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
    x = symbols[0]
    result = expr_to_holonomic(func,x=x,x0=x0,initcond=False)
    # SymPy's from_hyper: initial conditions of hyperexpand(func) at x0,
    # x0+1, ... without limits, or of the unexpanded function from x0 = 1.
    from .. import hyperexpand, oo
    simp = hyperexpand(func)
    if simp in (oo, -oo):
        return result
    order = result.annihilator.order
    if str(simp.func) == 'hyper':
        x0, evaluate = 1, evalf
    else:
        evaluate = False
    x0 = _expr(x0)
    y0 = _find_conditions(simp, x, x0, order, evaluate)
    while not y0:
        x0 += 1
        y0 = _find_conditions(simp, x, x0, order, evaluate)
    result.x0, result.y0 = x0, y0
    return result


def _find_conditions(func, x, x0, order, evalf=False):
    # SymPy's _find_conditions with use_limit=False.
    from .. import nan
    y0 = []
    for _ in range(order):
        value = func.subs(x, x0)
        if evalf:
            value = value.evalf()
        if value.is_finite is False or value is nan or value == nan:
            return None
        y0.append(value)
        func = func.diff(x)
    return y0


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


class BaseHolonomicError(Exception):
    def new(self, *args):
        raise NotImplementedError("abstract base class")


class NotPowerSeriesError(BaseHolonomicError):
    def __init__(self, holonomic, x0):
        self.holonomic = holonomic
        self.x0 = x0

    def __str__(self):
        return 'A Power Series does not exists for %s about %s.' % (self.holonomic, self.x0)


class NotHyperSeriesError(BaseHolonomicError):
    def __init__(self, holonomic, x0):
        self.holonomic = holonomic
        self.x0 = x0

    def __str__(self):
        return 'Power series expansion of %s about %s is not hypergeometric' % (self.holonomic, self.x0)


def _singular_combination(y1, y2, product):
    # Reference __add__/__mul__ with singular initial conditions: an
    # ordinary list becomes {0: [y0[i]/i!]}.
    from .. import factorial
    def as_dict(y):
        return y if isinstance(y, dict) else {_expr(0): [j/factorial(i) for i, j in enumerate(y)]}
    y1, y2 = as_dict(y1), as_dict(y2)
    y0 = {}
    if product:
        for i in y1:
            for j in y2:
                k = min(len(y1[i]), len(y2[j]))
                c = [sum((y1[i][b]*y2[j][a - b] for b in range(a + 1)), _expr(0)) for a in range(k)]
                y0[i + j] = c if i + j not in y0 else [a + b for a, b in zip(c, y0[i + j])]
        return y0
    for i in y1:
        y0[i] = [a + b for a, b in zip(y1[i], y2[i])] if i in y2 else y1[i]
    for i in y2:
        if i not in y1:
            y0[i] = y2[i]
    return y0


def _extend_y0(holonomic, n):
    """Reference _extend_y0: more initial values from the equation at an
    ordinary point."""
    from .. import cancel
    if holonomic.annihilator.is_singular(holonomic.x0) or holonomic.is_singularics() == True:
        return holonomic.y0
    x = holonomic.x
    a = holonomic.annihilator.order
    coefficients = holonomic.annihilator.coefficients
    y0 = holonomic.y0
    if len(y0) < a or n <= len(y0):
        return y0
    list_red = [cancel(-coefficients[i]/coefficients[a]) for i in range(a)]
    y1 = list(y0[:min(len(y0), a)])
    for _ in range(n - a):
        sol = _expr(0)
        for value, b in zip(y1, list_red):
            r = b.subs(x, holonomic.x0)
            if not getattr(r, 'is_finite', True):
                return y0
            sol += value*r
        y1.append(sol)
        derived = [cancel(list_red[0].diff(x))]
        for i, j in enumerate(list_red[1:]):
            derived.append(cancel(j.diff(x) + list_red[i]))
        derived.append(list_red[-1])
        list_red = derived
    return list(y0) + y1[len(y0):]


def _hyper_to_meijerg(func):
    """Reference _hyper_to_meijerg."""
    from .. import hyperexpand, gamma, meijerg, Integer
    ap, bq = func.ap, func.bq
    if any(i <= 0 and int(i) == i for i in ap):
        return hyperexpand(func)
    z = func.args[2]
    an = tuple(1 - i for i in ap)
    bm = (Integer(0),)
    bmq = tuple(1 - i for i in bq)
    k = Integer(1)
    for i in bq:
        k = k*gamma(i)
    for i in ap:
        k = k/gamma(i)
    return k*meijerg(an, (), bm, bmq, -z)
