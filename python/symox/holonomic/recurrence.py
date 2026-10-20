"""Polynomial Ore algebras and polynomial-coefficient recurrences.

Coefficients use Symox's native exact polynomial domains.  Multiplication
implements the Ore relation, rather than treating the generator as a symbol.
"""
import math
import operator


def _expr(value):
    from .. import sympify
    return value.as_expr() if hasattr(value, 'as_expr') else sympify(value)


class _ParameterizedPolynomialRing:
    """QQ(parameters)[x], stored in the native multivariate rational field.

    Native domain parents currently disallow nested composite storage. Validate
    polynomial dependence on x here, retaining native exact coefficient data.
    """
    def __init__(self,domain,x):
        self.domain,self.symbols = domain,(x,)
        self._field = domain.domain.frac_field(*domain.symbols,x)
        self.zero,self.one = self.convert(0),self.convert(1)

    def convert(self,value):
        from .. import cancel
        value = cancel(_expr(value))
        if not value.is_polynomial(self.symbols[0]):
            raise ValueError('Ore coefficients must be polynomial in the independent variable')
        return self._field.convert(value)

    def __eq__(self,other):
        return isinstance(other,_ParameterizedPolynomialRing) and (self.domain,self.symbols) == (other.domain,other.symbols)

    def __repr__(self):
        return '%s[%s]' % (self.domain,self.symbols[0])


def _base(domain,x):
    from .. import QQ,ZZ
    try:
        return domain.old_poly_ring(x)
    except NotImplementedError:
        if domain.is_FractionField and domain.domain in (QQ,ZZ):
            return _ParameterizedPolynomialRing(domain,x)
        raise


class _OreAlgebra:
    differential = False

    def __init__(self, base, generator=None):
        from .. import Symbol
        self.base = base
        symbols = tuple(base.symbols)
        if len(symbols) != 1:
            raise ValueError('Ore algebras require a univariate polynomial ring')
        self.x = symbols[0]
        # Display label only: noncommutativity is implemented by the Ore product.
        self.gen_symbol = Symbol(str(generator or ('Dx' if self.differential else 'Sn')))
        self.operator = self.operator_class([0, 1], self)
        if self.differential:
            self.derivative_operator = self.operator
        else:
            self.shift_operator = self.operator

    def __eq__(self, other):
        return type(self) is type(other) and self.base == other.base and self.gen_symbol == other.gen_symbol

    def __repr__(self):
        return '%s Operator Algebra in %s over %s' % ('Differential' if self.differential else 'Recurrence', self.gen_symbol, self.base)

    def __call__(self, value):
        if isinstance(value, _OreOperator):
            if value.parent != self:
                raise ValueError('operators belong to different algebras')
            return value
        return self.operator_class([value], self)


class _OreOperator:
    def __init__(self, list_of_poly, parent):
        self.parent = parent
        values = list(list_of_poly)
        while len(values) > 1 and _expr(values[-1]) == 0:
            values.pop()
        self.listofpoly = [parent.base.convert(_expr(c)) for c in values or [0]]
        self.order = len(self.listofpoly) - 1

    @property
    def coefficients(self):
        return [_expr(p) for p in self.listofpoly]

    def __repr__(self):
        g = str(self.parent.gen_symbol)
        return ' + '.join('(%s)%s' % (c, '' if i == 0 else g if i == 1 else g+'**'+str(i))
                          for i, c in enumerate(self.coefficients) if c != 0) or '0'

    def __eq__(self, other):
        try:
            other = self.parent(other)
        except (TypeError, ValueError):
            return False
        return self.coefficients == other.coefficients

    def __add__(self, other):
        other = self.parent(other)
        a, b = self.coefficients, other.coefficients
        return type(self)([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                           for i in range(max(len(a), len(b)))], self.parent)

    __radd__ = __add__

    def __neg__(self):
        return type(self)([-c for c in self.coefficients], self.parent)

    def __sub__(self, other):
        return self + -self.parent(other)

    def __rsub__(self, other):
        return self.parent(other) + -self

    def __mul__(self, other):
        from .. import Integer, expand
        other = self.parent(other)
        x = self.parent.x
        result = [Integer(0)] * (self.order + other.order + 1)
        for i, a in enumerate(self.coefficients):
            for j, b in enumerate(other.coefficients):
                if self.parent.differential:
                    for k in range(i + 1):
                        result[i+j-k] += a * math.comb(i, k) * b.diff(x, k)
                else:
                    result[i+j] += a * b.subs(x, x+i)
        return type(self)([expand(c) for c in result], self.parent)

    def __rmul__(self, other):
        return self.parent(other) * self

    def __pow__(self, power):
        power = operator.index(power)
        if power < 0:
            raise ValueError('Ore operator powers must be nonnegative')
        result, base = self.parent(1), self
        while power:
            if power & 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def __truediv__(self, scalar):
        from .. import cancel
        if isinstance(scalar, _OreOperator):
            raise NotImplementedError('operator division requires a side and a remainder')
        scalar = _expr(scalar)
        if scalar == 0:
            raise ZeroDivisionError('division by zero')
        return type(self)([cancel(c / scalar) for c in self.coefficients], self.parent)


class RecurrenceOperator(_OreOperator):
    """An operator sum p_i(n) S**i, where S f(n) = f(n+1)."""


class RecurrenceOperatorAlgebra(_OreAlgebra):
    operator_class = RecurrenceOperator


def RecurrenceOperators(base, generator=None):
    parent = RecurrenceOperatorAlgebra(base, generator)
    return parent, parent.shift_operator


class HolonomicSequence:
    def __init__(self, recurrence, u0=None):
        if not isinstance(recurrence, RecurrenceOperator):
            raise TypeError('expected a recurrence operator')
        self.recurrence = recurrence
        self.u0 = [_expr(v) for v in (u0 or [])]
        self.n = recurrence.parent.x

    def __repr__(self):
        return 'HolonomicSequence(%r, %r)' % (self.recurrence, self.u0)

    def __eq__(self, other):
        return isinstance(other, HolonomicSequence) and self.recurrence == other.recurrence and self.u0 == other.u0

    def eval(self, n):
        from .. import cancel
        n = operator.index(n)
        if n < 0:
            raise ValueError('sequence indices must be nonnegative')
        values = list(self.u0)
        r = self.recurrence.order
        if len(values) < r:
            raise ValueError('insufficient initial values')
        if r == 0 and n >= len(values):
            raise ValueError('order-zero recurrences do not specify new terms')
        coeff = self.recurrence.coefficients
        while len(values) <= n:
            start = len(values) - r
            lead = coeff[-1].subs(self.n, start)
            if lead == 0:
                raise ValueError('recurrence is singular at index %s; supply more initial values' % start)
            values.append(cancel(-sum(coeff[i].subs(self.n, start)*values[start+i] for i in range(r))/lead))
        return values[n]


def _install_scalar_dispatch():
    from .._core import Expr
    from functools import wraps
    if getattr(Expr, '_ore_dispatch_installed', False):
        return
    for name in ('__add__', '__sub__', '__mul__', '__truediv__'):
        original = getattr(Expr, name)
        def factory(original, name):
            @wraps(original)
            def operation(self, other):
                if isinstance(other, _OreOperator) or getattr(other,'_symox_holonomic_function',False):
                    if name == '__add__':
                        return other.__radd__(self)
                    if name == '__sub__':
                        return other.__rsub__(self)
                    if name == '__mul__':
                        return other.__rmul__(self)
                    raise NotImplementedError('scalar divided by an Ore operator')
                return original(self, other)
            return operation
        setattr(Expr, name, factory(original, name))
    Expr._ore_dispatch_installed = True


_install_scalar_dispatch()
