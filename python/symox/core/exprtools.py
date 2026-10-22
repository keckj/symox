"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'factor_nc factor_terms gcd_terms')

import symox as _S


def decompose_power(expr):
    expr=_S.sympify(expr);base,exponent=expr.as_base_exp()
    if exponent.is_Rational:return base**_S.Rational(1,exponent.q),int(exponent.p)
    coefficient,tail=exponent.as_coeff_Mul()
    if coefficient.is_Rational:return base**(tail/coefficient.q),int(coefficient.p)
    return expr,1

def decompose_power_rat(expr):
    expr=_S.sympify(expr);base,exponent=expr.as_base_exp()
    return (base,exponent) if exponent.is_Rational else decompose_power(expr)


from collections import defaultdict

from symox import S, I, Integer, Rational, Number, Basic, Expr, Mul, sympify
from symox.core.sorting import ordered
from symox.core.numbers import equal_valued
from symox.core.mul import _keep_coeff

SYMPY_INTS = (int,)


class NonCommutativeExpression(ValueError):
    """Raised when a commutative term is required."""


def _isnumber(i):
    return isinstance(i, (SYMPY_INTS, float)) or i.is_Number


# Factors and Term follow sympy/core/exprtools.py (SymPy 1.14) line for line:
# exponents may be negative or symbolic, and gcd/lcm/div/normal compare them
# exactly as the reference does.
class Factors:

    __slots__ = ('factors', 'gens')

    def __init__(self, factors=None):  # Factors
        if isinstance(factors, (SYMPY_INTS, float)):
            factors = S(factors)
        if isinstance(factors, Factors):
            factors = factors.factors.copy()
        elif factors in (None, S.One):
            factors = {}
        elif factors is S.Zero or factors == 0:
            factors = {S.Zero: S.One}
        elif isinstance(factors, Number):
            n = factors
            factors = {}
            if n < 0:
                factors[S.NegativeOne] = S.One
                n = -n
            if n is not S.One:
                if n.is_Float or n.is_Integer or n is S.Infinity:
                    factors[n] = S.One
                elif n.is_Rational:
                    # since we're processing Numbers, the denominator is
                    # stored with a negative exponent; all other factors
                    # are left .
                    if n.p != 1:
                        factors[Integer(n.p)] = S.One
                    factors[Integer(n.q)] = S.NegativeOne
                else:
                    raise ValueError('Expected Float|Rational|Integer, not %s' % n)
        elif isinstance(factors, Basic) and not factors.args:
            factors = {factors: S.One}
        elif isinstance(factors, Expr):
            c, nc = factors.args_cnc()
            i = c.count(I)
            for _ in range(i):
                c.remove(I)
            factors = dict(Mul._from_args(c).as_powers_dict())
            # Handle all rational Coefficients
            for f in list(factors.keys()):
                if isinstance(f, Rational) and not isinstance(f, Integer):
                    p, q = Integer(f.p), Integer(f.q)
                    factors[p] = (factors[p] if p in factors else S.Zero) + factors[f]
                    factors[q] = (factors[q] if q in factors else S.Zero) - factors[f]
                    factors.pop(f)
            if i:
                factors[I] = factors.get(I, S.Zero) + i
            if nc:
                factors[Mul(*nc, evaluate=False)] = S.One
        else:
            factors = factors.copy()  # /!\ should be dict-like

            # tidy up -/+1 and I exponents if Rational

            handle = [k for k in factors if k is I or k in (-1, 1)]
            if handle:
                i1 = S.One
                for k in handle:
                    if not _isnumber(factors[k]):
                        continue
                    i1 *= k**factors.pop(k)
                if i1 is not S.One:
                    for a in i1.args if i1.is_Mul else [i1]:  # at worst, -1.0*I*(-1)**e
                        if a is S.NegativeOne:
                            factors[a] = S.One
                        elif a is I:
                            factors[I] = S.One
                        elif a.is_Pow:
                            factors[a.base] = factors.get(a.base, S.Zero) + a.exp
                        elif equal_valued(a, 1):
                            factors[a] = S.One
                        elif equal_valued(a, -1):
                            factors[-a] = S.One
                            factors[S.NegativeOne] = S.One
                        else:
                            raise ValueError('unexpected factor in i1: %s' % a)

        self.factors = factors
        keys = getattr(factors, 'keys', None)
        if keys is None:
            raise TypeError('expecting Expr or dictionary')
        self.gens = frozenset(keys())

    def __hash__(self):  # Factors
        keys = tuple(ordered(self.factors.keys()))
        values = [self.factors[k] for k in keys]
        return hash((keys, values))

    def __repr__(self):  # Factors
        return "Factors({%s})" % ', '.join(
            ['%s: %s' % (k, v) for k, v in ordered(self.factors.items())])

    @property
    def is_zero(self):  # Factors
        f = self.factors
        return len(f) == 1 and S.Zero in f

    @property
    def is_one(self):  # Factors
        return not self.factors

    def as_expr(self):  # Factors

        args = []
        for factor, exp in self.factors.items():
            if exp != 1:
                if isinstance(exp, Integer):
                    b, e = factor.as_base_exp()
                    e = _keep_coeff(exp, e)
                    args.append(b**e)
                else:
                    args.append(factor**exp)
            else:
                args.append(factor)
        return Mul(*args)

    def mul(self, other):  # Factors
        if not isinstance(other, Factors):
            other = Factors(other)
        if any(f.is_zero for f in (self, other)):
            return Factors(S.Zero)
        factors = dict(self.factors)

        for factor, exp in other.factors.items():
            if factor in factors:
                exp = factors[factor] + exp

                if not exp:
                    del factors[factor]
                    continue

            factors[factor] = exp

        return Factors(factors)

    def normal(self, other):
        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                return (Factors(), Factors(S.Zero))
            if self.is_zero:
                return (Factors(S.Zero), Factors())

        self_factors = dict(self.factors)
        other_factors = dict(other.factors)

        for factor, self_exp in self.factors.items():
            try:
                other_exp = other.factors[factor]
            except KeyError:
                continue

            exp = self_exp - other_exp

            if not exp:
                del self_factors[factor]
                del other_factors[factor]
            elif _isnumber(exp):
                if exp > 0:
                    self_factors[factor] = exp
                    del other_factors[factor]
                else:
                    del self_factors[factor]
                    other_factors[factor] = -exp
            else:
                r = self_exp.extract_additively(other_exp)
                if r is not None:
                    if r:
                        self_factors[factor] = r
                        del other_factors[factor]
                    else:  # should be handled already
                        del self_factors[factor]
                        del other_factors[factor]
                else:
                    sc, sa = self_exp.as_coeff_Add()
                    if sc:
                        oc, oa = other_exp.as_coeff_Add()
                        diff = sc - oc
                        if diff > 0:
                            self_factors[factor] -= oc
                            other_exp = oa
                        elif diff < 0:
                            self_factors[factor] -= sc
                            other_factors[factor] -= sc
                            other_exp = oa - diff
                        else:
                            self_factors[factor] = sa
                            other_exp = oa
                    if other_exp:
                        other_factors[factor] = other_exp
                    else:
                        del other_factors[factor]

        return Factors(self_factors), Factors(other_factors)

    def div(self, other):  # Factors
        quo, rem = dict(self.factors), {}

        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                raise ZeroDivisionError
            if self.is_zero:
                return (Factors(S.Zero), Factors())

        for factor, exp in other.factors.items():
            if factor in quo:
                d = quo[factor] - exp
                if _isnumber(d):
                    if d <= 0:
                        del quo[factor]

                    if d >= 0:
                        if d:
                            quo[factor] = d

                        continue

                    exp = -d

                else:
                    r = quo[factor].extract_additively(exp)
                    if r is not None:
                        if r:
                            quo[factor] = r
                        else:  # should be handled already
                            del quo[factor]
                    else:
                        other_exp = exp
                        sc, sa = quo[factor].as_coeff_Add()
                        if sc:
                            oc, oa = other_exp.as_coeff_Add()
                            diff = sc - oc
                            if diff > 0:
                                quo[factor] -= oc
                                other_exp = oa
                            elif diff < 0:
                                quo[factor] -= sc
                                other_exp = oa - diff
                            else:
                                quo[factor] = sa
                                other_exp = oa
                        if other_exp:
                            rem[factor] = other_exp
                        else:
                            assert factor not in rem
                    continue

            rem[factor] = exp

        return Factors(quo), Factors(rem)

    def quo(self, other):  # Factors
        return self.div(other)[0]

    def rem(self, other):  # Factors
        return self.div(other)[1]

    def pow(self, other):  # Factors
        if isinstance(other, Factors):
            other = other.as_expr()
            if other.is_Integer:
                other = int(other)
        if isinstance(other, SYMPY_INTS) and other >= 0:
            factors = {}

            if other:
                for factor, exp in self.factors.items():
                    factors[factor] = exp*other

            return Factors(factors)
        else:
            raise ValueError("expected non-negative integer, got %s" % other)

    def gcd(self, other):  # Factors
        if not isinstance(other, Factors):
            other = Factors(other)
            if other.is_zero:
                return Factors(self.factors)

        factors = {}

        for factor, exp in self.factors.items():
            factor, exp = sympify(factor), sympify(exp)
            if factor in other.factors:
                lt = (exp - other.factors[factor]).is_negative
                if lt == True:
                    factors[factor] = exp
                elif lt == False:
                    factors[factor] = other.factors[factor]

        return Factors(factors)

    def lcm(self, other):  # Factors
        if not isinstance(other, Factors):
            other = Factors(other)
            if any(f.is_zero for f in (self, other)):
                return Factors(S.Zero)

        factors = dict(self.factors)

        for factor, exp in other.factors.items():
            if factor in factors:
                exp = max(exp, factors[factor])

            factors[factor] = exp

        return Factors(factors)

    def __mul__(self, other):  # Factors
        return self.mul(other)

    def __divmod__(self, other):  # Factors
        return self.div(other)

    def __truediv__(self, other):  # Factors
        return self.quo(other)

    def __mod__(self, other):  # Factors
        return self.rem(other)

    def __pow__(self, other):  # Factors
        return self.pow(other)

    def __eq__(self, other):  # Factors
        if not isinstance(other, Factors):
            other = Factors(other)
        return self.factors == other.factors

    def __ne__(self, other):  # Factors
        return not self == other



def _coeff_gcd(a, b):
    # Rational.gcd of the reference: gcd of numerators over lcm of denominators.
    from symox import gcd
    return a.gcd(b) if hasattr(a, 'gcd') else gcd(a, b)


def _coeff_lcm(a, b):
    from symox import lcm
    return a.lcm(b) if hasattr(a, 'lcm') else lcm(a, b)

class Term:

    __slots__ = ('coeff', 'numer', 'denom')

    def __init__(self, term, numer=None, denom=None):  # Term
        if numer is None and denom is None:
            if not term.is_commutative:
                raise NonCommutativeExpression(
                    'commutative expression expected')

            coeff, factors = term.as_coeff_mul()
            numer, denom = defaultdict(int), defaultdict(int)

            for factor in factors:
                base, exp = decompose_power(factor)

                if base.is_Add:
                    cont, base = base.primitive()
                    coeff *= cont**exp

                if exp > 0:
                    numer[base] += exp
                else:
                    denom[base] += -exp

            numer = Factors(numer)
            denom = Factors(denom)
        else:
            coeff = term

            if numer is None:
                numer = Factors()

            if denom is None:
                denom = Factors()

        self.coeff = coeff
        self.numer = numer
        self.denom = denom

    def __hash__(self):  # Term
        return hash((self.coeff, self.numer, self.denom))

    def __repr__(self):  # Term
        return "Term(%s, %s, %s)" % (self.coeff, self.numer, self.denom)

    def as_expr(self):  # Term
        return self.coeff*(self.numer.as_expr()/self.denom.as_expr())

    def mul(self, other):  # Term
        coeff = self.coeff*other.coeff
        numer = self.numer.mul(other.numer)
        denom = self.denom.mul(other.denom)

        numer, denom = numer.normal(denom)

        return Term(coeff, numer, denom)

    def inv(self):  # Term
        return Term(1/self.coeff, self.denom, self.numer)

    def quo(self, other):  # Term
        return self.mul(other.inv())

    def pow(self, other):  # Term
        if other < 0:
            return self.inv().pow(-other)
        else:
            return Term(self.coeff ** other,
                        self.numer.pow(other),
                        self.denom.pow(other))

    def gcd(self, other):  # Term
        return Term(_coeff_gcd(self.coeff, other.coeff),
                    self.numer.gcd(other.numer),
                    self.denom.gcd(other.denom))

    def lcm(self, other):  # Term
        return Term(_coeff_lcm(self.coeff, other.coeff),
                    self.numer.lcm(other.numer),
                    self.denom.lcm(other.denom))

    def __mul__(self, other):  # Term
        if isinstance(other, Term):
            return self.mul(other)
        else:
            return NotImplemented

    def __truediv__(self, other):  # Term
        if isinstance(other, Term):
            return self.quo(other)
        else:
            return NotImplemented

    def __pow__(self, other):  # Term
        if isinstance(other, SYMPY_INTS):
            return self.pow(other)
        else:
            return NotImplemented

    def __eq__(self, other):  # Term
        return (self.coeff == other.coeff and
                self.numer == other.numer and
                self.denom == other.denom)

    def __ne__(self, other):  # Term
        return not self == other

__all__ += ['decompose_power','decompose_power_rat','Factors','Term']
