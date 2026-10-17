"""Monomials as exponent tuples (``sympy.polys.monomials``).

``Monomial`` and the ``monomial_*`` helpers are plain tuple arithmetic.
``itermonomials`` builds products with the native ``Mul`` and, like SymPy,
yields the total-degree case from a set, so only the collection is specified,
not its iteration order.
"""
from itertools import combinations_with_replacement, product

__all__ = ["Monomial", "itermonomials", "monomial_count"]


def _is_sequence(value):
    return isinstance(value, (list, tuple)) or (
        hasattr(value, "__getitem__") and hasattr(value, "__len__") and not isinstance(value, str)
        and not hasattr(value, "is_Number"))


def itermonomials(variables, max_degrees, min_degrees=None):
    """Generate monomials in ``variables`` bounded by total or per-variable degrees."""
    from .. import Mul, S
    if _is_sequence(max_degrees):
        n = len(variables)
        if len(max_degrees) != n:
            raise ValueError('Argument sizes do not match')
        if min_degrees is None:
            min_degrees = [0]*n
        elif not _is_sequence(min_degrees):
            raise ValueError('min_degrees is not a list')
        else:
            if len(min_degrees) != n:
                raise ValueError('Argument sizes do not match')
            if any(i < 0 for i in min_degrees):
                raise ValueError("min_degrees cannot contain negative numbers")
        if any(min_degrees[i] > max_degrees[i] for i in range(n)):
            raise ValueError('min_degrees[i] must be <= max_degrees[i] for all i')
        power_lists = []
        for var, min_d, max_d in zip(variables, min_degrees, max_degrees):
            power_lists.append([var**i for i in range(min_d, max_d + 1)])
        for powers in product(*power_lists):
            yield Mul(*powers)
        return
    max_degree = max_degrees
    if max_degree < 0:
        raise ValueError("max_degrees cannot be negative")
    if min_degrees is None:
        min_degree = 0
    else:
        if min_degrees < 0:
            raise ValueError("min_degrees cannot be negative")
        min_degree = min_degrees
    if min_degree > max_degree:
        return
    if not variables or max_degree == 0:
        yield S.One
        return
    variables = list(variables) + [S.One]
    if all(getattr(v, "is_commutative", True) for v in variables):
        it = combinations_with_replacement(variables, int(max_degree))
    else:
        it = product(variables, repeat=int(max_degree))
    monomials_set = set()
    d = max_degree - min_degree
    for item in it:
        count = 0
        for variable in item:
            if variable == 1:
                count += 1
                if d < count:
                    break
        else:
            monomials_set.add(Mul(*item))
    yield from monomials_set


def monomial_count(V, N):
    """Number of monomials of total degree at most ``N`` in ``V`` variables."""
    from .. import factorial
    return factorial(V + N) / factorial(V) / factorial(N)


def monomial_mul(A, B):
    return tuple([a + b for a, b in zip(A, B)])


def monomial_ldiv(A, B):
    return tuple([a - b for a, b in zip(A, B)])


def monomial_div(A, B):
    C = monomial_ldiv(A, B)
    if all(c >= 0 for c in C):
        return tuple(C)
    return None


def monomial_pow(A, n):
    return tuple([a*n for a in A])


def monomial_gcd(A, B):
    return tuple([min(a, b) for a, b in zip(A, B)])


def monomial_lcm(A, B):
    return tuple([max(a, b) for a, b in zip(A, B)])


def monomial_divides(A, B):
    return all(a <= b for a, b in zip(A, B))


def monomial_max(*monoms):
    M = list(monoms[0])
    for N in monoms[1:]:
        for i, n in enumerate(N):
            M[i] = max(M[i], n)
    return tuple(M)


def monomial_min(*monoms):
    M = list(monoms[0])
    for N in monoms[1:]:
        for i, n in enumerate(N):
            M[i] = min(M[i], n)
    return tuple(M)


def monomial_deg(M):
    return sum(M)


def _exponents(other):
    from .. import Tuple
    if isinstance(other, Monomial):
        return other.exponents
    if isinstance(other, (tuple, Tuple)):
        return tuple(other)
    return None


class Monomial:
    """A monomial as a tuple of exponents, optionally tied to generators."""

    __slots__ = ('exponents', 'gens')

    def __init__(self, monom, gens=None):
        if not isinstance(monom, (tuple, list, Monomial)) and not (
                hasattr(monom, "__iter__") and not hasattr(monom, "is_Number")):
            from .. import Poly, sympify
            expr = sympify(monom)
            p = Poly(expr, *gens) if gens else Poly(expr)
            terms = p.terms()
            if len(terms) == 1 and terms[0][1] == 1:
                monom, gens = terms[0][0], p.gens
            else:
                raise ValueError("Expected a monomial got {}".format(monom))
        self.exponents = tuple(map(int, monom))
        self.gens = tuple(gens) if isinstance(gens, list) else gens

    def __getstate__(self):
        return (self.exponents, self.gens)

    def __setstate__(self, state):
        self.exponents, self.gens = state

    def rebuild(self, exponents, gens=None):
        return self.__class__(exponents, gens or self.gens)

    def __len__(self):
        return len(self.exponents)

    def __iter__(self):
        return iter(self.exponents)

    def __getitem__(self, item):
        return self.exponents[item]

    def __hash__(self):
        return hash((self.__class__.__name__, self.exponents, self.gens))

    def __str__(self):
        if self.gens:
            return "*".join(["%s**%s" % (gen, exp) for gen, exp in zip(self.gens, self.exponents)])
        return "%s(%s)" % (self.__class__.__name__, self.exponents)

    def as_expr(self, *gens):
        """Convert the monomial to an expression in ``gens`` (or its own gens)."""
        from .. import Mul
        gens = gens or self.gens
        if not gens:
            raise ValueError("Cannot convert %s to an expression without generators" % self)
        return Mul(*[gen**exp for gen, exp in zip(gens, self.exponents)])

    def __eq__(self, other):
        exponents = _exponents(other)
        if exponents is None:
            return False
        return self.exponents == exponents

    def __ne__(self, other):
        return not self == other

    def __mul__(self, other):
        exponents = _exponents(other)
        if exponents is None:
            raise NotImplementedError
        return self.rebuild(monomial_mul(self.exponents, exponents))

    def __truediv__(self, other):
        exponents = _exponents(other)
        if exponents is None:
            raise NotImplementedError
        result = monomial_div(self.exponents, exponents)
        if result is not None:
            return self.rebuild(result)
        from .polyerrors import ExactQuotientFailed
        raise ExactQuotientFailed("%s does not divide %s" % (Monomial(other), self))

    __floordiv__ = __truediv__

    def __pow__(self, other):
        n = int(other)
        if n < 0:
            raise ValueError("a non-negative integer expected, got %s" % other)
        return self.rebuild(monomial_pow(self.exponents, n))

    def gcd(self, other):
        """Greatest common divisor of monomials."""
        exponents = _exponents(other)
        if exponents is None:
            raise TypeError("an instance of Monomial class expected, got %s" % other)
        return self.rebuild(monomial_gcd(self.exponents, exponents))

    def lcm(self, other):
        """Least common multiple of monomials."""
        exponents = _exponents(other)
        if exponents is None:
            raise TypeError("an instance of Monomial class expected, got %s" % other)
        return self.rebuild(monomial_lcm(self.exponents, exponents))
