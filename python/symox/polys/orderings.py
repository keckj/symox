"""Monomial orderings as sort-key callables (``sympy.polys.orderings``).

An ordering maps an exponent tuple to a key whose natural Python ordering is
the monomial ordering. Native polynomial methods recognise the six named
orderings through their ``alias`` and sort exponent vectors without calling
back into Python; any other callable is used as a Python sort key.
"""
__all__ = ["lex", "grlex", "grevlex", "ilex", "igrlex", "igrevlex"]


def _iterable(value):
    return isinstance(value, (tuple, list))


class MonomialOrder:
    """Base class for monomial orderings."""

    alias = None
    is_global = None
    is_default = False

    def __repr__(self):
        return self.__class__.__name__ + "()"

    def __str__(self):
        return self.alias

    def __call__(self, monomial):
        raise NotImplementedError

    def __eq__(self, other):
        return self.__class__ == other.__class__

    def __hash__(self):
        return hash(self.__class__)

    def __ne__(self, other):
        return not (self == other)

    def __reduce__(self):
        return (_named_order, (str(self),))


class LexOrder(MonomialOrder):
    """Lexicographic order of monomials."""

    alias = 'lex'
    is_global = True
    is_default = True

    def __call__(self, monomial):
        return monomial


class GradedLexOrder(MonomialOrder):
    """Graded lexicographic order of monomials."""

    alias = 'grlex'
    is_global = True

    def __call__(self, monomial):
        return (sum(monomial), monomial)


class ReversedGradedLexOrder(MonomialOrder):
    """Reversed graded lexicographic order of monomials."""

    alias = 'grevlex'
    is_global = True

    def __call__(self, monomial):
        return (sum(monomial), tuple(reversed([-m for m in monomial])))


class ProductOrder(MonomialOrder):
    """A product order built from ``(order, selector)`` pairs."""

    def __init__(self, *args):
        self.args = args

    def __call__(self, monomial):
        return tuple(O(lamda(monomial)) for (O, lamda) in self.args)

    def __repr__(self):
        return self.__class__.__name__ + '(' + ", ".join(repr(x[0]) for x in self.args) + ')'

    def __str__(self):
        return self.__class__.__name__ + '(' + ", ".join(str(x[0]) for x in self.args) + ')'

    def __eq__(self, other):
        return isinstance(other, ProductOrder) and self.args == other.args

    def __hash__(self):
        return hash((self.__class__, self.args))

    def __reduce__(self):
        return (ProductOrder, self.args)

    @property
    def is_global(self):
        if all(o.is_global is True for o, _ in self.args):
            return True
        if all(o.is_global is False for o, _ in self.args):
            return False
        return None


class InverseOrder(MonomialOrder):
    """The inverse ``iO`` of an order ``O``: ``A >_iO B`` iff ``B >_O A``."""

    def __init__(self, O):
        self.O = O

    def __str__(self):
        return "i" + str(self.O)

    def __call__(self, monomial):
        def inv(l):
            if _iterable(l):
                return tuple(inv(x) for x in l)
            return -l
        return inv(self.O(monomial))

    @property
    def is_global(self):
        if self.O.is_global is True:
            return False
        if self.O.is_global is False:
            return True
        return None

    def __eq__(self, other):
        return isinstance(other, InverseOrder) and other.O == self.O

    def __hash__(self):
        return hash((self.__class__, self.O))

    def __reduce__(self):
        if str(self) in _monomial_key:
            return (_named_order, (str(self),))
        return (InverseOrder, (self.O,))


lex = LexOrder()
grlex = GradedLexOrder()
grevlex = ReversedGradedLexOrder()
ilex = InverseOrder(lex)
igrlex = InverseOrder(grlex)
igrevlex = InverseOrder(grevlex)

_monomial_key = {
    'lex': lex,
    'grlex': grlex,
    'grevlex': grevlex,
    'ilex': ilex,
    'igrlex': igrlex,
    'igrevlex': igrevlex,
}


def _named_order(name):
    return _monomial_key[name]


def monomial_key(order=None, gens=None):
    """Return a sort key implementing ``order`` on exponent tuples.

    ``order`` is a name (``'lex'``, ``'grlex'``, ``'grevlex'`` or an inverse
    ``'i...'`` order), a Symbol with such a name, or any callable. With
    ``gens`` the key accepts expressions instead of exponent tuples.
    """
    from .. import Symbol
    if order is None:
        order = lex
    if isinstance(order, Symbol):
        order = str(order)
    if isinstance(order, str):
        try:
            order = _monomial_key[order]
        except KeyError:
            raise ValueError("supported monomial orderings are 'lex', 'grlex' and 'grevlex', got %r" % order)
    if hasattr(order, '__call__'):
        if gens is not None:
            def _order(expr):
                from .. import Poly
                return order(Poly(expr, *gens).degree_list())
            return _order
        return order
    raise ValueError("monomial ordering specification must be a string or a callable, got %s" % order)


class _ItemGetter:
    """Helper picking exponents of selected generators."""

    def __init__(self, seq):
        self.seq = tuple(seq)

    def __call__(self, m):
        return tuple(m[idx] for idx in self.seq)

    def __eq__(self, other):
        return isinstance(other, _ItemGetter) and self.seq == other.seq

    def __hash__(self):
        return hash(self.seq)


def build_product_order(arg, gens):
    """Build a product monomial order on ``gens`` from ``(order, *vars)`` items."""
    gens2idx = {g: i for i, g in enumerate(gens)}
    order = []
    for expr in arg:
        name, var = expr[0], expr[1:]
        order.append((monomial_key(name), _ItemGetter(gens2idx[g] for g in var)))
    return ProductOrder(*order)
