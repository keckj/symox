"""Finite free algebra over ZZ or QQ modulo a monic univariate polynomial."""
from ..domains._modules import Ideal, QuotientRing, QuotientElement


class _MonicIdeal(Ideal):
    def __init__(self, ring, modulus):
        self.ring = ring
        self.gens = (ring.convert(modulus.as_expr()),)
        self._generator = self.gens[0]
        self._basis = None


class ExtensionElement(QuotientElement):
    def __init__(self, rep, ext):
        super().__init__(ext, ext.base_ideal.reduce_element(rep))

    @property
    def ext(self):
        return self.parent()

    @property
    def rep(self):
        from ... import Poly
        return Poly(self.ext.ring.to_sympy(self.data), self.ext.symbol, domain=self.ext.domain)

    def _symox_(self):
        from ... import CoercionFailed
        expression = self.ext.ring.to_sympy(self.data)
        if self.ext.symbol in expression.free_symbols:
            raise CoercionFailed('an abstract extension element has no canonical scalar embedding')
        return expression


class MonogenicFiniteExtension(QuotientRing):
    _symox_monogenic_extension = True

    def __init__(self, modulus):
        from ... import Poly
        if not isinstance(modulus, Poly) or not modulus.is_univariate:
            raise TypeError('a monogenic extension requires a univariate Poly')
        if not modulus.domain.is_ZZ and not modulus.domain.is_QQ:
            raise NotImplementedError('finite extensions currently require ZZ or QQ coefficients')
        if modulus.is_zero:
            raise ValueError('the defining polynomial must be nonzero')
        if modulus.LC() != 1:
            modulus = modulus.exquo_ground(modulus.LC())
        degree = int(modulus.degree())
        if degree < 1:
            raise ValueError('the defining polynomial must have positive degree')
        if degree > 4096:
            raise NotImplementedError('finite extension rank is limited to 4096')
        self.modulus, self.rank = modulus, degree
        self.domain, self.symbol = modulus.domain, modulus.gens[0]
        ring = self.domain.poly_ring(self.symbol)
        super().__init__(ring, _MonicIdeal(ring, modulus))

    @property
    def generator(self):
        return self(self.ring.gens[0])

    @property
    def basis(self):
        return tuple(self.generator**power for power in range(self.rank))

    @property
    def is_Field(self):
        return self.domain.is_Field and self.modulus.is_irreducible

    @property
    def is_Exact(self):
        return self.domain.is_Exact

    def characteristic(self):
        return self.domain.characteristic()

    def convert(self, value):
        from ... import Poly
        if isinstance(value, ExtensionElement):
            if value.parent() != self:
                raise ValueError('incompatible finite extensions')
            return value
        if isinstance(value, Poly):
            value = value.as_expr()
        return ExtensionElement(value, self)

    __call__ = convert

    def revert(self, value):
        from ... import CoercionFailed
        from ..polyerrors import NotReversible
        value = self.convert(value)
        field_ring = self.domain.get_field().poly_ring(self.symbol)
        try:
            inverse = field_ring.invert(field_ring.convert(self.ring.to_sympy(value.data)),
                                        field_ring.convert(self.modulus.as_expr()))
            result = self(field_ring.to_sympy(inverse))
        except CoercionFailed:
            raise NotReversible('inverse coefficients do not belong to the base ring') from None
        if result*value != self.one:
            raise NotReversible('the extension element is not a unit')
        return result
