"""Scalar coefficient results paired with their native domain."""

import operator


class DomainScalar:
    def __init__(self, element, domain):
        self.element = domain.convert(element)
        self.domain = domain

    def to_sympy(self):
        return self.domain.to_sympy(self.element)

    _symox_ = to_sympy

    def convert_to(self, domain):
        return type(self)(domain.convert(self.element), domain)

    def unify(self, other):
        if not isinstance(other, DomainScalar):
            raise TypeError('Expected a DomainScalar')
        domain = self.domain.unify(other.domain)
        return self.convert_to(domain), other.convert_to(domain)

    def __eq__(self, other):
        if not isinstance(other, DomainScalar):
            return NotImplemented
        return self.domain == other.domain and self.element == other.element

    def __bool__(self):
        return not self.domain.is_zero(self.element)

    def __repr__(self):
        return 'DomainScalar(%r, %s)' % (self.element, self.domain)

    def __neg__(self):
        return type(self)(-self.element, self.domain)

    def _binary(self, other, operation):
        if not isinstance(other, DomainScalar):
            return NotImplemented
        first, second = self.unify(other)
        return type(self)(operation(first.element, second.element), first.domain)

    def __add__(self, other):
        return self._binary(other, operator.add)

    def __sub__(self, other):
        return self._binary(other, operator.sub)

    def __mul__(self, other):
        return self._binary(other, operator.mul)

    def __pow__(self, exponent):
        return type(self)(self.element**operator.index(exponent), self.domain)
