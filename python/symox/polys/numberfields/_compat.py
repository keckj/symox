"""Small reference helpers used by the ported number-field modules.

``IntegerPowerable`` and ``dup_clear_denoms`` follow
sympy/polys/polyutils.py and sympy/polys/densetools.py; ``ANP`` keeps the
dense coefficient list and modulus of an algebraic number, as the reference
``ANP`` does, for the module conversions.
"""
from ... import ilcm


class IntegerPowerable:
    def __pow__(self, e, modulo=None):
        if e < 2:
            try:
                if e == 1:
                    return self
                elif e == 0:
                    return self._zeroth_power()
            except NotImplementedError:
                pass
            return NotImplemented
        else:
            bits = [int(d) for d in reversed(bin(e)[2:])]
            n = len(bits)
            p = self
            first = True
            for i in range(n):
                if bits[i]:
                    if first:
                        r = p
                        first = False
                    else:
                        r *= p
                        if modulo is not None:
                            r %= modulo
                if i < n - 1:
                    p *= p
                    if modulo is not None:
                        p %= modulo
            return r

    def _first_power(self):
        raise NotImplementedError

    def _zeroth_power(self):
        raise NotImplementedError


def dup_clear_denoms(f, K0, K1=None, convert=False):
    """Clear the denominators of a dense rational list, as the reference."""
    from ..domains import ZZ
    common = 1
    for c in f:
        common = ilcm(common, int(K0.denom(c)))
    common = K0.convert(common)
    if common != K0.one:
        f = [c * common for c in f]
    if not convert:
        return common, f
    return ZZ.convert(common), [ZZ.convert(c) for c in f]


class ANP:
    """Dense representation of an algebraic number modulo ``mod``."""

    def __init__(self, rep, mod, dom):
        self._rep = [dom.convert(c) for c in rep]
        while self._rep and not self._rep[0]:
            self._rep.pop(0)
        self.mod = [dom.convert(c) for c in mod]
        self.dom = dom

    def to_list(self):
        return list(self._rep)

    def mod_to_list(self):
        return list(self.mod)

    @property
    def rep(self):
        return list(self._rep)

    def __eq__(self, other):
        return (isinstance(other, ANP) and self._rep == other._rep
                and self.mod == other.mod)

    def __hash__(self):
        return hash((type(self).__name__, tuple(self._rep), tuple(self.mod)))

    def __repr__(self):
        return '%s(%s, %s, %s)' % (type(self).__name__, self._rep, self.mod, self.dom)


def public(obj):
    """The reference ``@public`` decorator; module ``__all__`` lists exports."""
    return obj
