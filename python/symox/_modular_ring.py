"""The ring Z/nZ that SymPy's ``GF(n)`` builds for moduli the native field
does not cover.

SymPy 1.14 accepts any positive modulus in ``GF(n)``: a composite ``n``
gives the ring of integers modulo ``n`` (``is_Field`` is then False) and the
nonsymmetric option changes only the element representatives. The native
``GF(p)`` covers symmetric prime fields below ``2**32``. This module ports
the reference ``FiniteField`` and ``ModularInteger`` element arithmetic for
the other moduli: composite moduli, primes above the native range and the
nonsymmetric representation. Polynomials over these domains are not
implemented (see docs/algorithmic-parity.md).
"""
import operator

from . import _core
from ._core import CoercionFailed, DomainError, NotInvertible


def _isprime(n):
    from .ntheory import isprime
    return bool(isprime(n))


def _invert(a, m):
    """``ZZ.invert(a, m)``: the inverse modulo ``m`` or ``NotInvertible``."""
    s, t, h = _gcdex(a, m)
    if h == 1:
        return s % m
    raise NotInvertible("zero divisor")


def _gcdex(a, b):
    # ZZ.gcdex (Python ints): s*a + t*b == h with h >= 0.
    if not b:
        return (1 if a >= 0 else -1, 0, abs(a))
    if not a:
        return (0, 1 if b >= 0 else -1, abs(b))
    x_sign = -1 if a < 0 else 1
    y_sign = -1 if b < 0 else 1
    a, b = abs(a), abs(b)
    x, y, r, s = 1, 0, 0, 1
    while b:
        c, q = a % b, a // b
        a, b, r, s, x, y = b, c, x - q*r, y - q*s, r, s
    return (x*x_sign, y*y_sign, a)


def _int_of(value):
    """The integer of a native ZZ/GF/QQ element or symbolic Integer, or None."""
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, ModularInteger):
        return None
    if isinstance(value, _core.Expr):
        if value.is_Integer:
            return int(value)
        return None
    if isinstance(value, _core.DomainElement):
        try:
            parent = value.parent()
        except Exception:
            return None
        if parent.is_ZZ:
            return int(value)
        return None
    try:
        import gmpy2
        if isinstance(value, type(gmpy2.mpz(0))):
            return int(value)
    except ImportError:
        pass
    return None


class ModularInteger:
    """An element of Z/nZ (SymPy's ``ModularInteger``)."""
    __slots__ = ('val',)
    mod = dom = sym = _parent = None

    def __init__(self, val):
        if isinstance(val, self.__class__):
            self.val = val.val % self.mod
        else:
            self.val = self._parent._convert_int(val) % self.mod

    def parent(self):
        return self._parent

    def modulus(self):
        return self.mod

    def __hash__(self):
        return hash((self.val, self.mod))

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.val)

    def __str__(self):
        return "%s mod %s" % (self.val, self.mod)

    def __int__(self):
        return int(self.val)

    def __index__(self):
        return int(self.val)

    def __getstate__(self):
        return {'val': self.val}

    def __setstate__(self, state):
        self.val = state['val']

    def __reduce__(self):
        return (_element, (self.mod, self.sym, self.val))

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.val)

    @classmethod
    def _get_val(cls, other):
        if isinstance(other, cls):
            return other.val
        return _int_of(other)

    def __add__(self, other):
        val = self._get_val(other)
        return NotImplemented if val is None else self.__class__(self.val + val)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        val = self._get_val(other)
        return NotImplemented if val is None else self.__class__(self.val - val)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        val = self._get_val(other)
        return NotImplemented if val is None else self.__class__(self.val * val)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        val = self._get_val(other)
        if val is None:
            return NotImplemented
        return self.__class__(self.val * _invert(val, self.mod))

    def __rtruediv__(self, other):
        return self.invert().__mul__(other)

    def __mod__(self, other):
        val = self._get_val(other)
        return NotImplemented if val is None else self.__class__(self.val % val)

    def __rmod__(self, other):
        val = self._get_val(other)
        return NotImplemented if val is None else self.__class__(val % self.val)

    def __pow__(self, exp):
        if not exp:
            return self.__class__(1)
        if exp < 0:
            val, exp = self.invert().val, -exp
        else:
            val = self.val
        return self.__class__(pow(val, int(exp), self.mod))

    def _compare(self, other, op):
        val = self._get_val(other)
        if val is None:
            return NotImplemented
        return op(self.val, val % self.mod)

    def _compare_deprecated(self, other, op):
        val = self._get_val(other)
        if val is None:
            return NotImplemented
        from .utilities.exceptions import sympy_deprecation_warning
        sympy_deprecation_warning(
            """Ordered comparisons with modular integers are deprecated.

            Use e.g. int(a) < int(b) instead of a < b.
            """,
            deprecated_since_version="1.13",
            active_deprecations_target="modularinteger-compare",
            stacklevel=4,
        )
        return op(self.val, val % self.mod)

    def __eq__(self, other):
        return self._compare(other, operator.eq)

    def __ne__(self, other):
        return self._compare(other, operator.ne)

    def __lt__(self, other):
        return self._compare_deprecated(other, operator.lt)

    def __le__(self, other):
        return self._compare_deprecated(other, operator.le)

    def __gt__(self, other):
        return self._compare_deprecated(other, operator.gt)

    def __ge__(self, other):
        return self._compare_deprecated(other, operator.ge)

    def __bool__(self):
        return bool(self.val)

    def invert(self):
        return self.__class__(_invert(self.val, self.mod))


_classes = {}
_domains = {}


def _element(mod, sym, val):
    return ModularRing(mod, sym)(val)


class ModularRing:
    """``GF(n)`` outside the native symmetric prime field range."""
    rep = 'FF'
    alias = 'FF'
    is_FiniteField = is_FF = True
    is_Numerical = True
    is_Exact = True
    is_Simple = True
    is_Ring = True
    is_Composite = False
    is_PID = False
    has_assoc_Ring = False
    has_assoc_Field = True

    def __new__(cls, mod, symmetric=True):
        key = (mod, bool(symmetric))
        domain = _domains.get(key)
        if domain is None:
            domain = object.__new__(cls)
            domain.mod = mod
            domain.sym = bool(symmetric)
            domain.dom = _core.ZZ
            name = ("SymmetricModularIntegerMod%s" if symmetric else "ModularIntegerMod%s") % mod
            domain.dtype = type(name, (ModularInteger,), {
                '__slots__': (), 'mod': mod, 'dom': _core.ZZ, 'sym': bool(symmetric), '_parent': domain})
            domain.dtype.__module__ = 'sympy.polys.domains.modularinteger'
            domain.zero = domain.dtype(0)
            domain.one = domain.dtype(1)
            _domains[key] = domain
        return domain

    def __getattr__(self, name):
        # The flags of the other domains.
        if name.startswith('is_') and name[3:4].isupper():
            return False
        raise AttributeError(name)

    def __reduce__(self):
        return (ModularRing, (self.mod, self.sym))

    @property
    def tp(self):
        return self.dtype

    @property
    def is_Field(self):
        return _isprime(self.mod)

    is_PID = is_Field

    def __str__(self):
        return 'GF(%s)' % self.mod

    __repr__ = __str__

    def __hash__(self):
        return hash(('FiniteField', self.mod, 'ZZ'))

    def __eq__(self, other):
        if isinstance(other, ModularRing):
            return self.mod == other.mod
        if isinstance(other, _core.Domain) and other.is_FiniteField:
            return other.characteristic() == self.mod
        return False

    def __ne__(self, other):
        return not self == other

    def __call__(self, *args):
        return self.dtype(*args)

    def _convert_int(self, value):
        integer = _int_of(value)
        if integer is None:
            raise CoercionFailed("expected an integer, got %s" % (value,))
        return integer

    def characteristic(self):
        return self.mod

    def get_field(self):
        return self

    def get_ring(self):
        raise DomainError('there is no ring associated with %s' % self)

    def to_int(self, a):
        aval = int(a)
        if self.sym and aval > self.mod // 2:
            aval -= self.mod
        return aval

    def to_sympy(self, a):
        return _core.Integer(self.to_int(a))

    def from_sympy(self, a):
        from . import sympify
        a = sympify(a)
        if a.is_Integer:
            return self.dtype(int(a))
        if getattr(a, 'is_Float', False) and a == int(a):
            return self.dtype(int(a))
        raise CoercionFailed("expected an integer, got %s" % a)

    def of_type(self, element):
        return isinstance(element, self.dtype)

    def __contains__(self, a):
        try:
            self.convert(a)
        except CoercionFailed:
            return False
        return True

    def convert(self, element, base=None):
        if base is not None:
            return self.convert_from(element, base)
        if self.of_type(element):
            return element
        if isinstance(element, ModularInteger):
            return self.dtype(int(element))
        integer = _int_of(element)
        if integer is not None:
            return self.dtype(integer)
        if isinstance(element, _core.DomainElement):
            parent = element.parent()
            if parent.is_FiniteField:
                return self.dtype(int(element))
            return self.convert_from(element, parent)
        if isinstance(element, _core.Expr):
            return self.from_sympy(element)
        raise CoercionFailed("Cannot convert %s of type %s to %s" % (element, type(element), self))

    def convert_from(self, element, base):
        if isinstance(element, ModularInteger) or getattr(base, 'is_FiniteField', False) or getattr(base, 'is_ZZ', False):
            return self.dtype(int(element))
        if getattr(base, 'is_QQ', False):
            q = base.denom(element) if hasattr(base, 'denom') else element.denominator
            if q == 1:
                return self.dtype(int(base.numer(element) if hasattr(base, 'numer') else element.numerator))
        raise CoercionFailed("Cannot convert %s of type %s from %s to %s" % (element, type(element), base, self))

    # Field arithmetic of the reference (``Field`` and ``Domain`` methods).
    def is_zero(self, a):
        return not a

    def is_one(self, a):
        return a == self.one

    def is_unit(self, a):
        return bool(a)

    def is_positive(self, a):
        return bool(a)

    def is_nonnegative(self, a):
        return True

    def is_negative(self, a):
        return False

    def is_nonpositive(self, a):
        return not a

    def abs(self, a):
        return abs(a) if not isinstance(a, ModularInteger) else a

    def neg(self, a):
        return -a

    def pos(self, a):
        return +a

    def add(self, a, b):
        return a + b

    def sub(self, a, b):
        return a - b

    def mul(self, a, b):
        return a * b

    def pow(self, a, b):
        return a ** b

    def exquo(self, a, b):
        return a / b

    def quo(self, a, b):
        return a / b

    def rem(self, a, b):
        return self.zero

    def div(self, a, b):
        return a / b, self.zero

    def revert(self, a):
        from .polys.polyerrors import NotReversible
        if a:
            return 1/a
        raise NotReversible('zero is not reversible')

    def gcd(self, a, b):
        return self.one

    def lcm(self, a, b):
        return a*b

    def numer(self, a):
        return a

    def denom(self, a):
        return self.one

    def is_square(self, a):
        raise NotImplementedError('quadratic residues modulo %s are not implemented' % self.mod)

    def exsqrt(self, a):
        raise NotImplementedError('square roots modulo %s are not implemented' % self.mod)

    def sqrt(self, a):
        raise NotImplementedError

    def poly_ring(self, *gens, **kwargs):
        raise NotImplementedError('polynomials over %s are not implemented' % self)

    old_poly_ring = poly_ring
    frac_field = old_frac_field = poly_ring


def GF(mod, symmetric=True):
    """``GF(n)``: the native prime field, or the reference ring Z/nZ."""
    if isinstance(mod, _core.Expr) and mod.is_Integer:
        mod = int(mod)
    elif not isinstance(mod, int):
        integer = _int_of(mod)
        if integer is None:
            raise ValueError("modulus must be a positive integer, got %s" % (mod,))
        mod = integer
    if mod <= 0:
        raise ValueError('modulus must be a positive integer, got %s' % mod)
    if symmetric and 2 <= mod < 1 << 32 and _isprime(mod):
        return _core.GF(mod)
    return ModularRing(mod, symmetric)
