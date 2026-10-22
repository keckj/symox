"""Polynomials over the ring Z/nZ of ``GF(n)`` outside the native path.

The native ``Poly`` covers symmetric prime moduli below ``2**32``. Composite
moduli, larger primes and ``symmetric=False`` give a ``ModularPoly``: a dense
reference representation over ``_modular_ring.ModularRing`` that runs the
ported reference dense routines (``polys/_dense``). The ``Poly`` methods
below follow the reference ``polytools.Poly`` code paths, so a zero divisor
raises ``NotInvertible`` where the reference does, and the reference ring
logic (``to_field``/``to_ring``) gives the same ``DomainError`` cases.
"""
from . import _core
from ._modular_ring import ModularRing, ModularInteger
from .polys._dense import densebasic as _db, densearith as _da, densetools as _dt
from .polys._dense import euclidtools as _de, sqfreetools as _ds, factortools as _df
from .polys.polyerrors import (DomainError, UnificationFailed, ExactQuotientFailed,
                               PolynomialError, CoercionFailed, OptionError)


def needs_ring(modulus, symmetric=True):
    """True when ``modulus`` is outside the native ``GF(p)`` polynomials."""
    if type(modulus) is int and symmetric is True:
        known = _NATIVE.get(modulus)
        if known is not None:
            return not known
        result = _needs_ring(modulus, symmetric)
        if len(_NATIVE) < 4096:
            _NATIVE[modulus] = not result
        return result
    return _needs_ring(modulus, symmetric)


_NATIVE = {}


def _needs_ring(modulus, symmetric):
    try:
        modulus = int(modulus)
    except (TypeError, ValueError):
        return False
    if modulus < 1:
        raise OptionError("'modulus' must a positive integer, got %s" % modulus)
    if symmetric is False:
        return True
    from .ntheory import isprime
    return not (2 <= modulus < 1 << 32 and isprime(modulus))


def ring_domain(domain=None, modulus=None, symmetric=True):
    """The ``ModularRing`` selected by the options, or None (native path)."""
    if isinstance(domain, ModularRing):
        return domain
    if modulus is not None and needs_ring(modulus, symmetric):
        return ModularRing(int(modulus), symmetric is not False)
    if isinstance(domain, str):
        import re
        match = re.match(r'^\s*(?:GF|FF)\s*\(\s*(\d+)\s*\)\s*$', domain)
        if match and needs_ring(int(match.group(1))):
            return ModularRing(int(match.group(1)))
    return None


def _sort_key(gens):
    return gens


class ModularPoly:
    """A ``Poly`` over ``GF(n)`` represented by the reference dense lists."""
    is_Poly = True
    is_commutative = True
    is_Matrix = False
    is_number = False
    _op_priority = 10.001

    @property
    def __class__(self):
        return _core.Poly

    # -- construction ---------------------------------------------------
    @classmethod
    def new(cls, rep, dom, gens):
        obj = object.__new__(cls)
        obj._rep = rep
        obj._dom = dom
        obj._gens = tuple(gens)
        return obj

    @classmethod
    def from_expr(cls, expr, gens, dom):
        from . import Poly
        native = expr if isinstance(expr, _core.Poly) else Poly(expr, *gens)
        gens = tuple(native.gens)
        terms = {}
        for monom, coeff in native.terms():
            terms[tuple(monom)] = dom.convert(coeff)
        return cls.new(_db.dmp_from_dict(terms, len(gens) - 1, dom), dom, gens)

    @classmethod
    def from_dict(cls, data, gens, dom):
        gens = tuple(gens)
        terms = {}
        for monom, coeff in data.items():
            if not isinstance(monom, tuple):
                monom = (monom,)
            terms[monom] = dom.convert(coeff)
        return cls.new(_db.dmp_from_dict(terms, len(gens) - 1, dom), dom, gens)

    @classmethod
    def from_list(cls, data, gens, dom):
        gens = tuple(gens)
        rep = _db.dmp_strip([dom.convert(c) for c in data], 0) if len(gens) == 1 else None
        if rep is None:
            raise NotImplementedError('multivariate from_list over %s' % dom)
        return cls.new(rep, dom, gens)

    def per(self, rep, gens=None, remove=None):
        if gens is None:
            gens = self._gens
        if remove is not None:
            gens = gens[:remove] + gens[remove + 1:]
            if not gens:
                return self._dom.to_sympy(rep)
        return ModularPoly.new(rep, self._dom, gens)

    @property
    def _u(self):
        return len(self._gens) - 1

    # -- basic properties -----------------------------------------------
    @property
    def gens(self):
        return self._gens

    @property
    def gen(self):
        return self._gens[0]

    @property
    def domain(self):
        return self._dom

    def get_domain(self):
        return self._dom

    def get_modulus(self):
        from . import Integer
        return Integer(self._dom.mod)

    @property
    def free_symbols(self):
        symbols = set()
        for gen in self._gens:
            symbols |= gen.free_symbols
        return symbols

    @property
    def free_symbols_in_domain(self):
        return set()

    @property
    def expr(self):
        return self.as_expr()

    @property
    def args(self):
        return (self.as_expr(),) + self._gens

    def as_dict(self, native=False):
        if native:
            return _db.dmp_to_dict(self._rep, self._u, self._dom)
        return {m: self._dom.to_sympy(c) for m, c in _db.dmp_to_dict(self._rep, self._u, self._dom).items()}

    def as_expr(self, *gens):
        from . import Add, Mul, Pow, Integer
        if gens:
            if len(gens) == 1 and isinstance(gens[0], dict):
                mapping = gens[0]
                gens = tuple(mapping.get(g, g) for g in self._gens)
            if len(gens) != len(self._gens):
                raise ValueError('expected %s generators, got %s' % (len(self._gens), len(gens)))
        else:
            gens = self._gens
        terms = []
        for monom, coeff in self.terms():
            factors = [coeff] + [Pow(g, e) if e != 1 else g for g, e in zip(gens, monom) if e]
            terms.append(Mul(*factors))
        return Add(*terms) if terms else Integer(0)

    def as_poly(self, *gens, **args):
        return self

    def terms(self, order=None):
        terms = _db.dmp_to_dict(self._rep, self._u, self._dom)
        if not terms:
            return [((0,) * len(self._gens), self._dom.to_sympy(self._dom.zero))]
        return [(m, self._dom.to_sympy(terms[m])) for m in sorted(terms, reverse=True)]

    def coeffs(self, order=None):
        return [c for _, c in self.terms()]

    def monoms(self, order=None):
        return [m for m, _ in self.terms()]

    def all_coeffs(self):
        if self._u:
            raise PolynomialError('multivariate polynomials not supported')
        return [self._dom.to_sympy(c) for c in self._rep] or [self._dom.to_sympy(self._dom.zero)]

    def all_terms(self):
        if self._u:
            raise PolynomialError('multivariate polynomials not supported')
        n = len(self._rep) - 1
        return [((n - i,), self._dom.to_sympy(c)) for i, c in enumerate(self._rep)]

    def as_list(self, native=False):
        return _db.dmp_to_tuple(self._rep, self._u) and self._to_list(self._rep, self._u)

    def _to_list(self, rep, u):
        if not u:
            return [self._dom.to_sympy(c) for c in rep]
        return [self._to_list(c, u - 1) for c in rep]

    def length(self):
        return len(_db.dmp_to_dict(self._rep, self._u, self._dom))

    def _gen_index(self, gen):
        if gen is None:
            return 0
        if isinstance(gen, int):
            if -len(self._gens) <= gen < len(self._gens):
                return gen % len(self._gens)
            raise ValueError('-%s <= gen < %s expected, got %s' % (len(self._gens), len(self._gens), gen))
        from . import sympify
        gen = sympify(gen)
        if gen in self._gens:
            return self._gens.index(gen)
        raise ValueError('a valid generator expected, got %s' % gen)

    def degree(self, gen=0):
        j = self._gen_index(gen)
        d = _db.dmp_degree_in(self._rep, j, self._u)
        if d < 0:
            from . import S
            return S.NegativeInfinity
        return d

    def degree_list(self):
        return tuple(self.degree(j) for j in range(len(self._gens)))

    def total_degree(self):
        terms = _db.dmp_to_dict(self._rep, self._u, self._dom)
        from . import S
        return max((sum(m) for m in terms), default=S.NegativeInfinity)

    def LC(self, order=None):
        return self._dom.to_sympy(_db.dmp_ground_LC(self._rep, self._u, self._dom))

    def TC(self):
        return self._dom.to_sympy(_db.dmp_ground_TC(self._rep, self._u, self._dom))

    def nth(self, *N):
        return self._dom.to_sympy(_db.dmp_ground_nth(self._rep, N, self._u, self._dom))

    def coeff_monomial(self, monom):
        from . import Poly, sympify
        monom = Poly(sympify(monom), *self._gens).monoms()
        if len(monom) != 1:
            raise ValueError('a monomial expected')
        return self.nth(*monom[0])

    @property
    def is_zero(self):
        return _db.dmp_zero_p(self._rep, self._u)

    @property
    def is_one(self):
        return _db.dmp_one_p(self._rep, self._u, self._dom)

    @property
    def is_ground(self):
        return _db.dmp_ground_p(self._rep, None, self._u)

    @property
    def is_monic(self):
        return self._dom.is_one(_db.dmp_ground_LC(self._rep, self._u, self._dom))

    @property
    def is_primitive(self):
        return self._dom.is_one(_dt.dmp_ground_content(self._rep, self._u, self._dom))

    @property
    def is_linear(self):
        return all(sum(m) <= 1 for m in _db.dmp_to_dict(self._rep, self._u, self._dom))

    @property
    def is_quadratic(self):
        return all(sum(m) <= 2 for m in _db.dmp_to_dict(self._rep, self._u, self._dom))

    @property
    def is_monomial(self):
        return len(_db.dmp_to_dict(self._rep, self._u, self._dom)) <= 1

    @property
    def is_homogeneous(self):
        degrees = {sum(m) for m in _db.dmp_to_dict(self._rep, self._u, self._dom)}
        return len(degrees) <= 1

    @property
    def is_univariate(self):
        return len(self._gens) == 1

    @property
    def is_multivariate(self):
        return len(self._gens) != 1

    @property
    def is_sqf(self):
        return _ds.dmp_sqf_p(self._rep, self._u, self._dom)

    @property
    def is_irreducible(self):
        return _df.dmp_irreducible_p(self._rep, self._u, self._dom)

    @property
    def is_cyclotomic(self):
        if self._u:
            return False
        return _df.dup_cyclotomic_p(self._rep, self._dom)

    def has_only_gens(self, *gens):
        indices = {self._gen_index(g) for g in gens}
        return all(not e or i in indices for m in self.monoms() for i, e in enumerate(m))

    # -- comparison and printing ----------------------------------------
    def __repr__(self):
        from . import Poly
        native = repr(Poly(self.as_expr(), *self._gens, domain='ZZ'))
        return native[:native.rindex(', domain=')] + ', modulus=%s)' % self._dom.mod

    __str__ = __repr__

    def _sympystr(self, printer):
        return repr(self)

    def _latex(self, printer):
        return r'\operatorname{Poly}{\left( %s, %s, domain=\mathbb{F}_{%s} \right)}' % (
            printer._print(self.as_expr()), ', '.join(printer._print(g) for g in self._gens), self._dom.mod)

    def __hash__(self):
        return hash(('Poly', _db.dmp_to_tuple(self._rep, self._u), self._gens, self._dom.mod))

    def __eq__(self, other):
        if not isinstance(other, ModularPoly):
            if isinstance(other, _core.Poly):
                return False
            try:
                other = ModularPoly.from_expr(other, self._gens, self._dom)
            except Exception:
                return False
        return (self._gens == other._gens and self._dom == other._dom
                and _db.dmp_to_tuple(self._rep, self._u) == _db.dmp_to_tuple(other._rep, other._u))

    def __ne__(self, other):
        return not self == other

    def eq(self, other, strict=False):
        return self == other

    def ne(self, other, strict=False):
        return not self == other

    def __bool__(self):
        return not self.is_zero

    def __getnewargs__(self):
        return ()

    def __reduce__(self):
        return (_rebuild, (self.as_expr(), self._gens, self._dom.mod, self._dom.sym))

    # -- unification ----------------------------------------------------
    def _unify(self, other):
        from . import sympify
        if not isinstance(other, ModularPoly):
            if isinstance(other, _core.Poly):
                raise UnificationFailed("Cannot unify %s with %s" % (self._dom, other.domain))
            other = sympify(other)
            if not other.is_Integer:
                # The reference converts the operand with the default domain.
                from . import Poly
                raise UnificationFailed("Cannot unify %s with %s" % (self._dom, Poly(other, *self._gens).domain))
            other = ModularPoly.from_expr(other, self._gens, self._dom)
        if other._dom != self._dom:
            raise UnificationFailed("Cannot unify %s with %s" % (self._dom, other._dom))
        if other._gens != self._gens:
            raise NotImplementedError('polynomials over %s with different generators' % self._dom)
        return self._dom, other

    def unify(self, other):
        _, other = self._unify(other)
        return self, other

    def _ground(self, value):
        return self._dom.convert(value)

    # -- arithmetic -----------------------------------------------------
    def _binary(self, other, func):
        try:
            _, other = self._unify(other)
        except (CoercionFailed, PolynomialError):
            return NotImplemented
        return self.per(func(self._rep, other._rep, self._u, self._dom))

    def add(self, other):
        return self._binary(other, _da.dmp_add)

    def sub(self, other):
        return self._binary(other, _da.dmp_sub)

    def mul(self, other):
        return self._binary(other, _da.dmp_mul)

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        return self.add(other)

    def __sub__(self, other):
        return self.sub(other)

    def __rsub__(self, other):
        return (-self).add(other)

    def __mul__(self, other):
        return self.mul(other)

    def __rmul__(self, other):
        return self.mul(other)

    def __neg__(self):
        return self.per(_da.dmp_neg(self._rep, self._u, self._dom))

    def neg(self):
        return -self

    def __pos__(self):
        return self

    def __pow__(self, n):
        from . import Integer
        if isinstance(n, (int, Integer)) and int(n) >= 0:
            return self.pow(int(n))
        return NotImplemented

    def pow(self, n):
        return self.per(_da.dmp_pow(self._rep, int(n), self._u, self._dom))

    def sqr(self):
        return self.per(_da.dmp_sqr(self._rep, self._u, self._dom))

    def add_ground(self, c):
        return self.per(_da.dmp_add_ground(self._rep, self._ground(c), self._u, self._dom))

    def sub_ground(self, c):
        return self.per(_da.dmp_sub_ground(self._rep, self._ground(c), self._u, self._dom))

    def mul_ground(self, c):
        return self.per(_da.dmp_mul_ground(self._rep, self._ground(c), self._u, self._dom))

    def quo_ground(self, c):
        return self.per(_da.dmp_quo_ground(self._rep, self._ground(c), self._u, self._dom))

    def exquo_ground(self, c):
        return self.per(_da.dmp_exquo_ground(self._rep, self._ground(c), self._u, self._dom))

    def abs(self):
        return self.per(_da.dmp_abs(self._rep, self._u, self._dom))

    def max_norm(self):
        return self._dom.to_sympy(_da.dmp_max_norm(self._rep, self._u, self._dom))

    def l1_norm(self):
        return self._dom.to_sympy(_da.dmp_l1_norm(self._rep, self._u, self._dom))

    def l2_norm_squared(self):
        return self._dom.to_sympy(_da.dmp_l2_norm_squared(self._rep, self._u, self._dom))

    # -- division (the reference ``auto`` ring logic) ---------------------
    def _to_ring(self, result):
        if not self._dom.is_Field:
            self._dom.get_ring()
        return result

    def div(self, other, auto=True):
        _, other = self._unify(other)
        q, r = _da.dmp_div(self._rep, other._rep, self._u, self._dom)
        if auto:
            self._to_ring(None)
        return self.per(q), self.per(r)

    def rem(self, other, auto=True):
        _, other = self._unify(other)
        r = _da.dmp_rem(self._rep, other._rep, self._u, self._dom)
        if auto:
            self._to_ring(None)
        return self.per(r)

    def quo(self, other, auto=True):
        _, other = self._unify(other)
        q = _da.dmp_quo(self._rep, other._rep, self._u, self._dom)
        if auto:
            self._to_ring(None)
        return self.per(q)

    def exquo(self, other, auto=True):
        _, other = self._unify(other)
        q, r = _da.dmp_div(self._rep, other._rep, self._u, self._dom)
        if not _db.dmp_zero_p(r, self._u):
            raise ExactQuotientFailed('%s does not divide %s' % (other.as_expr(), self.as_expr()))
        if auto:
            self._to_ring(None)
        return self.per(q)

    def __truediv__(self, other):
        return self.as_expr() / (other.as_expr() if isinstance(other, ModularPoly) else other)

    def __floordiv__(self, other):
        return self.quo(other)

    def __mod__(self, other):
        return self.rem(other)

    def __divmod__(self, other):
        return self.div(other)

    def prem(self, other):
        _, other = self._unify(other)
        return self.per(_de.dmp_prem(self._rep, other._rep, self._u, self._dom))

    def pquo(self, other):
        _, other = self._unify(other)
        return self.per(_de.dmp_pquo(self._rep, other._rep, self._u, self._dom))

    def pexquo(self, other):
        _, other = self._unify(other)
        try:
            return self.per(_de.dmp_pexquo(self._rep, other._rep, self._u, self._dom))
        except ExactQuotientFailed as exc:
            raise ExactQuotientFailed('%s does not divide %s' % (other.as_expr(), self.as_expr()))

    def pdiv(self, other):
        _, other = self._unify(other)
        q, r = _da.dmp_pdiv(self._rep, other._rep, self._u, self._dom)
        return self.per(q), self.per(r)

    # -- gcd and friends --------------------------------------------------
    def gcd(self, other):
        _, other = self._unify(other)
        return self.per(_de.dmp_gcd(self._rep, other._rep, self._u, self._dom))

    def lcm(self, other):
        _, other = self._unify(other)
        return self.per(_de.dmp_lcm(self._rep, other._rep, self._u, self._dom))

    def cofactors(self, other):
        _, other = self._unify(other)
        h, cff, cfg = _de.dmp_inner_gcd(self._rep, other._rep, self._u, self._dom)
        return self.per(h), self.per(cff), self.per(cfg)

    def half_gcdex(self, other, auto=True):
        _, other = self._unify(other)
        if self._u:
            raise ValueError('univariate polynomial expected')
        s, h = _de.dup_half_gcdex(self._rep, other._rep, self._dom)
        return self.per(s), self.per(h)

    def gcdex(self, other, auto=True):
        _, other = self._unify(other)
        if self._u:
            raise ValueError('univariate polynomial expected')
        if not self._dom.is_Field:
            raise DomainError('ground domain must be a field')
        s, t, h = _de.dup_gcdex(self._rep, other._rep, self._dom)
        return self.per(s), self.per(t), self.per(h)

    def invert(self, other, auto=True):
        _, other = self._unify(other)
        if self._u:
            raise ValueError('univariate polynomial expected')
        return self.per(_de.dup_invert(self._rep, other._rep, self._dom))

    def subresultants(self, other):
        _, other = self._unify(other)
        return [self.per(r) for r in _de.dmp_subresultants(self._rep, other._rep, self._u, self._dom)]

    def resultant(self, other, includePRS=False):
        _, other = self._unify(other)
        if includePRS:
            res, R = _de.dmp_resultant(self._rep, other._rep, self._u, self._dom, includePRS=True)
            return self.per(res, remove=0), [self.per(r) for r in R]
        return self.per(_de.dmp_resultant(self._rep, other._rep, self._u, self._dom), remove=0)

    def discriminant(self):
        return self.per(_de.dmp_discriminant(self._rep, self._u, self._dom), remove=0)

    def cancel(self, other, include=True):
        _, other = self._unify(other)
        dom = self._dom
        F, G = self._rep, other._rep
        if dom.is_Field and dom.has_assoc_Ring:
            dom = dom.get_ring()
        cp, cq, p, q = _de.dmp_cancel(F, G, self._u, self._dom, include=False)
        cp, cq = self._dom.to_sympy(cp), self._dom.to_sympy(cq)
        p, q = self.per(p), self.per(q)
        if include:
            return p.mul_ground(cp / cq if cq != 1 else cp) if False else (cp / cq, p, q)
        return cp, cq, p, q

    def content(self):
        return self._dom.to_sympy(_dt.dmp_ground_content(self._rep, self._u, self._dom))

    def primitive(self):
        cont, result = _dt.dmp_ground_primitive(self._rep, self._u, self._dom)
        return self._dom.to_sympy(cont), self.per(result)

    def monic(self, auto=True):
        return self.per(_dt.dmp_ground_monic(self._rep, self._u, self._dom))

    def terms_gcd(self):
        J, result = _db.dmp_terms_gcd(self._rep, self._u, self._dom)
        return J, self.per(result)

    def trunc(self, p):
        return self.per(_dt.dmp_ground_trunc(self._rep, self._ground(p), self._u, self._dom))

    def clear_denoms(self, convert=False):
        coeff, result = _dt.dmp_clear_denoms(self._rep, self._u, self._dom)
        return self._dom.to_sympy(coeff), self.per(result)

    # -- calculus and composition ------------------------------------------
    def diff(self, *specs, **kwargs):
        if not specs:
            specs = (self._gens[0],)
        rep = self._rep
        for spec in specs:
            if isinstance(spec, tuple):
                gen, m = spec
            else:
                gen, m = spec, 1
            rep = _dt.dmp_diff_in(rep, int(m), self._gen_index(gen), self._u, self._dom)
        return self.per(rep)

    _eval_derivative = diff

    def integrate(self, *specs, **args):
        if not specs:
            specs = (self._gens[0],)
        rep = self._rep
        for spec in specs:
            if isinstance(spec, tuple):
                gen, m = spec
            else:
                gen, m = spec, 1
            rep = _dt.dmp_integrate_in(rep, int(m), self._gen_index(gen), self._u, self._dom)
        return self.per(rep)

    def eval(self, x, a=None, auto=True):
        if a is None:
            if isinstance(x, dict):
                result = self
                for gen, value in x.items():
                    result = result.eval(gen, value) if isinstance(result, ModularPoly) else result
                return result
            if isinstance(x, (tuple, list)):
                result = self
                for value in x:
                    result = result.eval(result.gens[0], value)
                return result
            j, a = 0, x
        else:
            j = self._gen_index(x)
        a = self._eval_point(a)
        return self.per(_dt.dmp_eval_in(self._rep, a, j, self._u, self._dom), remove=j)

    def _eval_point(self, a):
        from . import sympify, Poly
        if isinstance(a, ModularInteger) and a.parent() == self._dom:
            return a
        a = sympify(a)
        if a.is_Integer:
            return self._dom.convert(a)
        if a.is_Rational:
            raise UnificationFailed("Cannot unify %s with QQ" % self._dom)
        gens = tuple(sorted(a.free_symbols, key=str))
        if gens:
            ground = Poly(a, *gens).domain
            raise UnificationFailed("Cannot unify %s with %s" % (self._dom, ground.unify(_core.ZZ) if False else
                                                               'ZZ[%s]' % ','.join(map(str, gens))))
        raise UnificationFailed("Cannot unify %s with EX" % self._dom)

    def __call__(self, *values):
        return self.eval(values) if len(values) > 1 else self.eval(values[0])

    def subs(self, *args, **kwargs):
        from . import sympify
        if len(args) == 2:
            old, new = sympify(args[0]), sympify(args[1])
            if old in self._gens and new.is_number:
                return self.eval(old, new)
        return self.as_expr().subs(*args, **kwargs)

    def _no_real_roots(self, *args, **kwargs):
        raise DomainError('isolation of real roots not supported over %s' % self._dom)

    intervals = refine_root = _no_real_roots

    def sturm(self, auto=True):
        raise DomainError('Cannot compute Sturm sequence over %s' % self._dom)

    def count_roots(self, inf=None, sup=None):
        raise DomainError('Cannot compute Sturm sequence over %s' % self._dom)

    def norm(self):
        raise DomainError('ground domain must be algebraic')

    sqf_norm = norm

    def replace(self, x, y=None, **_ignore):
        if y is None:
            x, y = self.gen, x
        from . import sympify
        x, y = sympify(x), sympify(y)
        if x in self._gens and y not in self._gens:
            gens = tuple(y if g == x else g for g in self._gens)
            return ModularPoly.new(self._rep, self._dom, gens)
        raise PolynomialError("Cannot replace %s with %s in %s" % (x, y, self))

    def reorder(self, *gens, **args):
        if not gens or tuple(gens) == self._gens:
            return self
        raise NotImplementedError('reordering polynomials over %s' % self._dom)

    def ltrim(self, gen):
        if self._gen_index(gen) == 0:
            return self
        raise NotImplementedError('ltrim over %s' % self._dom)

    def compose(self, other):
        _, other = self._unify(other)
        return self.per(_dt.dmp_compose(self._rep, other._rep, self._u, self._dom))

    def decompose(self):
        if self._u:
            raise PolynomialError('univariate polynomial expected')
        return [self.per(g) for g in _dt.dup_decompose(self._rep, self._dom)]

    def shift(self, a):
        if self._u:
            raise PolynomialError('univariate polynomial expected')
        return self.per(_dt.dup_shift(self._rep, self._ground(a), self._dom))

    def revert(self, n):
        if self._u:
            raise PolynomialError('univariate polynomial expected')
        return self.per(_dt.dup_revert(self._rep, int(n), self._dom))

    def slice(self, x, m, n=None):
        if n is None:
            j, m, n = 0, x, m
        else:
            j = self._gen_index(x)
        return self.per(_db.dmp_slice_in(self._rep, int(m), int(n), j, self._u, self._dom))

    def homogenize(self, s):
        from . import Poly
        total = self.total_degree()
        terms = {}
        for monom, coeff in _db.dmp_to_dict(self._rep, self._u, self._dom).items():
            terms[monom + (total - sum(monom),)] = coeff
        if s in self._gens:
            raise NotImplementedError('homogenize with an existing generator over %s' % self._dom)
        return ModularPoly.from_dict(terms, self._gens + (s,), self._dom)

    # -- square-free and factorization -------------------------------------
    def sqf_part(self):
        return self.per(_ds.dmp_sqf_part(self._rep, self._u, self._dom))

    def sqf_list(self, all=False):
        coeff, factors = _ds.dmp_sqf_list(self._rep, self._u, self._dom, all=all)
        return self._dom.to_sympy(coeff), [(self.per(g), k) for g, k in factors]

    def sqf_list_include(self, all=False):
        factors = _ds.dmp_sqf_list_include(self._rep, self._u, self._dom, all=all)
        return [(self.per(g), k) for g, k in factors]

    def factor_list(self):
        coeff, factors = _df.dmp_factor_list(self._rep, self._u, self._dom)
        return self._dom.to_sympy(coeff), [(self.per(g), k) for g, k in factors]

    def factor_list_include(self):
        factors = _df.dmp_factor_list_include(self._rep, self._u, self._dom)
        return [(self.per(g), k) for g, k in factors]

    def ground_roots(self):
        if self.is_multivariate:
            from .polys.polyerrors import MultivariatePolynomialError
            raise MultivariatePolynomialError("Cannot compute ground roots of %s" % self)
        roots = {}
        _, factors = self.factor_list()
        for factor, k in factors:
            if factor.degree() == 1:
                rep = factor._rep
                roots[self._dom.to_sympy(-rep[1] / rep[0])] = k
        return roots

    # -- domains -----------------------------------------------------------
    def set_domain(self, domain):
        from . import Poly
        if isinstance(domain, ModularRing):
            return ModularPoly.from_expr(self.as_expr(), self._gens, domain)
        return Poly(self.as_expr(), *self._gens, domain=domain)

    def set_modulus(self, modulus):
        from . import Poly
        dom = ring_domain(modulus=modulus)
        if dom is not None:
            return ModularPoly.from_expr(self.as_expr(), self._gens, dom)
        return Poly(self.as_expr(), *self._gens, modulus=modulus)

    def retract(self, field=None):
        from . import Poly
        return Poly(self.as_expr(), *self._gens)

    def to_field(self):
        return self

    def to_exact(self):
        return self

    def to_ring(self):
        if self._dom.is_Field:
            raise DomainError('there is no ring associated with %s' % self._dom)
        self._dom.get_ring()

    def inject(self, front=False):
        return self

    def exclude(self):
        return self

    def termwise(self, func, *gens, **args):
        from . import Poly
        terms = {}
        for monom, coeff in self.terms():
            result = func(monom, coeff) if True else None
            if isinstance(result, tuple):
                monom, result = result
            if result:
                terms[monom] = terms.get(monom, 0) + result
        return Poly(terms, *(gens or self._gens), **args)

    def __getattr__(self, name):
        if name.startswith('__') or name.startswith('_'):
            raise AttributeError(name)
        if hasattr(_core.Poly, name):
            raise NotImplementedError('Poly.%s over %s is not implemented' % (name, self._dom))
        raise AttributeError("'Poly' object has no attribute '%s'" % name)


def _rebuild(expr, gens, mod, sym):
    return ModularPoly.from_expr(expr, gens, ModularRing(mod, sym))


def construct(rep, gens, domain, args):
    """``Poly(rep, *gens, ...)`` for a ``ModularRing`` domain."""
    from . import sympify, Poly
    from .utilities.iterables import iterable
    if isinstance(rep, ModularPoly):
        if not gens or tuple(gens) == rep.gens:
            return rep.set_domain(domain) if rep.domain != domain else rep
        rep = rep.as_expr()
    if isinstance(rep, dict):
        if not gens:
            raise _gens_required()
        if len(gens) == 1 and gens[0] is not None and isinstance(gens[0], (tuple, list)):
            gens = tuple(gens[0])
        return ModularPoly.from_dict(rep, gens, domain)
    if not isinstance(rep, _core.Poly) and iterable(rep, exclude=str):
        if not gens:
            raise _gens_required()
        return ModularPoly.from_list(list(rep), gens, domain)
    if not isinstance(rep, _core.Poly):
        rep = sympify(rep)
        if getattr(rep, 'is_Equality', False):
            rep = rep.lhs - rep.rhs
        native = Poly(rep, *gens, **{k: v for k, v in args.items() if k in ('expand',)})
    else:
        native = rep if not gens else Poly(rep.as_expr(), *gens)
    return ModularPoly.from_expr(native, native.gens, domain)


def _gens_required():
    from .polys.polyerrors import GeneratorsNeeded
    return GeneratorsNeeded("generators are required")
