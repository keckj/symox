"""Sparse ring constructors using native exact domain arithmetic."""
import inspect
import operator
from .orderings import lex, monomial_key

__all__ = ['ring', 'xring', 'vring', 'sring']


def _symbols(value):
    from .. import symbols, sympify
    from .._core import Expr
    if isinstance(value, str):
        return tuple(symbols(value, seq=True)) if value else ()
    if isinstance(value, Expr):
        return (value,)
    return tuple(sympify(v) for v in value)


def _domain(value):
    from .. import _core, Poly, Dummy
    if isinstance(value, _core.Domain):
        return value
    return Poly(0, Dummy('domain'), domain=value).domain


class PolyRing:
    def __init__(self, symbols, domain, order=lex):
        self.symbols = _symbols(symbols)
        self.ngens = len(self.symbols)
        self.domain = _domain(domain)
        self.order = monomial_key(order)
        if not getattr(self.order,'is_global',False):
            raise NotImplementedError('polynomial arithmetic requires a global monomial order')
        if len(set(self.symbols)) != self.ngens:
            raise ValueError('duplicate generators')
        self._native = self.domain.poly_ring(*self.symbols, order=str(self.order)) if self.ngens else self.domain
        self.gens = tuple(self(s) for s in self.symbols)
        self.zero, self.one = self(0), self(1)
        self.zero_monom = (0,)*self.ngens
        self.dtype = self.__call__
        for symbol, gen in zip(self.symbols, self.gens):
            name = getattr(symbol, 'name', None)
            if name and not hasattr(self, name):
                setattr(self, name, gen)

    def __eq__(self, other):
        return isinstance(other, PolyRing) and (self.symbols,self.domain,self.order)==(other.symbols,other.domain,other.order)

    def __hash__(self):
        return hash((self.symbols,self.domain,self.order))

    def __repr__(self):
        return 'Polynomial ring in %s over %s with %s order' % (', '.join(map(str,self.symbols)),self.domain,self.order)

    def __call__(self, value=0):
        if isinstance(value, dict):
            return self.from_dict(value)
        if isinstance(value, PolyElement):
            if value.ring == self:
                return value
            value = value.as_expr()
        return PolyElement(self, self._native.convert(value))

    from_expr = __call__
    ground_new = __call__

    def from_dict(self, terms):
        from .. import Integer
        result = Integer(0)
        for powers, coeff in terms.items():
            powers = tuple(operator.index(p) for p in powers)
            if len(powers) != self.ngens or any(p < 0 for p in powers):
                raise ValueError('invalid monomial exponent tuple')
            term = self.domain.convert(coeff).as_expr()
            for x,n in zip(self.symbols,powers):
                term *= x**n
            result += term
        return self(result)

    def index(self, gen):
        if isinstance(gen, int):
            if not 0 <= gen < self.ngens:
                raise ValueError('generator index out of range')
            return gen
        value = gen.as_expr() if isinstance(gen, PolyElement) else gen
        try:
            return self.symbols.index(value)
        except ValueError:
            raise ValueError('expected a ring generator') from None

    def to_domain(self):
        return self._native

    def to_field(self):
        from .fields import FracField
        return FracField(self.symbols,self.domain,self.order)

    def is_element(self, value):
        return isinstance(value, PolyElement) and value.ring == self


class PolyElement:
    # Expressions on the left defer to the element, as the reference
    # ``Expr`` binary methods return ``NotImplemented`` for it.
    _op_priority = 10.1

    def __init__(self, ring, value):
        self.ring, self._value = ring, value

    def parent(self):
        return self.ring

    def as_expr(self, *symbols):
        value = self._value.as_expr()
        if symbols:
            if len(symbols) != self.ring.ngens:
                raise ValueError('wrong number of generators')
            value = value.xreplace(dict(zip(self.ring.symbols,symbols)))
        return value

    _symox_ = as_expr

    def to_sympy(self):
        return self.as_expr().to_sympy()

    def __repr__(self):
        return _poly_element_str(self)

    def __hash__(self):
        return hash((self.ring,self._value))

    def __eq__(self, other):
        try:
            return self._value == self.ring(other)._value
        except (TypeError,ValueError):
            return False

    def __bool__(self):
        return self._value != self.ring._native.zero

    def __add__(self, other):
        if not isinstance(other, (PolyElement, int)) and hasattr(other, 'field'):
            return NotImplemented
        try:
            return self.ring(self._value+self.ring(other)._value)
        except TypeError:
            return NotImplemented
    __radd__ = __add__

    def __neg__(self):
        return self.ring(-self._value)

    def __sub__(self, other):
        return self + (-self.ring(other))

    def __rsub__(self, other):
        return self.ring(other)-self

    def __mul__(self, other):
        if hasattr(other, 'field'):
            return NotImplemented
        try:
            return self.ring(self._value*self.ring(other)._value)
        except TypeError:
            return NotImplemented
    __rmul__ = __mul__

    def __pow__(self, n):
        n = operator.index(n)
        if n < 0:
            raise ValueError('polynomial powers must be nonnegative')
        return self.ring(self._value**n)

    def __truediv__(self, other):
        return self.exquo(other)

    def __rtruediv__(self, other):
        return self.ring(other).exquo(self)

    def __floordiv__(self, other):
        return self.div(other)[0]

    def __mod__(self, other):
        return self.div(other)[1]

    def __divmod__(self, other):
        return self.div(other)

    def _poly(self):
        from .. import Poly
        return Poly(self.as_expr(), *self.ring.symbols, domain=self.ring.domain)

    def to_dict(self):
        if not self.ring.ngens:
            return {(): self._value} if self else {}
        return {m:self.ring.domain.convert(c) for m,c in self._poly().as_dict().items()}

    def __len__(self):
        return len(self.to_dict())

    def __iter__(self):
        return iter(self.to_dict())

    def __getitem__(self, monom):
        return self.to_dict().get(monom,self.ring.domain.zero)

    def items(self):
        return self.to_dict().items()

    def terms(self, order=None):
        return sorted(self.items(), key=lambda item: monomial_key(order or self.ring.order)(item[0]), reverse=True)

    def monoms(self, order=None):
        return [m for m,c in self.terms(order)]

    def coeffs(self, order=None):
        return [c for m,c in self.terms(order)]

    def degree(self, x=None):
        from .. import oo
        if not self:
            return -oo
        i=self.ring.index(x) if x is not None else 0
        return max(m[i] for m in self)

    def degrees(self):
        return tuple(self.degree(i) for i in range(self.ring.ngens))

    def total_degree(self):
        return max(map(sum,self), default=0)

    @property
    def is_zero(self):
        return not bool(self)

    @property
    def is_one(self):
        return self == 1

    @property
    def LC(self):
        return self.coeffs()[0] if self else self.ring.domain.zero

    @property
    def LM(self):
        return self.monoms()[0] if self else self.ring.zero_monom

    def div(self, other):
        from .. import reduced
        many = isinstance(other,(list,tuple))
        divisors = list(other) if many else [other]
        divisors = [self.ring(g) for g in divisors]
        if any(not g for g in divisors):
            raise ZeroDivisionError('polynomial division by zero')
        if not self.ring.ngens:
            remainder=self
            quotients=[]
            for g in divisors:
                q,r=self.ring._native.div(remainder._value,g._value)
                quotients.append(self.ring(q));remainder=self.ring(r)
            return (quotients if many else quotients[0]),remainder
        q,r = reduced(self.as_expr(), [g.as_expr() for g in divisors], *self.ring.symbols,
                      domain=self.ring.domain, order=str(self.ring.order))
        q = [self.ring(v) for v in q]
        return (q if many else q[0]),self.ring(r)

    def rem(self, other):
        return self.div(other)[1]

    def quo(self, other):
        return self.div(other)[0]

    def exquo(self, other):
        return self.ring(self.ring._native.exquo(self._value,self.ring(other)._value))

    def gcd(self, other):
        return self.ring(self.ring._native.gcd(self._value,self.ring(other)._value))

    def lcm(self, other):
        return self.ring(self.ring._native.lcm(self._value,self.ring(other)._value))

    def diff(self, x):
        return self.ring(self.as_expr().diff(self.ring.symbols[self.ring.index(x)]))

    def evaluate(self, x, a=None):
        assignments = x if a is None else [(x,a)]
        substitutions = {self.ring.symbols[self.ring.index(g)]: v.as_expr() if isinstance(v,PolyElement) else v for g,v in assignments}
        remaining = tuple(g for g in self.ring.symbols if g not in substitutions)
        value = self.as_expr().subs(substitutions, simultaneous=True)
        return PolyRing(remaining,self.ring.domain,self.ring.order)(value) if remaining else self.ring.domain.convert(value)

    def drop(self, *gens):
        omitted={self.ring.index(g) for g in gens}
        if any(any(m[i] for i in omitted) for m in self):
            raise ValueError('cannot drop an active generator')
        return PolyRing(tuple(s for i,s in enumerate(self.ring.symbols) if i not in omitted),self.ring.domain,self.ring.order)(self.as_expr())

    def set_ring(self, ring):
        return ring(self.as_expr())


def ring(symbols, domain, order=lex):
    R=PolyRing(symbols,domain,order)
    return (R,)+R.gens


def xring(symbols, domain, order=lex):
    R=PolyRing(symbols,domain,order)
    return R,R.gens


def vring(symbols, domain, order=lex):
    R=PolyRing(symbols,domain,order)
    inspect.currentframe().f_back.f_globals.update({str(s):g for s,g in zip(R.symbols,R.gens)})
    return R


def _parenthesize(item, level, strict):
    """``StrPrinter.parenthesize`` over expressions and ring elements."""
    from ..printing.str import sstr
    text, prec = _element_text(item)
    if prec < level or (not strict and prec <= level):
        return "(%s)" % text
    return text


def _element_text(item):
    from ..printing.str import sstr
    from ..printing.precedence import precedence, PRECEDENCE
    if isinstance(item, PolyElement):
        return _poly_element_str(item), _poly_element_precedence(item)
    from .fields import FracElement
    if isinstance(item, _Coefficient):
        domain, coeff = item
        if getattr(domain, 'is_ZZ', False) or getattr(domain, 'is_QQ', False) \
                or getattr(domain, 'is_RR', False) or getattr(domain, 'is_CC', False):
            # Ground elements print natively (``1/2``, ``2.5``).
            return str(coeff), precedence(domain.to_sympy(coeff))
        return _element_text(domain.to_sympy(coeff))
    if isinstance(item, FracElement):
        prec = (_poly_element_precedence(item.numer) if item.denom == 1
                else PRECEDENCE["Mul"])
        return repr(item), prec
    return sstr(item), precedence(item)


class _Coefficient(tuple):
    """A ``(domain, element)`` pair printed like the reference raw coefficient."""


def _poly_element_precedence(item):
    """The reference ``precedence_PolyElement``."""
    from ..printing.precedence import precedence, PRECEDENCE
    terms = item.terms()
    zm = item.ring.zero_monom
    if len(terms) == 1 and terms[0][1] == item.ring.domain.one and sum(terms[0][0]) == 1:
        return PRECEDENCE["Atom"]
    if not terms or (len(terms) == 1 and terms[0][0] == zm):
        coeff = terms[0][1] if terms else item.ring.domain.zero
        return precedence(item.ring.domain.to_sympy(coeff))
    if len(terms) == 1:
        return PRECEDENCE["Mul"]
    return PRECEDENCE["Add"]


def _poly_element_str(self):
    """The reference ``PolyElement.str`` with ``StrPrinter``: terms in the
    ring order, non-atomic generators parenthesized (``x*(exp(x))``)."""
    from ..printing.precedence import PRECEDENCE
    from ..printing.str import sstr
    ring = self.ring
    domain = ring.domain
    terms = self.terms()
    if not terms:
        return _element_text(_Coefficient((domain, domain.zero)))[0]
    prec_mul, prec_atom = PRECEDENCE["Mul"], PRECEDENCE["Atom"]
    symbols, zm = ring.symbols, ring.zero_monom
    sexpvs = []
    for expv, coeff in terms:
        negative = domain.is_negative(coeff)
        sexpvs.append(" - " if negative else " + ")
        if expv == zm:
            scoeff = _element_text(_Coefficient((domain, coeff)))[0]
            if negative and scoeff.startswith("-"):
                scoeff = scoeff[1:]
        else:
            if negative:
                coeff = -coeff
            if coeff != domain.one:
                scoeff = _parenthesize(_Coefficient((domain, coeff)), prec_mul, True)
            else:
                scoeff = ''
        sexpv = []
        for i in range(ring.ngens):
            exp = expv[i]
            if not exp:
                continue
            symbol = _parenthesize(symbols[i], prec_atom, True)
            if exp != 1:
                sexpv.append("%s**%s" % (symbol, exp))
            else:
                sexpv.append(symbol)
        if scoeff:
            sexpv = [scoeff] + sexpv
        sexpvs.append("*".join(sexpv))
    if sexpvs[0] in (" + ", " - "):
        head = sexpvs.pop(0)
        if head == " - ":
            sexpvs.insert(0, "-")
    return "".join(sexpvs)


def _inferred(exprs, symbols, options, fractions=False):
    from .. import Poly, ZZ, QQ, sympify
    from .._core import Expr
    single = isinstance(exprs,(Expr,str,int,float)) or not hasattr(exprs,'__iter__')
    values=[sympify(exprs)] if single else [sympify(e) for e in exprs]
    options=dict(options)
    order=options.pop('order',lex)
    given=options.pop('domain',None)
    field=options.pop('field',False)
    extension=options.pop('extension',False)
    if options:
        raise NotImplementedError('unsupported construction options: %s' % ', '.join(options))
    symbols=tuple(symbols)
    if len(symbols)==1 and isinstance(symbols[0],(tuple,list)):
        symbols=tuple(symbols[0])
    if not symbols and extension:
        symbols=tuple(sorted(set().union(*(v.free_symbols for v in values)),key=str))
    elif not symbols:
        # The reference ``_parallel_dict_from_expr`` infers the generators,
        # including non-symbol ones such as ``exp(1/x)`` or ``x**(1/3)``.
        from .polytools import parallel_poly_from_expr
        parts=[]
        for v in values:
            parts.extend(v.as_numer_denom() if fractions else (v,))
        kw={}
        if given is not None:kw['domain']=given
        if extension:kw['extension']=extension
        try:
            polys,_=parallel_poly_from_expr(parts,**kw)
        except Exception:
            polys=None
        if polys:
            symbols=tuple(polys[0].gens)
            if given is None:
                domain=polys[0].domain
                if field and not domain.is_Field:
                    domain=domain.get_field()
                return values,symbols,domain,order,single
        else:
            symbols=tuple(sorted(set().union(*(v.free_symbols for v in values)),key=str))
    domain=_domain(given) if given is not None else ZZ
    if given is None:
        for v in values:
            parts=v.as_numer_denom() if fractions else (v,)
            for part in parts:
                if symbols:
                    d=Poly(part,*symbols,**({'extension':extension} if extension else {})).domain
                else:
                    from .constructor import construct_domain
                    d=construct_domain(part,**({'extension':extension} if extension else {}))[0]
                domain=domain.unify(d)
        if field and not domain.is_Field:
            domain=domain.get_field()
    return values,symbols,domain,order,single


def sring(exprs, *symbols, **options):
    values,symbols,domain,order,single=_inferred(exprs,symbols,options)
    R=PolyRing(symbols,domain,order)
    result=[R(v) for v in values]
    return R,result[0] if single else result
