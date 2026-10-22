"""The reference display of ground and inexact composite domain elements.

The reference elements are Python ``int`` (ZZ), ``PythonMPQ`` (QQ),
mpmath ``mpf``/``mpc`` in the field context (RR, CC) and ``PolyElement`` /
``FracElement`` over RR or CC. Other domains keep the native display.
"""
_native_repr = None
_native_str = None


def _context(ground):
    import mpmath
    ctx = mpmath.MPContext()
    ctx.prec = ground.precision
    return ctx


def _scalar(ground, value):
    """The mpmath value of a sympy number in ``ground``'s context."""
    from ... import Float
    ctx = _context(ground)
    if ground.is_ComplexField:
        re, im = value.as_real_imag()
        return ctx.make_mpc((Float(re, precision=ground.precision)._mpf_,
                             Float(im, precision=ground.precision)._mpf_))
    return ctx.make_mpf(Float(value, precision=ground.precision)._mpf_)


def _poly_str(ground, symbols, expr):
    from ... import Poly
    if expr == 0:
        return str(_scalar(ground, expr))
    # EX keeps each Float at its own precision for any RealField.
    poly = Poly(expr, *symbols, domain='EX')
    parts = []
    for monom, coeff in poly.terms():
        value = _scalar(ground, coeff)
        negative = not ground.is_ComplexField and value < 0
        parts.append(' - ' if negative else ' + ')
        if not any(monom):
            text = str(value)
            if negative and text.startswith('-'):
                text = text[1:]
        else:
            if negative:
                value = -value
            text = '' if value == 1 else str(value)
        factors = [text] if text else []
        for symbol, exponent in zip(symbols, monom):
            if not exponent:
                continue
            name = str(symbol)
            if not (symbol.is_Symbol or symbol.is_Atom):
                name = '(%s)' % name
            factors.append(name if exponent == 1 else '%s**%s' % (name, exponent))
        parts.append('*'.join(factors))
    head = parts.pop(0)
    if head == ' - ':
        parts.insert(0, '-')
    return ''.join(parts)


def _poly_level(ground, symbols, expr):
    """The reference precedence of a PolyElement: 2 atom, 1 term, 0 sum."""
    from ... import Add
    if expr.is_number or expr in symbols:
        return 2
    return 0 if isinstance(expr, Add) else 1


def _inexact_composite(element, domain):
    ground = domain.domain
    symbols = domain.symbols
    expr = element.as_expr()
    if domain.is_PolynomialRing:
        return _poly_str(ground, symbols, expr)
    numer, denom = expr.as_numer_denom()
    numer, denom = numer.expand(), denom.expand()
    if denom == 1:
        return _poly_str(ground, symbols, numer)
    top = _poly_str(ground, symbols, numer)
    if _poly_level(ground, symbols, numer) == 0:
        top = '(%s)' % top
    bottom = _poly_str(ground, symbols, denom)
    if _poly_level(ground, symbols, denom) <= 1:
        bottom = '(%s)' % bottom
    return '%s/%s' % (top, bottom)


def _display(element, use_repr):
    domain = element.parent()
    if domain.is_ZZ:
        return str(int(element))
    if domain.is_QQ:
        numer, denom = int(element.numerator), int(element.denominator)
        if use_repr:
            return 'MPQ(%d,%d)' % (numer, denom)
        return str(numer) if denom == 1 else '%d/%d' % (numer, denom)
    if domain.is_RealField or domain.is_ComplexField:
        value = _scalar(domain, element.as_expr())
        return repr(value) if use_repr else str(value)
    if domain.is_PolynomialRing or domain.is_FractionField:
        ground = domain.domain
        if ground.is_RealField or ground.is_ComplexField:
            return _inexact_composite(element, domain)
    return None


def _repr(self):
    text = _display(self, True)
    return _native_repr(self) if text is None else text


def _str(self):
    text = _display(self, False)
    return _native_str(self) if text is None else text


def install():
    from ... import _core
    global _native_repr, _native_str
    if _native_repr is None:
        _native_repr = _core.DomainElement.__repr__
        _native_str = _core.DomainElement.__str__
    _core.DomainElement.__repr__ = _repr
    _core.DomainElement.__str__ = _str
