"""Number-field methods of ``AlgebraicField``, as the reference
sympy/polys/domains/algebraicfield.py.

The native field keeps its primitive element and modulus; this module adds
the original generators (``orig_ext``), construction from a dense
coefficient list (``K([1, 2])``), ``to_alg_num`` and the arithmetic
invariants computed by the ported Round Two algorithm: ``maximal_order``,
``integral_basis``, ``discriminant``, ``primes_above`` and
``galois_group``. Results are cached per field object, as the reference
caches them on the field instance.
"""

_native_algebraic_field = None
_ff_installed = False
_state = {}


def _cache(K):
    entry = _state.get(id(K))
    if entry is None or entry[0] is not K:
        entry = (K, {'orig_ext': None, 'maximal_order': None,
                     'discriminant': None, 'nilradicals': {}, 'minpoly': None})
        _state[id(K)] = entry
    return entry[1]


def _symbolic_root_field(self, T, root, alias):
    """The reference ``QQ.algebraic_field((T, theta))``: a root named by a
    Symbol of the polynomial ``T``. The native field needs a numeric root;
    any root of an irreducible ``T`` generates an isomorphic field, so the
    last ``CRootOf`` is used (as ``alg_field_from_poly``) under the alias
    ``theta``, keeping ``T`` itself as the defining polynomial."""
    from ... import Poly, CRootOf
    if not isinstance(T, Poly):
        T = Poly(T)
    if not T.is_univariate or not (T.domain.is_ZZ or T.domain.is_QQ):
        raise NotImplementedError('a symbolic algebraic generator needs a univariate rational polynomial')
    if not T.is_irreducible:
        raise NotImplementedError('a symbolic algebraic generator needs an irreducible polynomial')
    K = _native_algebraic_field(self, CRootOf(T.as_expr(), T.degree() - 1), alias=alias if alias is not None else root)
    state = _cache(K)
    state['orig_ext'] = (root,)
    state['minpoly'] = T
    return K


def defining_poly(K):
    """The defining polynomial of the field (``K.ext.minpoly_of_element()``),
    or the polynomial given with a symbolic generator."""
    stored = _cache(K).get('minpoly')
    return stored if stored is not None else K.ext.minpoly_of_element()


def algebraic_field(self, *extension, alias=None):
    if (len(extension) == 1 and isinstance(extension[0], tuple) and len(extension[0]) == 2
            and getattr(self, 'is_QQ', False) and getattr(extension[0][1], 'is_Symbol', False)):
        return _symbolic_root_field(self, extension[0][0], extension[0][1], alias)
    K = _native_algebraic_field(self, *extension, alias=alias)
    if getattr(K, 'is_AlgebraicField', False) and getattr(self, 'is_QQ', False):
        if len(extension) == 1 and isinstance(extension[0], tuple):
            orig = extension[0][1:]
        else:
            orig = extension
        from ... import sympify
        _cache(K)['orig_ext'] = tuple(sympify(e) for e in orig)
    return K


def orig_ext(self):
    if not self.is_AlgebraicField:
        raise AttributeError('orig_ext')
    stored = _cache(self)['orig_ext']
    if stored is not None:
        return stored
    ext = self.ext
    return (ext if getattr(ext, 'alias', None) is not None else ext.root,)


def from_dense_list(K, coefficients):
    """The element whose dense coefficients, from the highest power of the
    primitive element down, are ``coefficients`` (the reference ``ANP``)."""
    dom = K.dom
    result = K.zero
    unit = K.unit
    for c in coefficients:
        result = result * unit + K.convert(dom.convert(c))
    return result


def to_alg_num(self, a):
    return self.ext.field_element([self.dom.to_sympy(c) for c in self.convert(a).to_list()])


def _do_round_two(self):
    from ..numberfields.basis import round_two
    state = _cache(self)
    ZK, dK = round_two(self, radicals=state['nilradicals'])
    state['maximal_order'] = ZK
    state['discriminant'] = dK


def maximal_order(self):
    state = _cache(self)
    if state['maximal_order'] is None:
        _do_round_two(self)
    return state['maximal_order']


def integral_basis(self, fmt=None):
    ZK = maximal_order(self)
    M = ZK.QQ_matrix
    n = M.shape[1]
    B = [from_dense_list(self, list(reversed(M[:, j].flat()))) for j in range(n)]
    if fmt == 'sympy':
        return [self.to_sympy(b) for b in B]
    elif fmt == 'alg':
        return [to_alg_num(self, b) for b in B]
    return B


def discriminant(self):
    state = _cache(self)
    if state['discriminant'] is None:
        _do_round_two(self)
    return state['discriminant']


def primes_above(self, p):
    from ..numberfields.primes import prime_decomp
    ZK = maximal_order(self)
    dK = discriminant(self)
    rad = _cache(self)['nilradicals'].get(p)
    return prime_decomp(p, ZK=ZK, dK=dK, radical=rad)


def galois_group(self, by_name=False, max_tries=30, randomize=False):
    return defining_poly(self).galois_group(
        by_name=by_name, max_tries=max_tries, randomize=randomize)


def _only_algebraic(function):
    def method(self, *args, **kwargs):
        if not getattr(self, 'is_AlgebraicField', False):
            raise AttributeError(function.__name__)
        return function(self, *args, **kwargs)
    method.__name__ = function.__name__
    method.__doc__ = function.__doc__
    return method


def _finite_field_attribute(name, value):
    native = getattr(__import__('symox')._core.Domain, name, None)

    def getter(self):
        if getattr(self, 'is_FiniteField', False):
            return value(self)
        if native is None:
            raise AttributeError(name)
        return native.__get__(self)
    return property(getter)


# The reference Gaussian elements: ``ZZ_I / x`` lands in QQ_I, and
# ``//``, ``%`` and ``divmod`` are ``K.div`` in the element's domain.
_native_ops = {}


def _gaussian_operand(K, other):
    if type(other) is int:
        return K.convert(other)
    parent = getattr(other, 'parent', None)
    if parent is not None and callable(parent):
        try:
            source = other.parent()
        except TypeError:
            return None
        if source == K:
            return other
        if getattr(source, 'is_ZZ', False) or getattr(source, 'is_ZZ_I', False) or \
                (getattr(source, 'is_QQ', False) and K.is_QQ_I):
            return K.convert(other)
    return None


def _zero_check(a, b, operation):
    if not b:
        if operation == 'truediv':
            raise ZeroDivisionError('%s / 0' % (a,))
        raise ZeroDivisionError('divmod(%s, 0)' % (a,))


def _gaussian_binary(name):
    native = _native_ops[name]
    reflected = name.startswith('__r')
    op = name.strip('_').lstrip('r') if reflected else name.strip('_')

    def method(self, other):
        K = self.parent()
        if K.is_ZZ_I or K.is_QQ_I:
            if op == 'truediv' and K.is_ZZ_I:
                from ... import QQ_I
                K = QQ_I
            o = _gaussian_operand(K, other)
            if o is not None:
                a = K.convert(self)
                a, b = (o, a) if reflected else (a, o)
                _zero_check(a, b, op)
                if op == 'truediv':
                    return K.exquo(a, b)
                q, r = K.div(a, b)
                return {'floordiv': q, 'mod': r, 'divmod': (q, r)}[op]
            return NotImplemented
        if native is None:
            return NotImplemented
        return native(self, other)
    method.__name__ = name
    return method


_native_exquo = None


def _exquo(self, a, b):
    try:
        return _native_exquo(self, a, b)
    except Exception as exc:
        from ... import ExactQuotientFailed
        if not isinstance(exc, ExactQuotientFailed):
            raise
        if not getattr(self, 'is_ZZ_I', False):
            # SymPy's ExactQuotientFailed(a, b, dom) text.
            try:
                text = '%s does not divide %s in %s' % (self.convert(b), self.convert(a), self)
            except Exception:
                raise exc from None
            raise ExactQuotientFailed(text) from None
        a, b = self.convert(a), self.convert(b)
        text = '(%d + %d*I) does not divide (%d + %d*I) in ZZ_I' % (
            int(b.x), int(b.y), int(a.x), int(a.y))
        raise ExactQuotientFailed(text) from None


def install():
    from ... import Domain
    from ... import _core
    global _native_exquo
    if not _native_ops:
        for name in ('__truediv__', '__rtruediv__', '__floordiv__', '__rfloordiv__',
                     '__mod__', '__rmod__', '__divmod__', '__rdivmod__'):
            _native_ops[name] = getattr(_core.DomainElement, name, None)
        for name in list(_native_ops):
            setattr(_core.DomainElement, name, _gaussian_binary(name))
        _native_exquo = Domain.exquo
        Domain.exquo = _exquo
    # The reference FiniteField attributes: ``mod`` (the characteristic),
    # ``dom`` (ZZ) and ``sym`` (symmetric representatives).
    global _ff_installed
    if not _ff_installed:
        _ff_installed = True
        Domain.mod = _finite_field_attribute('mod', lambda K: int(K.characteristic()))
        Domain.dom = _finite_field_attribute('dom', lambda K: __import__('symox').ZZ)
        Domain.sym = _finite_field_attribute('sym', lambda K: True)
    global _native_algebraic_field
    if _native_algebraic_field is None:
        _native_algebraic_field = Domain.algebraic_field
        Domain.algebraic_field = algebraic_field
    Domain.orig_ext = property(orig_ext)
    for function in (to_alg_num, maximal_order, integral_basis, discriminant,
                     primes_above, galois_group):
        setattr(Domain, function.__name__, _only_algebraic(function))
