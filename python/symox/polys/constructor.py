"""Infer native coefficient domains and convert scalar containers."""


def _composite_generators(values, symbols):
    """The symbols in the order of the native ``_construct_composite``
    generators (the reference ``_sort_gens``: ``x``, ``y``, ``z`` first),
    or ``None`` when those generators are not exactly these symbols."""
    from .. import _core
    if not symbols:
        return None
    try:
        found = _core._composite_generators(values)
    except Exception:
        return None
    if found is None or len(found[0]) != len(symbols) or set(found[0]) != symbols:
        return None
    return list(found[0])


def _infer_domain(values, *, field=False, extension=False, **options):
    from .. import (QQ, ZZ, QQ_I, ZZ_I, EXRAW, RealField, ComplexField,
                     Symbol, CoercionFailed, default_sort_key)
    if options:
        raise NotImplementedError('Domain inference options: ' + ', '.join(options))
    if extension is not True and extension not in (False, None):
        generators = extension if isinstance(extension, (list, tuple)) else (extension,)
        return QQ.algebraic_field(*generators)
    symbols = {s for v in values for s in v.free_symbols if isinstance(s, Symbol)}
    # The reference sorts the generators with ``_sort_gens``; the order
    # comes from the native ``_construct_composite`` generators.
    found = _composite_generators(values, symbols)
    symbols = found if found is not None else sorted(symbols, key=default_sort_key)
    if not symbols and all(v.is_Rational for v in values):
        return QQ if field or any(not v.is_Integer for v in values) else ZZ
    if not symbols:
        floats = [a for v in values for a in v.atoms() if getattr(a, 'is_Float', False)]
        if floats:
            precision = max(a._prec for a in floats)
            for domain in (RealField(precision), ComplexField(precision)):
                try:
                    for value in values:
                        domain.convert(value)
                except CoercionFailed:
                    continue
                return domain
            return EXRAW
        for domain in ((QQ_I,) if field else (ZZ_I, QQ_I)):
            try:
                for value in values:
                    domain.convert(value)
            except CoercionFailed:
                continue
            return domain
    if symbols:
        # Try the smallest exact polynomial ring before its fraction field.
        # Conversion itself certifies membership, including coefficient types.
        if field:
            domain = ZZ.frac_field(*symbols)
            try:
                for value in values:
                    domain.convert(value)
            except CoercionFailed:
                return EXRAW
            return domain
        for ground in (ZZ, QQ):
            domain = ground.poly_ring(*symbols)
            try:
                for value in values:
                    domain.convert(value)
            except CoercionFailed:
                continue
            return domain
        domain = ZZ.frac_field(*symbols)
        try:
            for value in values:
                domain.convert(value)
        except CoercionFailed:
            return EXRAW
        return domain
    if extension is True:
        # The reference ``_construct_algebraic``: the generators are the
        # leaves of the Add/Mul trees of the coefficients, so that
        # ``[1 + sqrt(2)]`` gives ``QQ<sqrt(2)>``.
        from ..core.sorting import ordered
        exts = set()

        def leaves(value):
            if value.is_Rational:
                return
            if value.is_Add or value.is_Mul:
                for arg in value.args:
                    leaves(arg)
            else:
                exts.add(value)
        for value in values:
            leaves(value)
        try:
            return QQ.algebraic_field(*ordered(exts))
        except (ValueError, NotImplementedError):
            pass
    return EXRAW


def _composite_domain(values, field):
    """The reference composite domain when a generator is not a Symbol.

    ``parallel_dict_from_basic`` also takes ``pi``, ``exp(a)`` or
    ``sin(x)`` as generators, so ``[pi, 2]`` gives ``ZZ[pi]``. Symbol-only
    generators keep the inference below, which certifies membership.
    """
    from .. import QQ, ZZ, Symbol, CoercionFailed
    from .. import _core
    found = _core._composite_generators(values)
    if found is None or all(isinstance(g, Symbol) for g in found[0]):
        return None
    gens, fractions = found
    for ground in (ZZ, QQ):
        domain = ground.frac_field(*gens) if fractions or field else ground.poly_ring(*gens)
        try:
            for value in values:
                domain.convert(value)
        except CoercionFailed:
            continue
        return domain
    return None


def _float_composite_domain(values, field):
    """The reference composite domain over RealField or ComplexField.

    ``_construct_composite`` takes the generators of the numerators and
    denominators, then a ComplexField ground when a Float meets a pure
    complex coefficient, else a RealField, at the largest Float precision.
    Floats are only accepted as coefficients; anything else keeps EX.
    """
    from .. import (Add, Mul, Dummy, Float, I, nsimplify, RealField, ComplexField,
                    CoercionFailed, _core)
    floats = []
    for value in values:
        for part in value.as_numer_denom():
            for term in Add.make_args(part.expand()):
                for factor in Mul.make_args(term):
                    if factor.is_Float:
                        floats.append(factor)
                    elif factor.has(Float):
                        return None
    if not floats:
        return None
    unit = Dummy('I')
    exact = []
    for value in values:
        # Only the shape matters here: any exact value of each Float will do.
        value = nsimplify(value, rational=True)
        exact.append(value.xreplace({I: unit}))
    found = _core._composite_generators(exact)
    if found is None:
        return None
    gens, fractions = found
    complexes = unit in gens
    gens = [g for g in gens if g != unit]
    if not gens:
        return None
    precision = max(f._prec for f in floats)
    ground = ComplexField(precision) if complexes else RealField(precision)
    domain = ground.frac_field(*gens) if fractions or field else ground.poly_ring(*gens)
    try:
        for value in values:
            domain.convert(value)
    except CoercionFailed:
        return None
    return domain


def _gaussian_domain(values, field):
    """The reference composite domain over ZZ_I or QQ_I, as ``ZZ_I[y]``
    or, for denominators or ``field=True``, ``ZZ_I(y)``."""
    from .. import QQ_I, ZZ_I, _core
    found = _core._gaussian_composite(values, field)
    if found is None:
        return None
    integral, ring, symbols = found
    ground = ZZ_I if integral else QQ_I
    return ground.poly_ring(*symbols) if ring else ground.frac_field(*symbols)


def _classify(coefficients, extension):
    """The reference ``_construct_simple`` classification.

    Returns ``'simple'`` when every coefficient is a rational, a Float, a
    pure complex number or (with ``extension=True``) an algebraic number,
    False for algebraic numbers mixed with Floats (EX), and None when a
    composite domain is needed.
    """
    from ..core.evalf import pure_complex
    floats = algebraics = False
    for coeff in coefficients:
        if coeff.is_Rational:
            continue
        if coeff.is_Float:
            if algebraics:
                return False
            floats = True
            continue
        parts = pure_complex(coeff)
        if parts:
            if not (parts[0].is_Rational and parts[1].is_Rational):
                floats = True
        elif extension is True and coeff.is_number and coeff.is_algebraic:
            if floats:
                return False
            algebraics = True
        else:
            return None
    return 'simple'


def construct_domain(obj, **args):
    """Return a native domain and coefficients with the input container shape.

    The control flow is the reference one: options go through the Poly
    option rules (an unknown option or ``extension=False`` raises
    OptionError, ``gens``, ``domain`` or ``greedy`` are accepted and
    ignored), simple coefficients give ZZ, QQ, ZZ_I, QQ_I, RR, CC or an
    algebraic field (only for ``extension=True``; any other extension
    value behaves like no extension), and the rest gives a composite
    domain unless ``composite=False``, otherwise EX.
    """
    from .. import EX, sympify, _core
    from .._poly_constructor import _options
    _, options = _options((), None, None, None, args)
    field = bool(options.get('field', False))
    composite = options.get('composite', None)
    extension = options.get('extension', None)
    mapping = isinstance(obj, dict)
    from ..utilities.iterables import iterable
    sequence = not mapping and iterable(obj)
    original = list(obj.values()) if mapping else list(obj) if sequence else [obj]
    coefficients = [sympify(value) for value in original]
    kind = _classify(coefficients, extension)
    if kind is False:
        domain = EX
    elif kind == 'simple':
        domain = _gaussian_domain(coefficients, field)
        if domain is None:
            domain = _infer_domain(coefficients, field=field, extension=extension is True)
    elif composite is False:
        domain = EX
    else:
        has_float = any(a.is_Float for v in coefficients for a in v.atoms())
        domain = None
        if not has_float:
            domain = _gaussian_domain(coefficients, field)
            if domain is None:
                domain = _composite_domain(coefficients, field)
            if domain is None and any(value.free_symbols for value in coefficients):
                domain = _infer_domain(coefficients, field=field)
        else:
            domain = _float_composite_domain(coefficients, field)
        if domain is None:
            domain = EX
    if domain.is_EXRAW:
        domain = EX
    elements = [domain.convert(value) for value in coefficients]
    if mapping:
        result = dict(zip(obj, elements))
    elif sequence:
        result = elements
    else:
        result = elements[0]
    return domain, result
