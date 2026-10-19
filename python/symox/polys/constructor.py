"""Infer native coefficient domains and convert scalar containers."""


def _infer_domain(values, *, field=False, extension=False, **options):
    from .. import (QQ, ZZ, QQ_I, ZZ_I, EXRAW, RealField, ComplexField,
                     Symbol, CoercionFailed, default_sort_key)
    if options:
        raise NotImplementedError('Domain inference options: ' + ', '.join(options))
    if extension is not True and extension not in (False, None):
        generators = extension if isinstance(extension, (list, tuple)) else (extension,)
        return QQ.algebraic_field(*generators)
    symbols = sorted({s for v in values for s in v.free_symbols if isinstance(s, Symbol)}, key=default_sort_key)
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
        try:
            return QQ.algebraic_field(*(v for v in values if not v.is_Rational))
        except (ValueError, NotImplementedError):
            pass
    return EXRAW


def construct_domain(obj, **args):
    """Return a native domain and coefficients with the input container shape.

    Inference covers implemented numerical and exact domains, symbolic
    polynomial rings and rational-function fields. General expressions use
    EX; ``extension=True`` enables the supported algebraic field models.
    """
    from .. import EX, sympify
    from .polyerrors import OptionError
    unknown = args.keys() - {'field', 'extension', 'composite'}
    if unknown:
        raise OptionError('Unknown domain inference options: ' + ', '.join(sorted(unknown)))
    field = args.get('field', False)
    composite = args.get('composite', None)
    extension = args.get('extension', None)
    if field not in (True, False) or composite not in (True, False, None):
        raise OptionError('field and composite must be boolean options')
    if extension is False:
        raise OptionError('extension=False is not an algebraic extension option')
    if extension is not None and extension is not True:
        raise NotImplementedError('construct_domain currently supports extension=True')
    mapping = isinstance(obj, dict)
    from ..utilities.iterables import iterable
    sequence = not mapping and iterable(obj)
    original = list(obj.values()) if mapping else list(obj) if sequence else [obj]
    coefficients = [sympify(value) for value in original]
    if composite is False and any(value.free_symbols for value in coefficients):
        domain = EX
    else:
        domain = _infer_domain(coefficients, field=field, extension=extension is True)
        if domain.is_EXRAW:
            if composite is True:
                raise NotImplementedError('Composite inference requires supported coefficient rings over native Symbols')
            domain = EX
    elements = [domain.convert(value) for value in coefficients]
    if mapping:
        result = dict(zip(obj, elements))
    elif sequence:
        result = elements
    else:
        result = elements[0]
    return domain, result
