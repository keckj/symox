"""Polynomial evaluation with explicit native coefficient domains."""
import operator

from .. import _core

_native_as_expr = _core.Poly.as_expr
_native_eval = _core.Poly.eval


def _level(polynomial, generator):
    from .. import sympify
    from .polyerrors import GeneratorsError
    try:
        index = operator.index(generator)
    except TypeError:
        generator = sympify(generator)
        try:
            return polynomial.gens.index(generator)
        except ValueError:
            raise GeneratorsError('Not a polynomial generator: %s' % generator) from None
    count = len(polynomial.gens)
    if not -count <= index < count:
        raise GeneratorsError('Polynomial generator index is out of range')
    return index % count


def as_expr(polynomial, *generators):
    if not generators:
        return _native_as_expr(polynomial)
    from .. import sympify
    from .polyerrors import GeneratorsError
    from ._poly_methods import _expr_from_dict
    if len(generators) == 1 and isinstance(generators[0], dict):
        replacements = list(polynomial.gens)
        for old, new in generators[0].items():
            old = sympify(old)
            if old not in polynomial.gens:
                raise GeneratorsError('Not a polynomial generator: %s' % old)
            replacements[polynomial.gens.index(old)] = sympify(new)
    else:
        if len(generators) != len(polynomial.gens):
            raise GeneratorsError('An expression value is required for every generator')
        replacements = [sympify(g) for g in generators]
    return _expr_from_dict(polynomial.as_dict(), replacements)


def _one(polynomial, index, value, auto):
    from .. import sympify
    from .constructor import construct_domain
    from .polyerrors import CoercionFailed, DomainError, UnificationFailed
    from ._poly_methods import _from_dict
    value = sympify(value)
    try:
        domain = polynomial.domain
    except NotImplementedError:
        if not auto:
            raise NotImplementedError('Retaining this legacy coefficient domain is unsupported') from None
        return _native_eval(polynomial, polynomial.gens[index], value)
    original_domain = domain
    try:
        point = domain.convert(value)
    except CoercionFailed:
        if not auto:
            raise DomainError('Cannot evaluate at %s in %s' % (value, domain)) from None
        point_domain, _ = construct_domain([value])
        try:
            domain = domain.unify_with_symbols(point_domain, polynomial.gens)
        except DomainError as error:
            raise UnificationFailed(str(error)) from error
        point = domain.convert(value)
    if all(generator.is_Symbol for generator in polynomial.gens):
        if len(polynomial.gens) == 1:
            # Scalar substitution and expansion stay in the Rust expression
            # engine; final conversion certifies the selected coefficient field.
            result = _native_as_expr(polynomial).subs(polynomial.gens[index], domain.to_sympy(point)).expand()
            return domain.to_sympy(domain.convert(result))
        representative = (polynomial if domain == original_domain else polynomial.set_domain(domain))
        result = _native_eval(representative, representative.gens[index], domain.to_sympy(point))
        try:
            same_domain = result.domain == domain
        except NotImplementedError:
            same_domain = False
        return result if same_domain else result.set_domain(domain)
    # Expression generators must retain their polynomial positions: substituting
    # into an expanded expression could lose a generator such as x + 1.
    terms, powers = {}, {0: domain.one}
    for monomial, coefficient in polynomial.as_dict().items():
        exponent = monomial[index]
        if exponent not in powers:
            powers[exponent] = point**exponent
        monomial = monomial[:index] + monomial[index+1:]
        term = domain.convert(coefficient)*powers[exponent]
        terms[monomial] = terms.get(monomial, domain.zero) + term
    generators = polynomial.gens[:index] + polynomial.gens[index+1:]
    if not generators:
        return domain.to_sympy(terms.get((), domain.zero))
    return _from_dict(polynomial, {m: domain.to_sympy(c) for m, c in terms.items()},
                      generators, domain=domain)


def eval(polynomial, x, a=None, auto=True):
    if a is not None:
        return _one(polynomial, _level(polynomial, x), a, auto)
    if isinstance(x, dict):
        result = polynomial
        for generator, value in x.items():
            if not isinstance(result, _core.Poly):
                from .polyerrors import GeneratorsError
                raise GeneratorsError('All polynomial generators have already been evaluated')
            result = _one(result, _level(result, generator), value, auto)
        return result
    if isinstance(x, (tuple, list)):
        if len(x) > len(polynomial.gens):
            raise ValueError('too many values provided')
        result = polynomial
        for value in x:
            result = _one(result, 0, value, auto)
        return result
    return _one(polynomial, 0, x, auto)


def evaluate_call(polynomial, *values):
    return eval(polynomial, values)


def install():
    _core.Poly.as_expr = as_expr
    _core.Poly.eval = eval
    _core.Poly.__call__ = evaluate_call
