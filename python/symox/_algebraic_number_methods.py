"""Structural AlgebraicNumber protocols and native scalar calculations.

The field generator, coordinates and alias remain attached to unchanged
numbers. Mathematical operations use the existing exact scalar embedding;
they do not invoke SymPy as a computational backend.
"""
from collections import defaultdict


def to_sympy(self):
    import sympy
    return sympy.AlgebraicNumber(self.root.to_sympy(),
                                 [c.to_sympy() for c in self.coeffs()],
                                 alias=None if self.alias is None else self.alias.to_sympy())


def class_key():
    return 5, 0, 'AlgebraicNumber'


def sort_key(self, order=None):
    from . import S, Tuple
    def key(value):
        return Tuple(*value).sort_key(order=order) if isinstance(value, tuple) else value.sort_key(order=order)
    return class_key(), (len(self.args), tuple(key(a) for a in self.args)), S.One.sort_key(), S.One


def compare(self, other):
    from . import AlgebraicNumber
    if not isinstance(other, AlgebraicNumber):
        a, b = self.class_key(), other.class_key()
        return (a > b)-(a < b)
    if len(self.args) != len(other.args):
        return (len(self.args) > len(other.args))-(len(self.args) < len(other.args))
    def compare_arg(a, b):
        if isinstance(a, tuple):
            if len(a) != len(b):
                return (len(a) > len(b))-(len(a) < len(b))
            for x, y in zip(a, b):
                result = x.compare(y)
                if result:
                    return result
            return 0
        return a.compare(b)
    for a, b in zip(self.args, other.args):
        result = compare_arg(a, b)
        if result:
            return result
    return 0


def _identity(self, *args, **kwargs):
    return self


def _one_self(self, *args, **kwargs):
    from . import S
    return S.One, self


def _self_one(self):
    from . import S
    return self, S.One


def _self_zero(self, *args, **kwargs):
    from . import S
    return self, S.Zero


def _none(self):
    return None


def as_coeff_add(self, *deps):
    from . import S
    return (self, ()) if deps and not self.has_free(*deps) else (S.Zero, (self,))


def as_coeff_mul(self, *deps, rational=True, **kwargs):
    from . import S
    if kwargs:
        raise TypeError('unsupported coefficient options: '+', '.join(kwargs))
    return (self, ()) if deps and not self.has_free(*deps) else (S.One, (self,))


def as_coeff_Add(self, rational=False, deps=None):
    from . import S
    return S.Zero, self


def args_cnc(self, cset=False, warn=True, split_1=True):
    return [{self} if cset else [self], []]


def as_coefficients_dict(self, *syms):
    from . import S
    return defaultdict(int, {S.One: self} if syms and not self.has_free(*syms) else {self: S.One})


def as_terms(self):
    value = complex(self)
    return [(self, ((value.real, value.imag), (), ()))], []


def as_ordered_terms(self, order=None, data=False):
    if data:
        return self.as_terms()
    return [self]


def as_ordered_factors(self, order=None):
    return [self]


def as_powers_dict(self):
    from . import S
    return defaultdict(int, {self: S.One})


def as_independent(self, *deps, **hint):
    from . import S, Symbol
    # Literal scalar symbols are screened against free_symbols by Expr's
    # contract; aliases and defining coordinates are structural metadata.
    other = [d for d in deps if not isinstance(d, Symbol)]
    dependent = bool(other) and self.has(*other)
    identity = S.Zero if hint.get('as_Add', False) else S.One
    return (identity, self) if dependent else (self, identity)


def coeff(self, x, n=1, right=False, _first=True):
    from . import S
    if n == 0:
        return S.Zero if self == x else self
    return S.One if n == 1 and self == x else S.Zero


def as_coeff_exponent(self, x):
    from . import S
    return (S.One, S.One) if self == x else (self, S.Zero)


def extract_additively(self, c):
    from . import S
    if c == 0:
        return self
    return S.Zero if self == c else None


def extract_multiplicatively(self, c):
    from . import S
    if c == 1:
        return self
    return S.One if self == c else None


def as_coefficient(self, expr):
    return self.extract_multiplicatively(expr)


def is_same(self, other, approx=None):
    from . import AlgebraicNumber
    if not isinstance(other, AlgebraicNumber) or len(self.args) != len(other.args):
        return False
    for a, b in zip(self.args, other.args):
        if isinstance(a, tuple):
            if len(a) != len(b) or not all(x.is_same(y, approx=approx) for x, y in zip(a, b)):
                return False
        elif not a.is_same(b, approx=approx):
            return False
    return True


def dummy_eq(self, other, symbol=None):
    return self == other


def count_ops(self, visual=False):
    from ._operation_counts import count_ops
    return count_ops(self.args, visual)


def _rebuild(self, root, coeffs, alias):
    from . import AlgebraicNumber
    if root == self.root and list(coeffs) == self.coeffs() and alias == self.alias:
        return self
    return AlgebraicNumber(root, coeffs, alias=alias)


def _map_arguments(self, transform):
    root = transform(self.root)
    coords = tuple(transform(c) for c in self.coeffs())
    alias = None if self.alias is None else transform(self.alias)
    return _rebuild(self, root, coords, alias)


def xreplace(self, rule):
    from . import AlgebraicNumber
    if self in rule:
        return rule[self]
    scalar = {a: b for a, b in rule.items() if not isinstance(a, AlgebraicNumber)}
    return _map_arguments(self, lambda a: a.xreplace(scalar))


def subs(self, *args, **kwargs):
    from . import AlgebraicNumber, sympify
    simultaneous = kwargs.pop('simultaneous', False)
    if kwargs:
        raise TypeError('unsupported substitution options: '+', '.join(kwargs))
    if len(args) == 2:
        pairs = [args]
    elif len(args) == 1:
        if hasattr(args[0], 'items'):
            from ._matrixexpr_methods import _walk
            pairs = [(old if isinstance(old, AlgebraicNumber) else sympify(old), new)
                     for old, new in args[0].items()]
            pairs.sort(key=lambda pair: (-sum(1 for _ in _walk(pair[0])), pair[0].sort_key()))
        else:
            pairs = list(args[0])
    else:
        raise TypeError('subs expects (old, new), a mapping, or a sequence of pairs')
    if simultaneous:
        for old, new in pairs:
            if self == old:
                return new if isinstance(new, AlgebraicNumber) else sympify(new)
        scalar = [(old, new) for old, new in pairs if not isinstance(old, AlgebraicNumber)]
        return _map_arguments(self, lambda a: a.subs(scalar, simultaneous=True))
    result = self
    for old, new in pairs:
        if result == old:
            result = new if isinstance(new, AlgebraicNumber) else sympify(new)
        elif isinstance(result, AlgebraicNumber):
            if not isinstance(old, AlgebraicNumber):
                result = _map_arguments(result, lambda a: a.subs(old, new))
        else:
            result = result.subs(old, new)
    return result


def matches(self, expr, repl_dict=None, old=False):
    from . import AlgebraicNumber
    if not isinstance(expr, AlgebraicNumber) or len(self.args) != len(expr.args):
        return None
    bindings = dict(repl_dict or {})
    for pattern, value in zip(self.args, expr.args):
        if isinstance(pattern, tuple):
            if pattern != value:
                return None
            continue
        result = pattern.xreplace(bindings).matches(value,old=old)
        if result is None:
            return None
        bindings.update(result)
    return bindings


def match(self, pattern, old=False):
    from . import AlgebraicNumber, Expr
    if isinstance(pattern, AlgebraicNumber):
        return pattern.matches(self, old=old)
    if isinstance(pattern, Expr) and pattern.is_Symbol and str(pattern) == pattern.name+'_':
        if any(self.has(excluded) for excluded in pattern.exclude):return None
        result = self.as_expr().match(pattern,old=old)
        return None if result is None else {key: self if key == pattern else value for key, value in result.items()}
    return None


def _kernel(name):
    def method(self, *args, **kwargs):
        from . import Expr
        original = self.as_expr()
        result = getattr(original, name)(*args, **kwargs)
        # Preserve the defining field and alias when the scalar operation
        # leaves the exact embedded value unchanged.
        if name in ('diff', 'integrate', 'invert', 'dir', 'round'):
            return result
        def restore(value):
            if isinstance(value, Expr):
                return self if value == original else value
            if isinstance(value, (list, tuple)):
                return type(value)(restore(v) for v in value)
            if isinstance(value, dict):
                return {restore(k): restore(v) for k, v in value.items()}
            return value
        return restore(result)
    method.__name__ = name
    return method


def expand(self, deep=True, **hints):
    if not deep:
        return self
    return _map_arguments(self, lambda a: a.expand(deep=True, **hints))


def apart(self, x=None, full=False, **options):
    from . import Symbol
    if x is not None and not isinstance(x, Symbol):
        raise ValueError('partial fractions require a scalar symbol generator')
    return self


def collect(self, syms, func=None, evaluate=True, exact=False, distribute_order_term=True):
    from . import S
    if evaluate:
        return self if func is None else func(self)
    return {S.One: self if func is None else func(self)}


def conjugate(self):
    from . import conjugate
    if self.is_real is True:
        return self
    return _rebuild(self, conjugate(self.root), self.coeffs(), self.alias)


def as_real_imag(self, deep=True, **hints):
    from . import S
    if self.is_real is True:
        return self, S.Zero
    return self.as_expr().as_real_imag(deep=deep, **hints)


def equals(self, other, failing_expression=False):
    from . import AlgebraicNumber
    if isinstance(other, AlgebraicNumber):
        other = other.as_expr()
    return self.as_expr().equals(other, failing_expression=failing_expression)


def is_constant(self, *wrt, **flags):
    return True


def _algebraic_query(self, *syms):
    if self in syms or not self.has_free(*syms):
        return True
    return None


def is_meromorphic(self, x, a):
    from . import Symbol
    if not isinstance(x, Symbol):
        raise TypeError('the meromorphic variable must be a Symbol')
    return True


def is_hypergeometric(self, k):
    return True


def series(self, x=None, x0=0, n=6, dir='+', logx=None, cdir=0):
    if dir not in ('+', '-'):
        raise ValueError("series direction must be '+' or '-'")
    return self


def lseries(self, *args, **kwargs):
    yield self


def leadterm(self, x, logx=None, cdir=0):
    from . import S
    return self, S.Zero


def taylor_term(self, n, x, *previous_terms):
    from . import S
    return self if n == 0 else S.Zero


def to_algebraic_integer(self):
    from . import AlgebraicNumber
    leading = self.minpoly.LC()
    return self if leading == 1 else AlgebraicNumber(leading*self.root, self.coeffs())


def fromiter(cls, args, **assumptions):
    return cls(*args, **assumptions)


def install():
    from . import AlgebraicNumber, S
    from ._matrixexpr_methods import atoms, find, count, has, has_free, has_xfree
    from ._matrixexpr_algebra import replace
    methods = {
        'to_sympy': to_sympy, '_sympy_': to_sympy,
        'class_key': staticmethod(class_key), 'sort_key': sort_key, 'compare': compare,
        'fromiter': classmethod(fromiter), 'atoms': atoms, 'find': find, 'count': count,
        'has': has, 'has_free': has_free, 'has_xfree': has_xfree,
        'is_same': is_same, 'dummy_eq': dummy_eq, 'count_ops': count_ops,
        'copy': _identity, 'doit': _identity, 'as_dummy': _identity,
        'transpose': _identity, 'conjugate': conjugate, 'adjoint': conjugate,
        'as_coeff_add': as_coeff_add, 'as_coeff_Add': as_coeff_Add,
        'as_coeff_mul': as_coeff_mul, 'as_coeff_Mul': _one_self,
        'as_content_primitive': _one_self, 'primitive': _one_self,
        'as_base_exp': _self_one, 'as_numer_denom': _self_one,
        'as_coeff_exponent': as_coeff_exponent, 'as_real_imag': as_real_imag,
        'as_independent': as_independent, 'args_cnc': args_cnc,
        'as_terms': as_terms, 'as_ordered_terms': as_ordered_terms,
        'as_ordered_factors': as_ordered_factors, 'as_powers_dict': as_powers_dict,
        'as_coefficients_dict': as_coefficients_dict, 'coeff': coeff,
        'could_extract_minus_sign': lambda self: False,
        'getO': _none, 'getn': _none, 'removeO': _identity,
        'extract_additively': extract_additively, 'extract_multiplicatively': extract_multiplicatively,
        'as_coefficient': as_coefficient, 'extract_branch_factor': _self_zero,
        'subs': subs, 'xreplace': xreplace, 'replace': replace,
        'match': match, 'matches': matches, 'expand': expand, 'equals': equals,
        'apart': apart, 'collect': collect,
        'is_constant': is_constant, 'is_polynomial': _algebraic_query,
        'is_rational_function': _algebraic_query, 'is_algebraic_expr': _algebraic_query,
        'is_hypergeometric': is_hypergeometric, 'is_meromorphic': is_meromorphic,
        'series': series, 'nseries': series, 'aseries': _identity,
        'lseries': lseries, 'as_leading_term': _identity, 'leadterm': leadterm,
        'taylor_term': taylor_term, 'fps': _identity, 'fourier_series': _identity,
        'to_algebraic_integer': to_algebraic_integer,
    }
    for name in ('cancel', 'combsimp', 'diff', 'dir', 'factor',
                 'gammasimp', 'integrate', 'invert', 'limit', 'normal', 'nsimplify',
                 'powsimp', 'radsimp', 'ratsimp', 'rcall', 'refine', 'rewrite',
                 'round', 'separate', 'simplify', 'together', 'trigsimp'):
        methods[name] = _kernel(name)
    for name, method in methods.items():
        setattr(AlgebraicNumber, name, method)
