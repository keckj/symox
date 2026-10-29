"""Contiguous hypergeometric operators backed by exact Euler polynomials.

Formula lookup uses a finite elementary catalogue and the native expansion
engine. Integer shift plans are bounded; singular shifts raise ValueError.
"""
from collections import Counter, defaultdict
from itertools import permutations
from symox._api_exports import export as _export
_export(globals(), 'symox', 'hyperexpand')


def _S():
    import symox
    return symox


def _params(values):
    S = _S()
    return tuple(S.expand(S.sympify(value)) for value in values)


def _mod1(value):
    # SymPy's _mod1: only the numeric part of the parameter is reduced.
    S = _S()
    value = S.sympify(value)
    if value.is_Number:
        return S.Mod(value, 1)
    c, rest = value.as_coeff_Add()
    return S.Mod(c, 1)+rest


class Hyper_Function:
    def __init__(self, ap, bq):
        self.ap, self.bq = _params(ap), _params(bq)

    @property
    def args(self):
        return self.ap, self.bq

    @property
    def sizes(self):
        return len(self.ap), len(self.bq)

    @property
    def gamma(self):
        return sum(a.is_integer is True and a.is_negative is True for a in self.ap)

    def __call__(self, z):
        return _S().hyper(self.ap, self.bq, z)

    def has(self, *symbols):
        # SymPy's Hyper_Function is an Expr with args (Tuple(*ap), Tuple(*bq)).
        Tuple = _S().Tuple
        return any(Tuple(*group).has(*symbols) for group in self.args)

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.args == other.args

    def __hash__(self):
        return hash((type(self), self.args))

    def __repr__(self):
        # The reference prints the Expr form, one Tuple per parameter group.
        def group(values):
            values = tuple(values)
            return '(%s,)' % values[0] if len(values) == 1 else '(%s)' % ', '.join(map(str, values))
        return '%s(%s)' % (type(self).__name__, ', '.join(group(g) for g in self.args))

    __str__ = __repr__

    def build_invariants(self):
        S = _S()
        return (self.gamma, *(tuple(sorted(Counter(_mod1(a) for a in group).items(),
                                          key=lambda pair: S.default_sort_key(pair[0]))) for group in self.args))

    def difficulty(self, func):
        try:
            pairs = [_pair_parameters(a, b) for a, b in zip(self.args, func.args)]
        except ValueError:
            return -1
        if self.gamma != func.gamma:
            return -1
        return sum(abs(int(_S().simplify(a-b))) for group in pairs for a, b in group)

    def _is_suitable_origin(self):
        # SymPy: no ai equals bj + n (n a non-negative integer), no ai is
        # zero and no bj is a non-positive integer.
        for a in self.ap:
            for b in self.bq:
                if (a - b).is_integer and (a - b).is_negative is False:
                    return False
        if any(a == 0 for a in self.ap):
            return False
        return not any(b.is_integer and b.is_nonpositive for b in self.bq)


class G_Function(Hyper_Function):
    def __init__(self, an, ap, bm, bq):
        self.an, self.ap, self.bm, self.bq = map(_params, (an, ap, bm, bq))

    @property
    def args(self):
        return self.an, self.ap, self.bm, self.bq

    @property
    def signature(self):
        return tuple(map(len, self.args))

    def __call__(self, z):
        return _S().meijerg(*self.args, z)

    def compute_buckets(self):
        S = _S()
        result = []
        for group, reverse in zip(self.args, (True, False, False, True)):
            buckets = defaultdict(list)
            for value in group:
                buckets[_mod1(value)].append(value)
            for values in buckets.values():
                origin = values[0]
                values.sort(key=lambda value: S.default_sort_key(value-origin), reverse=reverse)
            result.append(dict(buckets))
        return tuple(result)


def _annihilator(func, z, d):
    S = _S()
    if isinstance(func, G_Function):
        sign = -1 if (len(func.ap)-len(func.bm)) % 2 else 1
        return S.expand(sign*z*S.prod(d-a+1 for a in func.an+func.ap)-S.prod(d-b for b in func.bm+func.bq))
    return S.expand(d*S.prod(d+b-1 for b in func.bq)-z*S.prod(d+a for a in func.ap))


_X = []


def _xvar():
    # SymPy's module-level generator ``_x = Dummy('x')``.
    if not _X:
        _X.append(_S().Dummy('x'))
    return _X[0]


class Operator:
    def __init__(self, polynomial=0, generator=None):
        S = _S()
        self._d = generator or _xvar()
        self._poly = S.Poly(polynomial, self._d)

    def apply(self, obj, op):
        coefficients = list(reversed(self._poly.all_coeffs()))
        value, result = obj, coefficients[0]*obj
        for coefficient in coefficients[1:]:
            value = op(value)
            result += coefficient*value
        return result


def _native():
    from symox import _core
    return _core


def _from_record(record):
    # Rebuild the public operator object of a native shift operator.
    kind, groups, index, coefficients = record
    cls = globals()[kind]
    obj = object.__new__(cls)
    _set_poly(obj, coefficients)
    if index is not None:
        names = ('_ap', '_bq') if len(groups) == 2 else ('_an', '_ap', '_bm', '_bq')
        for name, group in zip(names, groups):
            setattr(obj, name, list(group))
        obj._i = _S().Integer(index)
    return obj


def _set_poly(obj, coefficients):
    S = _S(); x = _xvar()
    obj._d = x
    obj._poly = S.Poly(S.Add(*[c*x**k for k, c in enumerate(coefficients)]), x)


def _shift(obj, kind, value):
    value = _S().sympify(value)
    _set_poly(obj, _native()._hyperexpand_operator(kind, [[value]], 0, 0)[3])


def _unshift(obj, kind, groups, i, z):
    # ``list.pop(i)`` semantics of the reference for the shifted group.
    import operator
    S = _S()
    groups = [S.sympify(group) for group in groups]
    i = S.sympify(i)
    names = ('_ap', '_bq') if len(groups) == 2 else ('_an', '_ap', '_bm', '_bq')
    for name, group in zip(names, groups):
        setattr(obj, name, group)
    obj._i = i
    position = {'UnShiftA': 0, 'UnShiftB': 1, 'MeijerUnShiftA': 2, 'MeijerUnShiftB': 0,
                'MeijerUnShiftC': 3, 'MeijerUnShiftD': 1}[kind]
    values = [list(group) for group in groups]
    index = operator.index(i)
    if index < 0:
        index += len(values[position])
    if not 0 <= index < len(values[position]):
        raise IndexError('pop index out of range')
    _set_poly(obj, _native()._hyperexpand_operator(kind, values, index, S.sympify(z))[3])


class MultOperator(Operator):
    def __init__(self, p):
        super().__init__(p)


class ShiftA(Operator):
    """ Increment an upper index. """
    def __init__(self, ai):
        _shift(self, 'ShiftA', ai)

    def __str__(self):
        return '<Increment upper %s.>' % (1/self._poly.all_coeffs()[0])


class ShiftB(Operator):
    """ Decrement a lower index. """
    def __init__(self, bi):
        _shift(self, 'ShiftB', bi)

    def __str__(self):
        return '<Decrement lower %s.>' % (1/self._poly.all_coeffs()[0] + 1)


class UnShiftA(Operator):
    """ Decrement an upper index. """
    def __init__(self, ap, bq, i, z):
        _unshift(self, 'UnShiftA', (ap, bq), i, z)

    def __str__(self):
        return '<Decrement upper index #%s of %s, %s.>' % (self._i, self._ap, self._bq)


class UnShiftB(Operator):
    """ Increment a lower index. """
    def __init__(self, ap, bq, i, z):
        _unshift(self, 'UnShiftB', (ap, bq), i, z)

    def __str__(self):
        return '<Increment lower index #%s of %s, %s.>' % (self._i, self._ap, self._bq)


class MeijerShiftA(Operator):
    """ Increment an upper b index. """
    def __init__(self, bi):
        _shift(self, 'MeijerShiftA', bi)

    def __str__(self):
        return '<Increment upper b=%s.>' % (self._poly.all_coeffs()[1])


class MeijerShiftB(Operator):
    """ Decrement an upper a index. """
    def __init__(self, bi):
        _shift(self, 'MeijerShiftB', bi)

    def __str__(self):
        return '<Decrement upper a=%s.>' % (1 - self._poly.all_coeffs()[1])


class MeijerShiftC(Operator):
    """ Increment a lower b index. """
    def __init__(self, bi):
        _shift(self, 'MeijerShiftC', bi)

    def __str__(self):
        return '<Increment lower b=%s.>' % (-self._poly.all_coeffs()[1])


class MeijerShiftD(Operator):
    """ Decrement a lower a index. """
    def __init__(self, bi):
        _shift(self, 'MeijerShiftD', bi)

    def __str__(self):
        return '<Decrement lower a=%s.>' % (self._poly.all_coeffs()[1] + 1)


class _MeijerInverse(Operator):
    _label = ''

    def __init__(self, an, ap, bm, bq, i, z):
        _unshift(self, type(self).__name__, (an, ap, bm, bq), i, z)

    def __str__(self):
        return '<%s index #%s of %s, %s, %s, %s.>' % (self._label, self._i, self._an,
                                                    self._ap, self._bm, self._bq)


class MeijerUnShiftA(_MeijerInverse):
    """ Decrement an upper b index. """
    _label = 'Decrement upper b'


class MeijerUnShiftB(_MeijerInverse):
    """ Increment an upper a index. """
    _label = 'Increment upper a'


class MeijerUnShiftC(_MeijerInverse):
    """ Decrement a lower b index. """
    _label = 'Decrement lower b'


class MeijerUnShiftD(_MeijerInverse):
    """ Increment a lower a index. """
    _label = 'Increment lower a'


class ReduceOrder(Operator):
    def __new__(cls, ai, bj):
        S = _S(); ai, bj = S.sympify(ai), S.sympify(bj)
        n = S.simplify(ai-bj)
        if not n.is_Integer or n < 0 or bj.is_integer is True and bj.is_nonpositive is True:
            return None
        if int(n) > 256:
            raise ValueError('Contiguous order reduction exceeds 256 factors')
        return super().__new__(cls)

    def __init__(self, ai, bj):
        S = _S(); ai, bj = S.sympify(ai), S.sympify(bj)
        d = _xvar()
        self._a, self._b = ai, bj
        super().__init__(S.prod((d+bj+k)/(bj+k) for k in range(int(ai-bj))), d)

    @classmethod
    def _meijer(cls, b, a, sign):
        S = _S(); b, a = S.sympify(b), S.sympify(a)
        n = S.simplify(b-a)
        if not n.is_Integer or n < 0:
            return None
        if int(n) > 256:
            raise ValueError('Contiguous order reduction exceeds 256 factors')
        obj = object.__new__(cls); d = _xvar()
        Operator.__init__(obj, S.prod(sign*d+a+k for k in range(int(n))), d)
        obj._a, obj._b = b, a
        return obj

    @classmethod
    def meijer_minus(cls, b, a):
        return cls._meijer(b, a, -1)

    @classmethod
    def meijer_plus(cls, a, b):
        return cls._meijer(1-_S().sympify(a), 1-_S().sympify(b), 1)


def _reduce_order(ap, bq, gen, key):
    # The reference greedy pairing over key-sorted parameters.
    ap = sorted(ap, key=key)
    bq = sorted(bq, key=key)
    nap, operators = [], []
    for a in ap:
        op = None
        for i in range(len(bq)):
            op = gen(a, bq[i])
            if op is not None:
                bq.pop(i)
                break
        if op is None:
            nap.append(a)
        else:
            operators.append(op)
    return nap, bq, operators


def reduce_order(func):
    S = _S()
    nap, nbq, operators = _reduce_order(func.ap, func.bq, ReduceOrder, S.default_sort_key)
    return Hyper_Function(nap, nbq), operators


def reduce_order_meijer(func):
    S = _S()
    nan, nbq, ops1 = _reduce_order(func.an, func.bq, ReduceOrder.meijer_plus,
                                   lambda x: S.default_sort_key(-x))
    nbm, nap, ops2 = _reduce_order(func.bm, func.ap, ReduceOrder.meijer_minus,
                                   S.default_sort_key)
    return G_Function(nan, nap, nbm, nbq), ops1 + ops2


def make_simp(z):
    S = _S()
    return lambda value: S.cancel(S.expand(value))


def make_derivative_operator(M, z):
    def derivative(C):
        return (z*C.diff(z)+C*M).applyfunc(make_simp(z))
    return derivative


def apply_operators(obj, ops, op):
    for operator in reversed(ops):
        obj = operator.apply(obj, op)
    return obj


def _pair_parameters(target, origin):
    S = _S()
    if len(target) != len(origin) or len(target) > 8:
        raise ValueError('Shift planning requires equal parameter counts, at most eight per group')
    best = None
    for permuted in permutations(origin):
        differences = [S.simplify(a-b) for a, b in zip(target, permuted)]
        if all(d.is_Integer for d in differences):
            score = sum(abs(int(d)) for d in differences)
            if best is None or score < best[0]:
                best = score, list(zip(target, permuted))
    if best is None:
        raise ValueError('Parameters must differ by integers within corresponding groups')
    return best[1]


def devise_plan(target, origin, z):
    # Native port of the reference bucket-by-bucket plan.
    S = _S()
    try:
        records = _native()._hyperexpand_devise_plan(
            list(target.ap), list(target.bq), list(origin.ap), list(origin.bq), S.sympify(z))
    except ValueError as error:
        if str(error) == 'not reachable':
            raise ValueError('%s not reachable from %s' % (target, origin)) from None
        raise
    return [_from_record(record) for record in records]


def devise_plan_meijer(fro, to, z):
    # Port of SymPy 1.14: inverse-shift when possible, shift otherwise, and
    # give up when no progress can be made.

    def try_shift(f, t, shifter, diff, counter):
        for idx, (a, b) in enumerate(zip(f, t)):
            if ((a - b).is_integer and (b - a)/diff > 0 and
                    all(a != x for x in counter)):
                sh = shifter(idx)
                f[idx] += diff
                return sh
        return None

    fan, fap, fbm, fbq = list(fro.an), list(fro.ap), list(fro.bm), list(fro.bq)
    steps = (
        (fan, to.an, lambda i: MeijerUnShiftB(fan, fap, fbm, fbq, i, z), 1, lambda: fbm + fbq),
        (fap, to.ap, lambda i: MeijerUnShiftD(fan, fap, fbm, fbq, i, z), 1, lambda: fbm + fbq),
        (fbm, to.bm, lambda i: MeijerUnShiftA(fan, fap, fbm, fbq, i, z), -1, lambda: fan + fap),
        (fbq, to.bq, lambda i: MeijerUnShiftC(fan, fap, fbm, fbq, i, z), -1, lambda: fan + fap),
        (fan, to.an, lambda i: MeijerShiftB(fan[i]), -1, list),
        (fap, to.ap, lambda i: MeijerShiftD(fap[i]), -1, list),
        (fbm, to.bm, lambda i: MeijerShiftA(fbm[i]), 1, list),
        (fbq, to.bq, lambda i: MeijerShiftC(fbq[i]), 1, list),
    )
    ops = []
    change = True
    while change:
        change = False
        for f, t, shifter, diff, counter in steps:
            op = try_shift(f, t, shifter, diff, counter())
            if op is not None:
                ops.append(op)
                change = True
                break
    if (fan != list(to.an) or fap != list(to.ap) or fbm != list(to.bm)
            or fbq != list(to.bq)):
        raise NotImplementedError('Could not devise plan.')
    ops.reverse()
    return ops


class Formula:
    def __init__(self, func, z, res, symbols, B=None, C=None, M=None):
        S = _S()
        self.func, self.z, self.symbols = func, S.sympify(z), list(symbols)
        if B is not None:
            self.B, self.C, self.M = B, C, M
            return
        d = S.Dummy('Euler_operator')
        polynomial = S.Poly(_annihilator(func, z, d), d)
        dimension = int(polynomial.degree())
        if not 1 <= dimension <= 32:
            raise ValueError('Formula differential basis requires order one to 32')
        basis = [func(z) if res is None else S.sympify(res)]
        for _ in range(dimension-1):
            basis.append(z*S.diff(basis[-1], z))
        self.B = S.Matrix(basis)
        self.C = S.Matrix([[1]+[0]*(dimension-1)])
        self.M = S.zeros(dimension)
        for k in range(dimension-1):
            self.M[k, k+1] = 1
        leading = polynomial.LC()
        for k in range(dimension):
            self.M[dimension-1, k] = S.cancel(-polynomial.nth(k)/leading)

    @property
    def closed_form(self):
        return (self.C*self.B)[0]

    def find_instantiations(self, func):
        S = _S()
        if self.func.sizes != func.sizes:
            raise TypeError('Parameter counts must agree')
        if not self.symbols:
            return [{}] if self.func.build_invariants() == func.build_invariants() else []
        equations = [a-b for left, right in zip(self.func.args, func.args) for a, b in zip(left, right)]
        return S.solve(equations, self.symbols, dict=True)


class FormulaCollection:
    def __init__(self):
        formulae = []; add_formulae(formulae)
        self.formulae = formulae

    def lookup_origin(self, func):
        S = _S(); z = S.Dummy('z')
        value = _S().hyperexpand(func(z))
        if not value.has(S.hyper):
            return Formula(func, z, value, [])
        best = None
        for formula in self.formulae:
            if formula.func.sizes != func.sizes:
                continue
            for substitutions in formula.find_instantiations(func):
                origin = Hyper_Function(*(tuple(a.subs(substitutions) for a in group) for group in formula.func.args))
                difficulty = origin.difficulty(func)
                if difficulty >= 0 and (best is None or difficulty < best[0]):
                    result = Formula(origin, formula.z, formula.closed_form.subs(substitutions), [])
                    best = difficulty, result
        return best[1] if best else None


class MeijerFormula(Formula):
    def __init__(self, an, ap, bm, bq, z, symbols, B, C, M, matcher):
        self.func, self.z = G_Function(an, ap, bm, bq), z
        self.symbols, self.B, self.C, self.M, self._matcher = list(symbols), B, C, M, matcher

    def try_instantiate(self, func):
        if self.func.signature != func.signature or self._matcher is None:
            return None
        matched = self._matcher(func)
        if matched is None:
            return None
        substitutions, target = matched
        return MeijerFormula(*target.args, self.z, [], self.B.subs(substitutions),
                             self.C.subs(substitutions), self.M.subs(substitutions), None)


class MeijerFormulaCollection:
    def __init__(self):
        formulae = []; add_meijerg_formulae(formulae)
        self.formulae = defaultdict(list)
        for formula in formulae:
            self.formulae[formula.func.signature].append(formula)

    def lookup_origin(self, func):
        for formula in self.formulae.get(func.signature, []):
            result = formula.try_instantiate(func)
            if result is not None:
                return result
        S = _S(); z = S.Dummy('z'); value = _S().hyperexpand(func(z))
        if not value.has(S.meijerg):
            generic = Formula(func, z, value, [])
            return MeijerFormula(*func.args, z, [], generic.B, generic.C, generic.M, None)
        return None


def add_formulae(formulae):
    S = _S(); z, a = S.Dummy('z'), S.Dummy('a')
    entries = [((), (), S.exp(z), []), ((a,), (), (1-z)**(-a), [a]),
               ((1, 1), (2,), -S.log(1-z)/z, []),
               ((S.Rational(1, 2),), (S.Rational(3, 2),), S.sqrt(S.pi)*S.erfi(S.sqrt(z))/(2*S.sqrt(z)), [])]
    for ap, bq, value, symbols in entries:
        formulae.append(Formula(Hyper_Function(ap, bq), z, value, symbols))


def add_meijerg_formulae(formulae):
    S = _S(); z, b = S.Dummy('z'), S.Dummy('b')
    def matcher(func):
        return ({b: func.bm[0]}, func)
    formulae.append(MeijerFormula([], [], [b], [], z, [b], S.Matrix([z**b*S.exp(-z)]),
                                 S.Matrix([[1]]), S.Matrix([[b-z]]), matcher))


def build_hypergeometric_formula(func):
    return Formula(func, _S().Dummy('z'), None, [])


def try_polynomial(func, z):
    S = _S()
    degrees = [-int(a) for a in func.ap if a.is_Integer and a <= 0]
    poles = [-int(b) for b in func.bq if b.is_Integer and b <= 0]
    if poles and (not degrees or min(degrees) > min(poles)):
        return S.oo
    if not degrees:
        return None
    degree = min(degrees)
    if degree > 10000:
        raise ValueError('Terminating hypergeometric series exceeds 10000 terms')
    term, result = S.S.One, S.S.One
    for k in range(degree):
        term = S.cancel(term*z*S.prod(a+k for a in func.ap)/((k+1)*S.prod(b+k for b in func.bq)))
        result += term
    return S.expand(result)


def try_shifted_sum(func, z):
    S = _S()
    uppers = [a for a in func.ap if a.is_Integer]
    lowers = [b for b in func.bq if b.is_Integer and b > 0]
    if len(uppers) != 1 or uppers[0] <= 0 or not lowers:
        return None
    r, k = int(uppers[0]), int(min(lowers))-1
    if max(r, k) > 256:
        raise ValueError('Shifted sum exceeds 256 terms')
    ap, bq = list(func.ap), list(func.bq)
    ap.remove(S.Integer(r)); bq.remove(S.Integer(k+1))
    ap, bq = [a-k for a in ap], [b-k for b in bq]
    factor = S.factorial(k)*S.prod(S.rf(b, k) for b in bq)/(z**k*S.prod(S.rf(a, k) for a in ap))
    ops = [ShiftA(a) for a in reversed(range(1, r))]+[MultOperator(factor)]
    correction = -sum((z**n*S.prod(S.rf(a, n) for a in ap)/(S.factorial(n)*S.prod(S.rf(b, n) for b in bq))
                       for n in range(k)), S.S.Zero)
    return Hyper_Function(ap, bq), ops, correction


def try_lerchphi(func):
    """Recognize the exact coefficient identity pFq(1,a;a+1;z)."""
    S = _S()
    if len(func.ap) == 2 and len(func.bq) == 1 and S.S.One in func.ap:
        a = func.ap[1] if func.ap[0] == 1 else func.ap[0]
        if S.simplify(func.bq[0]-a-1) == 0:
            z = S.Dummy('z')
            return Formula(func, z, a*S.lerchphi(z, 1, a), [])
    return None


def hyperexpand_special(ap, bq, z):
    value = _S().hyperexpand(_S().hyper(ap, bq, z))
    return value if not value.has(_S().hyper) else None


def debug(*args):
    from ..polys.polyconfig import query
    try:
        enabled = query('SYMOX_HYPEREXPAND_DEBUG')
    except (KeyError, ValueError):
        enabled = False
    if enabled:
        print(*args)
__all__ += ['Hyper_Function', 'G_Function', 'Operator', 'MultOperator', 'ShiftA', 'ShiftB',
    'UnShiftA', 'UnShiftB', 'MeijerShiftA', 'MeijerShiftB', 'MeijerShiftC', 'MeijerShiftD',
    'MeijerUnShiftA', 'MeijerUnShiftB', 'MeijerUnShiftC', 'MeijerUnShiftD', 'ReduceOrder',
    'reduce_order', 'reduce_order_meijer', 'make_simp', 'make_derivative_operator',
    'apply_operators', 'devise_plan', 'devise_plan_meijer', 'Formula', 'FormulaCollection',
    'MeijerFormula', 'MeijerFormulaCollection', 'add_formulae', 'add_meijerg_formulae',
    'build_hypergeometric_formula', 'try_polynomial', 'try_shifted_sum', 'try_lerchphi',
    'hyperexpand_special', 'debug']
