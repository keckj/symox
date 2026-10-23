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
    S = _S()
    return S.simplify(S.Mod(value, 1))


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
        return any(value.has(*symbols) for group in self.args for value in group)

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
        return not any(a.is_integer is True and a.is_nonpositive is True for a in self.ap+self.bq)


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


class Operator:
    def __init__(self, polynomial=0, generator=None):
        S = _S()
        self._d = generator or S.Dummy('Euler_operator')
        self._poly = S.Poly(polynomial, self._d)

    def apply(self, obj, op):
        coefficients = list(reversed(self._poly.all_coeffs()))
        value, result = obj, coefficients[0]*obj
        for coefficient in coefficients[1:]:
            value = op(value)
            result += coefficient*value
        return result


class MultOperator(Operator):
    def __init__(self, p):
        super().__init__(p)


class ShiftA(Operator):
    def __init__(self, ai):
        S = _S()
        ai = S.sympify(ai)
        if ai == 0:
            raise ValueError('Cannot increment a zero upper parameter')
        d = S.Dummy('Euler_operator')
        super().__init__(1+d/ai, d)


class ShiftB(Operator):
    def __init__(self, bi):
        S = _S()
        bi = S.sympify(bi)
        if bi == 1:
            raise ValueError('Cannot decrement a unit lower parameter')
        d = S.Dummy('Euler_operator')
        super().__init__(1+d/(bi-1), d)


def _inverse_affine(target, z, slope, intercept):
    S = _S()
    d = S.Dummy('Euler_operator')
    polynomial = S.Poly(_annihilator(target, z, d), d)
    root = -intercept/slope
    remainder = S.simplify(polynomial.eval(root))
    if remainder == 0:
        raise ValueError('Singular contiguous shift: parameters cancel')
    quotient = S.div(polynomial, S.Poly(slope*d+intercept, d))[0]
    return -quotient.as_expr()/remainder, d


class UnShiftA(Operator):
    def __init__(self, ap, bq, i, z):
        ap, bq = list(_params(ap)), list(_params(bq))
        ap[int(i)] -= 1
        ai = ap[int(i)]
        if ai == 0:
            raise ValueError('Cannot decrement a unit upper parameter')
        polynomial, d = _inverse_affine(Hyper_Function(ap, bq), z, 1/ai, 1)
        super().__init__(polynomial, d)


class UnShiftB(Operator):
    def __init__(self, ap, bq, i, z):
        ap, bq = list(_params(ap)), list(_params(bq))
        bq[int(i)] += 1
        bi = bq[int(i)]-1
        if bi == 0:
            raise ValueError('Cannot increment a zero lower parameter')
        polynomial, d = _inverse_affine(Hyper_Function(ap, bq), z, 1/bi, 1)
        super().__init__(polynomial, d)


class MeijerShiftA(Operator):
    def __init__(self, bi):
        S = _S(); d = S.Dummy('Euler_operator')
        super().__init__(S.sympify(bi)-d, d)


class MeijerShiftB(Operator):
    def __init__(self, ai):
        S = _S(); d = S.Dummy('Euler_operator')
        super().__init__(d+1-S.sympify(ai), d)


class MeijerShiftC(Operator):
    def __init__(self, bi):
        S = _S(); d = S.Dummy('Euler_operator')
        super().__init__(d-S.sympify(bi), d)


class MeijerShiftD(Operator):
    def __init__(self, ai):
        S = _S(); d = S.Dummy('Euler_operator')
        super().__init__(S.sympify(ai)-1-d, d)


class _MeijerInverse(Operator):
    group = 0
    change = 1

    def __init__(self, an, ap, bm, bq, i, z):
        groups = [list(_params(group)) for group in (an, ap, bm, bq)]
        groups[self.group][int(i)] += self.change
        parameter = groups[self.group][int(i)]
        # Invert the forward affine operator evaluated on the new function.
        slope, intercept = ((1, 1-parameter), (-1, parameter-1),
                            (-1, parameter), (1, -parameter))[self.group]
        polynomial, d = _inverse_affine(G_Function(*groups), z, _S().sympify(slope), intercept)
        super().__init__(polynomial, d)


class MeijerUnShiftA(_MeijerInverse):
    group, change = 2, -1


class MeijerUnShiftB(_MeijerInverse):
    group, change = 0, 1


class MeijerUnShiftC(_MeijerInverse):
    group, change = 3, -1


class MeijerUnShiftD(_MeijerInverse):
    group, change = 1, 1


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
        d = S.Dummy('Euler_operator')
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
        obj = object.__new__(cls); d = S.Dummy('Euler_operator')
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
    S = _S()
    pairs = [_pair_parameters(a, b) for a, b in zip(target.args, origin.args)]
    source = tuple(tuple(b for a, b in group) for group in pairs)
    destination = tuple(tuple(a for a, b in group) for group in pairs)
    if sum(abs(int(S.simplify(a-b))) for group in pairs for a, b in group) > 256:
        raise ValueError('Shift plan exceeds 256 steps')
    failed = set()
    def search(state):
        if state == destination:
            return []
        if state in failed:
            return None
        if len(failed) >= 10000:
            raise ValueError('Shift planning exceeds 10000 states')
        for group in range(2):
            for i, value in enumerate(state[group]):
                delta = int(S.simplify(destination[group][i]-value))
                if not delta:
                    continue
                try:
                    if group == 0:
                        operator = ShiftA(value) if delta > 0 else UnShiftA(*state, i, z)
                    else:
                        operator = UnShiftB(*state, i, z) if delta > 0 else ShiftB(value)
                except ValueError:
                    continue
                following = [list(values) for values in state]
                following[group][i] += 1 if delta > 0 else -1
                suffix = search(tuple(tuple(values) for values in following))
                if suffix is not None:
                    return [operator]+suffix
        failed.add(state)
        return None
    result = search(source)
    if result is None:
        raise ValueError('No nonsingular monotone contiguous shift plan exists')
    return list(reversed(result))


def devise_plan_meijer(fro, to, z):
    S = _S()
    if fro.signature != to.signature:
        raise ValueError('Meijer signatures must agree')
    current, operators = [list(group) for group in fro.args], []
    for group, (source, target) in enumerate(zip(fro.args, to.args)):
        for i, (a, b) in enumerate(zip(source, target)):
            delta = S.simplify(b-a)
            if not delta.is_Integer or abs(int(delta)) > 256:
                raise ValueError('Meijer shifts require integer steps of magnitude at most 256')
            for _ in range(abs(int(delta))):
                value = current[group][i]
                if delta > 0:
                    operator = (MeijerUnShiftB(*current, i, z) if group == 0 else
                                MeijerUnShiftD(*current, i, z) if group == 1 else
                                MeijerShiftA(value) if group == 2 else MeijerShiftC(value))
                else:
                    operator = (MeijerShiftB(value) if group == 0 else
                                MeijerShiftD(value) if group == 1 else
                                MeijerUnShiftA(*current, i, z) if group == 2 else MeijerUnShiftC(*current, i, z))
                operators.append(operator)
                current[group][i] += 1 if delta > 0 else -1
    return list(reversed(operators))


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
