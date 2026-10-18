"""Scalar algebra over ordered matrix monomials."""


def linear_terms(expression, x=None):
    """Yield (scalar coefficient, ordered matrix monomial).

    No commutation of matrix factors is permitted. Variable-dependent matrix
    factors require entrywise algorithms rather than this linear view.
    """
    from . import MatMul, S
    for term in expression.args if expression.is_MatAdd else (expression,):
        if term.is_MatMul:
            scalar, matrices = term.as_coeff_matrices()
            matrix = MatMul(*matrices).doit()
        else:
            scalar, matrix = S.One, term
        if x is not None and x in matrix.free_symbols:
            raise NotImplementedError('matrix monomials must be independent of the scalar variable')
        yield scalar, matrix


def matrix_sum(terms, shape):
    from . import MatAdd, ZeroMatrix
    terms = [term for term in terms if not getattr(term, 'is_ZeroMatrix', False) and term != 0]
    return MatAdd(*terms).doit() if terms else ZeroMatrix(*shape)


def coefficient(expression, x, n):
    """A matrix Taylor coefficient, unlike the structural Expr.coeff view."""
    return matrix_sum([c.expand().coeff(x, n)*m for c, m in linear_terms(expression, x)], expression.shape)


def as_coeff_exponent(self, x):
    from . import S
    parts = [(c.as_coeff_exponent(x), matrix) for c, matrix in linear_terms(self, x)]
    exponents = {power for (_, power), _ in parts}
    if len(exponents) != 1:
        return self, S.Zero
    power = exponents.pop()
    return matrix_sum([c*m for (c, _), m in parts], self.shape), power


def as_poly(self, *gens, **args):
    from . import Poly
    from .polys.polyerrors import PolynomialError, GeneratorsNeeded
    if self.is_commutative is False:
        return None
    try:
        return Poly(self, *gens, **args)
    except (TypeError, PolynomialError, GeneratorsNeeded):
        # Poly has scalar coefficients; a matrix polynomial needs a distinct
        # coefficient backend. This is Expr.as_poly's failure convention.
        return None


def extract_branch_factor(self, allow_half=False):
    from . import S, MatMul
    if self.is_MatMul:
        scalar, matrices = self.as_coeff_matrices()
        remaining, power = scalar.extract_branch_factor(allow_half=allow_half)
        return MatMul(remaining, *matrices).doit(), power
    return self, S.Zero


def extract_multiplicatively(self, c):
    from . import MatrixExpr, MatMul, S, sympify
    if self == c:
        return S.One
    if not isinstance(c, MatrixExpr):
        c = sympify(c)
        scalar, matrices = self.as_coeff_matrices() if self.is_MatMul else (S.One, [self])
        rest = scalar.extract_multiplicatively(c)
        return None if rest is None else MatMul(rest, *matrices).doit()
    if not self.is_MatMul:
        return None
    target = c.args if c.is_MatMul else (c,)
    for index in range(len(self.args)-len(target)+1):
        if self.args[index:index+len(target)] == target:
            rest = self.args[:index]+self.args[index+len(target):]
            return MatMul(*rest).doit() if rest else S.One
    return None


def as_coefficient(self, expr):
    return self.extract_multiplicatively(expr)


def extract_additively(self, c):
    from . import MatrixExpr, S
    if c == 0:
        return self
    if self == c:
        return S.Zero
    if not isinstance(c, MatrixExpr) or self.shape != c.shape:
        return None
    terms = {}
    for scalar, matrix in linear_terms(self):
        terms[matrix] = terms.get(matrix, S.Zero)+scalar
    for scalar, matrix in linear_terms(c):
        if matrix not in terms:
            return None
        rest = terms[matrix].extract_additively(scalar)
        if rest is None:
            return None
        terms[matrix] = rest
    return matrix_sum([scalar*matrix for matrix, scalar in terms.items() if scalar != 0], self.shape)


def collect(self, syms, func=None, evaluate=True, exact=False, distribute_order_term=True):
    from . import MatrixExpr, collect as scalar_collect, S
    variables = tuple(syms) if isinstance(syms, (tuple, list, set)) else (syms,)
    if any(isinstance(v, MatrixExpr) for v in variables):
        raise NotImplementedError('matrix collect currently requires scalar generators')
    grouped = {}
    for scalar, matrix in linear_terms(self):
        values = scalar_collect(scalar, variables, evaluate=False, exact=exact,
                                distribute_order_term=distribute_order_term)
        for power, coefficient_ in values.items():
            part = coefficient_*matrix
            grouped[power] = grouped[power]+part if power in grouped else part
    if func is not None:
        grouped = {power: func(value) for power, value in grouped.items()}
    if not evaluate:
        return grouped
    return matrix_sum([power*value for power, value in grouped.items()], self.shape)


def leadterm(self, x, logx=None, cdir=0):
    from . import S
    terms, lowest = [], None
    for scalar, matrix in linear_terms(self, x):
        coefficient_, power = scalar.leadterm(x,logx=logx,cdir=cdir)
        if lowest is None or (power-lowest).is_negative is True:
            lowest, terms = power, [coefficient_*matrix]
        elif power == lowest:
            terms.append(coefficient_*matrix)
        elif (power-lowest).is_positive is not True:
            raise ValueError('cannot order symbolic matrix leading exponents')
    return matrix_sum(terms, self.shape), S.Zero if lowest is None else lowest


def as_leading_term(self, *symbols, logx=None, cdir=0):
    result = self
    for symbol in symbols:
        c, p = result.leadterm(symbol, logx=logx, cdir=cdir)
        result = c*symbol**p
    return result


def lseries(self, x=None, x0=0, dir='+', logx=None, cdir=0):
    from ._expr_analysis import _variable
    from . import Dummy, S
    x = _variable(self, x)
    small = Dummy('t', positive=True)
    change = x0+small if dir == '+' else x0-small
    translated = self.subs(x, change)
    finite_degree = 0
    for scalar, _ in linear_terms(translated, small):
        if scalar.is_polynomial(small) is not True:
            finite_degree = None
            break
        from . import Poly
        finite_degree = max(finite_degree, int(Poly(scalar, small).degree()) if scalar != 0 else 0)
    for n in range(257):
        c = coefficient(translated.series(small, 0, n+1).removeO().doit(), small, n)
        if not getattr(c, 'is_ZeroMatrix', False):
            yield c*(x-x0 if dir == '+' else x0-x)**n
        if finite_degree is not None and n >= finite_degree:
            return
    raise NotImplementedError('matrix lseries is limited to order 256')


def is_meromorphic(self, x, a):
    from . import Symbol
    if not isinstance(x, Symbol):
        raise TypeError('meromorphic variable must be a Symbol')
    try:
        parts = list(linear_terms(self, x))
    except NotImplementedError:
        return None
    results = [scalar.is_meromorphic(x, a) for scalar, _ in parts]
    return False if False in results else (None if None in results else True)


def is_hypergeometric(self, k):
    from . import hypersimp, cancel
    try:
        parts = list(linear_terms(self, k))
    except NotImplementedError:
        return None
    ratios = []
    for scalar, _ in parts:
        if scalar == 0:
            continue
        ratio = hypersimp(scalar, k)
        if ratio is None:
            return False
        ratios.append(ratio)
    return all(cancel(r-ratios[0]) == 0 for r in ratios[1:])


def invert(self, g, *gens, **args):
    from . import invert as scalar_invert, rem
    matrix = self.as_explicit()
    inverse_determinant = scalar_invert(matrix.det(), g, *gens, **args)
    return matrix.adjugate().applyfunc(lambda entry: rem(entry*inverse_determinant, g, *gens, **args)).as_immutable()


def round_(self, n=None):
    return self.as_explicit().applyfunc(lambda entry: entry.round(n)).as_immutable()


def rcall(self, *args):
    from ._matrixexpr_methods import _map_scalars
    return _map_scalars(self, lambda e: e.rcall(*args))


def replace(self, query, value, map=False, simultaneous=True, exact=None):
    from ._expression_replace import replace as replace_expression
    return replace_expression(self,query,value,map=map,simultaneous=simultaneous,exact=exact)


def install():
    from . import _core
    from ._operation_counts import count_ops
    from .series.matrix import matrix_fps, matrix_fourier_series
    for name, method in {
        'as_coeff_exponent': as_coeff_exponent, 'as_poly': as_poly,
        'extract_branch_factor': extract_branch_factor,
        'extract_multiplicatively': extract_multiplicatively,
        'as_coefficient': as_coefficient, 'extract_additively': extract_additively,
        'collect': collect, 'leadterm': leadterm, 'as_leading_term': as_leading_term,
        'lseries': lseries, 'is_meromorphic': is_meromorphic,
        'is_hypergeometric': is_hypergeometric, 'invert': invert, 'round': round_,
        'rcall': rcall, 'replace': replace,
        'count_ops': count_ops,
        'fps': matrix_fps, 'fourier_series': matrix_fourier_series,
    }.items():
        setattr(_core._MatrixExpr, name, method)
