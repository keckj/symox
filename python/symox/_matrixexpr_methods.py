"""Structural protocols for the native, noncommutative matrix DAG."""
from collections import Counter, defaultdict
from ._basic import Basic
from ._kind import matrix_tag

_native_has = None


class _NameAtom(str):
    """A structural matrix name; unlike a scalar Symbol it is never free."""
    args = ()
    is_Atom = True
    is_Number = is_Symbol = is_Add = is_Mul = is_Pow = is_Function = False

    def to_sympy(self):
        from sympy.core.symbol import Str
        return Str(self)


class _TupleView(tuple):
    """Basic-style view of native constructor tuples during inspection."""
    is_Tuple = True
    is_Atom = is_Number = is_Symbol = is_Add = is_Mul = is_Pow = is_Function = False

    @property
    def args(self):
        return tuple(self)

    @property
    def func(self):
        from . import Tuple
        return Tuple

    def to_sympy(self):
        from sympy import Tuple
        return Tuple(*(a.to_sympy() if hasattr(a, 'to_sympy') else a for a in self))


class ExprCondPair(_TupleView):
    """Presentation pair for the native flat Piecewise argument storage."""
    def __new__(cls, expr, cond):
        from . import true, false, Piecewise, piecewise_fold, ITE
        from .core.basic import as_Basic
        from .logic.boolalg import Boolean
        expr = as_Basic(expr)
        if cond == True:
            cond = true
        elif cond == False:
            cond = false
        elif getattr(cond, 'has', lambda *_: False)(Piecewise):
            cond = piecewise_fold(cond)
            if isinstance(cond, Piecewise):
                cond = cond.rewrite(ITE)
        if not isinstance(cond, Boolean):
            raise TypeError('Second argument must be a Boolean, not `%s`' % getattr(cond, 'class_name', type(cond).__name__))
        return tuple.__new__(cls,(expr,cond))

    @property
    def func(self):return type(self)

    @property
    def expr(self):return self[0]

    @property
    def cond(self):return self[1]

    @property
    def is_commutative(self):return self.expr.is_commutative

    free_symbols = Basic.free_symbols
    subs = Basic.subs
    xreplace = Basic.xreplace
    has = Basic.has
    atoms = Basic.atoms
    doit = Basic.doit

    def __eq__(self, other):
        return isinstance(other, ExprCondPair) and tuple.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(('ExprCondPair', self.expr, self.cond))

    def __getnewargs__(self):return self.expr, self.cond

    def to_sympy(self):
        from sympy.functions.elementary.piecewise import ExprCondPair
        return ExprCondPair(*(value.to_sympy() if hasattr(value,'to_sympy') else value for value in self))


_ExprCondPairView = ExprCondPair
Basic.register(ExprCondPair)
ExprCondPair.__module__ = 'symox.functions.elementary.piecewise'


def _children(value):
    if isinstance(value, (tuple, list)):
        return tuple(value)
    if getattr(value, '_symox_blockgrid', False):
        return tuple(value)
    if getattr(value,'is_Piecewise',False):
        arguments=value.args
        if arguments and isinstance(arguments[0],_ExprCondPairView):return arguments
        return tuple(_ExprCondPairView(*arguments[i:i+2]) for i in range(0,len(arguments),2))
    if getattr(value,'is_Order',False):
        arguments=value.args
        return (arguments[0],)+tuple(_TupleView(arguments[i:i+2]) for i in range(1,len(arguments),2))
    if matrix_tag(value) == 'MatrixSymbol':
        name, rows, cols = value.args
        return _NameAtom(name), rows, cols
    return getattr(value, 'args', ())


def _walk(value):
    # Occurrences are retained: count() counts repeated subexpressions.
    if isinstance(value, (tuple, list)):
        value = _TupleView(value)
    yield value
    for child in _children(value):
        yield from _walk(child)


def _query(value, query, matching=True):
    from ._structural_queries import matcher
    return matcher(query,matching=matching)(value)


def atoms(self, *types):
    types = tuple(t if isinstance(t, type) else type(t) for t in types)
    return {v for v in _walk(self)
            if (any(isinstance(v, t) for t in types) if types else not _children(v))}


def find(self, query, group=False):
    from ._structural_queries import find as find_in_graph
    return find_in_graph(self,query,group=group)


def count(self, query):
    from ._structural_queries import count as count_in_graph
    return count_in_graph(self,query)


def _fragment(value, pattern):
    from . import MatrixExpr
    if not isinstance(value, MatrixExpr) or not isinstance(pattern, MatrixExpr):
        return False
    if value._kind_tag != pattern._kind_tag:
        return False
    if value.is_MatMul:
        args = pattern.args
        return any(value.args[i:i+len(args)] == args for i in range(len(value.args)-len(args)+1))
    if value.is_MatAdd:
        return set(pattern.args).issubset(value.args)
    return False


def has(self, *patterns):
    from . import MatrixExpr, AlgebraicNumber
    from ._core import Expr
    from .combinatorics import Permutation
    from ._simptools import _head_name
    for pattern in patterns:
        predicate = isinstance(pattern, type) or _head_name(pattern) or (
            callable(pattern) and not isinstance(pattern, (Expr, AlgebraicNumber, Permutation)))
        if not predicate and isinstance(self, MatrixExpr) and not isinstance(pattern, (AlgebraicNumber, Permutation)):
            if _native_has(self, pattern):
                return True
            continue
        for value in _walk(self):
            if _query(value, pattern, matching=False):
                return True
            if (not predicate and isinstance(value, Expr) and
                    not isinstance(pattern, (AlgebraicNumber, Permutation)) and value.has(pattern)):
                return True
    return False


def has_free(self, *patterns):
    from ._core import Expr
    from .combinatorics import Permutation
    permutations = tuple(p for p in patterns if isinstance(p, Permutation))
    if permutations:
        # A permutation is an atomic constant and cannot be a bound variable.
        if any(_query(value, p, matching=False) for value in _walk(self) for p in permutations):
            return True
        patterns = tuple(p for p in patterns if not isinstance(p, Permutation))
        if not patterns:
            return False
    if any(_query(self, p, matching=False) or _fragment(self, p) for p in patterns):
        return True
    for child in _children(self):
        if isinstance(child, Expr):
            if child.has_free(*[p for p in patterns if not getattr(p, 'is_Matrix', False)]):
                return True
        elif hasattr(child, 'has_free'):
            if child.has_free(*patterns):
                return True
        elif any(_query(v, p, matching=False) for v in _walk(child) for p in patterns):
            return True
    return False


def has_xfree(self, s):
    if not isinstance(s, set):
        raise TypeError('has_xfree requires a set')
    from ._core import Expr
    if self in s:
        return True
    for child in _children(self):
        if isinstance(child, Expr):
            scalar = {p for p in s if not getattr(p, 'is_Matrix', False)}
            if child.has_xfree(scalar):
                return True
        elif hasattr(child, 'has_xfree'):
            if child.has_xfree(s):
                return True
        elif any(v in s for v in _walk(child)):
            return True
    return False


def copy(self):
    return self


def is_same(self, other, approx=None):
    from . import MatrixExpr
    from ._core import Expr
    def same(a, b):
        if isinstance(a, MatrixExpr) or isinstance(b, MatrixExpr):
            return (isinstance(a, MatrixExpr) and isinstance(b, MatrixExpr) and
                    a._kind_tag == b._kind_tag and len(a.args) == len(b.args) and
                    all(same(x, y) for x, y in zip(a.args, b.args)))
        if isinstance(a, Expr) and isinstance(b, Expr):
            return a.is_same(b, approx=approx) if approx is not None else a.is_same(b)
        return type(a) is type(b) and a == b
    return same(self, other)


def compare(self, other):
    if self == other:
        return 0
    a, b = self.sort_key(), other.sort_key()
    return (a > b) - (a < b)


def _map_scalars(value, transform, memo=None, preserve_indices=False):
    from . import MatrixExpr
    from ._core import Expr
    from . import _core
    if memo is None:
        memo = {}
    if isinstance(value, Expr):
        return transform(value)
    if isinstance(value, MatrixExpr):
        if value in memo:
            return memo[value]
        protected = set()
        if value._kind_tag == 'Inverse':
            protected.add(1)
        if preserve_indices:
            protected.update({
                'MatrixSymbol': (0, 1, 2), 'Identity': (0,), 'ZeroMatrix': (0, 1),
                'OneMatrix': (0, 1), 'FunctionMatrix': (0, 1), 'MatrixSlice': (1, 2),
                'MatPow': (1,), 'MatrixPermute': (2,),
            }.get(value._kind_tag, ()))
        args = tuple(a if i in protected else _map_scalars(a, transform, memo, preserve_indices)
                     for i, a in enumerate(value.args))
        result = value if args == value.args else value.func(*args)
        memo[value] = result
        return result
    if getattr(value, '_symox_blockgrid', False):
        return _core._BlockGrid(value.rows, value.cols,
                               [_map_scalars(a, transform, memo, preserve_indices) for a in value])
    if isinstance(value, (tuple, list)):
        parts = [_map_scalars(a, transform, memo, preserve_indices) for a in value]
        # Retain normalized slice metadata if no bound changed.
        return value if tuple(parts) == tuple(value) else type(value)(parts)
    if getattr(value, 'is_Matrix', False):
        return value.applyfunc(transform)
    return value


def as_dummy(self):
    return _map_scalars(self, lambda e: e if e.is_Atom else e.as_dummy())


def dummy_eq(self, other, symbol=None):
    from . import Dummy, Symbol, MatrixExpr
    if not isinstance(other, MatrixExpr):
        return False
    left, right = [], []
    for expression, found in [(self, left), (other, right)]:
        for value in _walk(expression):
            if isinstance(value, Dummy) and value not in found:
                found.append(value)
    if symbol is not None and not right and symbol in other.free_symbols:
        right.append(symbol)
    if len(left) != len(right):
        return False
    # Use fresh shared dummies, avoiding collisions with actual scalar symbols.
    replacements = [Dummy('d') for _ in left]
    return self.xreplace(dict(zip(left, replacements))) == other.xreplace(dict(zip(right, replacements)))


def matches(self, expr, repl_dict=None, old=False):
    from . import MatrixExpr
    from ._core import Expr
    if getattr(self, 'is_MatMul', False) and getattr(expr, 'is_MatMul', False):
        # The reference MatMul inherits Mul.matches: the commutative factors
        # match first (extending repl_dict), then the ordered matrix factors.
        from ._matching import _mul
        return _mul(self, expr, repl_dict, old)
    bindings = dict(repl_dict or {})
    budget = 10000
    def visit(pattern, value):
        nonlocal bindings, budget
        budget -= 1
        if budget < 0:
            raise NotImplementedError('matrix pattern matching exceeded its search budget')
        if isinstance(pattern, MatrixExpr):
            if not (isinstance(value, MatrixExpr) and pattern._kind_tag == value._kind_tag and
                    len(pattern.args) == len(value.args)):
                return False
            if pattern.is_MatAdd:
                def assign(patterns, values):
                    nonlocal bindings
                    if not patterns:
                        return True
                    before = dict(bindings)
                    for i, candidate in enumerate(values):
                        bindings = dict(before)
                        if visit(patterns[0], candidate) and assign(patterns[1:], values[:i]+values[i+1:]):
                            return True
                    bindings = before
                    return False
                return assign(pattern.args, value.args)
            return all(visit(a, b) for a, b in zip(pattern.args, value.args))
        if isinstance(pattern, Expr) and isinstance(value, Expr):
            result = pattern.xreplace(bindings).matches(value,old=old)
            if result is None:
                return False
            bindings.update(result)
            return True
        if isinstance(pattern, (tuple, list)):
            return (isinstance(value, (tuple, list)) and len(pattern) == len(value) and
                    all(visit(a, b) for a, b in zip(pattern, value)))
        return pattern == value
    return bindings if visit(self, expr) else None


def match(self, pattern, old=False):
    if getattr(pattern, 'is_MatrixExpr', False):
        return pattern.matches(self, old=old)
    if getattr(pattern, 'is_Wild', False):
        # Wild.matches: any expression free of the excluded atoms and
        # satisfying the properties binds to the wildcard.
        # (an excluded matrix symbol is stored as a native scalar node)
        free = self.free_symbols
        if any(x == self or self.has(x) or any(x == f for f in free)
               for x in getattr(pattern, 'exclude', ())):
            return None
        if not all(f(self) for f in getattr(pattern, 'properties', ())):
            return None
        return {pattern: self}
    return None


def _identity(self, *args, **kwargs):
    return self


def getO(self):
    orders = [v for v in _walk(self) if getattr(v, 'is_Order', False)]
    return None if not orders else sum(orders[1:], orders[0])


def getn(self):
    order = getO(self)
    return None if order is None else order.getn()


def removeO(self):
    return _map_scalars(self, lambda e: e.removeO())


def _scalar_method(name):
    def method(self, *args, **kwargs):
        if name == 'rewrite' and len(args) == 1:
            from ._matrix_permutations import rewrite_permutation
            result = rewrite_permutation(self, args[-1], deep=kwargs.get('deep', True))
            if result is not None:
                return result
        return _map_scalars(self, lambda e: getattr(e, name)(*args, **kwargs),
                            preserve_indices=name in ('evalf', 'n'))
    method.__name__ = name
    return method


def simplify(self, **kwargs):
    return _map_scalars(self, lambda e: e.simplify(**kwargs)).doit()


def as_coeff_add(self, *deps):
    from . import S, MatAdd
    if deps and not self.has_free(*deps):
        return self, ()
    if self.is_MatAdd:
        independent = [a for a in self.args if deps and not a.has_free(*deps)]
        dependent = tuple(a for a in self.args if a not in independent)
        return (MatAdd(*independent).doit() if independent else S.Zero), dependent
    return S.Zero, (self,)


def as_numer_denom(self):
    from . import S
    return self, S.One


def as_powers_dict(self):
    from . import S
    if self.is_MatMul:
        result = defaultdict(int)
        for factor in self.args:
            for base, exponent in factor.as_powers_dict().items():
                result[base] += exponent
        return result
    return defaultdict(int, {self: S.One})


def as_ordered_factors(self, order=None):
    if order is not None and self.is_MatMul:
        commutative,noncommutative=self.args_cnc()
        return sorted(commutative,key=lambda factor:factor.sort_key(order=order))+noncommutative
    return list(self.args) if self.is_MatMul else [self]


def args_cnc(self, cset=False, warn=True, split_1=True):
    from . import S
    commutative, noncommutative = [], []
    if self.is_MatMul:
        for a in self.args:
            (noncommutative if getattr(a, 'is_Matrix', False) else commutative).append(a)
    else:
        noncommutative.append(self)
    if split_1 and commutative and commutative[0].is_negative:
        commutative[:1] = [S.NegativeOne, -commutative[0]]
    return [set(commutative) if cset else commutative, noncommutative]


def as_real_imag(self, deep=True, **hints):
    from . import I, S
    conjugate = self.conjugate()
    if conjugate == self:
        return self, S.Zero
    return (self+conjugate)/2, (self-conjugate)/(2*I)


def primitive(self):
    from . import S
    return S.One, self


def as_content_primitive(self, radical=False, clear=True):
    from . import S, MatMul, MatAdd
    if self.is_MatMul:
        scalar, matrices = self.as_coeff_matrices()
        content, rest = scalar.as_content_primitive(radical=radical, clear=clear)
        return content, MatMul(rest, *matrices).doit()
    if self.is_MatAdd:
        parts = [a.as_content_primitive(radical=radical, clear=clear) for a in self.args]
        contents = [c for c, _ in parts]
        from . import Rational
        from math import gcd, lcm
        numerator, denominator = 0, 1
        for c in contents:
            p, q = c.as_numer_denom()
            numerator, denominator = gcd(numerator, int(p)), lcm(denominator, int(q))
        content = Rational(numerator, denominator)
        if not clear and any(c.is_Integer for c in contents):
            content = S.One
        return content, MatAdd(*(c/content*p for c, p in parts)).doit()
    return S.One, self


def could_extract_minus_sign(self):
    if self.is_MatMul:
        return self.as_coeff_matrices()[0].could_extract_minus_sign()
    if self.is_MatAdd:
        negative = sum(a.could_extract_minus_sign() for a in self.args)
        return negative > len(self.args)/2
    return False


def valid_index(self, i, j):
    from . import sympify
    from ._core import Expr
    try:
        i, j = sympify(i), sympify(j)
    except (TypeError, ValueError):
        return False
    if not isinstance(i, Expr) or not isinstance(j, Expr):
        return False
    return all(getattr(test, 'is_BooleanFalse', False) is False and test != False
               for test in (i >= -self.rows, i < self.rows, j >= -self.cols, j < self.cols))


def _zero(value):
    return value == 0 or getattr(value, 'is_ZeroMatrix', False)


def equals(self, other, failing_expression=False):
    from . import MatrixExpr, DenseMatrix as Matrix
    if not isinstance(other, MatrixExpr):
        if isinstance(other, Matrix):
            other = MatrixExpr(other)
        else:
            return False
    if self == other:
        return True
    if self.shape != other.shape:
        return False
    difference = (self-other).doit()
    if _zero(difference):
        return True
    try:
        rows, cols = int(self.rows), int(self.cols)
    except (TypeError, ValueError):
        return difference if failing_expression else None
    if rows*cols > 256:
        return difference if failing_expression else None
    uncertain = None
    for row in range(rows):
        for col in range(cols):
            entry = difference[row, col]
            result = entry.equals(0)
            if result is False:
                return False
            if result is None:
                uncertain = entry
    if uncertain is None:
        return True
    return uncertain if failing_expression else None


def as_coefficients_dict(self, *syms):
    from . import S
    result = defaultdict(int)
    terms = self.args if self.is_MatAdd else (self,)
    for term in terms:
        if syms and not term.has_free(*syms):
            result[S.One] += term
        else:
            coefficient, monomial = term.as_coeff_Mul()
            result[monomial] += coefficient
    return result


def as_independent(self, *deps, **hint):
    from . import S
    unknown = set(hint).difference({'as_Add'})
    if unknown:
        raise TypeError('unsupported as_independent hints: '+', '.join(sorted(unknown)))
    additive = hint.get('as_Add', self.is_MatAdd)
    identity = S.Zero if additive else S.One
    return (identity, self) if self.has_free(*deps) else (self, identity)


def as_terms(self):
    from . import S
    from ._simptools import _decompose_power
    terms, all_generators = [], set()
    for term in (self.args if self.is_MatAdd else (self,)):
        coefficient, powers, nc = complex(1), defaultdict(int), []
        for factor in (term.args if term.is_MatMul else (term,)):
            if getattr(factor, 'is_Matrix', False):
                nc.append(factor)
            elif factor.is_number:
                coefficient *= complex(factor)
            else:
                base, exponent = _decompose_power(factor)
                powers[base] += exponent
                all_generators.add(base)
        terms.append((term, (coefficient, powers, tuple(nc))))
    generators = sorted(all_generators, key=lambda e: e.sort_key())
    return [(term, ((coefficient.real, coefficient.imag),
                    tuple(powers.get(g, S.Zero) for g in generators), nc))
            for term, (coefficient, powers, nc) in terms], generators


def as_ordered_terms(self, order=None, data=False):
    if order is not None:
        from ._expr_methods import as_ordered_terms as ordered
        return ordered(self,order=order,data=data)
    terms, generators = self.as_terms()
    terms.sort(key=lambda item: (tuple(-e for e in item[1][1]),
                                  tuple(a.sort_key() for a in item[1][2])))
    return (terms, generators) if data else [term for term, _ in terms]


def coeff(self, x, n=1, right=False, _first=True):
    import operator
    from . import S, MatAdd, MatMul, MatrixExpr, sympify
    n = operator.index(n)
    if not isinstance(x, MatrixExpr):
        x = sympify(x)
    if n == 0:
        independent = [a for a in (self.args if self.is_MatAdd else (self,)) if not a.has(x)]
        return MatAdd(*independent).doit() if independent else S.Zero
    target = x**n if n != 1 else x
    if self == target:
        return S.One
    if self.is_MatAdd:
        parts = [term.coeff(x, n, right=right) for term in self.args]
        parts = [p for p in parts if not _zero(p)]
        if not parts:
            return S.Zero
        result = parts[0]
        for part in parts[1:]:
            result += part
        return result
    if self.is_MatMul and isinstance(target, MatrixExpr):
        wanted = target.args if target.is_MatMul else (target,)
        positions = [i for i in range(len(self.args)-len(wanted)+1)
                     if self.args[i:i+len(wanted)] == wanted]
        if positions:
            position = positions[0] if right else positions[-1]
            tail = self.args[position+len(wanted):] if right else self.args[:position]
            return MatMul(*tail).doit() if tail else S.One
    return S.Zero


def _analytic_method(name):
    def method(self, *args, **kwargs):
        result = _map_scalars(self, lambda e: getattr(e, name)(*args, **kwargs), preserve_indices=True)
        return result.doit(deep=False)
    method.__name__ = name
    return method


def taylor_term(self, n, x, *previous_terms):
    import operator
    from . import factorial
    n = operator.index(n)
    if n < 0:
        raise ValueError('Taylor term order must be nonnegative')
    return self.diff(x, n).subs(x, 0) * (x**n/factorial(n))


def dir_(self, x, cdir):
    current = self
    for order in range(1, 257):
        current = current.diff(x)
        coefficient = current.subs(x, 0)
        if not _zero(coefficient) or _zero(current):
            return coefficient*cdir**order
    raise NotImplementedError('matrix direction requires more than 256 derivatives')


def integrate(self, *args, **kwargs):
    from . import integrate as scalar_integrate, MatAdd, MatMul, FunctionMatrix, Lambda, Symbol
    from ._core import Expr
    if not args:
        candidates = {v for v in self.free_symbols if isinstance(v, Symbol)}
        if len(candidates) != 1:
            raise ValueError('specify an integration variable for a matrix expression')
        args = (next(iter(candidates)),)
    variables = [a[0] if isinstance(a, (tuple, list)) else a for a in args]
    def independent(matrix):
        # Matrix symbols are constant entries. A varying shape does not
        # define a matrix-valued integral on a fixed vector space.
        return not any(variable in matrix.free_symbols for variable in variables)
    if self.is_MatAdd:
        return MatAdd(*(term.integrate(*args, **kwargs) for term in self.args)).doit()
    if self.is_FunctionMatrix:
        function = self.lamda
        if any(v in function.variables for v in variables):
            raise ValueError('cannot integrate a FunctionMatrix over a bound row or column index')
        body = scalar_integrate(function.expr, *args, **kwargs)
        return FunctionMatrix(self.rows, self.cols, Lambda(function.variables, body))
    if self.is_MatMul:
        coefficient, matrices = self.as_coeff_matrices()
        if all(independent(a) for a in matrices):
            return MatMul(scalar_integrate(coefficient, *args, **kwargs), *matrices).doit()
    elif independent(self):
        from . import S
        return scalar_integrate(S.One, *args, **kwargs)*self
    try:
        rows, cols = int(self.rows), int(self.cols)
    except (TypeError, ValueError):
        raise NotImplementedError('integration of dependent matrix factors requires concrete dimensions') from None
    if rows*cols > 256:
        raise NotImplementedError('entrywise matrix integration is limited to 256 entries')
    return self.as_explicit().applyfunc(lambda e: scalar_integrate(e, *args, **kwargs))


def _algebraic_query(name):
    def method(self, *symbols):
        from . import Symbol
        if self.is_MatrixSymbol:
            return True
        symbols = symbols or tuple(self.free_symbols)
        if not set(symbols).intersection(self.free_symbols):
            return True
        if self.is_MatAdd or self.is_MatMul:
            for argument in self.args:
                selected = symbols if getattr(argument, 'is_Matrix', False) else tuple(s for s in symbols if isinstance(s, Symbol))
                if not selected:
                    continue
                result = getattr(argument, name)(*selected)
                if result is not True:
                    return False
            return True
        return None
    method.__name__ = name
    return method


def is_constant(self, *wrt, **flags):
    from . import Symbol
    if self.is_MatrixSymbol:
        return True
    wrt = tuple(wrt) or tuple(v for v in self.free_symbols if isinstance(v, Symbol))
    if not wrt or not self.free_symbols.intersection(wrt):
        return True
    if self.is_MatAdd or self.is_MatMul:
        results = [a.is_constant(*wrt, **flags) for a in self.args]
        return False if False in results else (None if None in results else True)
    if self.is_FunctionMatrix:
        return self.lamda.expr.is_constant(*wrt, **flags)
    if self.is_MatPow or self._kind_tag in ('Transpose', 'Adjoint', 'Inverse', 'MatrixConjugate'):
        return all(a.is_constant(*wrt, **flags) for a in self.args)
    return None


def _diagonal_length(self):
    if self._kind_tag=='DiagonalOf':
        return self.shape[0]
    if self._kind_tag!='DiagonalMatrix':
        raise AttributeError('diagonal_length is available for diagonal matrix expressions')
    rows,cols=self.shape
    if rows.is_Integer and cols.is_Integer:
        return min(rows,cols)
    if rows.is_Integer:return rows
    if cols.is_Integer:return cols
    if rows==cols:return rows
    try:return min(rows,cols)
    except TypeError:return None


def _block_quadrants(self, message):
    from . import ShapeError
    if getattr(self, 'blockshape', None) != (2, 2) or matrix_tag(self) not in ('BlockMatrix', 'BlockDiagMatrix'):
        raise ShapeError(message)
    grid = self.blocks
    return grid[0, 0], grid[0, 1], grid[1, 0], grid[1, 1]


_native_schur = None


def schur(self, mat='A', generalized=False):
    """SymPy's ``BlockMatrix.schur``; the generalized complement uses the
    Moore-Penrose inverse ``(X.T*X)**(-1)*X.T`` of the chosen block."""
    from . import NonInvertibleMatrixError
    A, B, C, D = _block_quadrants(self, 'Schur Complement can only be calculated for 2x2 block matrices')
    d = {'A': A, 'B': B, 'C': C, 'D': D}
    if mat not in d or (not generalized and mat in 'AD'):
        return _native_schur(self, mat)
    try:
        inv = (d[mat].T*d[mat]).inv()*d[mat].T if generalized else d[mat].inv()
    except NonInvertibleMatrixError:
        raise NonInvertibleMatrixError('The given matrix is not invertible. Please set generalized=True \
            to compute the generalized Schur Complement which uses Moore-Penrose Inverse')
    if mat == 'A':
        return D - C * inv * B
    if mat == 'B':
        return C - D * inv * A
    if mat == 'C':
        return B - A * inv * D
    return A - B * inv * C


def LUdecomposition(self):
    """SymPy's block LU decomposition ``(L, U)`` of a 2x2 block matrix."""
    from . import S, ZeroMatrix, BlockMatrix
    from . import NonInvertibleMatrixError
    A, B, C, D = _block_quadrants(self, 'Block LU decomposition is supported only for 2x2 block matrices')
    try:
        A = A**S.Half
        AI = A.I
    except NonInvertibleMatrixError:
        raise NonInvertibleMatrixError('Block LU decomposition cannot be calculated when\
                    "A" is singular')
    Z = ZeroMatrix(*B.shape)
    Q = self.schur()**S.Half
    L = BlockMatrix([[A, Z], [C*AI, Q]])
    U = BlockMatrix([[A, AI*B], [Z.T, Q]])
    return L, U


def _is_kronecker(value):
    return isinstance(value, _core_matrix_expr()) and matrix_tag(value) == 'KroneckerProduct'


def _core_matrix_expr():
    from . import _core
    return _core._MatrixExpr


def has_matching_shape(self, other):
    """Whether multiplication can be brought inside two Kronecker products."""
    return (_is_kronecker(self) and _is_kronecker(other)
            and self.cols == other.rows
            and len(self.args) == len(other.args)
            and all(a.cols == b.rows for (a, b) in zip(self.args, other.args)))


_native_structurally_equal = None


def structurally_equal(self, other):
    if _is_kronecker(self):
        return (_is_kronecker(other)
                and self.shape == other.shape
                and len(self.args) == len(other.args)
                and all(a.shape == b.shape for (a, b) in zip(self.args, other.args)))
    return _native_structurally_equal(self, other)


_native_inverse = None


def _inverse(self):
    # SymPy's MatPow._eval_inverse: the inverse of B**e is B**(-e), so that
    # A**(1/2)*(A**(1/2)).I cancels.
    if matrix_tag(self) == 'MatPow' and self.is_square is not False:
        from . import MatPow
        base, exp = self.args
        return MatPow(base, -exp)
    return _native_inverse(self)


_native_str = None
_printing = set()


def _str(self):
    # SymPy's ``Printable.__str__`` is ``sstr``: the Python printer orders
    # the terms of a MatAdd as SymPy does (``A - B*D**(-1)*C``).  The native
    # string is the fallback for a node the Python printer would print
    # through ``repr`` of itself.
    key = id(self)
    if key in _printing:
        return _native_str(self)
    _printing.add(key)
    try:
        from .printing.str import sstr
        return sstr(self)
    except Exception:
        return _native_str(self)
    finally:
        _printing.discard(key)


def install():
    from . import _core
    global _native_schur, _native_structurally_equal, _native_inverse, _native_str
    if _native_str is None:
        _native_str = _core._MatrixExpr.__str__
        _core._MatrixExpr.__str__ = _str
        _core._MatrixExpr.__repr__ = _str
    if _native_schur is None:
        _native_schur = _core._MatrixExpr.schur
        _native_structurally_equal = _core._MatrixExpr.structurally_equal
        _native_inverse = _core._MatrixExpr.inverse
        _core._MatrixExpr.I = property(_inverse)
        _core._MatrixExpr.inv = _inverse
        _core._MatrixExpr.inverse = _inverse
    global _native_has
    if _native_has is None:
        _native_has = _core._MatrixExpr.has
    _core._MatrixExpr.diagonal_length=property(_diagonal_length)
    methods = {
        'atoms': atoms, 'find': find, 'count': count, 'has': has,
        'has_free': has_free, 'has_xfree': has_xfree, 'copy': copy,
        'is_same': is_same, 'compare': compare, 'as_dummy': as_dummy,
        'dummy_eq': dummy_eq, 'matches': matches, 'match': match,
        'canonicalize': _identity, 'as_expr': _identity,
        'getO': getO, 'getn': getn, 'removeO': removeO, 'simplify': simplify,
        'as_coeff_add': as_coeff_add, 'as_numer_denom': as_numer_denom,
        'as_powers_dict': as_powers_dict, 'as_ordered_factors': as_ordered_factors,
        'args_cnc': args_cnc, 'as_real_imag': as_real_imag,
        'primitive': primitive, 'as_content_primitive': as_content_primitive,
        'could_extract_minus_sign': could_extract_minus_sign, 'valid_index': valid_index,
        'equals': equals, 'as_coefficients_dict': as_coefficients_dict,
        'as_independent': as_independent, 'as_terms': as_terms,
        'as_ordered_terms': as_ordered_terms, 'coeff': coeff, 'normal': _identity,
        'taylor_term': taylor_term, 'dir': dir_, 'integrate': integrate,
        'is_constant': is_constant,
        'schur': schur, 'LUdecomposition': LUdecomposition,
        'has_matching_shape': has_matching_shape, 'structurally_equal': structurally_equal,
    }
    for name in ('evalf', 'n', 'nsimplify', 'trigsimp', 'powsimp', 'combsimp',
                 'gammasimp', 'radsimp', 'ratsimp', 'cancel', 'factor', 'together',
                 'refine', 'rewrite', 'apart', 'separate'):
        methods[name] = _scalar_method(name)
    for name in ('series', 'nseries', 'aseries', 'limit'):
        methods[name] = _analytic_method(name)
    for name in ('is_polynomial', 'is_rational_function', 'is_algebraic_expr'):
        methods[name] = _algebraic_query(name)
    for name, method in methods.items():
        setattr(_core._MatrixExpr, name, method)
    from ._matrix_index_sums import from_index_summation
    _core._MatrixExpr.from_index_summation = staticmethod(from_index_summation)
    class FromIter:
        def __get__(self, instance, owner):
            from . import MatrixExpr
            constructor = instance.func if instance is not None else MatrixExpr
            return constructor.fromiter
    _core._MatrixExpr.fromiter = FromIter()
