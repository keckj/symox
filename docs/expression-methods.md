# Expression method options

[README](../README.md) · [Compatibility](compatibility.md)

All 97 measured `Expr` method names are exposed. Name presence does not imply
complete signature or behavioral parity; the following supported options have
dedicated implementation and regression tests.

## Coefficients and structural queries

`expr.as_coefficients_dict()` separates numeric coefficients from multiplicative
generators. It returns a `defaultdict(int)`, so missing generators have coefficient
zero. Supplying `*generators` treats factors independent of those generators as
coefficients, including symbolic parameters. Public arguments of held sums are
retained when extracting coefficients.

`coeff(x, n=1, right=False)` selects explicit multiplicative factors without
expansion or factoring. The residual coefficient may itself depend on `x`:
`(x*(x+1)).coeff(x)` returns `x+1`. A power `x**2` does not contain the explicit
factor `x`. Composite factors and exact numerical factors are accepted;
`coeff(1)` selects terms whose rational coefficient is one, and `coeff(-1)`
selects explicit negative signs. At power zero it selects independent terms,
with the standard special handling of additive generators. The power must be
an integer. `right` is accepted for commutative scalar factors; noncommutative
scalar symbols remain unsupported. Numeric coefficient powers above magnitude
100,000 raise `NotImplementedError` to avoid evaluating enormous integers.

```python
from symox import symbols, sin, Derivative, Function, Symbol

x, y, a = symbols('x y a')
coefficients = (3*a*x+x+2).as_coefficients_dict(x)
assert coefficients[x] == 3*a+1
assert coefficients[1] == 2
assert coefficients[y] == 0
assert sin(x).has(sin)
assert Derivative(Function('f')(x), x).has(Derivative)
```

`has` supports literal expressions, classes, built-in function heads, native
constructor heads and undefined function heads. Undefined functions retain
their argument-count and assumption metadata, even when their names collide
with built-in functions. Explicit foreign undefined-function selectors are
also recognized through the native conversion interface.

`find(query, group=False)` and `count(query)` accept scalar expression patterns,
classes, function heads and callable predicates. `group=True` returns a mapping
from matching values to their occurrence counts. These adapters share selectors
with matrix expression queries, including undefined-function metadata.

```python
f = Function('f')
expression = sin(x)+f(x)
assert expression.find(sin) == {sin(x)}
assert expression.find(Symbol, group=True) == {x: 2}
assert expression.count(f) == 1
```

Queries distinguish native node identity from mathematical equality when
walking bound-variable structures. Alpha-equivalent substitutions therefore
retain the names occurring in each original subtree. Repeated subgraphs have
their occurrence weights propagated without being unfolded. Callable predicates
should be pure: they are evaluated once per distinct graph node. Structural
tuple and condition metadata can be returned as immutable presentation views
whose entries remain native values. Queries are limited to 100,000 distinct
nodes.

## Wildcards and legacy matching

`Wild(name, exclude=(), properties=(), **assumptions)` creates a native wildcard
with its own type flag, constructor head, assumption facts and immutable exclusion
tuple, including their original order and repetitions. Conversion to/from explicitly
supplied SymPy objects and native pickle round trips retain this metadata.
`atoms(Wild)`, `find(Wild)` and `has(Wild)` select wildcards, while `Symbol`
selectors continue to include them. `symbols(..., cls=Wild)` uses the requested
constructor; custom callable constructors are also respected.

Predicate properties on wildcards remain unsupported. Nonempty `properties`
now raises `NotImplementedError` instead of silently accepting ineffective
constraints. Empty properties are supported. As in SymPy, `srepr` identifies
the wildcard and its assumptions but does not serialize its exclusion tuple;
use pickle or an explicit constructor to retain exclusions.

`match(pattern, old=False)` and `pattern.matches(value, repl_dict=None, old=False)`
accept the legacy matching flag. With `old=True`, a failed native match can solve
one wildcard that occurs linearly in the matching equation. The proposed binding
must satisfy its exclusions and reproduce the target expression symbolically.
This fallback does not search nonlinear roots or solve several wildcards at once.
Existing bindings passed to `matches` are copied, and matching still has the
native search limitations described under replacement below.

```python
from symox import Wild

a = Wild('a')
assert (2/x).match(a*x, old=True) == {a: 2/x**2}
```

## Operation counts

`count_ops(expr, visual=False)` and `expr.count_ops(visual=False)` share an
implementation. Numeric mode returns an integer. Visual mode returns a native
symbolic sum of labels such as `ADD`, `SUB`, `NEG`, `MUL`, `DIV`, `POW`, `SIN`
and `FUNC_F`. Native matrix products also expose matrix operation labels such
as `MATMUL`. Lists, tuples, dictionaries, dense/sparse matrices and arrays are
counted through their entries.

```python
from symox import count_ops, Symbol

expression = x*x+sin(x)
assert expression.count_ops() == 3
assert count_ops(expression, visual=True) == Symbol('ADD')+Symbol('POW')+Symbol('SIN')
```

Conditions are counted when supplied directly. As in SymPy's scalar counting
conventions, conditions nested inside a scalar expression can have different
labels. Counts describe symbolic operations and are not timing predictions.

Repeated graph nodes are processed once and tree-occurrence multiplicities are
propagated with arbitrary-size integers. This handles compact shared products
whose expanded operation count is exponential. At most 100,000 distinct graph
instructions are retained. Fraction normalization is limited to 10,000 unfolded
multiplicative factors and numeric integer powers with exponent magnitude at
most 100,000; unsupported normalization raises `NotImplementedError` rather
than evaluating a huge held expression.

## Numeric comparisons

Numeric atoms are compared at their stored exact values, including the binary
values of floats. Ordered relations and `Eq`/`Ne` use that comparison before any
rounded subtraction. SymPy 1.14 can round mixed float/rational operands in its
relational constructors; native comparisons keep their stored values distinct.
`Min` and `Max` retain the selected original value, so tiny exact rationals and
adjacent large integers do not collapse through machine rounding. Unresolved
nearby constant expressions can remain symbolic. Explicit nonreal arguments and
NaN are rejected by the public extrema constructors.

Non-strict sign bounds also simplify symbolic arguments. For a nonnegative `n`,
`Max(0, n)` returns `n` and `Min(0, n)` returns zero. Lattice absorption handles
forms such as `Min(x, Max(x, y))`, including nested argument subsets. These rules
also simplify clipped symbolic matrix dimensions.

Arbitrary-precision binary comparisons use sign, leading-bit positions and
bounded mantissa alignment. They do not allocate a denominator whose size is
proportional to the stored exponent just to compare against another value.

## Numerical recognition

`expr.nsimplify(constants=(), tolerance=None, full=False)` now forwards method
options to the native numerical recognition interface. A supplied constant
basis is supported, for example `value.nsimplify(constants=[pi])`. Exact inputs
retain their exact representation. Explicit tolerances and `full=True` for
irrational recognition remain unsupported.

The global `nsimplify(value, rational=True, rational_conversion='base10')`
converts float atoms to exact rationals, including symbolic coefficients and
exponents. Decimal conversion recognizes nearby small fractions at the displayed
precision and otherwise retains the decimal data. `rational_conversion='exact'`
preserves the stored binary value. Neither path rounds through a machine float,
so arbitrary-precision values and tiny nonzero numbers are retained.

```python
from symox import Float, nsimplify, Rational

assert nsimplify(Float('.1'), rational=True) == Rational(1, 10)
assert nsimplify(Float('.1'), rational=True, rational_conversion='exact') == Rational(3602879701896397, 36028797018963968)
assert nsimplify(Float('3.141592653589793'), rational=True, tolerance=Rational(1, 1000)) == Rational(355, 113)
```

With `rational=True`, a positive tolerance below one bounds the denominator by
`ceil(1/tolerance)`; a tolerance of at least one rounds on that tolerance grid
before bounding the denominator. These are rational approximation controls;
they do not enable the unfinished irrational identification search. Float
rationalization is limited to 100,000 precision bits and decimal exponent
magnitude 100,000. Unsupported sizes raise `NotImplementedError` before exact
denominator allocation. `rational=False` selects the existing native numerical
recognition path, which can return radicals and supplied constant multiples.

## Ordering, leading terms and assumptions

`nfloat(expr, n=15, exponent=False, dkeys=False)` converts numeric parts to
native floats at the requested decimal precision. Exact exponents are retained
at every nesting level by default; `exponent=True` converts them too. Undefined
function arguments stay exact, while symbolic built-in function arguments are
converted. Lambda signatures and asymptotic expansion points retain their
exact metadata.

```python
from symox import nfloat, Float, sqrt, sin

assert nfloat(sqrt(x)) == sqrt(x)
assert nfloat(sin(x/2)) == sin(Float('0.5')*x)
assert nfloat(x**2, exponent=True) == x**Float(2)
```

Lists, tuples, sets and dictionaries retain their container type. Dictionary
keys are preserved unless `dkeys=True`. Concrete matrices and arrays convert
entries, while symbolic matrices retain exact dimensions, indices and matrix
power exponents. Decimal precision must lie between 1 and 10,000. Symbolic
floating-conversion caches and structural traversal are limited to 100,000
distinct nodes; cyclic Python containers are rejected.

`as_ordered_terms(order=None, data=False)` accepts `lex`, `grlex`, `grevlex`,
their inverse `i...` orders, `rev-...` variants, and monomial-key callables.
It can return coefficient/monomial data together with the ordered generators.
The default presentation keeps SymPy's convention of putting the constant
first in `1-x`; explicit `order='lex'` follows the monomial order instead.
`as_ordered_factors`, recursive `sort_key` and `default_sort_key` propagate
explicit monomial orders to nested scalar expressions, matrix coefficients,
constructor tuples and algebraic-number arguments. Commutative factors are
sorted while matrix multiplication order is retained. Repeated native nodes
share their presentation keys; bound-variable names remain distinct in
alpha-equivalent expressions. Explicit recursive keys are limited to 128 nested
levels and 100,000 distinct nodes. The default uses the native cached ordering.

`leadterm(x, logx=None, cdir=0)` returns the coefficient and exponent of the
original variable on the selected approach ray. `as_leading_term(*symbols,
logx=None, cdir=0)` applies leading-term extraction successively to its variables;
without variables it returns the original expression. Finite complex directions
are supported, including the principal branches of logarithms and fractional
powers. `logx` substitutes a formal logarithm in the returned coefficient.

```python
from symox import exp, log, I

assert (exp(x)-1).leadterm(x, cdir=-1) == (1, 1)
assert log(x).leadterm(x, logx=y, cdir=I) == (y, 0)
assert ((x+x*x)*(y+y*y)).as_leading_term(x, y) == x*y
```

The native finite-center extension `leadterm(x, x0=point, dir='+')` is retained.
Its coefficient multiplies `(x-point)**exponent`. Essential singularities and
unsupported expansions raise an error; infinite centers use `aseries`.

`refine(assumption=...)` accepts the standard method keyword and delegates to
native assumption refinement. The previous spelling `assumptions=...` remains
an alias.

The `collect` method accepts coefficient callbacks and its evaluation, matching
and order-distribution controls positionally as well as by keyword, using the
same native collector as the global function.

## Replacement

`replace(query, value, map=False, simultaneous=True, exact=None)` supports native
expression patterns, classes, built-in and undefined function heads, and callable
predicates. Function-head replacements receive the matching node's arguments;
predicate replacements receive the node itself. Pattern callbacks receive the
wildcard bindings as named arguments.

`map=True` also returns the substitutions performed. With `exact=True`, matches
with a zero wildcard binding are skipped; this is the default for patterns
containing multiple distinct wildcards. `simultaneous=True` avoids applying the
rule twice when rebuilding a parent collapses it to an already processed child.
Returned replacement expressions are not recursively traversed again.

```python
from symox import Wild, cos

a = Wild('a')
result, changes = sin(x).replace(sin(a), cos(a), map=True)
assert result == cos(x)
assert changes == {sin(x): cos(x)}
assert (4*x).replace(2*x, y) == 4*x
```

Matching uses the native matcher; arbitrary algebraic rearrangements are not
fully equivalent to SymPy's matching search. In particular, symbolic fixed
factors must be present in multiplicative matches. Coefficient patterns should
give wildcards explicit generator exclusions when they represent independent
coefficients. Linked variable and shape metadata are validated when rebuilding
constructors. Compact shared graphs are processed once per distinct node, so
callbacks should be pure; the graph is limited to 100,000 distinct nodes.

## Substitution

Scalar expressions and concrete matrices accept `subs(old, new)` or one
mapping/iterable of pairs. Mappings, including `UserDict` and read-only mapping
proxies, replace larger patterns before their constituent symbols. Ordered
pair iterables preserve their sequential order. Each pair must contain exactly
two values; outer and inner iterators can be consumed once.

`simultaneous=True` applies all rules to the original expression. For example,
`(x/y).subs({x: 0, y: 0}, simultaneous=True)` is `nan`, whereas the sequential
substitution first produces zero. Matrix methods apply the same native scalar
binder, held-expression and dimension validation to every entry. Sparse matrices
retain sparse storage when implicit zeros stay zero. Empty matrices return an
empty result without visiting a scalar entry.

Substitution-priority size estimates memoize shared native graphs and inspect
held differentiation counts without expanding them.

## Differentiation

`diff(expression, *variables, evaluate=True)` and scalar
`expression.diff(*variables, evaluate=True)` share native differentiation
semantics. `evaluate=False` constructs a native held `Derivative`; its `doit()`
performs the calculation later. Orders may be Python or native integers, or
be supplied as `(variable, order)`, `[variable, order]` and native `Tuple` pairs.
Zero orders leave the expression unchanged. Negative, nonintegral and repeated
standalone order arguments are rejected before differentiation. Matrix elements,
undefined function applications and held derivatives are valid scalar variables.
The textual `sympify` parser accepts the same count/pair forms, variable inference
and evaluation flag. Orders outside the native 32-bit count representation are
rejected instead of truncating. Nested held derivatives preserve the total order
by retaining adjacent count groups when merging them would overflow that field.
Held derivative `free_symbols` come from the body; an unused differentiation
variable is not added to that set. Free-symbol, containment and substitution
identity queries read compact variable/count pairs without enumerating the
derivative order.

```python
held = sin(x).diff(x, 2, evaluate=False)
assert held == Derivative(sin(x), (x, 2))
assert held.doit() == -sin(x)
```

Without an explicit variable, scalar expressions and explicit matrices infer a
single free variable; multivariate inputs require a choice. Constant explicit
matrices differentiate to zero matrices. Matrix and array Jacobians retain the
finite-array implementation described in [arrays](arrays.md).

Held derivatives of whole matrices or arrays require a future lazy array
expression representation and raise `NotImplementedError`. Individual scalar
entries can be held explicitly with `applyfunc`. Evaluated derivatives use native
canonical simplification; `simplify=False` remains unsupported for evaluated
calls. It is accepted for held scalar derivatives, which perform no calculation.
