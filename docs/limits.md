# Limits and asymptotic orders

[README](../README.md) · [Compatibility](compatibility.md) · [Series](series.md) · [Integration](integration.md)

`limit(e, z, z0, dir='+')`, `Limit` and `Order` (`O`) follow SymPy 1.14's
`sympy.series.limits`, `gruntz` and `order`. The limit engine, the Gruntz
algorithm and the `Order` class are native Rust
(`crates/symox-core/src/limits/`, `crates/symox-core/src/series/order.rs`);
a Python front end (`limit` in `python/symox/__init__.py`) adds the
heuristics that SymPy's `Limit.doit` applies around Gruntz, branch-cut
sides and oscillating limits. Sequence limits (`limit_seq`) are described
in [Series and expression analysis](series.md#pade-approximants-and-sequence-limits).

## The engine

The native `limit` works on one side at a time. An expression free of `z`
is returned unchanged and `z0 = nan` gives `nan`. For each side it tries:

1. a leading term from the native series engine (fast path);
2. otherwise the point is mapped to `t -> oo` (`z = z0 + 1/t` or
   `z0 - 1/t`, `z = -t` at `-oo`) and `Gruntz::limitinf` runs: the
   most-rapidly-varying subexpressions (mrv set, compared through the
   limits of `log f/log g`), rewriting in a single variable `w -> 0+`,
   the leading term `c0*w**e0` (`mrv_leadterm`) and the sign of `e0`; when
   `t` itself is in the mrv set the variable moves up with `t -> exp(t)`.
   Recursion is limited to depth 48 and results are cached. A `Piecewise`
   is replaced by its eventually true branch, and exact sign tests,
   `Abs` powers, irrational exponents and cancelling exponentials have
   dedicated steps (`gruntz/`, `exact_sign.rs`);
3. when Gruntz declines: bounded factors and oscillation
   (`bounded_affine`, which turns `sin`/`cos` of a growing argument and
   rounding remainders into accumulation bounds), values of `harmonic` and
   `subfactorial`, a squeeze for negative bases, `gamma` and `loggamma`
   at infinite arguments (`gamma_heuristic.rs`).

A `nan` computed from a finite expression counts as undecided, not as the
value. `O` terms take the limit of their expression, and an `O` at another
point is dropped.

The Python `limit` then adds, in this order: directed infinities such as
`I*oo` (substitute `z -> cdir*z`, then the limit at `oo`); a held `Limit`
when the expression already contains `oo`, `zoo` or `nan`; `zoo` for a base
tending to `-oo` raised to an exponent tending to `oo`; the branch-cut
rewrite below; the native engine; `arg` through a positive dummy; a
`powsimp` merge of powers of one base; local germs of `SingularityFunction`;
one-sided limits of `Piecewise`; and the oscillating walker. If nothing
applies, the result is the unevaluated `Limit`.

```python
from symox import symbols, limit, sin, cos, exp, log, sqrt, oo, E, I

x = symbols('x')
p = symbols('p', positive=True)
assert limit(sin(x)/x, x, 0) == 1
assert limit(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x, oo) == -1
assert limit(x**x, x, 0) == 1
assert limit(log(x)**3/x**p, x, oo) == 0
assert limit((1 + 1/x)**x, x, oo) == E
limit(cos(x)**(1/x**2), x, 0)     # exp(-1/2)
```

## Directions

`dir` is `'+'` (default), `'-'` or `'+-'`; anything else raises
`ValueError`. At `oo` the direction is normalized to `'-'` and at `-oo` to
`'+'`, so `Limit(1/x, x, oo).args[3]` is `-`. With `'+-'` both sides are
computed: equal values are returned, two infinite values give `zoo`, and
different values raise `ValueError`:

```python
from symox import floor, Abs
limit(1/x, x, 0)              # oo
limit(1/x, x, 0, '-')         # -oo
limit(1/x, x, 0, '+-')        # zoo
limit(1/x**2, x, 0, '+-')     # oo
limit(floor(x), x, 0, '-')    # -1
limit(Abs(x)/x, x, 0, '+-')
# ValueError: The limit does not exist since left hand limit = -1 and right hand limit = 1
```

A one-sided limit onto a branch cut is rewritten to the side it approaches
from (`branch_side_rewrite` in `python/symox/calculus/_branch_limit.py`,
checked numerically before use):

```python
limit(log(I*x - 1), x, 0, '-')    # -I*pi
limit(log(I*x - 1), x, 0, '+')    # I*pi
limit(sqrt(-1 + I*x), x, 0, '-')  # -I
```

## Oscillating limits

When the native engine leaves the limit undecided, `oscillating_limit`
(`python/symox/calculus/_accum_limit.py`) evaluates the leading terms of
oscillating subexpressions as `AccumBounds`, as SymPy does through
`sin(oo)` and similar values. The native `bounded_affine` route is tried
first; the rest of the walker is Python.

```python
from symox import AccumBounds, frac
assert limit(sin(x), x, oo) == AccumBounds(-1, 1)
assert limit(sin(1/x), x, 0) == AccumBounds(-1, 1)
assert limit(1/x + sin(x), x, oo) == AccumBounds(-1, 1)
assert limit(exp(sin(x)), x, oo) == AccumBounds(exp(-1), E)
assert limit(frac(x), x, oo) == AccumBounds(0, 1)
limit(x*sin(x), x, oo)        # AccumBounds(-oo, oo); SymPy: oo*sign(AccumBounds(-1, 1))
```

## `Limit` and `gruntz`

`Limit(e, z, z0, dir='+')` is a held native node with `.args`, `.var` and
`.dir`; `doit()` runs the same front end as `limit`. `gruntz(e, z, z0,
dir='+')` in `symox.series.gruntz` calls the native engine without the
other heuristics; when it declines, the oscillating walker gives the
accumulation bounds that SymPy's `mrv_leadterm` returns, and otherwise
`PoleError("The native limit engine could not determine this limit")` is
raised. A non-symbol variable raises `NotImplementedError("Second argument
must be a Symbol")`, as in SymPy.

```python
from symox import Limit
from symox.series.gruntz import gruntz
Limit(sin(x)/x, x, 0)              # Limit(sin(x)/x, x, 0, dir='+')
assert Limit(sin(x)/x, x, 0).doit() == 1
assert gruntz(exp(-x)*x**5, x, oo) == 0
assert gruntz(sin(x), x, oo) == AccumBounds(-1, 1)
```

## Errors and undecided results

- An undecided limit returns the unevaluated `Limit`, never a guessed value
  (`limit(f(x), x, 0)` for an undefined `f`, `limit(exp(x), x, I*oo)`).
- A `'+-'` limit whose sides differ raises `ValueError`.
- A point that depends on the variable raises SymPy's
  `NotImplementedError("Limits approaching a variable point are not
  supported (x -> x + 1)")` (`limit(x, x, x + 1)`).
- A non-symbol variable follows SymPy's `Limit.doit`: `limit(x**2, x**2, 0)`
  is 0, an expression free of the variable is returned unchanged
  (`limit(x, x**2, 0)` is `x`), and otherwise `TypeError("x**2 should be of
  symbol type")` is raised, also by `Limit(x**2 + 1, x**2, 0).doit()`.
- `dir` other than `'+'`, `'-'`, `'+-'` raises SymPy's `ValueError`
  (`direction must be one of '+', '-' or '+-', not ...`), and a `dir` that is
  neither a string nor a Symbol raises its `TypeError`.
- `gruntz` raises `PoleError` where `limit` would hold.

## `Order`

`O(expr, *args)` is SymPy's `Order`, implemented in Rust (`series/order.rs`
with `order_point.rs`, `order_diff.rs`, `order_log.rs`, `order_dominant.rs`).
The expression is reduced to its leading monomial at the point (default 0),
the variables are sorted, sums absorb terms contained in an order, and
products multiply orders. Points other than 0 are given as `(x, x0)`;
`O(nan)` is `nan`.

```python
from symox import O, symbols
y = symbols('y')
assert O(x) + O(x**2) == O(x)
assert O(x)*O(x**2) == O(x**3)
assert x + x**2 + O(x) == O(x)
assert O(sin(x)) == O(x)
O(1/x + 1/x**2, (x, oo))           # O(1/x, (x, oo))
O((x - 1)**2 + (x - 1), (x, 1))    # O(x - 1, (x, 1))
O(x*y + x**2*y)                    # O(x*y, x, y)
O(x) + O(y, (y, oo))               # O(y, (y, oo)) + O(x)
O(x).subs(x, x - 1)                # O(x - 1, (x, 1))
(1 + x + O(x**2)).removeO()        # x + 1
```

Containment (`O(a).contains(b)`, `b in O(a)`) follows SymPy's
`Order.contains` (`python/symox/series/_order_contains.py`), which decides
with the limit of the ratio and can return `None`. An undecided
containment keeps both terms in a sum instead of dropping one.

```python
assert O(x).contains(O(x**2)) is True
assert O(x**2).contains(O(x)) is False
assert x**2 in O(x)
assert O(x*y).contains(x**2*y**2) is True
```

Several variables at different points (`O(x*y, (x, 0), (y, oo))`) raise
`NotImplementedError` in the source, as in SymPy; the installed binary
of 2026-10-30 predates that check and builds the term.

## Limits and differences from SymPy

On the limit, Gruntz, series and Order statements of SymPy 1.14's test files,
the last recorded run (2026-10-30,
`audits/performance/2026-10-30-limgaps/README.md`) gave 1118 of 1293
statements equal to SymPy's output text, and 90 distinct native Gruntz
declines, 68 of them on series of `sin`, `cos`, `tan` or `floor` in the
mrv variable. Those declines go on to the oscillating front end; on that
replay 19 were decided by the native `bounded_affine`, 38 by the Python
walker and 29 stayed held.

Intentional deviations (validated):

- `limit(factorial(x)/x, x, -oo)` and `limit(gamma(x)/x, x, -oo)` stay
  held; SymPy returns 0.

Deviations pending validation:

- `limit(Piecewise((x, x < 0), (x + 1, True)), x, 0, '+-')` raises
  `ValueError` (the sides are 0 and 1); SymPy returns 1.
- `limit(atan(a*x)/x, x, oo)` is 0; SymPy holds it.
- `limit(subfactorial(x)/factorial(x), x, oo)` is held; SymPy returns 0
  (the value is `exp(-1)`).
- `limit(Abs(x + y)**2, x, 0)` is `y**2*sign(y)**2`; SymPy gives `y**2`.

Other differences:

- `limit(exp(x), x, I*oo)` is held; SymPy raises `NotImplementedError`.
- Series gaps reached through limits and `O`: the series of
  `floor(1/sin(p))` raises `PoleError`, `floor(1/p + p**7)` and
  `floor(1/x**2 + x**7)` are not expanded as in SymPy (the latter raises
  `ValueError`), `(x*frac(1/x)).series()` raises `PoleError` where SymPy
  raises `ValueError`, `O(x*frac(1/x))` is not reduced to `O(x)`, and
  `O(exp(1/x)*sin(1/x) + x)` keeps both terms. See `TODO.md` and
  [the parity notes](algorithmic-parity.md#limgaps-wave-36-limit-and-order-gaps-of-limitsrest-gruntz-declines-2026-10-30).
- The oscillating walker is Python; it is not part of the bare Rust library.
