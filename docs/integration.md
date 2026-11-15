# Integration

[README](../README.md) · [Compatibility](compatibility.md) · [Limits](limits.md)

`integrate(f, *limits, **hints)` and `Integral(f, *limits)` follow SymPy
1.14's strategy order, hints and result classes. The rule-based stages, the
Meijer G method and the rational, Risch and PRDE kernels run in Rust
(`crates/symox-core/src/integrate*`); the other reference drivers are Python
ports in `python/symox/integrals/`. SymPy is not imported at runtime.

## Indefinite integration pipeline

`integrate(f, x)` first applies `Integral.doit`'s rewrites
(`piecewise::doit_rewrite`): `Abs` and `sign` become `Piecewise` when `x` is
real (or the bounds are finite and real), `Min` and `Max` become
`Piecewise`, and a nested `Piecewise` is folded to the top and integrated
branch by branch. Any other integrand goes through these stages in order;
the first one that returns an antiderivative wins:

1. Native rules (`integ` in `integrate.rs`), in the style of SymPy's
   `manualintegrate`: linearity, powers, a table of elementary and special
   functions with linear arguments, rational functions with a denominator
   of degree at most two (or a power of a linear or quadratic factor, or
   pieces of an `apart` decomposition), powers of `sin`/`cos`, integration
   by parts in LIATE order, cyclic `exp`/trigonometric products,
   u-substitution and `DiracDelta` rules. When this result contains an
   `atan` of an `exp`/`log` tower, the Risch result is preferred.
2. `risch_integrate(f, x, separate_integral=True)`. A part proven
   non-elementary is integrated again without Risch, and stays a
   `NonElementaryIntegral` if that fails.
3. A term loop over the summands, each with its `x`-free coefficient taken
   out (`as_independent`): SymPy's `_eval_integral` stages `ratint`,
   `trigintegrate`, `deltaintegrate`, `singularityintegrate`,
   `risch_integrate` and `heurisch`; then `meijerint_indefinite`; then
   `manualintegrate`; then the native rules again. A single term that fails
   is retried once after rewriting products of `sin`/`cos` as sums. In the
   bare Rust library the native rules replace the Python stages.
4. Otherwise the result is an unevaluated `Integral`, with the `x`-free
   factors moved in front of it.

```python
from symox import symbols, integrate, exp, sin, cos, log, sqrt, erf, pi

x, n = symbols('x n')
integrate(x/(x**3 - 1), x)
# log(x - 1)/3 - log(x**2 + x + 1)/6 + sqrt(3)*atan(2*sqrt(3)*x/3 + sqrt(3)/3)/3
integrate(x**n, x)       # Piecewise((x**(n + 1)/(n + 1), Ne(n, -1)), (log(x), True))
integrate(exp(-x**2), x)            # sqrt(pi)*erf(x)/2
integrate(x*exp(x)*sin(x), x)
# x*exp(x)*sin(x)/2 - x*exp(x)*cos(x)/2 + exp(x)*cos(x)/2
integrate(sin(x)**3*cos(x)**2, x)   # cos(x)**5/5 - cos(x)**3/3
```

### Rational functions

`ratint` (`rationaltools.py`) is Bronstein's algorithm as in SymPy: Hermite
reduction for the rational part, then the Lazard-Rioboo-Trager (LRT)
resultant for the logarithmic part, converted to real `log`/`atan` terms by
`log_to_real` and `log_to_atan`, or kept as a `RootSum` when the field is
not real or the roots are not found. Over QQ and ZZ (coefficients up to
4096 bits) the Rust driver (`integrate/ratint_driver.rs` with the Hermite,
LRT, `log_to_atan` and `log_to_real` kernels) does the whole computation;
elsewhere the Python driver calls those native kernels when they accept the
input. A selector in `rationaltools.py` takes the native route only when
the expression operators and `atan` are the stock ones (the vector and
differential geometry modules install their own). A refused native route
runs the full Python algorithm, never a "no result".

### Risch, heurisch, manual and Meijer G

`risch_integrate` (`risch.py`, `rde.py`, `prde.py`) is SymPy's transcendental
Risch algorithm for towers of `exp` and `log`. Polynomial reduction,
`splitfactor` and the derivation run in Rust; in the parametric Risch
differential equation, `prde_normal_denom`, `prde_special_denom`, the
linear constraints, `constant_system`, `prde_spde`, the no-cancellation
cases, `param_rischDE`, `limited_integrate` and its reduction,
`parametric_log_deriv_heu`, `is_deriv_k` and the base case of
`is_log_deriv_k_t_radical_in_field` run in Rust (`integrate/prde*.rs`),
called from the Python driver. The non-base cases of the last one, and the
drivers, remain Python.

`heurisch` (`heurisch.py`, the Risch-Norman heuristic) and `manualintegrate`
(`manualintegrate.py`, rule tree with `integral_steps`) are Python ports;
the cartesian coefficient preparation of `heurisch` is in Rust
(`integrate/cartesian.rs`). The Meijer G method is entirely in Rust
(`integrate/meijerint/`: table, conditions, indefinite, definite and inverse
Mellin transforms); `symox.integrals.meijerint` exposes
`meijerint_indefinite`, `meijerint_definite` and `meijerint_inversion`,
which return `None` when the method does not apply:
`meijerint_definite(exp(-x)*x**a, x, 0, oo)` is
`(gamma(a + 1), re(a) > -1)`.

### Not found and proven non-elementary

An unevaluated `Integral` means "no method found an antiderivative"; it is
not a proof. A `NonElementaryIntegral` (a subclass of `Integral`, in
`symox.integrals.risch`) means the Risch algorithm proved that no
elementary antiderivative exists. It prints like an `Integral` and keeps its
class through `subs`, `doit`, `simplify` and `expand`; an equal `Integral`
built elsewhere is a plain `Integral`.

```python
from symox.integrals.risch import risch_integrate, NonElementaryIntegral

r = risch_integrate(exp(-x**2), x)
assert isinstance(r, NonElementaryIntegral)       # Integral(exp(-x**2), x)
assert isinstance(integrate(exp(x**2), x, risch=True), NonElementaryIntegral)
assert not isinstance(integrate(sin(x)/log(x), x), NonElementaryIntegral)
```

## Options

- `risch=True`, `heurisch=True`, `manual=True` or `meijerg=True` runs only
  that method; more than one raises `ValueError("At most one of manual,
  meijerg, risch, heurisch can be True")`. `False` removes the method from
  the pipeline.
- `conds='piecewise'` (default), `'separate'` or `'none'` choose how
  convergence conditions of definite integrals are returned; another value
  raises `ValueError`.
- `Integral.doit(final=False)` omits the continuity correction that the
  final round adds to each `atan(c*tan(a))`:
  `sign(c)*pi*floor((a - pi/2)/pi)`. As in SymPy, `final` is a `doit` hint
  only.
- As in SymPy, any other keyword of `integrate` (including `final`) and any
  keyword of `Integral(...)` raise `TypeError("Basic.__new__() got an
  unexpected keyword argument '...'")`. `Integral.doit` ignores hints it
  does not use, and `doit(integrals=False)` returns the integral unchanged.

```python
integrate(1/(2 + cos(x)), x)
# 2*sqrt(3)*(atan(sqrt(3)*tan(x/2)/3) + pi*floor((x/2 - pi/2)/pi))/3
Integral(1/(2 + cos(x)), x).doit(final=False)
# 2*sqrt(3)*atan(sqrt(3)*tan(x/2)/3)/3
```

## Definite integrals

The integration variable is replaced by a dummy that is positive, negative
or real according to the bounds (not when `DiracDelta` is present), and
reversed bounds are swapped with a sign change. Improper integrals with an
infinite bound try `meijerint_definite` first and keep its convergence
condition; `P(x)*exp(-a*x**q)` over `(0, oo)` or `(-oo, oo)` uses the gamma
moment formula under its condition; rational integrands are analysed for
poles in the interval (`rational_definite.rs`, `rational_conditions.rs`).
Otherwise the antiderivative is evaluated at the bounds, by substitution or
by a one-sided limit at an infinite or singular endpoint. Jumps at interior
poles and logarithmic branch points are corrected with one-sided limits; a
non-integrable pole gives `nan`, as in SymPy. A `DiracDelta` on an
endpoint keeps the integral held.

With symbolic parameters the condition becomes a `Piecewise` whose last
branch is the held integral:

```python
from symox import symbols, integrate, exp, sin, cos, tan, pi, oo, Integral, S

x, a = symbols('x a')
p = symbols('p', positive=True)
integrate(exp(-a*x), (x, 0, oo))
# Piecewise((1/a, Abs(arg(a)) < pi/2), (Integral(exp(-a*x), (x, 0, oo)), True))
integrate(exp(-a*x), (x, 0, oo), conds='separate')   # (1/a, Abs(arg(a)) < pi/2)
integrate(exp(-a*x), (x, 0, oo), conds='none')       # 1/a
assert integrate(exp(-p*x), (x, 0, oo)) == 1/p
assert integrate(sin(x)/x, (x, 0, oo)) == pi/2
assert integrate(tan(x), (x, 0, pi)) is S.NaN
Integral(1/x, (x, -1, 2)).principal_value()          # log(2)
```

An integral over an infinite interval whose antiderivative oscillates gives
accumulation bounds (`oscillating_integral` in
`python/symox/calculus/_accum_limit.py`), as in SymPy:

```python
from symox import AccumBounds
assert integrate(sin(x), (x, 0, oo)) == AccumBounds(0, 2)
assert integrate(cos(x), (x, 0, oo)) == AccumBounds(-1, 1)
```

## Multiple, line and special integrands

Several limits are integrated innermost first, so a bound may depend on an
outer variable. With one argument and several free symbols, `integrate`
raises `ValueError("specify dummy variables for ...")`.

```python
from symox import Curve, line_integrate
y, t = symbols('y t')
assert integrate(x*y, (y, 0, x), (x, 0, 1)) == S(1)/8
assert integrate(x*y, (x, 0, 1), (y, 0, x)) == x**2/4
assert integrate(x*y, x, y) == x**2*y**2/4
line_integrate(x + y, Curve([t, t**2], (t, 0, 1)), [x, y])
# -1/12 - asinh(2)/64 + 67*sqrt(5)/96
```

An `Eq` integrand is integrated on both sides, `integrate(1 + x + O(x**3),
x)` is `x + x**2/2 + O(x**4)`, vector integrands are integrated
componentwise and `_eval_Integral` hooks are honored. A `Sum` integrand
whose limits depend on an integration variable raises `ValueError`.

## The `Integral` object

`Integral` has SymPy's `function`, `limits`, `variables`, `free_symbols`,
`doit(**hints)` (same hints as `integrate`), differentiation under the
integral sign, `transform`, `as_sum` and `principal_value`
(`_integral_extras.py`).

```python
I2 = Integral(x**2, (x, 0, 2))
assert I2.as_sum(4) == S(21)/8                     # midpoint by default
assert I2.as_sum(4, method='trapezoid') == S(11)/4
Integral(x, (x, 0, 1)).transform(x, 2*t)           # Integral(4*t, (t, 0, 1/2))
```

`as_sum` accepts `left`, `right`, `midpoint` and `trapezoid`; a
multidimensional integral raises `NotImplementedError` and an infinite
interval `ValueError`, as in SymPy. Integral transforms (Laplace, Fourier,
Mellin, Hankel, sine, cosine and their inverses) are exported from
`symox.integrals.transforms` and `symox.integrals.laplace` and are not
described here.

## Numerical evaluation

`Integral(...).evalf(n)` and `N` first evaluate the result of `doit()` at
the requested precision when `doit` finds a closed form. Otherwise an
integral with one `(x, a, b)` limit and no free symbol is computed by
numerical quadrature, as SymPy's `evalf_integral`
(`python/symox/_quadrature.py`): mpmath's `quadts` (tanh-sinh, with its error
estimate), or `quadosc` for `evalf(quad='osc')`, which never calls `doit`.
The integrand is evaluated by the native evaluator at each node. The working
precision grows with SymPy's schedule until the estimated accuracy reaches
the target, up to `maxn` digits (default 100). symox has no native
quadrature; the mpmath quadrature is used as in SymPy.

Differences: an integral whose estimated accuracy stays below the target
at `maxn` stays unevaluated (SymPy returns a Float of the lower accuracy);
so does a zero quadrature result, which SymPy returns as a Float without a
correct digit. Integrals with several limits or free symbols stay
unevaluated, as in SymPy.

```python
from symox import N, log
N(Integral(exp(-x**2), (x, 0, oo)), 30)    # 0.886226925452758013649083741671
Integral(sin(x)/log(x), (x, 2, 3)).evalf()  # 0.674380603436754
Integral(sin(x)/x**2, (x, 1, oo)).evalf(20, quad='osc')
Integral(sin(x)/x, (x, 0, 1)).evalf(quad='osc')
# ValueError: quadosc requires an infinite integration interval
```

## Limits and differences from SymPy

- Result forms can differ from SymPy's where another stage or `simplify`
  candidate wins. On SymPy 1.14's integration test files, the last replay
  (2026-10-29, parity notes section "integsuite") had failing tests
  (failures and errors) in test_integrals 51/205, test_risch 8/36,
  test_heurisch 3/33, test_meijerint 13/33, test_manual 0/42 and
  test_rationaltools 0/12, with no wrong result identified. These files
  have not been replayed since the later integration waves.
- Numerical quadrature goes through mpmath (above).
- `integrate(1/(x**8 + 16), x)` returns a `RootSum`; SymPy 1.14 returns `0`,
  which is wrong. This deviation is pending validation.
- Still in Python: the rational driver outside QQ/ZZ, general `log_to_real`,
  the Risch, PRDE and heurisch drivers, `manualintegrate`, the non-base
  cases of `is_log_deriv_k_t_radical_in_field`. Open gaps (the `log_to_atan`
  form over EX, `RootSum` with `real=False`, a `TypeError` over
  `QQ<sqrt(2)>`) are tracked in `TODO.md` and
  [the parity notes](algorithmic-parity.md).
