# Series and expression analysis

[README](../README.md) · [Compatibility](compatibility.md)

See also [Limits and asymptotic orders](limits.md) and [Integration](integration.md).

The public `Expr` method surface measured against SymPy 1.14 is present in
full. The methods below execute native symbolic operations; the bounds and
unsupported cases describe their current behavioral scope.

## Formal power series

`fps(f, x=None, x0=0, dir=1, hyper=True, order=4, rational=True, full=False)`
and `expr.fps(...)` return a lazy `FormalPowerSeries` for supported ordinary
integer-power expansions. Exact finite coefficients come from the native
series engine without importing SymPy; `FormalPowerSeries` itself is a lazy
Python wrapper around those expansions.

```python
from symox import symbols, sin, exp, fps

x = symbols('x')
s = sin(x).fps(x)
s.truncate(8)               # x - x**3/6 + x**5/120 - x**7/5040 + O(x**8)
s.polynomial(6)             # x - x**3/6 + x**5/120
s[3]                       # -x**3/6
s.ak.coeff(3)               # -1/6
fps(exp(x), x).product(s, n=5).truncate()
```

The result exposes `function`, `x`, `x0`, `dir`, coefficient and monomial
sequences, indexing, finite slices and iteration. `truncate(None)` returns
a lazy term iterator. `product`, `compose`, `inverse`, differentiation and
integration use exact native expressions and series. Composition requires
the inner constant to equal the expansion point; inversion requires a
provably nonzero constant coefficient.

Orders are bounded at 256. Generic symbolic coefficient formulas and SymPy's
complete formal-series expression hierarchy are not implemented. The
algorithm-selection flags are accepted as hints to the native recurrence
engine; disabling both `hyper` and `rational` retains the input. Nonordinary
expansions such as `sqrt(x)`, `log(x)`, `1/x` and `exp(1/x)` at zero retain the
original expression. As in SymPy's `_compute_fps`, a factor `x**n` with a
symbolic exponent free of `x` is separated at zero from the right: the rest is
expanded and each term multiplied by `x**n`, so
`fps(x**n*sin(x**2), x).truncate(8)` is
`-x**(n + 6)/6 + x**(n + 2) + O(x**(n + 8))`. Symbolic coefficients such as
`fps(exp(a*x), x)` are expanded (SymPy 1.14 returns the input there).
Unsupported native series operations raise explicitly.
At nonzero centers, coefficients are checked against Taylor expansion;
this avoids known SymPy 1.14 `fps` coefficient errors, for example with
`x**2 + 1` expanded about one.

## Fourier series

`fourier_series(f, limits=None, finite=True)` and
`expr.fourier_series(limits=None)` construct exact sine/cosine expansions.
The default interval is `(-pi, pi)` for a univariate expression.

```python
from symox import symbols, pi

x = symbols('x')
s = x.fourier_series((x, -pi, pi))
s.truncate(3)               # 2*sin(x) - sin(2*x) + 2*sin(3*x)/3
s[1]                       # 2*sin(x)
s.a0                       # 0
s.sigma_approximation(4)
s.shiftx(1).truncate(3)
```

Coefficients use native exact definite integration. Supported affine real
boundaries in `Abs`, `sign`, `Heaviside` and `Piecewise` split the interval
before integration. Trigonometric polynomials have finite support and return
`FiniteFourierSeries`. Other coefficients are computed lazily. `truncate(n)`
counts nonzero terms, including the constant term when nonzero;
`truncate(None)` returns an iterator. Sigma approximation uses harmonics below
`n` weighted by `sinc(k*pi/n)`. Shifts and scalings of values and arguments
preserve the Fourier expansion's coefficients and transform its basis.

The coefficient search is bounded at harmonic 4096, with at most 128 interval
boundaries. Unsupported integrals and undecidable branch ordering raise
`NotImplementedError`; nonfinite integrals raise `ValueError`. The result
objects expose the implemented series protocol rather than all SymPy `Basic`
methods. Missing harmonics in finite sigma approximations contribute zero,
avoiding SymPy 1.14's incidental `IndexError`. Functions provably constant
inside the interval terminate immediately, including cases where SymPy's
nonzero-term search does not terminate.

## Asymptotic expansions

`expr.aseries(x=None, n=6, bound=0, hir=False)` expands at positive infinity.
It supports native inverse-power, logarithmic and Puiseux expansions, plus
comparable decaying exponential scales. `hir=True` keeps coefficients in
their hierarchical form; otherwise coefficients expand recursively until
the dominant coefficient's remainder determines the result's accuracy.

```python
from symox import sin, exp, O, oo

f = sin(1/x + exp(-x)) - sin(1/x)
f.aseries(x, n=4)
# (-1/(2*x**2) + 1 + O(x**(-4), (x, oo)))*exp(-x)
f.aseries(x, n=3, hir=True)

g = 1/(1 + exp(-x) + exp(-x*x))
assert g.aseries(x, n=3) == 1-exp(-x)+exp(-2*x)+O(exp(-3*x), (x, oo))
```

The order is limited to 256 and `bound` to 32. A positive `bound` also splits
supported growing exponential factors using an inverse-power expansion of
their exponent's coefficient; an exact residual exponent preserves the
original value. The factorization is held to prevent immediate recombination.
Limits of exponent ratios select the fastest comparable exponential class.
Slower exponentials remain coefficient functions; rational leading ratios
retain their exact residual exponents. For example, `exp(-x*x-x)` keeps its
factor `exp(-x)` when expanding in the class of `exp(-x*x)`. Hierarchical
geometric expansions are verified with exact rational remainder identities.
Undecidable comparisons, irrational leading ratios and general Gruntz
most-rapidly-varying normalization remain unsupported.
Growing exponentials are preserved when the native series engine keeps them
exact. Polar-branch identities used here apply along the positive real
direction at infinity.

## Pade approximants and sequence limits

`approximants(coefficients, X=Symbol('x'), simplify=False)` generates consecutive
Pade approximants through a formal continued fraction. Coefficients can be
exact numbers or symbolic expressions. Exact cancellation terminates the
generator when the coefficient list already describes a rational function;
it does not emit a spurious zero or NaN after termination. Input is bounded at
256 coefficients and is copied before computation.

```python
from symox import approximants, limit_seq, factorial, cancel

assert cancel(list(approximants([1]*8, x))[-1] - 1/(1-x)) == 0
n = symbols('n', integer=True, positive=True)
assert limit_seq((-1)**n/n, n) == 0
assert limit_seq(2**n/factorial(n), n) == 0
```

`limit_seq(expr, n=None, trials=5)` is a Python port of SymPy's
`sympy.series.limitseq` (Kauers' algorithm): up to `trials` rounds of a
Stolz-Cesaro difference step on numerator and denominator, followed by their
`dominant` terms, which `symox.series.limitseq` also exports. Each round asks
the native limit engine at `oo`; an unevaluated native `Limit` counts as
undecided. Gamma-type terms bypass that engine and use `simplify(combsimp(...))`
and an exact ratio test, which also decides powers of nonreal numbers.
`fibonacci` and `lucas` are rewritten in Binet form; `factorial` and
`subfactorial` in `gamma`. Even/odd subsequences handle alternating signs and
trigonometric terms; different comparable limits return an `AccumBounds`
interval. Nested `dominant` calls are limited to depth 4. Undecided cases
return `None`, which is not a proof that the limit does not exist; some
sequences that SymPy decides still return `None` here.

`AccumBounds(min, max)` is an expression node ported from SymPy's
`AccumulationBounds`: it appears inside sums, products and powers, and its
arithmetic, powers, comparisons, intersection and union follow the reference
methods. It is the same class as `AccumulationBounds` and collapses equal
endpoints to a scalar. Elementary functions map bounds to bounds, and
`sin(oo)`, `cos(oo)` and `frac(oo)` evaluate to bounds as in SymPy.

The compatible `symox.series` module also exports the existing native
`Limit`, `limit`, `Order`, `O`, `residue` and `difference_delta` implementations.
Their `limits`, `order`, `residues` and `limitseq` module paths are available.
`gruntz(expr, variable, point, dir='+')` uses the Rust exact limit engine,
including its most-rapidly-varying comparison algorithm, and raises `PoleError`
when that engine cannot determine a result. Ordinary `limit` retains a held
`Limit` for unresolved cases. When the native engine leaves an oscillating
limit unevaluated, a Python term-by-term evaluation can return accumulation
bounds, for example `limit(sin(x), x, oo) == AccumBounds(-1, 1)`; one-sided
limits onto a branch cut, such as `limit(log(I*x - 1), x, 0, '-') == -I*pi`,
are rewritten to the non-principal side first. Series expansion and limit
solving have distinct algorithmic scopes; the restrictions on `aseries` above
still apply.

## Local regularity, branches and constructors

`is_meromorphic(x, a)` is a three-valued local test. It handles arithmetic,
integer and noninteger powers, elementary exponential/logarithmic and
trigonometric functions, supported inverse functions and gamma/error functions.
Unknown functions and undecidable cancellations return `None`. This is a
local test, not a global singularity finder.

`exp_polar(z)` is a native symbolic head distinct from an undefined function
of the same name. It preserves winding information, differentiates, parses,
pickles and converts explicitly to/from SymPy. `extract_branch_factor` returns
`(residual, turns)` for the `exp_polar(2*pi*I*turns)` factor, with
`allow_half=True` permitting half turns. Numerical evaluation and numerical
compilation of polar exponentials remain unsupported; they are not silently
replaced by ordinary exponentials.

`as_dummy()` renames bound variables in native integrals, sums, products and
lambdas, strips their assumptions and preserves free endpoint occurrences.
`canonical_variables` exposes the collision-free mapping. On a symbol,
`as_dummy()` returns a fresh `Dummy` with default commutativity assumptions.

`fromiter(iterable, **options)` rebuilds the concrete constructor through
classes or expression instances, including `Add`, `Mul`, `Pow`, `Tuple` and
function heads. `Expr.fromiter` constructs a generic native expression node.
Native DAG leaves normalize Python numbers into symbolic numbers; arbitrary
raw Python objects are not retained as generic expression arguments.
