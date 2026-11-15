# Scalar real and complex analysis

[README](../README.md)

The analysis adapters use native scalar expressions, equation solvers,
inequality solvers and symbolic sets. They do not import SymPy during calculations.

## Singularities

`singularities(expression, symbol, domain=None)` collects zeros of negative-power
bases, logarithm arguments and elementary reciprocal-function denominators.
Inverse hyperbolic branch endpoints are included for `atanh` and `acoth`.
The default domain is real for a symbol declared real, and complex otherwise.
An explicit domain can restrict the result to a real interval or another native set.

```python
from symox import symbols, singularities, S, FiniteSet, tan, Interval, pi

x = symbols('x')
assert singularities(1/(x*x-1), x, S.Reals) == FiniteSet(-1, 1)
assert singularities(tan(x), x, Interval(0, pi)) == FiniteSet(pi/2)
```

Branch points are not exhaustively classified. Unsupported special functions,
undefined functions, held derivatives and powers with undecidable exponent signs
raise `NotImplementedError`; an empty result is not fabricated for those cases.
Unresolved equations can remain `ConditionSet` objects.

The native equation solver represents complete complex solution families for
affine arguments of `sinh`, `cosh`, `tanh` and `coth` using integer `ImageSet`
lattices. In particular, hyperbolic denominator poles include their remote
imaginary translates. Real argument/range restrictions are checked separately.
Unsupported nonlinear hyperbolic equations retain a `ConditionSet`.

## Real continuity and stationary points

`symox.calculus.util.continuous_domain(f, x, domain)` supports elementary real
expressions involving polynomials, rational functions, exponentials, trigonometric
functions, logarithms, powers and supported inverse functions. It removes poles
and intersects real branch inequalities using native set operations. Bounds must
be independent of the analysis variable, and the domain must be known to be real.

Analysis uses a temporary real variable to respect the supplied domain even when
the original symbol has no real assumption. Powers follow the native principal
complex branch: negative inputs to noninteger rational powers are excluded from
the real domain. Real powers require a known sign for the exponent. Unsupported
complex intermediate values, including absolute values of expressions not known
real, raise `NotImplementedError`.

`stationary_points(f, x, domain=S.Reals)` solves the native first derivative
equal to zero within the real continuous domain. Finite trigonometric intervals
produce finite sets of points. Nondifferentiable sign/absolute-value branches
and unknown continuity domains are unsupported.

```python
from symox import stationary_points, sin

assert stationary_points(x*x, x) == FiniteSet(0)
assert stationary_points(sin(x), x, Interval(0, 2*pi)) == FiniteSet(pi/2, 3*pi/2)
```

## Periods

`symox.calculus.util.lcim(numbers)` computes the least common integral multiple
of commensurable finite real constants. Exact rational ratios identify a shared
irrational scale, so `[2*pi, 3*pi, pi/2]` gives `6*pi`, while
`[sqrt(2), sqrt(3)]` returns `None`. Rational numerators use integer LCM and
denominators use integer GCD. Incommensurable or unproved numeric ratios return
`None`; empty input returns one.

`periodicity(f, x, check=False)` is a Python port of SymPy 1.14's algorithm
over a real dummy variable. After `simplify` (with SymPy's deep `trigsimp`
step), it handles trigonometric functions of affine arguments, `Abs` (half the
period of `Abs(sin(...))`), exponentials through the periods of their real and
imaginary parts, powers, products and sums combined with the internal `lcim`,
`Mod` and compositions found by `decompogen`. Relationals use `lhs - rhs`;
`Piecewise` returns `None`. Constants return zero. `None` means this search
found no period; it is not a general proof of aperiodicity.

```python
from symox import periodicity, cos

assert periodicity(sin(2*x), x) == pi
assert periodicity(sin(x)**2, x) == 2*pi
assert periodicity(sin(x)+cos(2*x), x, check=True) == 2*pi
```

As in SymPy, returned periods need not be the smallest period: the power rule
takes the base's period, so `sin(x)**2` gives `2*pi`. `check=True` verifies
the candidate with `equals` after shifting, and raises `NotImplementedError`
with SymPy's message when that identity is not recognized.

## Monotonicity

`symox.calculus.singularities.monotonicity_helper(expression, predicate,
interval=S.Reals, symbol=None)` checks a custom predicate on the native first
derivative over the requested interval. For example, `lambda d: d > 0` tests
strict positivity of every derivative value, including finite endpoints; it
does not infer strict monotonicity from isolated zeros. The native proof uses
the same real continuity and connected-domain checks as the monotonicity
adapters. Empty intervals are vacuous and singleton intervals evaluate the
predicate at that point. When the native proof declines (disconnected or
non-interval domains, parameters not known to be real, bounds depending on
the variable, functions without a native continuity rule such as `Abs`) or is
undecided, the result is SymPy's: the interval is compared with the
`solveset` solution of the derivative predicate over the reals. A
multivariate expression without `symbol` raises SymPy's
`NotImplementedError`, and an undefined function keeps an explicit
`NotImplementedError` where SymPy overflows the recursion limit. The native
proof can be more decisive than that comparison: with `lambda d: d > 0`,
`x**3` over `S.Reals` gives `False` (the derivative vanishes at zero) and,
with `lambda d: d < 0`, `1/x` over `S.Reals` gives `False`, where SymPy
returns `None`. For a parameter `a`, real or not, `a*x` gives `None`, while SymPy
returns `False` because its `solveset` of the parameter-only condition is
empty.
Conjunctions, disjunctions and negations of derivative relations use native
set intersection, union and complement.

`is_increasing`, `is_decreasing`, their `is_strictly_...` variants and
`is_monotonic` accept `(expression, interval=S.Reals, symbol=None)`. A single
free variable can be inferred. The native proof below runs first. When it
declines or is undecided (parameters not known real, bounds depending on the
variable, disconnected or non-interval domains), the result is SymPy's
derivative-only answer: the interval is compared with the `solveset` solution
of the derivative predicate over the reals, or, for `is_monotonic`, with the
turning points. Results are `True`, `False`, or `None` when neither decides.
Boolean conditions, sets and tuples are rejected as scalar calculus inputs.

```python
from symox import is_increasing, is_strictly_increasing, is_decreasing, is_monotonic, sqrt

assert is_strictly_increasing(x**3, symbol=x) is True
assert is_strictly_increasing(sqrt(x), Interval(0, S.Infinity), x) is True
assert is_decreasing(1/x, Interval.open(0, S.Infinity), x) is True
assert is_monotonic(x*x, symbol=x) is False
```

These checks prove derivative signs using native inequalities and set operations.
The function must be real, defined and continuous throughout the requested
interval; a pole inside it disproves that requirement. Derivative continuity is
checked on the interval interior, so closed endpoints can have an infinite
derivative, as with `sqrt(x)` at zero. Finite isolated derivative zeros do not
prevent strict monotonicity. A constant is monotonic, and empty or singleton
domains are handled by the usual vacuous definitions.

These mathematical cases can differ from derivative-only helper results in
SymPy 1.14: `x**3` is strictly increasing despite its stationary point, a constant
is monotonic, and `1/x` is not a function defined throughout `S.Reals`.

The native proof domain is a connected real interval, all reals, or a
singleton/empty set, with parameters known to be real. Other domains get the
SymPy fallback above, which inherits SymPy's derivative-only semantics and
does not check continuity. Unknown functions and inequality families beyond
the native solver raise `NotImplementedError`; unresolved derivative domains
or zero-set classifications can return `None`. This is a proof procedure for
supported expressions, not a complete decision algorithm.

`is_convex(f, x, domain=S.Reals)` uses the same interval machinery to prove that
the first derivative is nondecreasing in the interior. The original function
must be continuous on the requested domain. This supports closed endpoints with
an infinite first derivative, for example `is_convex(-sqrt(x), x,
domain=Interval(0, S.Infinity))`. One variable must be supplied. Multivariate
convexity and nonsmooth proofs such as `Abs(x)` remain unsupported.

## Infimum and supremum

`minimum(f, x, domain=S.Reals)` and `maximum(f, x, domain=S.Reals)` combine native
critical-point equations, finite function values and one-sided endpoint limits.
They use the real continuous parts of the requested domain. Open endpoints and
infinite bounds contribute limits, so the returned infimum or supremum need not
be attained by the function.

```python
from symox import minimum, maximum, exp

assert minimum(x*x+2*x+5, x) == 4
assert maximum(sin(x), x) == 1
assert minimum(exp(x), x) == 0
assert minimum(1/x, x) == -S.Infinity
assert maximum(x, x, Interval.open(0, 1)) == 1
```

Domains can be finite unions of explicit real intervals and finite sets. A known
period reduces an unbounded interval, or an interval spanning at least a period,
before critical points are enumerated. Formal derivative poles and zeros of
absolute-value arguments are included as candidate points when the original
function remains continuous there. An infinite endpoint limit can prove the
corresponding bound without resolving all stationary points.

Critical-point sets must otherwise be finite and complete. Unresolved equations,
symbolic signs of infinite limits and unsupported domain decompositions raise
`NotImplementedError`. Other parameters must be known real. An empty requested
domain, or one with no real function values, raises `ValueError`. Domain
components, candidate points and distinct candidate values are limited to 2,048.

## Function images

`symox.calculus.util.function_range(f, x, domain)` returns an image set built
from continuous real components. It combines their exact bounds with an
attainment check and joins overlapping image intervals. Finite input components
produce finite sets. Empty input or an empty real continuous domain returns the
empty set. As in SymPy, a whole-line image is `Interval(-oo, oo)` rather than
`S.Reals`.

```python
from symox.calculus.util import function_range
from symox import cos, Union

assert function_range(exp(x), x, S.Reals) == Interval.open(0, S.Infinity)
assert function_range(1/x, x, S.Reals) == Union(Interval.open(-S.Infinity, 0), Interval.open(0, S.Infinity))
assert function_range(cos(x), x, Interval.open(0, 2*pi)) == Interval.Ropen(-1, 1)
```

Endpoint attainment is checked on the original component. For a whole open
period, `cos(x)` omits its maximum, while a period with either endpoint included
attains it. This also preserves the gap at zero in the image of `1/x`. A compact
continuous input interval has attained finite extrema. When an endpoint's
attainment cannot be proved with native sets and equations, the function raises
`NotImplementedError` instead of guessing its inclusion.

## Nonempty finite-set domains

`not_empty_in(set_expression, x)` finds the real values of `x` for which at
least one defined real member remains in a finite set or its intersection
with intervals, finite sets, unions and complements.

```python
from symox import FiniteSet, not_empty_in, log

assert not_empty_in(FiniteSet(x/2).intersect(Interval(0, 1)), x) == Interval(0, 2)
assert not_empty_in(FiniteSet(0, x).intersect(Interval(0, 1)), x) == S.Reals
assert not_empty_in(FiniteSet(log(x)), x) == Interval.open(0, S.Infinity)
```

A constant member can prove nonemptiness on its own. Reciprocal members retain
their excluded poles. Functions such as `floor`, `sign` and `acot` remain
defined at their real jumps; continuity is not required for membership.
Unsupported membership equations or inequalities, unresolved real coefficients
and multivariate requests raise `NotImplementedError`.
