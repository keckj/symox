# Holonomic functions

[README](../README.md) · [Differential equations](differential-equations.md) ·
[Special functions](special-functions.md)

`symox.holonomic` supplies the nine exports in SymPy 1.14's corresponding
package. This is an API presence count. The implemented algorithms and their
limits are described below; it is not a claim of complete behavioral parity.
The representation follows the [SymPy holonomic documentation](https://docs.sympy.org/latest/modules/holonomic/about.html).
All coefficients, scalar derivatives and exact linear algebra use Symox's
native engine. The holonomic algorithms themselves are Python control flow in
`python/symox/holonomic` over those native polynomial and domain kernels.
There is no runtime import or delegation to SymPy. The `holonomicerrors`
module provides SymPy's exception classes (`NotHolonomicError`,
`NotPowerSeriesError`, `NotHyperSeriesError`, ...) and messages.

## Differential and shift operators

`DifferentialOperators(base, generator)` returns an algebra and its derivative
operator. `DifferentialOperator(coefficients, algebra)` constructs
`sum(p_i(x) D**i)` from coefficients in increasing order. The Ore product uses
the Leibniz rule, including derivatives of coefficients:

```python
from symox import QQ, Symbol
from symox.holonomic import DifferentialOperators, HolonomicFunction

x = Symbol('x')
R, D = DifferentialOperators(QQ.old_poly_ring(x), 'Dx')
assert D*x == x*D + 1
assert D**2*x == x*D**2 + 2*D

sine = HolonomicFunction(D**2 + 1, x, 0, [0, 1])
assert sine.series(6, order=False) == x - x**3/6 + x**5/120
```

`RecurrenceOperators` and `RecurrenceOperator` use the shift relation
`S*p(n) = p(n+1)*S`. `HolonomicSequence(operator, u0)` stores the recurrence
and initial values; `eval(n)` propagates them exactly. A vanishing recurrence
leading coefficient requires additional initial values.

```python
from symox import ZZ
from symox.holonomic import RecurrenceOperators, HolonomicSequence

n = Symbol('n', integer=True)
R, S = RecurrenceOperators(ZZ.old_poly_ring(n), 'Sn')
assert S*n == (n + 1)*S
fib = HolonomicSequence(S**2 - S - 1, [0, 1])
assert fib.eval(10) == 55
```

Both algebras implement addition, subtraction, composition products,
nonnegative integer powers and exact coefficient division where the result
remains polynomial. Incompatible algebras and negative operator powers are
rejected. Operator division with a remainder is not implemented. The displayed
generator is a label; operator noncommutativity is implemented by these classes,
not by an ordinary scalar symbol.

## Functions, closure and conversion

`HolonomicFunction(annihilator, x, x0=0, y0=None)` stores a positive-order
homogeneous polynomial differential equation and ordinary derivative values
`[f(x0), f'(x0), ...]`, or Frobenius conditions `{p: [c0, c1, ...]}` at a
singular point (the coefficients of `x**p`, `x**(p+1)`, ...). `is_singularics()`
is True for a dictionary, False for a list and None without initial values.
A missing initial list leaves the solution unspecified. Zero and order-zero
differential equations are rejected as function constructors.

Addition and multiplication build direct-sum and tensor-product differential
systems. Exact native nullspaces find a cyclic relation for their output.
Differentiation and rational composition use the same construction. This is
an exact closure algorithm, not sampling or fitting to numerical values.
Scalar rational multiplication/division, integration, `shift_x`, `change_x`
and coefficient-domain unification are supported. Sums and products with
conditions at different points move them to a common point as in SymPy.
Definite integration from a removable singularity takes the antiderivative's
lower-limit value as a limit; undefined values raise or leave the initial
conditions unavailable instead of producing `nan`.

Powers follow SymPy's `__pow__`: a first-order equation accepts any exponent,
including negative, rational and symbolic ones; higher orders use
square-and-multiply for nonnegative integers, raise `NotHolonomicError` for
negative powers and `TypeError` for rational or symbolic powers. Division by a
holonomic function raises `NotImplementedError`, because the class is not
closed under division.

```python
from symox import exp, sin
from symox.holonomic import expr_to_holonomic

h = expr_to_holonomic(exp(x)*sin(x), x)
assert h.series(4, order=False) == x + x**2 + x**3/3
assert h.to_sequence()[0][0].eval(4) == 0
```

`expr_to_holonomic` recognizes rational functions and expressions with rational
logarithmic derivative, rational powers, sums/products, exponentials with
rational exponent, trigonometric and hyperbolic functions, logarithms,
`erf`/`erfc`, cylindrical Bessel functions, Airy functions, generalized
hypergeometric functions and Meijer G functions. Their argument substitutions
must be rational in the independent variable. A power `f**c` with `c` free of
`x` goes through `__pow__`, so `1/sin(x)` raises `NotHolonomicError` and
`sqrt(sin(x))` raises `TypeError`. Algebraic functions at a branch point and
rational functions with one simple indicial root receive Frobenius conditions,
for example `{1/2: [1]}` for `sqrt(x)*exp(x)`; removable singularities such as
`sin(x)/x` at zero take their values by limits. Unsupported conversions,
including `Ci`, raise `NotImplementedError`; they do not return a fabricated
annihilator.

`from_hyper` and `from_meijerg` construct the standard Euler-operator equations
and transform rational arguments. `initcond=False` skips initial derivative
evaluation. `evalf=True` evaluates available initial values numerically.
Singular special-function values can leave initial conditions unavailable;
choose a suitable `x0` or supply derivative values explicitly.

Coefficient rings use native rational/integer arithmetic. Conversion can infer
rational-function parameters, for example `expr_to_holonomic(sin(a*x), x)`.
Parameterized polynomial coefficients are validated in `x` and stored in a
native multivariate rational field, avoiding unsupported nested domain storage.
Float inputs select `RR` or `RR(parameters)` as in SymPy; the annihilator is
computed exactly with markers in place of the Floats, which are then restored.
General algebraic coefficient extensions are not implemented here. Annihilators
can differ from SymPy's by a constant factor because cyclic vectors are
normalized differently; the initial conditions agree.

## Series, expressions and numeric continuation

`series(n, coefficient=False, order=True)` propagates Taylor coefficients from
the differential equation. `coefficient=True` returns coefficients, not
derivatives. It handles ordinary points and analytic branches at singular
points when the supplied derivative values determine the missing terms.
Known converted expressions can supply an otherwise undetermined Taylor
derivative. Inconsistent numeric initial values raise an error (SymPy instead
returns a series with an undetermined constant). At singular points, `series`
and `to_sequence` follow SymPy's indicial equation and Frobenius step, including
fractional exponents and dictionary initial data; coefficients that the data
do not determine remain `C_i` constants. Logarithmic Frobenius branches are not
constructed, and a finite derivative list at a singular point can be
insufficient to identify a branch.

Closure operations inspect nonnegative integer roots of the indicial
polynomial and retain the extra initial derivatives they require. For example,
forming `x*sine` from the raw equation for `sine` retains `[0, 0, 2]` at zero;
the first two zero derivatives alone would leave its amplitude undetermined.

`to_sequence()` derives a polynomial recurrence for the Taylor coefficients.
Its default result is a list of `(HolonomicSequence, recurrence_start)` pairs,
or `(HolonomicSequence, exponent, recurrence_start)` triples at singular points;
the coefficients are relative to `x-x0`. `lb=False` returns the sequences.

`to_hyper()` is SymPy's algorithm on that recurrence. When it raises
`NotHyperSeriesError` but a known converted expression is a sum of
exponential, trigonometric/hyperbolic and polynomial terms, Symox still returns
the correct sum of hypergeometric terms (an intentional divergence; SymPy
raises). `to_meijerg()` converts each hypergeometric term, or preserves an
existing Meijer G expression. `to_expr()` retains a known converted expression;
otherwise it is `hyperexpand(to_hyper()).simplify()`, so Bessel equations give
`besselj` results. If that fails for ordinary initial values, the native ODE
solver is tried and its constants are fitted to the initial values;
reconstruction is then limited to that solver's supported equations. A general
solution that is undefined at the initial point raises `NotImplementedError`.

`evalf(points, method='RK4', h=0.05, derivatives=False)` follows the companion
system along straight complex segments using fourth-order Runge-Kutta or Euler
steps. As in SymPy, a list of points is the mesh itself: one step from each
point to the next, starting at `x0`, and `h` is unused. A single endpoint is
reached exactly with steps no longer than `h`. It requires complete initial
derivatives at an ordinary point and numeric coefficients along the path. It
uses machine precision, has a 100,000-step limit, and provides no adaptive
error bound or automatic singularity avoidance.
`change_ics` uses exact derivatives of a known expression or this numeric
continuation. General branch continuation around singularities is not provided.

The focused tests compare annihilators and Taylor coefficients with SymPy,
check exact Ore identities, test independent initial conditions, and verify
real/complex numerical continuation and execution with SymPy imports blocked.
