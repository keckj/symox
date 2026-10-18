# Holonomic functions

[README](../README.md) · [Differential equations](differential-equations.md) ·
[Special functions](special-functions.md)

`symox.holonomic` supplies the nine exports in SymPy 1.14's corresponding
package. This is an API presence count. The implemented algorithms and their
limits are described below; it is not a claim of complete behavioral parity.
The representation follows the [SymPy holonomic documentation](https://docs.sympy.org/latest/modules/holonomic/about.html).
All coefficients, scalar derivatives and exact linear algebra use Symox's
native engine. There is no runtime import or delegation to SymPy.

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
`[f(x0), f'(x0), ...]`. A missing initial list leaves the solution unspecified.
Zero and order-zero differential equations are rejected as function constructors.

Addition and multiplication build direct-sum and tensor-product differential
systems. Exact native nullspaces find a cyclic relation for their output.
Differentiation and rational composition use the same construction. This is
an exact closure algorithm, not sampling or fitting to numerical values.
Nonnegative integer function powers, scalar rational multiplication/division,
integration, `shift_x`, `change_x` and coefficient-domain unification are
supported. General division and negative powers of D-finite functions are
rejected because the class is not closed under those operations.

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
error functions, cylindrical Bessel functions, Airy functions, generalized
hypergeometric functions and Meijer G functions. Their argument substitutions
must be rational in the independent variable. Unsupported conversions raise
`NotImplementedError`; they do not return a fabricated annihilator.

`from_hyper` and `from_meijerg` construct the standard Euler-operator equations
and transform rational arguments. `initcond=False` skips initial derivative
evaluation. `evalf=True` evaluates available initial values numerically.
Singular special-function values can leave initial conditions unavailable;
choose a suitable `x0` or supply derivative values explicitly.

Coefficient rings use native rational/integer arithmetic. Conversion can infer
rational-function parameters, for example `expr_to_holonomic(sin(a*x), x)`.
Parameterized polynomial coefficients are validated in `x` and stored in a
native multivariate rational field, avoiding unsupported nested domain storage.
General algebraic coefficient extensions are not implemented here.

## Series, expressions and numeric continuation

`series(n, coefficient=False, order=True)` propagates Taylor coefficients from
the differential equation. `coefficient=True` returns coefficients, not
derivatives. It handles ordinary points and analytic branches at singular
points when the supplied derivative values determine the missing terms.
Known converted expressions can supply an otherwise undetermined Taylor
derivative. Inconsistent numeric initial values raise an error. Fractional
Frobenius powers, logarithmic Frobenius branches and dictionary initial data
are not implemented. A finite derivative list at a singular point can be
insufficient to identify a branch.

Closure operations inspect nonnegative integer roots of the indicial
polynomial and retain the extra initial derivatives they require. For example,
forming `x*sine` from the raw equation for `sine` retains `[0, 0, 2]` at zero;
the first two zero derivatives alone would leave its amplitude undetermined.

`to_sequence()` derives a polynomial recurrence for the Taylor coefficients.
Its default result is a list of `(HolonomicSequence, recurrence_start)` pairs;
the coefficients are relative to `x-x0`. `lb=False` returns the sequences.

`to_expr()` retains a known converted expression. Explicitly changed initial
values discard that expression and invoke the native ODE solver when an
expression is requested. Reconstruction from a raw annihilator is limited to
that solver's supported equations. `to_hyper()` reconstructs exponentials,
trigonometric/hyperbolic functions and their polynomial sums/prefactors, or
preserves an existing hypergeometric representation. General hypergeometric
reconstruction and general Meijer G reconstruction are not implemented;
`to_meijerg()` can preserve an existing Meijer G expression.

`evalf(points, method='RK4', h=0.05, derivatives=False)` follows the companion
system along straight complex segments using fixed-step fourth-order Runge–Kutta
or Euler integration. It requires complete initial derivatives at an ordinary
point and numeric coefficients along the path. It uses machine precision,
has a 100,000-step limit, and provides no adaptive error bound or automatic
singularity avoidance. Decreasing `h` improves accuracy on regular paths.
`change_ics` uses exact derivatives of a known expression or this numeric
continuation. General branch continuation around singularities is not provided.

The focused tests compare annihilators and Taylor coefficients with SymPy,
check exact Ore identities, test independent initial conditions, and verify
real/complex numerical continuation and execution with SymPy imports blocked.
