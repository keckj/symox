# Splines, variational calculus and integer equations

[README](../README.md)

## Piecewise expressions and splines

`piecewise_fold` distributes nested piecewise expressions through scalar
operations. `piecewise_exclusive` subtracts preceding branch conditions so that
the resulting conditions are disjoint. It appends a `nan` branch for uncovered
points unless `skip_nan=True`; `deep=False` limits the transformation to the
outer expression. Boolean simplification follows the native logic engine's
supported reductions.

`bspline_basis` and `bspline_basis_set` use the Cox–de Boor recurrence over exact
knots, returning piecewise polynomials. Repeated knots are supported; knots must
be real, comparable and nondecreasing. `interpolating_spline` builds clamped
knots, solves the exact collocation system and returns the interval polynomials.
Interpolation requires strictly increasing, concrete comparable abscissas, equal
numbers of abscissas/ordinates and at least `degree+1` data points. The degree
must be positive. Symbolically unordered knots and periodic boundary conditions
are outside this implementation.

## Euler equations and line integrals

`euler_equations(L, funcs=(), vars=())` extracts dependent functions and their
derivative jets, takes algebraic partial derivatives with temporary native
symbols, then applies the Euler–Lagrange total derivatives. It supports several
dependent functions, independent variables and higher-order derivatives.
Dependent function arguments must agree with the independent variables.

```python
from symox import Function, symbols, euler_equations
x = symbols('x')
f = Function('f')(x)
euler_equations(f.diff(x)**2/2 - f**2/2, f, x)
# [Eq(-f(x) - Derivative(f(x), (x, 2)), 0)]
```

`line_integrate(field, curve, vars)` substitutes a `geometry.Curve` into a scalar
field and integrates against arc length. The parameter must be distinct from
the field's coordinate symbols. Evaluation depends on the native integration
engine; unsupported integrals retain that engine's held form.

## Singularity functions

`SingularityFunction(x,a,n)` represents the truncated power for integral order
`n>=0` and Dirac delta derivatives for negative integral orders. It supports
exact point evaluation, substitution, argument differentiation and rewrites to
`Piecewise`, `Heaviside` or `DiracDelta`. The order-zero truncated power uses the
right-continuous value one at its origin. Nonintegral orders are rejected.

`singularityintegrate(f,x)` integrates a single singularity factor, independent
multipliers, polynomial multipliers and sums. Multiplication by delta
derivatives uses the distributional product rule. Unrelated or unsupported
integrands return `None`, as this is a specialized integration helper.
General products of distributions are not defined by this interface.

## Polynomial and integer systems

`factor_system` decomposes rational polynomial equations into combinations of
irreducible factors, discards inconsistent branches using native Gröbner bases,
and removes duplicate/redundant factor systems. It is not a comprehensive-Gröbner algorithm for stratifying external symbolic parameters.

`solve_triangulated` uses the native exact polynomial system solver, including
lexicographic elimination. Positive-dimensional systems, unsupported radical
solutions and algebraic extension towers raise explicitly; it never presents a
partial solution list as complete.

`substitution` combines provided solution dictionaries with native symbolic
solving, resolves chained values, checks residuals and excludes vanishing
specified denominators. It returns a `FiniteSet` of ordered native tuples.
It inherits native solver limits; full transcendental set-valued parameter
families from SymPy's substitution machinery are not implemented.

`diophantine` supplies complete parametrizations for rational linear equations
in any number of integer variables, integer univariate roots, factored supported
families, binary product equations, positive diagonal binary quadratic equations,
normalized Pell equations `x**2-D*y**2=1` and Pythagorean triples. Linear
parametrizations use unimodular integer column transformations. Finite positive
quadratic cases use bounds derived by completing squares. Pell solutions use
the continued-fraction fundamental unit. Returned parameters carry integer
assumptions. General indefinite quadratics, generalized Pell right-hand sides,
arbitrary higher-degree integer equations and SymPy's complete Diophantine
helper API remain unsupported and raise `NotImplementedError`.
