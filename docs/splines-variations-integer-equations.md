# Splines, variational calculus and integer equations

[README](../README.md)

## Piecewise expressions and splines

`piecewise_fold(expr, evaluate=True)` is SymPy's algorithm: it moves nested
`Piecewise` expressions outward through the surrounding operation and collects
branches with identical conditions. Conditions on a bound variable stay inside
`Sum` and `Integral`. The fold runs natively in Rust
(`crates/symox-core/src/simptools/piecewise_fold_public.rs`); the Python
version of the same algorithm is used when global evaluation is off or one of
the stock constructors it relies on has been replaced. `piecewise_exclusive`
is SymPy's: each condition excludes the earlier ones and is simplified with
`Boolean.simplify`, so the resulting conditions are disjoint. It appends a
`nan` branch for uncovered points unless `skip_nan=True`; `deep=False` limits
the transformation to the outer expression. Boolean simplification follows the
native logic engine's supported reductions.

`bspline_basis` and `bspline_basis_set` use the Cox-de Boor recurrence over exact
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
symbols, then applies the Euler-Lagrange total derivatives. It supports several
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

`SingularityFunction(x,a,n)` uses SymPy's constructor: it represents the
truncated power for nonnegative orders, including non-integer ones, and Dirac
delta derivatives for the orders -1 to -4. Orders below -4, imaginary orders
and arguments `x - a` known to be nonreal raise `ValueError`, and `nan`
propagates. It supports exact point evaluation, substitution, argument
differentiation and rewrites to `Piecewise`, `Heaviside` or `DiracDelta`
anywhere in an expression. The order-zero truncated power uses the
right-continuous value one at its origin.

`singularityintegrate(f,x)` integrates a single singularity factor of integral
order, independent
multipliers, polynomial multipliers and sums. Multiplication by delta
derivatives uses the distributional product rule. Unrelated or unsupported
integrands return `None`, as this is a specialized integration helper.
General products of distributions are not defined by this interface.

## Polynomial and integer systems

`factor_system`, `factor_system_cond`, `factor_system_bool` and
`factor_system_poly` follow SymPy 1.14's algorithm
(`python/symox/_public_solvers.py`): each equation is factored, one
irreducible factor is chosen from every equation, and systems that strictly
contain another are removed. `factor_system` drops the degenerate systems that
only hold under a condition on constant (parameter) factors;
`factor_system_cond` keeps them and `factor_system_bool` returns the result as
an `Or` of `And`s of equations. The decomposition does not check the
consistency of each system and is not a comprehensive-Gröbner stratification
of symbolic parameters.

`solve_triangulated(polys, *gens, extension=False)` is SymPy's Gianni-Mora
triangulated solver over a lexicographic Gröbner basis. As in SymPy, without
`extension` only roots in the ground domain are returned, so the list is not
the complete complex solution set (`[x**2 + y**2 - 1, x - y]` gives `[]`).
With `extension=True` every complex root is lifted through algebraic fields
built from `CRootOf` generators.

`substitution` (`symox.solvers.solveset.substitution`) is a port of SymPy
1.14's substitution method (`python/symox/solvers/_substitution.py`): each
equation is solved for one unknown with `solveset_real` and
`solveset_complex` and substituted into the others. Periodic solutions are
kept as `ImageSet` families, and unsolvable systems give a `ConditionSet`.
It inherits the limits of the native `solveset`.

`diophantine` applies SymPy's preprocessing (the numerator, with factors
independent of the variables removed) and supplies complete parametrizations
for rational linear equations in any number of integer variables (SymPy's
`Linear.solve` form and parameter signs), integer univariate roots, factored
supported families, binary product equations, parabolic binary quadratics
(`b**2 == 4*a*c`, SymPy's case), positive diagonal binary quadratic equations,
normalized Pell equations `x**2-D*y**2=1`, Pythagorean triples and SymPy's
sums of equal even powers `x_1**e + ... + x_k**e = n` with even `e > 3`.
Finite positive quadratic cases use bounds derived by completing squares.
Pell solutions use the continued-fraction fundamental unit. Returned
parameters carry integer assumptions. General indefinite quadratics,
generalized Pell right-hand sides, homogeneous ternary quadratics (where
SymPy's two-parameter answer is not the complete solution set), arbitrary
higher-degree integer equations and SymPy's complete Diophantine helper API
remain unsupported and raise `NotImplementedError`.
