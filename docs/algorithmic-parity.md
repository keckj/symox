# Algorithmic parity

[README](../README.md)

SymPy 1.14.0 is the current comparison reference. Mathematical capability at
least equivalent to that reference is a project requirement. Import availability
and benchmarks on selected examples do not establish this equivalence.

The [versioned API audit](sympy-parity-audit.md) records the traversal of all
310 module paths, including constructor and public method signatures. Run its
resumable driver rather than treating name coverage as behavioral evidence.

## Acceptance criteria

For each advertised operation, compare accepted domains and options, the meaning
of the result, exactness or requested numerical accuracy, and failure behavior.
Solver comparisons must compare the entire returned solution set, including
real-domain filtering, multiplicities where appropriate, excluded poles and
branch conventions. Checking that each returned solution satisfies an equation
does not detect missing solutions.

An unsupported method or exhausted resource budget must remain distinguishable
from a proof of no solutions. An undecidable exact residual must not be rejected
solely because evaluation in machine precision produces a nonzero value.

Complexity comparisons must include input size, coefficient size, expression
growth and output size. Rust implementation and faster timings alone do not
prove a better asymptotic algorithm. Performance optimization must preserve the
existing mathematical result and fallback behavior.

## Intentional corrections validated on 2026-10-22

The user validated these results, where SymPy 1.14 is mathematically wrong.
Symox keeps the correct answer and does not reproduce the reference. Each
case is tested with the symox value and detailed in its section below.

- Periodic inequalities on a bounded domain:
  `solveset(sin(x) > 0, x, Interval(2*pi, 3*pi))` is
  `Interval.open(2*pi, 3*pi)` (SymPy `EmptySet`), and
  `solveset(sin(x) >= 0, x, Interval(0, 2*pi))` includes `2*pi`.
- `integrate(Piecewise((1, sin(x) > 0), (0, True)), (x, 0, 4*pi))` is
  `2*pi` (SymPy `pi`); `simplify` keeps the condition `sin(x) > 0` and
  `simplify_univariate` keeps the point `x = 0`.
- `is_increasing(a*x, symbol=x)` with a real `a` is `None` (SymPy `False`,
  wrong whenever `a >= 0`).
- Noncommutative `solve`: `solve(B*A - C, B)` is `[C*A**-1]` (SymPy
  `[A**-1*C]`), `solve(A*B*C - 1, B)` is `[(C*A)**-1]` (SymPy
  `[(A*C)**-1]`), `solve(A*B - B, B)` is `[0]` (SymPy `[]`), and
  `solve(A*B - B*A, A)` raises `NotImplementedError` instead of SymPy's
  `[]`, since `A = B` is a solution. `solve(x**2*A - A, x)` and
  `integrate(A*x, x)` are computed where SymPy raises.
- Noncommutative `nonlinsolve` follows `solve`: `nonlinsolve([A*B - C], [B])`
  is `{(A**-1*C,)}` (SymPy `{(C*A**-1,)}`, which does not satisfy the
  equation), `nonlinsolve([A*x - B], [x])` is `EmptySet` (SymPy gives the
  commutative `x` the noncommutative value `B*A**-1`),
  `nonlinsolve([A*B - B*A], [A])` raises `NotImplementedError` (SymPy
  `{(Complexes,)}`) and systems of several equations raise
  `NotImplementedError` (SymPy crashes with `AttributeError`).
- `srepr` of an applied undefined function with assumptions keeps its facts
  (`Function('g', ..., commutative=False, ...)(Symbol('x'))`) so that it
  round-trips through `sympify`, pickling and `to_sympy`; SymPy prints the
  lossy `Function('g')(Symbol('x'))`.
- For an integer unknown `x**5 - x - 1` has no solution, and
  `x**3*exp(x) - 1` returns three valid Lambert representatives.

Decided alignment: `2nd_power_series_ordinary` follows SymPy and refuses
equations with a right-hand side (`ValueError`), so `classify_ode` lists the
same hints as the reference.

## Wave-4 intentional corrections pending validation, 2026-10-22

In these 15 cases SymPy 1.14 is wrong and symox keeps the correct value.
They still need the user's validation. Each case was checked against a
definition, mpmath at 30 to 40 digits, an equilibrium condition or a
normalization integral, not against the reference. The tests are in
`tests/test_wave4_h_validation_parity.py`.

- `DiagonalOf(MatrixSymbol('x', n, 3)).diagonal_length` is `Min(3, n)`
  for an integer `n`. SymPy returns 3, which is wrong for `n < 3`. A
  2 x 3 matrix has a diagonal of length 2.
- `hn1(nu, z).rewrite(hankel1)` and `hn2(nu, z).rewrite(hankel2)` use the
  order `nu + 1/2`, since `h1_n(z) = sqrt(pi/(2*z))*H1_{n+1/2}(z)`.
  At `nu = 1, z = 2`, the definition `jn + I*yn`, the closed form
  `-exp(I*z)*(z + I)/z**2` and the symox rewrite all give
  `0.43539777... - 0.35061200...*I`. SymPy's `hankel1(nu, z)` form gives
  `0.5111 - 0.0949*I`. For `hn2(2, 7/5)` the true value is
  `0.11334028 + 1.57276046*I`, and SymPy's form gives `0.2196 + 1.0830*I`.
- `Cable.solve()` with loads at fixed points `Z(2, 7.26)` (3 downward) and
  `X(4, 6)` (8 downward) between `A(0, 10)` and `B(10, 10)` raises
  `ValueError`. With vertical loads only, the horizontal tension is
  constant. Equilibrium at `Z` needs `H = 3/(1.37 - 0.63) = 4.054`, and
  equilibrium at `X` needs `H = 8/(0.63 + 2/3) = 6.17`, so no
  equilibrium exists. SymPy returns tensions whose horizontal components
  are 5.255, 4.054 and 5.700 in the three segments. `plot_tension`
  follows.
- `Cable(('A', 0, 40), ('B', 100, 20))` with a distributed load of 850 and
  `solve(58.58)` gives `R_A = (-36465, 49793)` and `R_B = (36465, 35207)`.
  The parabola gives `H = w*(x0**2 - (100 - x0)**2)/40 = 36465` and
  vertical reactions `w*x0` and `w*(100 - x0)`, which sum to the load
  85000. SymPy's reactions are not balanced: horizontally 36465 against
  44400, vertically 92661 against 85000.
- `gaussian_conj(3, 2, 5)` is `(-5/4, 25/4, 5*sqrt(2)/4)`. ABCD
  propagation of `q = s + I*z_r` through a thin lens gives
  `z_r_out = 6.25` and `m = sqrt(z_r_out/z_r) = 1.7678`, which matches
  `m = 1/sqrt((1 - s/f)**2 + (z_r/f)**2)`. SymPy uses `1 - (s/f)**2` and
  gives `z_r_out = 5/2` and `m = sqrt(5)/2`. `s_out` is the same in both.
- CDF support branches. Each symox CDF matches the integral of the
  density inside the support to 1e-20 and is 0 or 1 outside it:
  - `cdf(BetaNoncentral('x', 1, 1, 1), evaluate=False)(2).doit()` is 1.
    SymPy gives `2*exp(1/2) = 3.297`.
  - `BoundedPareto` has a `x >= H` branch. SymPy gives 1.11375 at
    `L = 1, H = 3, x = 10`.
  - `PowerFunction` has a `z >= b` branch. SymPy gives 1.7778 at
    `a = 1, b = 4, z = 5`.
  - `Reciprocal` has a `x >= b` branch. SymPy gives 1.209 at
    `a = 1, b = 5, x = 7`.
  - `Levy` is 0 for `z < mu`. SymPy's `erfc` form is complex there:
    `1 - 1.6504*I` at `mu = 1, c = 2, z = 0`.
  - `LogitNormal` is 0 for `z < 0` and 1 for `z >= 1`. SymPy gives
    `15.16 - 4.63*I` at `z = 2`.
- `marginal_distribution(MultivariateBeta('C', a1, a2), C[0])(x)` is the
  `Beta(a1, a2)` density `x**(a1 - 1)*(1 - x)**(a2 - 1)/beta(a1, a2)`.
  SymPy integrates `x2` over `[0, 1]` instead of the simplex. Its result
  integrates to 2 at `a1 = 2, a2 = 3`.
- `density(GaussianUnitaryEnsemble('U', 2))(X)` is
  `2*exp(-Trace(X**2))/pi**2`. Over the four real parameters,
  `Tr H**2 = h11**2 + h22**2 + 2*(a**2 + b**2)`, and
  `exp(-Tr H**2)` integrates to `pi*(pi/2) = pi**2/2`. SymPy's
  `1/(2*pi**2)` integrates to 1/4.
- `ImplicitRegion((x, y, z), (x + y + z)**4).regular_point()` raises
  `NotImplementedError`. Every partial derivative is `4*(x + y + z)**3`,
  which vanishes everywhere on the surface, so no regular point exists.
  SymPy returns `(-10, -10, 20)`, which is a singular point. The circle
  case returns a valid point.

## Solver work on 2026-10-18

The native univariate `solve` pipeline now falls back from radical formulae to
exact indexed root objects for rational-coefficient polynomials. Factorization
separates repeated factors; real-root intervals determine real/complex and sign
facts. Polynomial residuals and polynomial poles at a root object can be checked
by exact reduction modulo its defining polynomial.

Native Lambert W dispatch covers affine logarithmic equations,
affine exponentials plus linear terms, and positive integral powers of affine
factors multiplied by affine exponentials. Rational interval arithmetic selects
the additional real `W_-1` branch when its domain can be established. Like
SymPy's `solve`, this is a finite principal-branch convention, not a claim to
enumerate every complex Lambert branch.
For these forms, complex `solveset` requests and unresolved parameter cases
retain a `ConditionSet` rather than presenting finite representatives as a
complete set. Numeric real cases use the certified real-branch selection.

The audit records two intentional mathematical improvements separately from
reference matches. For an integer unknown, `x**5-x-1` has no solution: its value
is odd in both residue classes modulo 2. SymPy 1.14.0 retains a nonintegral root
object in this case. For `x**3*exp(x)-1`, symox returns three valid Lambert
representatives where the reference returns one. Their validity follows from
`W(z)*exp(W(z))=z` and the exact three cube roots of `1/27`.

Elementary trigonometric inversion now retains complete complex integer
lattices for affine arguments, including cotangent, secant and cosecant.
Real-domain restriction solves the imaginary-part equation for the integer
index and preserves an exact intersection when it cannot decide. Real range
bounds apply only to a provably real affine argument. For example,
`solveset(exp(I*x)-1, x, S.Reals)` retains every `2*pi*n` solution.
The [scoped solver review](../audits/sympy-1.14.0/solver-review.json) records
these invariants, exact whole-lattice comparisons and remaining complexity
obligations. It separately records finite-root corrections for reference
results containing `nan` or infinity; those are not reference-output matches.

Lambert W numerical evaluation uses native multiprecision Halley iteration,
including complex integer branches. Exact candidate checking tries symbolic
proof first; a numerical rejection requires a nonzero residual stable at 160
and 320 bits, above the equation-scaled rounding threshold. Undecidable residuals
are retained. Approximate candidate checking uses the declared input precision
and an equation-relative residual scale.

Polynomial `solve` and `roots` now honor cubic/quartic radical controls, including
independently enabled quartics. Disabling formulas in `solve` retains exact
indexed roots for rational-coefficient polynomials. `roots(strict=True)` checks
completeness after filters and predicates and raises when its retained radical
subset is incomplete. The univariate `solve` flag `incomplete` (default `True`)
returns the partial `roots` result when `all_roots` does not support the
coefficient domain, and `incomplete=False` raises `NotImplementedError` there;
systems always ask for complete roots. These options do not establish all
solver-option parity.

Multivariate polynomial product equations retain the union of all factor zero
sets, with every other equation applied to each branch. This prevents pivot
division from discarding a whole component. Native factorization can recover
expanded products before a potentially singular pivot, and a proved nonzero
linear pivot is preferred where available. For example, `x*y-z=0` has the full
polynomial parametrization `z=x*y`, including zero values of either parameter.
Branch exploration can grow as the product of factor counts before duplicate
elimination; comparative solver complexity and general parameter strata remain
open. Inferred/unordered unknowns use mapping output, while ordered explicit
unknowns retain tuple columns, following the reference interface conventions.

Scalar `solve` also attempts rational-function coefficient identities before
ordinary algebraic solving, preserving the distinction between a bare equation
and a one-equation list. Full parameter strata remain unproved.

Overdetermined systems follow SymPy's `_solve_system` (2026-10-22). As in
`solve`, an equation with free symbols but none of the unknowns is dropped
(`solve([x*y - 2, y - 1], [x])` is `{x: 2/y}`), unless its assumptions prove
it nonzero, and the Float rationalization flag only looks at the kept
equations (`solve([x - 1, 1.5*y - 2], [x], dict=True)` is `[{x: 1}]`). An
unknown that is not a symbol (`f(x)`, `Derivative(f(x), x)`) counts through
`has`, as SymPy swaps it for a Dummy before this test; otherwise the Floats
of `classify_ode` equations reached the native solver and the inverse
fallback recursed on `1.0*(...)` factorizations. A
system linear or polynomial in the unknowns is solved over the coefficient
domain of its parameters: a residual equation free of the unknowns whose
expanded numerator is nonzero makes the system generically inconsistent, so
`solve([1.5*x*y - 2, x + y - 1], [x])`, `solve([x - y, x - 1], [x])` and
`solve([x - 1, x - sin(y)], [x])` are `[]`, as in the reference, instead of
keeping a candidate valid only on a parameter slice. For a system that is not
polynomial in the unknowns, the residual is rejected only when it has no
function and no radical (SymPy's `checksol`), so `[exp(x) + x - y, x - 1]` is
`[]`, while an undecided residual such as `E - sin(y)` keeps the candidate.
Known differences: SymPy treats `sqrt(a)` as a generator independent of `a`
and returns `[]` for `solve([x**2 - a, x - sqrt(a)], [x])`, while symox keeps
the valid `[(sqrt(a),)]` (divergence not yet validated by the user). The
non-polynomial residual rule was compared on a few cases only (`[sin(x) - y,
x - 1]` keeps `x = 1` in both). SymPy's `_invert` preprocessing of each
equation is ported: independent terms move to the right-hand side, two-term
sums of powers are combined by `powsimp`, and `exp`, `log`, `tan`, `asin`,
`acos`, `atan`, negative integer powers and radicals of a base with more
than one unknown are inverted (not `sin` or `cos`, like the reference). An
equation whose inverted left side is an unknown in which it is not linear
gives the tuple output of the reference instead of a dictionary
(`solve([exp(2*x) - 1, y - 1], [x, y])` is `[(0, 1)]`, while
`solve([exp(x + y) - 2, x - y], [x, y])` is `{x: log(2)/2, y: log(2)/2}`).
With `rational=False`, a connected component of the system (equations
sharing unknowns) that is polynomial and has a Float coefficient keeps
Float values, like the reference RR domain: `solve([0.25*x**2 - 1, x - 2],
[x], rational=False)` is `[(2.00000000000000,)]`. Dictionary keys and the
rows of a solution list follow the reference canonical order (`ordered`
keys, then rows sorted by `default_sort_key`). A single expression target such as `f(x)`, `exp(x)` or `x**2` can be
recast as a fresh algebraic atom only with rational work equations and exact
reconstruction of the original. Original targets and derivative coefficients
are restored in the output; these are algebraic solving conventions, rather
than an assertion about all analytic ranges or differential equations.
Relational `solve` inputs now reach the native inequality engine with every
supplied condition retained. In particular, a mixed equation/inequality system
cannot silently discard the inequality. This inherits that engine's supported
domains and remaining multivariate restrictions.

Native complex square roots recover their small component from the product
identity, preserving genuine tiny imaginary inputs near the negative real axis.
Exact Cardano-pair norm/sign proofs prevent imaginary rounding noise from
changing nested square-root branches; real and imaginary cancellation are
tracked separately. Reciprocal rational-power product pairs also have exact
conjugacy proofs when their squared moduli are rational. This fixes held
numerical evaluations of real cubic roots and lost quartic solutions in the
checked solver.
The [scoped radical review](../audits/sympy-1.14.0/radical-review.json) records
the invariants and independent numerical regressions. General accuracy bounds,
precision options and comparative bit complexity remain separate obligations.

For ordinary scalar symbol substitutions, `evalf(subs=...)` now holds arithmetic
until the multiprecision evaluator runs. This prevents premature floating-point
cancellation such as `x+y-z` with values `1e16`, `1`, `1e16` becoming zero.
Sequential substitution semantics, `.n()`, `N` and scalar matrix entries use
the same native path. Binding expressions, algebraic substitution keys, Boolean
control flow and the Python special-function adapter retain their existing paths;
the change does not certify those paths or all numerical accuracy options.

Run the focused comparison audit with both libraries installed:

```sh
python benches/audit_solver_parity.py --out /tmp/solver-parity.json
```

This audit covers the cases it names. It is not a global parity certificate.

## Algebraic fields and composite domains on 2026-10-21

Algebraic number fields are no longer limited to quadratic, cubic and the
previously supported towers. `QQ.algebraic_field`, `QQ.alg_field_from_poly`,
`QQ.cyclotomic_field` and `extension=` accept any finite list of algebraic
generators, including fourth roots such as `2**(1/4)`, nested radicals,
`CRootOf` instances and arbitrary irreducible polynomials. Minimal polynomials
follow the reference `_minpoly_compose` composition; when the composition
leaves several candidate factors, the factor is selected by a certified
outward-rounded enclosure rather than by a floating-point tolerance. Primitive
elements, their coordinates and the embedding of each generator are computed
exactly. Membership (`to_number_field`, `field_isomorphism`, conversion into a
field) is decided by factoring the minimal polynomial over the field and
identifying the matching factor with a certified enclosure; an element that
cannot be identified is reported as such and never coerced to an approximate
value.

The field elements propagate through the polynomial algorithms. Univariate and
multivariate `gcd`, `sqf_list` and Trager `factor` over number fields are
line-by-line ports of `dmp_ff_prs_gcd`, `dmp_sqf_list`, `dmp_sqf_norm` and
`dmp_ext_factor`. Division, `gcdex`/`half_gcdex`/`invert` (also over
`QQ(a)`), `cancel` and `groebner` use the same exact arithmetic, with the
reference domain-error messages. For example, `factor(x**4 - 2*y**4,
extension=2**(1/4))` and `gcd(x**4 - 2, x - 2**(1/4), extension=True)` match
SymPy.

`ZZ_I`, `QQ_I` and `EX` are now polynomial coefficient backends with the
reference ring algorithms (Gaussian factorization, EX dense `dmp` algorithms
including `gcd`, `sqf`, `roots` and `nroots`). Composite domain inference
selects them as SymPy does: `Poly(x**2 - sqrt(2)*y, x, extension=True)` is
over `EX`, `Poly(x + y).eval(x, I)` is over `ZZ_I`, and transcendental
coefficients such as `ZZ[sin(a)]` or `QQ(pi)` support `sturm`, `count_roots`,
`groebner`, `factor` and `sqf` of a `Poly`. Unification keeps the receiver's
generator order and applies the reference ground rule: a fraction field
unified with a polynomial ring over a non-field ground is taken over the ring,
so `QQ(a)` with `ZZ[pi]` gives `ZZ(a,pi)` for `Domain.unify` and `Poly`
arithmetic alike.

Complex roots are isolated with the reference certified rectangles, and
`CRootOf` orders complex roots exactly and evaluates by rounded refinement
(commits `Port the reference root isolation with certified complex rectangles`
and `Port CRootOf with exact complex ordering and rounded evaluation`).

Differential tests: `tests/test_composite_domain_parity.py`,
`test_general_algebraic_fields.py`, `test_degree_six_eight_fields.py`,
`test_algebraic_fields.py` and `test_implicit_domain_rings.py`. SageMath's
`sage/interfaces/sympy.py` doctests now have 3 failures out of 319 (311
passing before this work): the 2 `ZZ[x][y]` cases and the Fresnel integral.
The `ZZ[x][y]` `_sage_` path (`domain`, `gens`, class hook) now matches the
reference; the doctests were not rerun in Sage for this change.

Remaining gaps in this area:

- Nested polynomial rings (`PolynomialRing(ZZ[x], 'y')`, `ZZ[x][y]`) and
  ring factories over `ZZ_I`, `QQ_I` and algebraic fields (`ZZ_I[x]`,
  `QQ_I[x]`, `QQ<sqrt(2)>[x]`) are implemented with the reference
  construction, element arithmetic, conversion, printing and `unify`, and
  `Poly` over nested rings such as `ZZ[x][y]` (differential tests in
  `tests/test_nested_rings_parity.py`). Remaining differences: domains and
  elements are `Domain`/`DomainElement` instances (they pass `isinstance`
  checks against `PolynomialRing` but `type(...).__name__` differs, and
  `.ring` and `has_assoc_Field` are missing); invalid elements raise
  `CoercionFailed` where the reference raises `ValueError`;
  `ZZ[x][y].convert(ZZ[x](x))` and `set_domain` between `ZZ[x][y]` and
  `ZZ[x,y]` convert where the reference raises `CoercionFailed`;
  `Poly` over `QQ<sqrt(2)>[x]` raises `NotImplementedError`; expression-level `factor` or `gcd` with a
  polynomial-ring `domain=` (e.g. `ZZ[x]`, `ZZ[x][y]`) raise
  `NotImplementedError` (the `Poly` methods work); adding a `ZZ[x][y]` and a
  `ZZ[x]` polynomial raises `TypeError` instead of `GeneratorsError`;
  `QQ<sqrt(2)>(sqrt(2))` is accepted where the reference rejects it.
- Gaussian composite grounds (differential tests in
  `tests/test_gaussian_domains_parity.py`, 148 cases with hardcoded SymPy
  1.14 outputs). `Poly` inference selects `ZZ_I[y]`, `QQ_I[y]`, `ZZ_I(y)`
  and `ZZ_I(N,pi)` as the reference does (`Poly(I*x*y, x)`,
  `Poly(I*x/2*y, x)`, `Poly(I*x/y, x)`, `Poly(-2*I*pi*n/N, n)`), with
  `domain='ZZ_I[y]'`, `domain=ZZ_I.frac_field(y)`, `Domain.unify` and
  `construct_domain` for ring grounds. A `Poly` over these domains keeps its
  coefficients in the reference ground normal form and matches the
  reference for printing (`(0 + 1*I)*y`), arithmetic, `div`/`rem`/`quo`,
  `gcd`, `lcm`, `cofactors` (methods and the top-level functions),
  `factor_list`, `sqf_list`, `sqf_part`, `content`, `primitive`, `monic`,
  `clear_denoms`, `invert`, `gcdex`, `cancel`, `is_irreducible` and
  `is_sqf`. `ZZ_I[x].is_unit(I)` is `False`, nested `is_unit` follows the
  reference (`a.is_ground` and a unit leading coefficient, only `+-1` in a
  non-field ground ring) and `canonical_unit` normalizes by the ground unit
  of the leading coefficient; a `domain='ZZ_I(y)'` string raises the
  reference `OptionError`, and `set_domain` from `ZZ_I[y]` to `ZZ[y]` raises
  the reference `CoercionFailed`. The fraction fields `ZZ_I(y)`,
  `QQ_I(y)` and `ZZ_I(N,pi)` have domain elements
  (`tests/test_gaussian_fractions_parity.py`): conversion, `construct_domain`
  (`[I/y]`, `[-2*I*pi/N]`, `[I*y]` with `field=True`), arithmetic, the
  reference cancellation to coprime numerator and denominator with integral
  Gaussian coefficients and a canonical denominator unit, the raw
  `FracElement` powers (`K(I/y)**-1` is `y/(0 + 1*I)`), `gcd`, `lcm`,
  `cofactors`, `numer`/`denom`, `is_unit`, signs and `canonical_unit`,
  printing (`(0 + 1*I)*y/(1 + 0*I)`) and the reference errors
  (`ValueError` for a non-rational expression, `ZeroDivisionError`);
  `quo_ground`/`exquo_ground` failures print the reference message. Remaining
  gaps: `sqf_list` of a multivariate polynomial over a composite Gaussian
  ground raises `ValueError` (not implemented); the SymPy 1.14 result is
  itself wrong in some of these cases (`Poly((x + I*y)**2*(x - z), x, z)`
  gives `-1, [(x - z, 1), (x + I*y, 2)]`, whose product is the negation of
  the input, and similarly over `QQ_I[y]`), so it is not reproduced;
  multivariate `sqf_list` over `ZZ_I` such as `I*(x + y)**3*(x - I)` raises
  the factorization certificate error, where the reference result
  `(1, [(x - I, 1), (I*x + I*y, 3)])` is also off by the unit `-1`;
  `K.convert(ZZ_I(1 + I))` succeeds where the reference raises
  `CoercionFailed: expected an integer`; `ZZ_I(y).from_sympy(Float(1.5))`
  raises `ValueError` where the reference gives `(3 + 0*I)/(2 + 0*I)`; elements are `DomainElement`
  objects, not `FracElement`; fraction fields over algebraic fields still
  have no domain elements. `.rep` prints the reference `DMP_Python` over
  `ZZ_I[y]` or `ZZ_I(y)`, but the representation is expression-based (algorithms run on a flattened Gaussian polynomial),
  so it is slower than the reference dense recursive one. The pre-existing
  `ZZ[x].canonical_unit(ZZ[x](-x))` repr is `ZZ[x](-1)` (reference `-1`).
- A `Poly` whose inferred coefficient domain is a fraction field over an
  algebraic ground is built over `EX`: the coefficients and the
  arithmetic stay exact, but the domain and domain-specific operations
  (division, `gcd`, `factor`) differ from the reference.
- Number factors of an unevaluated product are coefficients, as in the
  reference `_dict_from_expr`: `Poly(Mul(-1, 2, x, evaluate=False), x)` is
  `Poly(-2*x, x, domain='ZZ')`.
- `resultant`, `subresultants`, `discriminant`, `gcd`, `lcm` and `cancel`
  with `extension=True` and parameter coefficients follow the reference: the
  option is ignored when a coefficient is not an algebraic number, and
  `parallel_poly_from_expr` infers one joint domain (`EX` for
  `[x**2 - sqrt(2)*y, x - y]`). Only error messages differ for an explicit
  `extension=sqrt(2)` with parameters.
- `RR.invert` reduces `1/a` modulo `b` with the reference `mpf_mod` rule
  (`RR.invert(4, 0)` is `0.25`, `RR.invert(-2, 0)` raises
  `ZeroDivisionError`); `Domain.convert` rejects strings, `None`, bytes and
  sequences with `CoercionFailed`, and `nan`/`oo`/`zoo` outside `EXRAW`.
  Python floats, `mpf`, `complex` and `mpc` take the reference route
  through `RR`/`CC` (53 bits): `QQ.convert(2.5)` is `5/2` and
  `QQ.convert(0.1)` is `1/10` (the limited denominator of
  `RR.to_rational`), `ZZ`, `GF(p)` and `ZZ[x]` accept integer-valued
  values only, `QQ_I` takes `CC` values and `ZZ_I` refuses both, algebraic
  fields read `RR` through `QQ`, and the refusals carry the reference
  message (`Cannot convert 2.5 of type <class 'mpmath.ctx_mp_python.mpf'>
  from RR to ZZ`). The SymPy `Float` route of `from_sympy` is the
  reference one too (`ZZ.convert(Float(3))` is `3`, `QQ_I.convert(Float(2.5))`
  is `5/2`). `ZZ(2.5)` truncates and `QQ(0.1)` is the exact `Fraction`, as
  the reference constructors; `QQ(1j)` and `QQ(mpf(2.5))` raise its
  `TypeError`. The `from_sympy` refusals of `ZZ`, `QQ`, `RR`, `CC`, `ZZ_I`
  and `QQ_I` use the reference messages (`expected an integer, got pi`,
  `sqrt(2) is not Gaussian`), including the reference `AttributeError` for
  sets and dicts in `RR`, `CC` and the Gaussian domains; `EX` sympifies a
  set or dict and `EXRAW` refuses non-`Expr` objects (`Interval`,
  `FiniteSet`) with `Expecting an Expr instance but found: ...`.
  Univariate `Poly` objects over `RR[y]`, `CC[y]` and `RR(y)` are inferred
  and computed like the reference (see `docs/domains.md`, composite domains
  over `RR` and `CC`, for the ported algorithms and the remaining gaps:
  several generators, `factor`, `ground_roots`, `apart` and the `to_ring`
  message).
  Ground elements print like the reference raw values (`2`, `MPQ(5,2)`,
  `mpf('2.5')`, `mpc(real='1.0', imag='2.0')`, in the field's mpmath
  context), and elements of `RR[x]`, `CC[x]`, `RR(x)` or `RealField(100)[x]`
  like the reference `PolyElement`/`FracElement` (`-x - 1.5`,
  `(0.0 + 1.0j)*x + (1.0 + 0.0j)`, `1.0/(x + 1.5)`); `CC.invert` raises the
  reference `TypeError` (`unsupported operand type(s) for %: 'mpc' and
  'mpc'`); `Domain.unify` with an inexact ground follows the reference
  composite rule (`ZZ[x].unify(RR)` is `RR[x]`, `RR[x].unify(CC)` is
  `CC[x]`), `RR[x].is_Exact` is False and `RR[x].get_exact()` is `QQ[x]`.
  `parallel_poly_from_expr` infers one joint domain when the separately
  inferred ones differ (`EX` for `[x + exp(y), x*y]`, `RR[y]` for
  `[x*y, x + 1.5]`). Tests: `tests/test_domains_inexact_wave3_parity.py`
  (105 of 131 probed calls identical, 53 before).
  Remaining: `ZZ.convert_from(QQ(1, 2), QQ)` keeps the message `expected an
  integer`; algebraic field elements print as `QQ<sqrt(2)>(sqrt(2))` (reference
  `ANP(...)`); domain elements are `DomainElement` (reference `PolyElement`)
  and `RR.frac_field(x).field` is missing; `EXRAW` elements are wrapped and
  do not compare equal to the raw `Expr`. Floats over several generators
  report the reference domain: `Poly(1.5*x*z + y, x, z).domain` is `RR[y]`,
  `I*x*z + 1.5*y` gives `CC[y]` (and `1.5*x + I` gives `CC`),
  `1.5*x*z + 1/y` builds `RR(y)` with native ring arithmetic, and `monic`/
  `to_field` move to `RR(y)`. `ground_roots` reads the linear factors of
  `factor_list` like the reference (`{1.0: 1, -1.0: 1}` for
  `1.5*x**2 - 1.5`) and raises the reference `TypeError` over `RR[y]`, as
  does `apart` over `RR[y]` once the cancelled denominator has degree above
  one. Tests: `tests/test_domains_inexact_wave3_parity.py` (135 of 149
  probed calls identical; 105 of the first 131 before).
  Plain univariate `RR` `invert`/`gcdex`/`half_gcdex` (methods and
  expression functions) run the reference `dup_half_gcdex` in 53-bit
  floating point (`dup_ff_div` drops the rounding-created leading term,
  `dup_gcdex` takes `dup_quo`, `invert` needs `h` exactly `1.0`), so
  `invert(1.5*x, x**2 + 1)` is `-0.666666666666667*x` and
  `invert(x**2 - 1.0, x - 1)` raises `NotInvertible`. `Poly.cancel` over
  plain `RR` computes the cofactors over QQ like `dup_inner_gcd` and
  applies the `dmp_cancel` field sign rule, returning a Float `cp/cq`
  (`include=True` keeps the sign in `p`). `content`/`primitive` over plain
  `RR` give `1.0` like the reference field content (also several
  generators and the `content`/`primitive` functions). `PolynomialRing.ring`
  and `FractionField.field` exist (`RR.poly_ring(x).ring`,
  `CC.frac_field(x).field`). Tests: `tests/test_domains_rrcc_parity.py`
  (82 of 157 probed calls identical; 49 before).
  Remaining: products, `gcd`, `lcm`, `div` and `resultant` of a
  several-generator `RR[y]` polynomial with an exact `ZZ[y]` one fail
  (`TypeError`; the reference unifies to `RR[y]`), and several-generator
  `RR[y]` polynomials keep the legacy representation, whose
  `factor_list`/`sqf_list`/`content`/`integrate` differ from the reference
  (`Poly(1.5*x*z + 1.5*y, x, z).content()` is `3/2`, reference `1.0`).
  `apart` over RR (`1.5/(x**2 - 1)`, reference
  `-0.75/(x + 1) + 0.75/(x - 1)`, or `1/(x - 1.5*y)`) computes exactly over
  QQ and returns rationals instead of Floats. `CC` has no polynomial
  backend for an explicit `domain=CC` (`Poly(x, x, domain=CC)` raises
  `NotImplementedError`), and the inferred `CC` polynomials have no
  `invert`/`gcdex` (`ValueError` where the reference returns the CC result
  or raises `NotInvertible`). Elements of `RR.poly_ring(x).ring` print
  `0`/`1.0*x` where the reference prints `0.0`/`x`, and the zero `RR` Poly
  prints `Poly(0, x, domain='RR')` (reference `Poly(0.0, ...)`).
  `factor`/`factor_list`/`sqf_list` of RR expressions print factors as
  `1.0*x - 1.0*y` where the reference prints `x - y`, and `factor(1.5*x -
  1.5*y)` does not pull the Float content out as the reference
  `together`/`gcd_terms` does (`1.5*(x - y)`).
- Negative bases with rational exponents follow the reference
  `Integer._eval_power`/`Rational._eval_power` and the `Mul.flatten` sign
  join: `(-2)**(1/6)` stays one power (also after `(-1)**(1/6)*2**(1/6)`),
  `(-2)**(-1/6)` is
  `-(-2)**(5/6)/2`, `(-2/3)**(1/6)` is `(-2)**(1/6)*3**(5/6)/3`,
  `(-2)**(-3/2)` is `sqrt(2)*I/4`, `I*(-1)**(1/3)` is `(-1)**(5/6)`, and
  `sqrt(sqrt(2)/2)` is `2**(3/4)/2`. With `extension=True`, `(-2)**(1/6)`
  now builds the reference field `QQ<(-2)**(1/6)>` (minpoly `x**6 + 2`), and
  `factor`/`gcd`/`sqf_list`/`cancel` over it match. Tests:
  `tests/test_extension_root_parity.py` (79 of 86 probed calls identical;
  36 before). Remaining: the reference per-prime extraction of
  `Integer._eval_power` is not ported (`(-12)**(1/6)` stays one power, the
  reference gives `(-3)**(1/6)*2**(1/3)`; `144**(1/3)` is `12**(2/3)`, the
  reference `2*18**(1/3)`); porting it makes `simplify` of the quartic
  root residual of `x**4 + x + 1` exhaust the multivariate Gaussian gcd,
  which must be bounded first. `factor(x**6 + 2, extension=(-2)**(1/6))` exhausts
  the factorization work limit (`ValueError`; the reference factors it in
  about one second), `nsimplify((-2)**(1/6))` returns the power unchanged
  where the reference rewrites it as `2**(1/6)*sqrt(3)/2 + 2**(1/6)*I/2`,
  and `construct_domain` elements print as `QQ<...>(...)` instead of the
  reference `ANP(...)` repr (general, not specific to this field).
- `construct_domain` follows the reference control flow: the options go
  through the `Poly` option rules (an unknown option, `extension=False` or a
  non-boolean `field` raise the reference `OptionError`; `gens`, `domain`
  and `greedy` are accepted), only `extension=True` enables algebraic
  inference (`extension=sqrt(2)` gives `ZZ` for `[1]` and `EX` for
  `[sqrt(2)]`), non-algebraic coefficients with `extension=True` give a
  composite domain (`ZZ[pi]` for `[pi]`) and algebraic numbers among
  composite generators give `EX` (`[pi, sqrt(2)]`). Composite elements print
  like the reference `PolyElement` (`pi`, `x*pi + 1`, `y*(sin(x))`,
  `-x/(y - 1)`). Floats with generators build the reference
  `RealField`/`ComplexField` composite at the largest Float precision
  (`[1.5*x]` gives `RR[x]`, `[1.5*x, pi]` gives `RR[x,pi]`, `[1.5*x/y]`
  gives `RR(x,y)`, `[1.5*x + I]` gives `CC[x]`), also with
  `composite=True`. Remaining: a Float outside a coefficient position (an
  exponent such as `x**1.5`, a function argument) keeps `EX`.
- `EX` root forms and `nroots` digits can differ from the reference. An
  audit of 56 `nroots`/`roots` calls (47 identical) found these intentional
  corrections, where symox returns the more accurate roots: the triple root
  of `(x - 1/3)**3` (reference `NoConvergence`), `x**2 - 2*x + 1 + 1/10**10`
  (`1.0 +- 1.0e-5*I`, reference imaginary parts off by about `1e-18`), the
  double root `-I` of `x**2 + 2*I*x - 1` (reference `~1e-13` real parts), the small root
  `1.0e-20` of `x**2 - 1e20*x + 1` (reference `0`) and the double roots of
  `x**4 - sqrt(2)*x**2 + 1/2` (reference splits them near `1e-8`).
  Float products and powers are expanded in the reference order
  (`Mul._expandsums` halves, `Pow._eval_expand_multinomial` multinomial
  terms, over the canonical argument order where two Floats compare by
  their mpf tuple, positive first), so Float coefficients have the
  reference bits and clustered roots such as those of
  `(x - 1)*(x - 1.0001)*(x - 1.0002)` are as accurate as the reference
  (they were off by `~2e-8` from a last-bit coefficient difference).
  `Rational(Float)` gives the exact binary ratio like the reference. An
  audit of 223 clustered and near-multiple root cases (clusters of 2-5
  roots at distances `1e-2..1e-12`, real and complex, Float and exact
  coefficients, degree up to 40, `n` in 15/30/50/100, `maxsteps` and
  `cleanup`) measured against high-precision roots of the exact
  coefficients finds no case less accurate than the reference.
  `roots` over `EX` factors first like the reference
  (`x**2 - (1 + sqrt(2))*x + sqrt(2)` gives `{1: 1, sqrt(2): 1}`) and
  quadratic roots use `factor_terms(sqrt(d)/A)`. Remaining: the reference
  `to_rational_coeffs` step for an irreducible `EX` polynomial is not
  ported; the terms of `sqrt(I) + sqrt(2)*(-1 + I)/2` (a root of
  `x**2 + (1 - I)*sqrt(2)*x - 2*I`) print in a different order; `Rational(1, 0)`
  raises `ZeroDivisionError` instead of giving `zoo`. Exact reduction of root objects (`CRootOf` and
  `RootSum` simplification into an algebraic field) is not ported.
- `RootSum` follows the reference constructor: the polynomial is kept as an
  expression with its generator (`args == (expr, Lambda, gen)`), pure
  factors are split, linear and quadratic factors of a rational function are
  summed explicitly (`_rational_case`), and it prints and `srepr`s like the
  reference. `doit`, differentiation and `evalf` go through the explicit
  roots as in the reference. `ratint` of a rational function with symbolic
  coefficients (`1/(x**3 + x + a)`) gives the reference `RootSum`. The
  supporting polynomial rules are those of the reference: `composite=False`
  over a polynomial or fraction ground gives `EX`, the content and gcd over
  `QQ[a]` follow `dup_content` and `dup_rr_prs_gcd`, an algebraic number
  whose minimal polynomial is the generator is zero, and `expand` distributes
  a rational coefficient over a power of a sum with a denominator term by
  term, as the `mul` hint of the reference does.

The ports of the integration methods also rely on these reference rules of
the polynomial layer:
- an expression is expanded deeply before it becomes a polynomial: function
  arguments are expanded and `exp` of a sum splits, so `Poly(exp(x*(x + 1)))`
  has the generators `exp(x)` and `exp(x**2)` and
  `cancel(exp(-x**2/(2*(2 - I)**2)))` is `exp(-3*x**2/50)*exp(-2*I*x**2/25)`;
  the factor named by "contains an element of the set of generators" is the
  first one in the argument order of its term
  (`tests/test_poly_deep_expansion.py`);
- `Poly` has the reference operator priority, so `expr op poly` reaches the
  reflected methods of `Poly`;
- a quotient over a legacy ring `ZZ[c]` that is not a polynomial is computed
  in `ZZ(c)` and retracted, `sqf_list_include` and `factor_list_include`
  accept `RR` and the legacy `ZZ[c]` domains, and domain strings such as
  `ZZ[pi]` or `QQ(log(x))` are parsed;
- `primitive` and `content` over a transcendental coefficient ring follow
  `dmp_ground_primitive`;
- `factor_list` and `sqf_list` of an expression with symbolic multiplicities
  keep the `int` or `Integer` type of each multiplicity and the `Mul.args`
  order of the reference (`tests/test_symbolic_multiplicity_factor_lists.py`).
Remaining differences: the order of the factors in the TypeError message
of `sqf_list` on non-integer symbolic exponents depends on the hash seed in
the reference; `ground_roots` over `RR` or `ZZ[sin(y)]` raises ValueError
(reference `{}` for `x**2 - 2.0` and `{1: 1, -1: 1}` for
`sin(y)*x**2 - sin(y)`); the ExactQuotientFailed message of `exquo` over
`ZZ[y]` differs in text.

## Frobenius series on 2026-10-21

`dsolve(..., hint='2nd_power_series_regular')` was audited against SymPy
1.14.0 on 37 equations: Bessel of orders 0, 1/3, 1/2, 1, 3/2 and 2,
Laguerre (plain, generalized and with integer parameter), Gauss
hypergeometric, Euler (distinct, double, complex and integer-differing
roots), roots differing by an integer with and without a logarithm,
complex, irrational and symbolic roots, a nonhomogeneous equation, and
the options `n=2, 3, 4, 8, 10` and `x0=1, -1`. Every case where SymPy
returns a result now prints the same string in Symox (series, constants,
grouping and `O` term); every case where SymPy raises gives the same
exception and message. Each basis element Symox returns is certified by
substituting it into the ODE (with `log` held as a symbol): the residual
vanishes through the computed order.

Fixes that made the default output match: `cse` no longer counts or
enters `Order` terms when it looks for repeated subexpressions (it used to
rewrite `O(x**6)` into `O(1)` and drop series); `terms_gcd(expand=False)`
keeps a polynomial over shifted generators such as `x - 1` unexpanded
(`clear_denoms` no longer rebuilds a ground-domain polynomial); and
multiplicative `subs` splits a product with `fraction`, as `Mul._eval_subs`
does, and keeps the residual coefficient as one factor, so
`(3*C1*(x + 5)/7).subs(3*C1, C1)` is `C1*(x + 5)/7`.

Deliberate extensions, which never change the reference part:

- When the indicial roots coincide or differ by an integer with an
  inconsistent recurrence, SymPy returns only the series of the larger
  root. Symox returns the same expression plus `C2` times the logarithmic
  Frobenius solution, which is set aside before constant simplification
  and appended afterwards: `sol.subs(C2, 0)` is exactly the reference.
- For non-polynomial `p` or `q`, SymPy calls `series(ser, x=x0, n=n+1)`,
  which leaves the expression unexpanded, and reads wrong coefficients
  (for example `-4*x` instead of `4*x` in `2F1(1, 2; 1/2; x)`, or only the
  constant term of `2F1(1/2, 1; 3/2; x)`). Symox expands in `x - x0`; its
  output is identical to SymPy's pipeline with only that call corrected.
- With irrational roots (`x**2*f'' + x*f' + (x - 2)*f`), SymPy computes the
  same raw series (`simplify=False` matches) but raises
  `NotImplementedError: Not sure of sign of 6 - x0` from `Order.contains`
  during constant simplification. Symox decides that comparison and
  returns the certified basis.

Differential tests: `tests/test_frobenius_parity.py` (hardcoded SymPy
outputs, extension and certification checks).

Remaining gaps in this area:

- Complex indicial roots (Euler `x**2*f'' + x*f' + f`, Bessel of imaginary
  order) and a symbolic order `nu` raise the reference `ValueError:
  Expected Expr or iterable but got None`; no complex or parametric
  Frobenius basis is produced.
- `x*(1 - x)*f'' + (2 - 3*x)*f' - f` raises the reference `IndexError` in
  the `2nd_hypergeometric` matcher of `classify_ode`, before the Frobenius
  hint is reached.
- Roots differing by an integer that is at least the requested number of
  terms (`x**2*f'' - 11*x*f' + 36*f`, double root 6 with `n=6`) give
  `O(1)`, as in SymPy; increase `n` to obtain the series.
- With `x0 != 0` the result keeps SymPy's `O(x**n)`, an order at 0 rather
  than at `x0`.
- The logarithmic solution is not reduced to SymPy's normal form for
  combinations of the two constants: it stays the Frobenius derivative
  itself, so `C2` may contain a multiple of the first series.

## Min/Max, polarify, Piecewise and relationals on 2026-10-21

Six small behaviors were compared against SymPy 1.14.0 and fixed:

- `Min`/`Max` keep `ordered` args, as in SymPy. The str, latex and pretty
  printers now sort the args by `default_sort_key`, as `_print_LatticeOp`
  and `_hprint_variadic_function` do, so `Min(x**2, x, sin(x))` prints
  `Min(x, x**2, sin(x))` while `.args` and `srepr` stay
  `(x, sin(x), x**2)`. Their derivative is the reference Heaviside sum
  (`Max(y, x).diff(x)` is `Heaviside(x - y)`), no longer an unevaluated
  `Derivative`. Functions with neither `nargs` nor `eval` (`exp_polar`,
  `meijerg`, `lerchphi`, `LeviCivita`, `carmichael`, the integral
  transforms) and `Min`/`Max` get the reference class key `(4, 0, name)`,
  which puts `exp_polar` first in products.
- `polarify` is a port of `_polarify`: function args are always visited,
  other nodes only visit their `Expr` args, so `meijerg` parameter tuples
  and `Piecewise` conditions keep their exact values; Integral limits and
  `ExprCondPair` follow the reference. `hyper` orders its parameters
  before unpolarifying them, as `hyper.__new__` does.
- `series` of a top-level `Piecewise` follows `Expr.series` with
  `Piecewise._eval_nseries` (branch-wise expansion, the conditions in the
  positive dummy, the `O` term added only when not absorbed), including
  `x0 != 0` and `dir='-'`. `Order.subs` now moves the limit point like
  `Order._eval_subs` (`O(x**3).subs(x, x - 2)` is
  `O((x - 2)**3, (x, 2))`, `O(x**3).subs(x, -x)` is `O(x**3)`).
- Strict and non-strict inequalities between identical sides are decided
  only when the sides are extended real or comparable numbers: `y < y`
  stays unevaluated for a plain or complex symbol, `r < r` is `False` for
  a real one, and `I < I`, `zoo < zoo`, `nan < nan` raise the reference
  `TypeError`.
- `Piecewise.as_expr_set_pairs` is a port, including the complex-domain
  `ValueError` and the multivariate `NotImplementedError` messages.
- `piecewise_fold` is a port of the SymPy 1.14 function: condition
  folding of nested `Piecewise` conditions through `ITE`/`to_nnf`/
  `simplify_logic`, shared-prefix grouping in `Add`/`Mul`, and the
  `evaluate=None` duplicate-branch merge.

Differential tests: `tests/test_misc_wave2_parity.py` (hardcoded SymPy
outputs).

Remaining gaps in this area:

- The noncommutative branches of `piecewise_fold` are not audited.
- The order of the `polarify` replacement mapping follows set iteration
  and may differ from SymPy's dict order; the mapping itself is equal.

## Periodic inequalities and Piecewise simplification on 2026-10-21

A probe of 90 cases (solveset, reduce_inequalities, `as_set`, Piecewise
integrals, simplify, rewrites, powers and series) went from 42 to 79
outputs equal to SymPy 1.14.0.

Periodic inequalities:

- `solve_univariate_inequality` has the reference periodic branch: over
  `S.Reals` an inequality of a periodic expression is solved on one
  period `[0, T)`, as SymPy does (`solveset(sin(x) > 0, x, S.Reals)` is
  `Interval.open(0, pi)`). A decision by the function range gives the
  whole domain or `EmptySet` (`reduce_inequalities(cos(x) + 2 > 0, x)` is
  `True`). The branch is in `python/symox/solvers/_periodic_inequality.py`
  and the native `solveset`, `solve_univariate_inequality` and
  `reduce_inequalities` call it through a host hook.
- `Relational.as_set` of a periodic relation raises the reference
  `NotImplementedError` unless the range decides it; then it returns
  `Interval(-oo, oo)` or `EmptySet` as SymPy does.
- `Piecewise._intervals` has the `err_on_Eq` messages. `piecewise_simplify`,
  `piecewise_simplify_arguments` and `piecewise_exclusive` are ports, and
  `simplify` routes `Piecewise`, `BooleanFunction` (through
  `simplify_univariate` and `simplify_logic`) and inequalities (the
  linear canonical form of `Relational._eval_simplify`) as SymPy does.
  An `Add` or `Mul` holding a `Piecewise` is folded and then simplified
  with `piecewise_simplify(doit=False)`.
- `Pow` of a `Piecewise` distributes over the branches
  (`Piecewise._eval_power`).
- `KroneckerDelta(i, j).rewrite(Piecewise)` is
  `Piecewise((0, Ne(i, j)), (1, True))`.
- `Piecewise.rewrite(ITE)` raises the reference `TypeError` with the
  class name of the first non-Boolean value (`Integer`, `Half`, `Add`);
  a value equal to 1 is accepted and the check stops at the first true
  condition, as in SymPy.

Deliberate divergences (SymPy's answer is wrong):

- On a bounded domain that is not inside `[0, T)` the inequality is
  solved on the domain itself. SymPy intersects the one-period answer
  with the domain: `solveset(sin(x) > 0, x, Interval(2*pi, 3*pi))` is
  `EmptySet` in SymPy and `Interval.open(2*pi, 3*pi)` here. On a
  half-unbounded domain the answer is a `ConditionSet` (or the
  `NotImplementedError` of `solve_univariate_inequality`), not a
  one-period set.
- `solveset(sin(x) >= 0, x, Interval(0, 2*pi))` includes `2*pi`, where
  `sin` vanishes; SymPy gives `Interval(0, pi)`.
- A periodic condition seen on the whole line (an antiderivative,
  `piecewise_simplify`, `simplify`) is kept unless its range decides it.
  SymPy rewrites `Piecewise((1, sin(x) > 0), (0, True)).simplify()` with
  `0 < x < pi`, which is false at `x = 2*pi + 1`.
- A definite integral with numeric limits solves a periodic condition on
  the integration interval:
  `integrate(Piecewise((1, sin(x) > 0), (0, True)), (x, 0, 4*pi))` is
  `2*pi`; SymPy gives `pi`.
- `simplify_univariate` tests each interval endpoint on the whole
  condition. SymPy rebinds the condition in its loop and can drop an
  endpoint; here `Or(x <= 0, x >= 1)` keeps `x = 0`.

Differential tests: `tests/test_piecewise_periodic_parity.py` (hardcoded
SymPy outputs; the divergence tests assert the correct values).

Remaining gaps in this area:

- `solvify` cannot solve `sin(x) + cos(x)` or `exp(sin(x))` on a period,
  so `solveset(sin(x) + cos(x) > 0, x, Interval(0, 2*pi))` stays a
  `ConditionSet`; SymPy gives the intervals.
- `reduce_inequalities([sin(x) > 0, x > 0, x < 2*pi], x)` drops the
  redundant `x < 2*pi` that SymPy keeps (form only).
- The series, `And`/`Or` and `Order.subs` gaps listed here before
  2026-10-22 are closed, see the next section.

## Complex expansion without limits and Heaviside(x, H0) on 2026-10-22

`tests/test_core_misc_wave3_parity.py` holds the reference outputs.

- `expand(complex=True)`, `expand_complex`, `as_real_imag`, `re` and `im`
  have no default size limit, as in SymPy. The guards of 4096 terms,
  integer powers through 128, 20 000 nodes and 1 000 000 numeric bits used
  to raise `NotImplementedError` where SymPy computes the result (for
  example `expand((a + I*b)**200, complex=True)`, `re((a + I*b)**140)`,
  `((x + I)**130).as_real_imag()` or `re` of a 2000-term sum). These guards
  remain only on `as_real_imag_bounded`, the opt-in entry point of the
  internal heuristics without a SymPy counterpart (exact assumption
  refinement, `solveset`'s nonzero proof, the complex split of
  `floor`/`ceiling` and the exact numeric parts of the domains); they raise
  and never return a partial result, and their callers already treat the
  error as "unknown". On a probe of 20 calls Symox now matches SymPy on 19
  (8 before) and is faster, e.g. `expand((a + b + I + x)**25,
  complex=True)` in 0.4 s against 8.1 s.
- Three guards remain on the public entry points, all raising
  `NotImplementedError`: an expression nesting deeper than 256 levels (it
  protects the native stack), the binomial expansion of `(a + I*b)**n` with
  a non-real base for `|n| > 2**17`, whose output alone exceeds 2 GiB, and
  an unevaluated numeric product or sum that multiplies out to more than
  `2**34` bits (a shared `Mul(v, v)` DAG nested 80 times). SymPy runs out of
  memory or never finishes on these; the native code could not be
  interrupted. A real base keeps the power whole at any
  exponent, as in SymPy (`Pow(2, 10**9, evaluate=False).as_real_imag()` is
  `(2**1000000000, 0)`).
- `Integer**Rational` whose integer part does not fit in 64 bits stays a
  held power; it used to drop the integer part (`Pow(2, Rational(10**20 + 1,
  2))` gave `sqrt(2)`, a wrong value). SymPy never finishes this call.
- `Heaviside(x, H0)` is a native two-argument node: `args` is `(x, H0)` and
  `(x, 1/2)` for the default, `pargs` drops a default `1/2`, `Heaviside(x,
  S.Half)` is `Heaviside(x)`, a Float `H0` is kept (`Heaviside(x, 0.5)`).
  Printing (str, srepr, latex, pretty, C code), `rewrite(Piecewise)`, the
  value `H0` at zero (also through `subs` and `xreplace`),
  `diff`, the
  ordering of `Heaviside(x)` against `Heaviside(x, 0)`, `is_real`,
  `is_finite`, `conjugate`, `as_real_imag`, `lambdify`, pickling and the
  conversions to and from SymPy follow the reference. `integrate` and
  `laplace_transform` ignore `H0` like the SymPy rules do
  (`integrate(Heaviside(x, 0), x) == x*Heaviside(x)`); an integral left
  unevaluated keeps its `H0`. `Heaviside(x, None)` raises `AttributeError`
  and a nonreal argument raises `ValueError`, as in SymPy. A probe of 96
  calls went from 39 to 89 outputs equal to SymPy 1.14 (see below for the
  rest; one more is a set printed in another order).
- A number times a noncommutative sum stays a product, as in
  `Mul.flatten` (`-(A + B)` is `Mul(-1, A + B)`, `-(A + B) + (A + B)` prints
  `A + B - (A + B)`), and `expand` distributes it. A probe of 21
  noncommutative calls went from 8 to 21 outputs equal to SymPy.

Remaining gaps in this area:

- `laplace_transform(Heaviside(t - 2)*t, t, s)` keeps the factored form
  `(2/s + s**(-2))*exp(-2*s)` where SymPy distributes it (equal values).
- The error message of `Heaviside(I, 0)` is shorter than the reference.
- `lambdify` of `Heaviside(x, 0)` returns the float `0.0` at zero where
  SymPy's generated code returns the integer `0` (equal values).

## Limits and series at jumps on 2026-10-22

`tests/test_heaviside_limit_parity.py` holds the reference outputs and the
corrections below.

- Wrong results fixed: `limit(Heaviside(x), x, 0)` returned `1/2` for every
  direction (SymPy: `1`, `0`, `ValueError` for `'+-'`), so did every
  expression whose limit went through `Heaviside` at its jump
  (`Heaviside(x)*exp(x)`, `Heaviside(x)**2` gave `1/4`, `Heaviside(1/x)` at
  `oo`, `Heaviside(x)*Heaviside(-x)` with `'+-'` gave `1/4`).
  `limit(Heaviside(x)/x, x, 0, '-')` gave `-oo` and `'+-'` gave `zoo` (both
  `0`/`ValueError`), `limit(Heaviside(x)*log(-x), x, 0, '-')` gave `-oo`
  (`0`), and `limit(DiracDelta(x), x, 0)` gave `DiracDelta(0)` (`0`).
- The series engine expands `Heaviside(u, H0)` and `DiracDelta(u, k)`
  one-sided: the leading coefficient of `u` on the side approached gives
  `1` or `0` (`H0` is never reached) and `DiracDelta` vanishes, which is
  what SymPy's positive dummy produces. A nonvanishing symbolic value is
  continuous (`limit(Heaviside(x - y), x, 0)` is `Heaviside(-y)`); an
  undecided sign (`Heaviside(a*x)`) stays a `Limit`. The gruntz rewrite does
  the same at infinity (`limit(Heaviside(x), x, oo)` is `1`).
- `floor`, `ceiling` and `frac` expand to the one-sided constant (at an
  integer, the side is the sign of the first non-constant term), so their
  limits at jumps, `series(floor(x), x)` (`0`), `series(frac(x), x)` (`x`)
  and `floor(x).as_leading_term(x)` (`0`) follow SymPy; they used to stay
  unevaluated or raise. An expression that is exactly zero on the side
  approached has limit and leading term `0`.
- `limit(..., '+-')` of a `Piecewise` combines the one-sided limits and
  raises `ValueError` when they differ.
- `Heaviside(u, H0).evalf(subs=...)` (one or two arguments) and `N` give the
  value: `0` at negative points (the integer, as SymPy), `H0` at zero.
- The derivative with respect to a symbol inside `H0` follows
  `Function.fdiff`: `Derivative(Heaviside(x, y), y)` for a symbol `H0` not
  in the first argument, `Subs(Derivative(Heaviside(u, _xi_2), _xi_2),
  _xi_2, H0)` times the chain factor otherwise.
- A probe of 254 calls went from 70 to 194 outputs equal to SymPy 1.14. The
  other 60 are the corrections below, `Abs(x)/x`-style series that print a
  trailing order term, and evalf misses.

Intentional differences (SymPy 1.14 is wrong; the correct value is kept):

- `limit(Heaviside(sin(x)), x, 0)`, `Heaviside(x**3 - x)` and
  `Heaviside(exp(x) - 1)` give `1/2` on both sides in SymPy (gruntz
  evaluates the argument at the point); Symox gives the one-sided values
  and raises `ValueError` for `'+-'`. `series(Heaviside(sin(x)), x)` is `1`
  (SymPy: `1/2 + x*DiracDelta(0) + ...`).
- `limit(Heaviside(a*x), x, 0)` is `1/2` in SymPy although the sign of `a`
  is unknown; Symox keeps the `Limit`.
- `limit(Heaviside(x), x, -1)` and `limit(Heaviside(x - 1), x, 0)` never
  finish in SymPy; Symox gives `0`.
- `limit(Piecewise((0, x < 0), (1, True)), x, 0, '-')` is `1` in SymPy (the
  branch at the point); Symox gives `0`, and `-1` for
  `Piecewise((x - 1, x < 0), (x + 1, x >= 0))`; `'+-'` raises `ValueError`
  where SymPy gives `1`.
- `nseries` of `Heaviside` or `DiracDelta` at the jump is a Taylor
  expansion through the discontinuity in SymPy (`1/2 + x*DiracDelta(0) +
  x**2*DiracDelta(0, 1)/2 + O(x**3)`); Symox returns the one-sided
  expansion, as `series` does in both.
- `as_leading_term` of `Heaviside`, `DiracDelta`, `sign` and `Abs` raises
  `NotImplementedError` in SymPy (no `_eval_as_leading_term`); Symox
  returns the one-sided leading term (`Heaviside(-x)` gives `0`).
- With `limit(DiracDelta(x), x, oo)` now `0`, a delta derivative supported
  inside the interval or at infinity integrates to its value:
  `integrate(DiracDelta(x - 1, 1), (x, 0, oo))` and the same with a
  positive support `p` give `0` (SymPy: `DiracDelta(oo)`), and
  `DiracDelta(x, 1)` over `(-oo, oo)` gives `0` as in SymPy. When a delta
  of the antiderivative may sit on a finite bound, the definite integral
  stays unevaluated instead of the endpoint difference `DiracDelta(0)`
  that SymPy returns (`DiracDelta(x, 1)` over `(0, oo)`, `(-oo, 0)` and
  `(0, 1)`, `x*DiracDelta(x, 2)`, `DiracDelta(x, 1)*Heaviside(x)`,
  `DiracDelta(x - a, 1)` over `(0, oo)`); `(0, 1)` used to give
  `-DiracDelta(0)`.

Remaining gaps in this area:

- `series(Heaviside(x)/x, x)` and `series(Abs(x)/x, x)` print `1/x +
  O(x**3)` and `1 + O(x**3)` where SymPy returns the exact `1/x` and `1`
  (equal values).
- `DiracDelta(x).evalf(subs={x: 0})` gives `DiracDelta(0)` and
  `Heaviside(x, y).evalf(subs={x: 0})` gives `Heaviside(0, y)`; SymPy
  returns the original expression.

## Piecewise series, And/Or patterns and Order.subs on 2026-10-22

A probe of 112 cases (series of sums and products holding a `Piecewise`,
the order term per branch, powers, `piecewise_fold`, limits, `And`/`Or`
simplification, `Heaviside(x, H0)`, periodic inequalities and
`Order.subs`) went from 48 to 97 outputs equal to SymPy 1.14.0.

- `series` of any expression holding a `Piecewise` in the expansion
  variable follows `Expr.series` (`python/symox/series/_piecewise_series.py`):
  the positive dummy (which often decides the conditions), the shift for
  `x0 != 0` and `dir='-'`, then `Add._eval_nseries` term by term,
  `Mul._eval_nseries` with the leading terms of the factors (a `Piecewise`
  contributes the leading term of its first branch true at 0, as
  `Piecewise._eval_as_leading_term`), and the order term handling of
  `Expr.series`. `nseries` at 0 from the right is `_eval_nseries` without
  the dummy (`(Piecewise((sin(x), x < 0), (cos(x), True))*x).nseries(x, 0, 4)`
  is `x + O(x**4)`).
- A branch is expanded with `_eval_nseries` semantics: a power `x**k` with
  `k > n` is `O(x**k)` (`Pow._eval_nseries`, `maxpow < 0`), so
  `Piecewise((x, x < 1), (x**5, True)).series(x, 0, 3)` is
  `Piecewise((x, x < 1), (O(x**5), True))`.
- `simplify` of `And`/`Or` ports `And._eval_simplify` (linear `Eq`
  substitution, `And(Eq(x, 1), y > x)` is `Eq(x, 1) & (y > 1)`),
  `Or._eval_simplify` (`a <= b` split into `Eq(a, b) | (a < b)`) and
  `_apply_patternbased_simplification` with the two- and three-term
  pattern tables (`python/symox/logic/_pattern_simplify.py`):
  `Or(x < 0, x < 1)` is `x < 1`, `Or(x > 0, x < 1)` is `True`,
  `Or(x > 1, x < -1)` is `Abs(x) > 1`.
- `Order.subs` with a non-affine replacement solves the moving point as
  `Order._eval_subs` does (`solveset` of `old = new(d)`, first solution,
  `limit` at the old point) through a host hook
  (`python/symox/series/_order_point.py`): `O(x**2).subs(x, 1/x)` is
  `O(x**(-2), (x, oo))`, `O(x**3).subs(x, exp(x))` is
  `O(exp(3*x), (x, -oo))`.
- `limit` of an expression holding a `Piecewise` at a finite point from
  one side substitutes `z0 + d` (or `z0 - d`) with a positive `d`; when
  every condition is decided, the limit of the selected branches is
  taken.

Deliberate divergence (SymPy's answer is wrong):

- `limit(Piecewise((sin(x), x < 0), (cos(x), True)) + x, x, 0, '-')` is
  `0` here, the left branch value; SymPy gives `1`, the value of the right
  branch.

Differential tests: `tests/test_piecewise_series_parity.py` (hardcoded
SymPy outputs).

Remaining gaps in this area:

- `solvify` cannot solve `sin(x) + cos(x)`, `sin(x) - cos(x)`,
  `sin(2*x) + cos(x)` or `exp(sin(x)) - 1` on a period, so the periodic
  inequalities and the Piecewise integrals with such conditions stay a
  `ConditionSet` or an `Integral`; SymPy gives the intervals and `pi`.
- `O(x*y).subs(y, x)` gives `O(x**2)`; SymPy raises `ValueError`.
- `Xor` simplification by patterns (`Xor._eval_simplify`) is not wired;
  its pattern table is ported.
- A function of an undecided `Piecewise` (`exp(Piecewise(...))`) still
  raises `NotImplementedError` in `series`; `Mul._eval_nseries` is ported
  only for products holding a `Piecewise`, with Order terms of the factors
  skipped in the term products.

## Set emptiness and the subset decision path on 2026-10-21

A probe of 50 cases (`is_empty`, `is_finite_set`, `is_subset`,
`is_superset`, `is_proper_subset` on unions, intervals, ranges and named
sets, and monotonicity on disconnected domains) went from 29 to 49 outputs
equal to SymPy 1.14.0. `tests/test_union_sets_parity.py` holds the
reference values.

- `Set.is_subset` is a port of the reference decision path in
  `python/symox/_set_subset.py`: the trivial cases (equality, emptiness,
  finiteness), the `_eval_is_subset`/`_eval_is_superset` hooks
  (`FiniteSet`, and `Naturals`, `Naturals0`, `Integers` through
  `Range(.., oo)`), the pairwise `is_subset_sets` rules (`Reals` counts as an
  `Interval`, as in SymPy), then the `self.intersect(other) == self`
  fallback. `Union.is_subset` is the fuzzy and of its components.
  `is_superset`, `is_proper_subset` and `is_proper_superset` delegate as in
  the reference and raise `ValueError` for a non-set argument.
- `Union.is_empty` is the fuzzy and of its components and
  `ProductSet.is_empty` the fuzzy or; `Interval.is_empty` is
  `fuzzy_bool(start > end)` (`>=` when a bound is open), so a symbolic
  `Interval(a, b)` is undecided instead of non-empty. `Intersection` and
  `Complement` have no emptiness rule, as in the reference. `is_finite_set`
  follows the reference for unions, intersections, complements, products,
  intervals (`measure.is_zero`) and ranges.
- `S.Reals.is_subset(Union(Interval(0, pi/2), Interval.Ropen(3*pi/2, 2*pi)))`
  is now False, so `is_increasing(sin(x), symbol=x)` is False as in SymPy.
- The interval intersection rule with `(-oo, oo)` returns the other
  interval even when its bounds are not comparable
  (`Interval(0, z).intersect(S.Reals)` is `Interval(0, z)`).
- Monotonicity: when the native proof declines (disconnected or non-real
  domain, parameters not known to be real, bounds depending on the
  variable) or is undecided, the reference `monotonicity_helper`
  (`interval.is_subset(solveset(predicate(derivative), x, S.Reals))`) and
  the reference `is_monotonic` (no turning point in the interval) give the
  answer. `is_increasing(x, FiniteSet(0, 1), x)`, on a union of intervals
  and with `symbol=x` on `x + y` are True; on `S.Complexes` it is False;
  `is_monotonic(a*x, symbol=x)` is True.

Remaining gaps:

- `Interval(a, b).is_subset(Interval(a, b + 1))` is True; SymPy gives None
  because its intersection keeps `Min(b, b + 1)` unevaluated while the
  native intersection simplifies it (the True answer is correct).
- `is_increasing(a*x, symbol=x)` with a real parameter `a` stays None;
  SymPy gives False (its `solveset(a >= 0, x, S.Reals)` decides it).
- `is_increasing(Function('f')(x), symbol=x)` raises `NotImplementedError`;
  SymPy overflows the recursion limit.
- The native intersection and complement simplification still uses the
  approximate Rust `is_subset` (partly addressed on 2026-10-22, see the
  `imageset` section below).

## Interval(-oo, oo), Reals and the union rules on 2026-10-22

A probe of 166 calls (set operations, `solveset`,
`solve_univariate_inequality`, `reduce_inequalities`, `imageset`,
`ConditionSet`, `Piecewise.as_expr_set_pairs`, `integrate` of a
`Piecewise`, `is_increasing`, `singularities`, `continuous_domain`,
`function_range`, symbolic `Range` in `FiniteSet`) went from 124 to 154
outputs equal (by `srepr`) to SymPy 1.14.0.
`tests/test_sets_reals_parity.py` holds the 152 agreed reference values.

- SymPy 1.14 does not turn `Interval(-oo, oo)` into `Reals` in
  `Interval.__new__`: `srepr(Interval(-oo, oo))` stays an `Interval`. Instead
  `Reals.__eq__` and `Reals.__hash__` defer to `Interval(-oo, oo)`. The native
  node of `Reals` now has the hash of `Interval(-oo, oo)` and both compare
  equal, also inside other expressions (`ProductSet(Interval(-oo, oo),
  Interval(-oo, oo)) == ProductSet(S.Reals, S.Reals)`). As in the reference,
  `Union`, `Intersection` and `FiniteSet` keep the first of the two forms,
  `Interval(-oo, oo).is_subset(S.Reals)` is True, and `is_increasing`,
  `continuous_domain`, `function_range` and `as_expr_set_pairs` accept
  `Interval(-oo, oo)` as the real line.
- `solve_univariate_inequality(x**2 + 2 > 0, x, relational=False)` gives
  `Reals`, and a symbolic `Range` stays an element of a `FiniteSet`, as in
  the reference; `Range.size` is the reference property (symbolic when
  decidable, `ValueError` otherwise).
- The native union simplification follows the reference `union_sets`
  rules instead of a generic subset absorption: intervals merge only when
  all their endpoints are comparable (`Union(Interval(-oo, a),
  Interval(a, oo))` stays a union), open finite endpoints contained in the
  other set are closed (`Union(S.Integers, Interval.open(0, 1))`), elements
  of a `FiniteSet` contained in the other set are dropped, `Integers`
  absorbs or is absorbed through its intersection, `Complexes` absorbs real
  sets, `ProductSet` uses `is_subset` and the common factor rule, and the
  named sets follow the reference table (`Union(S.Reals, Range(0, 5))` stays
  a union, as in SymPy).

Remaining gaps: see the next section.

## imageset, image-set lattices and products on 2026-10-22

The probe was extended to 285 calls (119 new ones on `imageset`,
intersections, complements and `ProductSet`). It went from 214 to 281
outputs equal (by `srepr`) to SymPy 1.14.0: the 166 previous calls from 154
to 166, the new ones from 60 to 115. `tests/test_sets_rest_parity.py` holds
the 127 newly agreed reference values.

- `imageset` is a port of the reference `imageset` and of the
  `_set_function` handlers (`python/symox/_set_extensions.py`): a `Lambda`
  over an `Interval` or `Reals` is evaluated from the singularities in the
  interval, the one-sided limits at open ends and the critical points of
  `solveset(diff(expr, x), x)` (`Interval(1, 3)` for `2*x + 1` on
  `[0, 1]`, `Interval(0, oo)` for `x**2` on the reals, a union of two rays
  for `1/x`), and stays unevaluated where the reference does (several free
  symbols, a symbolic bound, a derivative whose zeros are not a finite set,
  such as `sin` or `Abs`). `Piecewise` maps are imaged piece by piece. The
  `FiniteSet`, `Union`, `Intersection` (when the map is invertible),
  `Range`, `Integers` (canonical `a*n + b mod a`), `Naturals` and generic
  handlers follow the reference, including the composition of nested image
  sets and the `Symbol('x')` variable of the generic handler. As in the
  reference, several base sets or a `Range` give an unevaluated `ImageSet`
  (the previous enumeration of finite products is gone), a constant map on a
  non-empty base is a `FiniteSet`, and a signature that does not match the
  base sets raises `ValueError`.
- The native intersection has the reference `ProductSet` rule (empty for
  different lengths, else the componentwise intersection), the
  `ConditionSet` rule (the condition set over the intersected base) and the
  linear case of the `ImageSet` rule: `{p*n + q}` against `Integers` or
  another such lattice with exact rational coefficients is solved as a
  linear diophantine equation and returned in the canonical form the
  reference `imageset` builds (`Intersection(S.Integers, ImageSet(Lambda(n,
  2*n), S.Integers))` is the image set, `2*n` with `3*n` gives `6*n`, `n/2`
  gives `Integers`). Through the `Integers` union rule
  `Union(S.Integers, ImageSet(Lambda(n, 2*n), S.Integers))` is `Integers`.
- The native `is_subset` decides `ProductSet` against `ProductSet` and an
  integer lattice against `Integers` by the reference fallback
  `self.intersect(other) == self` (True or None), so a product contained in
  another is absorbed by a union, and the complement has the reference rule
  `(A x B) - (C x D) = ((A - C) x B) U (A x (B - D))`.
- Unevaluated `Union` and `Intersection` order their arguments as the
  reference `ordered(args, Set._infimum_key)`: by the numeric infimum (`oo`
  when it is not comparable), ties in the default order.

Intentional correction:

- `imageset(Lambda(x, tan(x)), Interval(0, 1))` is `Interval(0, tan(1))`.
  SymPy 1.14 keeps the image set because its `solveset` leaves
  `tan(x)**2 + 1 = 0` as a `ConditionSet`; the equation has no solution (the
  native `solveset` gives `EmptySet`), so `tan` is monotonic on the interval
  and the interval is its exact image.

Remaining gaps: see the next section (the `Reals` and `Interval` cases of
the `ImageSet` rule and the `Range` rules were added there).

## Range, image-set and complex intersections on 2026-10-22

A probe of 218 calls (intersections, complements, symmetric differences
and membership with `Range`, `Interval`, `ImageSet`, `FiniteSet` with
symbols, `ProductSet`, `ConditionSet`, `Naturals`, `Integers`, `Reals`,
`Complexes`) went from 147 to 211 outputs equal (by `srepr`) to SymPy
1.14.0. `tests/test_sets_intersection_parity.py` holds the 64 newly agreed
reference values.

- The native intersection has the reference `Range` rules: `Range` with
  `Range` (the common arithmetic progression by the extended gcd, as a
  `Range` with a positive step, including half-infinite and decreasing
  ranges), `Range` with an `Interval` or `Reals` whose bounds are numbers
  (the interval trimmed to its integers: `Intersection(Range(0, 10),
  Interval(2, 5))` is `Range(2, 6)`, irrational and float bounds included),
  `Range` with `Naturals` and `Naturals0` (through `Interval(inf, oo)`), and
  `Range` with `Integers` or `Rationals` (the range). A symbolic range or
  interval stays unevaluated, as in the reference.
- `Range.contains` of an integer expression in a numeric range of unit step
  is the reference relation (`And(n >= 0, n <= 9)` for `Range(0, 10)`,
  `n >= 0` for `Range(0, oo)`, `Eq` for one element).
- The Python `Intersection` applies the reference `intersection_sets
  (ImageSet, Set)` handler (`python/symox/_set_intersection.py`) to what the
  native rules leave unevaluated: the diophantine equation `f(n) = g(m)`
  against `Integers` or another integer image set (`n**2` with `m**2 + 1`
  gives `{1}`), the real part of the map against `Reals` (the image set
  itself for `2*n` and `n*pi`, `EmptySet` for `n + I`, `{0}` for `x + I*x`
  over the reals, with the denominators excluded), and the inverted bounds
  of an `Interval`. `Complexes` against a set the reference subset decision
  proves real gives that set's intersection with `Reals`
  (`Intersection(S.Complexes, ImageSet(Lambda(n, 2*n), S.Integers))` is the
  image set). `ImageSet.base_sets` and `base_pset` exist as in the reference.
- `Complement` applies the first step of the reference `Complement.reduce`
  to what the native rules leave unevaluated: `EmptySet` when the reference
  subset decision (`python/symox/_set_subset.py`) proves `A <= B`
  (`Complement(Range(0, 5), Interval(-1, 10))`).

Remaining gaps (updated by the follow-up below):

- **Intentional correction**: when the inverted bounds are not real, the
  reference intersects the real zeros of the map with the interval, so
  `Intersection(Interval(0, 1), ImageSet(Lambda(n, n*I + 1), S.Integers))`
  is `EmptySet` in SymPy 1.14 although `n = 0` gives `1`. Symox intersects
  the real part of the image set with the interval and gives `{1}`.

## Inversion, diophantine and membership follow-up on 2026-10-22

A probe of 101 calls (`invert_complex`/`invert_real` on symbols with
assumptions, binary diophantine equations, image-set intersections,
`Range`, `ImageSet` and `ProductSet` membership) went from 30 to 96 outputs
equal (by `srepr`) to SymPy 1.14.0, and a probe of 31 native `solveset` and
intersection calls on assumed symbols went from 26 to 29.
`tests/test_sets_contains_parity.py` holds the agreed reference values.

- `solveset` follows the reference handling of the symbol's assumptions: a
  symbol with assumptions (other than a plain `real=True` symbol over a real
  domain) is replaced by a `Dummy` that is `real` when the domain is a
  subset of `Reals` and `complex` otherwise, the equation is solved, and the
  symbol is put back. The assumptions no longer filter the solutions, so
  `invert_complex(n*I, 1, n)` with an integer `n` is `(n, {-I})` and
  `solveset(2*k - 3, k, S.Reals)` with an integer `k` is `{3/2}`, as in
  SymPy.
- `diophantine` applies the reference preprocessing (numerator,
  `factor_terms`, then the factor independent of the variables is dropped,
  so `2*pi*n - pi*m` is `2*n - m`), factors a repeated factor, solves the
  linear equation with the reference `Linear.solve` (`base_solution_linear`
  and the extended gcd, variables in the default sort order, the reference
  parameter signs), and has the reference parabolic binary quadratic case
  (`discriminant == 0`: `n**2 - m`, `n**2 + n - 2*m`, `(n - m)**2 - n - m`,
  `4*n**2 + 4*n*m + m**2 + n + 1`). `Intersection(S.Integers,
  ImageSet(Lambda(n, n**2), S.Integers))`, `n**2` with `2*n`, and
  `2*n*pi` with `n*pi` are the reference image sets.
- `Range.contains` with any step is the reference rule: the relation
  `And(x >= inf, x <= sup)` when `(ref - x) % step` is zero (for an integer
  expression whose non-constant terms are multiples of the step), `false`
  when it is a nonzero integer, `Eq` for one element, unevaluated otherwise.
- `ImageSet.contains` returns only a decided value; an undecided membership
  stays `Contains(x, ImageSet(...))` as in the reference.
- `ProductSet.contains` of a tuple of the right length is the conjunction of
  the component memberships, `false` for a tuple of another length or a
  non-tuple number, and unevaluated for a symbol.
- **Intentional corrections** (SymPy 1.14 loses solutions through the
  assumption filtering of its image-set intersection):
  `solveset(exp(I*x) - 1, x, S.Reals)` is `ImageSet(Lambda(n, 2*n*pi),
  S.Integers)` (SymPy `{0}`), and `solveset(exp(I*x) - I, x, Interval(0,
  10))` is `{pi/2, 5*pi/2}` (SymPy `{pi/2}`).

Remaining gaps:

- `invert_real`/`invert_complex` are `(x, solveset(f - y, x, domain))`
  rather than a port of the reference `_invert`, so the partial inversions
  the reference returns differ: `invert_complex(n**2, 2, n)` gives
  `(n, {-sqrt(2), sqrt(2)})` where SymPy gives `(n**2, {2})`; the same holds
  for `invert_real(r**2, -1, r)` and `invert_complex(k**2, 4, k)`.
- `invert_complex(exp(r), -1, r)` gives the image set
  `2*n*I*pi + I*pi` where SymPy gives `I*(2*n*pi + pi)` (same set,
  different form).
- The image-set and complex intersection rules run in the Python
  constructors only. The probed native `solveset` paths already agree with
  SymPy, but intersections built inside native code do not use those rules.
- The native pairwise intersection still absorbs one argument into the other
  when the approximate Rust `is_subset` proves inclusion, and the native
  complement still starts with that `is_subset`, where the reference only has
  specific handlers. Native `solveset` relies on this absorption for image
  sets against `Reals`, so removing it needs the image-set handlers in native
  code first.

## Solver targets, solve_linear and matrix structure on 2026-10-22

A probe of 182 calls (`solve`, `solveset`, `linsolve`, `nonlinsolve` with
expression, number, `Indexed`, function and derivative targets,
`solve_linear`, `Matrix.diag`, `extract`, the structural matrix predicates
and `.kind` of matrices and matrix expressions) went from 112 to 178 outputs
equal (by `str` and exception type) to SymPy 1.14.0.
`tests/test_wrong_linalg_wave4_parity.py` holds the reference values.

- `solve` follows the reference `recast_to_symbols`: an unknown that is not
  a `Symbol` (an expression such as `x + 2`, a number, an `Indexed`) is
  replaced by a numbered `Dummy` in `ordered` order, the system is solved,
  and the keys and values are put back, the dict keys in `ordered` order.
  A number target also replaces its negative, as `Number._eval_subs` does,
  so `solve(x - 2, 2)` is `[x]` and `solve(sqrt(2) - 1, 1)` is `[sqrt(2)]`.
  `solveset` and `nonlinsolve` recast an expression target the same way, so
  `solveset(x + 2 - y, x + 2)` is `{y}` and `nonlinsolve([f(x)**2 - 1],
  [f(x)])` is `{(-1,), (1,)}`.
- `solve_linear` is a port of the reference: a symbol is only isolated when
  the derivative of the numerator survives differentiation by its free
  symbols, so `solve_linear(y*cos(x)**2 + y*sin(x)**2 - y)` is `(0, 1)`;
  derivatives of the symbols are evaluated, a candidate that zeroes a
  denominator is rejected, numeric integrals in the solution are evaluated,
  and a non-expression argument raises `AttributeError`.
- `Matrix.diag` has the reference `strict`, `unpack`, `rows`, `cols` and
  `cls` options: a single list is unpacked, a non-strict flat list is a
  column and ragged rows are padded with zeros, `rows`/`cols` pad the
  result and must be integers, and a too-small size raises `ValueError`.
  The module-level `diag` is `strict=True, unpack=False`, and
  `ImmutableMatrix.diag` and `SparseMatrix.diag` share the same code.
- `extract` accepts boolean masks and negative indices and raises
  `IndexError` out of range, for dense, immutable and sparse matrices.
- `is_symmetric` simplifies the differences with `simplify` by default (or
  the given function, or none with `simplify=False`), so
  `Matrix([[1, (x + 1)**2], [x**2 + 2*x + 1, 2]]).is_symmetric()` is True.
  `is_anti_symmetric` and `is_diagonal` are the reference fuzzy
  conjunctions, so an off-diagonal entry that is zero but not structurally
  zero gives `None`.
- `.kind` of every matrix expression (`MatrixSymbol`, `MatMul`, `MatAdd`,
  `Identity`, `ZeroMatrix`, `Transpose`, `Inverse`, powers) and of sparse
  matrices is `MatrixKind(NumberKind)`. The internal dispatch tag that
  used the name `kind` is now `_kind_tag`.

Intentional divergence:

- Function and derivative targets (`f(x)`, `f(x).diff(x)`) stay on the
  native `solve` path instead of being swapped for dummies. The probed
  outputs already agree with SymPy, and the native path is the one the ODE
  solvers rely on.
- `diag(MatrixSymbol('A', 2, 2))` places the explicit entries of the matrix
  expression, as symox did before; SymPy raises `AttributeError` because a
  matrix expression has no `tolist`.

Remaining gaps:

- The general `subs` does not let a number replace its negative
  (`(x - 2).subs(2, y)` is `x - 2` where SymPy gives `x - y`); only the
  solver recast applies that rule.
- `Indexed` is a Python object outside the native expression tree, so
  `atoms(Indexed)` of a sum is empty and `solve(eqs, eqs.atoms(Indexed))`
  gives `[]`. An explicit list of `Indexed` targets works.
- `solve(Integral(x, x) - y, y)` gives `[x**2/2]` where SymPy keeps
  `[Integral(x, x)]`: the default simplification evaluates the integral
  (`simplify=False` gives the reference value).
- `type(...)` of an immutable dense matrix is the native `Matrix` class, so
  `type(ImmutableMatrix.diag(1, 2)).__name__` is `Matrix` instead of
  `ImmutableDenseMatrix`; `isinstance(..., ImmutableMatrix)` holds.

## Special numbers, Indexed sums and ConditionSet.subs on 2026-10-22

There are 301 probes in three groups:

- 96 grid probes for genocchi, bernoulli, euler, andre and harmonic: 16
  differed from SymPy before, 1 after.
- 129 probes for rf/ff on `Poly` and `Ei` with polar arguments: 58 before,
  1 after.
- 76 probes for Sum/Product with Indexed terms, limit_seq, piecewise double
  integrals and ConditionSet.subs: 43 before, 19 after (four are the
  intentional ConditionSet divergence below).

The tests are in `tests/test_wrong_calculus_wave4_parity.py`.

- `genocchi(n, x)` used the wrong sign on every constant term except
  `G_1`, since `G_k(0) = -G_k`. This made every polynomial and every value
  at an integer or rational point wrong for `n >= 2`. Examples:
  `genocchi(4, x) = -4*x**3 + 6*x**2 - 1` and `genocchi(29, 2) = -58`.
- `euler(-1) = pi/2` and
  `euler(-1, x) = digamma((x + 1)/2) - digamma(x/2)`.
  `euler(0[, x]) = 1`, and an odd positive symbolic `n` gives 0.
- In `rf`/`ff` at minus infinity, `ff(-oo, k)` for an integer `k > 0` now
  follows the parity of `k` like `rf`: `ff(-oo, 2) = oo` and
  `ff(-oo, 3) = -oo`. For `k < 0`, both functions give `oo` at `±oo`. The
  rule "a negative integer `x` with a non-integer `k` gives 0" now applies
  only to `rf`.
- `rf(Poly, k)` and `ff(Poly, k)` with an integer `k` multiply shifted
  copies of the polynomial and keep the `Poly` type and domain. A negative
  `k` gives the reciprocal expression. A multivariate `Poly` raises
  `ValueError` ("only defined for polynomials on one generator").
- `Ei(z).evalf()` adds `I*pi` when `z` is polar and `z/exp_polar(I*pi)` is
  positive, like `Ei._eval_evalf`. For example,
  `Ei(polar_lift(-1)).n() = -0.2193... + 3.1415...*I`, and
  `exp_polar(-I*pi)` gives `-I*pi`. Precision follows the requested
  number of digits.
- `Sum`, `Product`, `summation`, `product` and `doit` of terms with symbolic
  `Indexed` indices (`f[n]`, `f[n, m]`, `f[2*n]`, `f[n + 1] - f[n]`) now
  expand the indices. Before, the flattened symbol `f[n]` was treated as a
  constant, so `Sum(f[n], (n, 0, 3)).doit()` returned `4*f[n]`. Reversed
  limits and symbolic upper bounds follow the reference.
- `ConditionSet.subs` of the bound symbol leaves the symbol and the
  condition unchanged and substitutes only in the base set, like
  `ConditionSet._eval_subs`. Examples:
  `ConditionSet(x, x < 1, {x, z}).subs(x, 0) = Union(ConditionSet(x, x < 1, {z}), {0})`
  and `ConditionSet(x, x < 1, Reals).subs(x, y)` unchanged.
- The gruntz `mrv` set of an unevaluated `Sum`, `Product`, `Integral`,
  `Limit`, `Subs` or `Order` is unknown, as in the reference, which raises.
  Before, the binder counted as a constant. Then
  `limit(Sum(2**m/m, (m, 1, n))/(2**n/n), n, oo)` returned 0, and so did
  `limit_seq` of that expression and of the nested-sum example
  `Sum(k**2*Sum(2**m/m, (m, 1, k)), (k, 1, n))/(2**n*n)`. The limit is now
  left unevaluated. limit_seq takes the Stolz-Cesaro difference path, which
  gives 2 for the first expression. The nested example gives 4 when it is
  computed first in a session.

Divergence kept, to be validated by the user: SymPy's
`ConditionSet._eval_subs` drops a substitution of a free symbol in the
condition whenever the replacement is not a symbol. For example,
`ConditionSet(x, x < y, Reals).subs(y, 1)` stays `ConditionSet(x, x < y, Reals)`,
and `ConditionSet(x, x > y, {1, 2, z}).subs(y, 1)` keeps `y`. Both are
wrong. Symox applies the substitution, which gives
`ConditionSet(x, x < 1, Reals)` and
`Union(ConditionSet(x, x > 1, {z}), {2})`. The capture case
`ConditionSet(x, x < y, Reals).subs(y, x)` gives `x < x` in both.

Remaining gaps:

- Piecewise double integrals whose condition involves an inner integration
  variable with finite real limits, such as
  `(x**2 + Piecewise((0, x/t < 0), (1, True))).integrate((t, 4/5, 1), (x, -1, 1))`.
  SymPy returns 1/3, and also 1 and 1/2 for related integrands. SymPy
  replaces such variables with dummies that inherit assumptions (here
  `t > 0`). Symox moves the inner integral into the Piecewise branches,
  leaves `t` in the conditions, and returns an unevaluated `Integral`. The
  result is still correct.
- limit_seq returns `None` where SymPy finds 1, 1/2, 1/3 or oo on
  `harmonic` terms (`harmonic(n)/log(n)`, `n*Sum(1/k**2, (k, 1, n))`); see
  the limit_seq section at the end.
- `Sum(x*f[n], (n, 0, 2)).doit()` is the expanded `f[0]*x + f[1]*x + f[2]*x`
  (SymPy `x*(f[0] + f[1] + f[2])`). Native function terms behave the same
  way.
- `Sum(f[n], ...).free_symbols` lacks the base `f`, because Indexed objects
  are flattened into symbols. `rf(Poly, k)` with a symbolic `k` cannot
  hold the `Poly` argument. Indexed objects with a shaped base print in
  flattened form.
- `harmonic(12, 1/3)` prints `12**(2/3)/12` (SymPy `18**(1/3)/6`). The
  values are equal, and the difference comes from radical normalization.
- `sqrt(exp_polar(I*pi)).n()` keeps a spurious real part of about 1e-26.

## Infinities, numeric comparison, Factors and Boolean.equals on 2026-10-22

A probe of 562 core calls (Add/Mul/Pow/division tables over `zoo`, `oo`,
`-oo`, `nan` and `I*oo`, `comp`, `all_close`, `check_assumptions`,
`failing_assumptions`, `unchanged`, `WildFunction` matching, bound function
variables, `Factors`/`Term` and `Boolean.equals`) went from 107 to 14
outputs different from SymPy 1.14.0 (by `str` and exception type), and a
second probe of 762 `Pow`/`Add`/`Mul` table calls went from 143 to 43.
`tests/test_wrong_core_wave4_parity.py` holds the reference values.

- `Add` follows `Add.flatten`: a second `zoo` gives `nan`, and a term
  `c*oo*s` is collected as `s` with the coefficient `oo`, so `I*oo + I` is
  `oo*I`, `x + x*oo` is `oo*x` and `I*oo - I*oo` is `nan`.
- `Pow` follows `Pow.__new__` and the `_eval_power` of `Zero`, `One`,
  `NegativeOne`, the integers, rationals and Floats, `Infinity`,
  `NegativeInfinity` and `ComplexInfinity`: any base to the power `zoo` is
  `nan`; `b**oo` is `oo` for `b > 1`, `0` for `-1 < b < 1`, `zoo` for a
  finite `b < -1`, `nan` otherwise; `b**-oo` is `(1/b)**oo` with `nan` for
  a negative non-integer; `0**-oo` is `zoo`; `(-oo)**e` is
  `(-1)**e*oo**e`; `oo**e` uses `re(e)`; `(I*oo)**(1/2)` is `oo*sqrt(I)`.
  Exponents such as `I*oo` are decided through the finite coefficient of
  `oo`, and bases such as `pi` through the numerical sign of `pi - 1`.
- `comp` is the reference function (non-Float numbers have precision -1,
  so `comp(1e-10, 0)` is `True`); `all_close`, `check_assumptions`,
  `failing_assumptions` and `unchanged` are the reference functions.
- `WildFunction('F').match(F)` and `F.matches(F)` are `{F_: F_}`.
- `iterfreeargs` and the native `free_symbols` of a definite `Integral`
  mask a non-symbol integration variable with a `Dummy`, so
  `Integral(f(x), (f(x), 1, 2)).free_symbols` is empty.
- `Factors` and `Term` are the reference classes of `exprtools.py`
  (`gcd`, `lcm`, `quo`, `rem`, `div`, `normal` with negative and symbolic
  exponents).
- `Boolean.equals` is the reference test: same atoms and
  `not satisfiable(Not(Equivalent(a, b)))`, so `(A >> B).equals(~B >> ~A)`
  is `True`; with relationals it raises `NotImplementedError`.

Divergences kept (not validated by the user yet):

- `zoo**oo` stays unevaluated; SymPy gives `0`, although `|zoo| = oo`.
- The indefinite `Integral(f(x), f(x)).free_symbols` is `{x}`; SymPy gives
  `{f(x)}`, a non-symbol.

Remaining gaps:

- `oo/0` and `-oo/0` with a Python `0` are `zoo`; SymPy's
  `Infinity.__truediv__` gives `oo` and `-oo` (`oo/S(0)` is `zoo` in both).
- `Mul` with `zoo` keeps finite non-real factors (`zoo*(1 + I)`,
  `zoo*(1 + oo*I)`) and `(1 + oo*I)**2` stays a power; SymPy gives `zoo`.
- Print order of `zoo` in sums (`zoo + x` versus SymPy `x + zoo`).
- `WildFunction` is a `Wild` symbol: `F.is_Function` is `False` and
  `F.free_symbols` is `{F_}` (SymPy `True` and `{F}`).
- `Sum` and `Product` with a function as the summation variable raise
  `ValueError`.
- `all_close((1, 2), [1, 2])` raises `AttributeError` (SymPy
  `SympifyError`).
- `Pow(I, -S.Half)` is `1/sqrt(I)` (SymPy `-I**(3/2)`) and `(-pi)**1.5`
  stays unevaluated (SymPy `(-1.0*I)*pi**1.5`).

## Geometry, vectors, combinatorics and holonomic numerics on 2026-10-22

A probe of 160 calls (`is_concyclic`/`is_collinear`, `contains` of rays,
lines and segments in 2D and 3D, `Plane.angle_between`, vector `equals`,
`express` and `projection`, `Cycle`, `Permutation.atoms`,
`Subset.subset_indices`, the polycyclic `Collector` and
`HolonomicFunction.evalf`) went from 87 to 152 outputs (str or exception
type) equal to SymPy 1.14.0. `tests/test_wrong_geomcomb_wave4_parity.py`
holds the reference values.

- `Point.is_collinear` and `Point.is_concyclic` are instance methods that
  include `self`, as in the reference (`Point(1, 1).is_concyclic(Point(1, 0),
  Point(0, 1), Point(-1, 0))` is `False`; the unbound form still works).
  `is_concyclic` uses the reference rank and rref test, so coplanar spatial
  points are accepted, and `Point.affine_rank` is available.
- `Ray.contains(Ray)` and `Line.contains(Ray)` require both points of the
  inner ray on the outer entity (a non-parallel ray from a point of the ray
  was accepted).
- `Plane.angle_between` of a spatial linear entity is the signed
  `asin(n.d/(|n|*|d|))`, and is `None` for other arguments, as in SymPy.
- Vector `equals` compares `(a - b).dot(a - b)` with zero, so vectors of
  rotated frames compare equal (`D.j.equals(C.k)`). `projection(...,
  scalar=True)` is the reference coefficient `dot(other)/dot(self)`
  (7/3, not the Euclidean scalar projection), and projecting onto the zero
  vector gives zero instead of `ValueError`.
- Vectors and dyadics print as in the reference: a unit coefficient is
  omitted, a coefficient is parenthesized only if its text has `(`, `-` or
  `+`, Floats are printed short, terms are grouped by system name then
  base order, and the zero dyadic prints `0`. `Permutation` repr is the
  reference cyclic form `Permutation(2)(0, 1)`.
- `Cycle` is a port of the reference class: it stores the listed points
  only, iterates over the array form (`list(Cycle(1, 2))` is `[0, 2, 1]`),
  has the reference `str`/`repr` and composition, and rejects dict
  arguments. `Permutation.atoms()` is the set of the array-form elements.
  `Subset.subset_indices` returns `[]` when an element is missing.
- `Collector.subword_index` returns letter positions `(low, high)` or
  `(-1, -1)`; `map_relation(x_i*x_j)` returns the conjugate relator value
  of `x_j**-1*x_i*x_j`; `pc_presentation` uses the reference insertion
  order and stores an identity value as `()`.
- `HolonomicFunction.evalf(points)` with a list takes one Euler or RK4 step
  from each point to the next, starting at `x0` (`h` is unused), as in the
  reference. A single endpoint is still reached with steps of length at most
  `h`.

Intentional divergences (SymPy 1.14 demonstrably wrong; to be validated):

- `Ray((0, 0), (4, 4)).contains(Ray((-1, -1), (0, 0)))` is `False`: the
  reference only checks collinearity and the sense, and returns `True` for a
  ray that starts outside.
- `Collector.collected_word` keeps the whole power relator: with the S4
  sequence `[(0 1 2 3), (0 2 1), (0 1)(2 3), (0 3)(1 2)]` the reference keeps
  only the first syllable of `x0**2 = x2*x3` and returns `x0*x1**2` for
  `x1*x0**3` (a different permutation); symox returns `x0*x1**2*x3`, and
  `exponent_vector((0 2)(1 3))` is `[0, 0, 1, 1]`, not `[0, 0, 1, 0]`.
- `HolonomicFunction.evalf(b)` at a single point integrates to exactly `b`.
  The reference builds the mesh with `int((b - a)/h)` steps and returns the
  value at the last mesh point, which is not `b` when `(b - a)/h` is not an
  integer in floating point.

Remaining gaps:

- `PermutationGroup.polycyclic_group` picks a different polycyclic sequence
  than the reference for `SymmetricGroup(4)` and `DihedralGroup(4)`; the
  reference sequence comes from `derived_series` with randomized normal
  closures and is not even stable between SymPy runs. Collector results that
  depend on it (`map_relation`, `exponent_vector`) differ for those groups,
  while they agree for an explicitly given sequence.
- `Collector.induced_pcgs` and `constructive_membership_test` keep the
  symox algorithms: the reference `induced_pcgs([(0 1 2 3)])` in S4 omits
  the square of the generator, so its membership test rejects
  `(0 2)(1 3)`; symox returns a complete sequence and a boolean (the
  reference returns the exponent list or `False`).

## Precision and incomplete results from the wave-4 replay on 2026-10-22

The tests are in `tests/test_precision_incomplete_wave4_parity.py`,
`tests/test_wave4_substitution_parity.py`,
`tests/test_wave4_index_biomech_control_parity.py` and
`tests/test_wave4_h_validation_parity.py`.

### expm1, log1p and cancellation

`expm1` and `log1p` are function classes with the reference `eval`,
`fdiff`, `expand_func`, `rewrite` and assumption handlers. Before, the
constant arguments were turned into `exp(a) - 1` and `log(1 + a)`. With a
Float, both lost every digit: `expm1(1e-99).evalf()` and
`expand_log(log1p(1e-99)).evalf()` gave 0. A symbolic application did not
evaluate after `subs`, and `diff` stayed unevaluated.

- `expm1` of a Float, or of a number with Float parts, is computed with
  `mpmath.expm1` at the precision of the Float. Otherwise `exp(a)` is
  tried, so `expm1(0) = 0` and `expm1(log(3)) = 2`. `expm1(2)` stays
  unevaluated, as in the reference, and `evalf` uses `mpmath.expm1`.
- `log1p` of a Rational is `log(a + 1)`. A Float is converted exactly,
  so `log1p(1e-99)` is `log(Rational(1e-99) + 1)`, as in the reference.
  `evalf` uses `mpmath.log1p`.
- Accuracy. Relative error against mpmath at 400 digits, for
  `t = 10**-k` with `k` in 10, 50, 99 and 200, at 15, 30, 60 and 120
  digits: every form stays below `10**(1 - dps)`. The forms are
  `expm1(t)`, `log1p(t)`, `expand_log(log1p(t))`, the same three on
  `Float(t)`, and `subs` into a symbolic application.
- SymPy is less accurate here. `log1p(1e-99).evalf()` is 0 in SymPy. It is
  `1e-99` in symox. SymPy also gives 0 for `log(1 + 10**-99)`,
  `log(1 + t) - t`, `cos(t) - 1`, `sinc(t) - 1` and `tan(t) - t` at
  `t = 10**-99`. Symox has full precision on these, and on `Si(t)`,
  `erf(t)`, `erfc(t) - 1`, `sqrt(1 + t) - 1` and `exp(t) - 1 - t`.
  `gamma(-1 + 10**-50)` raises at 15 and 30 digits in both, like the
  reference. At 60 digits the symox error is `2e-21` (SymPy `4e-15`).
- `lambdify` with NumPy or math evaluates `expm1` and `log1p` through
  `numpy.expm1` and `numpy.log1p`. Before, it raised "cannot compile".
  `ccode` and `cxxcode` print `expm1`/`log1p`, and `cxxcode` now prints
  `std::expm1` like the reference.

### lambdify with mpmath

`lambdify(..., 'mpmath')` and `modules=['mpmath']`, with optional dict
namespaces, now follow SymPy. The expression is printed with mpmath
calls and evaluated in the namespace of `from mpmath import *` and
`MPMATH_TRANSLATIONS`. So the callable computes with `mpf` at the
current `mp.dps`. Before, the native evaluator returned a double: at
`mp.dps = 40`, `lambdify(x, sin(x), 'mpmath')(1)` had only 16 correct
digits.

- Floats are printed as their exact binary mantissa and Rationals as
  `mpf(p)/mpf(q)`. `pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio` and
  `oo` are mpmath constants. Square roots are printed with `sqrt`.
  `Piecewise`, `sinc`, `Heaviside`, `Min`/`Max` and matrices are printed
  as in the reference (`mpmath.matrix`).
- 80 expression and namespace probes give the same values and types as
  SymPy, at 38 digits.
- A list with `math`, `numpy` or `sympy` keeps the native evaluator, as
  before. SymPy also computes doubles there. An expression that the
  printer does not cover also stays on the native evaluator.
- `MpmathPrinter` prints Floats as `mpmath.mpf((sign, man, exp, bc))`,
  Rationals without extra parentheses and roots as `mpmath.sqrt`. It uses
  `mpmath.inf`, `mpmath.euler` and similar constants. 12 of 14 probes
  print the same text as the reference. For `LambertW`, `Si` and
  `conjugate`, symox prints `mpmath.lambertw`, `mpmath.si` and
  `mpmath.conj`, where SymPy raises.

Remaining gaps:

- The native `evalf` of `log(1 + r)` loses every digit when the Rational
  `r` is below about `10**-230`, and returns 0. SymPy already returns 0
  from `10**-99`. `log1p(r)` follows `log(r + 1)`, as in the reference, so
  it is affected too.
- `series`, `limit` and `integrate` do not handle `expm1` and `log1p` (a
  native gap for function classes). SymPy gives
  `x + x**2/2 + x**3/6 + O(x**4)`, `limit(expm1(x)/x, x, 0) = 1` and
  `integrate(expm1(x), x) = -x*exp(x) + x*expm1(x) + exp(x)`. `pycode`
  raises where SymPy prints `math.expm1(x)`.
- `sympy.codegen.scipy_nodes` (`cosm1`, `powm1`) has no symox module.

### Index methods, musculotendon constants and gbt

- `get_contraction_structure` and `get_indices` follow
  `sympy.tensor.index_methods`. A product, power or function whose factors
  hold contractions gets a nested entry, so
  `get_contraction_structure(x[i, i]*y[j, j])` has the key
  `x[i, i]*y[j, j]` with `[{(i,): {x[i, i]}}, {(j,): {y[j, j]}}]`, and
  `A[j, j]**A[i, i]` is handled the same way. Dummy indices keep their first
  occurrence order. An index that occurs three times gives no outer index,
  as in the reference, instead of an error. 19 probes, including every
  docstring example, give the same keys and entries as SymPy.
- `MusculotendonDeGroote2016.p` lists the 31 reference symbols in the
  reference order: the six muscle parameters, the activation constants, then
  `c_0:4_fl_T`, `c_0:2_fl_M_pas`, `c_0:12_fl_M_act` and `c_0:4_fv_M`.
  The tendon curve constants are listed even for a rigid tendon. Default
  symbols carry no assumptions. `MusculotendonDeGroote2016.with_defaults`
  exists, and the explicit formulations accept a nonzero fiber damping and
  neglect it in the force balance, like SymPy. The implicit formulations,
  which SymPy does not implement, are kept.
- `gbt` converts `alpha` with `Rational(alpha).limit_denominator(1000)` and
  does not simplify, so `gbt(tf, T, 0.5)` gives the exact
  `T/(2*(L + R*T/2))`. `bilinear`, `forward_diff` and `backward_diff`
  give the same `srepr` as the reference on four transfer functions.

Remaining gaps:

- `ccode`/`fcode` of a contraction with `assign_to` and ranged `Idx` print
  a single assignment instead of the reference summation loops (fixed in the c3 section).
- The biomechanics curves are evaluated (`doit`) in the force and in
  `rhs`, where SymPy keeps unevaluated curve objects. The values are equal.
  For formulation 2 without defaults, SymPy leaves a `Derivative`
  unevaluated in `rhs`, where symox gives a number.

### substitution and nonlinsolve with transcendental equations

- `substitution` is a port of the reference algorithm. It runs a real pass
  and a complex pass, keeps the `ImageSet` families and tracks complements,
  intersections and excluded denominators. So
  `substitution([exp(x) - sin(y), y**2 - 4], [x, y])` keeps the
  `2*n*I*pi` families and gives the reference output.
- `nonlinsolve` falls back to the reference path (unrad, the polynomial and
  non-polynomial split, Groebner, then `substitution`) when the native
  parametric solver rejects a system that is not polynomial. Polynomial
  systems keep the native path.
- `ImageSet` objects inside a `Tuple` have `lamda`, `base_set`, `base_sets`
  and `base_pset`.
- 38 exp, log and trigonometric systems, each through both functions, went
  from 74 to 11 results that differ from the reference. All the reference
  docstring examples of `substitution` give the reference output.
- Symox is deliberately stricter in three cases where the reference is
  wrong. These are intentional corrections:
  - `[sin(x)**2 - 1, y - 1]`: the reference returns `(x, 1)`. Symox
    returns the families `2*n*pi + pi/2` and `2*n*pi + 3*pi/2`.
  - `[exp(y) - x, x**2 + 1]`: the reference returns `(exp(y), y)`. Symox
    returns `(I, 2*n*I*pi + I*pi/2)` and `(-I, 2*n*I*pi - I*pi/2)`.
  - `[2**x - 8, y - x]`: the reference returns `(y, y)`. Symox returns
    `(3, 3)`.
- A row that leaves a variable free is dropped when an equation still
  restricts that variable. If a solver call was unsolved, a `ConditionSet`
  is returned instead. Rows that differ only in the names of their
  parameters are merged. A finite row that no free row covers is kept,
  where the reference drops it.

The complex `solveset` gaps found here (`exp(x) - exp(-x)`,
`exp(2*x) - 3*exp(x) + 2`, `2**x - 8` and `sin(x) + cos(x)`) are fixed; see
"Complex `solveset` of exponential and trigonometric equations" at the end.

## Mechanics, vector, control and continuum mechanics on 2026-10-22

A probe of 149 calls (`Body`, `RigidBody`, `Particle`, `dynamicsymbols`,
`msubs`, `CoordinateSym`, the wrapping geometries and pathways, the margins,
`pole_zero_numerical_data`, `StateSpace.dsolve`, `TransferFunctionMatrix`
indexing, `Beam3D` and `Cable`) went from 104 to 14 outputs different from
SymPy 1.14.0 (by `str` and exception type).
`tests/test_wrong_mechanics_wave4_parity.py` holds the reference values.

- `RigidBody` and `Body` without inertia carry the symbolic central inertia
  `B_ixx, B_iyy, B_izz, B_ixy, B_iyz, B_izx`, so `kinetic_energy` keeps the
  rotational term. Default masses are plain `Symbol(name + '_mass')`.
  `Body` with a mass and no inertia is a particle (`is_rigidbody` is
  `False`) and uses the particle kinetic energy and parallel axis.
- `dynamicsymbols` passes its assumptions to the undefined functions.
- `msubs` is the reference algorithm (`_crawl`, `_sub_func`): it does not
  substitute inside derivatives. `smart=True` uses `_smart_subs` with
  `_fraction_decomp` and the `tan` rewrite.
- `CoordinateSym(name, frame, index)` returns the coordinate of `frame`, so
  `CoordinateSym('a_y', A, 1) == A[1]`.
- `Vector` and `Dyadic` print like `_sympystr`: unit coefficients are
  dropped, signs are pulled out (`- F*sin(q)*N.x`), and sums are wrapped in
  parentheses.
- `WrappingCylinder` and `WrappingSphere` are the reference classes:
  `point_on_surface`, the `ValueError` messages, `geodesic_end_vectors`,
  and a cylinder `geodesic_length` through `_directional_atan` (no `Mod`).
  `WrappingPathway.to_loads` applies `force*pA_force`, `force*pB_force` and
  `-force*(pA_force + pB_force)` at the axis point.
- `phase_margin` and `gain_margin` follow the reference: `solveset` over
  `Interval(0, oo, left_open=True)`, `(arg*180/pi + 180) % 360` folded
  below 180, `-180` with no crossover, and `oo` with no phase crossover.
  Non-SISO systems, time delays and extra symbols raise the reference
  errors.
- `pole_zero_numerical_data` returns `(zeros, poles)` like `_poly_roots`:
  exact roots through `all_roots` over QQ/ZZ, otherwise `nroots`, with
  real roots as `float` and others as `complex`.
- `StateSpace.dsolve` returns the output `C*x + D*u` and has the
  reference `ValueError` and `ShapeError` checks.
- `TransferFunctionMatrix.__getitem__` has the matrix semantics:
  `tfm[:, 0]` is a column `TransferFunctionMatrix`, a pair of integers is
  an entry and a flat index is row major.
- `Beam3D` is the reference Timoshenko formulation: `dsolve` of the shear
  and bending equations, then `linsolve` for the integration constants
  with the slope and deflection boundary conditions. Load magnitudes go
  into `load_vector` and `moment_load_vector`.
- `Cable.apply_length` sets `length`, raises `ValueError` below the support
  distance, and a distributed solve keeps the applied length. Supports keep
  their input order and `reaction_loads` start at zero.

Divergences kept (not validated by the user yet):

- `Beam3D(..., variable=y)` uses `y`; the reference passes the variable to
  `Beam.__init__` and then uses `x` anyway.
- `Cable.apply_length` compares with the support distance
  `sqrt(dx**2 + dy**2)`; the reference uses `dx**2 - dy**2`.
- `Cable.solve(x0, y0)` with distributed loads checks `y0` against the
  supports and raises `ValueError` if they disagree; the reference ignores
  the second argument. Distributed reactions stay balanced (the reference
  uses the left tension at both supports).
- `WrappingSphere.geodesic_length` with a symbolic angle gives
  `r*acos(cos(q))`; the reference raises `ValueError` because it cannot
  prove that the points are on the sphere.
- Exact forms differ but are equal for the margins (for example
  `20*log(3)/log(10)`), `StateSpace.dsolve` and `Beam3D` deflections, and
  for the order of `Cable.tension` keys.

Remaining gaps:

- Frame coordinates are `Dummy` atoms: they print as `_A_x`, `is_real` is
  `True` (SymPy `None`) and they have no `.frame` attribute.
- `LinearPathway(...).to_loads(F)[0].vector.simplify()` has another
  simplified form.
- `solveset` leaves `log(Abs(G(I*w)))` equations as a `ConditionSet`, and
  algebraic roots inside `Interval(0, oo)` as an unevaluated
  `Intersection`. The margins solve `|G|**2 - 1` instead and decide
  membership numerically at 50 digits; ties that cannot be decided raise
  `NotImplementedError`.
- `Beam3D` plotting and the `max_*` methods are not ported.

## Quantum states, Wigner symbols and second quantization on 2026-10-22

Three probes compare outputs with SymPy 1.14.0 by `str` and exception
type. `tests/test_wrong_quantum_wave4_parity.py` holds the reference values.

- Quantum states, spin representation and the Wigner/CG grid: 161 of 235
  calls matched before, 234 after.
- Rotation, tensor-product and spin siblings: 18 of 34 before, 30 after.
- `secondquant`: 6 of 88 before, 76 after (85 when term order in sums is
  ignored).

- `TensorProduct` follows `TensorProduct.__new__`/`flatten`:
  - The commutative parts are pulled out with `args_cnc`, so
    `TensorProduct(x*A, y*B)` is `x*y*AxB`.
  - A single argument returns itself.
  - Factors that are sums, products or powers are printed in parentheses.
  - Products of tensor products of the same length combine factorwise:
    `(A*C)x(B*D)`.
  - `represent` of a tensor product of kets is the Kronecker product.
- `IntQubit` follows `IntQubitState._eval_args`:
  - A single integer is converted to binary.
  - With two arguments, the second one is `nqubits` when it is larger
    than 1.
  - `IntQubit(n, nqubits=k)` is supported.
  - A `Qubit` argument is accepted.
  - `IntQubit(2, 1)` raises `ValueError`.
  - It prints as `|n>`, and `as_int()` is a method.
  - `IntQubitBra` is its dual.
  - Gates keep the qubit class.
  - `superposition_basis` and the Grover functions use
    `IntQubit(i, nqubits=n)`, so an oracle on `IntQubit(k)` marks the
    intended state.
- Spin states follow `SpinState._represent_base`: `represent(state,
  basis=...)` is the column of `Rotation.D(j, m, mp, ...).doit()` with the
  Euler angles of each axis pair.
  - The default basis is the state's own axis.
  - A bra is represented by the conjugate transpose of its ket.
  - `JzKetCoupled` and the other coupled kets give the direct-sum column
    of dimension `prod(2*ji + 1)` and print `|j,m,j1=..,j2=..>`.
  - Spin operators are represented only in the `Jz` basis; any other basis
    raises `NotImplementedError`.
  - `WignerD.doit` is the reference `_eval_wignerd`.
  - `qapply(Rotation*JzKet(j, m))` sums `D(j, m, mp)*|j,mp>` with `m` as
    the first index.
- `Wigner3j(j1, m1, j2, m2, j3, m3).doit()` passes its arguments to
  `wigner_3j(j1, j2, j3, m1, m2, m3)`; `Wigner3j(6, 0, 4, 0, 2, 0).doit()`
  is `sqrt(715)/143`.
  - The `wigner` functions accept Floats with half-integer values.
  - `wigner_9j(..., prec=n)` passes `prec` to `racah` and returns the
    unsimplified sum.
- `enumerate_states(state, start, num)` and `(state, index_list)` are the
  reference function: `XKet('foo')` gives `|foo_1>`, `|foo_2>`, and so on.
  The other states give `[]`, and a wrong number of arguments raises
  `TypeError`.
- `diff(Wavefunction(...), x)` differentiates the expression and keeps the
  limits. `OperatorExpr.diff` applies the product rule. `Wavefunction.limits`
  is a dict.
- `TimeDepKet`/`TimeDepBra`: `label` excludes the time, `time` is the last
  argument, and they print `|psi;t>` and `<psi;t|`.
- Sums of states print with the native order and signs: `|0>/2 + |1>/2`
  and `-|2>`.
- `secondquant` uses the Fermi vacuum of the reference:
  - The fermion operators carry the Fermi-level and quasi-particle
    properties of the reference.
  - `contraction` introduces the `_a`/`_i` dummies, and a fermion with a
    boson raises `ContractionAppliesOnlyToFermions`.
  - `NO` sorts its operators with the sign and is 0 on a Pauli violation.
  - `NO(...).doit()` splits general `Dummy` indices into the parts above
    and below the Fermi level, and raises
    `SubstitutionOfAmbigousOperatorFailed` for other ambiguous indices.
  - `wicks`, `evaluate_deltas`, `PermutationOperator.get_permuted`,
    `simplify_index_permutations` and the `SecondQuantizationError`
    exceptions are the reference ones.
  - `wicks(F(0)*Fd(0), keep_only_fully_contracted=True)` is
    `KroneckerDelta(0, _a)` (it was 1, under an empty-vacuum convention).
  - Operators print as `f+(p)`, and as `CreateFermion(p)` inside products.

Divergences kept (not validated by the user yet):

- `wigner_9j(1, 1, 1, 1, 1, 1, 1, 1, 0, prec=64)` ends in `...556`, the
  correctly rounded value of `1/18`; SymPy prints `...555`.
- `represent(Rotation(...), j=...)` evaluates the `WignerD` entries; SymPy
  keeps them unevaluated.
- Bosons: `NO(...).doit()`, `wicks` and `contraction` still give the
  correct bosonic results, where SymPy raises `NotImplementedError` or
  `ContractionAppliesOnlyToFermions`.
- `Fd(p)*Fd(p)` and `NO(Fd(p)*Fd(p))` are 0; SymPy leaves
  `CreateFermion(_p)**2`.
- `wicks(3)` is 3; SymPy raises `AttributeError`.

Remaining gaps:

- `qapply(TensorProduct(Jz, Jz)*TensorProduct(JzKet(1, 1), JzKet(1, -1)))`
  raises `TypeError` (operator expressions are not hashable as tensor
  factors); SymPy gives `-hbar**2*|1,1>x|1,-1>`.
- `Bra('a')*Ket('b')` prints `InnerProduct(Bra(a), Ket(b))`; SymPy prints
  `<a|b>`.
- `TensorProduct` distributes over sums when it is constructed; SymPy keeps
  `(A + B)xC` until `expand(tensorproduct=True)`.
- Coupling or uncoupling more than two spins (`couple`/`uncouple` with
  three `jn`) is not checked against the reference.
- Symbols with the `above_fermi`/`below_fermi` assumptions cannot be
  created (unknown assumption in the core). Only the dummies made by
  `secondquant` itself have a Fermi level, so a `KroneckerDelta` built by
  the user does not vanish across the Fermi level.
- `substitute_dummies`, the index ordering of `AntiSymmetricTensor` and
  Fermi-sea `FKet` (`fermi_level != 0`, holes) are not ported.
- Sums of operator products print `-1*NO(...)` and `+ -...`; SymPy prints
  `-NO(...)` and `- ...`.

## Remaining blocking gaps

### Noncommutative symbols

`Symbol('A', commutative=False)`, `Dummy`, `Wild` and
`Function('g', commutative=False)` are accepted, and an undefined function
applied to a noncommutative argument is noncommutative, as in the reference.
A feature bit marks expressions that contain such an atom, so commutative
expressions keep their native fast paths and their outputs unchanged.

Supported with reference behaviour (`tests/test_noncommutative_parity.py`
holds the reference outputs):

- `is_commutative` of symbols, products, sums, powers and functions, and the
  fact closure (`is_real`, `is_complex` and `is_zero` are False).
- `Mul` canonicalisation as in `Mul.flatten`: commutative factors are
  collected in front, noncommutative factors keep their order, and adjacent
  equal bases combine (`A*A**-1 == 1`, `A**x*A**y == A**(x + y)`).
  `A*B - B*A` stays nonzero, `(A*B)**2` stays unexpanded and
  `(2*A)**2 == 4*A**2`. A number times a noncommutative sum stays a
  product (`-(A + B)`, `2*(A + B)`) and only `expand` distributes it.
- `args`, `args_cnc`, `as_coeff_Mul`, `as_ordered_factors`, `coeff`,
  `as_independent`, `has` on subproducts, `str`, `srepr`, `latex`, equality,
  hashing, pickling, `subs` of a contiguous subproduct (a port of
  `Mul._eval_subs`), `xreplace`, `expand` (multinomials keep the factor
  order; `exp(A + B)` is not split), `diff`, `series` and `nc_simplify`.
- `factor` and `factor_nc` (ports of `factor_nc` and `_mask_nc`),
  `gcd_terms` and `factor_terms` (the reference term masks), `together`,
  `cancel` and `apart` (the reference `PolynomialError` branches), `ratsimp`
  (raises `PolynomialError`), `collect` (raises `AttributeError` for a
  noncommutative pattern), `combsimp` and `gammasimp` without gamma
  factors, and `simplify` (the reference pipeline, starting with
  `nc_simplify`). Commutative subexpressions still go to the native code.
- `Poly`, `linsolve` and the native expression-to-polynomial conversion
  raise `PolynomialError('non-commutative expressions are not supported')`.
- `solve` for commutative unknowns replaces each equation by the
  coefficient equations of its distinct noncommutative words, so a
  commutative unknown only takes commutative values
  (`solve(A*x - B, x) == []`, `solve(x*A - A, x) == [1]`). A single
  noncommutative unknown is solved from a one-sided linear equation
  without commuting factors (`solve(A*B - C, B) == [A**-1*C]`).
  `nonlinsolve` of one equation in one unknown uses the same solver.
- `adjoint` and `transpose` (ports of the reference classes and of the
  `_eval_adjoint`/`_eval_transpose` rules of `Mul`, `Add`, `Pow`, `exp`,
  `Piecewise` and `conjugate`): `adjoint(A*B) == adjoint(B)*adjoint(A)`,
  hermitian and antihermitian symbols, `conjugate(adjoint(A)) ==
  transpose(A)`, `subs`/`xreplace` re-evaluating the nodes, `Expr.adjoint`,
  `Expr.transpose` and `Matrix.H`/`Matrix.adjoint` of noncommutative entries
  (`tests/test_noncommutative_rest_parity.py`).
- `radsimp` reproduces the final numerator/denominator step of the
  reference (`radsimp(A/(1 + sqrt(2))) == -A + sqrt(2)*A`).
- `Symbol`, `Dummy` and `Function` with inconsistent assumptions (including
  `commutative=False` with `real=True`) raise `InconsistentAssumptions`, a
  `ValueError` subclass, as in the reference. `srepr(sin)` prints `sin` and
  `srepr(Function('q'))` prints `Function('q')`.
- `limit`, `nsimplify`, `evalf` and the matrix products, determinants,
  traces and transposes of noncommutative entries match the reference on
  the probed cases.

Stages of `simplify` without a noncommutative reference algorithm
(`trigsimp`, `exptrigsimp`, `hyperexpand`, `besselsimp`, `logcombine`,
`expand_log`, `radsimp`) run on a linear mask: the noncommutative part of
each term becomes a commutative Dummy, and the result is kept only when it
is still linear in these Dummies with none in a denominator. Otherwise the
stage is skipped, so no factor is commuted.

Remaining noncommutative gaps:

- `solve` and `nonlinsolve` differ from the reference where the reference
  result is wrong or where Symox can decide more (see the intentional
  corrections at the top): `solve(A*B + B, A)` is `[-1]` (reference `[]`),
  polynomial equations of higher degree in a commutative unknown are solved
  (`solve(x**2*A - A, x) == [-1, 1]`, `solve(exp(x)*A - A, x) == [0]`)
  where the reference raises `PolynomialError`. Two-sided equations in a
  noncommutative unknown and systems in noncommutative unknowns raise
  `NotImplementedError`.
- `solveset(A*x - B, x)` returns the reference's `{B*A**-1}`, inconsistent
  with `solve(A*x - B, x) == []`; `solveset` has no noncommutative branch.
- `integrate(A*x, x)` returns `x**2*A/2`; the reference raises
  `ValueError`.
- `radsimp(x*A/(1 + sqrt(2)))` is `-x*(-sqrt(2)*A + A)` (reference
  `x*(-A + sqrt(2)*A)`, the same value): `-(A + B)` distributes over a
  noncommutative `Add` in Symox, where the reference keeps
  `Mul(-1, A + B)`, so `could_extract_minus_sign` breaks the tie
  differently (and `(A + B) - (A + B)` is 0 where the reference keeps
  `A + B - (A + B)`).
- `adjoint` and `transpose` of a function of a noncommutative argument stay
  unevaluated (`adjoint(sin(A))`; the reference
  `transpose(sin(conjugate(A)))` through `Function._eval_conjugate`).
  `latex` and `pretty` print the nodes as functions
  (`\operatorname{adjoint}{\left(A \right)}`, reference `A^{\dagger}`).
- `srepr` of a symbol lists its assumptions sorted
  (`Symbol('A', antihermitian=True, commutative=False)`); the reference
  keeps the keyword order. `srepr` of an undefined function head with
  assumptions lists all derived facts (reference `Function('g')`).
- The message of `InconsistentAssumptions` lists the given facts where the
  reference lists the facts deduced up to the conflict.
- `Matrix([[A, 0], [0, B]]).inv()` returns `diag(A**-1, B**-1)`; the
  reference raises `PolynomialError`.
- `Expr.equals` of noncommutative operands returns True only for a zero
  difference and None otherwise.
- Noncommutative `Wild` matching beyond the probed cases, the
  `assumptions` system (`ask`), matrix algorithms of noncommutative entries
  beyond the probed ones, the printers other than `str`, `srepr` and
  `latex`, and the noncommutative
  user classes with `is_commutative = None` (below) are not audited against
  the reference. Polynomial simplification still treats noncommutative host
  factors of such user classes as commuting.

None of the cross-call reuse paths below (exact arithmetic, monomial `Add`
construction, substitution, expansion, numerical evaluation and the NumPy
float64 kernels) accepts an expression containing a host node of a user
class: such nodes dispatch their hooks (`_eval_subs`, `_eval_evalf`,
`_eval_power`, assumption handlers...) on every call and own Python instances
that must stay collectable. Monomial `Add` reuse also requires each argument
to be exactly the nominal SymPy class native wrapping gives its node, so user
subclasses and other wrappers keep full constructor dispatch.

The [public performance runner](../benches/README.md) checks accuracy and result
quality before ranking timings. Performance work on 2026-10-20 corrected
`nroots`' formerly ignored precision: rational-coefficient roots now use native
arbitrary-precision Aberth corrections after exact square-free decomposition.
The public function and `Poly.nroots` accept positional/keyword `maxsteps` and
`cleanup`; exhausted iteration raises mpmath's `NoConvergence` rather than
silently returning machine-precision estimates. Starting guesses may use
machine arithmetic; the computed answers do not. A Newton-polygon initializer
also covers root magnitudes outside binary64. Tests cover 15-100 digits,
all roots of the selected polynomials, multiplicities, zero roots and extreme
scales. These tests do not certify every condition number or coefficient domain.
General algebraic, transcendental and complex coefficient domains in `nroots`,
all argument/error edge cases, certified complex enclosures and comparative
worst-case convergence remain open obligations.

The same work corrects an overly conservative expansion bound for univariate
polynomials: the degree-plus-one support bound replaces exponential counts of
distributive products when the input grammar establishes that bound. Other
grammars retain their existing checked fallback. Reuse of completed expansion,
substitution and numerical calculations is bounded and includes the relevant
precision, options, physical binding identities and numeric-query policy.
Cached expansion replays its work reservation and retains distinct term limits.
Exact number-theory caches clone returned factor maps; public mutations cannot
alter cached answers. Rational dense polynomial construction avoids an
expression round trip, retaining exact domains and the established fallback
for other coefficient domains. These are scoped preservation arguments, not a
refresh of the separate full API or downstream qualification records.

`expand` has no default size limit (2026-10-22), as SymPy's has none. A
default term limit of 100 000 and a work budget of 5 000 000 elementary
products used to make SymPy-compatible calls fail, e.g. `expand((x + y +
z + w + v)**40)` (135 751 terms) or `expand((x + y)**3000)`, and the internal
callers that fell back to the unexpanded input on that error (simplification
helpers, `collect`, `integrate`, the inverse Laplace transform, the domain
construction of expressions) silently changed their results. The default
budget is now unbounded, for scalar and matrix expressions, and DiracDelta
expansion no longer stops at 100 000 roots. The `max_terms` hint is a
symox-only, opt-in resource guard: with it, exceeding the term limit or the
work budget raises `ValueError` and never returns a partial expansion. Only
an arithmetic overflow of the term counts (exponents near `2**64`) still
fails without a hint. `expand(complex=True)`, `expand_complex`,
`as_real_imag`, `re` and `im` have no default size limit either since the
same date, see the section "Complex expansion without limits and
Heaviside(x, H0) on 2026-10-22".

The small integer response caches for `primepi` and `partition` apply only to
one positional built-in Python integer in `[-10**6,10**6]` while evaluation is
enabled. This is cache eligibility, not a mathematical domain limit: larger
inputs, other accepted types, keywords and held evaluation retain the original
factory. Existing public head objects, aliases, signatures and pickle identity
are preserved. Each cache holds at most 32 immutable computed results; separate
thread-local evaluation contexts cannot retrieve an evaluated value while
requesting held output.

The versioned per-function ledger retains its own native snapshot. The
downstream candidate does not renew those historical certificates. Current
upstream and regression results are recorded separately in the
[integration report](downstream-integration.md) and its candidate-specific
validation evidence.

The [performance fix report](algorithm-performance.md) also records exact
Mersenne perfect-power rejection, bounded reuse of canonical exact monomial
`Add` inputs and a NumPy float64 elementary-expression path. These fast paths
retain their established fallback for inputs outside their eligibility rules.
NumPy compilation is charged to `lambdify`, array results remain independent of
their inputs, and unsupported printer output keeps the native evaluator.
Matching selected NumPy results or reusing identical exact computations does
not establish comparative complexity or full backend/API parity.

Follow-up frontend optimizations reuse bounded immutable lexical tokens rather
than evaluated expressions: namespaces, user transformations, callbacks and
evaluation options remain active on every parse. Scalar expansion retains the
native option and budget checks while bypassing unrelated container dispatch.
The exact involution shortcut for permutation powers first verifies every
image, applies to arbitrary signed integer exponents and retains subclass
constructor/inversion behavior. These are scoped preservation arguments;
timings remain separate from the open mathematical parity obligations below.

The downstream fixes also cover exact Python `Fraction` conversion of native
integers and rationals, two-dimensional NumPy/matrix arithmetic, external
`_sympy_` conversion hooks, and relational constructor identity and instance
checks. Cyclic conversion hooks raise an exception before exhausting the native
stack. These protocols do not establish general Python subclass compatibility:
custom symbolic subclasses, their metadata and container type predicates remain
open obligations. Dense matrix subclasses can be constructed and override Python
methods; preservation of those subclasses across every matrix operation is not
yet certified.

Python subclasses of `Symbol`, `Dummy`, `Function`, `Expr`, `AtomicExpr`,
`Number`, `Basic` mixed with `Expr`, and `Singleton` classes are now real
subclasses of the native expression type. Each instance is a native node whose
host handle names its class and is weakly bound to the Python instance;
expressions containing such nodes keep those instances alive, so instance
attributes survive `args`, `free_symbols`, `subs`, `xreplace`, arithmetic and
matrices. Equality and hashing follow `_hashable_content` when a class defines
it, construction through `Expr.__new__(cls, *args)`, `Basic.__new__`,
`Symbol.__xnew__` and Singleton registration follow the reference, and pickling
uses the reference `__getnewargs__`/`__reduce_ex__` protocol. A class
`_eval_power` is consulted by `**` and `Pow` once the generic exponent rules
leave the power unevaluated; `Add` and `Mul` fold their `Number` terms through
the arithmetic of a `Number` class that defines it, in the reference flatten
order, and `_as_mpf_val` gives `float()` and, as the reference
`Number._eval_evalf`, `evalf` (an infinite or undefined value evaluates to the
`oo`, `-oo` or `nan` singleton).

Native algorithms dispatch to the class hooks where the reference does:
`_eval_subs` (before the native rewrite of the node), `_eval_is_*` assumption
handlers (combined with the native deduction rules for `Add`, `Mul`, `Pow`,
`Abs` and the fact closure, and not cached after a raise), `fdiff` and
`_eval_derivative`, `_eval_evalf`, `_eval_rewrite_as_*`, `_eval_expand_*`
including user hints, `_eval_simplify`, overridden `func` and `free_symbols`,
`_eval_power` and `_op_priority`. An exception raised by a hook propagates as
the original exception object. Canonical argument order compares host class
names as the reference `ordering_of_classes` does. A class with
`is_commutative = None` (or a false or undetermined commutativity handler) is
noncommutative: `Mul` keeps such factors after the commutative part in their
written order, `Add` and printers keep the reference order, and powers,
expansion, substitution and differentiation respect it. Expressions without
such nodes take the unchanged native path. `Wild` subclasses are hosted and
keep their class through matching. Host instances stay alive exactly while an
expression or a Python object references them, reference cycles through
instance attributes are collected, and nodes rebuilt by a native operation go
through the class constructor again.

Remaining gaps: hook exceptions raised inside a native algorithm are reported
when the algorithm returns to Python rather than at the raising point, so
partial work of that call is discarded instead of interrupted; the
`frac`, `numer`, `denom` and `modulus` expand hints use the native expansion;
polynomial-based simplification (`cancel`, `together`, `factor`, `ratsimp`,
`radsimp`, `simplify` and `nc_simplify`) treats noncommutative host factors as
commuting; native `subs` does not match a pattern against a noncommutative
subsequence of a product; native tree hashes do not include a host
`_hashable_content` (this only causes hash collisions, equality is exact);
the reference shares one fact base between symbols with equal assumptions,
which can hide a `Symbol` subclass `_eval_is_*` handler after a plain-symbol
query, and caches `Function` applications so `super().__new__` may return a
shared object - neither aliasing is reproduced; a `Wild` subclass `matches`
override is not consulted by the native matcher; `Basic.__new__` sympifies its
arguments, whereas the reference stores raw Python values; nested Python
`_print_<Class>` printer methods are not consulted by native printers;
native operators on a tree that already contains such a `Number` do not refold
its coefficient (`(w*x)*3`).

Native nodes are instances of nominal classes that follow the reference
hierarchy. Native results equal to `oo`, `-oo`, `nan`, `zoo`, the Boolean atoms
and the small rationals are the `S` singleton objects. Function heads such as
`sin` and `Abs`, the Boolean gates and `Max`/`Min` are classes. They support
`isinstance`, `issubclass`, `singledispatch` registration, `LatticeOp`
`make_args`/`zero`/`identity`, pickling and user subclasses that evaluate,
differentiate and substitute. The abstract bases are importable from their
reference modules, for example `TrigonometricFunction`, `HyperbolicFunction`,
`ExpBase` and `MinMaxBase`. Relations have `reversed`, `reversedsign`,
`negated`, `canonical`, `weak`, `strict`, `lts` and `gts`. Boolean atoms have
`negated`, `canonical` and `_op_priority`. Arithmetic on a Boolean atom
raises the reference `TypeError`. Arithmetic and ordering on another non-`Expr`
operand (a relation, a Boolean function, a predicate, a set or a tuple) are
not supported, as in the reference. Evaluated inequalities reject non-real,
NaN and non-`Expr` sides with the reference messages. Sign facts decide an
ordering only when the difference is extended real, and a numeric difference
of symbolic sides orders them only when both sides are extended real, so
`s < s + 1` stays unevaluated for a complex `s` (`is_ge`). Sums, products and powers
of finite terms are finite, so `1 + I` is not extended real. Scalar and matrix
`kind` follow the reference rules, and `srepr` of a `Float` names its
precision. The `kind` of a tree with shared subexpressions is computed once
for each distinct node. A nominal class compares equal to its public class
and hashes like it, so `type(expr) == Eq` holds and dictionaries keyed by
public classes find `type(expr)`. `exp(x).exp` is `x` and `exp(x).base` is
`E`. Applications of user `Function` subclasses to numbers are numbers, as in
the reference, while `Function('f')(1)` is not. The `Basic` class flags such
as `is_symbol`, `is_NumberSymbol`, `is_Not` and `is_scalar` are defined. These
classes have remaining gaps:

- Head constructors expose an `evaluate` keyword in their signature, which the
  reference does not show.
- A user subclass of `sin` evaluates at a `Float` argument (`S2(1.0)`), where
  the reference keeps it unevaluated.
- The private kind classes are `_ScalarKind` instances, not the reference
  `_NumberKind` and `_BooleanKind` classes. Set kinds are not implemented.
- `is_complex` of a finite non-real sum such as `1 + I` is now True (see
  the core gaps section at the end). The finiteness of `I**(1/3)` is
  unknown.
- Noncommutative symbols are implemented; see the noncommutative section
  above for their remaining gaps.
- `Mul.is_polar` of nonpolar factors is unknown. The reference gives False.
- Binary float exponents are `i128` values bounded by `2**100`, so
  `exp(2**63)`, `exp(-2**63)`, `exp(2**70)` and `Float('1e999999999999999999')**3`
  evaluate, compare and print as in the reference (Torch value ranges compare
  `exp(9223372036854775807)`). Decimal conversion of huge exponents uses a
  200-bit logarithm estimate instead of exact powers of ten. Doubles that
  overflow or underflow (`exp(Float(710))`, `cosh(Float(1000))`) are
  re-evaluated with multiprecision, and `Max`/`Min` compare such values at 96
  bits. Beyond `2**100` (for example `exp(Float('1e100'))`, which the
  reference evaluates) results stay symbolic, and the argument reduction of
  `sin`/`cos` is not attempted above `2**(2**30)`.
- `Eijk` is a class, whereas the reference `Eijk` is a function that returns
  a `LeviCivita` node.
- `printer`, `str`, `precedence`, `codeprinter`, `c` and `cxx` are ports of
  the reference modules, so downstream subclasses of `C99CodePrinter` and
  `CXX11CodePrinter` (Devito's) override `_print_<Class>`, `type_mappings` and
  precedence as in the reference. The other code printers keep the earlier
  implementation and are not subclasses of the ported `CodePrinter`; `fcode`
  and `rust_code` print natively. Those printers reject the augmented
  assignment statements (`x *= y`) that the reference prints. The C printers
  extend the reference in one case: a `DotProduct` of vectors with integer
  lengths prints as the sum of element products, where the reference raises.
  Float digits come from a port of mpmath's `to_str`, without importing
  mpmath.
- `codegen.ast` nodes are not `Basic`: their equality, hashing and `func`
  follow the reference, but `atoms`, `xreplace` and pattern matching on
  statement trees do not.
- `IndexedBase`, `Indexed` and `Idx` are not `Basic` and become symbols in
  arithmetic, so `Idx` ranges are lost, `assign_to` does not generate
  contraction loops over repeated indices, and Devito's
  `IndexedBase.__new__(cls, ...)` subclasses fail.
- The public `Add`, `Mul`, `Basic`, `Symbol`, `Function` and `Eq` are routing
  classes distinct from the nominal node classes, so `cls.__base__` of a user
  subclass differs from the reference (Devito recurses through
  `cls.__base__.__new__`).
- Arithmetic with a string operand sympifies the string, where the
  reference raises `TypeError`.
- The strict-mode error message of the code printers names `_core.Matrix` for
  explicit matrices.
- Measured on 2026-10-20 after integrating the subclass, host-hook,
  performance and singleton-class work: Devito's `tests/test_symbolics.py`
  passes 27 of 136 tests (before: none, the import failed); most remaining
  failures need instance attributes on `Relational` and `Add` subclasses
  (`_subdomain`, `.d`) and the items above. Torch's `test/test_sympy_utils.py`
  had 34 failures (reference: none): 33 were the `exp(2**63)` exponent range,
  now closed (1 failure, 235 passed, 146621 subtests), and 1 needs the user
  `SingletonInt._eval_Eq` hook (`Eq(SingletonInt(1), SingletonInt(2))` stays
  unevaluated instead of `False`).

N-dimensional array subclasses follow the reference construction protocol
(`__new__` through `cls._new`, then `__init__`, with `_array`, `_shape`,
`_rank`, `_loop_size` and `_handle_ndarray_creation_inputs`), and arithmetic,
slicing, `applyfunc`, `diff`, `reshape`, `transpose` and the array operations
keep the subclass. Dense matrix subclasses select result classes with the
reference `classof` and `_class_priority` rules for joins, stacking, element
products, solvers, inserts and arithmetic. `subs` and `xreplace` accept
function classes as in `Function._eval_subs`: applications whose argument
count is accepted by the replacement are rebuilt bottom up, `exp` replaced by a
non-function becomes a power, `cot` and `CRootOf` keep their head, unordered
rules use the reference node-count ordering, and `xreplace` ignores
function-class keys. Unknown `subs` options are ignored hints, as in the
reference; only `hack2=True` remains explicitly unsupported. A proxy object
forwarding attribute access to an expression is returned by `sympify`
unchanged and converts to that expression in native operations. In
`lambdify`, a user namespace `sinc` receives the argument printed by the
selected reference printer (`x/pi` for the NumPy family) and a namespace
`sqrt` evaluates `b**(1/2)` and `1/b**(1/2)`, as the reference prints them.
A power of a power with a real inner exponent of magnitude below one combines
(`sqrt(sqrt(t))` is `t**(1/4)`), as in `Pow._eval_power`. Remaining gaps:
the other branches of `Pow._eval_power` (an inner exponent `-1` on a negative
base, an even inner exponent on an imaginary base, a base with nonnegative
real part and inner exponent below two in magnitude, and the branch factor of
a half-integer outer exponent) are not reproduced, so
`(1/n)**(1/2)`, `(n**(3/2))**(1/2)` for negative `n` and
`((1 + I*r)**(3/2))**y` stay nested; `lambdify` with `math` or `mpmath` evaluates natively, so complex
arguments to `math` functions do not raise and `mpmath` results are double
precision values rather than `mpc`/`mpf`; array `shape` is a Python tuple rather than `Tuple`; `Derivative`,
`integrate` and `Indexed` on arrays, symbolic array indices (a native
`Piecewise` rather than `Indexed`) and array plus matrix arithmetic (accepted
natively, rejected by the reference) differ; lambdify of integer arrays and
matrices returns float64 where the reference keeps int64; global
`transpose`/`adjoint`/`conjugate` and `QRsolve` return mutable matrices where
the reference returns immutable ones; `ImmutableMatrix` subclasses are not
preserved; `Matrix.D` is missing; function-class substitution into
`Integral`, `Sum` and `Product` raises a `TypeError` with a different message;
simultaneous substitution with a `Lambda` replacement succeeds where the
reference raises `IndexError`; `Lambda.nargs` is a raw native set node. SymPy 1.14 builds a
dense array from a sparse one out of its stored flat positions; symox keeps
the entries.

Physical quantities print through `StrPrinter`: its `abbrev` setting, and the
class default `StrPrinter._default_settings['abbrev']` that Lcapy sets, print
`V` instead of `volt` in `str`, `repr` and `sstr`, as in the reference.
`StrPrinter` has its own `_default_settings`. `has(Quantity)`,
`atoms(Quantity)` and `isinstance` recognise the atom that represents a
quantity inside an expression, and `get_dimensional_dependencies` of a
dimension name expression (`'length/time'`) returns the reference powers.
Remaining gaps: a quantity is represented by a positive `Symbol` of the same
name, so a user `Symbol('volt', positive=True)` is also recognised as a
quantity; when two quantities in one expression share an abbreviation, or an
abbreviation equals another symbol's name, the full name is printed; quantities
sort with ordinary symbols (`V*volt` where the reference prints `volt*V`, and
`volt/ampere` has its arguments in another order); `get_dimensional_expr`
returns base dimensions rather than the named dimensions; `base_dims` order and
`Dimension` symbols differ.

Rational bounded-denominator approximation uses the same exact standard-library
`Fraction.limit_denominator` kernel as the reference, including its tie breaking
and integer conversion of the bound. Scalar arithmetic and expansion preserve
the canonical Python objects for zero, one, negative one and one half. Scalar
`Iterable` and `Sized` predicates distinguish tuples from ordinary expressions.
Tuple results use a native-expression subtype with sequence slots, including
results of substitution, evaluation and serialization. Tuple indexing and
length use constant-time native access; iteration materializes arguments once,
then visits them in linear time. Scalar iteration raises immediately. Other
container protocols and general symbolic subclass preservation remain open.

The OctSymPy bridge now has native Boolean category queries, original Symbol
assumptions, explicit matrix representations, MatrixBase stacking, entrywise
global integration and simultaneous scalar substitution in matrix dimensions
and coefficients. Global matrix integration returns immutable output, while
the matrix method retains the receiver's mutability. Public condition pairs
validate their expression/condition arguments and expose the reference
`(expr, cond)` signature, typed equality, reconstruction and pickle behavior.
`Basic` is now a structural container and category, rather than an alias for
scalar `Expr`. Virtual registrations include native scalar nodes, symbolic
matrices and immutable explicit matrices. Sympification retains an existing
AlgebraicNumber's exact field metadata, as in the reference. These category
checks do not change the concrete native Python MRO or certify every Basic
method. Algebra involving
direct formal Boolean/Basic applications remains unsupported by the scalar
kernel; mixed simultaneous matrix/scalar substitution remains explicit
`NotImplementedError`. General nominal inheritance, printer settings and
matrix substitution patterns remain open. Boolean-valued Piecewise conditions
with an explicit default convert to ITE before native branch storage; scalar
rewriting uses the same conversion. Exhaustiveness proofs without a default
remain an open obligation.

Native scalar nodes are now instances of nominal classes named after their
SymPy classes (`Add`, `Integer`, `sin`, `Derivative`...), generated from the
reference class graph with the same names and MRO. These classes carry no
state and are subtypes of the native expression type. Downstream converters
that dispatch on `type(expr)`, its string, or MRO base names (SageMath's
`_sage_` methods, Cadabra's `SympyBridge`) therefore see SymPy's classes.
Under the downstream routing facade they also report SymPy's module paths.
Attributes patched onto public classes after import, including `Basic`, are
forwarded to these classes, so later patches follow the MRO with no stale
cache. Native-written attributes stay authoritative. The same work adds
`Poly.rep` as a dense `DMP` view, SageMath's
`_integer_`/`_rational_`/`_polynomial_` element hooks for ZZ, QQ and
polynomial-ring domain elements, `RealInterval`/`ComplexInterval` with
`CRootOf._get_interval`, `Printer._as_ordered_terms`, and a persistent
`Symbol.name` string (`x.name is x.name`, as in SymPy). Constant-coefficient
linear ODEs take their characteristic roots from exact `roots`, as SymPy does
(`X'' + m**2*X = 0` gives `exp(+-I*m*tau)`, not `exp(+-sqrt(-m**2)*tau)`).
The gamma, beta, harmonic and combinatorial families (`gamma`, `lowergamma`,
`uppergamma`, `polygamma`, `loggamma`, `digamma`, `trigamma`, `factorial`,
`subfactorial`, `factorial2`, `rf`, `ff`, `binomial`, `beta`, `harmonic`) now
rewrite with SymPy's rules, conditions and hints (`piecewise`, `limitvar`),
following SymPy's bottom-up traversal and pattern selection.

Measured with the routing facade against the unmodified upstream selections:
SageMath's `sage/interfaces/sympy.py` doctests went from 182 to 8 failures out
of 319, and Cadabra's `tests/sympy_cdb.cdb` passes with all 12 success
markers. The remaining SageMath failures are open gaps: complex root
isolation (3), polynomial-ring domains over polynomial rings such as
`ZZ[x][y]` (2), `CRootOf` arguments stored as a coefficient tuple instead of
SymPy's `(PurePoly, index)` (2, also visible in `str` and `srepr`), and the
Fresnel (Meijer G) integral of `sin(x**2)` (1). Related rewrite gaps remain:
`subfactorial(x).rewrite('factorial')` keeps an unevaluated `Sum` because the
alternating `summation` with a symbolic bound is not closed, and the error
function family (`erf` to `uppergamma`, `expint` or `tractable`) is not
rewritten. Other
differences found along the way: `sstr(expr, order=...)` ignores the order
setting, the matcher can choose a different valid assignment for ambiguous
patterns, `Poly.intervals` uses Sturm isolation rather than continued
fractions, and numeric evaluation of complex `CRootOf` stays unevaluated.
Tuples, sets, matrices, `Poly`, matrix expressions and `Lambda` are not yet
instances of nominal SymPy classes. The nominal classes derive from the host
tree type of user subclasses (described above), so an expression containing
user subclass instances reports its SymPy class and keeps those instances
alive. Nominal classes stand only for symox's own public classes: a user
class reusing a SymPy name (Torch's `Max`, Devito's `Derivative`) does not
see the attributes patched onto the node class of that name; user subclasses
see those patched onto their nearest SymPy ancestor (`Symbol`, `Function`,
`Expr`) and its bases, as their SymPy MRO would give.

The explicit `NonElementaryIntegral` constructor builds a held native integral
and retains its marker for reconstruction, pickle and the covered root
structural operations. Constructing this marker does not prove that its
integrand is non-elementary. `integrate` produces it only from a completed
proof of the ported `risch_integrate` (exp and log towers, as in the
reference): with `risch=True`, in the Risch stage of the default chain on
the whole integrand, and in the Risch stage of the term loop. A computation
that no method decides stays a plain `Integral`. The marker is kept by node
identity: an equal `Integral` built elsewhere is not marked, the dummy
substitution of a definite integral keeps it, and other rebuilds of a
composite native expression may drop it.

Boolean system constraints are converted to their equation meaning before
linear solving: true is an identity, and false proves inconsistency. The public
`solve` interface reports an empty assignment list for an identity or when no
requested variable is assigned, including dict/set output conventions. This
interface result is distinct from a mathematical proof that the solution set
is empty. Undecided parameter constraints retain their established core paths.

Native distribution integration now has an exact rule for a single simple
real affine Dirac delta multiplied by a locally analytic weight. It substitutes
the exact support point in the weight and integrates through Heaviside, with
the signed slope retained in the primitive. A substitution budget rejects
expensive forced coefficient growth before evaluation. The new rule declines
inexact input, singular or discontinuous weights and products of deltas; all
established integration methods remain available. Derivatives
`DiracDelta(a*x+c, n)` of literal order `n` integrate by repeated parts
against an analytic weight `g`, giving
`sum_j binomial(n, j)*(-1)**j*g^(j)(x0)` terms with the scale
`1/(a**n*Abs(a))`. SymPy 1.14 keeps only one term, omits the binomial
coefficient and the scale (for example `x*DiracDelta(x, 2)` gives
`-DiracDelta(x)` instead of `-2*DiracDelta(x)`), so these results differ from
the reference by design where the reference is wrong. Where the reference
returns `DiracDelta(0)`, `zoo` or `DiracDelta(oo)` for a derivative at an
endpoint, the integral stays unevaluated.
This does not certify general distribution integration or floating root
precision. Direct Heaviside/DiracDelta constructors now validate provably
nonreal arguments, and DiracDelta validates its derivative order. The same
validation after arbitrary substitution and reconstruction remains open.

The `diracdelta=True` expansion hint uses a complete-root guard for square-free
rational polynomials: the native solver must supply as many distinct exact
roots as the polynomial degree, and each root's reality must be decided before
any support is discarded. Multiple roots and undecided domains retain the
original delta. Affine expansion additionally requires a proved real nonzero
slope and a real support point. Generic parameter strata, higher derivative
orders and nonpolynomial root sets remain open. No new asymptotic superiority
claim is made for root isolation or expression growth.

Two distribution differences are intentional and are not reference matches.
For `DiracDelta(x**5-x-1)`, SymPy 1.14's incomplete radical-root dictionary
causes expansion to return zero. Symox retains the exact indexed real root;
an independent Sturm count and minimal-polynomial check certify its support.
Finite interval endpoints use the declared symmetric half-weight convention,
including at the origin. SymPy returns zero for some weighted deltas at a
finite origin endpoint while returning half weight after translation or at
an infinite endpoint. Translation and a normalized Gaussian regularization
are checked separately; this convention is not a universal definition of
distribution integration at discontinuous boundaries.

Improper definite integrals follow the reference order in `Integral.doit`.
The variable is replaced by a dummy that carries the sign the bounds imply.
For extended real bounds with an infinite end, `meijerint_definite` (now a
port of the reference stages, conditions and `hyperexpand`) is tried before
any antiderivative. Its convergence condition is kept: `conds='piecewise'`
returns `Piecewise((value, cond), (Integral, True))`, `'separate'` returns
the pair and `'none'` drops it. A parameter range where the integral
diverges therefore never receives the analytic continuation as its value.
`meijerg=True` integrates only through Meijer G (an unevaluated Integral
otherwise) and `meijerg=False` skips it. The `risch`, `heurisch` and `manual`
options follow the reference `_eval_integral` (see the integration methods
paragraph below). An improper
bound limit that is not settled keeps the Integral unevaluated. An
antiderivative jump at an interior pole or log branch point is removed
through one-sided limits (`Expr._eval_interval`), so for example `tan` over
`(0, pi)` gives nan as in the reference.

Piecewise integrands follow the reference `Integral.doit` per limit. Abs and
sign are rewritten over a real dummy (indefinite: for a real variable;
definite: for real finite bounds), Min and Max become Piecewise, and an
integrand containing a Piecewise is folded. A Piecewise integrand goes
through `Piecewise._eval_integral` (indefinite) or `piecewise_integrate` and
the split evaluation of `Piecewise._eval_interval` (definite); the Meijer G
stage is not tried on it, as in the reference. A condition holding a
Piecewise is folded before the ITE rewrite, as `ExprCondPair` does.
`Piecewise.piecewise_integrate` and `Undefined` are available. On the
reference `test_piecewise.py`, assertion runs go from 522 passing / 232
failing to 642 / 103, with no wrong value among the failures. Remaining
differences:
- an interval the inequality solver cannot decide (for example the periodic
  `sin(x) >= 0` of `Abs(sin(x))` over `(0, 2*pi)`) keeps the integral
  unevaluated, where the reference returns 4;
- a piece whose antiderivative is held gives a held definite integral;
- `1/(x+a)` over `(0, 1)` with an unassumed `a` stays unevaluated (the
  reference gives `-log(a) + log(a + 1)`);
- condition forms with ITE, Min and Max differ in some results
  (`test_piecewise_integrate1ca/cb`, `test_piecewise_integrate2`,
  issue 6900), and the Python `piecewise_fold` is not the native port (it
  orders the folded pieces differently).

Integration methods. `risch_integrate` (with `prde`, `rde`),
`manualintegrate`, `heurisch`/`heurisch_wrapper`, `ratint`,
`trigintegrate` and `deltaintegrate` are ports of the reference modules,
and the native integrator calls them in the reference order of
`_eval_integral`:
- `risch=True` runs only `risch_integrate`: a proven non-elementary part is
  a `NonElementaryIntegral`, NotImplementedError (algebraic or trigonometric
  integrands) leaves the factored `Integral`, as `Integral.doit` gives it,
  and a definite integral raises ValueError, as in the reference;
- `manual=True` and `heurisch=True` run `_eval_integral` with only that
  method (the reference port, including the `Order` terms, `Piecewise`
  and `sincos_to_sum` recursion); a True option turns the others off, two
  True options raise ValueError;
- a False option removes the method from the default chain: after the
  native rules, `risch_integrate` on the whole integrand
  (`result + i.doit(risch=False)`, or `NonElementaryIntegral(f).doit(
  risch=False)` when nothing was integrated), then per term the linear power
  shortcut, `ratint`, `trigintegrate`, `deltaintegrate`,
  `singularityintegrate`, `risch_integrate`, `heurisch`, `meijerint` and
  `manualintegrate` (whose leftover integrals are done with
  `manual=False`), then the native term rules;
- the undone result is `Integral(f, x).factor()`, so
  `integrate(y*exp(x)/log(x), x)` is `y*NonElementaryIntegral(exp(x)/log(x),
  x)` as in the reference; a limit `(x, b)` with one bound is factored the
  same way, so `Integral(c*g, (x, b)).doit()` gives `c*Integral(g, (x, b))`
  and integrals that differ by a constant factor share one held integral.
  Gap: the reference then factors that product again
  (`factored_function = function.factor()`), and its `factor` expands the
  held integral over the terms of its integrand
  (`AddWithLimits._eval_expand_basic`): for
  `Integral(3*(3*f(x) - 2*g(x))*exp(x)/13, x).doit()` the reference returns
  `3*(Integral(3*f(x)*exp(x), x) + Integral(-2*g(x)*exp(x), x))/13` where
  symox keeps one held integral of the factored integrand; symox
  `Integral.expand` does not split the integrand yet;
- `Derivative.doit` does its result again when it changed and still holds
  a derivative, as the reference, so the integrals created by
  differentiating under the integral sign are done (and factored). The
  substituted PDE of `checkpdesol` with the `1st_linear_constant_coeff_Integral`
  hint then cancels in `simplify` as in the reference, but that `doit`
  runs `manualintegrate` on held integrands of `Subs` derivatives, about
  8 s per check in both libraries (SymPy 1.14: 7.5 s);
- an `Integral` integrand (not a `NonElementaryIntegral`) is denested into
  the limits, as the reference constructor does: `integrate(Integral(y, x),
  x)` is `x**2*y/2`, and `risch=True` raises ValueError for an inner definite
  limit;
- with a True option, the exceptions of the selected port (ValueError,
  PolynomialError, TypeError, ...) reach the caller as in the reference.
The ports rely on these core rules of the reference, which are also
followed: a symbol is the power `y**1` for an old power of it in `subs`, so
`y.subs(1/y, u)` is `1/u` (manualintegrate finds `u = 1/y` with it); a
`Derivative` whose variable is replaced by a non-symbol is held in
`Subs(Derivative(f(_xi_0), _xi_0), _xi_0, value)`, and one whose new
variable is already present is held in `Subs(expr, old, new)`
(`tests/test_subs_symbol_powers.py`); the first derivative of a power is
`self*(b'*n/b)` built with binary products, so `diff(sqrt(-x), x)` is
`sqrt(-x)/(2*x)`; heurisch keeps the repeated unknowns of the reference
`poly_coeffs` (the reference ring maps each to its last generator).
`log(a).is_zero` is `(a - 1).is_zero` as in `log._eval_is_zero`, so
`(y + 1)**x` with a zero `y` integrates to `x` in every method.
The expansion rules used by the ports follow the reference:
`expand(trig=True)` and `expand_trig` expand `sec`, `csc`, `sech` and `csch`
through their partner (`1/cos(...)` expanded), `tan`, `cot`, `tanh` and
`coth` of sums and integer multiples with the symmetric-function and
binomial formulas, `sinh` and `cosh` with the sum formula, and
`expand_mul` rebuilds `n/d` after splitting a fraction, so
`expand(x**2*y/(2*x*(x*y/2 + x - y/(2*x))))` is `x*y/(x*y + 2*x - y/x)`
(`tests/test_expand_trig_reference.py`).
Differential corpus: the 294 integration examples of `test_risch.py`,
`test_manual.py` and `test_heurisch.py`, run with SymPy 1.14 and symox as
the literal test call and as `integrate` in the default, `risch=True`,
`manual=True` and `heurisch=True` modes (1470 runs, 20 s each, the
timeouts rerun in a fresh process). Measured on 2026-10-22: the printed
form is identical in 1374 runs and equal by differentiation or
simplification in 72, 20 runs give another form, 3 stay unevaluated (all
`integrate(f(z), z)`, below), no symox run exceeds 20 s (the reference
exceeds it once); no symox result fails the derivative check where the
reference passes it, and the 10 `NonElementaryIntegral` results are proven
in both libraries (none in only one). The cases fixed in this measurement
are `tests/test_integrate_rest_parity.py`. The agreeing fast cases are `tests/test_integrate_methods_corpus.py`
(1072 cases; `heurisch(sin(x)*cos(x))` is left out, its reference form
depends on the hash seed), and the hand-picked ones are
`tests/test_integrate_methods.py`. Two corpus inputs are right here and
wrong in the reference: `integrate(f(z), z)` with `z = z(t)` (the reference
treats `f(z(t))` as a constant and returns a Meijer G sum; symox leaves the
`Integral`) and `x*exp(a*x + 1)` (below). `Integral.doit` adds the reference floor terms to every antiderivative,
indefinite or evaluated between bounds: an `atan` whose argument is linear
in `tan(a)` gets `sign(c)*pi*floor((a - pi/2)/pi)` (`floor(a/pi)` for
`cot(a)`), so `integrate(1/(cos(x) + 2), x)` is continuous and its integral
over `(0, 2*pi)` is `2*sqrt(3)*pi/3` (it was 0 before, a wrong value). The
reference adds the terms a second time to the antiderivative of a
Piecewise, which makes `integrate(Piecewise((1/(cos(x) + 2), x > 0), (0,
True)), (x, -1, 7))` wrong there (`4*sqrt(3)*pi/3 + ...`); symox adds them
once (an intentional difference). On an exp/log tower where the native
rules give an `atan` of the tower, the reference `risch_integrate` runs
first and its `RootSum` is returned (`exp(x)/(1 + exp(2*x))`,
`1/(x*(1 + log(x)**2))`). Remaining gaps:
- the cases that exceeded 20 s in the previous measurement
  (`atan(x)*log(x)`, `test_issue_10847` and its slow variant,
  `heurisch=True` on `test_issue_6799` and on `x**(4/3)*log(x)`) now finish
  within the limit; `atan(x)*log(x)` takes about 5 s (reference 8 s);
- `integrate(exp(10*x)*sin(exp(x)), x)` finishes here and times out in the
  reference;
- `solve_lin_sys` over an inexact ground (`RR(k10)`, as `heurisch` builds
  it for a Float rate) row reduces with the arithmetic of the domain like
  the reference `_solve_lin_sys_component`, instead of `linsolve` over
  expressions, where the float rounding kept a spurious common factor
  (`(k10 - c)/(k10**2 - 2*c*k10 + c**2)`, whose denominator `solve` then
  split into two complex roots depending on `s`, giving conditions on the
  integration variable). `integrate(exp(k10*s)*exp(-1.66666666666667*s), s)`
  is now the reference `Piecewise` with the generic case under
  `Ne(k10, 1.66666666666667)` and the degenerate rate last
  (`tests/test_transforms_wave3_parity.py`). The top-level `cancel` of an
  expression with 53-bit Float coefficients follows the reference
  `dmp_cancel` over `RR` (`tests/test_forms_rest_parity.py`): `RR` has no
  associated ring, so no denominators are cleared; the gcd is computed over
  `QQ` after `to_rational(limit=True)` of each coefficient, the cofactors
  come back through `RR.from_QQ`, and a negative denominator leading
  coefficient is moved to the numerator. `cancel(x/(0.3*x + 0.1))` is
  `1.0*x/(0.3*x + 0.1)` and the generic piece of the integral above is
  `1.0*exp(k10*s)/(1.0*k10*exp(c*s) - c*exp(c*s))`, where symox used to print
  `2251799813685248.0*...`. Intentional divergence: the degenerate piece is
  float noise in the reference (`-300239975158033.0*exp(-3.33066907387547e-15*s)`)
  and is now the exact `s` in symox. Remaining: `Poly.cancel` over `RR` still
  clears denominators (`3.6e16*x`) and `include=False` returns the integer
  `1` where the reference has `1.0`; floats with another precision keep the
  legacy cancellation; the definite integral of the float Piecewise of
  `test_issue_15498` is not rechecked here;
- `manual=True` on the rational function of case 153 leaves an `Integral`
  where the reference finishes (manualintegrate leftovers); `manual=True` on
  `atan(x)*log(x)` gives the reference result, about twice as slowly;
- the native rules also run under the reference `risch_integrate` base
  case (`integrate(..., risch=False)`), so rational parts keep the native
  partial-fraction form: `risch=True` on
  `(4*x**4 + ...)/(x**6 + ...)` gives `8*x/(4*x**2 + 4) + 3*atan(x) -
  1/(x**2 + 1) - 3/(x + 1)` (reference `(-x**2 + x - 4)/(x**3 + x**2 + x +
  1) + 3*atan(x)`);
- an unevaluated sum is kept as one integral: `(f(x) - f(-x))/x` gives
  `-Integral((f(-x) - f(x))/x, x)` (reference `-Integral(f(-x)/x, x) -
  Integral(-f(x)/x, x)`), likewise for a Bessel sum with `risch=True`, and the
  `NonElementaryIntegral` of `(exp(x) - exp(-x))/x` holds
  `(exp(2*x) - 1)*exp(-x)/x` (reference `(1 - exp(-2*x))*exp(x)/x`);
- held radicands are factored: `risch=True` on `1/sqrt(4 - 4*x**2)` gives
  `Integral(1/sqrt(-(x - 1)*(x + 1)), x)/2` (reference
  `Integral(1/sqrt(1 - x**2), x)/2`), likewise `3/sqrt(1 - 9*x**2)`,
  `1/sqrt(4*x**2 - 4)` and `heurisch(1/sqrt(3*x**2 - 4), x, hints=[])`;
- the variable order of a held mixed `Derivative` under `risch=True` and
  `heurisch=True` on `f(x, y).diff(y).diff(x)` differs;
- `sqrt(a*x + b)/(x + c)` with real symbols gives an unsimplified
  equivalent form, and `exp(-I*x**2/2)` over `(0, 1)` gives
  `sqrt(2)*sqrt(pi)*exp(-I*pi/4)*erf(sqrt(2)*exp(I*pi/4)/2)/2` (reference
  `(1 - I)*sqrt(pi)*erf(1/2 + I/2)/2` expanded), equal values;
- `ratint(1/(x**3 + x + 1), x)` finishes with a sum of logarithms and
  arctangents with radical coefficients, equal in value to the reference
  but in a larger form (about 80 KB printed against 17 KB) and about five
  times slower. The cause is the numeric radical canonicalization of
  products: `Integer**Rational` keeps `1922**(2/3)` where the reference
  extracts `31*124**(1/3)`, `sqrt(3*2**(2/3))` gives
  `sqrt(3)*sqrt(2**(2/3))` (reference `2**(1/3)*sqrt(3)`, also visible in
  the arctangent coefficient of `ratint(1/(x**3 + 2), x)`), and `Mul` does
  not combine numeric radicals by the order-dependent gcd rule of
  `Mul.flatten`. A port of these rules matched the reference on isolated
  products but made the `EX` division of `log_to_atan` on this integrand run
  for more than ten minutes, so it is not applied;
- `pdsolve(x*u(x, y).diff(x) - y*u(x, y).diff(y) + y**2*u(x, y) - g(x, y))`
  runs for more than 100 s on a tree that merges the integration options
  into the second integration wave (the reference takes 27 s). The time goes
  to `heurisch` on `g(x, eta/x)*exp(-eta**2/(2*x**2))/x`, whose ansatz has
  the reference size (1717 coefficients): the walk `find_non_syms` over
  `raw_numer` is slow. The slowdown remains with the RootSum, polynomial and
  expansion changes of this branch reverted, so it comes from the merged
  base, not verified against the later `heurisch` fix of the second wave;
- `test_issue_6799` over `(-pi, pi)` gives `pi*cos(n*phi)` in the default
  chain (reference the unsimplified sum
  `-(-pi*n*cos(n*phi)/2 - sin(n*phi)/2)/n + (pi*n*cos(n*phi)/2 -
  sin(n*phi)/2)/n`, equal value);
- the native rules run before the reference chain, so an integrand they
  integrate keeps the native form: `x*exp(x) + exp(x**2)` gives
  `(x - 1)*exp(x) + sqrt(pi)*erfi(x)/2` (reference `x*exp(x) - exp(x) +
  ...`), and a definite integral that stays unevaluated is not split into
  the reference sum of one-sided integrals
  (`y*Integral((x - 1)*(x + 1)*exp(x)/log(x), (x, 1, 2))`);
- `integrate(x + 1/log(x), (x, 1, 2))` gives `li(2) + oo` (reference `oo`):
  numeric values of special functions such as `li(2)` have no sign here,
  where the reference decides it by evaluation;
- an exception other than NotImplementedError inside a ported method makes
  that stage give no result in the default chain, where the reference lets
  some of them out of `integrate`;
- `deltaintegrate` is not used on the terms of the default chain that hold
  a `DiracDelta`: the native delta rules (above) stand in for it, so the
  singular products the reference assigns a value (`DiracDelta(x)/x` gives
  `zoo*Heaviside(x)`) stay unevaluated; `heurisch=True` keeps the reference
  stage;
- after the reference chain fails, the native term rules can still find an
  antiderivative the reference does not give;
- `final=False` is accepted and ignored, and with `manual=True` or
  `heurisch=True` the native definite stages (rational antiderivatives,
  gamma moments) are still tried before the antiderivative;
- `ratint(1/(x**5 + x + 1))` and similar `RootSum` results are slow (also
  in the reference);
- a derivative of total order above one is not passed through the
  reference `factor_terms(signsimp(...))`, and a product is differentiated
  n times by repetition instead of the general Leibniz rule of
  `Mul._eval_derivative_n_times`: `diff(x**n, x, 2)` is
  `n**2*x**n/x**2 - n*x**n/x**2` (reference `n*x**n*(n - 1)/x**2`), equal
  values in another form;
- `heurisch_wrapper(sqrt(x**2/((y - x)*(y + x))), x)` gives the reference
  antiderivative with the radicand written `x**2/(-x**2 + y**2)` (reference
  `-x**2/(x**2 - y**2)`).
The ports depend on `cancel`: float coefficients outside the double range
(`1e-300*u`) no longer underflow to zero or overflow to `oo`, but the forms of
float `cancel` results still differ from the reference RR domain
(`1/(0.2*x + 0.4)` is printed with scaled integer coefficients).

`meijerint_indefinite` is ported (`_meijerint_indefinite_1`, the splitting
points, the hyperbolic rewrite and the choice of the simplest result) and
is the last stage of indefinite integration unless `meijerg=False`.
`meijerg=True` uses only it. The public `meijerint_indefinite` raises the
errors the reference lets out, with their classes: the inverse Mellin
transform errors of the recursive rewrite other than
`IntegralTransformError` (`NotImplementedError('Gammas partially over the
strip.')` for `sin(x)**3/x**3`, a second `MellinTransformStripError`),
`Unexpected form...`, the `_get_coeff_exp` and `meijerg` `ValueError`s and
the `hyperexpand` errors. Inside `integrate` they give no result and the
next method runs; the reference does that for `NotImplementedError` only,
so an integrand where the method raises another class would differ (none
is known).

Polynomial times `exp(a*x)` antiderivatives use the Risch form
`num*exp(a*x)/den` and keep `Piecewise((.., Ne(den, 0)), (value at a = 0,
True))` when the coefficient may vanish. Polynomial times sin or cos, and
products with both an exponential and a trigonometric factor, keep the cases
of `heurisch_wrapper`: the parameter values that zero the denominator are
solved for, joined when there are several, and integrated separately, for
example `Eq(a, 0) & Eq(b, 0)`, `Eq(a, -I*b)`, `Eq(a, I*b)` for
`exp(a*x)*sin(b*x)`. The cyclic rule no longer divides by a vanishing
`a**2 + b**2` (it returned `zoo` terms for `exp(I*x)*sin(x)`); those products
are integrated in exponential form. Known differences:
- the special-case antiderivatives come from the exponential form, not the
  heurisch form (equal up to a constant);
- for `x*exp(a*x + 1)` the reference special case is `x**2/2`, which is
  wrong; the native value is `E*x**2/2`;
- for `x*sin((a + 1)*x)` heurisch fails in the reference and the result has
  no `Ne(a, -1)` case; the native result keeps it.

The incomplete gamma functions, `Ei` and `expint` apply the reference branch
factors for arguments carrying `exp_polar(2*I*pi*n)`, and `expand_log`
recurses and repeats until the expression is stable, as the reference does.
`solve` takes even powers out of the radical of a quadratic root
(`a**2 + b**2 = 0` gives `-I*b, I*b`, like `roots_quadratic._sqrt`).
Like the reference `_solve`, it leaves unsimplified the solutions of a
polynomial in one generator of the unknown that has more than two roots,
unless `simplify` is given; the roots of each factor of a product are decided
separately. For an exact polynomial of degree 3 and above in one generator
the native solver takes the `roots` formulas (Cardano, Ferrari, binomial,
with the `cubics`, `quartics` and `quintics` flags), and `all_roots` when they
are incomplete: linear, quadratic and binomial factors give radicals
(`CRootOf._roots_trivial`), the others give `CRootOf`. A non-real binomial
root is expanded in complex form like `roots_binomial` (`x**6 + x**3 + 1`).
The generator is inverted for every root through `solve(g - t)`, and, like
the reference `decompose_power`, `exp(p*x/q + c)` with `|p| > 1` is a power
of `exp(x/q)`, so `exp(2*x) - 2` and `exp(3*x + 1) - 2` give every branch
(the `I*pi` solutions were missing); `b**(p*x/q + c)` with `b` free of the
unknown is decomposed the same way, also with several such powers of one
base (`2**(3*x) - 8` has three solutions, `2**(2*x) - 3*2**x + 2` gives
`[0, 1]`, where the native solver found one solution or failed). A positive
real root of a power of `b` is also inverted on its own so that an exact
rational exponent is kept (`2**x - 8` gives `3`, as the reference simplify
gives it). The sign analysis of `solve_univariate_inequality` takes the
roots of the numerator as exact zeros: with the Cardano forms it no longer
closes the intervals of `x**3 - 3*x + 1 > 0` at the roots (the reference
raises `NotImplementedError` for this inequality). A product is solved factor by factor,
each with its own simplification (`(x**3 - 2)*(x**2 - a**2/(a + 1))` gives
`a*sqrt(1/(a + 1))`), the union is filtered by the denominators and the
product is not simplified again. Like the reference `_solve`, the union is
not checked against the whole product: only the roots of the simple
denominators (`_simple_dens`, the base of a power with a numeric exponent)
and the assumptions of the unknown are rejected, so `solve(x*log(x), x)` is
`[0, 1]`, `solve((x - 1)*log(x - 1), x)` is `[1, 2]` and
`solve(x*log(x)/(x - 1), x)` is `[0]`. Like the reference `roots` with
`incomplete=True` (the default of `solve`), a polynomial whose roots are
incomplete and whose domain `all_roots` does not support (symbolic,
algebraic or transcendental coefficients) gives the partial result of
`roots`: `solve(x**5 + a*x + 1, x)` is `[]` and
`solve((x - 1)*(x**5 + a*x + 1), x)` is `[1]`; `incomplete=False` still
raises `NotImplementedError`, since no root was proven absent, and a
polynomial over `QQ` or the Gaussian rationals still gives `CRootOf`. Like
the reference, a perfect-power rational base `n**k` or `1/n**k` of a
generator is rewritten as a power of `n` before solving, unless the
equation is `b**g = b**h` with the same base (`4**x - 4` gives `[1]`), so
`8**x - 4`, `9**x - 3`, `4**x - 2**x` and `(1/4)**x - 4` give the
reference solutions. The roots of a power of `exp` with at most two roots,
and of a power of a rational base, have their logarithms expanded like the
reference simplify (`exp(2*x) - 2` gives `log(2)/2`). When every
generator is a power `b**g` of a base other than `E`, the equation is
solved like the reference `_tsolve` fallback (the reference inversion of
such a generator always fails on its Dummy right hand side): a single power
term `c*b**g + d` with `-d/c = b**h` gives only `g = h` (`2**(2*x) - 2`
gives `[1/2]`), and a sum is factored after `powdenest` and solved factor by
factor (`4**(2*x) - 3*4**x + 2` gives `[0, 1/2, I*pi/log(2)]`, the three
reference solutions; the fourth value `1/2 + I*pi/log(2)` is not returned
by the reference). Floats are converted like the reference `solve`: unless
`rational=False`, each equation is converted with `nsimplify(f,
rational=True)` and, when `rational` is not given, the solutions are mapped
through `nfloat(exponent=False)`, deduplicated and sorted
(`4.0**x - 8` gives `[1.50000000000000, 1.5 + 4.53236014182719*I]`,
`exp(2.0*x) - 3` gives two solutions). Like the reference `solve`,
hyperbolic functions of the unknowns are rewritten through `exp` before
solving (`cosh(x) - 2` gives `[log(2 - sqrt(3)), log(sqrt(3) + 2)]`). The
generator of `exp(a)` follows `ExpBase.as_numer_denom`: when `a` is
negative, or could extract a minus sign while `-a` is not known negative,
the numerator is taken with `exp(-a)` (`exp(-x) + 1` gives `[I*pi]`,
`exp(-x/2) + 1` gives `[2*I*pi]`). The checksol numerical step evaluates
`val.n(18).n(12, chop=True)` like the reference, so an exactly cancelling
real part is chopped (`log(x) + log(x - 1)` gives `[1/2 + sqrt(5)/2]`), and
`log(u)` is known positive or negative from the sign of `u - 1`, so
solutions are filtered by the sign assumptions of the unknown. Two
dependent terms in powers of different bases are combined like the
reference `_invert` (`powsimp` of their quotient), and a dependent product
of powers through `powsimp(powdenest(...))` (`3**x - 2**x` gives `[0]`,
`exp(x)*2**x - 3` gives `[log(3)/(log(2) + 1)]`). The logarithmic Lambert
form `k*(p*x + q)**n*log(c*x + d) = R` with `p*x + q` proportional to
`c*x + d` gives `(exp(W(n*S)/n) - d)/c` with `S = R/(k*(p/c)**n)`, over the
branches of the existing Lambert convention (`x**2*log(x) - 1` gives
`[exp(LambertW(2)/2)]`, `x*log(x) + 1/4` gives both real branches), and
`f(x)**g(x) = c` is solved like the reference `_tsolve`: the solutions of
`g*log(f) - log(c)`, plus the checked candidates `f = +-r`, `g = k` for
each exact rational `c = r**k` (`c = 1` gives `g = 0` and `f = +-1`), so
`x**x - 4` gives `[2]` and `x**x - 1` gives `[0, 1]`. Like the reference,
each Lambert solution is replaced by a simpler numeric form when it agrees
with it (to 200 bits here): a rational, `exp(r)` or `sqrt(r)`, as
`nsimplify` (mpmath `identify`) finds them (`x*exp(x) - 2*exp(2)` gives
`[2]`, `x*log(x) - 3*exp(3)` gives `[exp(3)]`). The real-branch domain of
`LambertW` is decided for arguments with logarithms of rationals by an
outward rational enclosure of `log`, so `x*log(x) + log(4)` keeps the
complex principal solution `exp(LambertW(-log(4)))`. Remaining
differences:
- symox keeps the `LambertW(z, -1)` branch, real for `-1/e < z < 0`, in
  equations where the reference loses it through its `bivariate_type`
  path with a symbolic `_u`, which it takes when the direct Lambert attempt
  fails (an `Add` such as `exp(x) - x - 2`, where `log(exp(x))` does not
  expand for complex `x`): `exp(x) - x - 2` gives
  `[-2 - LambertW(-exp(-2)), -2 - LambertW(-exp(-2), -1)]` where the
  reference has only the first solution (also `(x + 2)*exp(-x) - 1`,
  `exp(x) - 3*x`, `x*exp(-x) - 1/4`, `exp(x)/x - 3`, `exp(x) - 2*x - 1`;
  both values are solutions, so the reference set is incomplete);
- the same `W_-1` branch is kept for a power of a non-`E` base, where the
  reference returns only the principal branch: `2**x - 4*x` gives
  `[4, -LambertW(-log(2)/4)/log(2)]` where the reference has only the second
  value (`2**4 = 16 = 4*4`), `3**x - 3*x` adds `1` and `2**x - x - 2` adds
  `2`; `x*2**x + 1/4` and `x**2 - 2**x` match the reference. The Lambert
  argument of an exponential factor `b**u` is `u*log(b)`, and
  `k*(a*x + b)**n + C*exp(c*x + d)` is rewritten as
  `k*(a*x + b)**n*exp(-c*x - d) + C` before the Lambert attempt;
- an `exp(c + k*x)` generator with a symbolic constant `c` is inverted like
  the reference split `exp(c)*exp(k*x)`, so `exp(a - x) + 1` gives
  `[log(-exp(a))]`. With a rational exponent `k` other than `+-1` the forms
  still differ: `exp(2*x + a) - 1` gives `[log(exp(-a))/2,
  log(-exp(-a/2))]` where the reference has `[log(-sqrt(exp(-a))),
  log(exp(-a))/2]` (equal values);
- `simplify` applies the reference step `shorter(expand_log(e, deep=True),
  logcombine(e))` after `futrig` when the expression has a logarithm, so the
  Lambert log forms `x**(2*x) - 4`, `x*log(x) - log(9)/2`,
  `x*log(x) + 2*log(2)` and `x**x - sqrt(2)` match. The step is not
  idempotent and the remaining power-base forms still differ:
  `4**x - 8` gives `[3/2, 3/2 + I*pi/log(2)]` where the reference prints
  `(log(8)/2 + I*pi)/log(2)` (also `(1/4)**x - 8`, `4**(-x) - 8`,
  `4**(2*x) - 8`, `4**x - 16`, `4**x + 8`, `8**x - 2`, `9**x - 27`,
  `2**(2*x) - 8`, `2**(3*x) - 8`, `(1/4)**x - 16`, equal values in the same
  order). The divergent step is the outer `simplify` of each solution:
  `simplify(3/2 + I*pi/log(2))` is `(log(8)/2 + I*pi)/log(2)` in the
  reference, reached by its sequence (`normal` bottom up gives
  `(3*log(2) + 2*I*pi)/(2*log(2))`, `Mul(*powsimp(e).as_content_primitive())`
  keeps it, then `factor_terms` and `shorter(expand_log, logcombine)`). The
  native size search keeps the input. Following the ported
  `reference_sequence` for numeric logarithmic expressions fixed
  `4**x - 16`, `4**x + 8` and `(1/4)**x - 16` without regressions in a
  56-call probe, but not the others: the native `powsimp` distributes the
  rational content into the sum (`(3*log(2)/2 + I*pi)/log(2)`, where the
  reference `powsimp` keeps the input), and the sequence then returns
  `3/2 + I*pi/log(2)`. Flattening the content product did not change that,
  so a later step also diverges; the change was not kept. `x**(3*x) - 8` gives `exp(LambertW(log(8)/3))` where the
  reference, whose polynomial generator is `x**x`, has
  `exp(LambertW(log(2)))` (equal value);
- the branch order and forms of a cubic in a power of `2` or `exp` differ:
  `2**(3*x) - 3*2**x + 1` and `exp(3*x) - 3*exp(x) + 1` list the roots in
  another order, and the native `log` takes `I*pi` out of a provably
  negative argument, so one root prints `log(A) + I*pi` where the reference
  keeps `log(-A)` (equal values);
- a rational base that is not a perfect power of an integer or of its
  reciprocal keeps another form: `(4/9)**x - 8` gives
  `log(2**(3/(2*log(2/3))))` where the reference has
  `log(2**(3/log(4/9)))` (equal value);
- `solve(x**3 - 3*x + 1, x, simplify=True)` differs in form: the native
  simplify search does not reach the reference forms of the Cardano roots
  (equal values);
- the native Pow does not fold `sqrt(2**(1/3)/2)` to `2**(2/3)/2`, so
  `exp(-3*x) - 2` gives `log(sqrt(2)*sqrt(2**(1/3))/2)` (equal value);
- for `2**(3*x) - 3*2**x + 1` the reference prints one solution with the
  logarithm of an exactly negative number in a form that its evalf puts at
  `-I*pi`; the native value is the principal one, `+I*pi`, so the two
  solutions differ by `2*I*pi/log(2)` (both are solutions; the other
  forms differ by the simplify gap above);
- a single logarithm of a rational function, `c - log(x) + log((y - 1)/(y +
  1))/2`, gives `(-x**2 - exp(2*c))/(x**2 - exp(2*c))` where the reference
  rewrites it to `1/tanh(c - log(x))` (equal value; this form predates the
  generator inversion above, which applies the substituted inverse only to
  `exp` and `b**a` generators and inverts other generators, such as
  logarithms, root by root as before).

Laplace and Mellin transforms fall back to the reference integration
backends (`_laplace_transform`, `_mellin_transform`) when no table rule
applies. The Fourier, sine and cosine transforms, the Hankel transform and
all their inverses are computed like the reference: first through the
integral (`_fourier_transform`, `_sine_cosine_transform`,
`_hankel_transform`), whose first Piecewise condition is the returned
convergence condition (`Eq(Abs(arg(w)), 0)`, `Abs(arg(a)) < pi/2`, `k > 0`,
...). The native closed form tables are used only when the native integral
is not computed. They read the half plane, the strip and the auxiliary
conditions from the conditions of the Meijer G integral, and they keep the
reference error behavior: `needeval=True` raises `IntegralTransformError`,
the Hankel transforms raise it instead of holding the result, and the
`noconds`/`simplify` defaults and the returned tuples are those of the
reference. The transform stays unevaluated when the integral does not
converge. A Fourier transform also stays unevaluated when the function is not
integrable on the real axis. A Mellin transform with no admissible strip
stays unevaluated as well. Examples that now match the reference:
`laplace_transform(cos(sqrt(t)), t, s)` (an equivalent form),
`exp(-a*t**2)` with `Abs(arg(a)) <= pi/2`, `erfc(sqrt(t))`, the Fourier
transform of `exp(-x**4)`, Hankel transforms of `exp(-a*x)/x`, and Mellin
transforms of `exp(-x)*log(x)`, `cos(x)*exp(-x)` and `log(x)/(1+x)`. The
results and conditions of the reference transform tests (sine and cosine
transforms of `t**(-a)`, `exp(-a*t)`, `1/(a**2+t**2)`,
`exp(-a*sqrt(t))*cos(a*sqrt(t))`, the Meijer G inverses, the Fourier
transforms of Gaussians, boxes and one sided exponentials, the Hankel
transforms of `1/r**m` and `r**nu*exp(-a*r)`) are matched;
`tests/test_transforms_remaining_parity.py` holds the reference outputs.
The `simplify=True` result step is `simplify(powdenest(piecewise_fold(e),
polar=True))`. Because the native `simplify` is lighter, the steps of the
reference `simplify` that matter for these results are applied as well: the
`expand_log`/`logcombine` step, the `gammasimp` step for gamma functions of
non-integer argument and the final `exptrigsimp` step (`cosh(a*t) -
sinh(a*t)` is `exp(-a*t)`).

Remaining transform differences:
- Equal values printed differently, because the native `simplify` is not
  the reference one: `sine_transform(log(t)/t)` (`-log(w**2)/4 -
  EulerGamma/2` inside the product), sums such as `exp(-t) + 1/sqrt(t)` or
  `exp(-t) + f(t)` (not put over one denominator),
  `cosine_transform(1/(a+t))` (the sign of the `Si` factor),
  `fourier_transform(x*exp(-a*x)*Heaviside(x))` (`1/(a**2*(1 +
  2*I*pi*k/a)**2)`), the indicator of `[0, 1]` (one fraction), the `2**(1 -
  m)` factor of `hankel_transform(1/r**m, r, k, nu)` and the arrangement of
  `laplace_transform(cos(sqrt(t)))`.
- A Piecewise result of a Fourier, sine, cosine or Hankel transform (and
  their inverses) goes through the reference `piecewise_simplify`, the
  `simplify` of a Piecewise: univariate conditions become real intervals,
  `Abs(w**2) < 1` is `(w > -1) & (w < 1)` like the reference, which
  assumes a real `w` there (sine and cosine transforms of `besselj(0, t)`
  and `sin(t)/t`, the inverse Fourier transform of `sin(k)/k`,
  `tests/test_transforms_wave3_parity.py`). The nested Piecewise of the
  cosine transform of `sin(t)/t` is still folded (its outer conditions are
  the reference ones).
- `hankel_transform(1, r, k, 0)` (and the inverse) raises
  `IntegralTransformError`: the integral of `r*besselj(0, k*r)` over
  `(0, oo)` diverges, and the reference returns a meaningless
  `Piecewise` of `nan`, `oo` and `besselj(1, zoo*k)` from the limits of the
  antiderivative at infinity; symox does not compute this integral.
- The native definite integral of `exp(-a*x)*sin(b*x)*Heaviside(x)` over the
  real line with positive `a`, `b` and a real `k` is not computed. The table
  gives the reference value with the condition `True` instead of the
  reference `Abs(arg(1 - 2*pi*k*exp_polar(-I*pi/2)/a)) < pi/2` (the
  reference `meijerint_definite` result). Without assumptions the native
  integral is computed with `(re(a) + im(b) > 0) & (re(a) - im(b) > 0)`
  where the reference holds the transform. The Laplace transform of the
  same function matches the reference (`-re(a) + Abs(im(b))`, `True`).
- `fourier_transform(exp(-a*x**2), x, k)` with an unrestricted `a` is held
  like the reference, with and without `noconds`: when the transform
  integral is computed only under conditions and the integral of the
  function is not computed, the reference error is final and the closed
  form table is not tried (it returned the value, invalid for
  `re(a) <= 0`, with no condition under `noconds=True`).
- `hankel_transform(1, r, k, 0)` diverges: native raises
  `IntegralTransformError`, where the reference returns a Piecewise with
  `nan` and `oo` pieces.
- Laplace table rules return `True` conditions with `Abs(k)` values, where
  the reference integral gives `Eq(Abs(arg(k)), 0)` or `k > 0`.
- Exponential rate conditions of the Laplace rules are `re(a) > 0`, where
  the reference has `Abs(arg(a)) < pi/2`.
- Some Laplace rule results differ where SymPy 1.14 is wrong:
  `1/(t+1)`, `exp(-t)/(t+1)` and `Heaviside(t-1)*sqrt(t)`.

The inverse Laplace transform now applies the reference rules for the
irrational functions of Bateman 5.3 (`_inverse_laplace_irrational`): once
after the simple rules and again as the last program rule. The expansion
step is the reference sequence `expand(deep=False)`, `expand_mul`, `expand`,
then `apart`. The cosine part of a quadratic denominator uses the first root
of `x**2 - b`, as in the reference. Factors match `(a*s**m + b)**n` through
`as_base_exp`. Zeros and poles are sorted by `(n, b != 0, b)`. An undecidable
order leaves the function to the later stages, where the reference raises
TypeError. Where the reference closed form fails `L[f](s) == F(s)` (checked
by quadrature), the native rule returns the correct value:
- 5.3.6, `sqrt(s)/((s-b)*(sqrt(s)+a))`: the middle term has `erf`, not `erfc`;
- 5.3.10, `1/(sqrt(s)+a)**2`: the value is
  `(2*a**2*t+1)*exp(a**2*t)*erfc(a*sqrt(t)) - 2*a*sqrt(t)/sqrt(pi)`;
- 5.3.16, `1/(sqrt(s)+a)**4`: the factor `a**3*t**(5/2)` is `a*t**(3/2)`;
- 5.3.14 and 5.3.15: the shift is divided by the coefficient of the cubed
  pole and the scale multiplies the whole value. The reference is right only
  for a unit coefficient;
- 5.3.14 requires the zero `sqrt(a*s)` (not `sqrt(a*s**2)`), and 5.3.2
  requires the pole `sqrt(a*s)` (not `sqrt(a*s+c)`). The reference returns
  wrong values for those shapes, so they stay unevaluated.

Not yet ported:
- the integration stage of the inverse Laplace transform
  (`_inverse_laplace_transform_integration` through the inverse Mellin
  transform and `meijerint_inversion`). Inputs that only it handles stay
  unevaluated, for example `a**b*(s+sqrt(s**2-a**2))**(-b)/sqrt(s**2-a**2)`;
- the `cos` patterns of the Mellin `process_conds`;
- the `hyperexpand` step of the reference `simplify`.

Scalar containment queries now visit explicit arithmetic arguments, including
coefficients and additive constants. They exclude the exponent one stored
implicitly in compact products, and accept `None` as a nonmatching scalar
query. Other structural containment patterns still require separate audits.

Mutable dense matrix subclasses now reconstruct real Python subclass instances
after the qualified arithmetic, slicing, transpose, substitution and linear
algebra operations. The `_new` protocol bypasses a custom public constructor.
Result-class precedence preserves the reference distinction between addition
and multiplication on equal-priority subclasses. Explicit mutability
conversions keep their target category. This is scoped protocol work: custom
Symbol/Function inheritance, immutable/sparse subclasses, subclass serialization
and unqualified matrix methods remain open.

The lightweight trigonometric engine used by native general simplification
memoizes shared subexpressions. Its rule rewrites, rule fixpoints and
trigsimp results are shared by every pass, nested call and candidate of one
`simplify` call. Keys compare structurally, and a hit is used only when the
stored key has the same presentation. Node sizes and head sets are cached by
node identity through weak references, so rejected rewrites of a large sum
are freed. An accepted rewrite of a sum takes back the original terms that
the pattern bindings rebuilt. Head discovery visits each distinct node once
and combines the same head bits; rewrite rules, acceptance conditions and
candidate ranking are retained. Local caches keep expression owners alive
to prevent reuse of identity keys. Exact identity and unchanged-expression
tests cover graphs whose expanded trees have more than `2**64` occurrences.
Public traversal multiplicities are unchanged.

The commutative matcher behind these rules (`pattern.rs`) indexes the terms
of a large sum (at least 8 terms) by value and by the factors left once a
pattern term's fixed factors are divided out. A pattern term whose value is
already determined by the bindings, or a product whose bare wildcard factor
is already bound, is then tried only on the terms that can match. A bound
lone wildcard whose value is no longer among the free terms fails at once.
The scan order, and so the first match returned, is unchanged. Patterns with
wildcard properties are not filtered, since their predicates may report
errors. With `SYMOX_CHECK_DETERMINED_FILTER=1` every skipped term is checked
to fail. Before this, a pattern such as `2*a*cos(b)**2 - a + c`, whose `-a`
matches every term, matched in quadratic time with a linear scan per step.
On a sum of 12,000 terms this made one `simplify` call of the `dsolve` case
in the ODE section run for more than 40 minutes.

These changes remove repeated work, not the inherent cost of expression
growth or pattern matching. The rule engine is not SymPy's Fu-based
`futrig`, and the separate public matching/Fu driver and the complexity of
the complete `simplify` pipeline remain separate obligations.

Single-expression `cancel` follows the reference preprocessing (`signsimp`,
then `factor_terms(radical=True)`), so `cancel(sqrt(4*R1*T*k + 4*R2*T*k))`
with positive symbols is `2*sqrt(T)*sqrt(k)*sqrt(R1 + R2)`. Native
`simplify` also accepts the reference state after its unconditional
`factor_terms(expr, sign=False)` step when the reference's later `cancel`,
`factor_terms` and `expand_power_exp(expand_mul)` comparisons keep it (so
the radical form above is kept while `n**3 - n` stays expanded), and it
applies the reference's unconditional `exptrigsimp` when the result has
trigonometric, hyperbolic or exponential functions
(`simplify(im(exp(2*I)/(1 - exp(4*I))))` is `1/(2*sin(2))`). The remaining
steps of the reference pipeline are still approximated by the native
candidate search.

Expression-domain `DomainMatrix` Gauss-Jordan reduction now uses the native
exact pivot simplifier. This exposes zero pivots such as
`a - (-a)**2/a` before raw-domain elimination divides by them. Public power
constructors retain their established forms, including interpolation formulas
whose factors must cancel at their own nodes. Both EX and EXRAW retain their
coefficient domain. The upstream Thevenin circuit regression reaches the
correct result. SymPy 1.14 EXRAW can produce an invalid inverse for the
regression matrix; the native result is checked against EX and the exact
identity matrix product. This is an intentional mathematical correction.
Matrix numeric zero queries use exact facts and symbolic pivot simplification,
instead of accepting a fixed machine-precision residual as nonzero. Closed
finite entries proved zero are canonicalized before default domain dispatch
and in default RREF output. Explicit host predicates retain their callback
path. Closed function values are not assumed to be independent polynomial
generators during default inversion: trigonometric relations can make a matrix
singular. Exact rank, RREF and inversion checks include scaled identities and
both scalar and two-dimensional singular matrices.

Scalar numerical assumption refinement now checks real/imaginary components
at 160 and 320 bits, with relative stability and scale checks. Unresolved
cancellation retains unknown; a rounded numerical zero is never itself a zero
proof. The strict exact-query policy remains independent. Multiprecision
inference shares precision-specific descendant values across each query, and
held sums may collect identical finite terms with exact coefficients.
Floating coefficients retain their established arithmetic evaluation path.
Tiny sine evaluation now preserves relative precision instead of truncating
the argument to an absolute fixed-point zero. Comparisons cover values as
small as `sin(10**-1000)` at 60 requested decimal digits. General interval
certification and all special-function accuracy contracts remain open.

The rejection of a proved singular numeric matrix by default/DM and LU
inversion is a separately recorded mathematical correction. SymPy 1.14's
formal-domain paths return an inverse for some of the independently proved
zero identities in these regressions. This is not a reference-output match.
General functional relations, parameter strata and domain inference are still
blocking obligations. Typed polynomial rings and fraction fields over Gaussian
integer/rational ground domains remain unsupported, even though the covered
scalar and matrix Gaussian computations are available.

Automatic polynomial division now retries inferred closed algebraic
coefficients with the existing exact number-field engine when the legacy
coefficient view cannot represent the quotient. The original division and
explicit domain paths remain available. Exact quotient/remainder identities
cover radical sums and complex algebraic coefficients. Transcendental fields
are not approximated by this fallback. Default inferred output-domain metadata
can still differ from the reference EX domain; broader coefficient-domain
inference remains an explicit gap. All four synthesis forms for the covered
Lcapy radical-impedance circuit now terminate with the original impedance.

Global `div`, `quo`, `rem` and `exquo` now respect the default disabling
of automatic field extension when a coefficient domain or field choice is
explicit. `auto=True` remains an explicit override, while `None` follows the
reference's absent-option preprocessing. Differential checks cover integer,
rational, polynomial-ring, finite and algebraic fields, including result-domain
metadata and exact-quotient errors. The previously incorrect default-extension
expectation in the ring regression is now requested explicitly with `auto=True`;
separate default-option tests preserve the reference ring behavior.

The LaTeX symbol printer now exposes the real parent `_print_Symbol` hook for
subclass delegation, preserving its `style` argument, exact symbol-name overrides
and super/subscript conventions. Token translation calls the existing Rust
translator directly instead of reparsing each component as a complete symbol.
Output modes wrap the final result once, retaining inline fraction defaults.
This qualifies the tested symbol hooks and printer settings; arbitrary recursive
printer overrides, matrix styles and every layout option remain separate gaps.

Implicit external conversion never invokes a native object's reference exporter.
Conversion chains are traversed iteratively, with retained owners and identity
checks for cycles; reentrant hooks are guarded separately. Matrix `nsimplify`
applies scalar recognition entrywise. These changes do not certify all
recognition options or arbitrary external object protocols. `exp` on a concrete
matrix dispatches to the native matrix function and returns immutable output;
general floating-point matrix eigenvalue support remains incomplete.

Numeric real skew-symmetric matrices of size two and three now use Rodrigues'
polynomial identity. Floating entries are lifted to their exact dyadic values
before forming the squared angle, then evaluated by the native multiprecision
engine and rounded to the input precision. The half-angle quadratic coefficient
avoids subtracting nearly equal cosines. Working precision accounts for large
and tiny inputs; rounded results must agree at two successive working
precisions. Resource exhaustion raises `NotImplementedError`. This numerical
stability check is not an interval accuracy certificate. Other matrices retain
the existing exact spectral path. General floating matrix functions and global
scalar trigonometric accuracy remain open obligations.

Public builtin function heads now retain callable identity, symbolic spelling
and serialization identity for structural replacement queries. Relational
aliases share their constructor class and signature. Neither change supplies
the complete nominal `FunctionClass` hierarchy required by some downstream
subclass declarations. Higher-priority external operands keep control of their
arithmetic before scalar conversion, preserving wrapper metadata such as units.
Finite permutation and scalar protocol tests are scoped differential checks;
the downstream report records the upstream suites that actually reproduce
their passing SymPy 1.14 reference cases.

Undefined function heads support instance queries that preserve their name,
assumptions and arity; held derivatives expose the public constructor for
reconstruction. Callable slots belong to native lambda, predicate and registered
differential-field subtypes, so constants and ordinary expressions are not
reported as callable, including after importing differential geometry.
For registered charts, symbolic atom identities preserve point evaluation
after reconstruction and pickle round trips. The full nominal function and
symbolic subclass hierarchy remain uncertified.
Ordinary undefined applications now use native expression subclasses specific
to their head, assumptions and arity, with an actual `UndefinedFunction`
metaclass ancestry. Reconstruction, traversal and pickle restore these classes;
application patterns remain distinct from constructor queries. The public
`Function('f')` factory remains a separate object, so full function-class
identity, `type(expr) is expr.func` and arbitrary symbolic subclasses remain
uncertified.
Tuple reconstruction accepts the reference keyword signature; raw
`Tuple(sympify=False)` items outside the native expression representation
remain unsupported.

Python `float` and `complex` conversion now prefer the existing adaptive native
precision evaluator before rounding to binary64. Float conversion rejects
genuine small imaginary components instead of treating them as rounding noise.
Exact numeric literals retain their direct single-rounding conversion, including
subnormal ties. Numerical evaluation memoizes shared physical nodes by working
precision and refinement context; repeated DAG branches are not expanded into
repeated evaluations at the same precision.
Functions outside that evaluator retain the existing binary64 numerical path;
this is not a certificate of all special-function precision behavior.
The default compiled numerical callable preserves its float64 numeric path and
uses generated native NumPy source for foreign operand protocols and complex
arrays. Array subclasses can retain units or masks through their own arithmetic.
Unsupported scalar conversions return `NotImplemented` from arithmetic so
foreign reflected operators can run; errors raised by conversion hooks still
propagate to callers.
Fraction and Decimal inputs still use the existing float64 numerical path;
preserving their exact operand arithmetic remains an open obligation.
Arbitrary module overrides and all numerical-callable domains remain separate
parity obligations.

The NumPy numerical-callable path also flattens nested formal arguments and
restores nested list/tuple and concrete matrix results, including empty shapes.
Applied functions and derivatives used as formal parameters are replaced
simultaneously by independent native dummy arguments before compilation.
Other module configurations retain their existing paths; arbitrary signatures,
all vectorized mixed shapes and every lambdify option remain separate obligations.

Downstream mechanics now exposes Kane's stored independent state vectors and
honors explicit inverse-method choices. Dependent and auxiliary Kane coordinates
remain explicitly unsupported. Matrix assignments lower to row-major scalar
assignments through the existing printer hooks. Canonical CSE ordering compares
expanded mathematical node counts and sort keys, with arbitrary-size counters
and no former 64-argument ordering cutoff; `order='none'` keeps its unsorted
traversal. These changes do not certify all printing or CSE options.

Symbolic Cholesky/LDL solves preserve compact exact factorization expressions
instead of expanding every intermediate fraction. A defined rational
specialization of full column rank can certify generic full rank for a
commutative rational-function matrix. Failed probes prove nothing and retain
the original exact rank computation, including its singular-matrix errors;
non-rational-function inputs also retain that rank path. Pivot and nonfinite
checks and the normal-equation fallback remain;
factor inspection visits shared physical nodes once. Rational matrices retain
the existing rank and cancellation path. Regression checks specialize dense
symbolic systems to exact rational and complex numerical inputs. These checks
and the PyDy pendulum case establish specific correctness and expression-growth
improvements, not comparative bit-complexity bounds for all matrix solvers.

Polynomial `gcd` and `lcm` accept native symbolic tuples and finite sets as
sequence inputs, preserving coefficient conversion errors. Nonempty reductions
start with the first input, so single negative terms and rational LCMs retain
the reference normalization. Sequence options and all coefficient domains still
require their own behavioral audit. Symbolic `Dict` values retain their mapping
identity during sympification and expose atom traversal. `nsimplify` takes the
exact rationalization path for expressions with free symbols, including mapping
keys and values, before attempting numerical recognition. General container
protocols and full irrational recognition remain separate obligations.

Permutation groups now expose `generate`, `generate_dimino`, `random`,
`conjugacy_class` and `conjugacy_classes`. Dimino generation builds incremental
subgroup cosets, retaining output-sized storage. A conjugacy class is an orbit
under generator conjugation and does not enumerate unrelated group elements.
Uniform sampling chooses one representative independently from each exact
stabilizer transversal; the unique transversal product is a bijection with the
group. Sampling therefore does not enumerate a factorial-sized group. Tests
compare generated sets, Dimino order for identical generators and conjugacy
partitions with the reference, and exhaust all transversal choices on a small
group. They do not certify all permutation-group methods or a comparative bound
for constructing every stabilizer chain.

- Accumulation bounds follow `sympy/calculus/accumulationbounds.py`:
  `AccumBounds` is a port of the reference class (arithmetic, powers,
  comparisons, `Min`/`Max`, `intersection`, `union`, `in`), and its whole
  test file passes (`tests/test_accumbounds_wave4_parity.py`). The `eval`
  branches that return bounds are reproduced: `sin`, `cos` at `oo`/`-oo`
  (and `cos` of a real infinity such as `oo*sign(s)` with a nonzero `s`) are
  `AccumBounds(-1, 1)`, `tan`, `cot`, `sec`, `csc` at `oo`/`-oo` are
  `AccumBounds(-oo, oo)`, `frac(oo)` is `AccumBounds(0, 1)`, `atanh(zoo)`
  and `asech(zoo)` are `I*AccumBounds(-pi/2, pi/2)`, `atan(zoo)` is
  `AccumBounds(-pi/2, pi/2)`, `sinh`, `cosh`, `tanh`, `coth`, `sech`,
  `csch` of `I*oo` are the trigonometric bounds (`sinh(I*oo)` is
  `I*AccumBounds(-1, 1)`), `exp(I*oo)` is `nan`,
  and `sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `exp`, `log`, `sqrt` of
  bounds are bounds. `subs` of an infinity or of bounds re-evaluates these
  applications as the reference rebuild does (`sin(x).subs(x, oo)`,
  `(x**2 + 2*x + 1).subs(x, AccumBounds(-1, 1))`).
  A limit the native engine leaves unevaluated and that holds an
  oscillating function (`sin`, `cos`, `tan`, `cot`, `sec`, `csc`, `frac`) of
  an argument tending to an infinity is evaluated term by term with the
  bounds arithmetic, as the reference leading terms do
  (`limit(sin(x), x, oo)`, `limit(1/(2 + sin(x)), x, oo)`,
  `limit(sin(1/x), x, 0)`), and a definite integral with an infinite
  endpoint subtracts such limits of the antiderivative, as
  `ExprWithLimits._eval_interval` does (also for `Limit(...).doit()` and
  `Integral(...).doit()`). A limit the native engine already
  answers is never replaced, and any `nan`, `oo - oo` or unevaluated part
  keeps the limit or the integral unevaluated.
  Intentional divergences, where the reference value is degenerate or an
  error: `limit(x*sin(x), x, oo)`, `limit(x**2*cos(x), x, oo)`,
  `limit(exp(x)*sin(x), x, oo)` and `integrate(x*sin(x), (x, 0, oo))` are
  `AccumBounds(-oo, oo)` (the reference returns
  `oo*sign(AccumBounds(-1, 1))`), `limit(sin(x)/(1 + exp(-x)), x, oo)` is
  `AccumBounds(-1, 1)` (the reference raises `TypeError`), and
  `limit(sqrt(x + sin(x)) - sqrt(x), x, oo)` stays `0` (the reference leaves
  it unevaluated).
  Remaining gaps: `sin(oo)` built inside the kernel (a native
  substitution or simplification that does not go through the public
  constructors or `subs`) stays unevaluated, and so does a held `Limit` or
  `Integral` nested in a larger expression on `doit()`
  (`(x + Integral(sin(x), (x, 0, oo))).doit()`); `cos(oo*sign(s))` is evaluated
  but `sin(oo*sign(s))` stays, as in the reference; `limit(sinc(x), x, oo)`
  and `limit((-1)**n, n, oo)` stay unevaluated (the reference gives `0` and
  `nan`); a limit at `zoo` and a limit of an `Integral` stay unevaluated
  where the reference raises `NotImplementedError`.

- The integration engine does not implement the full collection of Risch,
  heurisch and Meijer G methods available in SymPy. Meijer G definite
  integration, its convergence conditions and `hyperexpand` are ported (see
  above), and so are the Meijer G indefinite integrals and the Piecewise
  integration of `Integral.doit`, and so are the `risch`, `heurisch` and
  `manual` methods with their options and the `NonElementaryIntegral`
  results of a Risch proof (see the integration methods paragraph for the
  remaining gaps).
  The integration corpus (`tools/integration_corpus.py`) now reports 1083
  identical and 68 equivalent results (1070 and 70 before the parameter
  work below); 31 cases stay unevaluated where the reference answers.
  The rational pole analysis of definite integrals decides the position of
  a pole with parameters from the assumptions (`-a` is outside `(0, 1)` for
  a positive `a`, so `1/(x+a)`, `1/(1+a*x)` and `1/x` on `(a, 1)` evaluate),
  the Gruntz sign of a parameter asks the assumptions as the reference does
  (`x**a` on `(1, oo)` is `oo` for a positive `a`), `log(a*x)` integrates to
  the reference `x*log(a*x) - x` (its value on `(0, 1)` holds for every
  nonzero `a`), and an indefinite integral left inside a Piecewise
  antiderivative becomes a definite integral, as in
  `ExprWithLimits._eval_interval` (`Heaviside(x - a)` on `(0, 1)`).
  The 31 open cases are of two kinds:
  - 14 oscillatory integrals on `(0, oo)` (`sin(x)`, `cos(x)`, `sin(s*x)`,
    `cos(s*x)`) that the reference returns as `AccumBounds`. They now
    evaluate like the reference (`sin(x)` is `AccumBounds(0, 2)`,
    `sin(s*x)` with a positive or nonzero `s` is
    `AccumBounds(-1, 1)/s + 1/s`; see the accumulation bounds paragraph
    below). Without assumptions or with a real or complex `s` both libraries
    keep them held. The one left held is `cos(s*x)` with a nonzero `s`, for
    which the reference returns `sin(oo*sign(s))/s`, an unevaluated
    function of an infinity;
  - 17 values the reference returns without a condition although they are
    wrong for some allowed parameter values, held on purpose: `1/(x+a)`,
    `1/(x-a)` and `1/(1+a*x)` on `(0, 1)` and `1/x` on `(a, 1)` without
    assumptions or with real, complex or nonzero `a` (the pole lies on the
    interval for some `a`: `-log(a)` at `a = -1` is not the value of a
    divergent integral), and `1/(x-a)` with a positive `a` (divergent for
    `0 < a <= 1`). With a positive `a`, the pole `-a` inside `(-a-1, 0)`
    makes `1/(x+a)**2` diverge to `oo` here, where the reference returns
    `-1 - 1/a`.

  The remaining numerical differences are reference errors: `laplace_transform(1/(t+1))`
  is `-exp(s)*Ei(-s)` (the reference gives `-exp(-s)*Ei(-s)`), and the
  Fourier transform of `1/(1+x**2)` keeps `Abs(k)`.
  Gamma moments `c*x**a*exp(delta-beta*x**q)` on `(0, infinity)` now retain the
  exact `re(a)>-1` condition for provably positive real `q`, or `re(a)<-1`
  for provably negative real `q`, with `re(beta)>0`. This includes rational
  and symbolic signed powers, complex parameters/rates and a free phase offset.
  `conds='piecewise'`, `'none'` and `'separate'` select actual result forms.
  A failed convergence test does not authorize analytic continuation as an
  integral value. Sums retain one original-integral fallback, preserving
  cancellation cases. The moment reduction keeps expression-valued exponents;
  even `2**100` does not require a degree-sized array or an i64 degree counter.
  Complex/sign-unknown powers, other convergence regions and the full strategy
  families remain separate obligations.
  The [scoped integration review](../audits/sympy-1.14.0/integration-review.json)
  records the Mellin/contour argument and sparse phase handling.
- `evalf(strict=True)` raises `PrecisionExhausted` with the reference message
  when a cancelling sum cannot be distinguished from zero. `maxn` bounds the
  extra working precision of each sum by `max(prec, int(maxn*log2(10)))` bits
  with SymPy's restart schedule (the extra precision is bounded by
  `min(maxprec, 2*working)` before each restart, and the last restart may
  overshoot by the missing accuracy, at most `2**20` bits in total). Without
  `strict` an explicit `maxn` returns the sum at the accuracy it reached, as a
  Float of lower precision. Without `strict`, a sum evaluated as a whole
  that exhausts its budget is returned as SymPy's `evalf_add` returns it:
  each part has precision `max(min(prec, accuracy), 1)`, and a part without a
  correct bit (an exact zero such as
  `-I*exp(3*I/5)/2 - exp(3*I/10)*sin(3/10) + I/2`) is the precision 1 Float
  equal to its absolute error bound (SymPy's scaled zero), printed without
  digits as `0.e-168 + 0.e-168*I`. The schedule starts from SymPy's `prec+4`
  target and compares the terms' working precision without their ten guard
  bits, so the reference example gives SymPy's exponents at 30 digits.
  `chop=True` still returns `0` and `strict=True` still raises. evalf without
  options first tries the native 4096-bit budget, which resolves
  cancellations that SymPy's default `maxn=100` reports as a bound, and
  reports an unresolved sum with SymPy's default `maxn=100` schedule.
  Remaining divergences: the bound is always the positive error bound, where
  SymPy returns its rounding residue when it is not exactly zero (sign and
  magnitude of noise, such as `-0.e-125` for `sin(pi/7)**2+cos(pi/7)**2-1`);
  the accuracy estimate is a rigorous rounding bound rather than SymPy's
  `add_terms` estimate, so exponents and achieved precisions can differ by a
  few bits, and for some sums by up to about 30 decimal digits of the bound;
  an exhausted sum nested in a numerical expression (`exp(z)`, `1 + I*z`,
  `3*z` with an unevaluated product) stays unevaluated, because native
  evaluation does not propagate SymPy's per-value accuracy (only a symbolic
  factor, as in `x*z`, is evaluated around it);
  `tests/test_evalf_cancel_parity.py` records the agreed cases.
  Branch cuts audit (same file, 123 expressions at 15, 30 and 100 digits,
  369 points: `log`, `sqrt`, symbolic powers, inverse trigonometric and
  hyperbolic functions, `LambertW` branches, `polylog`, `zeta` near 1,
  `hyper`/`meijerg` at the unit argument, `loggamma`, Bessel functions of
  negative arguments, on exact cuts and 10**-30 above and below them): every
  point has SymPy's branch and digits correct to the requested precision,
  357 part by part and 9 relative to the larger part where SymPy itself is
  only accurate in that sense (`hyper([1, 1], [2], 2)` and its neighbours).
  Native evaluation is more accurate than SymPy on tiny parts near a cut
  (`sqrt(-2 + I/10**30)`, `acosh(1/2 + I/10**30)`) and evaluates
  `zeta(1 + 1/10**20)` and `hyper([1, 1], [2], 1 - 1/10**20)` at 15 digits
  where SymPy rounds the argument onto the pole and fails. Remaining gap:
  `meijerg([[1, 1], []], [[1], [0]], -1)` stays unevaluated where SymPy
  returns `-oo` (the logarithmic singularity).
  `verbose` prints nothing; `quad` is ignored without an Integral (as in
  SymPy) and refused with `NotImplementedError` for an Integral. The former
  10,000-decimal-digit cap has been removed; native precision still uses a bounded machine counter with
  space reserved for working precision. Binary exponents use an `i128`
  counter bounded by `2**100`; out-of-range arithmetic must retain an exact expression
  rather than wrap to an unrelated number. Explicit rational output retains
  allocation budgets even though Float literals have compact storage at much
  larger exponents. Removing policy caps does not certify all numerical
  algorithms, counter ranges or precision options.
- Complex algebraic roots now have the reference certified isolating
  rectangles, exact complex ordering and rounded refinement; see the
  2026-10-21 section above for the remaining algebraic-domain gaps.
- Special functions that SymPy evaluates through mpmath (Bessel, Airy,
  gamma family, error and exponential integrals, zeta family, elliptic
  integrals, Lambert W, `hyper`, `appellf1`, `meijerg`, orthogonal
  polynomials and the non-integer Bernoulli/harmonic/Euler extensions)
  now evaluate at any precision through the same mpmath functions, with
  certified agreement between two guard precisions (up to `2*prec + 512`
  bits); see [the special-function guide](special-functions.md#multiprecision-evaluation).
  `meijerg([[],[]],[[1/3],[]],2/3).evalf(50)` gives 50 correct digits, and
  `exp_polar` branches follow SymPy's `meijerg._eval_evalf`. mpmath errors
  (`erfinv(2.7)`, `zeta(1.0)`, `loggamma(-2.0)`) propagate as `ValueError`;
  assumption queries do not raise. Float arguments auto-evaluate under
  SymPy's `_should_evalf` rule at the largest argument precision, and negative
  Float bases with fractional exponents follow `Float._eval_power`
  (`(-1.3)**(1/3)` keeps a `(-1)**(1/3)` factor; `sqrt(-1.3)` has no real
  noise). Remaining divergences: digits are certified, so ill-conditioned
  arguments can be more accurate than SymPy (`elliptic_k(1 - 10**-30)`,
  `(-1.3)**(-1/3)` in the last digit) and slowly converging mpmath cases such
  as `hyper([2,3,4],[1],1/2)` take a few times longer; a non-finite mpmath
  result stays unevaluated instead of becoming `oo`/`nan`; the bare Rust
  library without the Python package keeps double precision kernels for these
  families; `marcumq` stays double precision (SymPy has no kernel).
- Elementary functions (exponential, logarithm, trigonometric, hyperbolic
  and their inverses) evaluate each component of a value at the requested
  relative accuracy for an exact argument, including near their zeros, on
  and near branch cuts and under cancellation (`log1p` forms, products such
  as `(1 - x)*(1 + x)`, the Hull, Fairgrieve and Tang complex `asin` and
  `acos`). Branch cut values follow mpmath, through which SymPy evaluates
  them, and the argument is evaluated with as many bits as the condition of
  the function demands. A Float argument is evaluated only after the
  function's exact rules return nothing, as in `Function.__new__`. The
  remaining exponent limit is the bounded `i128` counter above. Analytic
  continuation is that of mpmath for numerical values and of SymPy's
  `exp_polar`/`unpolarify` rules for symbolic ones; there is no audit of
  every multivalued function's continuation outside these paths.
- Ordinary differential equations, ODE systems and first-order PDEs run
  SymPy 1.14's own methods (`ode.py`, `single.py`, `riccati.py`,
  `nonhomogeneous.py`, `hypergeometric.py`, `systems.py`, `pde.py`) on the
  native algorithms, with every hint, option, `ics` handling and the
  Airy/Bessel families; see
  [the differential equation guide](differential-equations.md). The native
  solvers run only when every reference method raises `NotImplementedError`,
  and their bounds (64 series coefficients, 64 functions in a system, eight
  in a matrix exponential) apply only to that fallback. On 315 equations
  drawn from SymPy's ODE tests, classification output matches for 296, the
  test's hint for 262 and the default solution for 258 (string comparison).
  The differences come from equivalent candidates chosen by `simplify`,
  `integrate` and `factor`, from Gaussian polynomial domains (`ZZ_I`/`QQ_I`
  inference) and from integrals that remain unevaluated. The native
  `simplify` is a candidate search followed by SymPy's
  `factor_terms(sign=False)` and `hollow_mul` steps, not SymPy's whole
  sequence of rewrites, so its form can differ. SymPy applies
  `factor_terms` unconditionally; the native step keeps it only when it is
  not larger and keeps a provably real Cardano root provably real. The
  `factorable` hint depends on that form (two corpus classifications).
  SymPy's `cancel(expand_log(expr, deep=True))` step is used only to
  recognize zero (`log(100)` against `2*log(10)`), not to choose a form.
  Orders above `10**6` raise `NotImplementedError` in `classify_ode`
  (SymPy does not return for them). The Lie group
  heuristic `abaco2_unique_general` is still slower than in SymPy on large
  expressions. For `(x + cos(f(x)))*f(x).diff(x) + f(x) + 2*x` it calls
  `simplify` on three expressions. The native candidate search now also
  runs SymPy's own sequence (the rational steps, `factor_terms`, Fu's
  `trigsimp`, the `powsimp`, `cancel` and `factor_terms` shorter steps and
  `exptrigsimp`) as one more candidate, so the factored powers of
  `x + cos(_y)` that SymPy keeps are no longer lost to the multiple-angle
  expansions. The first two calls take 0.3 s and 1.3 s (SymPy: 2.3 s and
  10.1 s) and return 142 and 150 operations, SymPy's sizes (before: 219
  and 176). The third input is then SymPy's (750 operations, before 930)
  and that call returns 641 operations (SymPy: 642) in about 77 s (before:
  210 s), where SymPy takes 71 s. On that input, the reference sequence
  itself takes about 2.5 s; most of the rest goes to the native `trigsimp`
  candidates of the input and of its expansion (a sum of about 12000
  terms), in the native trigonometric rule engine. That engine now interns
  the nodes it visits (equal sub-expressions share one allocation, so the
  memo lookups compare pointers instead of whole sums), skips the
  sub-expressions that lack a function head of a rule before any lookup,
  matches each sum term only against the pattern terms whose heads it
  has, remembers the terms that failed a pattern term from empty bindings
  and no longer re-adds the terms of a sum when the pattern has no fixed
  term (before: 97 s and 62 s for the two candidates). A default `dsolve`
  of this equation takes about 78 s in a release build (before: 210 s,
  and 395 s earlier), where SymPy takes 85 s.
  SymPy then raises `NotImplementedError`, while symox returns its native
  fallback answer `Eq(x**2 + x*f(x) + sin(f(x)), C1)`. Before the
  trigonometric memo and the matcher indexes it took more than 40 minutes.
  The pattern matcher still rebuilds the terms of a sum (`args()`) for each
  match, so that third call stays a little slower than SymPy's; this
  remains an open gap. Deliberate
  extensions are the logarithmic Frobenius basis for repeated or resonant
  indicial roots, the corrected expansion of non-polynomial Frobenius
  coefficients, irrational indicial roots (see the Frobenius section), `ics` that name the generated constants, exact
  characteristic multiplicities in the fallback, the native PDE
  characteristic fallback, and `pdsolve` names that do not capture the
  equation's `xi`, `eta` or `F`. Nonlinear first-order PDEs, Charpit's
  method and second-order PDEs are not implemented, as in SymPy.
- Complex algebraic roots do not yet have certified complex isolating boxes;
  complex root ordering and numerical refinement retain the limitations
  described in [the root guide](algebraic-roots-quaternions.md).
- General transcendental solving, all solver options, symbolic parameter case
  distinctions and full set-valued solution-family coverage still need a
  versioned behavioral audit. Supported Lambert families do not establish parity
  for arbitrary mixtures of powers, exponentials and logarithms.
- Native `Wild(properties=...)` now checks ordered Python predicates at binding
  time, after exclusions, and continues commutative backtracking on rejection.
  Free wildcards can absorb groups or neutral elements; canonical reciprocal
  coefficients are matched without a floating-point conversion. Predicate
  errors propagate through matching, replacement, discovery and collection.
  Owned predicate tuples participate in symbol identity and survive native
  expression and matrix serialization. This covers ChemPy's scalar coefficient
  constraint, not a global matching certificate. General commutative matching
  can require exponential partition search; comparative complexity, repeated
  side-effectful predicates, exotic unhashable callable arithmetic and the
  complete nominal subclass hierarchy remain uncertified.
- Retaining an undecidable exact residual preserves possible solutions but is
  not an exact proof that the candidate is valid. A general exact algebraic and
  branch-aware verification engine remains necessary.
- `PurePoly` now delegates arithmetic and polynomial algorithms to native
  `Poly` backends, with generator-independent unification. Characteristic
  polynomials over coefficient domains without a complete native `Poly`
  backend retain their existing data representation; those domains and all
  constructor options remain an open parity obligation. Reference comparisons
  through inexact coefficient rounding or between `Poly` and `PurePoly` also
  require care: their equality relation need not be transitive.
- Finite permutation-group order and membership use exact incremental
  stabilizer chains, with no former 100,000-element enumeration cutoff.
  `schreier_sims_incremental(slp_dict=True)` now returns exact reconstruction
  witnesses in earlier strong generators. Remaining group/Basic protocols
  are incomplete. Scoped invariant reviews retain the outstanding comparative
  complexity obligations, including witness storage, in the versioned audit.
  Public `Piecewise.args` now exposes reconstructible expression/condition
  pairs with `expr`/`cond` accessors; constructor and structural consumers use
  those pairs. The public constructor is also a class supporting native
  instance predicates, `fromiter`, unbound methods and class queries in
  `has`/`replace`. Held evaluation discards false branches and branches after
  the first true condition, and collapses a single true branch. Nominal
  `Function`/`Expr` inheritance and every `ExprCondPair`/`Basic` protocol
  remain separate obligations.

These gaps prevent a claim of global algorithmic equivalence with SymPy.

## Wrong-result fixes: optics, Gaunt, Poisson process, Maple, codegen (WP-G14)

Fixed to match SymPy 1.14:

- `physics.optics.utils` is a faithful port: `deviation` sign, `fresnel_coefficients`
  returns two values under total internal reflection, `refraction_angle`
  returns `Ray3D` and raises `ValueError` on total internal reflection, and
  lens/mirror formulas handle `oo`.
- `real_gaunt` uses the SymPy U-matrix convention (checked identical to the
  reference on 859 index combinations).
- `maple_code(..., assign_to=y)` emits `:=`; the printer is a port of the
  reference `MapleCodePrinter`. Sparse matrices keep their sparse type through
  `sympify`, and `S.BooleanTrue`/`S.BooleanFalse` are available.
- `codegen.ast.Element` coerces indices and strides to `Tuple`, is no longer an
  arithmetic expression, and `ccode` honours strides and offset.
  `Node.attr_params` returns `None` for a missing attribute.
- `PoissonProcess` conditional probabilities with random-event or `Contains`
  conditions are computed exactly (finite enumeration with a lumped tail
  state) instead of being returned inside a `Piecewise`. `state_space` is
  `Naturals0`; `__add__` and `split` follow the reference, including errors.

Deliberate divergences where the reference is mathematically wrong (symox keeps
the exact value; `X` is `PoissonProcess('X', 3)`):

| Query | SymPy | symox |
|---|---|---|
| `P(Eq(X(t1),2), Ne(X(t1),4) & Contains(t1, [2,4)))` | `1-36e^-6` | `1-54e^-6` |
| `P(X(2)<=3, X(1)>1)` | `5e^-3` | `45e^-3/(2(e^3-4))` |
| `P(X(2)<=3, X(1)>=1)` | `27e^-3/2` | `48e^-3/(e^3-1)` |
| `P(X(3)>2, X(1)<2)` | `e^-6` | `1-23e^-6/2` |
| `P(X(4)<3, Eq(X(1),1))` | `e^-9` | `10e^-9` |
| `P(X(2)<3, X(2)>1)` | `0` | `18/(e^6-7)` |
| `P((X(1)<2) & (X(2)>3))` | `4e^-6` | `4e^-3-77e^-6/2` |
| `P(Eq(X(t1),1) \| Eq(X(t2),0))`, disjoint `Contains` | `e^-6+3e^-3` | `3e^-3+e^-6-3e^-9` |
| `E(X(3), Eq(X(1),2))` | `9` | `8` |

Remaining gaps:

- `fcode` of codegen nodes (`Element`, `Assignment`, `Declaration`) is not
  supported by the native Fortran printer; `fcode(Element)` now raises
  `TypeError` instead of emitting the wrong `x[i, j]`.
- Reprs of other codegen nodes (`Variable`, `Assignment`, `Node`,
  `FunctionCall`, `For`) differ from the reference `Token` repr, and
  `str(none)` is `None` instead of `NoneToken()`.
- A symox `Tuple` now compares equal to the Python `tuple` of its items
  and shares its hash (see the core gaps section at the end).
- Unevaluated process probabilities display internal atoms (`_X(t1)`).
- Outside the process path, joint probabilities of independent discrete
  variables over multivariate relationals (for example
  `P(Eq(A,2) & Eq(B,3))` with plain `Poisson` variables) still return a
  wrong `Piecewise` with dummies outside the `Sum`.

## Wave 4: polys domains and number fields (WP-C4, WP-C6)

Implemented, with SymPy 1.14 values in
`tests/test_missing_polys_wave4_parity.py`:

- `numberfields.modules`, `basis` (`round_two`), `primes` (`prime_decomp`,
  `prime_valuation`, `PrimeIdeal`), `utilities` and `exceptions` are ports of
  the reference pure-Python modules on symox domains, with exact integer HNF
  (`_normalforms.py`) and the `DomainMatrix` helpers they need (`_solve`,
  `flat`, `__setitem__`, `from_rep`, ZZ/scalar division to QQ). The
  reference module API (`Module(col)`, `element_from_rational`, `represent`,
  `find_min_poly`, `kernel(modulus)`) replaces the former custom one.
- `AlgebraicField`: `maximal_order`, `integral_basis(fmt)`, `discriminant`,
  `primes_above`, `galois_group`, `orig_ext`, `to_alg_num` and construction
  from dense coefficient lists (`K([1, 2])`).
- Finite fields: elements print `SymmetricModularIntegerMod<p>(v)` /
  `v mod p`, compare with integers modulo p, and `GF(p).mod`, `.dom`, `.sym`
  exist. `domain=GF(p)` (domain object or string) is the `modulus=p` option
  for `factor`, `factor_list`, `sqf`, `sqf_list`, `gcd`, `lcm`, `cancel`,
  `div`, `quo`, `rem`, `gcdex` and `invert`. `cancel` modulo p now cancels the
  gcd and makes the denominator monic, as `Poly.cancel` does.
- Gaussian domains: `ZZ_I / x` lands in `QQ_I`; `//`, `%` and `divmod` are
  `K.div` in `ZZ_I` and `QQ_I`; `ZZ_I.exquo` failures carry the reference
  message.

Intentional divergences:

- `prime_valuation`: the reference loop (Cohen 4.8.17) reconverts a
  non-integral intermediate matrix to ZZ and raises `CoercionFailed` on
  integral ideals at ramified primes or primes of residue degree > 1, such
  as `2*ZK` in `QQ(sqrt(2))`. symox keeps the reference result whenever it
  succeeds. Otherwise it returns the valuation given by ideal containment,
  bounded by the norm.
- `GF(p**k)` and `GF(n)` for composite `n` are now the reference ring
  `Z/nZ` (see the core gaps section); polynomials over it: see the section
  on polynomials over `Z/nZ` at the end.

Remaining unsupported:

- Galois groups of degree 5 and 6 (for example `QQ.cyclotomic_field(7)`)
  raise `NotImplementedError`. The resolvent tables of
  `galois_resolvents` are not ported.
- `to_sympy` of fields generated by an aliased `AlgebraicNumber`
  (`alg_field_from_poly(.., alias)`, `cyclotomic_field`) returns the root
  expression, not the alias. The native `AlgebraicNumber` does not survive
  arithmetic. Field elements print natively, not as `ANP(...)`.
- `QQ.algebraic_field((T, theta))` with a symbolic `theta` is supported
  since c-doc (2026-10-23) for an irreducible rational `T`; see the c-doc
  section at the end.
- `GF(p, symmetric=False)` is the Python `Z/nZ` domain with the reference
  `ModularIntegerMod{p}` elements.
- The reference `ModularInteger % ModularInteger` is not supported.
- `ModuleHomomorphism.kernel()` over ZZ with a non-integral rational
  nullspace raises `CoercionFailed`, like the reference.
- The audited `DomainMatrix`/`DDM`/`SDM`/`DFM` method examples (`fflu`,
  `lll`, `qr`, `clear_denoms`, `primitive`, `charpoly_factor_*`, ...)
  remain open.

## Special-function evaluation and rewrites (wave 4, missing_special)

Implemented, following the SymPy 1.14 `eval` and `_eval_rewrite_as_*`
methods:

- Values at infinity: `erfi(-I*oo) = -I`, `polygamma(0, -oo)` and
  `polygamma(0, +-I*oo)` give `oo`, `loggamma(+-I*oo) = zoo`,
  `uppergamma(a, oo) = 0`, `subfactorial`, `lucas`, `fibonacci` and `bell` at
  `oo`, `airybi`/`airybiprime(-oo) = 0`, `bessely(nu, +-I*oo)` and
  `besselk(nu, oo or +-I*oo)`.
- `bernoulli(n) = 0` for odd `n > 1`, including symbolic odd `n`.
- `conjugate(Ynm(n, m, theta, phi)) = (-1)**m*Ynm(n, -m, theta, phi)` (so
  `Ynm_c` evaluates as in SymPy).
- `SingularityFunction` uses SymPy's constructor (nan propagation, real
  arguments, orders down to -4, non-integer nonnegative orders), and its
  Piecewise, Heaviside and DiracDelta rewrites apply anywhere in an
  expression.
- `Heaviside.fdiff` and `DiracDelta.fdiff`, with SymPy's
  `ArgumentIndexError` message.
- `doit()` on a relation rebuilds it from its arguments, so
  `Eq(2, 5, evaluate=False).doit()` is `False`.
- Rewrite targets of Ei, expint, li, Li, Si, Ci, Shi, Chi, erfc to erf, all
  four Airy functions, the Bessel family (besselj, bessely, besseli, besselk,
  jn, yn, with SymPy's integer-order guard on yn), betainc and
  betainc_regularized, zeta, dirichlet_eta, polylog, elliptic_k/e/f/pi, and
  catalan.
- `elliptic_pi(1, z, m)`, `elliptic_pi(n, z, 0)` and `elliptic_pi(n, z, n)`
  evaluate in closed form.
- `stieltjes` and `riemann_xi` use SymPy's `eval`: `stieltjes(n)` keeps a
  single argument, negative or non-integer numeric orders give `zoo`, and
  `riemann_xi(s)` evaluates whenever `zeta(s)` does.
- `hyper(ap, bq, z, evaluate=False)` sorts its parameters without cancelling
  them.

Intentional divergences, where symox keeps the correct value:

- `polygamma(-1, oo)` is `oo` in symox. SymPy returns 0, but the order -1
  polygamma is an antiderivative of the digamma function and diverges.
- `besselj(1, -I*oo)` is `-oo*I` in symox, from
  `besselj(1, -I*y) = -I*besseli(1, y)`. SymPy returns nan.
- `Ci(-I*oo)` raises `UnboundLocalError` in SymPy. Symox leaves it
  unevaluated.
- `diff(Heaviside(x - a), x, k)` gives `DiracDelta(x - a, k - 1)`, whereas
  SymPy sometimes normalises the argument to `a - x` with a sign. The two
  forms are equal.

Remaining unsupported items. These raise or stay unevaluated and never
return a wrong value:

- Series of the Airy and complete or incomplete elliptic functions raise
  `NotImplementedError`.
- `limit(harmonic(n, m), n, oo)` stays an unevaluated `Limit`
  (`expand_func(harmonic(n - k))` is now expanded, see the c-final section).
- Derivatives of `expint` in its order, `elliptic_pi` in its characteristic
  and `assoc_laguerre` in its parameter stay as `Derivative`.
- `dirichlet_eta(s, a)` (two arguments) is not implemented.
- `motzkin.is_motzkin`, `Integer.factors()` and
  `partition(q)` for a negative integer symbol (SymPy gives 0) are missing.
- Printing differences remain: the term order of
  `li(z).rewrite(meijerg)`, and `4*(1 - m)` versus `4 - 4*m` in
  `elliptic_f(z, m).diff(m)`.
- Native functions given `nan` sometimes return `nan` where SymPy stays
  unevaluated.
- Some infinite complex arguments depend on how `exp(I*oo)` is handled in
  the core: `expint(0, +-I*oo)`, `lowergamma`/`uppergamma(2, +-I*oo)` and
  `elliptic_pi(+-I*oo, 0)`.
- The following still differ at infinite degree or argument, or in their
  error type:
  - legendre and chebyshevt at infinite degree;
  - `DiracDelta(-oo)`;
  - the error types of laguerre, partition, bell and factorial2;
  - `beta(1, +-oo)`;
  - `Chi(-I*oo)`;
  - the odd-part extraction in `elliptic_e(-I*oo, 1)`.

## Wave 4: core matching, series, sums and fancy sets (missing_core)

Implemented as SymPy 1.14 ports (tests in
`tests/test_missing_core_wave4_parity.py`):

- `Basic.matches`/`Basic.match`. `Wild`, `Pow`, `Add`, `Mul` and `MatMul`
  `matches` accept an old-match dictionary and the `old` flag. `MatMul`
  matching goes through the commutative-part matcher. `WildFunction` keeps
  `nargs` per instance.
- `series(n=None)` lazy generator. `series`/`nseries` with `logx`.
  `Add.extract_leading_order`.
- `Abs` series where the sign of the leading coefficient is not known: when
  the native engine refuses, the expansion is retried in a positive dummy, as
  `Expr.series` does. `Abs(y*x).series(x)` gives `x*Abs(y)`.
- `N`/`evalf` of `Sum(f, (k, a, oo))` via `hypsum` (fixed-point summation or
  Richardson/Shanks extrapolation). Euler-Maclaurin is the fallback.
  Divergence raises SymPy's `ValueError`. `Sum.euler_maclaurin`.
- `expand_power_base` (the full `Pow._eval_expand_power_base`, including
  `force` and `deep`).
- `Range`: indexing, slicing, `reversed`, iteration and the canonical stop of
  symbolic ranges. `ImageSet`: one-variable membership (solving over the base
  set) and iteration.
- `ComplexRegion` now returns `CartesianComplexRegion`/`PolarComplexRegion`
  (polar angles normalized). Added `from_real`, plus the reference union and
  intersection rules (two regions of the same kind, or a region and a real
  set).
- Bordered `hessian(f, vars, constraints)`, `Trace`/`Determinant` simplify
  (`doit` of an unhandled head), and `FunctionMatrix` with string lambdas or a
  `KroneckerDelta` head. Sparse `RL`/`CL`.
- `solve_triangulated` follows the reference algorithm: without `extension`
  only the ground-domain roots are returned, honouring `domain=`.

Intentional divergences (SymPy is wrong):

- Symbolic `Range` canonicalization: when the stop has a rational offset,
  test `((end - stop)/step).is_negative` rather than `(end - stop).is_negative`.
  SymPy drops an element for negative steps: it gives `Range(n, n - 20, -3)`
  6 elements instead of 7.
- `Range` with a finite start and an infinite stop:
  - Integer indices test `(n - i).is_extended_nonnegative`. SymPy raises
    `IndexError` for `Range(1, oo)[5]`.
  - Slices are indexed from the start. SymPy's mapping onto the reversed
    range is off by one: `Range(1, oo)[0:3]` gives `Range(0)` and
    `Range(1, oo)[2:5]` gives `Range(3, 5)`.
  - SymPy's `ValueError` for an ambiguous unbounded re-stride (`[::2]`) is
    kept.

Remaining gaps:

- `solve_triangulated(..., extension=True)` raises `NotImplementedError`. The
  reference lifts roots into `CRootOf` algebraic-field towers, and the native
  polynomial domains cannot yet evaluate polynomials at `CRootOf` generators.
- The `Abs` series of an argument whose leading coefficient has an unknown
  sign, and which does not factor out of `Abs` in a positive dummy, still
  raises `NotImplementedError`. Example: `Abs(x*y + x**2)`; SymPy returns
  `sign(y)*(x*y + x**2)`. `Abs(y*x).nseries(x)` gives `x*Abs(y)`, where SymPy
  gives the equal form `x*y*sign(y)`.
- The default `expand()` and the native `expand_power_base_deep` (used by
  meijerint) do not yet split `(2*x)**y` the way the public
  `expand_power_base` does.
- Other `ComplexRegion` operations raise `NotImplementedError`: complements,
  mixed polar/Cartesian unions, and unions with non-real sets. Membership of
  float points also differs from SymPy, which raises `TypeError` (the
  message differs). `ImageSet` membership for maps of several variables is
  still decided only by the native predicate.

## Wave 4: DomainMatrix, DDM, SDM and DFM API (WP-C5)

`symox.polys.matrices` now ships the reference `domainmatrix`, `ddm`, `sdm`,
`dense`, `rref`, `lll`, `normalforms`, `linsolve`, `dfm`/`_dfm` and
`domainscalar` modules, ported over the native domain elements. The dense
polynomial helpers they need (`dup_mul`, `dup_transform`,
`dup_clear_denoms`, `dup_factor_list`, the factor sorting and Tarjan's
strongly connected components) are in `_dense_poly.py`. This covers
`cancel_denom*`, `clear_denoms*`, `content`, `primitive`, `eval_poly*`,
`charpoly_factor_*`, `charpoly_base`/`charpoly_berk`, `adj_poly_det`,
`fflu`, `lu`, `lu_solve`, `qr`, `lll`/`lll_transform`,
`nullspace(divide_last=...)`, `scc`, `solve_den*`, the `to_/from_` dok, dod,
flat and flat_nz constructors, `from_dict_sympy`, and the DDM/SDM helper
functions (`sdm_irref`, `ddm_irref_den`, `sdm_nullspace_from_rref`,
`sdm_particular_from_rref`, ...). All 85 audited polys.matrices records match
the reference, against 28 before. The 984 docstring outputs of the vendored
modules also match. Algorithms keep the reference domains, fraction-free
paths and storage formats. Like the reference, dense and sparse matrices
compare unequal.

Intentional divergences:

- Entries are converted to native elements of the matrix domain on
  construction: a Python `1` in a QQ matrix is stored as `QQ(1)`. The
  reference keeps the raw `int`. Algebraic and EXRAW elements print natively,
  for example `QQ<sqrt(2)>(3)`, not `ANP(...)`.
- `DDM.nullspace` and `SDM.nullspace` over ZZ are fraction-free and exact.
  The reference returns floats there, for example `-2.0`.
- `rowspace` returns the pivot rows of the transpose. The reference keeps
  the first `rank` rows, which can be dependent, for example the zero row of
  `[[0, 0], [1, 2]]`.
- `M**-n` returns `inv(M)**n`, kept from earlier symox releases. The reference
  raises `NotImplementedError('Negative powers')`.
- For `rref` and `inv` over EX and EXRAW, the expression kernel is used, so
  that `sin(x)**2 + cos(x)**2 - 1` is recognized as a zero pivot.
- Factors are sorted with `sort_key` over EX and EXRAW. With the reference
  `<` comparison, `charpoly` raises `TypeError` for symbolic EXRAW entries.
- Over polynomial and fraction-field domains, `eval_poly`, `adj_det` and
  `solve_den_charpoly` return the result. Over EX, `eval_poly_mul` does too.
  In these cases the reference raises `TypeError`, from
  `DomainMatrix + PolyElement` or from `'method' object is not iterable`.

Remaining unsupported:

- `DFM` requires python-flint, as in the reference without flint:
  `to_dfm()` raises `NotImplementedError`.
- `sympy.polys.matrices.eigen` is not shipped. It needs
  `agca.extensions.FiniteExtension`.

## Missing sympy.geometry features (wave 4, WP-C8) on 2026-10-22

Implemented, with SymPy 1.14 values hard-coded in
`tests/test_missing_geometry_wave4_parity.py` (probe: 27/236 matching before,
234/236 after; the two remaining lines are listed below):

- Area moments: `second_moment_of_area`, `polar_second_moment_of_area`,
  `section_modulus`, `first_moment_of_area` and `cut_section` for polygons,
  triangles and regular polygons, plus the ellipse and circle versions.
- `reflect` for points (exact formula), linear entities, polygons, regular
  polygons (SymPy's `closing_angle` rotation with a negative radius), circles,
  axis-aligned ellipses, parabolas and curves. `Ellipse.evolute`,
  `director_circle`, `auxiliary_circle` and `normal_lines` (including `prec`)
  are also implemented.
- `plot_interval` and `random_point(seed)`. Lines, rays, segments, ellipses
  and planes draw the same numbers as SymPy, using `random.Random(seed)` or the
  global `random` module.
- `Plane.perpendicular_plane`, `parallel_plane`, `are_concurrent`,
  `is_coplanar`, `arbitrary_point(u, v)`, `parameter_value` and
  `random_point`. Also added: `Point3D.direction_ratio`, `direction_cosine`
  and `are_collinear`, and the Point helpers `canberra_distance`, `is_zero`,
  `is_nonzero`, `origin`, `project` and `orthogonal_direction`. Points can
  now have more than three dimensions.
- `Line3D.equation` now returns a `Tuple`. Other additions: `Line(Segment)`,
  `Ray(angle=...)`, `coefficients`, `bisectors`, `parameter_value`,
  `Ray.xdirection`, `ydirection` and `zdirection`, and `Ray2D.closing_angle`.
- `Polygon(c, r, n=..., rot=...)` and `Triangle(sss=/asa=/sas=)` (`None` when
  SSS has no solution, `IndexError` for parallel ASA rays). Also added:
  `Polygon.bisectors(prec)`, the Piecewise `Polygon.arbitrary_point`,
  `Triangle.bisectors()` as a method, and `RegularPolygon.inradius`,
  `circumradius`, `length` and `spin`. `RegularPolygon` now reduces the
  rotation modulo `2*pi/n` and accepts a negative radius. `Polygon.angles`
  keeps SymPy's key order.
- Exact conic intersections. Ellipse/ellipse, parabola/ellipse and
  parabola/parabola intersections use a resultant elimination, then an exact
  membership check of each root, returned in SymPy's `ordered` order. If a
  root is missing or its reality cannot be decided, the result is
  `NotImplementedError`, never a partial list.
- When a Point is passed to a constructor, its coordinates are no longer
  re-evaluated, so floats survive as in SymPy. `Line(p, slope=s)` keeps a
  float slope.

Intentional divergences:

- `Curve.reflect` returns the reflected curve. SymPy 1.14 returns the curve
  unchanged, because it only maps Point atoms, which is mathematically wrong.
- `Point3D.reflect` keeps the z coordinate when the line lies in the xy-plane.
- `Plane.parameter_value(p, t)`: when the system solver gives up, symox
  solves one trigonometric component at a time and checks every candidate
  exactly. The result is a valid parameter, but it can differ from SymPy's
  by a multiple of `2*pi`; for example, `Plane((1,2,3),(1,1,1))` at `5*pi/4`
  gives `-3*pi/4`.
- The bisector ray points of non-convex polygons equal SymPy's in value, but
  their radical forms can differ.

Remaining unsupported (they raise like SymPy, never a wrong value):

- `Ellipse.reflect` across an oblique line raises `NotImplementedError` with
  the reflected equation, as SymPy does. `Circle.evolute` raises
  `NotImplementedError`, as SymPy does.
- Conic intersections whose eliminant cannot be fully solved by `roots`, or
  whose roots have undecidable reality, raise `NotImplementedError`.

## Wave 4: statistics joint events and symbolic queries (WP-C10)

Implemented (tests in `tests/test_missing_stats_wave4_parity.py`):

- Joint events over independent probability spaces. `P` of `&` and `|`
  events, including discrete x discrete and discrete x continuous events,
  now factors over independent blocks. `P(Eq(A, 2) & Eq(B, 3))` with
  independent Poisson variables used to return a wrong `Piecewise`; it is
  now the product of the two probabilities. `E` and `P` drop condition
  blocks that are independent of the queried expression. `E` of a product
  over independent blocks factors.
- `piecewise_fold` no longer moves a `Piecewise` out of a `Sum`, `Product`
  or `Integral` when a condition uses the bound variable. Before, a
  `simplify` of such sums was wrong.
- `rewrite` (Integral, Sum, Probability, Expectation) and `expand` of
  `Expectation`, `Variance`, `Covariance`, `Probability`, `Moment` and
  `CentralMoment`, plus their `doit`, `evaluate_integral` and
  `evaluate_sum`, follow the reference.
- A random symbol used as a condition means that it is observed, as in
  SymPy's `Eq(Y, Y.symbol)`. This applies to `given`, `E`, `P`, `density`
  and `variance`: `E(X + Y, Y) = Y`, and `density(X + Y, Y)` is the shifted
  normal. A component of a joint distribution raises
  `NotImplementedError`. `given(expr, cond)` returns `expr` when `cond` is
  independent of it.
- `given(Tuple(...), cond)` is supported when at most one entry depends on
  the event and that entry shares no space with the others. Other cases
  raise `NotImplementedError`.
- `RandomSymbol.subs` on distribution parameters rebuilds the variable.
  Before, `Normal('X', m, 1).subs(m, 2)` silently kept the mean `m`.
- `IdealSoliton` and `RobustSoliton` accept a symbolic `k` (and a symbolic
  `c`). `density(.).dict` is the held `Density(...)`, and `.subs(...).doit()`
  evaluates it. The ideal soliton mass function is symbolic. `E` and `P` of
  a symbolic-size soliton raise `NotImplementedError` (SymPy raises
  `TypeError`).
- `coskewness` with a condition (the audited doctest) now runs in a fraction
  of a second; before, it took more than 60 s. The cause was the factored
  expectation of independent blocks.
- `DiscreteMarkovChain`:
  - Float rows are checked exactly, as in the reference.
  - Symbol state labels are decided by index.
  - A `MatrixSymbol` of explicit size is used entrywise, both as `trans_probs`
    and through `TransitionMatrixOf`.
  - A `MatrixSymbol` of symbolic size has state space `Range(0, n)` and
    answers `P(Eq(Y[a], b), Eq(Y[c], d))` as `(T**(a - c))[d, b]`.
  - `stationary_distribution(condition_set=True)` works for a `MatrixSymbol`.
  - Added `limiting_distribution` and `fixed_row_vector`.
- Added `ContinuousMarkovChain.limiting_distribution()`.

Probe results against SymPy 1.14 (scratch probes):

- Joint and conditional `P`/`E`, 24 cases: every case now matches the
  reference value, or is correct where the reference raises or is wrong.
  The top-priority `Piecewise` result was wrong before the fix.
- Observed-variable and tuple `given`, 25 cases: 22 raised or were wrong
  before. After the fix, 11 match. Of the 14 textual differences:
  - 2 are equivalent forms and 2 are `_X` dummy printing.
  - In 2 cases both sides raise.
  - In 4 cases symox returns a correct value where SymPy raises.
  - 3 are the divergences listed below, and 1 is unsupported (the shared
    `Eq(X + Y, 3)` tuple).
- Markov chain doctest variants, 20 cases: 17 failed before and 9 differ
  after. Of those 9, 7 are open (listed below), 1 is the printing of
  `(T**5)[1, 2]` and 1 is an equivalent form.
- Soliton doctests: these raised before. They now match the reference for
  `IdealSoliton`, and for `RobustSoliton` everywhere except the
  normalization divergence.

Intentional divergences:

- `P(A | B)` for non-disjoint independent discrete events: SymPy adds the
  two probabilities. symox uses inclusion-exclusion.
- `Variance(X).doit()` and `Covariance(X, Y).doit()` evaluate where the
  reference keeps them unevaluated.
- `P(X > Y, Y)` evaluates to `1/2 - erf(sqrt(2)*Y/2)/2`, where SymPy returns
  it unevaluated. `P(Y > 1, Y)` is the indicator `Piecewise`, where SymPy
  returns a `BernoulliDistribution`.
- `given(Tuple(X, Y), X > 0)`: the entries are independent, so
  `independent` is True. SymPy reports them as dependent through the shared
  conditional space. `dependent(A, A + B)` for discrete variables is True;
  SymPy raises `NotImplementedError`.
- `RobustSoliton` normalizes over the support `1..k`, so the mass sums to 1
  (Luby). The reference `Z` sums the spike terms up to `round(k/R)`, beyond
  `k`, and its masses do not sum to 1. For example, `k = 10` with
  `c = 0.05` sums to `0.975`. Numeric values therefore differ from the
  reference doctest.

Remaining unsupported:

- `Str` states of a Markov chain raise `NotImplementedError`. The native
  `Eq` does not accept `Str`.
- Markov queries at symbolic times (`P(Eq(Y[a], b), Eq(Y[c], d))` on an
  explicit matrix) raise `NotImplementedError`. SymPy builds a symbolic
  matrix power.
- `stationary_distribution()` and `limiting_distribution()` of a
  `MatrixSymbol` chain without `condition_set=True` raise
  `NotImplementedError`. The native `Lambda` rejects matrix-symbol
  signatures.
- Not implemented and still raising `AttributeError`:
  `DiscreteMarkovChain.canonical_form`, `communication_classes`,
  `decompose` and `fundamental_matrix`.
- `given(Tuple(...), cond)` with an event shared by several entries, such
  as `Eq(X + Y, 3)`.
- The native summation of `k*3**k/k!` over `k >= 0` stays unevaluated, so
  `Expectation(Poisson).rewrite(Sum).doit()` is a `Sum`. `E` itself is
  exact. Sums of an expression times an indicator `Piecewise` of the
  summation variable stay unevaluated; SymPy also leaves them unevaluated.

## Wave 4: permutation, finitely presented and polycyclic groups (WP-C7)

`PermutationGroup` now inherits the reference algorithms: the methods of
`sympy.combinatorics.perm_groups` are ported and only the constructor,
`order`, `contains`, `orbit(s)`, `elements`, `make_perm`, equality and
hashing keep the native fast paths (the native stabilizer chain gives the
order and membership). This covers `center`, `centralizer`,
`normal_closure`, `derived_series`, `lower_central_series`,
`composition_series`, `is_primitive`, `minimal_block(s)`,
`sylow_subgroup`, `subgroup_search`, `schreier_sims_incremental`,
`schreier_sims_random`, `coset_factor`, `coset_rank`, `coset_unrank`,
`presentation`, `strong_presentation`, `abelian_invariants`, `baseswap`,
`commutator`, `equals`, `index`, `is_alt_sym`, `is_elementary`,
`is_normal`, `is_subgroup`, `orbit_rep`, `pointwise_stabilizer`,
`schreier_vector`, `copy` and the reference `str`/`repr`. The modules
`util`, `testutil`, `free_groups`, `fp_groups`, `coset_table`,
`homomorphisms`, `rewritingsystem(_fsm)`, `pc_groups`, `tensor_can` and the
named groups (`DirectProduct`, `AbelianGroup`, ...) are ports of the
reference. Base, strong generators, basic orbits and transversals are the
reference ones. `Permutation.rmul` accepts sequences (`__rmul__`) and
`Partition.sort_key` is available.

The reference docstring examples of the combinatorics group modules went
from 105 to 225 of 249 objects passing (1145 to 1534 of 1588 examples), and
a probe of 98 group calls from 22 to 96 values equal to SymPy 1.14 (the two
others are `random_stab`, random, and an error message of
`Permutation(3)(0, 1)` resizing). Expected values are in
`tests/test_missing_groups_wave4_parity.py`.

`tensor_can.canonicalize` is the reference Butler-Portugal algorithm, so
the previous limit of 100000 orbit states is gone, and
`testutil.graph_certificate` is the reference certificate (vertices other
than `0, ..., n - 1` are renumbered in sorted order, where the reference
raises `IndexError`). `fp_groups` computes abelian invariants with an exact
integer Smith normal form.

Randomized outputs: `normal_closure`, `sylow_subgroup`, `random_stab`,
`schreier_sims_random`, `polycyclic_group` and the presentations use
`random_pr` and `randrange` like the reference, so their generators differ
from run to run in both libraries. The tests check orders, group equality,
normality, `_verify_bsgs` and generator counts. `schreier_sims_random` is
Monte Carlo in both libraries (with `consec_succ=8`, about 1 % of runs give
an incomplete BSGS).

Intentional divergences (SymPy 1.14 demonstrably wrong):

- `Collector.constructive_membership_test` returns the exponent vector of a
  member. The reference compares the sifted permutation with `1`, which is
  never equal, and returns `False` for every element.
- `Collector.induced_pcgs` and the membership test invert the leading
  exponent modulo the relative order (the reference multiplies by the float
  `e**-1`), sift elements of the last depth and add the power of each new
  element. The reference returns `[(0 1 2 3)]` for the subgroup generated by
  `(0 1 2 3)` in S4 (missing `(0 2)(1 3)`) and one generator for the Sylow
  3-subgroup of S8 (order 9). A non-invertible leading exponent (relative
  order 4) raises `NotImplementedError`. `induced_pcgs` no longer empties
  the list it receives.
- `Collector.collected_word` keeps the whole power relator (see the
  geometry and combinatorics section above); `Collector(pcgs, [], orders)`
  builds the tail series of the sequence and checks it, where the reference
  raises `IndexError`.
- `schreier_sims_random` of a group with a one-point base (`CyclicGroup(7)`)
  returns the BSGS; the reference raises `IndexError`.
- `schreier_sims_incremental` of identity-only or empty generators returns
  the empty BSGS (`([], [])`, with `{}` for `slp_dict=True`); the reference
  raises `IndexError`.

The earlier note that `induced_pcgs` and `constructive_membership_test`
keep symox algorithms is superseded by the ported and corrected reference
algorithms above.

Remaining gaps:

- `tetrahedron.pgroup` has other generators than the reference, so the
  `PermutationGroup.generate` example lists the elements in another order.
- `group_numbers._holder_formula` is not ported, and the `permutations`
  module examples that depend on `Permutation` printing modes,
  `print_cyclic`, error messages and `AppliedPermutation` still differ
  (24 objects of the docstring probe).
- `coset_enumeration_r` of `B(2, 4)` takes about as long as in the
  reference (4 s); free group words are pure Python.
- `Collector` and `PolycyclicGroup` printing recurses without end, as in
  the reference.

## Polynomial systems with quartic radical roots (wave 4 performance regression)

`solve([x**2 + y**2 - 4, y - x**4 + 1], [x, y])` ran for more than 200 s
(SymPy 1.14: 12 s). Root cause, in order of cost:

- Commit 927d768 removed the size limits of `expand(complex=True)`.
  `roots_binomial` expands `sqrt` of each root of the decomposed quartic with
  it; the bounded expansion used to fail and keep the unexpanded root (the
  archived 2.8 ms timing came from that form).
- The complex split ran its sign queries in exact mode only, so `sqrt(s)` with
  `s` a positive radical sum was not proved nonnegative and was expanded
  through `atan2`, `Piecewise` and multinomial squares (18 MB per root).
  SymPy decides these signs with its certified numerical evaluation; the split
  now follows the caller's numeric refinement policy, so `sqrt(s)` stays whole
  as in the reference.
- Univariate `solve` verified the roots of a numerator polynomial in the
  unknown by numerical substitution. The reference `_solve` sets
  `check = False` for them (only denominators are checked); symox does the
  same now.
- The system check (`checksol`) evaluated exact-zero residuals with the
  4096-bit cancellation budget, recomputing `ln 2` at every precision. The
  check now caps the cancellation at its own precision (a larger cancellation
  is already within its tolerance) before its chopped evaluation, and `ln 2`
  is cached per precision.
- With the numeric split, `arg` of radical complex numbers (for example
  `arg((1 + I*sqrt(-1 + 10*10**(1/3)))**(-3))`, used by `phase_margin`)
  needed `extended_real` of `sqrt(c)`, which the exact assumption system
  leaves open; the `arg`/`atan2` evaluator now accepts a certified real
  number as extended real, and these values match the reference text.

Timings (symox before / after, SymPy 1.14), same solution sets checked at
50 digits:

| System | before | after | SymPy |
|---|---|---|---|
| `x**2+y**2-4, y-x**4+1` | >60 s | 4.1 s | 12.0 s |
| `x**2-y**2-1, y-x**4+2` | >60 s | 8.9 s | 12.0 s |
| `x**2+y**2-4, y-x**4` | >60 s | 3.4 s | 12.2 s |
| other `data/cases3.json`-style degree 2-4 systems (9) | 0.7 ms - 3.4 s | 0.8 ms - 0.1 s | 0.1 - 68 s |

Intentional divergence: for `x**2 + y**2 - 4, y - x**3 + 1` SymPy returns
`[]`, which is wrong; symox returns the 6 `CRootOf` solutions (zero residuals).

Remaining gaps:

- symox solves these systems by substitution (`y = x**4 - 1` expanded at the
  radical `x` roots) instead of the reference Groebner path, so the printed
  result is much larger (about 1 MB for the first system, SymPy 16 kB) while
  being equal.
- Most of the remaining time is the numerical system check on those large
  expressions.
- Not caused by this change: `limit_seq(Sum(2**m/m, (m, 1, n))/(2**n/n), n)`
  returns `None` in about a quarter of fresh processes (the reference: `2`),
  because the Stolz-Cesaro path depends on the process-random `Dummy` ordering.

## Wave 4: sympy.physics.quantum and secondquant (WP-C11, missing_quantum)

Regression tests: `tests/test_missing_quantum_wave4_parity.py` and
`tests/test_missing_quantum_secondquant_wave4_parity.py`. Both hard-code
SymPy 1.14 outputs. Probe score: 42/126 before, 114/126 after.

Implemented (ports of the SymPy 1.14 algorithms):

- `Operator`/`Ket`/`Bra` `.label`, `dual_class`, `Operator.inverse`/`inv`.
  `UnitaryOperator`: `U*Dagger(U)` and `Dagger(U)*U` give 1.
- `Commutator`/`AntiCommutator` `doit` and `expand(commutator=True)`, and
  `Dagger` of both.
- `qapply` inside `TensorProduct` factors. `Bra*Ket` and `InnerProduct`
  printing.
- `hbar` is a `NumberSymbol`: printing order and `evalf` at any precision
  match SymPy.
- `represent` follows SymPy's dispatch: `_represent_<Basis>`, then the
  default basis, then `rep_innerproduct`/`rep_expectation`. Products go
  through the index/unities path and `integrate_result`. Also ported:
  `rep_innerproduct`, `rep_expectation`, `integrate_result`, `get_basis`
  and `enumerate_states`. `XKet`/`PxKet` default labels and
  `.position`/`.momentum` were added.
- `CG`/`Wigner3j`/`6j`/`9j` keyword arguments, attributes and `is_symbolic`.
  `WignerD` is a held function (6 arguments, `doit`, `simplify` evaluates
  it). Also added: `Rotation` printing and `inverse`, and `Rotation.D`/`d`.
- `couple`/`uncouple` for any number of spins, ported from SymPy, with
  explicit `jcoupling`. Coupled states take `(j, m, jn, *jcoupling)`, use
  SymPy's validation and errors, and provide `coupling`, `coupled_jn`,
  `coupled_n`, `uncoupled_class`/`coupled_class` and SymPy's printing.
- Hilbert spaces. Tensor products flatten and merge equal neighbours into
  powers, giving `C(2)**3`. Direct sums flatten and print with `+`.
  `FockSpace` prints `F`. The new attributes are `.spaces`, `.base` and
  `.exp`.
- `DifferentialOperator.variables`, `len(Qubit)` and `Qubit.dimension`.
- New module `sympy.physics.quantum.pauli`. It provides `SigmaX/Y/Z`,
  `SigmaMinus/Plus` (powers, `Dagger`, commutator and anticommutator
  tables, `represent`, `qapply` on `SigmaZKet`), `SigmaZKet`/`SigmaZBra`
  and `qsimplify_pauli`.
- Secondquant, as ports:
  - `apply_operator` on `B`/`Bd`/`F`/`Fd`.
  - Fock-state `up`/`down` with Fermi-sea semantics (holes, Pauli and
    hole-count zeros, sorting sign, Kronecker guards on symbolic indices).
  - `VarBosonicBasis`/`FixedBosonicBasis`.
  - `Commutator.eval`/`doit`, including `wicks=True`.
  - `AntiSymmetricTensor` canonical sorting, sign, `symbol`/`upper`/`lower`.
  - Fermi-level-aware deltas in `evaluate_deltas`/`contraction`/`wicks`.
  - `wicks` expansion and the full `substitute_dummies` algorithm.
  - `above_fermi`/`below_fermi` symbol assumptions.

Intentional divergences (same value, or SymPy is wrong):

- Operator sums print their terms in a different order (for example,
  `B*[A,C] + [A,B]*C`).
- Operator coefficients are simplified, so `exp(a)*exp(b)` prints combined
  (`Rotation` applied to states). `TensorProduct` distributes over sums
  when built (`AxC + BxC` rather than `(A + B)xC`). `qsimplify_pauli`
  returns expanded sums such as `I/2 + I*SigmaZ()/2`.
- Some results are an equal unitary form: `U*Dagger(U)**2` gives
  `Dagger(U)` where SymPy gives `U**(-1)`. `simplify(WignerD(...))` can
  print `(-1)**(1/4)` for `exp(I*pi/4)`.
- `Ket.dual_class()` returns symox's own class, so the module path in its
  repr differs.
- Bra-ket products of orthogonal states such as `SigmaZBra(0)*SigmaZKet(1)`
  evaluate to the Kronecker delta immediately. SymPy keeps `<0|1>` until
  `qapply`.
- The SymPy Pauli rules that ignore the spin name are wrong:
  - `Commutator(SigmaZ('a'), SigmaMinus('b')).doit()`: SymPy returns
    `-2*SigmaMinus()`; symox returns 0.
  - Anticommutators of differently named spins: symox keeps
    `A*B + B*A` instead of SymPy's same-spin value.
- `represent(SigmaX(), format='numpy')` returns the numpy matrix; SymPy
  raises `NotImplementedError`.
- `conjugate(<b|k>)` is held outside `Dagger`.
- `srepr` and the type of `above_fermi`/`below_fermi` symbols differ.
- Secondquant:
  - `fermi_level` is stored on each state, and `Dagger` keeps it.
  - `up`/`down` accept plain Python integers (SymPy raises
    `AttributeError`).
  - The printed argument order of an unevaluated `Commutator` uses an
    approximation of `sort_key`.

Remaining unsupported (they raise, never a wrong value):

- `couple`/`uncouple` with symbolic spins raise `NotImplementedError`. SymPy
  returns a held `Sum`, which symox cannot build over operator expressions.
- `XOp` represented in `PxKet` and `PxOp` in `XKet` raise
  `NotImplementedError`. SymPy builds a one-argument `DifferentialOperator`
  there.
- `represent` of products such as `XOp()*XKet()` returns a held `Integral`
  of a `DiracDelta` product. The value is correct, but integrate does not
  apply delta sifting there.
- Secondquant:
  - The core `KroneckerDelta` does not vanish across the Fermi level when
    built and has no `is_only_above_fermi`/`is_only_below_fermi`.
  - LaTeX and `srepr` of operator products print `OrderedProduct`.
  - `Sum` over operator expressions (`matrix_rep` in the oscillator
    example).
  - `NO.subs` with operator keys works only at the top level.

## Wave 4: control, Quaternion and tensor helpers

Implemented:

- `sympy.physics.control.lti` is the reference 1.14 module on symox: `Series`,
  `Parallel` and `Feedback` with `StateSpace` operands, `doit` and
  `rewrite(TransferFunction)`, `MIMOSeries`, `MIMOParallel`, `MIMOFeedback`,
  `TransferFunctionMatrix.elem_poles`, `elem_zeros`, `from_Matrix`,
  `eval_frequency`, `StateSpace.controllable_subspace`,
  `observable_subspace` and `TransferFunction.from_zpk`. 58 of 61 triage
  probes match the reference output exactly (24 before), 35 of 41 reference
  `test_lti.py` tests pass (none collected before).
- `Quaternion` is the reference 1.14 class as a symox host expression:
  `from_euler`, `to_euler`, `axis`, `angle`, `arc_coplanar`,
  `vector_coplanar`, `parallel`, `orthogonal`, `mensor`, `set_norm`, `add`,
  `mul`, `index_vector`, `pow_cos_sin`, `to_Matrix`, `from_Matrix`,
  `from_axis_angle`, `str`/`latex`/top-level `pprint`, `conjugate`, `diff`,
  `integrate`/`Integral` and pickling (which keeps `real_field` and `norm`).
  The top-level `conjugate`, `diff`, `integrate` and `Integral` now call
  the `_eval_conjugate`, `diff` and `_eval_Integral` methods of host
  expression classes. 170 of 171 probe lines match (the other is the
  printing order of a `free_symbols` set); 19 of 20 reference
  `test_quaternion.py` tests pass.
- Tensor: `Idx` bounds (infinite, symbolic, unordered tuples, integer
  labels returned as such) and the reference error messages; a label not
  known to be an integer is now rejected as in the reference;
  `IndexedBase` assumptions; `shape` of an `Integral`; `TensMul.split`,
  `Tensor.split`, `matches`, repeated-index errors; `PartialDerivative`
  dummy names and symbol index orders in `replace_with_arrays`; the
  `from_*`/`conv_*` array expression converter modules. 178 of 199
  reference doctest probes match (142 before).

Intentional divergence:

- `Quaternion.from_rotation_matrix` of a rotation through `pi` (for example
  `diag(1, -1, -1)`) takes the relative signs of the vector components from
  the symmetric part of the matrix and returns a quaternion whose rotation
  matrix is the input. The reference returns the zero quaternion.
- `PartialDerivative(H(i, j), A(j)).replace_with_arrays(...)` contracts
  the derivative index with the right axis (`[2, y]` for
  `H = [[x, y], [x*y, x**2]]`, `A = [x, y]`); the reference gives
  `[x + 1, 0]`, which disagrees with its own manual contraction.

Remaining gaps:

- (Resolved by the core gaps section at the end: `Tuple(1, 2) == (1, 2)`
  is True and `(2 + I).is_complex` is True, so `Quaternion` uses the native
  `is_complex`.)
- A `Quaternion` nested inside another expression does not pretty print;
  `pprint` of control objects raises `TypeError` as before.
- Control: `poles`/`elem_poles` list conjugate roots in another order;
  `evalf` keeps `-1.0*s`; Floats inside matrices and tuples print in full;
  `StateSpace` to `TransferFunction` negates numerator and denominator and
  `MIMOFeedback.doit()` returns an expanded form (same functions, core
  matrix `solve`/`simplify` forms); `TransferFunctionMatrix` from a symox
  `Tuple(Tuple(TF))` fails because the core `Tuple` cannot hold systems;
  `StateSpace.dsolve` does not produce the reference `Piecewise`;
  `phase_margin` keeps the previous symox algorithm, as the reference one
  needs relational truth values that symox does not provide.
- Tensor: array expressions are still the simplified non-`Basic` symox
  implementation (no `C.T.T` simplification, other canonical forms,
  `permutedims`/`tensorproduct` of matrix symbols evaluate to explicit
  arrays, `convert_array_to_matrix` returns an equal explicit matrix,
  `convert_indexed_to_array` of a trace gives `Trace(M)`, diagonal cases
  and `sort_args_by_name` raise `NotImplementedError`). `Integral` of an
  array or a matrix raises `TypeError` (native conversion). `IndexedBase`
  and `Indexed` are still `codegen.ast` classes with another repr.

## Run-to-run determinism of Dummy-dependent results (wave 4)

Fresh Dummy indices start from a seed made of the clock, the process id and
a random value, so that they stay unique across processes and pickles.
Dummy hashes used that index. Add and Mul store their operands in hash
order, so the storage order changed from run to run. `difference_delta`
split a product on its first stored factor. SymPy splits on `args[0]`, in
canonical order. For `n*Sum(2**m/m, (m, 1, n))`, symox then gave
`Sum(2**m/m, (m, 1, n + 1)) + 2**(n + 1)*n/(n + 1)` in about a quarter of
the runs, instead of the reference `2**(n + 1) + Sum(2**m/m, (m, 1, n))`.
Both values are equal, but the Stolz-Cesaro iteration of `limit_seq` does
not converge on the first form, so
`limit_seq(Sum(2**m/m, (m, 1, n))/(2**n/n), n)` returned `None` instead
of 2.

Fixes:

- `difference_delta` splits a product on its first factor in canonical
  order, as the reference does.
- A Dummy hashes the offset of its index from the first fresh index of the
  process. Fresh Dummies of a given program therefore hash, and are stored,
  in the same way in every run. The indices stay unique across processes,
  pickles keep them, and the origin does not change after a fork.

Remaining risks and gaps:

- Rust code that walks the stored operand order (`Expr::args`) where SymPy
  walks `args` in canonical order can still diverge from the reference. It
  is now reproducible, but it can change when an unrelated Dummy is created
  earlier in the run. `traverse::ordered_args` gives the reference order.

## Wave 4: code generation and printers (WP-C14)

Implemented:

- `fcode`, `FCodePrinter`, `rust_code` and `RustCodePrinter` are now ports
  of the SymPy 1.14 printers, replacing the native expression printers.
  This covers Matrix `assign_to`, Piecewise statements, `user_functions`,
  F77/F95 line wrapping, `strict`, the standard checks of `merge` and the
  other `FFunction` intrinsics, and the Rust `TypeCast`, `float_floor` and
  `float_ceiling` overrides.
- `glsl_code`, `jscode`, `julia_code`, `octave_code` and `rcode` are ports
  of the reference printers. Matrix `assign_to` and inline matrices print
  as in SymPy.
- `cfunctions` follow the reference for `fdiff`, `rewrite`, `evalf` and
  printing. `isnan` and `isinf` are boolean functions with the reference
  `eval`. `pycode` and `NumPyPrinter` know the reference math names, and
  unsupported names raise `PrintMethodNotImplementedError`.
- `codegen.ast`:
  - `Print`, `FunctionPrototype`, `FunctionDefinition`, `Declaration`, `For`
    and `While` coerce their fields as SymPy does.
  - `Variable.deduced`, `Type.from_expr` and `cast_check` follow SymPy,
    including the integer and float limit messages.
  - `CodeBlock` gained `left_hand_sides`, `right_hand_sides`,
    `topological_sort` and `cse`.
- `fnodes`:
  - New nodes: `Module`, `Program`, `Subroutine`, `SubroutineCall`, `use`,
    `use_rename`, `Extent`, `dimension`, `array`, `bind_C`, `sum_` and
    `product_`.
  - New attributes: `intent_*`, `allocatable`, `pure`, `elemental`,
    `assumed_extent` and `assumed_size`.
  - `literal_dp` and `literal_sp` give plain `Float` values, as in SymPy.
- `cnodes` `sizeof` and `alignof` return function calls, as in SymPy.

Intentional divergences:

- `fcode(Do(..., concurrent=True))` prints the Fortran 2008 form
  `do concurrent (i = 1:n)`. The reference prints `do concurrent i = 1, n`,
  which gfortran rejects.
- In F77, SymPy raises `F77 init./parameter statem. req. multiple lines.`
  for every `Declaration`, because its `NoneToken` is truthy. Symox prints
  scalar declarations without a value. It still raises for values,
  `value_const` and dimensions, so it never drops a dimension silently.
  `GoTo` without an expression prints `go to 10`; SymPy raises.
- `N(isnan(x).subs(x, Rational(3, 7)))` returns `False`; SymPy raises
  `AttributeError` (`BooleanFalse` has no `evalf`).
- `ccode(Sqrt(x))` raises `PrintMethodNotImplementedError`; SymPy raises
  `KeyError`.
- For the `fnodes.size` example, SymPy prints `sum(a**2)*1d0/size(a)` and
  symox prints `sum(a**2)/size(a)`. Symox builds `sum_` and `size` as
  applied functions named after the intrinsic, so `fcode` knows them as
  calls. The two expressions are equivalent.
- `RustCodePrinter().doprint(FunctionDefinition(...))` keeps symox's typed
  Rust function printer (integer and float arithmetic, loops, returns);
  SymPy raises `PrintMethodNotImplementedError`.
- `rust_code(x + 1)` prints the integer literal as SymPy does (`x + 1`).
  Rust rejects this in `f64` arithmetic. The earlier native printer wrote
  `1.0`.
- Unsupported-object messages name symox classes, for example
  `symox._core.Matrix` and `symox.printing.lambdarepr.NumPyPrinter`.

Remaining gaps:

- Codegen nodes are not `Basic`. The reprs of `Token`, `Variable`,
  `Attribute` and `Type` differ (12 of 215 `ast` doctest examples).
  `str(none)` is `None` and not `NoneToken()`. `intc` is a 32-bit
  `SignedIntType`; SymPy uses `IntBaseType('intc')`.
- `Idx` arithmetic (`i + 1` raises `TypeError`) and contraction loops are
  not supported (both fixed in the c3 section). A Function applied to a Matrix in julia and octave is not
  supported either. These account for the remaining glsl, jscode, julia,
  octave, rcode and codeprinter doctest differences.
- Symox evalf turns `N(-x)` into `-1.0*x`; SymPy keeps `-x`. The Fortran
  printer emulates SymPy's rule for function arguments, but the core
  difference remains.
- Not covered by this work package: smtlib, `codegen.rewriting`,
  `SeriesApprox` and the `algorithms` and `matrix_nodes` doctests (17/18 and
  15/17).

Doctest probe, matching examples before and after:

| Module      | Before  | After   |
|-------------|---------|---------|
| ast         | 163/215 | 203/215 |
| fnodes      | 89/107  | 106/107 |
| cnodes      | 19/24   | 24/24   |
| glsl        | 21/36   | 34/36   |
| jscode      | 17/21   | 19/21   |
| julia       | 26/34   | 31/34   |
| octave      | 24/34   | 31/34   |
| rcode       | 17/23   | 21/23   |
| codeprinter | 59/70   | 64/70   |

The cfunctions probes went from 1/14 to 14/14 (semantics) and from 40/187 to
183/187 (printing). The four printing differences left are the two
corrections above and two error messages that name a class at another path.
`tests/test_missing_codegen_wave4_parity.py` hard-codes 145 reference cases.


## Mechanics, vector and continuum mechanics (WP-C12, wave 4)

`symox.physics.vector`, `symox.physics.mechanics`,
`symox.physics.continuum_mechanics.beam`, `truss` and `arch` are now the
reference SymPy 1.14 sources, rewritten to import symox. This covers:

- the joint arguments `parent_point`, `child_point`, `parent_axis`,
  `child_axis`, `parent_interframe` and `child_interframe` of all six joints;
- `System.form_eoms`, `JointsMethod`, mass matrices, forcing and `rhs`;
- Beam hinges, influence lines, `Beam3D` torsion and the `max_*` methods.

The probe of these modules went from 137 of 244 to 215 of 244 matching
outputs.

Supporting changes:

- Frame coordinates are `CoordinateSym` host symbols, so they print as `A_x`,
  `is_real` is None and they have `.frame`. `_hashable_content` holds the
  frame and the index. The class is marked as a host class, so that content
  is part of the native identity. As a result, the coordinates of two frames
  that share a name stay distinct: `A[0] + B[0]` does not merge, and
  `(A[0]*B[0]).diff(A[0])` is `B[0]`.
- Limits of `SingularityFunction` use its one-sided germ. `lambdify` prints
  it as a `Piecewise`, so it works with numpy, math and mpmath.
- `integrate` re-evaluates `SingularityFunction` atoms in its result. A
  definite integral such as `integrate(SF(x, 0, 1), (x, 2, x))` therefore
  ends as `SF(x, 0, 2)/2 - 2`, not with an unevaluated `SF(2, 0, 2)`.
- `Piecewise(...).rewrite(Piecewise)` also rewrites `SingularityFunction`
  inside the pieces.
- `solve` ignores the `domain` flag, as the reference does (`Beam.max_*`
  passes it). A single univariate `Piecewise` equation that the native
  solver rejects falls back to the reference algorithm. Each piece is solved
  under its condition, with the earlier conditions negated. Undecided
  conditions give `Piecewise((sol, cond), (nan, True))`. Any other
  unsupported input still raises `NotImplementedError`.
- Matrix stacking skips 0x0 blocks.
- `Basic.is_same` exists. It compares mutable and immutable matrices as
  different classes, as the reference does.
- `factor`, `together` and `apart` map over lists and tuples, as the
  reference does.
- `vlatex` and `mlatex` print nested dynamic symbols and their time
  derivatives as the reference does (`q_{1}`, `\dot{q}_{1}^{2}`), not only
  at the top level, and print matrices. The general `latex(Matrix)` still
  raises `TypeError` on main; that gap belongs to the printing package.
- `plot` accepts a degenerate range such as `(x, 1, 1)`. The matplotlib
  backend draws the `markers`, `annotations`, `rectangles` and `fill`
  decorations used by Beam, Cable, Truss and Arch drawings.

`Cable` stays native, because the reference point-load and distributed-load
equilibria are wrong (see above). Two methods were added to it:

- `change_support`: the reference compares the first character of each load
  label with the support coordinates, which raises `TypeError` as soon as a
  point load exists. symox checks that every load position stays strictly
  between the supports.
- `draw`: it builds the reference plot with the same curves, support
  rectangles, load annotations and limits. The parabola base line can differ
  from the reference in the last printed digit (`-7.49552913752914` against
  `...915`), because the lowest point comes from the symox equilibrium.

Remaining differences:

- Plot objects: the reprs of `Beam.plot_*`, `Beam3D.plot_*` and `draw` differ
  from the reference `MatplotlibBackend` repr, because the symox plotting
  module has its own `Plot`.
- Form only:
  - `Beam3D` doctest: `collect(simplify(dy), x)` is mathematically equal to
    the reference expression (their difference simplifies to 0), but its
    structure differs.
  - `Beam3D.max_deflection` prints `-(-10 + a)**3` where the reference prints
    `(10 - a)**3`.
  - `singularityintegrate(x*SF(x, 0, 1), x)` gives `SF(x, 0, 3)/3`, where the
    reference gives `x**3*SF(x, 0, 0)/3`. The two are equal.
- Set and term ordering of `args_cnc`, `atoms` and `find_dynamicsymbols`
  outputs.
- Series: the `aseries` order term and the `PoleError` message of `nseries`
  on frame coordinates.
- `ObstacleSetPathway`: a `sqrt` form differs.
- `Beam.max_*` with symbolic `E` and `I` on the reference doctest beam
  exceeds 150 s in both libraries.

## Wave 4: limit_seq follows the reference algorithm

`limit_seq` is now a port of `sympy.series.limitseq`: `dominant`, the
`trials` loop of `_limit_seq` (difference step on the expanded numerator
and denominator, `gammasimp`, then the dominant terms of both), the odd and
even `Dummy` indices for alternating signs and trigonometric terms,
accumulation bounds, the termwise sum and the `Abs` fallback. The input is
rewritten as in the reference (`fibonacci` in `GoldenRatio`, `factorial`
and `subfactorial` in `gamma`). `dominant` is exported from
`symox.series.limitseq`.

The reference `_limit_inf` is `Limit(...).doit(deep=False)`. Symox uses the
native limit, with these differences:

- An unevaluated native `Limit` counts as undecided (`None`), and the
  algorithm continues. The reference returns the unevaluated `Limit` in
  that case (for example `catalan(n + 1)/catalan(n)`, and `F(n)` for an
  undefined `F`). Symox returns `None` there.
- Terms with a Sum or Integral give `None`, as the reference gruntz raises
  on them. A Product is first evaluated with `doit()`. The reference does
  not terminate on `Product(1 + 1/k, (k, 1, n))/n`. Symox gives 1.
- The native gruntz engine does not know the asymptotics of Gamma-type
  functions, and it returned wrong values on them (`0` for
  `factorial(n)*exp(n)/(n**n*sqrt(n))`). It is not used on Gamma-type terms.
  `simplify(combsimp(...))` is tried first. Then an exact ratio test proves
  0, oo or -oo (`5**n/factorial(n)`, `n**n/factorial(n)`). Any other case
  gives `None`.
- The ratio test also decides the powers of non-real numbers, such as
  `n**7*(I/3)**n`.
- With an odd or even index, `(-b)**e` becomes `b**e` or `-b**e` for a
  positive number `b`. The reference Pow constructor does this when the
  exponent is the index itself, and its gruntz handles the other parity
  exponents.
- `dominant` calls are nested at most 4 deep. The reference recurses until
  gruntz decides each ratio. Symox leaves more ratios undecided and would
  otherwise recurse until a RecursionError.
- When the odd subsequence has a limit but the even one does not, the
  result is `None`. In that case the reference raises an AttributeError.
- Plain Python numbers are sympified (the reference raises an
  AttributeError on `limit_seq(5, n)`).

Divergence kept on purpose: `lucas` is rewritten in its Binet form. The
reference does not rewrite `lucas`. It returns an unevaluated `Limit` for
`lucas(n + 1)/lucas(n)`, and a wrong 0 for `lucas(n)/fibonacci(n)`. Symox
gives `GoldenRatio` and `2*GoldenRatio - 1` (that is, `sqrt(5)`).

Probe against SymPy 1.14 (the cases of `test_limitseq.py`, the docstrings
and 120 more sequences): 131 of 159 match, up from 111 with the previous
implementation. The nested example
`limit_seq(Sum(k**2*Sum(2**m/m, (m, 1, k)), (k, 1, n))/(2**n*n), n)`
gives 4, and
`Sum(binomial(3*k, k)*binomial(5*k, k), (k, 1, n))/(binomial(3*n, n)*binomial(5*n, n))`
gives 84375/83351. Both results are the same in 20 fresh processes with
different `PYTHONHASHSEED` values. Each remaining difference is a `None` where the
reference finds a value. None of them is a wrong value.

Remaining gaps (all `None` in symox):

- `harmonic` at infinity: the native limit does not expand
  `harmonic(n)`, `harmonic(n, s)` or `polygamma` at `oo`. This affects
  `harmonic(n)/log(n)`, `harmonic(n) - log(n)`, `harmonic(2*n) - harmonic(n)`,
  `Sum(1/k, (k, 1, n))` and `Sum(1/k**2, (k, 1, n))` (that is, `harmonic(n, 2)`),
  `Sum(1/x, (x, 1, y)) - log(y)`, and the reference tests
  `(harmonic(n)*Sum(harmonic(k), (k, 1, n)))/(n*harmonic(2*n)**2)` and
  `Sum(harmonic(k)**2/k, (k, 1, 2*n))/harmonic(n)**3`.
- Gamma asymptotics (Stirling): `factorial(n)*exp(n)/(n**n*sqrt(n))`,
  `binomial(2*n, n)*sqrt(n)/4**n`, `binomial(2*n, n)/4**n`,
  `catalan(n)*n**(3/2)/4**n`, `subfactorial(n)/factorial(n)` (`exp(-1)`),
  `log(factorial(n))/(n*log(n))`, `1/gamma(n + 1/2)` (`gammasimp` does not
  reduce `gamma(n + 1/2)/gamma(n + 3/2)`), and `Sum(1/factorial(k), (k, 0, n))`.
- `floor(n/2)/n`, and the closed form of `Sum(m**k, (m, 1, n))` over
  `n**(k + 1)` with a symbolic `k`.
- The native limit crashes the process (segmentation fault, deep native
  recursion) on some large `Abs` quotients of powers `(n + a)**(n + b)`.
  This was found in the ratio test, which is now not used for such terms.

## Core gaps: Tuple equality, complex assumptions and GF(n) (wave 4)

Implemented:

- `Tuple` follows the reference `Tuple.__eq__`/`__hash__`: it equals the
  Python `tuple` of its items (item by item, so `Tuple(1.0, 2) != (1, 2)`
  as in the reference) and hashes as `hash(self.args)`, so a `Tuple` and its
  tuple share dict and set slots. The native hash uses CPython's tuple
  algorithm over the item hashes. `Tuple.__lt__`/`__le__` compare the args,
  and `index`, `count`, `tuple_count`, `in`, `+` and `*` with Python tuples
  and integers follow the reference. Expr equality is unchanged
  (Expr==Expr about 127 ns before and after; `hash(Tuple(1, x))` goes from
  90 ns to 105 ns; `Tuple(1, x) == (1, x)` from 2.5 us, unequal, to 0.8 us).
- The assumption engine ports the reference `Add`, `Mul`, `Pow` and `exp`
  handlers for `complex`, `algebraic`, `extended_real`, `imaginary`,
  `hermitian` and `antihermitian` (`_fuzzy_group`, `Mul._eval_real_imag`,
  `Mul._eval_herm_antiherm`, `Add._eval_is_imaginary`,
  `exp._eval_is_extended_real`, Lindemann for `exp._eval_is_algebraic`),
  plus `Add._eval_is_finite`, `Mul._eval_is_infinite` and
  `exp._eval_is_finite`. `(2 + I).is_complex`, `(pi*I).is_transcendental`,
  `(2*I).is_antihermitian`, `(oo + I).is_infinite` and `exp(I).is_finite`
  are True. Facts that imply a real value (`real`, `zero`, `rational`,
  `integer`, the signs) are False for numbers these handlers prove
  non-real, and for any expression whose native `is_real` is False
  (`(r + I).is_zero` is False). `(oo + I).is_extended_positive` was True and
  is now False. `Quaternion` no longer carries its local `is_complex`
  rule.
- `GF(n)` for composite `n`, `n` above `2**32` and `symmetric=False` is the
  reference `Z/nZ` domain (`python/symox/_modular_ring.py`): `is_Field` is
  `isprime(n)`, elements are `SymmetricModularIntegerMod{n}` /
  `ModularIntegerMod{n}`, non-units raise `NotInvertible('zero divisor')`,
  `get_ring` raises `DomainError`, and `GF(0)`/`GF(-4)` raise the reference
  `ValueError`. Symmetric prime moduli below `2**32` keep the native field.
  `lfsr_sequence` over `GF(4)` works as in the reference.

Remaining gaps:

- `is_square` and `sqrt` on `GF(9)` raise `NotImplementedError`
  (polynomials over `GF(n)`: see the section on polynomials over `Z/nZ`).
- `sqrt(2 + I).is_extended_real` is unknown (the reference proves it False
  through `arg(2 + I)/(2*pi)` not being an integer), and
  `(-2)**(1/3)` is not proven non-imaginary. A few products of non-real
  factors such as `(2 + I)*(3 - I)` are decided (False) where the
  reference is unknown.
- Ordering a `Tuple` against a plain Python tuple raises an error of
  another type than the reference. A `Tuple` cannot hold `Permutation` or
  matrix objects (native conversion), as before.

## Complex `solveset` of exponential and trigonometric equations on 2026-10-22

The native scalar `solveset` over `Complexes` returned incomplete finite sets
for transcendental equations, for example `{0, I*pi}` for `exp(x) - exp(-x)`.
A probe of 78 equations (exp, log, rational bases, trig, hyperbolic,
parameters, products and quotients, in `Complexes`, `Reals` and intervals)
went from 40 to 61 outputs identical to SymPy 1.14. The others are the same
sets in another form, or the divergences listed below.

What was implemented:

- `_solve_as_poly` for exponential generators. The equation becomes a
  polynomial in `t = exp(g*x)`, where `g` is the gcd of the exponents. Each
  root `k` is inverted as in `invert_complex`, to
  `ImageSet(Lambda(n, (I*(2*n*pi + arg(k)) + log(Abs(k)) - c)/g), Integers)`.
  The constant offset `c` of the generator is kept, so `exp(x + 1) - 2` gives
  `2*n*I*pi - 1 + log(2)`. The lattice is canonicalized as in
  `_set_function(Lambda, Integers)`: the sign of `n` is chosen with the fewest
  negatives, and the offset is reduced modulo the period.
- `_solve_as_rational`: the solutions of the numerator minus the roots of the
  denominator.
- `_solve_trig1`: trig and hyperbolic equations are rewritten in `exp` and
  solved for `exp(I*x)` or `exp(x)`. The poles are removed.
- Products are split factor by factor.
- When a step is undecided (no complete roots, a symbolic multiple of `x` in a
  trig equation, nested exponentials, `x**x`, mixed bases), the result is a
  `ConditionSet`. It is never the finite roots of `solve`.
- `nonlinsolve` and `substitution` inherit the results.
  `[exp(x) - sin(y), y**2 - 4]` now matches SymPy.
  `[2**x - 8, y - x]` and `[exp(y) - x, x**2 + 1]` are complete.

Intentional divergences:

- `b**x` with a positive rational base is solved through `exp(x*log(b))`.
  `2**x - 8` gives `2*n*I*pi/log(2) + 3`, where SymPy returns a `ConditionSet`.
- `2**x - 3**x` is a `ConditionSet`. SymPy returns the incomplete `{0}`.
- `exp(2*x) - exp(x)` has no spurious `nan`.
- `log(x) + log(x - 1)` keeps only `1/2 + sqrt(5)/2`. SymPy also keeps
  `1/2 - sqrt(5)/2`, which is not a solution.
- Some unions are printed in another order (`sin(x)*cos(x)`) or with another
  representative (`cos(x) - 2`). Another case is
  `exp(x)/(exp(x) - 2) - 3`, which gives a plain ImageSet where SymPy gives a
  `Complement`. These are the same sets.
- In `nonlinsolve([2**x - 8, y - x])` the real solution `(3, 3)` is listed
  next to the lattice that contains it.

What remains:

- `exp(exp(x)) - 1` is a `ConditionSet`. SymPy gives a `ProductSet` ImageSet.
- Trig equations with a symbolic coefficient of `x` stay a `ConditionSet`.
- `2**x*3**x - 6` (mixed bases) is a `ConditionSet`. SymPy also gives one.
- `solve(sin(x) + cos(x))` raises `NotImplementedError`. This predates the
  change; SymPy gives `[-pi/4]`.

Tests: `tests/test_solveset_complex_wave4_parity.py`.

## Canonical order and equality of Add and Mul over host classes

Root cause of `x + AccumBounds(1, 2) != Add(AccumBounds(1, 2), x)` (seen
after an import-order change): `AccumBounds.__radd__` returns
`Add(self, other, evaluate=False)`, a held node, while `Add(B, x)` was
canonical with the public arguments `(x, AccumBounds(1, 2))`. Native equality
compared the held storage with the canonical storage, which follows the node
hash order; the hash of a host node includes its class identity, so the
result depended on the import and registration order. The same accident made
209 of the 870 ordered pairs of plain symbols compare
`Add(a_j, a_i, evaluate=False) == a_i + a_j` as True where SymPy says False.

Fixes:

- A held Add or Mul equals a canonical one only when their public arguments
  agree in order (SymPy compares `args`). Equality is only stricter, so equal
  nodes keep equal hashes; canonical-canonical comparisons are unchanged.
- As in SymPy `Add.flatten` and `Mul.flatten`, AccumBounds is the
  coefficient: it leads the arguments of a canonical Add or Mul, and Numbers
  and other bounds fold into it through the bounds arithmetic
  (`Add(B(1, 2), 1, x).args == (B(2, 3), x)`, `(x + B(1, 2) + 1)`,
  `2*(x*B(1, 2))`, `(B(1, 2) + x)*2`). The folding runs only once bounds have
  been built, and in the native operators only for nodes holding host terms.

Remaining differences:

- Folding happens in the `Add`/`Mul` constructors and the arithmetic
  operators; other native rebuilds (`subs`, `expand`, `xreplace`) do not fold
  a Number into bounds.
- `2*(B(1, 2)*x + y)` gives `x*AccumBounds(2, 4) + 2*y`; SymPy keeps the
  distributed `2*x*AccumBounds(1, 2)` unfolded. The values are equal.
- Not native Add/Mul nodes, so outside this fix: `x + Ket('k')` is a quantum
  `OperatorExpr` whose `args` differ from SymPy and which differs from
  `Add(Ket('k'), x)` (an `OrderedProduct` wrapper); `codegen.ast.Indexed`
  sorts before symbols (`(x + A[x]).args == (A[x], x)`, SymPy `(x, A[x])`).
  Both construction orders agree for `Indexed`.

## Polynomials over Z/nZ on 2026-10-22

`Poly` with `modulus=n` or `domain=GF(n)` for composite `n`, `n` above
`2**32` or `symmetric=False` is `python/symox/_modular_poly.py`
(`ModularPoly`, reported as a `Poly`; `GF(n)` is reported as a `Domain`).
It runs the reference dense algorithms (`densearith`, `densetools`,
`euclidtools`, `sqfreetools`, `galoistools`, `factortools`, ported in
`python/symox/polys/_dense/`) over the reference `Z/nZ` ring, so results and
errors follow the reference: ring division (`div`, `rem`, `gcd` through the
PRS) and its `DomainError('there is no ring associated with GF(9)')`,
`NotInvertible('zero divisor')`, `ExactQuotientFailed`, `UnificationFailed`
with `QQ` or symbolic points, the `OptionError` for `modulus` below 1, the
reference printing (`Poly(x**2 - 4*x - 2, x, modulus=9)`). Arithmetic,
`div`/`rem`/`quo`/`exquo`/`prem`, `gcd`/`lcm`/`cofactors`/`gcdex`/`invert`,
`resultant`/`subresultants`/`discriminant`, `eval`, `degree`, `LC`, `monic`,
`diff`/`integrate`, `compose`/`decompose`/`shift`, `sqf_*`, `factor_list`
and the expression functions (`factor`, `sqf`, `gcd`, `rem`, ...) are
covered. Symmetric prime moduli below `2**32` keep the native path. Tests:
`tests/test_composite_modulus_parity.py`.

Divergence: the reference equal-degree factorization is randomized, so over a
composite modulus `factor_list(x**2 - 1, modulus=9)` returns factors or
raises `NotInvertible` from run to run; symox seeds it (seed 0) and is
reproducible.

Remaining gaps:

- `groebner` over a composite modulus raises the reference `TypeError` only
  when a single non-monic generator decides it (lex order); other inputs
  raise `NotImplementedError` (the reference runs Buchberger/F5B over the
  ring and returns a basis or raises `TypeError`). `groebner` over a prime
  modulus above `2**32` or with `symmetric=False` raises
  `NotImplementedError`.
- `roots`, `real_roots`, `nroots` and `all_roots` of such a polynomial raise
  `NotImplementedError` (the reference returns the integer roots of the
  symmetric representative). `solve_poly_system` raises another error.
- `rep` and `per(rep)` (the reference `DMP_Python`), `LM`/`EM`/`LT`/`ET`
  (`Monomial` objects), `eject`, `nth_power_roots_poly` and other methods
  not listed above raise `NotImplementedError`. `srepr` shows the modulus
  form and `pretty()` raises `TypeError`.

## Limits at oo of gamma-like functions (Stirling, wave 4)

Root cause of the wrong `0` for `limit(factorial(n)*exp(n)/(n**n*sqrt(n)), n, oo)`
(reference `sqrt(2)*sqrt(pi)`): gruntz took the mrv set of the arguments of
`gamma`/`factorial` and so treated them as being in the class of `n`, and the
series engine had no expansion of `loggamma` at an infinite argument.

Changes:
- gruntz `prepare` applies the reference tractable rewrites: `gamma(a)` ->
  `exp(loggamma(a))` when `a -> +oo`; `factorial`, `binomial`, `rf`, `ff` ->
  gamma ratios; `harmonic` -> `polygamma`; `fibonacci` -> powers of the
  golden ratio. `loggamma` near a non-finite, non-`+oo` limit is unknown.
- mrv refuses (unknown, so `Limit` is returned) a function that grows or
  decays beyond its argument's class (`gamma`, `Ei`, `erfi`, `besseli`,
  `lucas`, `catalan`, ...) when no tractable rewrite applies, instead of
  returning a wrong value. Before: `Ei(n)*exp(-n)*n`, `erfi(n)*exp(-n)`,
  `fibonacci(n)/GoldenRatio**n`, `besseli(0, n)*exp(-n)*sqrt(n)` gave `0`.
- The series engine expands `loggamma` (Stirling series with Bernoulli
  terms) and `polygamma(k, z)`/`digamma`/`trigamma` (A&S 6.3.18, 6.4.10) at
  an argument going to `+oo`; other functions of the gamma family there are
  still a pole.
- The mrv leading term writes `x` as `log(w)/k` when `log(w) = k*x` (and
  `exp(alpha*log(w) + beta)` as `w**alpha*exp(beta)`), like the reference
  `logx` expansion, so cancellations between `x` and `log(w)` are seen; a
  coefficient that vanishes only after the `log(w)` substitution is
  subtracted and the next term taken. A product of exponentials that the
  canonical `Mul` re-splits after the rewrite (no progress) is merged into
  one exponential before recursing.

Matches the reference: the Stirling forms, `gamma(n + 1/2)/(gamma(n)*sqrt(n))`,
`gamma(n)/gamma(n + 1/3)*n**(1/3)`, `factorial(n)**(1/n)/n`,
`binomial(n, 2)/n**2`, `factorial(n + 1)/factorial(n) - n`,
`loggamma(n) - (n - 1/2)*log(n) + n` (`log(2)/2 + log(pi)/2`).

Divergences, symox correct: `rf(n, n)*factorial(n)/factorial(2*n)` is `1/2`
(reference `oo`, wrong); `ff(2*n, n)*factorial(n)/factorial(2*n)` and
`ff(n, n)/factorial(n)` are `1` (reference unevaluated or raises);
`gamma(n + a)/(gamma(n)*n**a)` is `1` (reference raises). The reference
`subfactorial(n)/factorial(n) -> 0` is wrong (true value `exp(-1)`); symox
returns `Limit`.

Remaining gaps (symox returns `Limit`, never a number): `binomial(2*n, n)`
and `factorial(2*n)/factorial(n)**2` forms with `4**n`, `harmonic`,
`polygamma` and `loggamma` differences, `uppergamma`, `erfc*exp`, `Ei`, `li`,
`Si`, `Ci`, `zeta`, `lucas`, `catalan`, `Shi`, `Chi`, `factorial2`,
`subfactorial`, `polylog(2, 1/n)*n`, `gamma(log(n))`, `loggamma(1/n)`,
`gamma(-n)`, `harmonic(n, 1/2)`; `limit((x + 1)**(x + 1)/x**x/x, x, oo)`
(reference `E`). `loggamma(1/x).series(x, 0, 2)` prints `O(x**2*log(x))`
where the reference prints `O(x**2)`.

Native segfault on large `Abs` power quotients (`(n + a)**(n + b)`) from the
`limit_seq` ratio test: not reproduced on this base (192 `Abs` quotient
variants and the `limit_seq` paths run in subprocesses); a subprocess
regression test guards it. `limit_seq((n**2 + a)**(n + b)/(n + a)**(2*n), n)`
returns `0` where the reference returns `None` (the value depends on `b`):
an open wrong value in the `limit_seq` Python layer, not in the native
`limit` (which stays unevaluated).

## Wrong results: multiple Piecewise integrals, Range bounds, Indexed atoms, exact lambdify (wave 4, final A)

Fixes (regressions in `tests/test_wrong_final_a_parity.py`):

- Multiple integrals. `piecewise_fold` no longer lifts a Piecewise condition
  that holds a bound variable out of its Integral, Sum or Product: the
  folded `Piecewise((Integral(x**2, (t, 4/5, 1)), x/t < 0), ...)` left `t`
  free outside the t-integral. When an inner definite limit stays
  unevaluated and more limits follow, its Piecewise conditions are reduced
  by the sign the limits give the variable (SymPy's signed dummy of
  `Integral.doit`: `x/t < 0` is `x < 0` for `t` in `[4/5, 1]`) and the
  limit is tried again, so
  `(x**2 + Piecewise((0, x/t < 0), (1, True))).integrate((t, 4/5, 1), (x, -1, 1))`
  is `1/3` as in SymPy. A single inner limit stays unevaluated as in SymPy.
- `Range._inf`, `_sup`, `inf`, `sup`, `is_empty`, `is_finite_set` and
  `__bool__` are SymPy's: a symbolic Range whose emptiness or direction is
  undecided raises `ValueError` (`Range(ii, jj, kp).inf` returned `jj - kp`,
  and `sup` `ii`), an empty Range raises `NotImplementedError` like
  `EmptySet.inf`. The native set infimum and supremum no longer answer for a
  symbolic or empty Range. The infimum and supremum of a Union holding such a
  Range fall back to SymPy's `Min`/`Max` of its parts:
  `Union(Range(ii, ii + 5*kp, kp), Interval(100, 200)).inf` is
  `Min(100, ii)`, and an undecided Range part raises its `ValueError`.
- `expr.atoms(Indexed)` lists the Indexed references of the expression, so
  `solve(eqs, eqs.atoms(Indexed))` is `{A[1]: 1, A[2]: 2}` instead of `[]`.
- `lambdify` of a polynomial with Integer coefficients and nonnegative
  integer powers evaluates `int` arguments with Python integer arithmetic,
  as SymPy's generated code does: `lambdify(x, x + 1)(2**60)` is
  `1152921504606846977` for the default, `math`, `mpmath`, `sympy` and
  `numpy` modules.

Intentional divergence:

- SymPy folds a condition on a bound variable out of its integral and can
  return a result where that variable is free
  (`integrate(Piecewise((t, x/t < 0), (1, True)), (t, -1, 1), (x, 0, 1))` is
  `Integral(Piecewise((0, x/t < 0), (2, True)), (x, 0, 1))`, true value
  `1/2`); symox keeps the double integral unevaluated.

Remaining gaps:

- `Indexed`, `IndexedBase` and `Idx` are still the `codegen.ast` dataclass
  references (and a plain `Idx` class), not SymPy's `sympy.tensor.indexed`
  Expr nodes: in arithmetic `A[1]` becomes the symbol named `A[1]`, so
  `solve` returns those symbols as keys (equal to `A[1]`, which hashes like
  its symbol so `sol[A[1]]` works), `isinstance(atom, Indexed)` is False
  for the atoms found, `fcode(A[i])` is
  `A[i]` where SymPy prints `A(i)`, and `ccode(A[1] + A[2])` prints instead
  of raising `IndexException`. Making them Expr nodes touches every code
  printer and the codegen contraction support.
- Other lambdify differences: the native evaluator returns Python `float`
  where SymPy's numpy namespace returns `numpy.float64` (equal values), and
  the compiled function has no Python source for `inspect.getsource`.
  Non-polynomial expressions with `int` arguments are evaluated in doubles.

## Dagger, random matrix products, nested Expectation, rational nonlinsolve

Wrong results fixed against SymPy 1.14:

- `Dagger(A)` of a noncommutative symbol returned `conjugate(A)` and kept the
  factor order of products. It now returns the held adjoint, printed
  `Dagger(A)` as in SymPy (whose `dagger` module renames `adjoint`), and
  `Dagger(A*B) == Dagger(B)*Dagger(A)`. Add, Pow, numbers, commutative symbols
  and matrices follow the reference. `TensorProduct.expand(tensorproduct=True)`
  distributes over Add factors.
- `A*X` with `A = MatrixSymbol(k, k)` and `X = RandomMatrixSymbol(k, 1)` was
  reordered as a scalar product. Native matrix-expression products now accept
  host objects that expose a native matrix expression (`_symox_matrixexpr`),
  and scalar products defer to them, so `A*X` stays `A*X` with shape `(k, 1)`
  and `a*X` stays a matrix expression.
- `ExpectationMatrix`, `VarianceMatrix` and `CrossCovarianceMatrix` are matrix
  expressions with the reference shapes, and their `expand` is a port of the
  reference (linearity, constant matrix factors, `A*Var(X)*A.T`, cross terms).
  `VarianceMatrix`/`CrossCovarianceMatrix` of a symbolic non-vector raise
  ShapeError. As before, matrix-distribution variables (`MatrixNormal`, ...)
  are still accepted and evaluated, an extension over SymPy, which raises
  ShapeError for them.
- `Expectation(...).doit()` did not evaluate nested Expectations; it now
  follows the reference: deep doit evaluates inner expectations
  (`Expectation(X + Expectation(Y + Expectation(2*X))).doit() == 3*mu + 1`),
  and `deep=False` keeps held inner terms
  (`mu + Expectation(Expectation(Y))`).
- `nonlinsolve([x**2 + 2/y - 2, x + y - 3], [x, y])` lost solutions. A
  rational polynomial system whose numerators form a zero-dimensional ideal is
  now solved like the reference (grevlex Groebner basis, FGLM to lex,
  `solve_poly_system`).

Intentional divergence:

- Roots that make a denominator zero are dropped. SymPy returns the spurious
  `(0, 0)` for `nonlinsolve([1/x + 1/y - 1, x - 2*y], [x, y])`; symox returns
  `{(3, 3/2)}`. When it cannot be decided whether a denominator vanishes, the
  numerator path is not used and the previous path runs.

Remaining gaps:

- `srepr` and pretty printing of `Dagger(A)` still show `adjoint` (native repr
  is not routed to host printing hooks); `latex(Dagger(A)**2)` prints
  `\operatorname{adjoint}^{2}` instead of `A^{\dagger 2}`.
- Type names differ: `TensorProductAtom` vs `TensorProduct`, `_MatrixExpr` vs
  `MatMul`/`MatAdd`, `Matrix` vs `MutableDenseMatrix`. The deprecated
  `tensor_product_simp(e**2)` returns a value where SymPy returns None.
- `CrossCovarianceMatrix(a*X, b*Y).expand()` prints `a*b*CrossCovarianceMatrix(X, Y)`;
  SymPy keeps an unflattened scalar and prints `(a*b)*CrossCovarianceMatrix(X, Y)`.
- `ExpectationMatrix((X + Y)*(X - Y).T).expand()` is equal to the reference
  but its MatAdd terms print in a different order.
- Printing of held Expectation nodes shows internal dummy names (`_Y`) and
  native factor order inside the held argument; the values are equal.
- `nonlinsolve([x**3 + 1/y - 2, x - y], [x, y])` and similar systems whose
  numerators are not handled by the new path still hang in the native
  parametric solver (already the case before this change).
- Not from this change, present on the base commit: `sqrt(-sqrt(2))` stays
  unevaluated (SymPy `2**(1/4)*I`), so `dsolve(f(x).diff(x, 3) +
  sqrt(2)*f(x).diff(x))` returns a `cosh(x*sqrt(-sqrt(2)))` term instead of
  `cos(2**(1/4)*x)` (ODE corpus cases `lin_const_coeff_hom_38`, `reducible_03`).

## Trig equations in `solve`, symbolic trig coefficients in `solveset` on 2026-10-22

A probe of 77 calls (sin/cos/tan combinations, shifted and multiple
angles, powers, hyperbolic, parameters, `dict=True`, `check=False`, and
`solveset` with a symbolic coefficient) went from 28 to 65 outputs identical
to SymPy 1.14.

What was implemented (`crates/symox-core/src/solve/trig.rs`):

- The angle expansion of the reference `solve`: when more than one distinct
  trig atom depends on the unknown, `TR10` expands `sin(x + c)` and `TR11`
  doubles the angles equal to twice the gcd of the arguments.
- The trig branch of the multiple-generator case of `_solve`. Two terms
  `A*sin(u)**p + B*cos(u)**p` become the numerator of `TR2i(A/B) + 1`
  (`sin(x) + cos(x)` gives `[-pi/4]`). Anything else is solved as
  `TR1(f).rewrite(tan)`, a rational function of `tan(x/2)`, then checked
  against the original equation. When the rewritten form is unsolved, the
  result is `NotImplementedError`, as in the reference.
- The `rewrite(exp)` fallback of `_tsolve` for trig equations whose unknown
  occurs only in trig functions of `r*x + c` with rational `r` (other
  equations, such as `y*sin(y/x)` in `odesimp`, are left unsolved without
  the rewrite, which only made them slow). A sum of `exp(r_i*I*x)` with
  rational `r_i` is solved as a polynomial in `exp(I*x/q)`; the lattice is adopted only when every
  logarithm of a root evaluates (`log(-1)`, `log(I)`), so
  `tan(x) - sin(x)` gives `[0, -pi, pi, 2*pi]` and `cos(x) + cos(3*x)`
  gives `[-3*pi/4, -pi/2, -pi/4, pi/4, pi/2, 3*pi/4]`.
- `sec(u) = r` and `csc(u) = r` are inverted through `acos(1/r)` and
  `asin(1/r)` when these evaluate to multiples of `pi` (`sec(x) - 2` gives
  `[-pi/3, pi/3]`, a double root such as `sec(x) + 1` gives `[pi]`);
  otherwise the equation is left to the other paths.
- `2*atan(r)` from the tan half-angle rewrite is evaluated at
  `radsimp(r)` when that gives a tabulated angle
  (`sin(x) + cos(x) - sqrt(2)` gives `[pi/4]`).
- `solveset` of a single trig or hyperbolic function of `a*x + b` with a
  symbolic `a` (the reference `_solve_trig`): the unit lattice is inverted
  without expansion, `ImageSet(Lambda(n, (2*n*pi + pi/2 - b)/a), Integers)`,
  inside `ConditionSet(x, Ne(numer(a), 0) & Ne(denom(a), 0), ...)`. This
  applies only when the function occurs once, linearly; otherwise the result
  stays a `ConditionSet`, as in the reference.

Intentional divergences:

- Candidates containing `oo`, `zoo` or `nan` are dropped, also with
  `check=False`. `solve(tan(x) - sin(x), check=False)` gives
  `[0, -3*pi/2, -pi, -pi/2, pi/2, pi, 3*pi/2, 2*pi]`, the unchecked roots of
  the exp polynomial, where SymPy gives `[0, zoo, -pi, pi, 2*pi]`. With the
  default check both give `[0, -pi, pi, 2*pi]`.
- `solve(sin(x/2) - cos(x))` gives `[-pi, pi/3, 5*pi/3]` (each one checked).
  SymPy raises `NotImplementedError`.
- `solve(cos(x)**2 - sin(x))` prints `sqrt(2 - 2*sqrt(5))/2` where SymPy
  prints `sqrt(2)*sqrt(1 - sqrt(5))/2`. These are the same values.

What remains:

- `solve(sin(3*x) - sin(x))` gives `[0, -3*pi/4, -pi/4, pi/4, 3*pi/4, pi]`.
  SymPy gives the same list without `pi`, which is also a solution: its
  selection of principal values along the tan and exp rewrites is not
  reproduced.
- `sec` and `csc` equations whose inverse does not evaluate to a multiple of
  `pi` (`sec(x) - 3`) raise `NotImplementedError`; SymPy gives logarithms
  such as `I*(log(3) - log(1 - 2*sqrt(2)*I))`.
- `solveset(exp(exp(x)) - 1, x)` and other nested exponentials stay a
  `ConditionSet`. SymPy inverts them over a `ProductSet(Integers, Integers)`
  ImageSet (`_invert_complex` on an ImageSet right-hand side).
- With a symbolic coefficient, `solveset` of several trig functions
  (`sin(a*x) + cos(a*x)`, the reference `_solve_trig1` with a symbolic
  `mu`) and of a symbolic right-hand side (`cos(a*x) - b`) stays a
  `ConditionSet`.

Tests: `tests/test_solve_trig_rest_parity.py`.

## Missing domains and features, batch c1 (integrals, sums, series, calculus, holonomic, stats; wave 4)

Implemented (regressions in `tests/test_missing_c1_parity.py`, values from
SymPy 1.14):

- `Integral.as_sum`, `Integral.transform` and `Integral.principal_value` are
  ports of the reference methods (`python/symox/integrals/_integral_extras.py`);
  their docstring examples match.
- `Order.contains` and `in` are the reference algorithm
  (`python/symox/series/_order_contains.py`): `O(1) in O(1, x)` is True, an
  undecided containment raises `TypeError` for `in` like the reference.
- A multivariate `Order` of a sum keeps its leading terms only (the
  reference `extract_leading_order` on the expanded expression, every point
  moved to the origin): `O(1 + x*y)` is `O(1, x, y)`, `O(x + y + x*y)` is
  `O(x + y, x, y)`, `O(1 + x*y, (x, oo), (y, oo))` is
  `O(x*y, (x, oo), (y, oo))`. A term whose ratio to another is not a
  monomial is kept.
- The reference `eval_sum_hyper` stage runs after the native summation on a
  single-limit `Sum` with an infinite bound that the native engine leaves
  held (`python/symox/concrete/_sum_hyper.py`): `Sum(x**k/factorial(k),
  (k, 0, oo))` is `exp(x)`, `summation(x**k/k, (k, 1, oo))` is the
  reference `Piecewise` with its convergence condition. `Sum(f, (i, oo, oo))`
  stays held.
- `limit` no longer returns `0` for a power of the variable with an
  exponent of unknown sign: gruntz treats the exponential of an argument
  whose limit is `oo*sign(b)` as infinite (reference mrv), and a leading
  coefficient that holds `w**s` with a symbolic `s` goes to the symbolic
  exponent algebra, which needs the sign of `s`. The mrv comparison drops
  the x-free factors of the ratio of logarithms, so `exp(x)/x**a` is still
  `oo` and `exp(-x)*exp(b*sqrt(x))` still `0`. `limit(x**b/x, x, oo)` and
  `limit(exp(b*x)*exp(-x), x, oo)` stay held (reference: `oo**(b - 1)` and
  NotImplementedError), and
  `limit_seq((n**2 + a)**(n + b)/(n + a)**(2*n), n)` is None like the
  reference instead of `0`.
- `periodicity` of a power whose base and exponent both vary is the common
  period of the two: `periodicity(sin(4*x)**cos(2*x), x)` is `pi`.
- `Idx` supports arithmetic, so `apply_finite_diff` on `x[i + j]` points
  works.
- Holonomic: differential operators print as `(1)*Dx**2`, a
  `HolonomicSequence` prints its initial values (`u(0) = 1`), the operator
  algebras print the reference text, and `to_sequence` at an ordinary point
  returns the reference `smallest_n` (`lower + degree`, raised past the
  integer roots of the leading coefficient) and `order + smallest_n`
  initial values.
- Stats: `density()` of a finite expression returns a `FiniteDensity` (so
  `.dict` works); `DiscreteRV(..., set=set(range(n)))` (a finite support of
  consecutive integers) is the equal `Range`; `DiscreteMarkovChain` has
  `communication_classes`, `decompose`, `canonical_form`,
  `number_of_states` and `index_of`; the `MultivariateNormal` density
  exponent is expanded like the reference.

Intentional divergences (equal values, other forms or stricter answers):

- `Integral(x**2, (x, 0, 1)).as_sum(n, 'left')` and
  `Integral(sqrt(x), (x, 0, 1)).as_sum(3, 'right')` come back in the native
  canonical form (`sqrt(3)/9 + sqrt(6)/9 + 1/3` for
  `sqrt(3)*(1 + sqrt(2) + sqrt(3))/9`).
- `apply_finite_diff` with indexed points prints its terms in another order.
- `cdf(Moyal(...))` simplifies to `erfc(...)` where the reference prints
  `1 - erf(...)`.
- `approximants` cancels its coefficients exactly, so the continued fraction
  of a rational series ends: for the binomial-row docstring example the
  generator stops after `-1/(t*(x + 1) - 1)` where the reference goes on to
  yield `nan` (an artefact of uncancelled coefficients).
- `Sum(1/x, (x, a, b)).euler_maclaurin()` keeps `Integral(1/x, (x, a, b))`
  where the reference writes `log(b) - log(a)`, which is wrong when
  `a < 0 < b`.

Remaining gaps:

- `fps(x**n*sin(x**2), x)`: the native series does not expand `x**n` with a
  symbolic `n`, and the native `Order` of `x**(n + 8)` lists `n` as a
  variable.
- `Product` docstring: `(1/4)**n*4**n` is not merged to `1` by the native
  `Mul`, and limits of `RisingFactorial`/`gamma` ratios with half-integer
  shifts (Wallis product) stay held.
- Holonomic: Frobenius initial conditions (a dict `y0`), `to_expr` of
  `besselj` (the initial values are not resolved), and `to_meijerg` of
  general annihilators raise NotImplementedError/ValueError.
- Stats: `BernoulliProcess.joint_distribution` and
  `JointDistributionHandmade` are missing (and discrete-time process values
  print as `B(1)`, not `B[1]`); `density(X).dict` of a `Binomial` or `Die`
  with a symbolic parameter (the reference held `Density` object);
  `ContinuousMarkovChain` with float generators (numeric matrix
  exponential) and Markov queries at symbolic times; `DiscreteMarkovChain`
  with `Str` states; `stationary_distribution` with a `MatrixSymbol`;
  multivariate `Laplace`; `MultivariateNormal` with `MatrixSymbol`
  parameters; `GaussianSymplecticEnsemble` density at a `MatrixSymbol`;
  Wiener and Gamma process queries under `Contains` conditions;
  `given(Tuple(X, Y), cond)`. Each raises NotImplementedError or stays held;
  `sympy.stats.Expectation` is handled separately.
- `limit(x**p*x**(-p - 1), x, oo)` with a positive `p` stays held
  (reference `0`).

## Missing callables, batch c2 (2026-10-23)

Ported from SymPy 1.14 and checked against it in
`tests/test_missing_c2_parity.py`: `symbols` with extra assumptions,
`Rational` from strings and `Rational.factors`, `entails`, the `motzkin`
helpers, `Set.kind`, `Set.powerset`, `ProductSet` iteration, `Predicate`
subclasses, `fu` names (`TR1`, ...), `sqrt_biquadratic_denest`,
`factor_system`, `ff(n, n)` for integer `n`, `partition(n)` for negative or
zero integer `n`, `loggamma(p/2)` for symbolic rational `p/2`, `Expint`
derivative in `nu`, the `nseries` of `sinc`, `elliptic_k`, `elliptic_e`,
`elliptic_f`, `elliptic_pi`, `airyai`, `airybi` and their derivatives,
`nsimplify(..., tolerance=..., full=True)` through `mpmath.identify`,
mpmath numbers in `sympify` at their own precision, the `frac`, `numer`,
`denom`, `modulus` and `exact` hints of `expand`, `refine` of a
`MatrixElement`, `opt_cse` (with `Unevaluated`), `diff(..., simplify=False)`,
the sum-of-even-powers solver of `diophantine`, and `1/(b)` printing of an
unevaluated power base. `harmonic`, `nonlinsolve` and `sysode` were already
correct.

Correct results that differ from SymPy:

- `series(sinc(x), x, 1, 3)` returns the Taylor expansion; SymPy returns
  `sin(x)/x` unexpanded.
- `series(elliptic_k(m), m, 1, 2)` raises NotImplementedError (logarithmic
  singularity); SymPy returns `nan`. Some SymPy RecursionErrors in special
  function series become NotImplementedError.
- `srepr(Symbol('x', Integer=True))` omits the extra assumption.

Remaining gaps:

- Multiple derivatives are not ended with `factor_terms(signsimp(...))`:
  `diff(sqrt((x + 1)**2 + x), (x, 5)).count_ops()` is 136 (SymPy 30); the
  values are equal. Applying it in `diff` broke the native undetermined
  coefficients ODE solver (10 corpus cases), so it was not kept.
- `(x**y).nseries(x, 0, 2)` raises NotImplementedError; SymPy returns `x**y`
  (the series engine cannot carry a symbolic power prefactor).
- `trigsimp(..., old=True, compare=True)` uses futrig instead of the legacy
  pattern-matching `_trigsimp` (`tan(x)**2` vs `cot(x)**(-2)`, and the
  `futrig:` comparison line is not printed).
- `diff(assoc_laguerre(n, a, x), a)` stays unevaluated on purpose: the
  reference `Sum(assoc_laguerre(_k, a, x)/(n - a), (_k, 0, n - 1))` is wrong
  (for n = 1 it gives `1/(1 - a)` instead of 1).
- `devise_plan_meijer(G_Function([0], [], [0], []), G_Function([1], [], [1], []), z)`
  raises "Singular contiguous shift"; SymPy gives
  `[MeijerUnShiftB, MeijerShiftA]`.
- `solve_triangulated(..., extension=True)` raises NotImplementedError: the
  extension towers over `QQ<CRootOf(...)>` are not ported (the polynomial
  constructor rejects a CRootOf coefficient whose generator symbol is also a
  polynomial generator).

## Powers of products with a negative coefficient (wave 4, pow-negbase)

`Pow(Mul, r)` for a non-integer rational or float `r` now ports SymPy's
`Mul._eval_power` -> `Pow._eval_expand_power_base` in the native constructor
(`expand_power_base_split` in `ops.rs`), for commutative and noncommutative
products. Factors that are polar or extended nonnegative (numeric sign proofs
included, e.g. `sin(1)`) get their own power, an even number of negative
factors becomes positive, a negative numeric coefficient sharing the residual
power with factors of unknown sign splits into `-1` and its absolute value,
and the residual (non-real factors, unknown signs, the odd `-1`, the
noncommutative product) keeps one power. Before, only non-numeric positive
factors were split, so `sqrt(-sqrt(2))` stayed unevaluated (now `2**(1/4)*I`),
as did `(-sqrt(2))**(1/3)` (now `(-1)**(1/3)*2**(1/6)`). `sqrt(n*m)` with
negative `n`, `m` is now `sqrt(-m)*sqrt(-n)`, and `sqrt(-2*x*A*B)` is
`sqrt(2)*sqrt(-x*A*B)`. `Pow._eval_is_extended_real` now also gives `False`
for a negative base to a non-integer rational power (`sqrt(n)`,
`sqrt(-sqrt(2) - 1)`). This also fixes the ODE regression
`dsolve(f' + sqrt(2)*f''', f)`, which now gives the reference
`C1 + C2*sin(2**(3/4)*x/2) + C3*cos(2**(3/4)*x/2)` deterministically.

Probes (101 bases times 17 exponents, plus facts): 1 232 of 1 669 matched
before, 1 588 match now. The 81 remaining differences are all numerically
equal (checked at sample points) and come from two form-only gaps:
- SymPy's `Mul.flatten` keeps nested products when it re-evaluates held
  powers, e.g. `(1.0*I)*p**0.5`, `(-n)*sqrt(x)`, `-I**(3/2)` from `I**(-1/2)`.
  Symox prints the flat product.
- Numeric radical products combine differently: SymPy gives
  `(-3)**(1/3)*2**(1/6)*sqrt(3)/3` for `(-sqrt(6)/3)**(1/3)`, symox gives
  `(-1)**(1/3)*3**(2/3)*6**(1/6)/3`. The cause is a pre-existing difference
  in how `Mul` combines numeric powers with different exponents, for example
  `6**(1/6)*(-1/3)**(1/3)`. It is not specific to `Pow`.

Remaining gap: `im`/`re` of a product holding an unevaluated
`Pow(-sqrt(2), 1/2, evaluate=False)` returns the correct but unsimplified
`-sqrt(2)*I*sqrt(-sqrt(2))/2`. SymPy re-evaluates held powers inside
`Mul.flatten` and gets `2**(3/4)/2`. Canonical inputs are not affected.

The split adds dependent numeric radicals to the bottom-up `normal` form that
the reference `simplify` state cancels, for example `12**(1/6)`, `12**(5/6)`,
`sqrt(3)` and `sqrt(229)` in the residual of a Ferrari root of
`x**4 + x + 1`. The native algebraic gcd behind that `cancel` does not reduce
the radical relations and blew up, so `solve_poly_system([x**4 + x + 1], x)`
stopped returning. That state now runs under the default factor budget. When
the budget runs out, simplify skips this one candidate and keeps the others
(trigsimp, the cancelled, expanded and together forms, powsimp).
`solve_poly_system` again returns in about 4 s, as before. Remaining gap: on
such inputs the native cancel is still slower than SymPy's, which treats the
radicals as independent generators.

## Matrices, polys and printing/codegen callables (c3) on 2026-10-22

The doctests of 41 callables were compared against SymPy 1.14. Before this
work, most of them had a few examples that did not match. Now 29 match on
every example. `MutableSparseMatrix` has no examples of its own; it is
checked through the `SparseRepMatrix` examples.

Implemented, printing and codegen:

- `Indexed` stays a flattened symbol in the expression tree. Making it a
  native node would need the native printers and core to change, so it is
  handled at the printer level:
  - The code printers print a flattened `Indexed` symbol through
    `_print_Indexed`, which gives the strided form (`A[i*n + j]`) for a
    shaped base.
  - `Idx` arithmetic (`y[i + 1]`) works.
  - `ccode`, `fcode`, `jscode`, `rcode`, `rust_code`, `glsl_code`,
    `julia_code` and `octave_code` with `assign_to` and ranged `Idx`
    generate the reference loops, including the summation loops of a
    contraction. This closes the earlier gaps on `Idx` arithmetic and
    contraction loops.
- `latex`:
  - tuples, lists, dicts and sets as containers;
  - matrices with `mat_str`, `mat_delim` and `decimal_separator`;
  - `long_frac_ratio`, in the native printer;
  - `emptyPrinter`.
- `codegen.rewriting` follows the reference:
  - `Optimization.cheapest`, and `optimize` uses it;
  - `exp2_opt`, `sinc_opts`, `log2_opt`, `log2const_opt` and
    `logsumexp_2terms_opt`;
  - `FuncMinusOneOptim`, used by `expm1_opt` and `cosm1_opt`;
  - `log1p_opt`;
  - `create_expand_pow_optimization`, with `UnevaluatedExpr`;
  - `matinv_opt`, which requires a full-rank `MatrixSymbol` (asked through
    `Q.fullrank`);
  - the `optims_c99`, `optims_numpy` and `optims_scipy` tuples.
- `SeriesApprox` is a port: it checks the parity of `n_point_checks`, and
  lowers the order while the relative tolerance holds.
- `MatrixSolve` has `matrix` and `vector`, and the numpy printer prints
  `numpy.linalg.solve`. `CodegenAST.has` looks inside native expression
  arguments, which makes `newtons_method` match.
- `utilities.codegen` is SymPy 1.14's module, with the imports changed:
  - matrix outputs are kept in an `_Exprs` container;
  - `get_datatype` covers every language;
  - `Argument` datatypes use `getattr` for `is_integer` and `is_real`.
- `autowrap`:
  - `ufuncify` passes `helpers` to the C ufunc wrapper as static functions;
  - with a non-numpy backend, it vectorizes the scalar wrapper.
- `smtlib_code`:
  - infers the types of symbols and of undefined functions;
  - gives the reference warnings for a type it cannot infer and for a
    non-Boolean expression;
  - raises TypeError for an unknown undefined function;
  - sorts its declarations;
  - prints `Integer` and `Rational` as SymPy does.
- `public` is the reference decorator.
- `plot_implicit` of an `And`, `Or`, `Not` or `Xor` of inequalities draws
  the filled region (`contourf`). The reference raises
  NotImplementedError when an equality is part of such a combination, and
  so does symox.

Implemented, matrices:

- `BlockMatrix.LUdecomposition` (BlockDiagMatrix).
- The generalized and off-diagonal Schur complements.
- `MatPow.inverse()` as `MatPow(base, -exp)`.
- `KroneckerProduct.has_matching_shape` and `structurally_equal`.
- `Matrix.has` with pattern classes.
- `M[i, j] = B` block assignment, for dense and sparse matrices.
- `SparseMatrix` creation following `_handle_creation_inputs`, with SymPy's
  messages.
- `func_name` of dense matrices.
- The product and integer-power rules of `combine_kronecker`.

Implemented, polys:

- `apart(full=True)`, ported in `polys.partfrac`.
- `assemble_partfrac_list` builds a `RootSum` for a non-linear factor.
- `RootSum.doit` inside larger expressions.
- `Poly(domain='R')` and the domain aliases.
- `Poly.per` with a DMP.
- `CRootOf.eval_approx` with exact zero parts.
- `which_all_roots`: it evaluates before taking `Abs`, which removes a
  timeout.
- `sfield`/`sring` infer non-Symbol generators, and ring and field elements
  print as in SymPy.

Intentional divergences:

- A flattened `Indexed` orders its strides and factors by symbol name.
  `A[i*n + j]` is printed where SymPy prints `A[n*i + j]`, and
  `ccode(A[i]*x + 1)` gives `A[i]*x + 1` where SymPy gives `x*A[i] + 1`.
- Generated code headers say "Code generated with SymPy" followed by the
  symox version.
- With a non-numpy backend, `ufuncify` returns a `numpy.vectorize` wrapper,
  so an array argument gives an array.
- With `helpers`, the numpy `ufuncify` returns the correct array. The
  reference raises a numpy TypeError in this environment.
- `logaddexp_opt`, `cosm1_opt` and `powm1_opt` produce mathematically equal
  closed forms (`log(exp(x) + exp(y))`, `-2*sin(x/2)**2`, `x**y - 1`)
  instead of the numpy/scipy function nodes.
- `combine_kronecker` leaves `KroneckerProduct(A, B) + KroneckerProduct(C,
  D)` unchanged. SymPy's `KroneckerProduct(A + C, B + D)` is not an
  identity. It also does not distribute non-integer powers.

What remains:

- `has_dups(Set())` is supported since c-doc (2026-10-23).
- `julia_code` and `octave_code` of a Function applied to a Matrix
  (`g(mat)`) are not supported: the core cannot build that expression.
- `MatrixKind`: `Integral` of a `MatrixSymbol` is supported since c-doc
  (2026-10-23).
- The native matrix-expression `str` differs from SymPy in four ways:
  - It orders sums differently: `-B*D**(-1)*C + A` where SymPy has
    `A - B*D**(-1)*C`.
  - It prints an `Inverse` of a product or power without parentheses
    (`C.T*C**(-1)`); `sstr` gets this right.
  - It prints rational `MatPow` exponents as `A**1/2`.
  - It prints a BlockMatrix as `BlockMatrix([...])` instead of
    `Matrix([...])`.

  These cause the remaining `schur` and `LUdecomposition` differences. The
  fix belongs in `to_str` of `crates/symox-core/src/matrixexpr.rs`.
- `galois_group` of degree 5 and 6 raises NotImplementedError. A port needs
  SymPy's galoisgroups, galois_resolvents and resolvent_lookup modules,
  plus the dense `dup_*` toolkits.
- `prime_decomp` and `primes_above` over
  `QQ.algebraic_field((T, theta))` with a symbolic `theta` work since
  c-doc (2026-10-23); see the c-doc section at the end.
- Unsupported `Poly` domains raise an error:
  - composite domains with inexact or non-generator coefficients (`'R[y]'`,
    `'Z[y]'` or `'Q(y)'` with `sqrt(3)`);
  - `CC`/`'C'`;
  - `RR_n` with n other than 53.

  The exception types match SymPy's, but the messages may be worded
  differently.
- Ring elements over `QQ<a>` print coefficients as `sqrt(2)`, where SymPy
  prints the raw `ANP(...)` repr.
- `Expr*elem` with an Expr that does not convert into the ring raises
  CoercionFailed, where SymPy raises TypeError.
- `ufuncify`: five doctest lines fail only in the reference environment.
  Symox matches the docstring there.

Tests: `tests/test_missing_c3_parity.py`,
`tests/test_missing_c3_matrices_parity.py`,
`tests/test_missing_c3_polys_parity.py`.

## Error classes and messages, category f (2026-10-23)

Invalid input now raises the reference exception class, and usually the
reference text. Covered: BlockMatrix (irregular blocks), ConditionSet
(non symbol-like dummy), banded (collision), Domain.exquo, Poly.exquo,
Poly.pexquo, exquo and pexquo ("... does not divide ..."),
all_roots/real_roots over EX or ZZ[pi] (NotImplementedError) and with a
non-Symbol generator (PolynomialError), ITE (`as_Boolean` TypeError),
ImmutableSparseMatrix assignment, Permutation call out of range
(IndexError), Permutation.resize, CantSympify and string parse failures
(SympifyError texts), padded_key, sift(binary=True), Ellipse.reflect,
diff, BaseDyadic.diff, integrate and Integral without variables, series
with a negative or infinite term count, linear_eq_to_matrix, linsolve and
linear_coeffs (nonlinear term and cross-term texts), solve_poly_system with
strict=True (UnsolvableFactorError), lpmin/lpmax (unbounded, empty and
inconsistent texts), testing.pytest raises/warns (Failed texts) and
aesara_function without aesara (ImportError). Regressions:
`tests/test_errtype_f_parity.py`.

Intentional divergences:

- `series(expr, x, x0, oo)` returns `expr` when it does not depend on `x`
  (the reference also does); a constant is never turned into an error.
- `solve_poly_system(..., strict=False)` still raises NotImplementedError
  when a factor is not solvable in radicals instead of returning a partial
  list.
- SympifyError also subclasses TypeError; the LP errors and CodeWrapError
  keep their extra ValueError/RuntimeError bases.
- `warns` re-emits warnings of other categories instead of raising them, and
  does not check the stack level.
- Messages naming a class or module show `symox.` paths
  (ImmutableSparseMatrix, warns_deprecated_sympy).
- autowrap and binary_function report the meson f2py command line.

Remaining gaps:

- Ge, Gt, Le, Lt `bool()`: the native text lacks the ": x < 1" suffix.
- ImmutableDenseMatrix assignment: native text "Cannot set values of an
  immutable matrix" (class matches).
- nsolve: the tolerance message lacks the mpmath residual numbers.
- The native ExactQuotientFailed is a DomainError/ValueError, not a
  BasePolynomialError; GF(p).exquo(1, 0) raises ZeroDivisionError instead of
  NotInvertible.
- `Derivative(x*y)` and `diff(Matrix([x*y]))` keep their own texts for a
  missing variable; `Integral(Eq(...))` names the equation, not its lhs.
- `ConditionSet((x, y), ...)` with a plain tuple dummy is rejected natively.

## Category-c batch 4: geometry, optics, vector, quantum, control

Reference algorithms ported (doctests now match SymPy 1.14): `Line` from an
equation (`linear_coeffs`), `RayTransferMatrix(A, B, C, D)`,
`CoordSys3D.orient_new`/`orient_new_quaternion` with a symbolic quaternion,
`Del.gradient` of a constant, `ImplicitRegion` regular/singular points,
multiplicity and rational parametrization (`vector/_implicitregion_ref.py`),
`Wavefunction` calls, `prob` and `norm`, `IntQubit(Qubit(...))`,
`integrate_result`/`represent` with symbolic delta roots, and the spin state
`rewrite` to `Jx`/`Jy`/`Jz` (coupled and uncoupled). `phase_margin` decides
`pm >= 180` from a 50-digit value when the native relational stays
symbolic, and raises TypeError when that value is ambiguous.

Deliberate choices:

- The `CoordSys3D` orthogonality check and the `QuaternionOrienter`
  unit-norm check were removed, as in the reference, which stores the
  matrix as given. The older tests that expected ValueError were updated.
- The integration term loop now calls the reference `deltaintegrate` on
  products the native rules decline, but only when every delta argument in
  `x` has slope +1 (or -1 with no derivative order) and at least one root
  is symbolic. Results that contain `zoo`, `nan` or a numeric
  `DiracDelta(c)` are refused. As a result, `x*DiracDelta(x - a)`
  integrates to `a*Heaviside(x - a)`. `DiracDelta(x)/x` stays unevaluated,
  where the reference returns `zoo*Heaviside(x)`. Numeric roots stay with
  the native preflighted rules.
- `ImplicitRegion.rational_parametrization` of a plane in three variables
  keeps the earlier exact linear solution, where the reference raises
  NotImplementedError.

Remaining gaps:

- Printing (printing owner): `pprint` and `pretty` of Python-level objects
  and containers fail. This covers MIMO systems, `TransferFunctionMatrix`
  (and `from_Matrix`), `jones_2_stokes` lists, the `Diagram`/`Diagram.hom`
  tuples and `DiagramGrid`. `DiagramGrid` also lays out as 3x1 instead of
  the reference 2x2, because `diagram_drawing` is not ported.
- Core/sets: `FiniteSet` of Python-level `Diagram`/`Object` instances
  (`Category`, `Diagram.subdiagram_from_objects`), `Lambda` with a `Matrix`
  body (`CoordinateSymbol`), `Integral` of an `NDimArray`/`Matrix`
  (`shape`, `ArrayKind`), and `Add` of Python `Vector` objects
  (`VectorKind`).
- `phase_margin`: the value is correct, but the native `arg` returns a
  different closed form.
- `couple`/`uncouple` with symbolic spins return the reference `Sum`
  since c-doc (2026-10-23); `rewrite` with symbolic spins still raises
  NotImplementedError.
- `simplify_gpgp` follows the reference since c-doc (2026-10-23), see the
  c-doc section at the end.
- Array expressions: the `sort_args_by_name`, `PermuteDims`,
  `convert_array_to_matrix` and `convert_indexed_to_array` doctests match
  since c-doc (2026-10-23); see the c-doc section at the end for the
  remaining canonicalization gaps.
- `ImplicitRegion.rational_parametrization` matches the reference since
  c-doc (2026-10-23), see the c-doc section at the end.
- Divergence: for three variables, `ImplicitRegion.regular_point` skips
  grid points where the gradient vanishes. The reference returns such a
  singular point, for example `(-10, -10, 20)` for `(x + y + z)**4`. This
  version raises NotImplementedError when every candidate is singular. In
  that case the reference fallback `list[self.singular_points()][0]`
  raises TypeError.

## Output forms of sympy.stats documentation examples (wave 4, stats-forms)

The audit listed 47 sympy.stats callables whose example output was equal to
SymPy's but written in another form. Re-running every example of their
docstrings (398 examples) against SymPy 1.14: 1 callable matched before
(`Normal`), 19 match now. The fixes follow the reference formulas exactly.
- CDF layout. A distribution can now carry the reference CDF expression,
  support Piecewise included (`_cdf_form`). It is ported for `Arcsin`,
  `Uniform`, `Kumaraswamy`, `UniformSum` (`(0, x < a), (F, x <= b), (1, True)`),
  `Dagum`, `Exponential`, `ChiSquared`, `Frechet`, `Gamma`/`Erlang`,
  `GammaInverse`, `Nakagami` (`(F, x >= lower)` or `x > 0`, then `(0, True)`),
  and the integrated form `Piecewise((F, x >= lower), (0, True))` of `Lomax`
  and `Benini`, which have no reference `_cdf`. The in-support `_cdf` is
  unchanged and still used for probabilities and sampling.
- `cdf(X, **options)` forwards the options. With options SymPy integrates the
  density, so `cdf(Gamma(k, theta), meijerg=True)` now gives
  `k*lowergamma(k, z/theta)/gamma(k + 1)` for `z >= 0`.
- `ExGaussian` CDF composes the two Gaussian CDFs in the scaled variable, as
  the reference does.
- Densities written as in the reference: `Gompertz`, `ShiftedGompertz` and
  `NormalGamma` keep separate `exp` factors; `Geometric` is `(1-p)**(k-1)*p`
  and `FlorySchulz` is `a**2*k*(1-a)**(k-1)` (no unsimplified `binomial`);
  `FDistribution`, `Benini` and `Nakagami` use the reference expressions;
  `Trapezoidal` uses half-open pieces `[a, b)`, `[b, c)`, `[c, d]`. Shapes
  where `b == a` or `c == d` still drop the empty piece.
- `GeneralizedMultivariateLogGamma(Omega)` sums over `Symbol('n',
  integer=True, negative=False)` with `gamma(n + 1)` in the denominator.
- `Wishart` and `MatrixGamma` densities fold the scalar into the explicit
  inverse scale matrix (`exp(Trace(-Sigma**-1*X/2))`, `exp(Trace(-Sigma**-1*X)/beta)`)
  with `Determinant(Sigma)**(-n/2)`. `MatrixNormal` drops explicit identity
  factors from the trace product.

Remaining form gaps, by root cause. All are numerically equal.
- Core printing (not changed here): multi-line `Matrix` inside `str`
  (`Wishart`, `MatrixGamma`, which are otherwise identical); pretty printing of
  `beta` as `B` (`FDistribution`), of `Sum` and `factorial` (`UniformSum` pdf).
- Core canonicalization (not changed here): `Rational**symbol` stays
  `(1/2)**n` instead of `2**(-n)` (`GeneralizedMultivariateLogGamma`,
  `Logarithmic` density); `beta(2, 2)` evaluates to `1/6` (`BetaBinomial`
  `density(X).dict`); `Trace` of the product keeps `(X - M)*(X - M).T`
  instead of `(X - M).T*(X - M)` (`MatrixNormal`, equal by cyclicity).
- Moments from closed forms plus `simplify`, where SymPy integrates and does
  not simplify: `ChiSquared` `moment(X, 3)` (factored vs expanded),
  `ExGaussian` `E`, `Rayleigh` and `Zeta` variance, `Nakagami` `E` and
  variance (`mu*gamma(mu)` vs `gamma(mu + 1)`), `H` (`entropy`, `log(2*pi)`
  vs split logs), `Logarithmic` `E`/variance, `GaussianInverse`
  `skewness(...).expand()` (symox leaves `exp*(cosh - sinh)` terms),
  `kurtosis`, `skewness` and `coskewness` (simplified quotients),
  `Probability.evaluate_integral`, `Multinomial` marginal (symox uses
  `p1 + p2 + p3 = 1`). Matching them needs SymPy's integration path and its
  absence of `simplify`, which would change many other results; left as gaps.
- Kept on purpose:
  - `LogCauchy` CDF keeps the `Piecewise((0, z < 0), ...)` guard. SymPy's
    `atan((log(z) - mu)/sigma)/pi + 1/2` is complex for `z < 0`.
  - `MatrixStudentT` stays exact; SymPy prints `pi**1.0` and float
    determinant powers.
  - `GaussianOrthogonalEnsemble` evaluates the normalising integral
    (`2*pi**(3/2)` for n = 2). SymPy leaves the `Integral` unevaluated.
  - `Circular*Ensemble` for n = 1 and `joint_eigen_distribution` give the
    evaluated density, where SymPy keeps an empty or held `Product`.
  - `MultivariateEwens` marginal substitutes `a1 = 2` inside the `Eq(a1, 2)`
    branch (`1/2` vs `1/factorial(a1)`).

Tests: tests/test_stats_forms_parity.py.

## Matrix-expression printing (c3 follow-up) on 2026-10-23

This closes the `str` gaps that the c3 section lists for matrix expressions.

Implemented:

- `str` and `repr` of a native matrix expression go through the Python
  `sstr`, as SymPy's `Printable.__str__` does. It orders MatAdd terms as
  SymPy does (`A - B*D**(-1)*C`). If the Python printer fails, the native
  string is used.
- The native `to_str` in `crates/symox-core/src/matrixexpr.rs` follows
  `StrPrinter`'s precedences:
  - `Inverse`, `Transpose` and `MatPow` put parentheses around a product,
    sum, power, Hadamard or Kronecker base: `(A*B)**(-1)`, `(A**(-1)).T`.
  - A `MatPow` exponent with precedence Pow or lower gets parentheses:
    `A**(1/2)`, `A**(2*x)`, `A**(x**2)`, `A**(-2)`. `A**sin(x)` and `A**x`
    get none.
  - A scalar is printed nested, so `A**0.5` drops the trailing zeros.
  - `_print_MatMul` sign rules: a negative real or negative imaginary
    coefficient becomes the leading sign (`-(2*I)*A`). A symbolic negative
    coefficient keeps its parentheses (`(-x)*A`).
  - `BlockMatrix` and `BlockDiagMatrix` print as SymPy's `Matrix([...])`
    table of blocks: right aligned columns, one row per line, zero blocks
    shown as `0`. An explicit matrix inside an expression uses the same
    table.
  - `MatrixSlice` drops a unit step, a zero start and a full stop
    (`X[:2, 1:]`), and puts parentheses around a compound parent.
  The native string still orders MatAdd terms by the internal order. It is
  used only as the fallback.
- Python `sstr`:
  - `sstr(Matrix)` no longer raises AttributeError (`_format_str`). It
    prints SymPy's right aligned table.
  - `sstr(BlockMatrix)` prints the block table instead of
    `BlockGrid(...)`.
  - `sstr` of a MatAdd works (no `is_Add` lookup on a native node).
  - `Adjoint`, `HadamardProduct`, `KroneckerProduct` and `HadamardPower`
    print under their SymPy class names, not `_MatrixExpr(...)`.
  - A slice stop at or past the dimension is omitted.
- The `schur`, `LUdecomposition`, `LDUdecomposition` and
  `UDLdecomposition` doctests of `BlockMatrix`/`BlockDiagMatrix` now match.
  On 82 probe expressions, `str`/`repr` now match SymPy on 80. Before this
  work, 17 differed.

What remains:

- The two remaining probe differences are about evaluation, not
  printing. SymPy's `*` and unary minus call `doit()`, so
  `Inverse(A*B)*Transpose(C + D)` becomes `B**(-1)*A**(-1)*(C.T + D.T)`
  and `-Inverse(A*B)` becomes `-B**(-1)*A**(-1)`. Symox keeps
  `(A*B)**(-1)*(C + D).T` and `-(A*B)**(-1)`.
- `latex` and `pretty` have no matrix-expression printers. They fall
  back to the `str` form (`A*B` instead of `A B`, `\left(A B\right)^{-1}`,
  or the 2-D pretty layout). A port needs SymPy's `LatexPrinter` and
  `PrettyPrinter` handlers for MatrixSymbol, MatMul, MatAdd, MatPow,
  Inverse, Transpose, Adjoint, BlockMatrix, Identity, ZeroMatrix,
  HadamardProduct and KroneckerProduct.
- `julia_code(g(mat))` and `octave_code(g(mat))` remain unsupported, so
  the `user_functions` docstring examples still fail. Building `g(mat)`
  raises TypeError ("cannot convert object of type Matrix to a symox
  expression"). SymPy holds an explicit matrix as a function argument. The
  native core cannot, and supporting it needs a matrix-valued argument
  node in the core.
- `conjugate(A*B)` evaluates to `Adjoint(A.T)*Adjoint(B.T)`, where SymPy
  keeps `Adjoint((A*B).T)`. This is evaluation, not printing.
- A symbolic slice stop, as in `A[0:1, 0:1]` with symbolic `n`, is stored
  as `Max(0, ceiling(Max(0, Min(1, n))))`, where SymPy stores `1`.

Tests: `tests/test_matexpr_print_parity.py`.

## Printing parity, audit batch print-a1

Implemented (ports of the reference code, with regressions in
`tests/test_print_a1_parity.py`):

- Printing helpers: `multiline_latex`, pretty-printed `print_mathml` XML,
  `rawlines`/`strlines`, `mathematica_code` of products with rational
  coefficients, the `Ordinal` exponent repr, and the codegen AST reprs
  (`Token`, `Type`, `intc`, `NoneToken`, `Variable`, `Comment`, `For`).
  `str(none)` is now `NoneToken()`, as in the reference, which replaces
  the earlier `None` divergence noted above. `strlines` and `rawlines` now
  return the reference text (`dedent(...)` blocks and raw literals), and
  `multiline_latex` accepts `terms_per_line=0` as the reference does.
- Computation fixes found through printing:
  - `SumApprox`, plus `reduce_order` and `reduce_order_meijer`, which now
    use the reference sort keys.
  - `hermite_prob` uses the reference `eval`: symbolic or negative orders
    stay unevaluated.
  - `Options` and `build_options`: `auto` is off when `domain` or `field`
    is given, and `extension` builds an algebraic field.
  - `decompose`, `best_origin`, `disambiguate`, `PropKB` clauses
    (`tell`/`retract` through CNF), `kernS`, `collect_abs`,
    `collect_sqrt` and `cse_release_variables`.
  - `function_range` returns `Interval(-oo, oo)` and not `Reals`.
  - `sympy_eqs_to_ring` returns the ring's domain, and `sring(...,
    field=True)` gives `ZZ(a)` rather than `QQ(a)`.
  - The public modules of `Tuple` and `AlgebraicNumber` match the
    reference.
- Core (Rust):
  - `conjugate` now passes through `beta`, `betainc`,
    `betainc_regularized`, `erfi`, `erf2`, `li` (unless the argument is
    extended-negative), the Fresnel and Airy functions,
    `lowergamma`/`uppergamma` (except at `z = 0` or `-oo`), `sech`,
    `csch`, and the orthogonal polynomials (the order is kept).
  - `ground_roots` keeps the `factor_list` order.
  - `Pow` of a `Pow` whose base is `1/q` reads it as `q**(-e)`, as the
    reference `as_base_exp` does. A symbolic exponent that is known to be
    an integer also combines: `((1/2)**s)**-1 = 2**s` and
    `(x**s)**k = x**(k*s)`.

Intentional divergences:

- `Subs(y, y, 0).has(y)`: the reference cached `has` result is wrong;
  symox returns the correct answer.
- Orders that depend on hashes (`args_cnc`/`atoms` sets, `itermonomials`,
  `denoms`, `decompose(separate=True)`) and the core canonical `Add`/`Mul`
  order (`Abs.aseries`, `gamma` series, `singularityintegrate`) are
  form-only.

Remaining unsupported:

- `UndefinedPredicate`: a native `Predicate` has different `.func` and
  `.name`.
- Order of models and mappings from `satisfiable` and `bool_map`; operator
  order in `devise_plan`.
- The type name of the native `Domain`.
- Alias printing for `AlgebraicField` (`alg_field_from_poly`,
  `cyclotomic_field`, `integral_basis`).
- `DDM`, `SDM` and `DomainMatrix` elements print as `MPQ`, not the
  flint-style mix.
- `Add`/`Mul` of a `MatrixSymbol` do not become `MatAdd`/`MatMul`.
- `ImageSet` intersected with a `FiniteSet`, and `solveset` with an
  imaginary symbol.
- Pretty printing of sets, `Order`, `hyper`, `meijerg` and `Integral`
  (covers `solveset`, `stationary_points`, `Integers`, `Naturals`, the
  ODE power-series records and `odesimp`).
- Last-digit `evalf` precision for `catalan`, `jn` and `yn`.
- `sring` with a float coefficient and a symbolic one panics. This
  predates this batch.

## Printing parity, batch a2 (physics, diffgeom, combinatorics, tensor, categories, geometry)

Implemented (tests in `tests/test_print_a2_parity.py`):
- Quantum: `Ket`/`Bra`/`TimeDepKet` labels (`|n,m>` only for spin states),
  state classes in `symox.physics.quantum.state`, gate names (`H(0)`,
  `C((0),Z(1))`), `CPHASE` built as `CGateS` as in the reference, operator
  powers (`A**(-1)`, `sqrt(A)`), and adjacent powers merge (`A*A.inv() == 1`).
- Diffgeom: base scalar fields print by coordinate name (`x**2`, LaTeX
  `\mathbf{x}`), and `CoordSystem.transformation` prints as a `Lambda`.
- Tensor: `TensMul`/`TensAdd` signs and coefficients, the `WildTensor*`
  module path, and `Indexed`/`IndexedBase` repr and LaTeX.
- Values corrected so that they print like the reference: optics
  `phase_retarder`, `mueller_matrix` and `polarizing_beam_splitter` (the
  reference formulas); `Cable.loads_position` lists; `transpose` of an
  `Array` stays an `Array`; `Polyhedron` generators and `copy`;
  `Segment.perpendicular_bisector`; `Polygon.intersection` ordering, using
  a geometry `sort_key` that mirrors `Basic.sort_key`; the keys of
  `DimensionSystem.get_dimensional_dependencies` (symbols only for base
  dimensions, not for tabled derived ones); and the categories,
  `Partition` and `Prufer` printers.
- `repr(Matrix)` strips trailing Float zeros, as the reference `sstr`
  does. `str(Matrix)` keeps full precision.

Form-only (core canonical `Add`/`Mul` order, left as is):
- `metric_to_Ricci_components`, `metric_to_Riemann_components`,
  `vectors_in_basis` and `get_contraction_structure`;
- `FiberForceLengthPassiveInverseDeGroote2016` (`log(x + 1)`), `psi_n` and
  `PIDController.to_expr`, whose reference extra parentheses come from an
  unevaluated `Pow`;
- the term order of `Commutator.expand` and `Operator.expand`, and the order
  of `TensAdd` terms.

Ordering only:
- `PIDController.poles` (order of `roots`);
- `find_dynamicsymbols`, `PermutationGroup.conjugacy_classes` and
  `PermutationGroup.generate` (set/generation order).

`Bra.subs` combines `<mm| - I*<mm|` to `(1 - I)*<mm|`. The value is equal.

Reference-wrong: `wigner_9j` at high precision. This version keeps the
correctly rounded last digit.

Remaining gaps:
- Xypic layout (`XypicDiagramDrawer`, `xypic_draw_diagram`): `DiagramGrid`
  is a simplified port.
- `MatplotlibBackend.append`/`extend`: `plot` returns a `Plot`, not a
  backend subclass.
- `Reshape.as_explicit` prints `ArrayElement(A, 0)` instead of `A[0]`.
- LaTeX of quantum states inside `FiniteSet` (sets printer).
- Diffgeom LaTeX of `TensorProduct`/`WedgeProduct` (`\otimes`, `\wedge`),
  and the unevaluated `Differential` of a scalar expression.

## Trigonometric and nested exponential solve, remaining items (2026-10-23)

Implemented:

- `solve(sec(x) - r)` and `solve(csc(x) - r)`, when the reference `acos`
  or `asin` of `1/r` does not evaluate. The equation is rewritten to a
  quadratic in `t = exp(I*x)`, as `_tsolve` does through
  `rewrite(exp)`, and each root gives `-I*log(t)`. At the top level the
  reference simplify splits a rational content `1/c` out as
  `I*(log(c) - log(p))` (for example `I*(log(3) - log(1 - 2*sqrt(2)*I))`).
  Under a generator such as `sec(x)**2 - 9` the plain form is kept, as in
  the reference.
- `solveset` of nested exponentials over the complexes (`exp(exp(x)) - 1`,
  `exp(exp(exp(x))) - 1`, `exp(2*exp(x)) - 3`, `exp(exp(x) + 1) - 1`,
  `exp(exp(2*x)) - 1`). Each level of `_invert_complex` for `exp` adds a
  `k` dummy and gives `ImageSet(Lambda(((_k1, _n),), I*(2*_k1*pi +
  arg(body)) + log(Abs(body))), ProductSet(Integers, Integers))`. An outer
  affine map stays held as `ImageSet(Lambda(x, ...), inner)`. Over the
  reals the empty set is kept.
- `solveset` with a symbolic coefficient (`_solve_trig1` with a symbolic
  `mu`, and `_invert_trig_hyp_complex` with a symbolic right-hand side):
  `sin(a*x) - b`, `cos(a*x) - b`, `tan(a*x) - b` and `sin(a*x + 1) - b`
  give `ConditionSet(x, Ne(a, 0), Union(ImageSet(Lambda(_n, (2*_n*pi +
  asin(b))/a), Integers), ...))`. Over the reals, `sin` and `cos` add the
  condition `(-1 <= b) & (b <= 1)`. When the trig function occurs more
  than once (`sin(a*x) - cos(a*x)`), the exponential polynomial path keeps
  the image sets unexpanded, intersects them with the domain when it is
  real, and guards them with `Ne(a, 0)`.
- `complete_roots` in `poly/systems.rs` now tries the radical certificate
  before the simplify-based residual check. The two checks are still both
  required to fail before a root is rejected.

Probe comparison (40 probes, reference SymPy 1.14): 8 matched before and
27 after. The tests are in `tests/test_solve_rest2_parity.py`.
`tests/test_solveset_complex_wave4_parity.py` no longer lists
`exp(exp(x)) - 1` as an undecided equation.

Intentional divergence:

- `solve(sin(3*x) - sin(x))` also returns `pi`, which is a genuine root
  that the reference loses.

Remaining gaps:

- `sec(2*x) - 3`, `csc(3*x) + 4` and other non-unit arguments: the
  reference solves a polynomial of degree 4 or 12 in `exp(I*x)`, with
  forms such as `-I*log(1/3 - 2*sqrt(2)*I/3)/2` and cube roots. Symox
  raises NotImplementedError. `sec(x/2) - 3` raises in both, but the
  message text differs.
- `solve(sec(x) - a, x)` prints one root as `-I*log((1 - sqrt(1 -
  a**2))/a)`. The reference prints the equal `-I*log(-(sqrt(1 - a**2) -
  1)/a)`.
- `solveset(exp(exp(x)) - y, x)` with a symbolic `y` stays a ConditionSet,
  because the nested inversion requires a value proved nonzero. The
  reference returns the `ProductSet(Integers, Integers)` image.
- The Union order of the symbolic `_solve_trig1` results (`sin(a*x) +
  cos(a*x)`, `sin(2*a*x) + sin(a*x)`) differs. The reference order follows
  the creation order of the `_n` dummies (`default_sort_key` compares their
  indices). The sets are equal.
- `solveset(cos(ap*x) - b, x, Reals)` with a positive `ap` goes through
  the numeric-coefficient path. It prints `2*_n*pi/ap + acos(b)/ap`, where
  the reference prints `(2*_n*pi + acos(b))/ap`. The sets are equal.
- `nonlinsolve([x**3 + 1/y - 2, x - y], [x, y])` still does not finish.
  The time goes into `normalize` (`simplify(expand(e))`) in
  `poly/systems.rs` `descend`/`complete_roots`. That call reaches
  `cancel`, and from there the algebraic-field multivariate gcd
  (`algebraic_multi::gcd_from`), on the substituted cubic radicals. The
  reference does not simplify there. It keeps the raw `roots` output for
  `y` and an EX-domain form for `x`. `[x**3 + 1/y - 2, x + y]` raises
  NotImplementedError. `[x**3 + 1/y, x - y]` prints `y` in the factored
  `sqrt(2)*(1 + I)/2` form, where the reference prints the expanded
  `sqrt(2)/2 + sqrt(2)*I/2`.

## Matrix-expression LaTeX and pretty printing on 2026-10-23

The native `latex` and `pretty` printers now port the SymPy 1.14 handlers
for matrix expressions. Before this work, symox printed these expressions
through their `str` form. The handlers cover:

- MatrixSymbol, with `mat_symbol_style="bold"` in `latex`.
- MatMul, MatAdd and MatPow.
- Inverse, Transpose, Adjoint (with the `adjoint_style` setting) and
  conjugate.
- HadamardProduct, HadamardPower and KroneckerProduct.
- BlockMatrix and BlockDiagMatrix.
- MatrixSlice and MatrixElement.
- Trace and Determinant.
- Identity, ZeroMatrix and OneMatrix.
- Explicit matrices, with the `mat_delim` and `mat_str` settings, the
  inline `smallmatrix` form and the `array` form for more than 10
  columns.
- ElementwiseApplyFunction, DotProduct and Lambda.

The parenthesization follows SymPy's `precedence_traditional`, and MatAdd
signs follow SymPy's rules. In pretty output a negative term is joined
with a space. In LaTeX, " - " is used only when the coefficient is a
negative Number, so `A + - i B` prints as in SymPy. The MatAdd arguments
are now sorted with `default_sort_key`, as SymPy's canonicalize does. The
native `str` order of MatAdd terms is therefore SymPy's `args` order. The
Python `str` printer still orders the terms as SymPy's `str` does.
When the expression tree of the terms exceeds 4096 nodes (for example
the shared DAG of `block_collapse(X**1000000000)`), the sort falls back to
the structural order. Computing `default_sort_key` on that tree is
exponential.

There are 126 main probes and 70 extra probes (settings, nesting, edge
cases). On the main probes, the outputs that differ from SymPy went from
121 to 2 for latex, from 94 to 2 for pretty ascii and from 117 to 2 for
pretty unicode. On the extra probes, they went from 67/45/56 to 2/2/2.
The `str` time of `A*B + C` is unchanged: 131 us, against 509 us for
SymPy.

The remaining differences come from evaluation or representation, not
from the printers:

- SymPy evaluates `(A*B)[0, 1]` to a Sum. Symox keeps the element.
- `(A + B)[:2, 1:]` with symbolic `n` keeps the normalized
  Max/Min/ceiling stop.
- `conjugate(A*B)` is distributed, as noted in the previous section.
- `X[:2, :2].T` is rewritten as a slice of the transpose.

Still unsupported:

- PermutationMatrix in cyclic notation (the `perm_cyclic` setting).
- FunctionMatrix, DiagMatrix, DiagonalOf and Permute. In `latex` and
  `pretty` they fall back to the `str` form.
- `mat_symbol_style` for `pretty`.
- `decimal_separator` inside matrix expressions.

Tests: `tests/test_matexpr_pretty_parity.py` (177 expressions in the 3
printers, plus 13 settings cases, 4 of them skipped as known divergences).

## Printing parity, audit batch print-a3

This batch follows up on the gaps of print-a1. On the same 90 records, the
probe went from 50 differing records to 42. Regression tests are in
`tests/test_print_a3_parity.py`.

Implemented:

- Pretty printing (ASCII and unicode) and latex of sets, ported from the
  reference `_print_FiniteSet`, `_print_Interval`, `_print_Union`,
  `_print_Intersection`, `_print_Complement`, `_print_ProductSet`,
  `_print_ImageSet`, `_print_ConditionSet` and `_print_Range`, and of the
  named sets (`Integers`, `Naturals`, `Naturals0`, `Reals`, `Complexes`,
  `Rationals`, `EmptySet`, `UniversalSet`). This covers the results of
  `solveset` and `stationary_points`. Several base sets of an `ImageSet`,
  which symox stores as one `ProductSet`, print as separate sets.
- Pretty printing of `Order`, `hyper`, `meijerg` and `Integral`, and
  latex of `Order`, `hyper`, `meijerg` and `Integral` (`\iint`-style
  signs, `\limits`, which is left out inline and with `itex`). This
  covers the ODE power-series records and `odesimp`.
- Elements of an algebraic field print as the reference `ANP(rep, mod, QQ)`.
  `EX` elements print as `EX(...)` and `EXRAW` elements as the raw
  expression, in both `str` and `repr`.
- `construct_domain(..., extension=True)` takes as generators the leaves of
  the `Add`/`Mul` trees of the coefficients, as the reference
  `_construct_algebraic` does: `[1 + sqrt(2), 1]` gives `QQ<sqrt(2)>`.

Intentional divergences:

- `DDM`, `SDM` and `DomainMatrix` over `QQ` (16 records) hold `MPQ`
  elements. The reference keeps the raw Python `int` given for `QQ`. Its
  results then mix types (`rref` returns floats, `inv` raises
  `TypeError`). The reference is wrong here, and the exact conversion
  stays, as in WP-C5.

Remaining unsupported:

- Alias printing for `AlgebraicField` (`alg_field_from_poly`,
  `cyclotomic_field`, `integral_basis`). It needs an aliased
  `AlgebraicNumber` atom in the core expression tree.
- The type name of the native `Domain` and `DomainElement` classes. They
  are single Rust classes and cannot be subclassed per domain.
- A `PolyElement` over an algebraic field prints its coefficients as
  expressions, not as `ANP`.
- `Add`/`Mul` of a `MatrixSymbol` do not become `MatAdd`/`MatMul`. Not
  attempted: it changes core construction, which the matrix and core
  caller tests depend on.
- Last-digit `evalf` precision for `catalan`, `jn` and `yn`. Not attempted
  in this batch.
- `FiniteSet(...).intersect(ImageSet(...))` stays unevaluated. The
  reference tests each element for membership.
- `solveset` over the reals with an imaginary symbol keeps the symbol in
  its `ConditionSet`. The reference uses a real dummy `_R`.
- Pretty printing of `Sum` and `Limit` still uses the `str` form.

## Cubic and quartic roots in rational `nonlinsolve` systems on 2026-10-23

`nonlinsolve([x**3 + 1/y - 2, x - y], [x, y])` did not finish. SymPy
clears the denominators, takes the grevlex Groebner basis of the numerators,
converts it to lex with FGLM and calls `solve_poly_system`. Its
`solve_generic` keeps the raw `roots(f.ltrim(gen))` output of the
univariate basis element over QQ. It neither simplifies nor re-verifies
these roots. Symox checked each root by substituting it into the
polynomial. When the radical certificate failed, it fell back to
`simplify(expand(residual))`, and for the complex cubic and quartic roots
that call does not end (it goes through `cancel` into the
algebraic-field multivariate gcd).

The native system solver (`poly/systems.rs`, `complete_roots`) now follows
the reference:

- Roots of a univariate polynomial with rational coefficients are the
  exact `roots` formulas. They are kept unsimplified and are not
  re-verified, as in SymPy. The value of the last lex variable now prints
  like the reference: `y = sqrt(2)/2 + sqrt(2)*I/2` in `[x**3 + 1/y, x - y]`,
  and the cubic-formula forms such as `-1/3 - 2/(9*(17/27 +
  sqrt(33)/9)**(1/3)) + (17/27 + sqrt(33)/9)**(1/3)`.
- The root `-c0/c1` of a linear polynomial whose leading coefficient is
  certified nonzero is exact by construction (SymPy's `roots_linear`), so it
  no longer needs a residual zero certificate. Other degrees with
  non-rational coefficients keep the certificate.
- A value already assigned is re-normalized only if it contains the variable
  just solved.

Denominator filtering is unchanged: roots where a denominator of the
original system vanishes are still dropped (intentional divergence, see the
rational nonlinsolve section above).

A probe of 46 rational, polynomial and parametric systems (2 and 3
unknowns, degrees 1 to 5, denominators, parameters) was run with a 60 s
timeout per call. It went from 22 outputs identical to SymPy 1.14, 1 hang
and 5 errors to 26 identical, 0 hangs and 3 errors. No identical output
became different. Newly solved: `[x**3 + 1/y - 2, x - y]` (0.2 s),
`[x**3 + 1/y - 2, x + y]` (it used to raise NotImplementedError, SymPy
returns the four quartic-formula solutions; symox takes 5 s) and
`[x**3 + 1/y - 3, x - y]`. The other differences are equal sets in another
printed form, or intentional.

Remaining gaps:

- The components solved after the last variable (the `x` of the examples
  above) print in another form. SymPy computes them with `roots` over the
  EX domain. There, `preprocess_roots` clears the denominators before
  `roots_linear` simplifies `-c0/c1`, which gives forms such as
  `(-(17 + 3*sqrt(33))**(1/3) + 4 + ...)/(3*(1 - sqrt(3)*I)*(17 +
  3*sqrt(33))**(1/3))`. Symox simplifies the substituted equation. The
  values are equal.
- With radicals of the quartic formula, several radicands are exactly
  negative reals. At 15 or 30 digits, `evalf` (in SymPy too) can land on
  the wrong side of the branch cut. A 30-digit substitution check then
  reports non-solutions that are correct: at 80 digits, through mpmath,
  they are exact. SymPy's own output for `[x**3 + 1/y - 1, x - y]` shows
  the same effect.
- Still unsupported (NotImplementedError, SymPy solves them):
  `[a*x + y - 1, x - y]` and `[x**2 + a/y, x - y]` with a free parameter
  `a`. Also `[x**2 + 1/y - 1, x - y**2]`, whose quintic factor SymPy
  returns as `CRootOf` through the substitution path.

## Equivalent-form parity, audit batch form-b2 (2026-10-23)

These are records where the value was already equal to SymPy's but the form
differed. Each fix ports the reference step that builds the result.

Implemented:

- Geometry:
  - `Ellipse.tangent_lines` uses SymPy's rise/run formula on the curve and
    its single-solution branch for outside points.
  - `perpendicular_line` follows SymPy's 2D and 3D constructions.
  - `Line3D.distance` for skew lines goes through the plane with the
    cross-product normal.
  - `Plane.distance` is the projection on the unit normal.
  - `direction_ratio` is componentwise in 3D. This also fixes
    `deviation` and `refraction_angle`.
- `Matrix.eigenvects(simplify=...)` and `Matrix.diagonalize` are SymPy
  ports: primitive eigenvectors, the same errors, and an immutable result
  for immutable input.
- DiracDelta expansion uses SymPy's `roots` algorithm. Two guards keep
  exactness:
  - When `roots` misses a root (a quintic), the exact CRootOf fallback
    runs. SymPy returns 0 there.
  - When the scale `|g'(r)|` may vanish (`DiracDelta((x**2 - 1)*y)`),
    the delta is kept whole. SymPy divides by `Abs(y)`, which is
    undefined at y = 0.
- `erf2inv` special values follow SymPy's `eval`.
- `factor(deep=True)` runs SymPy's bottom_up `_try_factor` followed by the
  partials pass.
- `together(deep=True)` also combines function arguments and power
  exponents.
- `unrad` canonicalizes the result as SymPy's `_canonical` does: it drops
  numeric factors and factor powers, then applies `could_extract_minus_sign`
  per factor, instead of making the leading coefficient positive.
- `poly_from_expr` passes `wrt` and `sort` through, so
  `horner(f, wrt=y)` now matches SymPy.
- `split_surds`, `_split_gcd` and `rad_rationalize` are direct ports.
- Vector: `curl`, `divergence` (`doit=False` holds the Derivatives),
  `gradient`, `laplacian` and `Vector.integrate` are ported.
- Physics:
  - Ported from SymPy: musculotendon, activation and the curve mixin;
    `E_nl_dirac`, `R_nl` and `wigner_d_small` (which returns an
    ImmutableMatrix); `Medium`; `TWave` (repr, superposition, speed);
    `conjugate_gauss_beams`; `geometric_conj_ab`.
  - Kept beyond SymPy 1.14: musculotendon formulations 3 and 4,
    `conjugate_gauss_beams(s_in=...)`, and the corrected Rayleigh range and
    magnification of `gaussian_conj`.
- Other modules:
  - Ported: `factorint(limit=...)` (trial division, Fermat, pm1 and rho
    rounds), `find_linear_recurrence`, `from_hyper` initial conditions,
    hyperexpand `_mod1`, and `lambdastr`.
  - `fnodes` `FunctionCall` and `sum_` are noncommutative.
  - `lambdify` with Python ints, an implemented function or a
    DeferredVector uses exact integer arithmetic.

Intentional divergences:

- The general case of `erf2inv(x, y)` evaluates to `erfinv(y + erf(x))`,
  where SymPy holds `erf2inv(x, y)`.
- The two DiracDelta guards above.
- diffgeom objects are evaluated eagerly, where SymPy keeps lazy objects.
  The values are equal.

Remaining form gaps (values are equal):

- Rust core:
  - `Matrix.exp` uses skew_exp or analytic_func, not P*exp(J)*P^-1.
  - `expand_power_base` splits only when every factor is positive. SymPy
    also pulls out the nonnegative factors. This also affects
    `polarify(...).expand()`.
  - `(-2)**i` with an odd `i` keeps the base.
  - `2*2**(-m)` merges to `2**(1 - m)`. This affects `hankel_transform`
    and `inverse_hankel_transform`.
  - `diff(..., n)` distributes the numeric coefficient.
  - `expand_func(gamma(p/q))` does not shift by the integer part, so
    `expand_func(integrate(fresnel integrand))` keeps a gamma ratio.
- Simplification and transforms:
  - `simplify` picks `-(x**2/y**2 + 1)/y` over `(-x**2 - y**2)/y**3`
    (`idiff(...).simplify()`, and a custom `measure`).
  - `trigsimp(method='old')` simplifies more than SymPy's `trigsimpold`.
  - `trigsimp_groebner` ignores `hints`.
  - `inverse_mellin_transform` and `laplace_transform` build their
    results differently.
- Solvers:
  - `solve` on the `solve_undetermined_coeffs` system returns
    `c - b**2/(4*a)`, not `(4*a*c - b**2)/(4*a)`.
  - `reduce_inequalities` does not add the `x < oo` bound for symbolic
    coefficients. It also raises NotImplementedError for `x*y > 1` with
    respect to `x`, where SymPy returns the relational unchanged.
  - The `1st_linear_constant_coeff` PDE: the exp product and Integral
    hint forms differ.
  - `solve_decomposition` ImageSet output differs in pretty printing and
    in the order of the union terms.
- Polys: `factor_list`, `sqf_list` and `div` raise NotImplementedError for
  `wrt`.
- Vector: top-level `integrate(vector, x)` raises TypeError.
- Physics:
  - Bra*Ket is evaluated at multiplication time. SymPy keeps an
    unevaluated InnerProduct (`OrthogonalKet`, `Qubit`).
  - The core multiplies a Float into a sum that contains an
    UnevaluatedExpr.
  - `-(-10 + c)**3` is printed where SymPy has `(10 - c)**3`.
  - `limit` and `expand` forms differ in the PIDController records.
  - `simplify` gives a different form for the `ObstacleSetPathway` result.
- Other modules:
  - The `Product.change_index` and `reverse_order` forms differ (Rust
    summation term order and simplify).
  - The core distributes the numeric content into the
    `find_linear_recurrence` generating function.
  - The `round_two` and `integral_basis` doc forms need
    `QQ.algebraic_field((T, theta))` and `alias` support in the native
    algebraic domain.
  - `lookup_origin` uses native closed forms, not SymPy's `HyperRep_*`
    table.
  - The holonomic `_normal` sign and scale differ.
  - `lambdify` on numpy int arrays returns floats.

Tests: `tests/test_forms_b2_parity.py` and `tests/test_forms_b2_parity_*.py`.

## Printing order of Add terms and Mul factors (print-order, 2026-10-23)

`str`, `latex` and `pretty` order the terms of a sum and the factors of a
product as SymPy 1.14 `Add.as_ordered_terms` and `Mul.as_ordered_factors`
do. The stored canonical order of `Add`/`Mul` arguments is unchanged.
`printing::order` holds one implementation of `as_terms` and of
`default_sort_key`, and the keys of one print are computed once and shared:
- `as_terms`: every numeric factor multiplies into the complex coefficient,
  where `complex(zoo)` is `nan+nanj`. The generators are sorted by
  `default_sort_key`, and the sort is stable. With an `Order` term, both
  groups are sorted in reverse and the `Order` terms follow the regular
  terms.
- The special case of SymPy is kept: a positive Number or NumberSymbol stays
  before a negated single factor (`2 - 3*x`, `pi - x`). `oo - x` prints as
  `-x + oo`.
- `Function.class_key` uses the `funcs` table. A held undefined head whose
  name is a SymPy class with a finite `nargs` orders with 10000, as SymPy
  does, which puts `f(1)` before `SingularityFunction(x, 1, 0)`. `Max` and
  `Min` use `Basic.class_key`, and `Heaviside(x)` sorts with the arguments
  `(x, 1/2)`.
- Held sums and products keep the `make_args` order of their held
  arguments, as before.
- Series coefficients of `gamma` and summation results hold the `zeta`
  function, not an undefined `zeta`. `qho_1d` uses the quantum `hbar`
  number symbol, which prints as `\hbar` in LaTeX and as the Planck glyph in
  pretty. A host atom's one-line `_pretty` form now reaches the native
  pretty printer.
- `OperatorExpr` (`Commutator.expand`, `Operator` powers) prints through the
  same `Add`/`Mul` order.

Probe: 1973 random expressions with seed fixed, against SymPy 1.14.

| | str | as_ordered_terms | as_ordered_factors | latex | str order-only |
|---|---|---|---|---|---|
| before | 182 | 175 | 183 | 1182 | 57 |
| after | 138 | 138 | 139 | 1176 | 17 |

The remaining `str` mismatches are not caused by ordering. They are
evaluation differences: `(1/2)**x` versus `2**(-x)`, the distribution of
`re`/`im`, and powers of `I`.

Remaining gaps:
- `re`, `im`, `Heaviside` and `DiracDelta` of noncommutative arguments are
  commutative in SymPy and noncommutative here. This explains the 17
  remaining order-only mismatches.
- `hbar.is_number` is False here and True in SymPy, which leaves
  `psi_n(n >= 1)` with a different factor order.
- One probe expression segfaults symox, where SymPy raises `ValueError`.
  It was skipped.
- Printing gaps unrelated to order: LaTeX `\zeta{\left(3 \right)}` versus
  `\zeta\left(3\right)`; `Order` in LaTeX and pretty, printed as
  `Order(x**3, x, 0)` instead of `O(x**3)`; LaTeX and pretty of
  `SingularityFunction`; LaTeX and pretty of `Commutator` expansions, which
  print the `OrderedProduct` structure.

## Special-function evalf digits

A probe (seed 20261023) compares `str(N(e, n))` with SymPy 1.14 for about
60 special functions, about 10 arguments each, at 15, 30 and 50 digits
(2796 cases). Before this work: 2641 exact, 5 last-digit, 97 larger
differences and 53 errors. After: 2680 exact, 5 last-digit, 85 larger and
26 errors. Every remaining difference was checked against mpmath at 3n+10
digits, or is listed below.

What was implemented:

- The gamma identities used for `beta`, `rf`, `ff`, `binomial` and
  `catalan` take Float arguments as their exact binary value. Before, the
  inner gamma of a Float was computed at 53 bits, so `beta(5/2, 0.75)`,
  `rf(5/2, 0.75)` and `ff(1, 2.5)` were wrong from about the 17th digit at
  30 and 50 digits.
- `beta(x, y)` with numerical arguments follows SymPy's `beta.doit`. A
  zero argument gives `zoo` (including `beta(0, 0)`, which gave `nan`).
  `y == 1` or `x == 1` gives `1/x` or `1/y`. `y == x + 1` gives
  `1/(x*y*catalan(x))`. A negative integer sum of two non-integers gives 0.
  Any other case stays unevaluated: `beta(2, 4)` and `beta(-1, 2)` are no
  longer evaluated to `1/20` and `zoo`.
- `N(lowergamma(a, z))` is `oo` at a pole of mpmath's `gammainc(a, 0, z)`:
  `a` a nonpositive integer, `z != 0` and `re(z) >= 0`. When `re(z) < 0`,
  mpmath's ValueError "gamma function pole" is raised, as in SymPy.
- SymPy's ValueErrors when the expression is built: `laguerre` with a
  non-integer degree, `hermite` with a negative degree and
  `assoc_legendre` with a negative degree or `|m| > n`. The messages are
  SymPy's.
- `assoc_laguerre` with a non-integer degree evaluates through mpmath's
  `laguerre`, as `jacobi` and `gegenbauer` already did.

Kept on purpose (symox is correctly rounded and SymPy is not; the tests
check these against mpmath):

- `catalan(I)` and `catalan(1/2 + I)`. SymPy evaluates the gamma factors
  separately.
- `jn`/`yn` with a complex Float argument, and `jn(1, 2.5)`/`yn(1, 2.5)`.
  SymPy's rewrite computes `1/z` at 53 bits.
- `Li(2.3)` and `catalan(2.3)`. They are evaluated when built, at 53 bits.
  SymPy loses about 9-10 ulp there to cancellation.
- `zeta(107/3, 4)` at 50 digits, and the last digit of `airyaiprime`,
  `airybiprime`, `bernoulli` and `euler` at large arguments.
- `jacobi`, `gegenbauer` and `assoc_laguerre` with a non-integer degree.
  SymPy truncates the degree when it builds the expression: for example
  `jacobi(1/3, a, b, x)` is 1 and `gegenbauer(-1/2, a, x)` is 0. Symox
  evaluates the analytic continuation, as mpmath does.
- `chebyshevu(-1/2, 3)` evaluates. SymPy raises a TypeError (a bug in its
  negative-degree branch).
- `lucas` with a non-integer argument uses the Binet form (see the `lucas`
  section above).

Remaining gaps:

- `elliptic_pi`. mpmath's `ellippi` time is erratic in the precision, and
  symox's precision refinement re-evaluates it, so 7 probes time out
  (10 s) where SymPy takes 1 to 15 s (SymPy times out on 1 other probe).
  The 7 differences at 30 and 50 digits are closer to mpmath than SymPy's,
  but `elliptic_pi(2, 7/2)` at 50 digits is still wrong from about the
  45th digit.
- `jacobi(n, ...)` and `assoc_laguerre(n, ...)` with a negative
  non-integer `n` raise a ValueError in SymPy, through order-dependent
  transformations. Symox evaluates them (6 + 6 probes).

Tests: `tests/test_evalf_digits_parity.py`. It has 37 SymPy digit
strings, 15 construction cases, 14 domain errors, and mpmath checks for
the cases kept on purpose.

## Small gaps: inequalities, wrt, vector integrate, lambdify, matrices, Sum/Limit printing, ANP

Implemented (SymPy 1.14 behaviour):

- `reduce_inequalities` with a relational in several symbols and one symbol
  of interest goes through a port of `_solve_inequality`: a safe linear
  inversion that divides only by a factor of known sign. Symbols of
  unknown realness are recast to real dummies as in the reference.
  Otherwise the input stays unchanged (`x*y > 1`), or NotImplementedError
  is raised when no solver applies.
- The polys `wrt` and `sort` options order the automatically detected
  generators, as `_sort_gens` does. This covers `factor_list`, `sqf_list`
  (through a port of `_generic_factor_list`), `div`, `rem`, `quo`, `pdiv`,
  `prem`, `pquo`, `exquo`, `pexquo`, `resultant`, `subresultants`,
  `compose`, `reduced`, `groebner`, `sqf_part`, `discriminant`,
  `decompose`, `sturm`, `gff_list` and `ground_roots`. Explicit
  generators ignore the option, as in the reference.
- Multivariate `decompose`, `gff_list` and `ground_roots` raise the
  reference errors.
- `integrate(vector, ...)` and `Integral(vector, ...)` work componentwise
  through `_eval_Integral`. A single term is evaluated; a sum keeps held
  Integrals, as in SymPy.
- `lambdify` keeps the integer dtypes of numpy arrays and scalars, and is
  exact on Python ints such as `2**60`.
- Pretty printing of `Sum` and `Limit` is ported from `_print_Sum` (with
  `asum`) and `_print_Limit`, in both ASCII and unicode. The `Limit`
  LaTeX is `\lim_{x \to 0^+}\left(...\right)` for an Add or Mul, and
  `\lim_{x \to 0^-} \frac{1}{x}` otherwise.
- A PolyElement over an algebraic field prints its coefficients in the
  reference `ANP(rep, mod, QQ)` form, for both str and repr. ANP has Atom
  precedence and the sign comes from the leading coefficient.
  `ANP * PolyElement` now defers to the ring element.
- `Add` and `Mul` of MatrixSymbols already build MatAdd and MatMul (same
  `func`, isinstance, str and shape as SymPy).

This closes the `wrt`, vector `integrate` and numpy-int `lambdify` items
in the forms B2 gap list above.

Intentional divergences:

- `Add(A, 1)` with a MatrixSymbol `A` raises TypeError. SymPy builds a
  meaningless scalar+matrix `Add` (`1 + A`).
- `type(e).__name__` of a matrix expression is the native class name.
  `.func` and isinstance follow SymPy.
- Error messages differ in wording but have the same exception type:
  - the ShapeError of `Add(A, C)` and `Mul(C, A)`;
  - the TypeError of `A + 1`;
  - the NotImplementedError of `reduce_inequalities(x**2*y > 1, [x])`.
- The ANP repr assumes the reference's python ground types (`MPQ(p,q)`),
  i.e. a SymPy without gmpy2 or flint.

Still unsupported or different:

- The pretty printing of `Integral` still falls back to str.
- In an Add that holds both a `Sum` and a `Limit`, the term order differs
  (symox prints the Limit first). This is a general printer ordering
  gap, not a Sum/Limit printer gap.
- The global `div`, `rem`, `quo` and `exquo` still raise
  NotImplementedError for `greedy=True`, which the reference accepts, and
  for `symmetric=True`, where the reference raises OptionError (this
  predates this change).

Tests: `tests/test_gaps_misc_parity.py`.

## Batches c1 and c2, follow-up (c-rest-a, 2026-10-23)

Ported from SymPy 1.14, regressions in `tests/test_c_rest_a_parity.py`:

- `devise_plan_meijer` is the reference loop (unshifts, then shifts, in
  the reference order); every docstring example gives the same operator
  list, and an unreachable target raises `NotImplementedError('Could not
  devise plan.')`.
- `(x**y).nseries(x, 0, n)` and `.series(...)` with an x-free symbolic
  exponent return `x**y` (at `oo`, `(1/x)**(-y)`), and with `logx` return
  `exp(logx*y)`. `limit` retries a product of powers after
  `powsimp(combine='exp')`, so `limit(x**p*x**(-p - 1), x, oo)` is `0`.
- `trigsimp(..., old=True)` is the reference legacy `trigsimp_old`
  (`_trigsimp`, `__trigsimp` and the `_trigpats` patterns, in
  `python/symox/simplify/_trigsimp_old.py`), with `deep`, `recursive`,
  `method` and `compare` (which prints the `futrig:` line).
- `fps(x**n*f(x), x)` with an x-free symbolic `n` at `x0 = 0` expands `f`
  and multiplies each term by `x**n`: `fps(x**n*sin(x**2), x).truncate(8)`
  is `-x**(n + 6)/6 + x**(n + 2) + O(x**(n + 8), (x, 0))`.
- The Wallis product limit of the `Product` docstring,
  `limit(simplify(W.rewrite(gamma)), n, oo)`, is `pi/2`.
- Stats: `density(X).dict` of a `Die` or `Binomial` with a symbolic size is
  the held `Density(...)` object, and `.subs(n, 4).doit()` gives the table.
  `cdf` of these is the reference `Lambda(_k, Sum(Piecewise(...)))`.
  `BernoulliProcess` values print as `B[1]`, and `joint_distribution`
  returns a `JointDistributionHandmade` (process values or indices) that
  can be called or `.pdf`-evaluated.

Reference-wrong, with checks:

- `(2*x**y).nseries(x, 0, 2)` is `O(x**2)` in SymPy, which is false at
  `y = 0` (the value is 2). Symox keeps `NotImplementedError`.
- `fps(x**(n + 1)*exp(x), x).truncate(3)`: SymPy drops `x**(n + 3)/2` but
  still claims `O(x**(n + 4))`. Symox returns
  `x**(n + 3)/2 + x**(n + 2) + x**(n + 1) + O(x**(n + 4), (x, 0))`.
- `cdf(Binomial('X', n, p, succ, fail))` with a symbolic `n` and
  `(succ, fail) != (1, 0)`: SymPy sums the count pmf from 0, which is not
  the cdf of the shifted values (for `n = 3`, `succ = 2`, `fail = 1`, its
  `cdf(3)` is 1, but `P(X <= 3) = 1/8`). Symox raises NotImplementedError.

Form-only divergences:

- An `fps` with a symbolic power prefactor prints `x**n` first in the sum.
- The symbolic `Binomial` cdf term prints `(1/2)**_ki*(1/2)**(-_ki + n)`
  where SymPy prints `2**(_ki - n)/2**_ki` (the `Mul` gap below).
- `E` of a symbolic `Die` or `Binomial` is the closed form (`n/2 + 1/2`);
  SymPy returns a `Sum`.
- `fps(x**n*f(x), x, 1)` is a Taylor expansion of the whole product; SymPy
  keeps `x**n*(...)`.

Remaining gaps:

- The native `Mul` does not multiply numeric bases that share a
  non-rational exponent (SymPy `num_exp`): `(1/4)**n*4**n` is not `1` and
  `2**n*3**n` is not `6**n`. It does not map a `1/q` base to `q` with the
  opposite exponent either (`x*(1/2)**x` is `x/2**x` in SymPy). The values
  are equal. Both rules were implemented in `finish_mul` and matched the
  reference on every probe, but the new factors reach the native inverse
  Mellin strip test in another order, so
  `meijerint_indefinite(sin(x)**3/x**3, x)` raised
  `MellinTransformStripError` instead of the reference
  `NotImplementedError('Gammas partially over the strip.')`. The change
  was reverted. That strip test needs to be made order-independent first.
- Random symbols in arithmetic print as their internal dummies
  (`N + N**2` is `_N**2 + _N`; SymPy `N**2 + N`). This predates this batch.
- `given(Tuple(X, Y), cond)` with dependent entries, multivariate
  `Laplace`, `MultivariateNormal` with `MatrixSymbol` parameters,
  `DiscreteMarkovChain` without transition data (symbolic
  `stationary_distribution`), `GaussianSymplecticEnsemble` density at a
  `MatrixSymbol`, the continuous Markov chain and Wiener/Gamma process
  queries, and the holonomic items of batch c1 (Frobenius `y0`, `to_expr`
  of `besselj`, general `to_meijerg`) still raise or stay held.
- `diff(..., (x, 5))` is still not ended with `factor_terms(signsimp(...))`
  (see batch c2), and `solve_triangulated(..., extension=True)` still
  raises NotImplementedError.

## Category-c batches 3 and 4, remaining gaps (c-rest-b) on 2026-10-23

All c3 and c4 doctest records were re-run against SymPy 1.14. Callables
whose doctests all match: c3 went from 33/41 to 33/41, c4 from 12/37 to
23/37.

Implemented:

- Pretty printing of Python-level values. The Python `PrettyPrinter` now
  composes Python objects with the reference `stringPict` and `prettyForm`
  (a verbatim copy of SymPy's `stringpict.py`). Native subexpressions come
  from the new `_core._pretty_form`, which returns the picture with its
  baseline and binding. This covers the reference `_print_seq`, list,
  tuple, dict, set and frozenset printers (keys and items sorted by
  `default_sort_key`), `TransferFunction`, `Series`, `Parallel`,
  `Feedback`, the MIMO variants, `TransferFunctionMatrix` and
  `DiagramGrid`. Objects with a `_pretty` method are tried first, as in
  the reference dispatch.
- Categories:
  - Category values have the reference `sort_key` (Object like a named
    atom, morphisms as Basic, `CompositeMorphism` through a `Tuple`).
  - `Diagram.premises`, `Diagram.conclusions` and `Diagram.hom` return
    sets that print like SymPy's `FiniteSet`/`EmptySet` and are ordered
    like the reference `Dict`.
  - `Morphism(A, B)` raises SymPy's NotImplementedError.
  - `DiagramGrid`, `ArrowStringDescription`, `XypicDiagramDrawer` and
    `xypic_draw_diagram` are now a verbatim port of
    `sympy/categories/diagram_drawing.py`. Groups use a frozenset with
    the FiniteSet sort key. The previous simplified layout differed (for
    example a 3x1 grid instead of SymPy's 2x2 triangle). All 153
    reference doctests of the module pass. The errors are SymPy's too:
    IndexError for an empty diagram, KeyError for groups that miss an
    object.
- `cancel(f, expand=False)` follows the reference: unexpanded factors
  become generators, and `dmp_cancel(include=True)` runs over fresh
  symbols. Before, the option was ignored and the result was expanded.
  This also fixes `TransferFunction.simplify`,
  `TransferFunctionMatrix.simplify` and `MIMOFeedback.doit`. For the
  same reason, `parallel_poly_from_expr(..., expand=False)` now finds
  the generators without expanding.
- `jones_2_stokes` uses the reference formula (`Abs(e)**2`), without the
  extra `trigsimp`.
- `Add(v1, v2)` of vectors returns a `VectorAdd`, so its kind is
  `VectorKind(NumberKind)`. `Add(N.i, 1)` still raises TypeError.

Intentional divergences:

- `combine_kronecker(A⊗B + C⊗D)`: SymPy returns `(A + C)⊗(B + D)`. That is
  wrong in general, so symox keeps the sum unchanged.
- `ufuncify`: the SymPy doctest fails in this environment (f2py "return
  arrays must be of ArrayType"), while symox's path works.
- `preview_diagram` loads `\xyoption{all}` and passes `filename`,
  `outputbuffer` and `dvioptions` through to `preview`. SymPy omits the
  option, so its LaTeX run fails on curved arrows and loops. SymPy also
  hands those three options to `DiagramGrid` as ignored hints, so
  `viewer="file"` cannot work there.

Still unsupported or different (none of them returns a wrong value):

- c3:
  - `galois_group` of degree 6 over an AlgebraicField raises
    NotImplementedError.
  - `prime_decomp` over `QQ.algebraic_field((T, theta))` and
    `has_dups(Set())`: resolved by c-doc (2026-10-23).
  - `julia_code` and `octave_code` of a user function applied to a
    Matrix, and `Integral` of a `MatrixSymbol` (`MatrixKind` doctest),
    raise TypeError or ValueError: native expressions cannot hold
    matrices. (`Integral` of a `MatrixSymbol` is resolved by c-doc,
    2026-10-23.)
- c4:
  - `FiniteSet` of category values and the `CoordinateSymbol` doctests:
    resolved by c-doc (2026-10-23).
  - `shape`/`ArrayKind` of `Integral(Array)` and the array-expression
    doctests: resolved by c-doc (2026-10-23).
  - `phase_margin` returns an equal closed form written differently,
    because of the native `arg`/`atan` simplification.
  - `simplify_gpgp(ps*qs*qs)`: resolved by c-doc (2026-10-23).
  - `couple` and `uncouple` with symbolic spins: resolved by c-doc
    (2026-10-23).
  - `ImplicitRegion.rational_parametrization`: resolved by c-doc
    (2026-10-23).

## Odd/even heads on noncommutative sums, DiracDelta real-value check (nc-segv, 2026-10-23)

`erf(B*A - 2*x)` with noncommutative `A`, `B` crashed the interpreter
(SIGSEGV, unbounded native recursion). The negation of a noncommutative sum
`z` stays `Mul(-1, z)`, as in SymPy, because a number is not distributed over
a noncommutative sum. The native sign test (`ops::could_extract_minus_sign`)
broke ties between terms with its own structural rule, and reported both `z`
and `Mul(-1, z)` as negative. The rewrite `f(z) -> -f(-z)` then flipped back
and forth. Implemented:

- Noncommutative sums and products now use SymPy's exact
  `could_extract_minus_sign` (term count, then `sort_key` of `z` against
  `-z`). That rule is antisymmetric, which fixes erf, erfc, erfi, erf2, Si,
  Shi, fresnels, fresnelc, besselj, besseli and the trigonometric and
  hyperbolic heads (60 probe cases went from a crash to SymPy's output).
  Commutative arguments keep the previous rule, and the printers keep the
  syntactic test.
- Si, Shi, fresnels and fresnelc use SymPy's `extract_multiplicatively(-1)`
  on a noncommutative sum: `Si(-A*B - x)` is `-Si(x + A*B)`, as in SymPy.
- DiracDelta and Heaviside port SymPy's `fuzzy_not(im(arg).is_zero)` check
  for noncommutative arguments. `im(X)` and `re(X)` of a noncommutative `X`
  are nonzero (SymPy's rule `extended_real -> commutative`), and products of
  such parts follow `Mul._eval_is_zero`. The error messages are SymPy's: the
  DiracDelta one is wrapped like `filldedent`.

Intentional divergences:

- When a noncommutative sum `z` could itself extract a minus sign (for
  example `-A*B - x`), SymPy 1.14 raises RecursionError for sin, cos, tan,
  erf, asin, atan, sinh, the Bessel functions, Abs, sinc and others. Symox
  returns the equivalent `-f(-(-x - A*B))` for an odd head, `f(-(-x - A*B))`
  for an even head, and `Abs(-x - A*B)`.

Remaining gaps:

- SymPy treats `im`, `re`, `Abs`, `sign`, `arg`, `DiracDelta`, `Heaviside`
  and `hyper` of a noncommutative argument as commutative, because they are
  extended real or because of their own assumption rules. Symox keeps them
  noncommutative. So factor order differs in products. For example, the
  DiracDelta message of probe 628 prints
  `atan(p)**(3/2)*im(t*sinh(y)*Derivative(g(x, y), x, y)*B)`, and SymPy
  prints `im(...)*atan(p)**(3/2)`.
- The general assumption system does not yet derive `extended_real=False`
  (and every fact that implies it) from `commutative=False` for
  noncommutative products and sums (`(x*A).is_extended_real` is None, while
  SymPy gives False). The real-value check uses its own port, so it does not
  depend on this. Sums such as `DiracDelta(A + I)` still do not raise (SymPy
  raises), because `Add._eval_is_zero` settles them through those facts.

## sympy.stats follow-up (stats-rest, 2026-10-23)

This batch handles the stats gaps listed for batches c1 and c2. Doctest
probes over ContinuousMarkovChain, DiscreteMarkovChain (with
`stationary_distribution`), GammaProcess, WienerProcess,
GaussianSymplecticEnsemble, Laplace, MultivariateNormal, density,
dependent, independent, Moyal, DiscreteRV and given went from 137/171 to
163/171 matches with SymPy 1.14. Tests are in
`tests/test_stats_rest_parity.py`.

Implemented:

- Random symbols in arithmetic are a `RandomAtom` Dummy that prints as
  its name in str, latex and pretty (`N**2 + N`, `N^{2} + N`), and that
  keeps `pspace` and `distribution`.
- `given(Tuple(...), cond)` with dependent entries conditions the joint
  space, as SymPy does.
  - Finite spaces filter and renormalize the joint outcomes.
  - Independent continuous variables conditioned on a linear `Eq` eliminate
    one variable on the hyperplane and normalize by an exact integral.
  - The new entries share one space, so `dependent` is True.
  - Results are closed forms. SymPy leaves `erfc` residue, and leaves an
    unevaluated `Integral` for `E(B)` of `given(Tuple(X, Y), Eq(2*X + Y, 1))`.
    The correct values are 2/5 and 1/5.
  - A continuous inequality condition raises SymPy's
    `NotImplementedError('Multivariate Inequalities not yet implemented')`.
- `Laplace(name, mu_list, sigma_matrix)` dispatches to
  `MultivariateLaplace`, as SymPy does.
- DiscreteMarkovChain without transition data:
  - The transition is `MatrixSymbol('_T', n, n)` over `Range(0, n)`.
  - A symbolic `stationary_distribution` is
    `Lambda((wm, T), Eq(wm*T, wm))`, or a `ConditionSet` with
    `condition_set=True`. The native `Lambda` takes only scalar
    signatures, so this is a small Python `MatrixLambda`.
  - A numeric chain uses SymPy's linsolve system, so a reducible chain
    gives the reference parametric row (`1 - tau0, tau0`).
- Markov chains at float or symbolic times:
  - The kernel is SymPy's transition function: `Q*exp(t*D)*Q**-1` (CTMC)
    or `Q*D**t*Q**-1` (DTMC) when diagonalizable, otherwise `exp(t*G)`.
  - Symbolic-time queries port `_symbolic_probability`: a `Piecewise` on
    `a > c`, with the held `Probability` otherwise. Gt, Ge, Lt and Le use the
    reference sum from `next_state + upper`.
  - The kernel entry at symbolic states is a `Piecewise` over the states,
    where SymPy has a symbolic matrix index. The native explicit matrix
    takes only integer indices.
  - The values after `subs` match the doctests (0.4002723175, 0.6832579186,
    0.3096, 0.64705).
- Wiener and Gamma processes under `Contains(t, interval)` conditions use
  the reference convention already implemented for Poisson. The value at t
  is the increment over the interval, and disjoint intervals are
  independent.
- `GaussianSymplecticEnsemble` density at an n x n `MatrixSymbol` is the
  normalized closed form, for example `32*exp(-2*Trace(X**2))/pi**3` for
  n = 2. SymPy holds the normalizer as `Integral(exp(-2*Trace(_H**2)), _H)`.

Divergences kept, where SymPy is wrong:

- `P((X(t1) > 2) | (X(t1) < 4), Contains(t1, Interval.Ropen(2, 4)))` for
  a Wiener process is 1. SymPy adds the two overlapping events and returns
  `-erf(1)/2 + erf(2)/2 + 1`, as in the Poisson `Or` divergences already
  listed.
- Poisson `P(X(2) <= 3, X(1) > 1)` is
  `45*exp(-3)/(2*(-4 + exp(3)))`, the exact conditional probability.
  SymPy returns `5*exp(-3)`.

Remaining gaps:

- `MultivariateNormal` (and T and Laplace) with `MatrixSymbol` parameters
  of symbolic dimension raise
  `NotImplementedError('multivariate distributions of symbolic dimension
  (MatrixSymbol parameters) are not supported')`. SymPy holds
  `(exp(matrix)/sqrt(...))[0, 0]`, but the native exp of an abstract
  matrix is not implemented.
- `Str` states of a Markov chain are rejected, because the native `Eq` does
  not take `Str`.
- The `str` printer does not parenthesize a compound matrix parent of a
  `MatrixElement`: it prints `T**5[1, 2]` where SymPy prints `(T**5)[1, 2]`.
  The fix is in Rust (`printer.rs`) and was not made in this batch.
- After `subs`, the held `Probability` branch of a symbolic-time query keeps
  `C(a)`, because the random atom does not carry its time through `subs`.
  SymPy prints `Probability(Eq(C(1), 2))`.
- Term order: a random symbol sorts as a plain symbol, so it prints
  `X + x` where SymPy prints `x + X`. In SymPy, RandomSymbol's class_key is
  `(5, 0, ...)`. The native print order keys every symbol as `(2, 0, ...)`.
- `simplify(cdf(Moyal)(z))` is `erfc(...)` where SymPy gives `1 - erf(...)`.
  This is a simplification form difference.


## Small remaining gaps: order, containers, division options, orthopolys

Implemented (regression tests in `tests/test_small_rest_parity.py`):
- Add ordering of `Sum`, `Product`, `Integral`, `Limit`, `Subs` and
  `Derivative` terms in str, LaTeX and pretty. As in `Expr.as_terms`, a
  closed `Integral`/`Sum`/`Product` that is a number is folded into the
  numeric coefficient, and a factor whose commutativity is unknown (`Limit`)
  goes to the noncommutative part.
- `re`, `im`, `Abs`, `arg` and `sign` of noncommutative arguments are
  commutative, as in SymPy. This removes the order-only mismatches of the
  previous probe.
- Host `NumberSymbol` subclasses (`hbar`) report `is_number` True and order
  as numeric coefficients (`x**2 + y + 2*hbar`).
- LaTeX `Limit`: `\lim_{x \to 0^+} \frac{1}{x}`, with `\left(...\right)`
  around Add/Mul bodies only. LaTeX `Subs` always uses `\substack{ ... }`.
- Pretty printing of Python lists, tuples (with the one-element trailing
  comma), dicts and sets (sorted by `default_sort_key`), frozensets, empty
  `set()`, and `Tuple` without a `Tuple` prefix. Before, these fell back to
  their str form. `Integral` pretty printing already matched in every probe
  case.
- `div`/`rem`/`quo`/`exquo`/`pdiv`/`prem`/...: `greedy` is accepted
  (`greedy=False` turns numeric factors into coefficients, so the domain
  becomes EX). `symmetric` is accepted. The requires/excludes checks of
  `polyoptions` are ported, in keyword order, with SymPy's `OptionError`
  messages (for example `'symmetric' option is only allowed together with
  'modulus'`).
- `jacobi` and `assoc_laguerre` with a negative number degree raise
  SymPy's `ValueError` (`Cannot generate Jacobi polynomial of index -1/2`,
  `The index n must be nonnegative integer (got -1/2)`). `laguerre` with a
  non-integer degree raises `Error: n should be an integer.`, and with a
  negative Float degree gives `exp(x)*laguerre(-n - 1, -x)`. The symmetric
  Jacobi branches (`a == b`, `b == -a`) come first, as in SymPy. A
  substitution that rebuilds such an application raises too. Held
  applications (`evaluate=False`) still evaluate numerically through
  mpmath.

Probe counts (mismatching output lines against SymPy 1.14, before is the
main build): order and LaTeX 19 -> 2, where both are the pretty `Sum` +
`Limit` form below; re/im and hbar 15 -> 0; division options 80 -> 0;
orthopolys 116 -> 26, where all 26 are the positive non-integer degree
divergence below; container pretty printing 50 -> 0.

Intentional divergences:
- `jacobi` and `assoc_laguerre` with a positive non-integer number degree:
  SymPy's `jacobi_poly`/`assoc_laguerre_poly` truncate the degree with
  `int(n)`, so `jacobi(1/2, a, b, x)` returns 1, which is mathematically
  wrong. Here the application stays unevaluated and evaluates numerically
  through mpmath.
- `gegenbauer` with a negative non-integer number degree: SymPy returns 0,
  which holds only for negative integer degrees. Here it stays
  unevaluated. This also applies to `jacobi(-1/2, a, a, x)`, which SymPy
  reduces through `gegenbauer`.

Remaining unsupported:
- `elliptic_pi` outside the native real kernel domain (`n >= 1` or
  `m >= 1`) evaluates through mpmath with the certified refinement loop,
  which calls `mpmath.ellippi` at several precisions. mpmath's RJ
  quadrature is very sensitive to the working precision:
  `elliptic_pi(2, 7/2).evalf(50)` takes about 120 s here and 0.3 s in
  SymPy, which makes a single call. Fixing this needs a native
  multiprecision Carlson RJ, or branch-exact transformations to the real
  domain.
- With a modulus, `symmetric=False` raises `NotImplementedError`, because
  there is no nonsymmetric finite-field display. `composite`, `sort`,
  `wrt` and `split` are still `NotImplementedError` in division.
- `Sum` and `Limit` pretty forms inside an Add still print their str form
  on this branch. The 2D forms come from the separate pretty printing port.

## Numeric bases sharing an exponent (wave 4, numbase)

Implemented in the Rust `Mul` (`ops.rs::finish_mul`), following SymPy 1.14
`Mul.flatten`:
- The `num_exp` rule: powers with a finite nonzero number base and a
  non-rational exponent are gathered by exponent when the base is positive or
  the exponent is an integer. The bases sharing an exactly equal exponent are
  multiplied and rebuilt as one power: `2**n*3**n` is `6**n`,
  `(1/4)**n*4**n` is `1`, `Float(2.5)**n*2**n` is `5.0**n`,
  `(-2)**k*(-3)**k` is `6**k` for an integer `k`, while `(-2)**n*3**n` stays
  as it is for a general `n`. The combined power is not merged again with
  the other factors, as in SymPy: `2**n*3**n*6**(-n)` is `1` only through
  the binary evaluation order.
- `Pow.as_base_exp` for a unit-numerator rational base inside a product:
  `x*(1/2)**n` is `x/2**n`, while `(1/2)**n` alone stays as it is.
- The gammasimp internal `factor` no longer expands each factor base before
  factoring. SymPy runs `together` first. Expanding first changed the
  result once `2**s*6**s` combined into `12**s`.
- The SymPy sequence candidate of the native `simplify` now applies
  `combsimp` to gamma and combinatorial functions before its `powsimp`,
  `cancel` and `factor_terms` steps, as SymPy does. With the combined
  numeric powers, this keeps `meijerint_indefinite(sin(x)**3/x**3, x)`
  raising SymPy's `NotImplementedError('Gammas partially over the strip.')`.

Probes against SymPy 1.14 (`str` of 86 products and sums of numeric-base
powers): 54 differed before, none after.

Consequence: `solveset(2**x*3**x - 6, x, Complexes)` is now the equation
`6**x - 6`, as in SymPy, and falls under the documented divergence for
`b**x` with a positive rational base: the complete lattice
`2*n*I*pi/log(6) + 1`, where SymPy returns a `ConditionSet`.

Remaining gaps: the native `simplify` is still a candidate search, not
SymPy's chain of rewrites. Forms of other expressions with combined numeric
powers can still differ where the native candidates win. For example, the
simplified Mellin integrand of `sin(x)**3/x**3` keeps `2**s` where SymPy
keeps `6**s/3**s`. That form already differed before this change.

## Holonomic, multiple derivatives, triangulated extension and parametric systems (holo-rest, 2026-10-23)

Ported from SymPy 1.14. The regressions are in `tests/test_holo_rest_parity.py`.

- Holonomic functions (`sympy/holonomic/holonomic.py`):
  - Frobenius initial conditions `y0 = {p: [c0, c1, ...]}` are accepted.
    `is_singularics()` returns True for a dict, False for a list and None
    without initial conditions.
  - Sums and products combine singular conditions at the same point.
    `diff`, `integrate`, `change_ics` and `__neg__` follow the reference,
    and `to_sequence` returns the reference 3-tuples `(sequence, p, n0)`.
  - `series` and `to_sequence` follow the reference at singular points: the
    indicial equation and `_frobenius`, with the `C_i` constants left
    undetermined.
  - `to_hyper` is the reference algorithm. `to_expr` is
    `hyperexpand(to_hyper()).simplify()`, so the Bessel equations give
    `besselj(1, x)`, `besselj(2, x)/2` and `besselj(0, x)`.
  - `to_meijerg` sums `_hyper_to_meijerg` over the hypergeometric terms.
  - `expr_to_holonomic` gives singular conditions for algebraic functions
    at a branch point: `sqrt(x**3 + x)` and `x**(3/2)*exp(x)`.
  - `NotPowerSeriesError` and `NotHyperSeriesError` exist, with the
    reference messages.
- A derivative of total order greater than 1 ends with
  `factor_terms(signsimp(...))`. This applies to `diff`, `Expr.diff`,
  `Derivative(..., evaluate=True)` and `Derivative.doit()`, including held
  derivatives nested in a larger expression.
  `diff(sqrt((x + 1)**2 + x), (x, 5)).count_ops()` is 30.
  - The earlier attempt broke the undetermined coefficients ODE solver.
    The cause was the native `expand_mul`: it did not distribute a held
    `-(a*sin(x) + b*cos(x))` term inside a sum. It now does for any numeric
    coefficient, not only for noncommutative sums.
- `solve_triangulated(..., extension=True)` is the reference algorithm. The
  roots come from `all_roots(radicals=False)`, and each partial solution
  carries its `QQ<CRootOf(...)>` field. This needed two fixes:
  - Polys whose coefficients are `CRootOf` objects named after a generator
    (`CRootOf(x**2 - 2, 0)` in a polynomial in `x`) are now built over
    placeholders in `Poly._from_dict` and `Poly.eval`. This also fixes
    `Poly(y**2 - CRootOf(x**2 - 2, 0), y, extension=True).all_roots()`.
  - The native `algebraic_field` rejects a rational generator. A rational
    root therefore keeps the current domain (the reference builds `QQ<0>`).
- `solve(sec(k*x) - r)` and `solve(csc(k*x) - r)` for `k = 2, 3`, when the
  inverse trig value has no closed form. The quadratic in `u = exp(I*k*x)`
  gives the `k` roots of `t**k = u`, as `-I*log(u)/k` and
  `-I*log(w*u**(1/k))`. `sec(2*x) - 3`, `2*sec(2*x) - 5` and
  `csc(3*x) + 4` print as in the reference.
- `nonlinsolve` of a polynomial system with external parameters, when the
  exact parametric solver cannot represent the parameter conditions, now
  returns the reference `_handle_zero_dimensional` answer: generic values
  of the parameters. `nonlinsolve([a*x**2 - 1, y - x], [x, y])` matches.
  `nonlinsolve([a*x], x)` is `{(0,)}` as in SymPy; it used to raise. The
  `a = 0` component is still obtained by including `a` among the
  unknowns. An inconsistent generic basis does not prove that there are no
  solutions, so it still raises NotImplementedError.

Probe comparison (SymPy 1.14; values hard-coded in the tests):

| Probes | Before | After |
| --- | --- | --- |
| Holonomic | 8/40 | 37/40 |
| Multiple derivatives | 4/13 | 13/13 |
| `solve_triangulated` | 0/12 | 6/12 (the 6 others differ only in the `extension=False` OptionError) |
| Item-4 solve probes | 3/11 | 6/11 |

Intentional divergences:

- `expr_to_holonomic(exp(x) + sin(x)).to_hyper()` returns the correct sum
  of hypergeometric terms. SymPy raises NotHyperSeriesError.
- `HolonomicFunction(x*Dx - 1, x, 0, [1]).series()` keeps the ValueError,
  because the initial conditions are inconsistent. SymPy returns `C_0*x`.
- Some products and sums give the annihilator as a scalar multiple of the
  reference one, because symox normalizes cyclic vectors differently. The
  initial conditions match: `(h + h)` is `(-1) + (2*x)*Dx` and SymPy gives
  `(-1/2) + (x)*Dx`.
- In some `to_sequence` results symox determines a value such as `u(1) = 1`
  where SymPy leaves the constant `C_1`.
- `solve_triangulated(..., extension=False)` returns the ground-domain
  roots. SymPy's `build_options` rejects the explicit `False`.
- The denominator-zero roots of `nonlinsolve` stay dropped. SymPy also
  returns `(0, 0)` for `[1/x + 1/y - 1, x - 2*y]` and `(0, 1)` for
  `[x/(y - 1) - 1, y**2 - 1]`, where a denominator vanishes.

Remaining gaps:

- Parametric `nonlinsolve` and `solve` with symbolic coefficients agree
  with SymPy in value, but not in form. Symox prints
  `(1 - sqrt(-4*a*b + 1))/(2*a)`; the `solve_poly_system` back-substitution
  in SymPy prints `-(-sqrt(-4*a*b + 1)/2 - 1/2)/a`.
- `sec(4*x) - 3` and `sec(x/2) - 3` raise NotImplementedError. SymPy
  solves the first; the second raises in both, with a different message.
- With a symbolic right-hand side, `sec(2*x) - a` is solved, but it prints
  in an unsimplified form.
- When the inverse trig value has a closed form (`csc(2*x) - 2`,
  `sec(2*x) + 1`), only the principal family is returned. SymPy also lists
  `-I*log(-sqrt(sqrt(3)/2 + I/2))`-type roots, and both `-pi/2` and `pi/2`.
- Holonomic, out of scope: the `erfc` rewrite, `integrate` of singular
  conditions returning `nan`, float domains
  (`expr_to_holonomic(1.2*cos(3.1*x))`), symbolic or negative powers, and
  the reference expand canonicalization of rational coefficients.

## elliptic_pi numerical evaluation (ellippi, 2026-10-23)

SymPy 1.14 evaluates `elliptic_pi` with one `mpmath.ellippi` call at the
target precision. When the Carlson `R_J` arguments leave the region where
Carlson's duplication is valid (`m > 1`, `n > 1`, most complex arguments),
mpmath integrates first along `[0, N + i*margin]` with `quadsubdiv`. That
quadrature passes close to the integrand's singular points. Symox's evalf
evaluates at two or more precisions (`refined` in `evalf_mp.rs`) and checks
that the results agree, so `elliptic_pi(2, 7/2).evalf(50)` took about 120 s.
SymPy took 0.3 s, but timed out on 7 of the evd probe evaluations. At 30
digits SymPy also returns wrong trailing digits, for example `...923255`
instead of `...933061` for `elliptic_pi(2, 7/2)`.

`symox._ellippi` follows mpmath's `ellippi` and `elliprj` step by step: the
same special values, the quasi-periodicity reduction, `sum_accurately` and
the Carlson tail. Only the initial integral changes:

- When mpmath's margin is `+i` or `-i`, the path is moved to the half-line
  `i*sigma*[0, oo)`. The region between the two paths contains no singular
  point and meets no branch cut. The integral is computed numerically on
  `[0, i*sigma*H]`. The rest is `exp(-3*i*sigma*pi/4)*R_J(H - i*sigma*w)`,
  where every argument has a positive real part.
- Otherwise mpmath's own segment is kept.

In both cases the segment is split until every piece is at least four times
its length away from the singular points. The `t = 0` endpoint singularity
is removed by `t = b*u**2`. Each piece is then integrated with cached
Gauss-Legendre rules and an error check. Where mpmath applies Carlson's
algorithm directly, symox calls it with the same arguments and gets the same
bits.

Results on the 342 evaluations of the elliptic_k, e, f and pi probe at 15,
30 and 50 digits:

- Symox takes 11 s in total, against 220 s for SymPy. SymPy times out (120
  s) on 11 elliptic_pi evaluations. Before the change, symox had 7 timeouts
  on the evd probe.
- Every value is within one unit in the last digit of a 150-digit reference.
  That reference agrees with mpmath at 100 digits to more than 60 digits,
  which is the limit of mpmath's own quadrature.
- SymPy returns wrong trailing digits on 23 elliptic_pi evaluations: about
  25 correct digits out of 30, and 35 to 42 out of 50. Symox keeps the
  correct value there.

Remaining gaps:

- The 50-digit `elliptic_pi(2/3, 2 + 3*I/2)` takes 0.5 s, against 0.3 s for
  SymPy, because `refined` evaluates twice.
- If a singular point lies on the integration path (bisection depth 64),
  the integral falls back to mpmath's `quadsubdiv`, with mpmath's speed and
  accuracy.

## Commutative heads of noncommutative arguments, nc assumptions (nc-rest, 2026-10-23)

Probe of 78 functions on 16 noncommutative arguments (1248 cases): 174
mismatches before this branch, 84 after merging feat-small-rest, 31 now
(29 SymPy RecursionError cases, and csc/csch below). Regression tests in
`tests/test_nc_rest_parity.py`, run in a child interpreter.

Implemented:

- Heads that SymPy declares commutative at class level are commutative for
  any argument (`decompose::commutative_head`): re, im, Abs, sign, arg,
  DiracDelta, Heaviside, Max, Min, hyper, meijerg, KroneckerDelta,
  LeviCivita and the ntheory heads (partition, totient, reduced_totient,
  mobius, divisor_sigma, udivisor_sigma, legendre/jacobi/kronecker_symbol,
  primenu, primeomega, primepi). Products, expansion and distribution of a
  rational over a sum treat them as scalars:
  `hyper((1,), (2,), A)*hyper((1,), (2,), B)*B`, `2*re(A) + 2*x`. This
  supersedes the nc-segv gap on the probe-628 DiracDelta message, which is
  now SymPy's.
- Every fact whose closure contains `commutative` (real, extended_real,
  complex, zero, nonzero, integer, rational, positive, imaginary,
  algebraic, ...) is False for a noncommutative symbol, sum, product or
  function: `(x*A).is_extended_real`, `(A + x).is_zero`,
  `(A*B - B*A).is_zero` are False. `KroneckerDelta(A, x)` is 0.
- `re`/`im`/`Abs`/`sign` of a noncommutative argument are nonzero, and
  re/im are not integers (`re._eval_is_algebraic`).
- A commutative sum of such parts follows Add's `_fuzzy_group` on
  integer/algebraic: `(im(A) + 1).is_zero` is False, `(x + im(A)).is_zero`
  is None. So `DiracDelta(A + I)` and `Heaviside(x*A - I)` raise SymPy's
  "Function defined only for Real Values" errors.
- `Eq`/`Ne` of a noncommutative difference stay unevaluated (SymPy's
  `is_eq` only returns False for a commutative difference).
- floor/ceiling of a complex noncommutative argument keep it whole except
  for integer terms (`floor(x - I*A*B)`, `1 + floor(A)`).
- `arg` uses `factor_terms(arg).as_coeff_Mul()` as SymPy does:
  `arg(-(x + A*B))`, `arg(x*(y + 1))` for `arg(2*x*y + 2*x)`.
- str `MatrixElement` parenthesizes its parent below atom precedence:
  `(M**5)[1, 2]`, `(M*N)[1, 2]`, `(M**(-1))[0, 1]`, `(-M*N)[0, 1]`
  (LaTeX and pretty already matched).

Intentional divergences:

- A noncommutative Pow keeps the native answers instead of False: SymPy
  decides them only through Pow's own handlers (`(A**2).is_positive` is
  None, `(A**-1).is_real` False), and Symox answers None for some facts
  SymPy settles to False there (`(A**2).is_integer`, `sqrt(A).is_real`).
- SymPy prints `csc(I*z)` as `(-I)*sinh(z)**(-1)` (Mul.flatten keeps the
  evaluated `I**-1` as a nested `Mul(-1, I)`); Symox prints
  `-I*sinh(z)**(-1)`. The values are equal; this belongs to Mul numeric
  combining.
- The SymPy RecursionError cases of the nc-segv section are unchanged.

Remaining gaps:

- Derived facts that need the closure of several facts are not computed
  for these heads: `re(A).is_noninteger`, `Abs(A).is_negative` and
  `Abs(x).is_negative` are None (SymPy True/False), and `sign(z).is_complex`
  is None (SymPy True, also for commutative `z`).
- `floor(x + 2*I)` (commutative) is split through `as_real_imag` into
  `I*(floor(im(x)) + 2) + floor(re(x))`; SymPy's RoundFunction.eval only
  splits integer and Gaussian-integer terms (`floor(x) + 2*I`).
- `MatrixElement(M + N, i, j)`, `MatrixElement(M.T, i, j)` and
  `MatrixElement(HadamardProduct(M, N), i, j)` evaluate on construction in
  Symox and stay unevaluated in SymPy. Conversely `(M*N)[i, j]` and
  `(M**5)[1, 2]` stay unevaluated in Symox, while SymPy expands them into a
  `Sum` over dummy indices.

## c-final: harmonic, loggamma, implicit solve, lambdify namespace, autowrap helpers

Implemented (SymPy 1.14 behaviour), closing the five undocumented "c"
objects of the final wave-4 triage:

- `expand_func(harmonic(u + p/q))` for a rational argument with `u >= 0`
  gives the reference closed form `q*Sum(1/(q*k + p), (k, 0, u)) +
  2*Sum(cos(2*pi*p*k/q)*log(sin(pi*k/q)), (k, 1, floor((q - 1)/2))) -
  pi*cot(pi*p/q)/2 - log(2*q)` (so `harmonic(1/3)` and `harmonic(25/7)`
  reach their values through `doit`), and `expand_func(harmonic(n - k))`
  pulls out negative integer offsets like positive ones. Negative rationals
  and orders `m != 1` stay unexpanded, as in SymPy. `harmonic(1/3)` and
  `harmonic(1/3).doit()` stay unevaluated in both.
- `expand_func(loggamma(n + p/q))` gives the reference sum
  `loggamma(p/q) - n*log(q) + Sum(log((k - 1)*q + p), (k, 1, n))`, with
  `pi*I*n - Sum(log(k*q - p), (k, 1, -n))` for negative `n`.
- `series(loggamma(z), x, 0, n)` at a zero of the argument expands
  `log(gamma(z))` like `loggamma._eval_nseries`
  (`-log(x) - EulerGamma*x + pi**2*x**2/12 - ...`); at a negative integer it
  still raises `PoleError` (with the message `loggamma(x - 1)` instead of
  SymPy's `Cannot expand loggamma(_x - 1) around 0`).
- `solve(eq)` without symbols reaches the masked, inverted and unrad
  fallbacks when the equation has one free symbol, so the radical equations
  `x**(1/3) - x**(1/5) + 1/7` (two `CRootOf(...)**15` roots),
  `(x**3 - 3*x**2)**(1/3) + 1 - x` (no solution: the unrad candidate 1/3
  fails the check) and `sqrt(x + 1) - sqrt(x) - 1` are solved exactly.
- `solve(..., implicit=True)` masks every Pow, function application and
  other non-invertible object that contains an unknown (re/im included,
  Piecewise traversed), solves the masked system and restores the masked
  objects, as SymPy's `_solve` does (`solve(x + exp(x), x, implicit=True)`
  is `[-exp(x)]`). `implicit=False` is plain solve. Non-Symbol unknowns with
  `implicit=True` raise `NotImplementedError`.
- `reduce_inequalities` with univariate inequalities in several symbols
  (`[0 <= x - 1, y < 3]`) reduces each symbol's group separately and
  returns their `And`, as `_reduce_inequalities` does, including symbols
  that are not in the requested list.
- `symox.testing.pytest.ignore_warnings` is re-exported as in SymPy, and
  the compiled lambdify function has a `__globals__` mapping assembled on
  first access from SymPy's module tables (`f.__globals__['sin'] is
  numpy.sin` by default, translations such as `ceiling -> ceil` for math,
  the dict modules taking precedence, the expression symbols, `builtins`
  and `range`); the Python-level results (user-function dicts, implemented
  functions, structured and mpmath paths) are rebuilt with that mapping as
  their globals. The native evaluator does not execute in that namespace,
  so modifying it does not change the function (SymPy's would), and a
  user-dict function receives the evaluator's float (`mysin(1.0)` where
  SymPy's generated code passes the int `1`).
- `autowrap(..., helpers=("name", expr, args))` accepts a single triple
  like SymPy, besides an iterable of triples.

Remaining gaps: the native `CompiledFunction` still has no `__name__`,
`__code__`, signature or source (`inspect.getsource` raises `TypeError`, see
the lambdify notes above); `limit(harmonic(n, m), n, oo)`; `solve(sqrt(x) -
1 + x)` without symbols returns the equivalent `(-1/2 + sqrt(5)/2)**2`
where SymPy prints `3/2 - sqrt(5)/2` (with `x` given both agree).

## Documented c objects of the final triage (c-doc, 2026-10-23)

Implemented (doc examples now match SymPy 1.14):

- `has_dups`: the abstract `Set(*args)` can be built (a bare node, as the
  reference `Basic`), and `has_dups` returns False for `dict`, `set`,
  `Dict` and `Set` without iterating, as the reference does.
- `nonlinsolve` of a polynomial system with Float coefficients: the
  rationalized system is solved exactly and its isolated solutions are
  evaluated (`{(-0.5, -0.5), (0.5, 0.5)}` for `x**2 - 0.25, y - x`, as
  the reference). Positive-dimensional Float systems still raise
  NotImplementedError ("exact coefficients").
- `ImplicitRegion.regular_point` / `rational_parametrization`: conic regular
  points use a private port of the reference parametric solution of a
  homogeneous ternary quadratic in normal form
  (`vector/_ternary_reference.py`), so the circle example gives
  `(4*t/(t**2 + 1), 4*t**2/(t**2 + 1) - 2)` and ellipses such as
  `x**2 + 2*y**2 - 3` get the reference point instead of
  NotImplementedError. The public `diophantine` still refuses this family,
  because the reference two-parameter answer is not the complete integer
  solution set. The port covers a positive `x**2` coefficient, the only
  case conic regular points use. The Float `reg_point` example works through
  the `nonlinsolve` change above.
- `simplify_gpgp`: the reference algorithm with `sorted_components`
  (`ps*qs*qs` gives `GammaMatrix(-L_0)*p(L_0)*q(L_1)*q(-L_1)`).
- `QQ.algebraic_field((T, theta))` with a Symbol `theta` (`prime_decomp`
  doctest): the field is generated by the last `CRootOf` of `T` under the
  alias `theta`, with `T` kept as the defining polynomial, so
  `primes_above`, `discriminant`, `galois_group` and the prime printing in
  the variable of `T` match. A reducible `T` raises NotImplementedError.
  As for other aliased fields, `to_sympy` writes the `CRootOf` instead of
  `theta`.
- `couple` / `uncouple` with symbolic spins return the reference held `Sum`
  of CG terms (scalar arithmetic and `subs` supported; `doit()` raises the
  reference `ValueError('Coefficients must be numerical')`; other
  arithmetic raises NotImplementedError). `TensorProduct.subs` now reaches
  its factors.
- `Category` / `Diagram.subdiagram_from_objects`: `FiniteSet` holds
  category objects and diagrams, so `K.commutative_diagrams ==
  FiniteSet(d)` is True. The printed `Diagram` repr inside the set still
  differs.
- `CoordinateSymbol`: `Lambda` accepts a Matrix body, so `CoordSystem`
  relations can be given; `x.coord_sys`, `x.index` and `rewrite(Pol)`
  match.

- `Integral` of an explicit array, Matrix or `MatrixSymbol` (`shape`,
  `ArrayKind`, `MatrixKind` doctests) is a held shaped integral with the
  reference `.kind` and `.shape`; `doit()` integrates explicit entries one
  by one and a `MatrixSymbol` constant in the variable (`x*A`), and keeps
  anything else unevaluated.
- Array expressions: `tensorproduct`, `tensorcontraction`,
  `tensordiagonal` and `permutedims` of `MatrixSymbol`s and array
  expressions build symbolic array expressions; `PermuteDims` composes
  nested permutations (`.expr`, `.permutation` as the reference);
  `sort_args_by_name` is the reference algorithm; `convert_array_to_matrix`
  turns contractions of rank-2 tensor products into matrix products,
  transposes, `Trace` or an `ArrayTensorProduct` of products;
  `convert_indexed_to_array` is the reference port for `Indexed`,
  `MatrixElement`, `Mul` and `Sum`, with `first_indices`.

Intentional divergences:

- `combine_kronecker(KroneckerProduct(A, B) + KroneckerProduct(B.T, A.T))`
  stays unchanged. The reference `KroneckerProduct(A + B.T, B + A.T)` drops
  the cross terms (for 1x1 matrices equal to 1 it gives 4 instead of 2).
- `simplify_gpgp(ps*ps)` gives `p(L_0)*p(-L_0)`; the reference multiplies a
  lone pair by 4 (the trace), although the slashed square is `p.p` times
  the identity, and does not do so when other factors remain
  (`ps*qs*qs`). `qs*ps*ps*qs` likewise gives `p.p*q.q` (reference
  `4*p.p*q.q`).
- `Sum.euler_maclaurin` and `approximants`: unchanged, see the existing
  notes (the reference `log(b) - log(a)` and the trailing `nan`).

Remaining unsupported (all raise; none returns a wrong value):

- `AlgebraicField.galois_group` of degree 5 and 6 (`C6` doctest) raises
  NotImplementedError: SymPy's `galoisgroups`, `galois_resolvents` and
  `resolvent_lookup` tables are not ported.
- `julia_code` / `octave_code` of a user function applied to a Matrix
  (`g(mat)`) raise TypeError: native expressions cannot hold an explicit
  matrix as a function argument.
- `DiscreteMarkovChain` with `Str` states raises NotImplementedError:
  `Str` is not an expression atom in symox, so `Eq(Str(...), Y[t])` cannot
  be built. Use symbols or integers as states.
- `MultivariateNormal` with `MatrixSymbol` parameters of symbolic dimension
  raises NotImplementedError: the reference density needs `exp` of an
  abstract matrix expression, which the native core does not implement.
- `shape(NewExpr(M))`: a user subclass of `Expr` cannot hold a Matrix
  argument (arguments are stored natively), so this raises TypeError.
- `convert_indexed_to_array`: summations that need the full reference
  canonicalization (single-axis contractions such as `Sum(A[i, j], (j, 0,
  d - 1))`, diagonals over several axes, `KroneckerDelta`, `Add`, `Pow`,
  functions) fall back to the matrix-summation path, which gives an equal
  form such as `Trace(A)*Trace(B)`, or raise NotImplementedError.
- Array expressions apply only part of the reference canonicalization:
  nested contractions, diagonals and permutations inside other operations
  are not flattened, and `convert_array_to_matrix` leaves a `PermuteDims`
  of a tensor product unconverted (an equal form, where the reference gives
  for example `ArrayTensorProduct(A.T, B)`).
- Printing: `Trace(A)*B` prints as `Expr(A)*B` (native printer), which
  affects some `convert_array_to_matrix` results containing a trace; the
  value is correct.

## Error types and messages of the category-f doc examples (f-final, 2026-10-23)

The 15 objects of the final triage whose SymPy 1.14 doc example raised a
different error type or message now raise SymPy's exception type with its
text. Where the message contains a class repr, the package name is `symox`
instead of `sympy`.

- `nseries` of a function whose argument tends to a logarithm at 0
  (`sin(log(x))`, `Abs(sin(log(x)))`, `cos(log(x))`) raises
  `PoleError: Cannot expand sin(log(x)) around 0`, naming the innermost
  function as `Function._eval_nseries` does. Other expansion points keep the
  native message (SymPy prints the substituted dummy there, `sin(1/_k)`).
- `Derivative` and `diff` with respect to a non-differentiable variable
  raise `ValueError("\nCan't calculate derivative wrt x*y.")`.
- `bool()` of an undecided relational raises
  `TypeError: cannot determine truth value of Relational: x < 1`
  (`1 < x < 2` reports `x > 1`, as in SymPy).
- Item assignment on an immutable dense matrix raises
  `TypeError: Cannot set values of <class 'symox.matrices.immutable.ImmutableDenseMatrix'>`.
  The Python class is now named `ImmutableDenseMatrix` (in module
  `symox.matrices.immutable`) with `ImmutableMatrix` as its alias, as in
  SymPy. Sparse matrices report `symox._sparse.ImmutableSparseMatrix`.
- `Predicate("P").register(...)` and `register_many` raise
  `TypeError: <class 'symox.assumptions.assume.PredicateMeta'> cannot be dispatched.`;
  the predicate metaclass is named `PredicateMeta` as in SymPy.
- `Rational` accepts only what SymPy accepts: integers (including bool),
  rationals, floats, `Float` and strings, in either slot. `Rational(pi)`,
  `Rational(x)`, `Rational(sqrt(2))`, `Rational(I)`, `Rational(oo)`,
  `Rational(nan)`, `Rational(pi, 2)` and `Rational(1, pi)` raise
  `TypeError: invalid input: pi` (previously `Rational(pi)` returned `pi`).
  `Rational(p, 0)` is `zoo` and `Rational(0, 0)` is `nan` (previously
  `ZeroDivisionError`), and `Rational(True)` is `1`.
- `banded` reports `"\nsequence does not fit an integral number of times in the matrix"`.
- `sympify` of a string wraps only parse failures (`SyntaxError`,
  `TokenError`) into `SympifyError`, like SymPy; evaluation errors of the
  parsed code propagate unchanged. Arithmetic on a Symox class
  (`sympify("O + 1")`, `-O`, `1 + O`) raises Python's `TypeError` text with
  the operand type `'type'`, as for SymPy's plain classes.
- `aesara_code` and `aesara_function` emit SymPy 1.14's
  `SymPyDeprecationWarning` before the `ImportError` raised when Aesara is
  absent.
- `nsolve` failures report the residual and the tolerance as mpmath's
  `findroot` does: `Could not find root within given tolerance. (r > tol)`
  followed by `Try another starting point or tweak arguments.`. The
  verification uses mpmath's criterion, the squared Euclidean norm of the
  residual vector against `tol`, and the default tolerance is mpmath's
  `eps * 2**10` at `prec + 20` bits, `2**(-prec - 9)` (`2.168...e-19` at 15
  digits; it was `2**(-prec - 8)`). A failed machine-precision Newton run
  for a system now falls back to the arbitrary-precision solver, which
  decides the result or the error. A Newton iteration that cannot decrease
  the residual stops and is judged by that final check, as in mpmath,
  instead of raising a separate message.
- `warns_deprecated_sympy` already reports SymPy's "DID NOT WARN" text; the
  class is pytest's `Failed`, which SymPy also uses when run under pytest.

Intentional divergences:

- The residual printed by a failed `nsolve` differs from SymPy's. For the
  doc example `nsolve(cos(x)*cosh(x) - 1, 3.14*100)` both solvers reach the
  same root `312.588469032184` (`verify=False` gives it in both), and the
  residual is the rounding noise of `cosh(312) ~ 1e135`. Symox iterates
  with 64 guard bits where mpmath uses 20, so the noise differs (Symox
  `5.0e+214`, SymPy `1.39e+230`). The tolerance printed is the same.

Remaining gaps:

- `nsolve(x**2 + 2, I)` gives `-5.2e-25 + 1.4142135623731*I`, while SymPy
  prints `1.4142135623731*I`. The real part is below the requested
  precision. mpmath's secant iterate happens to round it to exactly zero,
  and Symox's stopping rule ends one iteration earlier.
- `nsolve` of a singular two-variable system raises `ValueError: Jacobian is
  numerically singular`, where mpmath raises `ZeroDivisionError: matrix is
  numerically singular`.
- The machine-precision Newton fast path for systems still verifies each
  squared residual against `1e-14`, not against mpmath's tolerance.

## Algebraic fields, root objects and inexact domains (algfield, 2026-10-23)

Probes against SymPy 1.14 for TODO section 1 (algebraic domains in
polynomial operations, polynomial rings and fraction fields, exact reduction
of radicals and root objects, numeric/expression domain differences, EX root
forms and nroots digits). Two probe sets went from 158/177 to 166/177
and to 57/59 matches. Regression tests: `tests/test_algfield_parity.py`, with
values hard-coded from SymPy 1.14.

Implemented:

- `gcd`, `lcm`, `div`, `cancel` and `factor` with `extension=True` accept
  `CRootOf` coefficients over the same symbol as the generator, as SymPy
  does. Explicit generators keep SymPy's `PolynomialError`.
- `CRootOf` nodes evaluate natively in `evalf` from their isolating
  enclosure. Their sign, `is_zero` and `equals` are decided exactly, so
  `simplify`, `Max`, `Min`, `sign`, `Abs`, `floor` and `Piecewise` see the
  same decisions as SymPy.
- `simplify` returns 0 at once for a numeric expression known to be zero,
  like SymPy's `simplify`.
- The extended signs (`is_extended_positive` and the other three) of a finite
  numeric non-atom delegate to the plain sign queries.
- Integer roots follow `Integer._eval_power` prime by prime: perfect
  powers are extracted, a remainder sharing a factor with the denominator
  gives a reduced radical, and the other primes share one radical. For example,
  `sqrt(8)`, `12**(2/3)` and `(2**6*3**4)**(1/4)` give SymPy's forms.
- A Float coefficient is a factor of `Term`, both in `exprtools` and in
  the `gcd_terms` used by `together`. `factor_terms`, `gcd_terms`,
  `together` and `factor` over RR therefore pull out the float content as SymPy
  does (`factor(2.0*x**2 - 2.0)` is `2.0*(x - 1)*(x + 1)`).
- `apart` of an expression with Float coefficients is a port of SymPy's
  `apart` with `apart_undetermined_coeffs` over RR/CC, without turning
  floats into binary rationals.

Intentional divergences:

- `nroots` of a polynomial with an exact double root gives the double root
  exactly. SymPy's floating iteration gives two nearby values.
- `minimal_polynomial` of sums and products of `CRootOf` gives a result,
  where SymPy 1.14 does not finish in a reasonable time. Symox `is_zero`
  gives None on a case where SymPy also does not finish.

Remaining gaps:

- `apart` with a non-field coefficient domain and `factor(...,
  extension=a**(1/3))` raise different error types than SymPy.
- `nsimplify((-2)**(1/6))`: SymPy gives the expanded complex radical.
- An explicit `Poly(..., domain=CC)` raises NotImplementedError, from the
  inexact backend in `crates/symox-py/src/domains.rs`. An inferred CC domain
  and `CC[y]` work.
- `nsimplify(CRootOf(x**2 - 2, 1) - sqrt(2))` is kept. SymPy gives 0
  through evalf chopping.
- `gcd(Float, Float)` without generators raises.
- `10**(2/3)*2**(1/3)` is not combined into `2*5**(2/3)`, because the Mul
  combination of radicals with different bases is missing.
- `apart(1/(x**2 + 1.0*I))` falls back to the native path and gives
  `1/(x**2 + I)`.
- `roots` dictionaries can differ in key order only.

## Solve options, Abs/Min/Max solving and periodic inequalities (solve-rest, 2026-10-23)

Probe of the open items of TODO section 3 against SymPy 1.14 (303 cases:
249 identical before, 273 after; periodic inequalities 11 -> 17 of 19;
`periodicity` 41 -> 42 of 42; Abs/Min/Max/Piecewise solveset 21 -> 35 of 36).

Implemented:

- `periodicity` is a port of SymPy's (`_periodicity`, `lcim`, the trig
  `_period`, `decompogen`, Abs/Mod/exp cases, `check`). So `sin(x)**2` has
  period `2*pi` as in SymPy, and the periodic inequality windows match
  (`sin(x)**2 > 1/4` on `[0, 2*pi]`). SymPy's `f.simplify()` turns
  `sin(x)*cos(x)` into `sin(2*x)/2`. Symox's simplify does not, so
  `periodicity` adds a deep `trigsimp` step under SymPy's 1.7 `count_ops`
  ratio. The gap itself is in `simplify`.
- `minsolve_linear_system` is SymPy's port (the quick heuristic and the
  combination search). `solve(..., particular=True[, quick=True])` builds the
  linear system and uses it, with SymPy's ValueError messages.
- `solve` flags: `exclude` (with and without symbols), `manual` for systems,
  and `numerical`/`minimal` (only used by checksol). Unknown flags that SymPy
  never reads (`positive=True`) are ignored instead of raising.
- `solve` rewrites `Abs` of the symbol as the Piecewise of `separatevars`,
  as SymPy does: `solve(Abs(xr) - 2, xr) == [-2, 2]`. The argument that is
  neither real nor imaginary raises SymPy's NotImplementedError.
- `solveset` over a real domain reproduces SymPy's Abs/Min/Max/Piecewise
  steps when the native solver returns a ConditionSet:
  - the symbol is made real and the other heads are rewritten as Piecewise
    and folded;
  - each `as_expr_set_pairs` piece is solved on its own set;
  - `p*Abs(q) + r` with `p` and `r` free of the symbol is inverted
    (`_invert_abs`, with the `Contains(rhs, Interval(0, oo))` condition);
  - the other equations are split on the sign of a polynomial `q`
    (`_solve_abs`).
  In the complex domain SymPy's "Absolute values cannot be inverted"
  ValueError is raised.
- A wrong native result is fixed. A univariate equation whose Abs arguments
  and remaining part are polynomial now follows `_solve_abs` before the
  native solver, which lost the pieces where a branch vanishes
  identically. For example, `solveset(x + Abs(x - 1) - 1, x, Reals)` was
  `{1}` and is now `Interval(-oo, 1)`, and `Abs(x) + Abs(x - 1) - 1` gives
  `Interval(0, 1)`.
- `solve_univariate_inequality` of a continuous Min/Max of polynomials
  uses SymPy's interval test between the zeros given by solveset
  (`Max(x, 1) > 2`).
- `reduce_inequalities` of a single non-rational univariate inequality
  works on a real symbol, as SymPy's `_reduce_inequalities` does. The
  infinite bound disappears: `sqrt(x) > 1` gives `1 < x` and `exp(x) < 1`
  gives `x < 0`.
- The "cannot be solved using solve_univariate_inequality" message is
  passed through `filldedent`, as in SymPy.

Intentional divergences (SymPy is wrong or less precise):

- `solveset(cos(2*x) < 0, x, Interval(0, 2*pi))`: SymPy returns only
  `(pi/4, 3*pi/4)`. Symox also returns `(5*pi/4, 7*pi/4)`. The same holds
  for `sin(2*x) > 0` and for `cos(x) > 1/2` at the endpoint `2*pi`.
- `solveset(Abs(x**2 - 1) - y, x, Reals)` keeps
  `Intersection({-sqrt(1 - y), sqrt(1 - y)}, Reals)`, because those values
  are not real for `y > 1`. SymPy lists them unconditionally.
- `minsolve_linear_system` (and `particular=True`) also tries a single
  nonzero value. SymPy's search stops at two
  (`range(n0 - 1, 1, -1)`), although its docstring promises the minimal
  number. For example, `solve([x + 2*y + 3*z - 6, 2*x + y - 3],
  particular=True)` gives `{x: 0, y: 3, z: 0}` (SymPy
  `{x: 3/2, y: 0, z: 3/2}`).
- `solveset(Abs(sin(x)) < 1/2, x, Interval(0, 2*pi))` is solved. SymPy
  leaves a ConditionSet because its `function_range` fails on the
  Piecewise derivative.
- `periodicity` keeps the ValueError for a variable that is not a Symbol.
  SymPy returns 0 for `periodicity(x, x*x)`.
- `periodicity((sin(x) + cos(x))**2, x)` gives `pi` where SymPy gives
  `2*pi`. Both are periods, and the simplify forms differ.
- Extra valid Lambert branches (`x**3*exp(x) - 1`, `exp(x) - x**3`,
  `exp(x) - 2*x - 1`, `3**x - 2*x - 1`), the ImageSet of
  `solveset(2**x - 4, x)` (SymPy gives a ConditionSet), and inequalities
  that SymPy refuses (`exp(x) > x + 1`, `x*exp(-x) < 1/4`,
  `log(x) > x - 2`).

Remaining gaps:

- `solveset(x**Rational(1, 3) + 2, x, Reals)`: SymPy's `invert_real` uses
  the real cube root and gives `{-8}`. Symox uses the principal branch
  and gives EmptySet.
- Radical solveset `sqrt(x + 1) - sqrt(x) - 1` and `sqrt(x) + sqrt(x - 1) - 1`
  (SymPy `{0}` and `{1}`) stay ConditionSets.
- `solveset(Mod(x, 3) - 1, x, Integers)` (SymPy `ImageSet(3*n + 1)`) stays a
  ConditionSet.
- `roots(x**3 - 3*x + 1, trig=True)` does not give the cos form.
- Forms that still differ:
  - `solve` of `x*exp(x**2) - 1`, `x**y - 2`, `x**a*exp(x) - 1`, `x**n - 1`;
  - `solve([x**2 + y**2 - 1, x - a])` gives `sqrt(-(a - 1)*(a + 1))`;
  - the Union order of `solveset(sin(x) - sin(2*x), x, Reals)`;
  - `Interval(-oo, oo)` is printed where SymPy prints `Reals`.
- The `reduce_inequalities` NotImplementedError message prints the original
  symbol, where SymPy prints its real Dummy (`_a`). The strict-`roots`
  message is not wrapped with `filldedent`.

## ODE corpus forms and errors (ode-forms, 2026-10-23)

ODE corpus: cls 296 -> 298, hint 263 -> 268, def 259 -> 263 (out of 315).
Regression tests: `tests/test_ode_forms_parity.py`.

Implemented:

- `subs` and `collect` no longer enter a `ComplexRootOf` (SymPy's
  `ComplexRootOf._eval_subs` returns self). `nonhomogeneous`'s
  `collect(sol, 0)` for a real root replaced the root index by a dummy and
  raised (euler_hom_08, lin_const_coeff_hom_33..37).
- `cse` neither counts nor enters a RootOf (`tree_cse._find_repeated`), so
  constantsimp no longer rebuilds `CRootOf(x0, i)`.
- `cancel(..., extension=True)` infers the generators of `exp(2*x)` and
  `exp(x)` through the same power decomposition as the other conversions
  (rational_riccati_01 hinted).
- `terms_gcd(..., clear=False)` keeps a polynomial over `QQ_I` with a
  composite generator instead of raising CoercionFailed
  (var_of_parameters_13/14).
- `1st_rational_riccati`: `solve_aux_eq` and `remove_redundant_sols` fall
  back to `EX` when the coefficients mix algebraic fields. SymPy uses a
  polynomial ring over `QQ<a>`, which Symox does not have
  (rational_riccati_08 now gives SymPy's solution).
- `simplify`: on a count_ops tie between `factor_terms(selected)` and the
  selection, the selection is kept when SymPy's own sequence ends on it
  (`-2*x*f(x) + x*exp(-x**2)` stays a sum). The `factorable` classification
  of dsolve depends on that form through `solve`; lie_group_08/12 no longer
  list `factorable`.
- `re`/`im` of a numeric argument take the numerical sign refinement of the
  other functions. It is inert for `CRootOf` until a CRootOf evaluates
  inside the sign query (see below).

Intentional divergences:

- 2nd_2F1_hyper_01..03 (default hint `2nd_power_series_regular`): SymPy's
  series are wrong (residual of order 1: for
  `x*(x - 1)*f'' + (3/2 - 2*x)*f' + 2*f`, SymPy gives
  `1 + 4*x/3 + ...`, while `1 - 4*x/3` is an exact solution). Symox keeps
  the correct series.
- `solve(7*log(y - 3)/2 - 7*log(y + 3)/2 - C1, y)` returns `[]` in both
  (SymPy's false empty result is reproduced, not fixed).

Remaining gaps:

- lin_const_coeff_hom_33..37: Symox keeps `sin(x*Abs(im(r)))` for a
  `CRootOf` r. SymPy proves `im(r) < 0` numerically and moves the sign into
  the constant. The sign query does not evaluate a CRootOf yet.
- lie_group_08/12 default: `integrate(x*(2*y - exp(-x**2))*exp(x**2), x)`
  gives `log(exp(-x**2))/2` (manualintegrate's form) where SymPy's risch
  step comes first and gives `-x**2/2`.
- dep_div_indep_04, undet_26: SymPy's `simplify` gives
  `(-x/sin(f(x)/x) + f(x))/x` and `sqrt(2)*sin(x + pi/4)` for the solved
  derivative, and `factorable` matches on that form; Symox keeps the sum.
- separable_26: `solve(1.0*log(y - 3) - 1.0*log(y + 3) - C1, y,
  rational=False)` raises; SymPy reaches `-3/tanh(C1/2)` through `_tsolve`'s
  `factor` step and a generic substitution of the independent part.
- var_of_parameters_13 (hinted): both solutions pass checkodesol. SymPy's
  keeps terms such as `x*(3*exp(I*x)/2 + 3*exp(-I*x)/2)` inside the
  coefficients of `sin(x)` and `cos(x)`; Symox's integration of the
  Wronskian quotients gives the shorter, equivalent form.
- rational_riccati_01 default: SymPy's `separable` solves the implicit log
  solution to `-sqrt(2)/tanh(C1 - sqrt(2)*x)`; Symox keeps the implicit
  equation: `solve(x + sqrt(2)*log(y - sqrt(2))/4 - sqrt(2)*log(y +
  sqrt(2))/4 - C1, y)` raises, the same gap as separable_26 (a common
  non-rational coefficient of the logarithms).
- Other differing forms not analysed in detail: euler_undet_03/06,
  fact_12/13/16, 1st_exact_11, almost_lin_01/05, undet_38/39,
  liouville_08, separable_13/14/21/27, rational_riccati_04/07/09/10,
  lie_group_06/07/09, 1st_homogeneous_coeff_best_03/08, dep_div_indep_02,
  linear_coeff_01. Reference timeouts or crashes (nothing to do): fact_06,
  fact_11, 1st_exact_05/07, 2nd_nonlinear_autonomous_conserved_01..05,
  separable_reduced_01/02/04/08/11, dep_div_indep_03/05.

## Wave 5 integration leftovers (integ-rest)

Done:

- `Integer**Rational` now follows SymPy's `Integer._eval_power`: exact
  root first, then `perfect_power`, then `factors(limit=2**15)` with the
  `out_int`/`out_rad`/`sqr_gcd` regrouping, and `n**(a*(q-1)/q)/n**a` for a
  negative exponent. For example, `1922**(2/3)` is `31*124**(1/3)`,
  `72**(1/3)` is `2*3**(2/3)`, `248**(1/6)` is `sqrt(2)*31**(1/6)` and
  `12**(-1/3)` is `18**(1/3)/6`, as in SymPy.
- Products of positive-integer radicals follow `Mul.flatten`: they are
  grouped by base, then by exponent (`comb_e`), then the common integer gcd
  of two radicals is extracted (`num_rat`), and the result is regrouped by
  exponent (`pnew`). For example, `sqrt(12)*18**(1/3)` is
  `2*3**(5/6)*6**(1/3)`, `2**(1/3)*6**(1/4)` is `2**(7/12)*3**(1/4)` and
  `1922**(1/3)*62**(2/3)` is `62*31**(1/3)`. Held powers of a Rational
  base take part as in `pnum_rat`, with the Rational gcd
  `igcd(p)/ilcm(q)`. A `1/q` base reads as `q**(-e)` through
  `Pow.as_base_exp`, and a negative base joins `neg1e`, so
  `Mul(Pow(1/2, 1/3, evaluate=False), 6**(1/3))` is `3**(1/3)`. Before
  this wave, 14 of the 23 integer probes and 6 of the 13 Rational-base
  probes in `tests/test_integ_rest_parity.py` differed from SymPy. Now
  none of the integer probes and 3 of the Rational-base probes differ.
- A power whose integer part cannot be materialized stays held on the
  base of a perfect power: `4**((10**20 + 1)/4)` is `2**((10**20 + 1)/2)`.
  SymPy tries to build `2**(5*10**19)` and never finishes.
- Integration corpus: the 18 MISSING and 10 DIFF cases were triaged again.
  Every one is an intentional hold or a reference error that is already
  documented: the definite integrals with parameters above, where SymPy's
  value is wrong for some parameter values; `cos(s*x)` on `(0, oo)` with a
  nonzero `s`, where SymPy returns `sin(oo*sign(s))/s`; and the
  `laplace_transform`/Fourier values. No case was changed.

Remaining gaps:

- SymPy's `Mul.flatten` result depends on the argument order when a held
  negative Rational-base power meets an integer radical:
  `Mul(6**(3/4), (-1/3)**(3/2))` is `-sqrt(3)*6**(3/4)*I/9`, but
  `Mul((-1/3)**(3/2), 6**(3/4))` is `-24**(1/4)*I/3`. Symox flattens in
  canonical order and gives the second form in both cases, so
  `(-sqrt(6)/3)**(3/2)` is `-24**(1/4)*I/3`, where SymPy's
  `_eval_expand_power_base` gives the first. The two forms are equal. The
  case is a strict xfail in `tests/test_pow_negative_base_parity.py`.
- `ratint(1/(x**3 + x + 1))` raises `PolynomialDivisionFailed` in
  `log_to_atan`, while SymPy returns a 16653-character closed form. The
  `1922**(2/3)` reduction is now in place, but it is not enough. The
  radicand that `roots_quadratic` passes through `simplify` gets a
  different form from symox's native `simplify` than from SymPy's
  pipeline. The resulting `EX` constant then mixes `1922**(1/3)` with
  `2**(1/3)*31**(2/3)`, so `cancel` cannot prove `(1/b)*b - 1 == 0`, and the
  division by a constant polynomial fails. A step-by-step emulation of
  SymPy's simplify pipeline with symox functions gives SymPy's
  intermediates. Fixing this needs SymPy's `simplify` pipeline (or its
  output form for this radicand) in the native simplify.

## Simplification and matrix probes (simp-mat, 2026-10-23)

Probe sets were generated from the assertions of SymPy 1.14's
`sympy/simplify/tests/*.py` and `sympy/matrices/tests/*.py`, comparing each
asserted expression with SymPy statement by statement.

- simplify tests: same 1544 -> 1573, different 194 -> 169, exceptions
  31 -> 27 (`test_simplify.py`: same 288 -> 317).
- matrix tests: same 1649 -> 1728, exceptions 145 -> 73.

Implemented:

- `nsimplify` follows SymPy: `tolerance=0` uses PSLQ through `identify`,
  `rational=True` uses `_real_to_rational` (with `nsimplify(fl,
  rational=False)` for each Float), infinite Floats map to `oo`/`-oo`, and
  Matrices are mapped entrywise. The decimal-fraction fallback keeps its
  explicit resource guard for exponents larger than 10**6.
- `radsimp` replaces a numeric radical denominator that `nsimplify` proves
  equal to a Number (SymPy's `d.is_number` branch), and `Expr.equals` ends
  with SymPy's `minimal_polynomial` test for numeric differences.
- `simplify` now runs these steps of SymPy's pipeline around the native
  search: `sign._eval_simplify`, `collect_abs`, the deep replacement that
  simplifies function arguments, the early return for function heads that
  the reference does not touch, then `sign` rewritten as `Abs`,
  `hyperexpand`, `kroneckersimp` and `besselsimp`. A `Piecewise` produced
  by that last `hyperexpand` is not simplified piecewise again, because the
  reference runs its Piecewise block before `hyperexpand`.
- Matrices: `is_positive_definite`, `is_positive_semidefinite`,
  `is_negative_definite`, `is_negative_semidefinite`, `is_indefinite`,
  `is_echelon` and `is_weakly/strongly_diagonally_dominant` are ports of
  SymPy's `eigen.py`/`reductions.py`, including the Cholesky and
  division-free elimination tests and `_find_reasonable_pivot`.
  `eigenvals`/`eigenvects` of numeric Float matrices use SymPy's mpmath
  path with its precision doubling and `PrecisionExhausted`, and
  `eigenvals` accepts `simplify` (bool or function) and `rational`.

Intentional divergence: `nsimplify` raises `ValueError` for a negative,
infinite or NaN tolerance, where SymPy accepts `-1` and fails with an
`AssertionError` for `inf`/`nan`.

Remaining gaps (probe evidence):

- simplify: `nthroot`; `nc_simplify` of noncommutative `MatMul` scalars;
  forms of `hypersimp`, issue 17141 and issue 22210 (derivative order);
  issue 11004 Float precision; `posify` with `Integral` and `Matrix`;
  `cos` has no `_eval_rewrite_as_exp` for `simplify_complex`;
  `log(exp(x), 2)` with `inverse=True`; some `KroneckerDelta` products.
- fu: `TR4` on `False`, missing `_TR11`/`_TR56`, `hyper_as_trig`, issue
  25590. hyperexpand: `lerchphi`/`meijerg` timeouts, `test_plan`,
  `meijerg` with Floats. powsimp (10) and radsimp (5) form differences;
  cse differences (91) were not triaged.
- matrices: eigenvalues of characteristic polynomials with irreducible
  factors of degree > 4 (SymPy returns `CRootOf`), Jordan forms of Float
  matrices, `cofactor`/`minors` index errors, Matrix hashing in
  `columnspace` comparisons, decomposition and solver form differences in
  `test_matrixbase.py` (54 different, 52 exceptions).

## Distribution integrals and Integral special integrands (integ-risch, 2026-10-23)

A probe set of 732 integrals was taken from SymPy 1.14's `test_risch.py`,
`test_heurisch.py`, `test_manual.py`, `test_integrals.py`,
`test_deltafunctions.py` and `test_singularityfunctions.py`. Symox gave no
wrong result on it before or after this wave. The Risch and heurisch
leftovers listed in the TODO already matched SymPy: `atan(exp(x))`, the
`RootSum` forms, partial fractions with `risch=True`, `atan(x)*log(x)`,
`x**(4/3)*log(x)` with `heurisch=True`, and `1/(cos(x) + 2)`.

Implemented:

- Delta integrals now use SymPy's `deltaintegrate` stage in more cases:
  - an order-0 delta of a rational polynomial of degree 2 or more, when
    the guarded expansion splits it into simple real roots. Examples:
    `DiracDelta(x**2 - 1)`, `cos(x)*DiracDelta(x**2 - 1)`,
    `DiracDelta(2*x**2 - 8)/x`, and `DiracDelta(x**2 - 3*x + 2)*x` on
    `(0, 3)`;
  - a unit-slope delta weighted by `Heaviside`, `sign` or `Abs`, as in
    `Heaviside(x - a)*DiracDelta(x - 1)`, which gives
    `Heaviside(1 - a)*Heaviside(x - 1)`;
  - a delta weighted by an undefined function: `f(x)*DiracDelta(x)` gives
    `f(0)*Heaviside(x)`, and `f(0)` over the real line.

  The stage is skipped for Float input, for integer powers above 1000,
  for deltas of order 1 or more (SymPy keeps only one term of
  `g(x)*DiracDelta(x - p, n)`), and when the result still holds an
  `Integral`. In that last case, `DiracDelta(x**2 - 1)/(x - 1)` stays held
  where SymPy returns `zoo`.
- `integrate` handles the special integrands that SymPy's `Integral`
  handles:
  - `Eq(l, r)` gives `Eq(integrate(l), integrate(r))`;
  - an integrand that is a `Sum` with finite limits is integrated term by
    term. A `Sum` with an infinite limit stays held. When a limit of the
    sum is an integration variable, `integrate` raises SymPy's
    `ValueError`;
  - a limit `(x, None, b)` or `(x, a, None)` keeps only the given
    endpoint, with strings sympified;
  - `O(g)` integrates to `O(G)`, where `G` is an antiderivative of `g`.

Probe before and after this wave:

| Category | Before | After |
|---|---|---|
| same | 620 | 634 |
| equivalent | 36 | 34 |
| differs (form, or a reference error or timeout) | 46 | 45 |
| symox unevaluated | 19 | 14 |
| symox exception only | 3 | 0 |
| SymPy exception only | 3 | 0 |
| symox wrong | 0 | 0 |

Intentional divergences, unchanged:

- For derivatives of deltas, symox applies the `1/(a**n*|a|)` scale of an
  affine argument. SymPy omits it: for `exp(x)*DiracDelta(2*x - 3, 1)` on
  the real line, SymPy gives `-exp(3/2)`, and the correct value is
  `-exp(3/2)/4`.
- The real root of `DiracDelta(x**5 - x - 1)` is kept as a `CRootOf`
  (SymPy gives 0). It was checked numerically against `exp(r)/|g'(r)|`.

Remaining gaps:

- Indefinite integrals of `g(x)*DiracDelta(x - p, n)` with `n >= 1` keep
  SymPy's single leading term. For example, `x*DiracDelta(x - 1, 1)` gives
  `-Heaviside(x - 1)`, without the `DiracDelta(x - 1)` term. The definite
  values are exact.
- `SingularityFunction` integrals over `(-oo, oo)` stay held, because
  `SingularityFunction` has no limit at infinity. SymPy gives `1` and
  `f(5)`.
- `Heaviside(x - 1)*DiracDelta(x - 1)` gives `Heaviside(x - 1)**2/2` from
  the native step, where SymPy gives `Heaviside(x - 1)/2`. The two are
  equal as distributions.
- `meijerint_definite` does not find these SymPy values:
  - `sqrt(log(1/x))` on `(0, 1)`;
  - the Frullani integrals of `test_integrals.py`;
  - `log(log(x))/x**n` on `(1, oo)`.
- Integrals over the unit cube with a parameter in an inner integrand,
  for example `1/(x + y)**2` on `(0, 1)`, stay held because of a possible
  pole.
- `integrate(1/(1 + exp(x - 1/2)/(1 + exp(x))), x)` takes more than 30 s.
- `DiracDelta(x)**2` over the real line stays held (SymPy gives
  `DiracDelta(0)`).
- `ratint(1/(x**3 + x + 1))` is still the gap described in the integ-rest
  section above. Two further findings:
  - SymPy's `roots_quadratic` calls `simplify` for an `EX` domain, but the
    native `roots_quadratic` in `poly/roots.rs` calls `factor` for every
    non-numeric domain.
  - Patching the Python `roots_quadratic` with an emulation of SymPy's
    simplify pipeline still fails in `log_to_atan`. The remaining
    difference is therefore in the normalization or `cancel` of the `EX`
    domain, not only in `simplify`.

## Solve and ODE leftovers (ode-solve, 2026-10-23)

Regression tests: `tests/test_ode_solve_parity.py`. This wave settles the
"Remaining gaps" of the solve-rest section above.

Implemented:

- `solveset` of a modular equation over the integers: ports of
  `_is_modular`, `_invert_modular` and `_solve_modular`
  (`python/symox/solvers/_solveset_steps.py`). `Mod(x, 3) - 1` gives
  `ImageSet(Lambda(n, 3*n + 1), Integers)`, `Mod(7**x, 41) - 15` gives
  `ImageSet(Lambda(n, 40*n + 3), Naturals0)`, `7 - Mod(x, 5)` and
  `Mod(x, 3) - I` are proven EmptySets, a domain that is not a subset of
  the integers keeps the ConditionSet. `pow(Integer, Integer, Integer)` and
  `pow(Expr, n, m)` (`Mod(x**2, 3)`) follow SymPy's three-argument pow.
- `solveset` radical step: when the native solver leaves a ConditionSet,
  SymPy's `_solve_radical` (unrad, then `checksol` on every candidate; an
  undecided candidate stays in a ConditionSet) runs if the
  symbol-dependent part is a sum, and the `_invert_real` branch for
  `g**(odd/even) = rhs` (rhs must be nonnegative) if it is a single power.
  `sqrt(x + 1) - sqrt(x) - 1` gives `{0}`, `sqrt(x) + sqrt(x - 1) - 1`
  gives `{1}`, and `sqrt(17*x - sqrt(x**2 - 5)) - 7` gives `{3}`.
- `Set.as_dummy` and `Set.dummy_eq`.
- `solve` fallback after the native NotImplementedError: SymPy's `_invert`
  and the `_tsolve` steps `factor(powdenest(...))`,
  `logcombine(lhs, force=True)`, `f(x)**e = rhs` with a symbol-free
  non-integer `e`, and the final posify retry
  (`python/symox/solvers/_tsolve_steps.py`). `x**y - 2` gives
  `[2**(1/y)]`, `x**n - 1` gives `[1]`. The separable/Riccati implicit
  solutions of rational_riccati_01, separable_26 and almost_lin_01 are now
  solved for `y`.
- `roots(..., trig=True)`: the trig branch of `roots_cubic` (cos/acos
  form when the discriminant is decidably positive); SymPy's
  `-x**3 + 2*x**2 + 3*x - 2` and `x**3 - 3*x + 1` forms are reproduced.

Intentional divergences:

- `solveset(x**Rational(1, 3) + 2, x, Reals)` stays EmptySet (SymPy:
  `{-8}`). SymPy's `invert_real` uses the real cube root there, but Symox
  and SymPy both evaluate `(-8)**(1/3)` to the principal value
  `1 + sqrt(3)*I`, so `-8` does not satisfy the equation, and SymPy's own
  Complexes result is EmptySet. The odd-denominator branch of
  `_invert_real` is therefore not ported.
- Complex radical candidates that SymPy leaves undecided in a ConditionSet
  can be decided by Symox's `checksol`: `x**(1/3) + sqrt(x) - 2` gives `{1}`
  where SymPy keeps the candidates `(-1 +- I)**6` (they are not solutions:
  `(-8*I)**(1/3) + sqrt(-8*I) - 2` is not zero).
- `exptrigsimp` no longer writes `(c*(1 + s*exp(x)))**m` as
  `c**m*(1 + s*exp(x))**m` for a non-integer `m` and a `c` that is not
  positive. SymPy does, and that is wrong:
  `simplify((-1 - sqrt(3)*I)*(exp(a) - 1)**(1/3)/2)` gives
  `(-1)**(1/3)*(-1 - sqrt(3)*I)*(1 - exp(a))**(1/3)/2`, which differs for
  `a = 2`. The solve fallback above exposed it: one cube root of
  almost_lin_01 was lost.

Remaining gaps:

- The tsolve steps give different but equal forms: rational_riccati_01
  `sqrt(2)*(-exp(2*sqrt(2)*C1) - exp(2*sqrt(2)*x))/(exp(2*sqrt(2)*C1) -
  exp(2*sqrt(2)*x))`, where SymPy writes `-sqrt(2)/tanh(sqrt(2)*(C1 - x))`
  (SymPy reaches it through the bivariate `_u` substitution and its
  `simplify`). almost_lin_01 gives `(exp(3*C1 + 3/x) - 1)**(1/3)`, where
  SymPy writes `(exp(3*C1*x + 3)**(1/x) - 1)**(1/3)`. The two forms agree
  for real `C1` and `x`, and Symox's form is the exact inverse of the
  equation.
- Lambert forms that need the log of a product of a positive unknown
  (`x*exp(x**2) - 1`, `log(x) + x**2`, `x**a*exp(x) - 1`) still raise
  NotImplementedError. SymPy reaches `_solve_lambert` with a posified
  symbol and `expand_log`.
- Radical solveset: `sqrt(x)*sqrt(x - 7) - 12` (a product, SymPy `{16}`),
  `sqrt(x) + sqrt(x + 1) + sqrt(1 - x) - 6*sqrt(5)/5` (the native Reals
  filter cannot decide that the Cardano roots with cube-root-of-unity
  factors are not real, and `roots` does not apply SymPy's
  `_integer_basis` scaling), and the `Complement(..., {-1/a})` form of
  `solveset(y - b/(1 + a*x), x)`.
- Symox has no public `Function.inverse()`. The solve fallback uses a
  table of SymPy's one-argument inverses.

## Precision near cuts and huge exponents (prec, wave 6, 2026-10-23)

Probe: 6225 evaluations (2075 expressions at 15, 50 and 200 digits) of every
multivalued elementary and special function symox evaluates (log, roots and
complex powers, inverse trigonometric and hyperbolic functions, LambertW
branches 0/-1/1/2, polylog, zeta/Hurwitz zeta/polygamma/dirichlet_eta,
Ei/li/Li/E1/expint/Si/Ci/Shi/Chi, gamma/loggamma/incomplete gamma, erf family,
hyper/meijerg, elliptic_k/e/f/pi, Bessel/Hankel/Airy) on and next to their
cuts (exact offsets such as `-2 + I/10**30`), near poles and zeros, and exact
cancellations. Each value was compared with SymPy's evalf and with mpmath on
the exact argument at 3n+60 digits.

- Relative accuracy of the complex value (n-3 digits required): 5947 correct
  in both, 188 correct in symox only (SymPy rounds the argument first), 29
  evaluated by symox only and correct (SymPy raises a pole error or returns
  `oo` for `atanh(1 + 10**-25)`, `zeta(1 - 10**-25)`, `gamma(-3 + 10**-25)`,
  ...), 30 left unevaluated by both, 22 inaccurate in both (18 identical
  strings). No case where SymPy is accurate and symox is not, no wrong branch,
  and no exception or timeout specific to symox.
- Componentwise, 7 cases give a tiny component with more printed digits than
  it holds. SymPy tracks the accuracy of each component and prints fewer
  digits (`0.e-21`, or a 26-digit mantissa). Examples are
  `(1 + 10**-25)**I` at 50 digits (32 correct digits in the imaginary part,
  SymPy prints 26), `(-1 + 10**-25)**(1+I)` at 15 digits (symox
  `-1.08e-30*I`, SymPy `0.e-21*I`, true `+4.3e-27*I`),
  `Li(2 + I/10**30)` at 15 digits (real part 5.2e-61 with 9 correct
  digits), and the cancelling sum
  `-I*exp(3*I/5)/2 - exp(3*I/10)*sin(3/10) + I/2 + 10**-25` at 50 digits.
  The complex value is accurate to n digits, but the small part's digits are
  not certified. Open gap: per-component accuracy tracking (SymPy's
  `complex_accuracy`) in the Pow, Li and Add paths of `evalf_mp`.

Huge exponents, fixed:
- Floats with 53-bit precision are stored as doubles, and their arithmetic
  overflowed or underflowed silently. `Float(2)**1024`, `Float(10)**400`,
  `Float(1e300)*Float(1e300)`, `Float(1)*10**400`, `Float(2)**1100.0` and
  `Float(2)**Rational(2049, 2)` gave `inf` or `oo`, and `Float(2.5)**-1100`
  gave `0.0`. These were wrong results. Add, mul, div and integer or real
  powers now detect a result or operand outside the normal double range and
  recompute it at 53 bits with a wide exponent. The results are SymPy's
  (`1.79769313486232e+308`, ...).
- `Float(b)**n` with an integer `n` beyond i64 (`Float(2)**(2**63)`,
  `Float(0.5)**(10**30)`, `Float(-2)**(2**63 + 1)`) is now evaluated like
  SymPy. Before, it stayed an unevaluated power.
- `Pow(-2, 10**20 + 1).evalf()` gave a spurious imaginary part
  (`... + 4.2e+30102999566398119492*I`). A negative real base to an integer
  power beyond i64 is now real, with the sign set by the parity of the
  exponent. SymPy does not finish these in 20 s.
- `gamma(x).evalf()` for large real `x` lost about log10(x log x) digits
  (`gamma(10**20/3)` was correct to 4 digits only). The argument is now
  evaluated with mag(x) + 16 extra bits, and Spouge's working precision
  includes mag(x). Integer arguments above 4096 use Spouge's sum instead of
  the quadratic exact product, so `factorial(10**6).evalf()` takes 3 ms
  instead of 169 s (SymPy: 3.8 s).
  Intentional correction: SymPy 1.14 rounds the argument to the target
  precision first and returns `3.11482135286148e+636286142112569524544` for
  `gamma(Rational(10**20, 3)).evalf()`. That is 52 orders of magnitude off.
  mpmath at 120 digits gives `2.70659425601764e+636286142112569524492`, which
  is symox's value. SymPy's `gamma(Rational(10**12, 3)).evalf()` is also
  wrong from its 6th digit.

Huge exponents, remaining (explicit, never a wrong number):
- Binary exponents stay bounded by `MAX_EXPONENT = 2**100`. Beyond it, symox
  keeps the expression unevaluated where mpmath (unbounded exponents)
  answers: `exp(10**30).evalf()`, `exp(2**100).evalf()`,
  `Float(2)**(2**101)`, `exp(Float('1e30'))`, `factorial(10**30).evalf()`.
  Lifting the bound needs big-integer exponents in `BigFloat`.
  `gamma(Rational(10**40, 3)).evalf()`, `gamma(Float(1e40))` and
  `factorial(10**40).evalf()` used to raise `OverflowError`, because the
  mpmath result has an exponent past i128. They now stay unevaluated too.
  So does `loggamma(10**40).evalf()` (`log(gamma(10**40))`, about 9.1e41).
  SymPy does not finish it in 60 s. Evaluating it needs a log-of-gamma path
  that avoids the huge intermediate.
- `Float(1)/Rational(1, 10**400)` gave `inf` (the divisor rounds to the
  double 0.0). It now gives SymPy's `1.00000000000000e+400`.
- Where SymPy hangs (`factorial(10**20).evalf()`, `Pow(2, 10**20).evalf()`,
  `loggamma(10**30).evalf()`, `Pow(Rational(3, 2), 2**80).evalf()`), symox
  answers. The values were checked against mpmath.

Tests: `tests/test_prec_parity.py`.

## Matrix leftovers (mat-rest)

Implemented (ports of SymPy 1.14):
- `eigenvals` and `eigenvects` of rational matrices whose characteristic
  polynomial has irreducible factors of degree > 4. `eigenvals` follows
  `_eigenvals_dict`/`_eigenvals_list` (`roots`, then `all_roots`, so
  `CRootOf`). `eigenvects` follows `_eigenvects_DOM`: nullspace over
  `QQ[lambda]/(minpoly)`, then substitution of each root. `jordan_form` of
  such a matrix raises SymPy's `MatrixError` ("... eigenvalues in CRootOf
  form").
- `jordan_form` of Float matrices: rationalize, compute, then restore floats
  at the input precision (`evalf(n=max_dps, chop)`). With `calc_transform`,
  the result is a list, as in SymPy.
- `inv(method='LU')`, `inverse_LU` and `inv(method='BLOCK')` run
  `_verify_invertible` (berkowitz det `.equals(0)`, then `rref`) on symbolic
  matrices. A singular symbolic matrix now raises `NonInvertibleMatrixError`
  instead of returning a matrix of `zoo`/`nan`. `BLOCK` is SymPy's
  `_inv_block`.
- `inv()` with noncommutative entries raises SymPy's `PolynomialError`
  (the default `DomainMatrix` path); `inv(method='LU')` works as in SymPy.
- `adjoint(f(A))` for a noncommutative `A` and a function `f`:
  `transpose(f(conjugate(A)))` (`Expr._eval_adjoint`).
- `minor`/`cofactor`/`minor_submatrix` index and squareness errors; boolean
  masks in `__getitem__` and `extract`; list keys (with repeats and masks) on
  `SparseMatrix`; `Matrix(2, 2)` `TypeError`; non-integer and negative
  dimensions of `eye`/`ones`/`zeros`/`Matrix`; `Matrix([[1, 2], [3]])`
  ("mismatched dimensions") and `Matrix([1, 2, [3]])` ("expecting list of
  lists") errors; `diag(..., strict=False)` of ragged lists.
- `is_Identity`, `is_upper_hessenberg`/`is_lower_hessenberg`, the Dirac
  conjugate `M.D`, and the `is_hermitian` fallback (`simplify(a - adjoint(b))`
  per entry pair when the native test is undecided).
- `solve(Matrix_equation, Matrix([a, b]))` (matrix of symbols).
- `sympy.matrices.utilities` (`dotprodsimp` context manager,
  `_get_intermediate_simp`) and `sympy.matrices.determinant`
  (`_find_reasonable_pivot`, `_find_reasonable_pivot_naive`).

Intentional divergences:
- `CRootOf` equality is structural: it depends on the generator name.
  SymPy's equality ignores it (PurePoly), and the generator SymPy prints
  (`x` or `_lambda`) depends on its cache.
- `inv(method='LU', iszerofunc=f)`: SymPy's `LUsolve` ignores `f`; symox
  applies it to the pivots, so a predicate that calls a pivot zero raises
  `NonInvertibleMatrixError` where SymPy divides by it.

Remaining gaps (explicit errors or different forms, never a wrong number):
- `inv(method='GE')` with noncommutative entries gives a different form
  (SymPy's result is also wrong for nc entries).
- `adjoint`/`transpose` of `Abs(A)` and `sqrt(A)` for noncommutative `A`.
- The deprecated mixin classes of `test_commonmatrix` (`_MinimalMatrix`,
  `ShapingOnlyMatrix` and the other `*OnlyMatrix`) do not exist.
- Class identity: `ImmutableMatrix` is `Matrix` with a flag, so
  `type(m) is ImmutableMatrix` and `Hashable` checks differ.
- `pinv` of some symbolic matrices times out.
- Non-integer matrix powers (`_matrix_pow_by_jordan_blocks`),
  `_eval_matrix_exp_jblock`/`_eval_matrix_log_jblock` helpers, `diff` by a
  matrix, `Mod` of a matrix, `match` on a matrix, `epsilon_eq`, `Point2D`
  entries, `DeferredVector` iteration, `Dagger`/`Operator` entries.
- Unsimplified forms differ in `cholesky`, `GramSchmidt`, `LDL`,
  `matrix_norm` (`Max`), `exp`, `refine` and `simplify(ratio=...)`, and the
  `expand` hint `numer`.
- `roots(multiple=True)` orders cubic roots differently (native `roots`).

Tests: `tests/test_mat_rest_parity.py`.

## simp-rest (wave 7): simplification leftovers

Implemented (statement-by-statement probe over `sympy/simplify/tests`):
- `cse`: SymPy's algorithm is now ported (`opt_cse` with `match_common_args`,
  then `tree_cse` rebuilt in `ordered` order, `basic_optimizations` with
  `sub_pre`/`sub_post` from the new `symox.simplify.cse_opts`, `ignore`,
  `postprocess`, `list=False`, sparse and immutable matrices, MatrixSymbol and
  MatrixExpr inputs). The former native shortcut produced the 91 differences
  noted by simp-mat. `symox.simplify._cse_diff` (`_forward_jacobian` and
  helpers) is ported too. `iterable(MatrixExpr)` is now False as in SymPy
  (`MatrixExpr._iterable = False`); it made `opt_cse` recurse forever.
- `Function.inverse(argindex=1)` on the builtin heads that define it in
  SymPy (`exp -> log`, `asin -> sin`, `erf -> erfinv`, `tanh -> atanh`, ...;
  `sin`, `cos` still have none). `inversecombine` is SymPy's (any head with an
  `inverse`, user `Function` subclasses included), so
  `simplify(f(g(x)), inverse=True)` gives `x` like SymPy.
- `trigsimp(..., inverse=True)` also cancels `asec(sec(x))` and
  `acsc(csc(x))` (SymPy accepts every inverse whose `g()(1)` is a
  TrigonometricFunction).
- `sech` and `csch` are hyperbolic functions for `hyper_as_trig`, `_osborne`
  and `_osbornei` (`sech -> sec`, `csch -> csc/I` and back), so trigsimp now
  reduces `tanh(x)**2 + sech(x)**2`, `cosh(x)*sech(x)`,
  `sinh(2*x)*csch(x)`, `csch(x)**2 - coth(x)**2` like SymPy (before, they were
  returned unchanged).
- `TR4` evaluates trigonometric functions of special angles that were built
  held (`evaluate(False)`); `_TR11`, `_TR56`, `_osborne`, `_osbornei` are
  importable from `symox.simplify.fu`.
- `f._eval_rewrite_as_exp(*args)` exists on the trigonometric and hyperbolic
  heads (used directly by SymPy's tests).
- `kroneckersimp` is SymPy's (Piecewise rewritten as KroneckerDelta, then
  cancellation to a fixed point), and `simplify` runs it before the native
  search when KroneckerDelta and Piecewise appear together:
  `simplify(KroneckerDelta(1, i) + Piecewise((1, Eq(j, 2)), (0, True)))` is
  `KroneckerDelta(1, i) + KroneckerDelta(2, j)` and
  `simplify(KroneckerDelta(0, j).rewrite(Piecewise)*KroneckerDelta(1, j))`
  is 0.

Remaining (documented, no wrong result found):
- `log(exp(x), 2, evaluate=False)` builds an undefined `Function('log')` with
  two arguments instead of a held `log`, so `simplify(..., inverse=True)`
  leaves it alone (SymPy: `x/log(2)`).
- `Derivative(sin(x), (x, d))` with a symbolic integer order `d` raises
  ValueError at construction (issue 22210); SymPy keeps it.
- Explicit-matrix expressions in `cse` print differently:
  `Inverse(ImmutableMatrix)` prints as `ImmutableMatrix(2, 2, [...])**(-1)`
  where SymPy prints `Matrix([...])**(-1)`, `MatrixSolve` is not a Basic, and
  the kalman_filter / matrix-matrix cases give equivalent but different
  replacements.
- `nc_simplify` of noncommutative matrix products (41 probe statements:
  MatrixExpr `doit` options unsupported), the noncommutative `simplify` of
  issue 25590 (`cos(x)*(2*sin(x)*B*A + cos(x))`, equivalent), hyperexpand
  timeouts (`test_prudnikov_1`, `test_lerchphi`, meijerg shift operators,
  `test_plan`, meijerg with Floats), trigsimp_groebner, powsimp_polar and the
  other powsimp/radsimp forms, nthroot, hypersimp form, posify with
  Integral/Matrix, the issue 17141 exponent form and the issue 11004 Float
  printing (`0.918937683105469` vs `0.9189`).
- `hyper_as_trig` results differ from SymPy only by the Dummy index.

Tests: `tests/test_simp_rest_parity.py`.

## Exponents past 2**100 and per-part accuracy (bigexp, wave 7, 2026-10-23)

Huge exponents, fixed:
- The binary exponent bound of `BigFloat` is now `MAX_EXPONENT = 2**120`,
  up from `2**100`. It stays an i128 field, so the common path is unchanged:
  the evalf and precision test files take 16.0 s against 17.2 s on the base
  commit. The bound leaves room below 2**127 for exponent sums and decimal
  scaling. SymPy's values are now returned for `exp(10**30).evalf()`
  (any precision), `exp(2**100).evalf()`, `exp(+-2**118).evalf()`,
  `Float(2)**(2**101)`, `Float(2)**(2**119)`, `exp(Float('1e30'))`, and for
  arithmetic on these values (`V**2`, `V**1000`, `1/V`, `V*V*V`, `sqrt(V)`,
  `V**0.5`, `log(V)`, comparisons, `_mpf_`).
- Numeric powers whose base has a wide exponent (`V**2`, `V**Rational(1, 2)`)
  are evaluated instead of being refused by the wide-exponent guard.
- `factorial(10**30).evalf()`: integer arguments beyond i64 go through
  Spouge's sum. mpmath at 120 digits agrees. SymPy does not finish in 20 s.
- The working precision of real and complex `x**y` now includes the size of
  ln(x). The absolute error of ln(x), multiplied by y, becomes the relative
  error of the result. Before this fix, the guard bits were too few for a
  base whose binary exponent is above about 2**20. Intentional corrections,
  checked against mpmath at 60 digits:
  - `(Float(2)**(2**40))**Rational(1, 3)` gives `4.31912023499076e+110328660180`.
    SymPy gives `4.31905934482588e+...`.
  - `V**Rational(1, 3)` with `V = exp(10**30).evalf()` gives
    `3.42789437569241e+144764827301083942550376306305`. SymPy rounds 1/3 to
    53 bits and moves the exponent by 8e15.
  - `Float('1e1000000', 30)**Rational(1, 3)` ends in `...356652`. SymPy's
    `...348497` is wrong.
- Intentional correction: `gamma(Rational(10**30, 3)).evalf()` gives
  `5.38088111291356e+9696194754459028578351281059261` (mpmath). SymPy rounds
  the argument first and gives `9.23792600424838e+...7643`.

Per-part accuracy, fixed:
- symox does not print fewer digits the way SymPy's `complex_accuracy` does.
  Instead, `eval_mp_root` recomputes a complex value whose parts differ in
  magnitude by `gap` bits at `prec + gap + 8` bits, so every printed digit of
  the small part is correct. The loop is bounded: past 512 extra bits, two
  evaluations are compared. Only the bits they agree on are kept, and SymPy's
  residue mechanism prints the shorter mantissa. A part with no agreeing bit
  becomes a zero of precision 1 at the error scale.
- Same strings as SymPy:
  - `((1 + 10**-25)**I).evalf(50)`
  - `Li(2 + I/10**30).evalf()`
  - `log(-1 - I/10**30).evalf(20)`
  - `asin(2 + I/10**30).evalf()`
- Intentional corrections, checked against mpmath. In these cases SymPy drops
  the small part to `0.e-n`, keeps a spurious one, or returns `oo`:
  - `((-1 + 10**-25)**(1 + I)).evalf()` gives `-0.0432139182637723 + 4.32139182637722e-27*I`.
    SymPy gives `+ 0.e-21*I`.
  - The cancelling sum `-I*exp(3*I/5)/2 - exp(3*I/10)*sin(3/10) + I/2 + 10**-25`
    at 50 digits gives `1.0000...e-25`, with an exactly zero imaginary part.
    SymPy gives `1.0e-25 + 0.e-80*I`.
  - `exp(I/10**200).evalf()` gives `1.0 + 1.0e-200*I`. SymPy gives `0.e-200*I`.
  - `sqrt(-1 + I/10**40).evalf()` gives `5.0e-41 + 1.0*I`. SymPy gives `0.e-41`.
  - `((1 + I/10**6)**(1/3)).evalf()` gives `3.33333333333272e-7*I`. SymPy
    prints only 10 digits.
  - `exp(I/10**8).evalf(30)` prints 30 correct digits of the imaginary part,
    where SymPy prints 23.
  - `atan(I*(1 + 10**-20)).evalf()` gives `1.5707963267949 + 23.3724245202204*I`.
    SymPy gives `oo`.

Remaining (explicit, never a wrong number):
- (Resolved by bigexp2, see "Arbitrary-size binary exponents" below.)
  Binary exponents beyond 2**120 still need big-integer exponents in
  `BigFloat`. Changing the field touches about 95 constructors and over 100
  field uses across the crates. The following stay unevaluated where SymPy
  answers:
  - `exp(10**40).evalf()`, `exp(10**100).evalf()`
  - `Float(2)**(2**121)`
  - `gamma(Rational(10**40, 3)).evalf()`, `gamma(Float(1e40))`
- `factorial(10**40).evalf()` and `loggamma(10**40).evalf()` also stay
  unevaluated. SymPy does not finish them.
- `sin(exp(10**30).evalf())` stays symbolic. SymPy does not finish it.
- Float literals with a decimal exponent above about 4e35 raise
  `ValueError`, where SymPy returns a value. Up to about 1.7e38, the binary
  exponent does not fit the 2**120 counter
  (`Float('1e10000000000000000000000000000000000000')`). Beyond that, the
  decimal exponent itself overflows i128 (`Float('1e10**39')`-sized
  literals). `Float('1e1000000000000000000000000000000000')` is now
  evaluated.

Tests: `tests/test_bigexp_parity.py`, `tests/test_prec_parity.py`.

## Solve and integration leftovers (solint-rest, wave 7, 2026-10-23)

Implemented:
- `_tsolve`'s bivariate block (`python/symox/solvers/bivariate.py`, a port
  of `sympy.solvers.bivariate`, called from `_tsolve_steps.py` after the
  native NotImplementedError): the `sym**spow` change of variable, the
  `eq_down` factoring, `_solve_lambert` with the even-degree `t`
  substitution, the nsimplify replacement and `bivariate_type`.
  `solve(x*exp(x**2) - 1, x)` gives `[exp(-LambertW(2)/2)]`,
  `solve(x**a*exp(x) - 1, x)` gives `[a*LambertW(1/a)]`, and `x**x - 2`,
  `x*log(x) - 1`, `log(x) + x**2`, `x**2 - 2**x`, `exp(x) + x**3` and
  `x - log(x)*y` give SymPy's LambertW roots. The posify retry tries the
  rewritten equation first (as SymPy does), then the original one.
- solveset `radical_step`: a product of radicals without a symbol-free
  coefficient goes to `unrad` unchanged (the reference inverter does not
  invert it). `sqrt(x)*sqrt(x - 7) - 12` gives `{16}`, as do
  `2*sqrt(x)*sqrt(x - 7) - 24` and the Reals domain.
- Rational equations with symbolic poles keep SymPy's
  `_result - singularities` form: `solveset(y - b/(1 + a*x), x)` gives
  `Complement({b/(a*y) - 1/a}, {-1/a})` (SymPy writes the element
  `(b/y - 1)/a`). This also fixes a wrong result: `solveset((x - a)/(x - b),
  x)` gave `{a}`, which is wrong for `a = b`. It now gives
  `Complement({a}, {b})`.
- Real domains: an exact algebraic root whose realness is undecided goes
  through `expand_complex` (as in `_solve_as_poly`), and a provably nonzero
  imaginary part drops it. Roots with logarithms stay undecided. `solveset(x**3 + x + 1, x, S.Reals)` now gives
  the single real Cardano root instead of a union with an undecided
  Intersection.
- `roots` applies `preprocess_roots`'s `_integer_basis` rescaling to rational
  polynomials. The roots of `3125*x**3 + 37100*x**2 + 100800*x - 82944` now
  have SymPy's form.
- `limit` never returns `nan` computed from an expression without `nan`.
  Before this change, `limit(Ei(-x) - Ei(-2*x), x, 0, '+')` gave `nan`, and
  so did `Ei(x) - Ei(2*x)`, `Ci(x) - Ci(2*x)` and `Chi(x) - Chi(2*x)`: the
  substitution of the singular value `Ei(0) = -oo` gave `-oo + oo`. SymPy
  gives `-log(2)` here. The limit now stays unevaluated.
- `solve(x**(1/3) + sqrt(x) + 1, x)` gives SymPy's `[]`.
  `bivariate_type` no longer raises PolynomialError on its dependent
  generators.

Intentional divergences:
- `solve_univariate_inequality((I*x + 1)**3 < 2, x)` raises
  NotImplementedError, as before this wave. SymPy 1.14 returns
  `Interval(-oo, oo)`, but at x = 1 the left side is `-2 + 2*I`, which has
  no order. Now that `solveset((I*x + 1)**3 - 2, x, S.Reals)` decides
  EmptySet (as SymPy does), an inequality with `I` is solved only when its
  expanded difference is real for a real variable.
- `solveset(sqrt(x)*(x + 1)**(1/3) - 2, x)` gives
  `{CRootOf(x**5 + 2*x**4 + x**3 - 64, 0)}`, which was checked numerically.
  SymPy raises PolynomialError.
- In the three-radical cubic case in Complexes, symox rejects the complex
  candidates after a numerical check. SymPy keeps them in a ConditionSet.
- On Reals, `solveset(y - b/(1 + a*x), x, S.Reals)` and
  `solveset(a/(x - 1) - 2, x, S.Reals)` keep the pole Complement. SymPy
  returns `Intersection({(b/y - 1)/a}, Reals)`, which keeps the pole
  `x = -1/a` when `b = 0` (and `{1}` for `a = 0`). Those values are not
  solutions.

Remaining:
- The limits above (log-singular leading terms of Ei/Ci/Chi/li/expint)
  need the `_eval_as_leading_term` of those functions. Until then, the
  definite integrals that depend on them stay unevaluated:
  - the Frullani integrals `(exp(-x) - exp(-2*x))/x` and
    `(cos(x) - cos(2*x))/x` on (0, oo);
  - `log(log(x))/x**2` on (1, oo);
  - `sqrt(log(1/x))` on (0, 1). Its limit at 0+ of
    `x*sqrt(log(1/x)) + sqrt(pi)*erfc(sqrt(log(1/x)))/2` stays held.

  The antiderivatives are identical to SymPy's. SymPy's `meijerint_definite`
  also fails, so these integrals go through the limits.
- `SingularityFunction` integrals over (-oo, oo) stay unevaluated.
- The real inverter of negative bases is not ported. On Reals, SymPy
  gives `{3}` for `(-5)**(x + 1) = 625` and `{42}` for
  `(-3)**(x - 3) = -3**39`. symox keeps the complex candidate
  (`log(625)/log(-5) - 1`, expanded) in an undecided Intersection.
  `solveset(Eq((-2)**x, 4), x, S.Reals)` gives EmptySet where SymPy gives
  `{2}`. That false "no solution" predates this wave and remains open.
- `ratint(1/(x**3 + x + 1))` still raises (xfail in
  `tests/test_integ_rest_parity.py`). The step-by-step comparison of
  `ratint_logpart` was not done.
- `integrate(1/(1 + exp(x - 1/2)/(1 + exp(x))), x)` still times out.
- Form difference: the second root of `(a*x + b*y)**2 + (a*x + b*y) - 1`
  for `x` is `(-b*y - sqrt(5)/2 - 1/2)/a`. SymPy gives
  `(-2*b*y - sqrt(5) - 1)/(2*a)`.

Tests: `tests/test_solint_rest_parity.py`.

## simp-rest2 (wave 8)

Implemented:
- The real inverter of negative bases (SymPy's `_invert_real`, the
  `integer_log` branch), in `negative_base_step`
  (`python/symox/solvers/_solveset_steps.py`). On a real domain, the
  symbol-free terms and factors are moved to the right side, `exp` and
  `log` are inverted, and `b**g(x) = r` with `b < 0` is solved as
  `g(x) = k` for the exact integer `k` with `b**k = r`. A product of
  finite factors is split first, as `_solveset` does. Results now match
  SymPy: `(-2)**x = 4` gives `{2}` (it gave a false EmptySet),
  `(-5)**(x + 1) = 625` gives `{3}`, `(-3)**(x - 3) = -3**39` gives
  `{42}`, `(-2)**x = -8` gives `{3}`, and `(-2)**x = 2` gives EmptySet.
  `(-2)**x` is real only on the integers, where it is `(-1)**x*2**x`, so
  EmptySet is correct there.
- A guard (`negative_base_guard`). On a real domain, a native result for
  an equation that still has a power `b**g(x)` whose constant base is not
  known to be positive becomes a ConditionSet. The principal logarithm
  misses the integer exponents of a negative base: `c**x = 4` has the
  solution 2 for `c = -2`, and `(-2)**x*2**x = 16` has the solution 2.
  SymPy keeps a ConditionSet for a base of unknown sign or a complex
  base (`c**x`, `r**x`, `I**x`, `(1 + I)**x`). Before, symox returned an
  incomplete Intersection, `{0}` or EmptySet.
- `solve` adds SymPy's `_tsolve` guess `e_rat`: `(-2)**x - 4` gives the
  integer solution 2 before the principal logarithm root (native
  `power_generator_tsolve`).
- `nthroot` is a port of SymPy's `nthroot` and `_nthroot_solve`.
- `sqrtdenest`: ports of `_sqrt_numeric_denest` and `_subsets` (issue
  19914).
- `minpoly`: ports of `_is_sum_surds`, `_separate_sq`, `_choose_factor`
  (over QQ) and `_minimal_polynomial_sq`.

Intentional corrections (each value checked by substitution):
- SymPy's `integer_log` raises ValueError for a non-integer right side or
  base. symox solves exactly:
  - `(-2)**x = 1/4` gives `{-2}`;
  - `(-1/2)**x = 4` gives `{-2}`;
  - `(-2/3)**x = 81/16` gives `{-4}`;
  - `(-E)**x = E**2` gives `{2}`;
  - `(-sqrt(2))**x = 4` gives `{4}`.
- SymPy raises "base must be 2 or more" for base -1. symox gives the
  lattices `2*Integers` (for `1`) and `2*Integers + 1` (for `-1`),
  intersected with a bounded domain (`{0, 2, 4, 6}` on `[0, 7]`).
- `(-2)**x*(x - 1) = 0` gives `{1}`. SymPy raises.
- `(-2)**x = 2**x` stays a ConditionSet. SymPy returns EmptySet, but the
  equation holds on every even integer.
- A symbolic or complex right side, a Float, a symbolic negative base, and
  equations with several negative-base powers stay ConditionSets where
  SymPy raises ValueError.
- `solve((-2)**xr - 4, xr)` for a real `xr` gives `[2]`. SymPy also keeps
  a non-real second root.

Remaining:
- `log(exp(x), 2, evaluate=False)`: the native `log` takes one argument,
  so the held two-argument form is a generic function. As a result,
  `simplify(..., inverse=True)` does not give `x/log(2)`.
- `Derivative(sin(x), (x, d))` with a symbolic integer order (issue
  22210): native derivative orders are machine integers, so construction
  raises ValueError.
- Still not done:
  - `nc_simplify` of noncommutative matrix products (and issue 25590);
  - the hyperexpand timeouts and forms;
  - `trigsimp_groebner`;
  - `powsimp_polar`;
  - `posify` with Integral/Matrix;
  - issues 17141 and 11004.
- `nthroot1`: the native `factor_list` stops with "factorization resource
  limit", and `nthroot` then returns its argument.
- hypersimp and simplify forms: `simplify((x + 5)/(16*x**2 + 28*x + 6))`
  keeps the expanded denominator. SymPy gives `(x + 5)/(2*(8*x**2 + 14*x
  + 3))`. The `factor_terms` step of native simplify keeps the reference
  form on a count_ops tie.
- `(I**x).is_finite` is None for a real finite `x` (SymPy: True). As a
  result, `I**x*(x - 2) = 0` on Reals stays a ConditionSet where SymPy
  gives `{2}`.
- SymPy rewrites `(r**2)**x` as `Abs(r)**(2*x)` for a real `r`. symox
  does not.

Tests: `tests/test_simp_rest2_parity.py`.

## Matrix leftovers, part 2 (mat-rest2, wave 8, 2026-10-23)

Ported from SymPy 1.14:
- `MatrixBase.pow` and `__pow__` for non-integer and symbolic exponents
  (`python/symox/_matrix_power.py`): Jordan-block powers, the
  recursion/cayley/multiply methods, `NonInvertibleMatrixError` for a
  singular matrix to a non-integer power, and an unevaluated `MatPow` when
  the Jordan form cannot be used. `M**S.Half` and `M**n` now match SymPy.
  The `_eval_matrix_exp_jblock` and `_eval_matrix_log_jblock` helpers are
  ported too.
- `pinv(method='ED')` uses SymPy's `_pinv_diagonalization`
  (`python/symox/_matrix_pinv.py`). The old native path hung on symbolic
  matrices.
- The deprecated mixin classes that `test_commonmatrix` imports:
  `symox.matrices.common`, `symox.matrices.matrices` and the `_legacy_*`
  modules, plus `matrices.exceptions` and `matrices.normalforms`.
- `diff(expr, A)` for a matrix or array variable, and the `[[x, y]]` and
  `([x, y], 2)` forms, go through `derive_by_array` as SymPy's
  `ArrayDerivative` does. The array keeps the mutability of the matrix.
  A scalar derived by a matrix stays an ImmutableMatrix
  (`_call_derive_scalar_by_matrix`).
- `Mod(Matrix, n)` applies `Mod` to each entry (`_eval_Mod`).
- `Float(10.0) % 1` is the exact `S.Zero`, as in `Float._new`. A Float
  or Rational modulus keeps a Float zero.
- `MatrixExpr.match` with a `Wild`, including `exclude`.
- `Float.epsilon_eq`.
- A geometry `Point` is accepted as a matrix row, because SymPy's
  `is_sequence` accepts it.
- `DeferredVector` is not iterable.
- `expand_mul`, `expand_func`, `expand_trig`, `expand_log` and
  `expand_power_exp` on a matrix return an ImmutableMatrix, as in SymPy.
- `Matrix._rep` is a sparse `DomainMatrix` over ZZ, QQ or EXRAW, and
  `_fromrep` builds a matrix from one. Both exist on the dense, immutable
  and sparse classes.

Remaining gaps:
- ImmutableMatrix is not a real class: `type(ImmutableMatrix(...))` is
  `symox._core.Matrix`, and the Hashable/type-name checks of the SymPy
  tests fail (41 probe statements).
- `Eq(ImmutableMatrix, ...)` raises "cannot convert Matrix", because
  relational nodes cannot hold matrices.
- Quantum `Dagger`/`Operator` entries become native noncommutative nodes,
  so the matrix adjoint gives `conjugate(...)` forms.
- `adjoint`/`transpose` of `Abs(A)` and `sqrt(A)` for a noncommutative `A`
  are not ported.
- Poly entries in a Matrix are converted to Expr without the
  "deprecated-non-expr-in-matrix" warning. This makes `test_smith_normal`
  crash in the probe harness, although `smith_normal_form` itself works.
- Derivative variables are not sorted as in SymPy (core Rust).
- Matrix expression derivatives of traces fail with "Shape should contain
  integers only" (`expressions/tests/test_derivatives.py`).
- A `_fromrep` with an algebraic-field domain keeps that domain on the
  dense classes only. The sparse classes recompute the domain (EXRAW).
- Unsimplified forms still differ from SymPy in `cholesky`,
  `GramSchmidt`, `LDLdecomposition`, `matrix_norm`, `exp`, `refine`,
  `simplify(ratio=...)` and the `numer` expand hint. The values are
  equal. The SymPy ED pinv form is also unsimplified, and
  `tests/test_sparse_matrices.py` now compares it after `simplify`.
- The order of the cubic roots of `roots(..., multiple=True)` differs.
- `DFT`, `GenericIdentity` and `_csrtodok` are missing.

Probe (all of `sympy/matrices/tests` and `expressions/tests`, mat-rest
harness): before same 3641, diff 301, exc 389, not reached 1801;
after same 4174, diff 320, exc 245, not reached 1393, with no statement
going from same to different. The new diffs are statements that used to
raise or were not reached. They are the type-name and form gaps above.

Tests: `tests/test_mat_rest2_parity.py`.

## Wave-8 replay of the documentation examples (replay, 2026-10-23)

The docstring audit (`tools/sympy_parity.py examples`, ordinary examples and
`--keywords`) was replayed on all 2916 documented objects. The triage is in
`audits/sympy-1.14.0/triage-wave8.json` and `.md`, compared with
`triage-wave4-final.json`.

Regressions found with the a2a0038 build (waves 5-7), all fixed:

- `FunctionMatrix`: `(A**2)[i, j]` of a symbolic-size matrix expression timed
  out. `MatPow(b, 2)` entries are now one row times one column, as in
  SymPy's `MatMul._entry` (`crates/symox-core/src/matrixexpr.rs`).
- `BinaryRelation`: `str()` of an applied user relation raised
  AttributeError. `__str__` is now `__repr__` (`Q.R(2, 1)`).
- `gaunt`: the `n(64)` example timed out in `perfect_power` on a large
  factorial quotient. Big integers are now trial-divided by the primes below
  1000. The multiplicity of the first prime factor bounds the candidate
  exponents, as in SymPy's `perfect_power(factor=True)`
  (`crates/symox-core/src/ntheory/factor.rs`).
- `PIDController.eval_frequency` (keyword example): `expand()` rebuilt the
  held `1/(0.99 + 0.2*I)` of `TransferFunction.to_expr` and evaluated it.
  SymPy's multinomial hint leaves `base**-1` alone, and no other hint
  applies, so that node is now returned unchanged
  (`crates/symox-core/src/expand.rs`).

Documented, not changed:

- `stats.Zeta` keyword example (b): `variance(X)` comes back over a common
  denominator. The polynomial-moment path simplifies its sum, while SymPy
  uses the moment-generating function. The value is the same.
- `solve` (e, new): `(root(r**15, 3) - root(r**15, 5) + Rational(1, 7)).n(2)`
  with `r = CRootOf(7*x**5 - 7*x**3 + 1, 1)` is exactly 0. SymPy reports
  `0.e-110`, while symox prints spurious digits. The cancellation budget
  and the residue logic of `eval_sum` (`crates/symox-core/src/evalf_mp.rs`)
  only apply at the root, but here `abs()` wraps the sum.

The original a2a0038 replay outputs were lost in a power cut. Both modes
were therefore re-run on the build that has the first three fixes. The
objects whose examples call `expand`, and the objects flagged in that run,
were re-run on the final build.

Counts of the still-differing ordinary examples, wave-4 final to wave 8: a 89
to 89, b 68 to 66, c 30 to 13, e 0 to 1, f 15 to 5, h 21 to 20, i 33 to 34,
and 256 to 228 in total. 28 objects are fixed since the last triage. Two
objects are `reference_timeout`: SymPy itself exceeded the 30 s budget on the
loaded machine while symox completed every example (`ode_sol_simplicity`
and the keyword examples of `ReferenceFrame.orient_space_fixed`). The ODE
corpus was not run, because its reference was still being regenerated after
the power cut.

Tests: `tests/test_replay_parity.py`.

## Arbitrary-size binary exponents (bigexp2, wave 8, 2026-10-23)

Fixed:
- The binary exponent of `BigFloat` is now an `Exponent`
  (`crates/symox-core/src/bigfloat/exponent.rs`): an `i128` up to 2**120 in
  absolute value, otherwise a boxed `BigInt`. The i128 path is unchanged
  for ordinary values. SymPy's values are now returned for
  `exp(10**40).evalf()`, `exp(10**100).evalf()`, `exp(10**1000).evalf()`,
  `Float(2)**(2**121)`, `Float(2)**(2**200)`, `Float(3)**(10**40)`,
  `exp(10**40 + I).evalf()`, Float literals such as `Float('1e10**50')`
  (decimal exponent above 4e35), and arithmetic, comparisons, hashing,
  `_mpf_` and `srepr` on these values.
- Materialized exponents are bounded at `MAX_EXPONENT_BITS = 65536` bits,
  which keeps decimal printing and repeated squaring bounded.
- `floor`, `ceiling` and `int` of a Float whose integer part needs more
  than 2**31 bits (`MAX_INTEGER_SHIFT`) raise
  `OverflowError: too many digits in integer`, as in SymPy. Before this,
  `floor(Float(2)**(2**100))` returned 1 and the same call at 2**121
  panicked in the core.
- `ln` carried too few guard bits after its square-root reduction: the
  last digits of `log(2).evalf(2500)` were wrong, and so was
  `Float(2)**(2**k)` for k >= 3000 (mantissa 2**53 - 1 instead of 1). The
  reduction now carries `reduce + 4` extra bits.
- `factorial(10**40).evalf()` = `2.68554168269385e+395657055180967481723488710810833949177076`
  and `loggamma(10**40).evalf()` = `9.11034037197618e+41` (mpmath). SymPy
  does not finish either in 30 s.

Intentional corrections, checked against mpmath at 60 digits:
- `gamma(Rational(10**40, 3)).evalf()` gives
  `3.55973154373079e+130295280877923619116846143926094265361665`. SymPy
  rounds the argument first and gives `2.99850408006855e+...539429`.
- `gamma(Float(1e40))` ends in `...626811080`. SymPy's `...626811120` is
  wrong in the exponent.
- `gamma(Float(1e100))` gives `5.54790221218382e+...443467`, the 53-bit
  rounding of the true value. SymPy gives `5.54790221218381e+...443567`.
- `exp(10**40).evalf()**Rational(1, 3)` gives
  `4.53613542926576e+1447648273010839425503763063055350274314`. SymPy
  rounds 1/3 to 53 bits and gives `9.34867898800170e+...839853`.

Benchmarks (load average about 25 during the runs, two interleaved runs
each, base commit a2a0038):
- The 15 evalf/precision test files: 41.0 s and 32.2 s on the base, 43.3 s
  and 32.8 s here.
- `test_performance_regressions.py` and `test_performance_validation.py`:
  20.3 s and 13.0 s on the base, 16.3 s and 17.2 s here.
- Common-path microbenchmark (best of 9): Float arithmetic loop 54.9 ms on
  the base, 56.2 ms here; `log(n).evalf(500)` 126.6 ms and 120.1 ms;
  `gamma(Rational(n, 5)).evalf()` 624 ms and 657 ms. The differences are
  within the noise of this load.

Remaining (explicit, never a wrong number):
- Binary exponents of 65536 bits or more stay unevaluated
  (`exp(2**70000).evalf()`, `Float(2)**(2**70000)`), where mpmath has no
  bound. `Float(2)**(2**70000)` takes a few seconds to give up, because it
  computes the logarithm before checking the bound.
- Printing a Float whose exponent is close to the bound is slow (0.3 s for
  `Float(2)**(2**8000)`, more than 2 minutes near 2**65000). SymPy cannot
  print these values at all past about 4300 exponent digits (Python's
  int-to-str limit raises ValueError).
- `sin(exp(10**30).evalf())` stays symbolic. SymPy does not finish it.

Tests: `tests/test_bigexp2_parity.py`. Older tests that expected the
2**120 bound to hold values now check the evaluated values
(`tests/test_bigexp_parity.py`, `tests/test_float_exponent_arithmetic.py`,
`tests/test_precision.py`, `tests/test_prec_parity.py`).

## solint-rest2 (wave 9)

Implemented:
- Precision of cancelled sums with CRootOf. The `Expr.evalf` wrapper used to
  replace each CRootOf by a Float at `prec + 64` bits before the native
  evaluation, which left a rounding residue: for
  `r = CRootOf(7*x**5 - 7*x**3 + 1, 1)`,
  `(root(r**15, 3) - root(r**15, 5) + Rational(1, 7)).n(2)` gave
  `-6.6e-24`. The wrapper now first tries the native evaluator, which
  refines the root through its enclosures, and keeps that value when no
  CRootOf is left. The native evaluator also propagates a cancelled nested
  sum: the enclosing function or power is evaluated again at `prec + 64`
  and each part keeps only the digits on which both evaluations agree.
  This gives a scaled zero (`0.e-111`, SymPy `-0.e-110`) for the sum and
  inside `re`, `Mul`, `Pow`, `pi*f`, `f*(f + 1)` and `sin(f)`; `Abs(f)` and
  `im(f)` give 0 as in SymPy, and `exp(f)`, `f + 1` and `cos(1 + f)` give
  their correct values. The sign and exponent of SymPy's zero are rounding
  noise and are not reproduced (`(f**2).n(2)` is `0.e-282` here and
  `6.1e-238` in SymPy, and `sqrt(f).n(2)` is `0.e-70` here and
  `5.6e-59*I` in SymPy).
- Log-singular series of `Ei`, `Ci`, `Chi`, `li` and `expint(1, z)` at a
  zero of the argument: `EulerGamma + log(u) + sum` with the convergent
  series (minus `I*pi` for `Ei` of an argument with a negative real
  part; `li(u)` is `Ei(log(u))` near `u = 1`). The leading terms are
  `log(x) + EulerGamma` as in SymPy. SymPy raises NotImplementedError for
  the leading terms of `li` and `expint`; symox returns
  `-log(x) - EulerGamma` for `expint(1, x)`. The argument is expanded
  until it is nonzero, so `Ei(2*x)`, `Ei(-x)` and `Ei(x**2)` are covered
  at every order. The limits `limit(Ei(-x) - Ei(-2*x), x, 0, '+')`, and the
  same for `Ci` and `Chi`, give `-log(2)` as in SymPy. The Frullani
  integrals `integrate((exp(-x) - exp(-2*x))/x, (x, 0, oo))` and
  `integrate((cos(x) - cos(2*x))/x, (x, 0, oo))` give `log(2)` as in SymPy.
- Definite integrals of SingularityFunction with an infinite endpoint.
  They use the antiderivative, with each `SingularityFunction` replaced by
  its value at the endpoint (0 at `-oo`; at `oo` 0 for negative orders, 1
  for order 0 and `(x - c)**n` for positive orders, then a limit). An
  endpoint at a singularity uses the one-sided value. The results match
  SymPy for negative orders: `integrate(SingularityFunction(x, 0, -1),
  (x, -oo, oo))` is 1, with `x**2*SingularityFunction(x, 2, -1)` it is
  4, and `SingularityFunction(x, 0, -2)` gives 0.
  Divergence: SymPy gives 0 for `integrate(SingularityFunction(x, a, 0),
  (x, -oo, oo))` and for `SingularityFunction(x, 1, 1)`, which is wrong.
  Those integrands are a unit step and a ramp, so both integrals diverge,
  and symox returns `oo`. Over `(0, oo)` symox gives
  `oo - SingularityFunction(0, a, 1)` where SymPy gives
  `-SingularityFunction(0, a, 1)`.
- `ratint(1/(x**3 + x + 1), x)`. `log_to_atan` divided by a constant EX
  polynomial with `Poly.div`. SymPy's EX elements are canonicalised by
  `cancel`, so the remainder is exactly 0 there. symox's EX arithmetic does
  not recognise `(a/c)*c - a` as 0 for these nested cube roots, and the
  division raised PolynomialDivisionFailed. A constant divisor now divides
  the coefficients directly. The result is real and has the same value as
  SymPy's (checked by differentiation and against the definite integral
  over (1, 2)). It is written in a different, larger radical form (82704
  characters, SymPy 16653).
- `solveset` of `sec`, `csc` and `cot` equations uses SymPy's
  `_invert_trig_hyp_complex` lattices when the inverse stays unevaluated:
  `2*n*pi +- asec(k)`, `2*n*pi + acsc(k)` and `2*n*pi + pi - acsc(k)`, and
  `n*pi + acot(k)`. `solveset(sec(2*x) - 3, x)` gives
  `Union(ImageSet(Lambda(_n, _n*pi - asec(3)/2), Integers),
  ImageSet(Lambda(_n, _n*pi + asec(3)/2), Integers))` as in SymPy, where it
  used to give the equivalent `acos(1/3)` form. Over the reals this form
  is used only for a real argument and a real `k` outside `(-1, 1)`, where
  every lattice point is real. The old answer was an `Intersection` with
  `Reals` and is now the plain union, as in SymPy. Equations whose inverse
  evaluates, such as `sec(x) - sqrt(2)`, keep the reduction to `cos`, which
  already matched SymPy.
- Systems with symbolic coefficients: 14 linear and nonlinear systems,
  among them `solve([a*x + b*y - 1, c*x - y], [x, y])`,
  `linsolve([a*x + b*y - c, d*x - y], [x, y])`,
  `solve([x*y - a, x + y - b], [x, y])` and
  `nonlinsolve([x*y - a, x - b*y], [x, y])`, give the same results as
  SymPy. `solve(sec(2*x) - 3, x)` already matched.

Remaining (explicit, never a wrong value):
- `integrate(log(log(x))/x**2, (x, 1, oo))` (SymPy `-EulerGamma`) and
  `integrate(sqrt(log(1/x)), (x, 0, 1))` (SymPy `sqrt(pi)/2`) stay
  unevaluated. The antiderivatives are found, and the limit at `x = 1` of
  `Ei(-log(x)) - log(log(x))/x` is `EulerGamma`. The native gruntz has no
  asymptotic expansion of `Ei`, `Ci`, `Si` or `erfc` at an infinite
  argument, so `limit(Ei(-log(x)), x, oo)`,
  `limit(erfc(sqrt(log(1/x))), x, 0, '+')` and even `limit(Si(x), x, oo)`
  stay as `Limit` objects. The base already had this gap. SymPy uses the
  tractable rewrites with `_eis` and `_erfs`.
- `limit(li(x) - li(2*x), x, 0, '+')` (SymPy 0) stays a `Limit`: near 0
  `li(x)` is `Ei(log(x))` with `log(x) -> -oo`, which needs the same
  asymptotic expansion. SymPy raises NotImplementedError for
  `li(x).as_leading_term(x)`, and so does symox.
- `integrate(1/(1 + exp(x - 1/2)/(1 + exp(x))), x)` still does not finish.
  SymPy solves it through risch in 0.2 s, while its manualintegrate also
  loops. In symox risch stops with CoercionFailed. The coefficients over
  `ZZ(exp(1/2))` contain `E`, because `exp(1/2)**2` collapses to `E`, and
  `infer_parameters`, `coefficient` and `checked_fraction_pair` in
  `crates/symox-core/src/poly/fraction.rs` treat `E` as a separate
  parameter. As a result `Poly(exp(1/2)/(1 + E + 2*exp(1/2)), t,
  domain='ZZ(exp(1/2))')` fails, while SymPy accepts it. The fix is to map
  `E` and `exp(r)` onto powers of the `exp(1/m)` generator in the fraction
  field conversion.

Tests: `tests/test_solint_rest2_parity.py`. The ratint test in
`tests/test_integ_rest_parity.py` no longer hits its xfail branch.

## simp-rest3 (wave 9, 2026-10-23)

Implemented:
- `(I**x).is_finite` is True for a real `x`: a power of a finite nonzero
  base to a finite exponent is finite (`finite_power` in
  `assumptions.rs`). `solveset(I**x*(x - 2), x, S.Reals)` gives `{2}`.
- `(r**2)**x` is `Abs(r)**(2*x)` for a real `r`, and also for a symbolic
  even exponent (Pow-of-Pow rule in `ops.rs`), as in SymPy.
- nthroot1: Zassenhaus' Hensel lifting and recombination multiply
  integers of the size of the factor bound p**l, and the fixed work
  allowance of the factor budget (2e7 units, about 5 ms) stopped them on
  the degree-20 polynomial of `nthroot` (coefficients of 351 to 2011
  bits). `zassenhaus` now extends the allowance by the cost of that stage,
  r*(n + 1)**2*bits(p**l)**2/64 units per lifting level, at most 100 times
  the scope's limit (`factor_budget::grant`); the attempt, degree,
  coefficient and storage limits are unchanged and still bound the
  exponential stages. The heuristic gcd also stops before an evaluation
  point whose images would pass the coefficient-bit or work limits and uses
  the modular gcd instead (`factor_budget::affordable`). The factorization
  takes 22 ms (1e-20 case) and 31 ms (1e-30 case), SymPy 92 ms, and
  `nthroot(expand_multinomial(q**5), 5) == q` for
  `q = 1 + sqrt(2) + sqrt(3) + 1/10**20` and `1/10**30`.
- `nc_simplify` of matrix expressions is a port of SymPy's
  (`python/symox/simplify/_nc_simplify.py`), for example
  `A*B*A*B*A*B*C*(A*B)**3*C` gives `((A*B)**3*C)**2` (issue 25590).
- `trigsimp_groebner` is a port of SymPy's
  (`python/symox/simplify/_trigsimp_groebner.py`), with `hints`,
  `polynomial` and `quick`. `trigsimp(method='groebner')` and
  `method='combined'` use SymPy's `groebnersimp` traversal. `groebner` and
  `ratsimpmodprime` accept the `gens` keyword (OptionError when both forms
  are given, as in SymPy). Intentional correction: where SymPy's
  `analyse_gens` raises ValueError (`1/2 is not an integer`, angles that
  are rational multiples such as `sin(beta/2)` with `sin(beta)`),
  `method='groebner'` and `'combined'` fall back to the native
  common-angle Groebner reduction, which proves
  `-6*sin(beta/2)**8 + 12*sin(beta/2)**6 - 6*sin(beta/2)**4 + 3*sin(beta)**4/8`
  is 0, as symox did before.
- powsimp of polar bases: `exp_polar(a)` has the base `exp_polar(1)` and
  the exponent `a`, and polar bases combine like nonnegative bases.
  `powsimp(exp_polar(x)*exp_polar(y)*p**x*p**y)` gives
  `(p*exp_polar(1))**(x + y)` for a polar `p`.
- posify of a Matrix, and `expand` of an Integral or a Sum of a sum with
  limits not known to be infinite (SymPy's
  `AddWithLimits._eval_expand_basic`): `expand(Integral(x + y, (y, 1, 3)))`
  is `Integral(x, (y, 1, 3)) + Integral(y, (y, 1, 3))`, and an infinite
  bound keeps the Integral.

Remaining (explicit):
- `simplify((x + 5)/(16*x**2 + 28*x + 6))` stays
  `(x + 5)/(16*x**2 + 28*x + 6)`; SymPy gives
  `(x + 5)/(2*(8*x**2 + 14*x + 3))` (same count_ops, its
  `shorter(together(expr), ...)` puts together first). A tie rule that
  takes SymPy's rational state in the native selection did not reach the
  output (later steps of the selection cancel again); it was reverted. The
  test is an xfail.
- test_lerchphi of the probe: the symbolic statements are correct and
  fast, but the `can_do` numerical checks of the same function pass the
  60 s probe timeout.
- `trigsimp_groebner(ex)` of SymPy's test can return
  `(-5*sin(x) + 5*cos(x) - 1)/(-8*sin(x) + 6*cos(x))`, the same quotient
  with both signs flipped. SymPy's test accepts both forms.
- `((1 + I)**x).is_zero` is None (SymPy: False), because `is_zero(1 + I)`
  is not decided by the exact queries.
- `MatrixExpr.doit(inv_expand=False)` is not supported natively;
  nc_simplify falls back to `doit()`.
- Issue 17141: `simplify(acos(I)**2)` gives a different form from SymPy's
  `(pi - 2*I*log(1 + sqrt(2)))**2/4`. It needs a change of the
  factor_terms tie rule of the simplify selection.
- Issue 11004: cancel converts Floats to 53 bits.
- `log(exp(x), 2, evaluate=False)` held: SymPy re-evaluates the rebuilt
  log in subs, doit, expand and simplify (`subs(x, 3)` gives
  `3/log(2)`, `simplify(e, inverse=True)` gives `x/log(2)`). Symox keeps
  the held form.
- `Derivative(sin(x), (x, d))` with a symbolic order (issue 22210): the
  native Derivative stores integer orders, so a symbolic order needs a
  core change.

Tests: `tests/test_simp_rest3_parity.py`.

## Remaining c/f objects (audit9, wave 9, 2026-10-23)

Replay of the 16 remaining c and f objects of `triage-wave8`: 16
observation_difference before, 5 examples_match and 11
observation_difference after (`audits/sympy-1.14.0/triage-wave9.md`).

Implemented:
- Galois groups of degree 5 and 6 (`galois_group`, `Poly.galois_group`,
  `AlgebraicField.galois_group`), ported from SymPy with its resolvent
  lookup table. Identical to SymPy on 29 test polynomials.
- Parameter derivatives of `jacobi` (a and b), `gegenbauer` (a) and
  `assoc_laguerre` (alpha), as Sum expressions.
- `harmonic(oo, m)` for symbolic m: zeta(m) when m - 1 is positive, oo when
  it is nonpositive.
- Matrix arguments of undefined functions and of Expr subclasses (code
  printers with matrix user functions, `shape` dispatch).
- `nsolve`: the secant and bisection solvers follow mpmath's findroot
  (start step, stop rule, step limits, midpoint result, error messages).
  Complex multiplication and division round as mpmath's mpc_mul/mpc_div.
- `SymPyDeprecationWarning`, `sympy_deprecation_warning`, `ignore_warnings`
  and `warns`/`warns_deprecated_sympy` ported (message, repr, once filter,
  stacklevel check).
- `lambdify` raises SymPy's ImportError for a missing module
  (tensorflow and the other string modules).

Intentional divergences (SymPy is wrong, category h):
- `diff(assoc_laguerre(n, a, x), a)`: SymPy gives
  Sum(L_k^a(x)/(n - a)); the correct sum divides by n - k (checked against
  the explicit polynomials for n = 2, 3, 5).
- `lambdify` of a Piecewise at x = 0: SymPy raises ZeroDivisionError because
  numpy evaluates every branch; symox returns the selected branch value.
- Unchanged from earlier waves: `Sum.euler_maclaurin`, `approximants`,
  `combine_kronecker`, `ImplicitRegion.regular_point`.

Remaining (unsupported, raised or left unevaluated, never a wrong value):
- `limit(harmonic(n, m + 1), n, oo)` stays unevaluated (SymPy uses the
  limit of the Hurwitz zeta rewrite).
- `DiscreteMarkovChain` with `Str` states and `MultivariateNormal` with
  MatrixSymbol parameters raise NotImplementedError.
- `lambdify`: `inspect.getsource` of the compiled function, numpy scalar
  result types (`np.float64` and `np.True_` versus Python floats).
- `nsolve`: a real part of about 1e-61 instead of 0 for
  `nsolve(x**2 + 2, I)`, and the last digits or the format of the residual
  in some "Could not find root" messages.

Tests: `tests/test_audit9_parity.py`.

## mat-rest3 (wave 9): dense matrix classes, matrix Eq, CSR helpers

Done:
- `MutableDenseMatrix` and `ImmutableDenseMatrix` are real classes. Both
  derive from the native `DenseMatrix` (module `symox.matrices.dense`).
  The native kernels allocate every result as one of the two classes,
  chosen by mutability, so `type()`, `isinstance`, `issubclass` and the
  class names match SymPy (`Matrix is MutableDenseMatrix`,
  `ImmutableMatrix is ImmutableDenseMatrix`, an immutable matrix is not a
  `MutableDenseMatrix`). Internal modules keep `isinstance` checks
  against the native base, so they still accept both classes.
- `eye`, `zeros` and `ones` (functions and static methods) accept SymPy's
  `cls=` option. `_eval_eye`, `_eval_zeros`, `_eval_ones`, `_eval_diag`,
  `_eval_jordan_block` and `_eval_wilkinson` exist on both dense classes,
  for `SpecialOnlyMatrix`-style callers.
- `Eq` and `Ne` with explicit matrices follow
  `_eval_is_eq(ImmutableDenseMatrix, ...)`: different shapes or a
  non-matrix side give `False`. Otherwise the entries decide: all zero
  differences give `True`, any proven non-zero difference gives `False`,
  and an undecided entry keeps the relation unevaluated (it re-evaluates
  after `subs`). A `MatrixSymbol` against an explicit matrix of the same
  shape stays unevaluated, as in SymPy.
- `symox.matrices.sparsetools._doktocsr` and `_csrtodok` are ported
  (including SymPy's ValueError for an all-zero matrix in `_doktocsr`).
- `Trace(X).arg` returns the traced matrix, as in SymPy.
- Internal checks that meant "any explicit dense matrix" now test the
  native base (`DenseMatrix`) instead of the public `Matrix`, so they still
  accept immutable matrices (stats, quantum, codegen, optics, simplex).

Remaining (documented gaps):
- `Eq(M, N).lhs` for explicit matrices is a matrix-expression proxy that
  prints as `ImmutableMatrix`; SymPy gives an `ImmutableDenseMatrix`.
- `Matrix([[a, 1]]).adjoint()` with a quantum `Operator('a')` gives
  `conjugate(OrderedProduct(Operator(a)))`; SymPy gives `Dagger(a)`.
- For a noncommutative `A`: `adjoint(Abs(A))` and `transpose(Abs(A))` stay
  as `transpose(Abs(A))`, where SymPy gives `Abs(A)` (Abs is real).
  `adjoint(sqrt(A))` stays `Dagger(sqrt(A))`, where SymPy expands it
  through `re`/`im`/`atan2`. `transpose(sqrt(A))` matches.
- `Trace(X*X).diff(X)` returns an explicit array of entries where SymPy
  gives `2*X.T`.
- Sparse `_rep.domain` of `SparseMatrix([[sqrt(2), 1], [0, 1]])` is
  `EXRAW` in both libraries for this input; an algebraic-field domain
  kept through `_fromrep` was not reproduced or changed here.
- `DFT`, `IDFT`, `GenericIdentity` and `GenericZeroMatrix` are still
  missing: they need new native matrix-expression node kinds.
- `Poly` entries in a matrix are converted to expressions
  (`Matrix([[Poly(x)]])` gives `Matrix([[x]])`), without SymPy's
  deprecation warning; `test_smith_normal` still crashes.
- `simplify(SparseMatrix(...))` returns a `MutableSparseMatrix`; SymPy
  returns an `ImmutableSparseMatrix` (dense matrices match).
- `BlockDiagMatrix`, `MatMul` and other matrix expressions still report
  the native `_MatrixExpr` class from `type()` and `__class__`.

Tests: `tests/test_mat_rest3_parity.py`.

## liminf (wave 10): limits of special functions at infinite arguments

Implemented:
- Gruntz preparation rewrites a special function whose argument tends to
  +-oo the way SymPy's `_eval_rewrite_as_tractable` does: an elementary
  part (exponentials, oscillating sin/cos, constants) times an auxiliary
  function that has an asymptotic series in 1/z (`_erfs`, `_eis` and the
  `_eval_aseries` of Si/Ci, expint and the Fresnel integrals). This covers
  erf, erfc, erfi, Ei, li and Li (also at 0+), Si, Ci, Shi, Chi,
  expint(nu, z), E1, uppergamma, lowergamma, fresnels and fresnelc
  (`series/asym.rs`, `limits/gruntz.rs`). An auxiliary function never
  escapes: a limit value that still contains one is left unevaluated.
- The bounded heuristic (oscillating sin/cos replaced by parameters) runs
  on the prepared expression, so `limit(Si(x) + erf(x), x, oo)` and the
  other Si/Ci/Fresnel combinations get the rewrite's sin/cos.
- Series of Si, Shi, erfi, fresnels and fresnelc at 0 and at finite points
  through the derivative integral (`ser_entire_integral`). The results
  match SymPy term for term, e.g. `limit((Si(x) - x)/x**3, x, 0) = -1/18`.
- `integrate(log(log(x))/x**2, (x, 1, oo)) = -EulerGamma`,
  `integrate(sqrt(log(1/x)), (x, 0, 1)) = sqrt(pi)/2` and
  `limit(li(x) - li(2*x), x, 0) = 0`, as in SymPy.
- Fraction-field coefficients follow `FracField._rebuild_expr`: a power
  `b**e` is read as `g**(e/eg)` for a generator `g = b**eg`, so over
  `ZZ(exp(1/2))`, `E` is `exp(1/2)**2`. This is only a fallback after the
  direct conversion fails. With it, risch sees a single transcendental in
  `integrate(1/(1 + exp(x - 1/2)/(1 + exp(x))), x)`, which used to never
  finish and now returns SymPy's result in a few hundredths of a second.
- Probe of 71 limits: 15 equal to SymPy before, 61 after. Of the other
  10, 9 are the corrections below and 1 is still unevaluated.

Intentional divergences (SymPy 1.14 is wrong or unevaluated; the values
were checked with mpmath at x = 200):
- `limit(expint(n, x), x, oo)` for n = 1, 2, 3 is 0. SymPy gives
  `expint(1, oo)`, `-oo*I` and `oo*I`.
- `limit(Chi(x)*x*exp(-x), x, oo)` and `limit(Shi(x)*x*exp(-x), x, oo)`
  are 1/2. SymPy gives 0.
- `limit(expint(a, x)*x*exp(x), x, oo) = 1` for a positive. SymPy gives
  `oo*sign(expint(a, 1/_w))`.
- SymPy leaves these unevaluated, and symox gives the value:
  `x*exp(x)*expint(1, x)` -> 1, `erfc(x)/erfc(2*x)` -> oo and
  `Shi(x) - Chi(x)` -> 0.
- Equivalent form: `limit(lowergamma(a, x), x, oo)` gives `gamma(a)`.
  SymPy gives `lowergamma(a, oo)`.

Remaining (left unevaluated, never a wrong value):
- AccumBounds limits of an oscillating rewrite times a growing factor:
  `limit(x*Ci(x), x, oo)`, `limit(x*(Si(x) - pi/2), x, oo)` (SymPy:
  AccumBounds(-1, 1)) and `limit(x*(fresnels(x) - 1/2), x, oo)` (SymPy:
  AccumBounds(-1, 1)/pi).
- Unchanged and not from this wave: `limit(x*sin(x), x, oo)` gives
  AccumBounds(-oo, oo). SymPy gives `oo*sign(AccumBounds(-1, 1))`.

Tests: `tests/test_liminf_parity.py`.

## Integration corpus residue (icorpus, wave 10)

Integration corpus (`tools/integration_corpus.py`) before -> after:
SAME 1134 -> 1161, EQUIV 46 -> 21, COND 2 -> 0, MISSING 18 -> 18,
DIFF 10 -> 10. Every change is EQUIV/COND -> SAME.

Fixed:
- `solve_univariate_inequality` with a generator that has sign assumptions
  (positive, negative, ...): SymPy keeps a real generator, but its root
  finding over `S.Reals` ignores the other assumptions. symox used them and
  turned `k**2 > 1/(4*pi**2)` for a positive `k` into
  `(k > 1/(2*pi)) | (k < 1/(2*pi))`. That is a wrong set: it contains
  `0 < k < 1/(2*pi)`. Fixed by solving over a plain real dummy, which also
  fixes `fourier_transform(sin(x)/x, x, k)` for a positive `k` (the two COND
  entries) and `reduce_inequalities`.
- Laplace `t**a` rule: SymPy applies it only when `n > -1` is True. For an
  open condition it falls back to integration, giving
  `gamma(a + 1)/(s*s**a)` (and `exp(-t)*t**a` likewise).
- `meijerint_definite` split points: a bare `x` was treated as a possible
  boolean and skipped, so the split point 0 was missing. On `(1, oo)`,
  `1/(x*(1+x))` now gives `log(2)` rather than a meijerg, and `x**a` gives
  `-1/(a + 1)`.

Remaining (documented causes):
- MISSING (17): definite integrals with parameters, such as `1/(x+a)`,
  `1/(x-a)` and `1/(1+a*x)` on (0, 1) and `1/x` on (a, 1). These are
  intentional holds (see the triage above): SymPy's unconditional value is
  wrong for some parameter values (a pole in the interval).
- MISSING (1): `cos(s*x)` on (0, oo) with a nonzero `s`. SymPy returns
  `sin(oo*sign(s))/s` from its Limit leading-term machinery, and symox
  leaves `limit(sin(s*x), x, oo)` unevaluated. This is a limit gap, not a
  wrong value.
- DIFF (10): `laplace_transform(1/(t+1))`. SymPy's `-exp(-s)*Ei(-s)` is
  wrong: at s = 1 the true value is e*E1(1) = 0.5963, which equals
  symox's `-exp(s)*Ei(-s)`. Intentional correction.
- EQUIV (21), equivalent forms only:
  - `mellin_transform(besselj(0, x))`: `2**s/2` versus `2**(s - 1)`. Native
    simplify candidate selection; simplify area.
  - `cos(x)/x` on (1, oo): a meijerg-reduction form.
  - `hankel_transform(exp(-a*x)/x)` with real symbols: `Abs(k)/(k*sqrt(..))`
    versus `1/(k*sqrt(a**2/k**2 + 1))`.
  - `Heaviside(t-a)*sin(t)`: `cos(a - pi/2)` versus the exp form.
  - The AccumBounds entries `sin(s*x)` and `cos(s*x)`: they print
    identically, so this is a comparison artifact.

Tests: `tests/test_icorpus_parity.py`.

## mat-rest4 (wave 10)

Implemented:
- Derivatives of matrix expressions by a `MatrixSymbol` now follow
  SymPy's matrix-line algorithm (`_eval_derivative_matrix_lines`,
  `_LeftRightArgs`, `_matrix_derivative`), in
  `python/symox/_matrix_derivatives.py`. It covers MatrixSymbol, MatAdd,
  MatMul, Transpose, Inverse, MatPow (integer exponents and 1x1 bases),
  Trace, scalar Add/Mul/Pow, and a scalar free of the variable (a
  ZeroMatrix). Rank-2 results are contracted as in SymPy, and rank-4
  results are `PermuteDims(ArrayTensorProduct(...))` / `ArrayAdd`, as in
  SymPy. Examples: `Trace(X*X).diff(X) == 2*X.T` and
  `(x*x.T*x).diff(x) == KroneckerProduct(I, x.T*x) + 2*x*x.T`. A vector
  derived by a vector uses the products of SymPy's array conversion
  (SymPy's old line contraction gives a wrong `x**2` term here). Any node
  the port does not cover falls back to the previous array-of-entries
  path, and never to a guessed matrix form. Probe
  `matrices/expressions/tests/test_derivatives.py`: same 13 -> 67 of 107.
- A native scalar node (Trace, Determinant, ...) used as a coefficient
  prints under its head (`Trace(A)*I`, not `Expr(A)*I`).
- `Eq(M, N).lhs` / `.rhs` of explicit matrices are `ImmutableDenseMatrix`.
- The entrywise public functions (`simplify`, `expand`, `factor`,
  `cancel`, `together`, `trigsimp`, `powsimp`, `N`, ...) applied to a dense
  or sparse matrix return the immutable class, as SymPy's sympify does
  (`simplify(SparseMatrix(...))` is an `ImmutableSparseMatrix`).

Remaining (documented gaps):
- Matrix derivatives that are still unsupported, and fall back to the
  array path or raise: HadamardProduct/HadamardPower (SymPy gives
  `DiagMatrix(...)`), ElementwiseApplyFunction, MatrixElement derivatives
  (KroneckerDelta), `X.diff(A)` for a different symbol (ArrayDerivative),
  `Trace(A)*A` (a rank-4 result with a scalar factor), and expressions
  with explicit matrices of symbolic size.
- Form differences in test_derivatives: a 1x1 factor in a product of
  powers is placed differently (`b*(c.T*a*x.T*b)**(-1/2)*c.T*a` against
  SymPy's `b*a.T*c*(...)**(-1/2)`; both are equal because the factor is
  1x1). The large `test_derivatives_of_complicated_matrix_expr` sum also
  orders its terms differently.
- Items carried over from mat-rest3 and not done here: the adjoint of an
  `Operator` entry (`Dagger(a)`), the adjoint/transpose of `Abs(A)` and
  `sqrt(A)` for a noncommutative `A`, the `DFT`/`IDFT`/`GenericIdentity`/
  `GenericZeroMatrix` nodes, the Poly-entry deprecation warning and
  `test_smith_normal`, and the `type()` of BlockDiagMatrix/MatMul/... (still
  the native `_MatrixExpr`).

Tests: `tests/test_mat_rest4_parity.py`. Older tests that pinned the
previous symox behaviour were updated to SymPy's results:
`Trace(A).diff(A)` is `Identity(2)` (not an array), a symbolic-size
`A.diff(A)` is a `PermuteDims` (it used to raise), and the entrywise
functions return immutable matrices (test_array_public_api,
test_differentiation_options, test_matrix_scalar_codegen,
test_dirac_expansion). The finite jacobian test (each entry of a rank-4
result against the derivative of the explicit matrix) checks the axis
order of `A.T.diff(A)` and the other rank-4 results.

## simp-rest4 (wave 10)

Implemented:
- Issue 11004, high precision Floats through the polys layer: the RR domain
  now carries its binary precision (the largest precision of the input
  coefficients, as SymPy's RealField does). cancel, together, factor,
  factor_list, sqf, gcd, lcm, div, quo, rem, apart, ratsimp, terms_gcd,
  simplify and Poly (also with `domain=RealField(136)`) keep the 103 and 136
  bit Floats SymPy returns instead of rounding to 53 bits. The
  double-precision fast paths only run for 53-bit coefficients. Poly prints
  an RR zero as 0.0 and `all_coeffs` returns RR zeros.
- `Add._eval_is_zero` ported (extended real, imaginary and `coeff*I` terms):
  `(r + I).is_zero` and `(1 + I*r).is_zero` are False, so the zero and
  finiteness queries of powers with complex bases follow. On a 252-line probe
  (18 complex bases, 14 exponents, is_zero/is_finite/is_nonzero/is_infinite)
  the differing lines fell from 105 to 13. SymPy 1.14 gives
  `((1 + I)**x).is_zero` None (x may be infinite), and symox matches that.

  The port runs after the sign and numeric paths, and only for a sum with a
  term not known to be real, so shared real DAGs stay linear.

Simplify probe (sympy/simplify/tests): same 1719 -> 1720, diff 128 -> 127,
no statement from same to different.

Intentional divergence:
- `(b**I).is_nonzero` for a nonreal complex base b (1 + I, 2*I, E + I, ...):
  symox gives False (the value is nonreal, and nonzero implies real), SymPy
  None.

Remaining (unsupported or a different form, never a wrong value):
- `I**I` and `(-I)**I` stay unevaluated (SymPy gives exp(-pi/2) and
  exp(pi/2), so is_nonzero is True there); `(exp(I)**ip).is_nonzero` is None
  instead of False.
- `simplify((x + 5)/(16*x**2 + 28*x + 6))`: the native size-based selection
  keeps the cancelled form on a count_ops tie, where SymPy's `shorter` picks
  the `together` form `(x + 5)/(2*(8*x**2 + 14*x + 3))` (strict xfail in
  `tests/test_simp_rest3_parity.py`). The same tie rule is behind issue
  17141 (`simplify(acos(I)**2)`).
- `log(exp(x), 2, evaluate=False)`: re-evaluation in subs/doit/expand and
  `simplify(..., inverse=True)` (SymPy: `x/log(2)`).
- Derivative with a symbolic integer order (issue 22210) and
  `MatrixExpr.doit(inv_expand=...)`.
- `srepr(Poly(...))` prints the expression form instead of SymPy's
  `Poly(expr, x)` (pre-existing).

Tests: `tests/test_simp_rest4_parity.py`.

## mat-rest5 (wave 11)

Implemented:

- A dense matrix holding quantum operators (`Matrix([[Operator('a'), 1]])`)
  stores each operator atom as a noncommutative native node that prints
  as the atom. `M.adjoint()`, `M.H` and `Dagger(M)` give
  `Matrix([[Dagger(a)], [1]])`, `M.T` and `M.conjugate()` match SymPy, and
  the entries compare equal to the operator expressions.
- `adjoint`/`transpose` of a commutative node over a noncommutative
  symbol follow `Expr._eval_adjoint/_eval_transpose`: `adjoint(Abs(A))`
  and `transpose(Abs(A))` are `Abs(A)`, `adjoint(sign(A))` is
  `sign(conjugate(A))`. Powers follow `Pow._eval_adjoint/_eval_transpose`
  (integer exponent, positive or complex base, then `expand_complex`).
- Matrix expression classes are real Python subclasses of the native
  node: `type(X*Y)` is `symox.matrices.expressions.matmul.MatMul`, and the
  same holds for MatAdd, MatPow, Transpose, Inverse, Adjoint,
  BlockMatrix, BlockDiagMatrix, Identity, ZeroMatrix, OneMatrix, the
  Hadamard, Kronecker and diagonal classes, MatrixSlice,
  PermutationMatrix and MatrixPermute (`type(e) is e.func`, pickling
  roundtrips). The native constructor allocates the registered class
  directly, so no Python wrapper is created per node.
- A Poly entry in a Matrix emits SymPy's `deprecated-non-expr-in-matrix`
  warning, and `RawMatrix` emits `deprecated-rawmatrix`. SymPy's
  `test_smith_normal` passes; normal forms accept a `RawMatrix`.
- `(Trace(A)*M).diff(A)` gives SymPy's `ArrayTensorProduct(I, M)` terms
  (also with a transpose, a second trace factor or a coefficient).

Probe on SymPy's matrices and matrices/expressions tests (statements
compared with SymPy 1.14): same 4266 -> 5320, diff 324 -> 322, exc
161 -> 157, not reached 1381 -> 333. No statement went from same to
different.

Intentional divergences (equivalent values):

- `adjoint(sqrt(A))` gives the same expansion as SymPy, but the modulus
  `(re(A)**2 + im(A)**2)**(1/4)` is not wrapped in `Dagger`, because the
  native `re(A)` is commutative.
- `type(adjoint(a))` is `DaggerAtom`, not `Dagger`, and
  `type(Matrix([[a]])[0])` is the native operator node, not `Operator`.
- A Poly entry is stored as its expression (the native matrix holds
  `Expr` entries): `Matrix([[Poly(x)]])[0]` is `x`, not a `Poly`. The
  warning is emitted once per Poly entry (SymPy: once per matrix).
- `(Trace(A)**2*X).diff(A)` is one term `(2*Trace(A))*I` where SymPy has
  two equal terms; `Trace(A.T*X)` is not rewritten to `Trace(A*X.T)`.

Remaining gaps:

- DFT, IDFT, GenericIdentity and GenericZeroMatrix are not implemented.
  Every matrix expression class is a native node kind, so they need new
  native kinds (DFT entries `exp(-2*I*pi*i*j/n)/sqrt(n)`, and a shapeless
  node for the generic identity and zero). `MatMul()` gives `1` where
  SymPy gives `GenericIdentity()`.
- `test_smith_normal_deprecated`: `Matrix == RawMatrix` is False because
  the native comparison does not defer to the other operand.
- Matrix derivatives that SymPy builds through `array_derive` and
  `ArrayDiagonal`: Hadamard products and powers (with the DiagMatrix
  conversion), element-wise functions, `x.T*HadamardProduct(X, A)*x`.
  They need the array-to-matrix conversion of diagonals.
- `X.diff(A)` with A not in X: SymPy keeps `Derivative(X, A)`, but held
  matrix derivatives are not implemented natively (TypeError).
- Derivatives by a matrix element: `X[i, j].diff(X[m, n])` needs
  `KroneckerDelta(i, m, (0, k - 1))` (the native KroneckerDelta has no
  range argument), and `(X*A)[i, j]` stays a matrix element instead of
  SymPy's `Sum`.

## limacc (wave 11): AccumBounds limits, symbolic oscillation, branch cuts

Implemented:
- Bounded affine limits at infinity: `x*Ci(x)`, `x*(Si(x) - pi/2)`,
  `x*(fresnels(x) - 1/2)` and their sums, multiples and affine images give
  SymPy's `AccumBounds` (`AccumBounds(-1, 1)`, `AccumBounds(-1, 1)/pi`, ...).
- `sin(s*x)` at `oo` for a nonzero real `s` gives `sin(oo*sign(s))` (and the
  held forms built from it), so `integrate(cos(s*x), (x, 0, oo))` is
  `sin(oo*sign(s))/s` like SymPy. A possibly zero `s` stays a `Limit`.
- Complex infinite limits: `limitinf` returns `c0*oo` for an imaginary
  `c0 = I*a` and `oo*sign(c0)` for another nonreal `c0`
  (`limit((1 + I)*x, x, oo)` is `oo*sign(1 + I)`).
- `cos(m*x)/x` with a complex `m` is no longer 0: only a real infinity keeps
  `sin` and `cos` bounded.
- The bounded-value heuristic refuses a base or log argument made of the
  bounded values alone, so `sin(x)**(1/x)` and `(cos(x) + 1)**(1/x)` at `oo`
  are no longer the wrong value 1.
- `limit()` returns `Limit` for an expression holding `oo`, `-oo`, `zoo` or
  `nan` (as `Limit.doit`): `limit(oo/x, x, oo)` was 0, `(1 + x)**oo` was 1.
- One-sided limits onto a branch cut (`calculus/_branch_limit.py`): log,
  non-integer powers, asin, acos, asec, acsc, atan and acot whose argument
  reaches the cut from the non-principal side are rewritten (log(a) as
  log(-a) + s*I*pi, a**r with the matching phase, inverse trigonometric
  nodes by their one-sided value), each replacement checked numerically on
  the approach path. At `+-oo`, a log or power whose argument approaches the
  negative real axis from below is rewritten the same way
  (`limit(log(x - I) - log(-x - I), x, oo)` is `I*pi`, was `-I*pi`).
- The gruntz `sign` recursion is depth-bounded (the gamma expressions of
  test_gruntz lines 116, 130 and 139 segfaulted).
- Series: an exact Heaviside factor no longer makes a product inexact
  (`series(Heaviside(x)/x)` is `1/x`).
- `evalf(subs=...)` of an expression with a DiracDelta substitutes exactly
  and keeps `DiracDelta(x)` when the substituted node is a DiracDelta (SymPy
  returns self); `O(x*y).subs(y, x)` raises SymPy's ValueError about
  repeated variables.

Probes:
- liminf's 71-limit probe: 61 -> 62 same.
- sympy/series/tests (test_limits.py, test_gruntz.py, 648 asserted values):
  same 462 -> 493, different 186 -> 155, crashes 3 -> 0; no
  statement from same to different.
- Corpora: ODE unchanged (cls 301, hint 275, def 270); integration SAME
  1161 -> 1162, MISSING 18 -> 17, DIFF 10.

Intentional divergences (symox is right, SymPy is loose or fails):
- `limit(x*(Si(x) - pi/2) + cos(x), x, oo)` is 0 (SymPy `AccumBounds(-2, 2)`).
- `limit(x*Ci(x)**2, x, oo)` and `limit(x*Ci(x)/log(x), x, oo)` are 0
  (SymPy raises NotImplementedError).
- `gruntz(y*I*x, x, oo)` gives `oo*I*sign(y)` (SymPy `oo*I*y`).

Remaining:
- Wrong values still given: issue 22334 (`(n + 1)**k/((n + 1)**(k + 1) -
  n**(k + 1))` is 0, not `1/(k + 1)`); issue 6364 (the leading constant
  `1 - sin(a)**2 - cos(a)**2` is not recognized as zero; SymPy simplifies
  the leading term of a sum); issue 22220 at `-pi` and `pi+`; issues 6366
  and 6843 (nan, `oo*sign(0)`); issue 18399 (`oo` instead of `zoo`);
  `frac(x)**(1/x)` at `oo` is 1 (SymPy `AccumBounds(0, 1)`);
  `acoth(exp(1/x))` at `0-` is `I*pi/2` (SymPy `-I*pi/2`).
- Missing answers: Bessel functions at infinity, floor/ceiling, asech/acsch,
  LambertW, Mod, `asin(exp(x))`, `x*Shi(x)*I`, the non-affine
  `x**2*Ci(x)` and `Ci(x)*x*sin(x)`.
- Forms: Floats against exact values (issues 5164, 5383),
  `(2*DiracDelta(x)).evalf(subs={x: 0})` keeps `2*DiracDelta(0)` (SymPy
  `2.0*DiracDelta(0)`).
- Piecewise leftovers: `solveset(exp(sin(x)) > 1, x, S.Reals)` needs the
  inversion in solve.rs (left untouched for the pending user work);
  `series(sin(Piecewise(...)))` and `series(exp(Piecewise(...)))` raise
  NotImplementedError (SymPy's sin result is itself wrong, its exp result is
  right).

Tests: `tests/test_limacc_parity.py`.

## intrest3 (wave 11): integration leftovers

A probe of 460 parametric definite integrals (about 95 integrands under the
none/real/pos/nz/cplx assumption sets, values checked against mpmath
quadrature at real samples including negatives and 0) drove this pass.

Wrong results fixed:
- Powers and products of sin/cos with a symbolic frequency
  (`sin(a*x)**2`, `sin(a*x)*cos(a*x)`, `cos(a*x)**2`, `sin(a*x)**3`,
  `sin(a*x)**4`, `tan(a*x)**2`, `sec(a*x)**2`, ...) dropped the
  `Ne(a, 0)` guard of `trigintegrate`, so `integrate(sin(n*x)**2,
  (x, 0, pi))` was `pi/2` at n = 0. The native rules now integrate in
  u = a*x + b and return `Piecewise((F(u)/a, Ne(a, 0)), (f(b)*x, True))`,
  SymPy's form for these cases, and the definite integral is
  `Piecewise((pi/2, Ne(n, 0)), (0, True))`.
- The fallback of `deltaintegrate` for a held simple delta returned
  `Heaviside(a*x - b)` (SymPy's code, reached in symox because the
  expansion is kept while a may vanish): now
  `Piecewise((Heaviside(a*x - b)/a, Ne(a, 0)), (x*DiracDelta(b), True))`.
- Improper rational integrals whose poles cannot be located on the interval
  (`1/x**2` on (a, oo) with a real or nonzero, `1/(x + a)**2` on (0, oo),
  `1/(a + x**2)` on (-oo, oo) and on (0, oo) for a not known positive):
  the Meijer G method sees no pole and returns an unconditional value
  (`1/a`, `pi*a/(2*polar_lift(a)**(3/2))`, ...). SymPy returns it; it is
  wrong for a < 0 (at a = -1, `1/(x**2 - 1)` has a non-integrable pole at
  x = 1 and SymPy's value is the finite `-I*pi/2`). An unconditional Meijer
  value is now held when `rational_definite::analyze` cannot decide the
  poles; a conditioned one keeps its conditions. Intentional holds: the
  integration corpus goes from SAME 1161 to 1157 (the four
  `1/(x**2+a) (x,0,oo)` none/real/cplx/nz entries become MISSING).

Values SymPy gives that symox now gives too:
- `x/(x**2 + a**2)` and `1/(x**2 + a**2)` on (0, 1) with a real: the roots
  `I*a`, `-I*a` are real only at a = 0, where they sit on the interval. The
  analysis returns `Degenerate`, and the generic value is kept only when it
  is not finite at a = 0 (it is `-log(a**2)/2 + log(a**2 + 1)/2` and
  `atan(1/a)/a`, both infinite there); with 0 outside the interval the
  value is plain.

Intentional corrections (SymPy wrong, checked numerically):
- `integrate(DiracDelta(a*x), (x, -1, 1))`: SymPy 0, symox
  `Piecewise((-Heaviside(-a)/a + Heaviside(a)/a, Ne(a, 0)),
  (2*DiracDelta(0), True))` (1/|a| for a != 0).
- `integrate(DiracDelta(x - a, 1), (x, 0, 1))`: SymPy gives
  `DiracDelta(a) - DiracDelta(1 - a)`, the opposite sign; symox holds.

Remaining (unsupported or a different form, never a wrong value):
- `atan(a*x)/x**2` on (1, oo) for a real: needs
  `limit(a*(log(x) - log(a**2*x**2 + 1)/2) - atan(a*x)/x, x, oo)` (limits,
  task limacc).
- `DiracDelta(a*x - b)` on (-oo, oo): SymPy expands to
  `Heaviside(x - b/a)/Abs(a)` and gives `1/Abs(a)`; symox keeps the delta
  while a may vanish, and for a nonzero `limit(Heaviside(a*x - b)/a, x, oo)`
  is not evaluated. The native affine delta antiderivative is
  `Heaviside(a*x + b)/a`, not SymPy's `Heaviside(x + b/a)/Abs(a)` (same
  value for real nonzero a).
- `Heaviside(x - a)` on (-1, 1) with a positive: SymPy gives a Meijer G
  value; symox returns the integral of the rewritten Piecewise.
- Forms: `sin(a*x)**2*cos(a*x)**2` and `sin(2*a*x + 1)**2` (SymPy uses
  other reduction formulas), `log(1 + a*x)` on (0, 1) (Ne guard where
  SymPy's form has none), `(exp(-a*x) - exp(-b*x))/x` positive
  (`-log(-a) + log(-b)` vs `-log(a) + log(b)`).
- `ratint(1/(x**3 + x + 1))`: 82 KB against SymPy's 17 KB (the
  log_to_real/atan path; 0.4 s against 7.6 s).
- Risch/heurisch item: the reviewed cases (`atan(exp(x))`,
  `exp(x)/(1 + exp(2*x))`, `1/(x*(1 + log(x)**2))`, risch=True partial
  fractions, `x**(4/3)*log(x)`, issue 10847) already agree with SymPy.

Probe counts (460 cases): before {SAME 391, EQUIV 18, MISSING 31, COND 6,
DIFF 2, EXTRA-WRONG 4, COMPARE-FAIL 8}, after {SAME 391, EQUIV 18,
MISSING 35, COND 3, DIFF 2, EXTRA-WRONG 4, COMPARE-FAIL 7}. The four
EXTRA-WRONG are `Heaviside(x - a)*Heaviside(b - x)` on (-oo, oo), where
symox's `(b - a)*Heaviside(b - a)` is right and SymPy holds (a quadrature
artifact of the narrow support); the two DIFF are the DiracDelta(a*x)
correction above. The six SAME lost are the holds above, the six gained
are the guarded `sin(n*x)**2` and the degenerate-pole cases.

Tests: `tests/test_intrest3_parity.py`.

## simp-rest5 (wave 11)

Simplify selection (measured, not kept): returning SymPy's own sequence
(`reference_sequence`, the chain of `shorter` calls whose ties keep the
earlier form) whenever it applies, skipping the later native tie rules and
powsimp/exptrigsimp passes for it, and making the native `count_ops` rebuild
the numerator and denominator the way `fraction` does (a number times a
single sum is distributed, so 2*(8*x**2 + 14*x + 3) counts as
16*x**2 + 28*x + 6) gives SymPy's forms for
`simplify((x + 5)/(16*x**2 + 28*x + 6))` and issue 17141. It raised the
simplify probe from same 1720 to 1730 with no regression, the ODE corpus
from cls/hint/def 301/275/270 to 303/279/274 and the integration corpus
from SAME 1161 to 1162. It also changed 7 tests of the suite: the form of
the Poisson process probability
`(9*t1**3 + 9*t1**2 + 6*t1 + 2)*exp(-3*t1)/2`, the period of
`(sin(x) + cos(x))**2` (2*pi instead of pi), one ODE symmetry search, and 4
nonlinsolve systems that became NotImplementedError. The wave budget left
no time to track those callers down, so the native selection is unchanged
and the two cases stay strict xfails.

Assumptions (Pow._eval_is_extended_real / _eval_is_extended_positive):
- An imaginary base to an imaginary power is real and nonzero when log(base)
  is imaginary: `I**I`, `(-I)**I`, `I**(2*I)`, `I**i`, `exp(I)**i` with i
  imaginary. They are positive except `I**I` itself: SymPy's positive test
  starts with a `base == exp` branch, so `(I**I).is_positive` is None while
  `is_real` is True; symox follows that.
- A product of an even number of imaginary factors with known-sign real
  factors is real and nonzero (`I*i`, `2*I*i`).

Held two-argument log (`log(exp(x), 2, evaluate=False)`): it keeps its two
arguments in str, args and srepr, and re-evaluates to `log(exp(x))/log(2)` in
subs, doit, expand, `func(*args)` and simplify (log._eval_simplify), so
`simplify(..., inverse=True)` gives `x/log(2)`.

`MatrixExpr.doit(inv_expand=False)` keeps `Inverse(A*B)` unexpanded
(Inverse.doit returns self); other doit hints still raise
NotImplementedError.

Simplify probe (sympy/simplify/tests) with the kept changes: same 1720 -> 1721
(test_simplify_function_inverse, the held log), diff 127 -> 126, no
statement from same to different. ODE corpus (cls/hint/def 301/275/270 of
315) and integration corpus (SAME 1161, EQUIV 21, DIFF 10, MISSING 18) are
unchanged.

Remaining (unsupported or a different form, never a wrong value):
- The simplify selection above: `simplify((x + 5)/(16*x**2 + 28*x + 6))`
  keeps the cancelled form and issue 17141 keeps the expanded square.
- Derivative with a symbolic integer order (issue 22210): the order of a
  native Derivative is a machine integer, so `Derivative(f(x), (x, n))` with
  a symbolic n raises ValueError where SymPy keeps it unevaluated.
- `arg(I**I)` is printed in a different form.
- `(2*I)**i` facts: symox gives None where SymPy proves positive/imaginary
  False.

Tests: `tests/test_simp_rest5_parity.py`.

## audit12 (wave 12): c lots c1, c2 and c4

Replay of the 113 objects of the c lots c1, c2 and c4 against SymPy 1.14
(`audits/sympy-1.14.0/triage-wave12.md`): c1 27/36 -> 27/36, c2 30/40 ->
31/40, c4 37/37 -> 37/37. The gaps listed for these lots (`to_meijerg`,
symbolic `fps`, stochastic processes, multiple derivatives, `nseries` of
`x**y`, `trigsimp(old=True)`, `devise_plan_meijer`, `solve_triangulated`
with extension, c4 pretty printing) were already closed by earlier waves.

Implemented (Rust, `expand`):

- The `power_base` hint is the reference `Pow._eval_expand_power_base` for
  a commutative product base: an integer exponent distributes over every
  factor, otherwise the nonnegative and polar factors get their own power
  and the rest keeps one (with the reference handling of negative factors).
  Before, a base was split only when every factor was positive, so
  `expand((2*y)**z)` stayed `(2*y)**z` (SymPy `2**z*y**z`).
- The hint runs after `mul` and `multinomial`, as in the reference sorted
  hint order, so `expand((x*(y + z))**x)` with positive symbols is
  `(x*y + x*z)**x`. A second arithmetic pass runs when it changed something.
- The factors of a product are visited as power nodes: `expand(a*(a*x)**n)`
  with a positive `x` is `a*a**n*x**n`.

Remaining (form only, no wrong value):

- `Derivative(..., evaluate=True)`, `diff` and `doit` do not sort the
  variables with `_sort_variable_count`: `diff(f(x, y), y, x)` is
  `Derivative(f(x, y), y, x)` (SymPy `Derivative(f(x, y), x, y)`).
- The chain rule through an undefined function keeps the `Subs` form
  (`Derivative(g(x), x)*Subs(Derivative(f(_xi_1), _xi_1), _xi_1, g(x))`);
  SymPy 1.14 writes `Derivative(f(g(x)), g(x))*Derivative(g(x), x)`.
- `fps(exp(a*x), x)` with a symbolic coefficient returns the correct series;
  SymPy returns `exp(a*x)` unchanged, because its `_transform_explike_DE`
  passes `a` to `linsolve` as an unknown (NonlinearError, caught).
- The other differences of the replay are order, set-order or already
  documented items (see the triage file).

Tests: `tests/test_audit12_parity.py`.

## mat-rest6 (wave 12)

Native node kinds:
- `DFT(n)` and `IDFT(n)` (`sympy.matrices.expressions.fourier`) are native
  nodes reporting their public class. Entries are `exp(-2*pi*I/n)**(i*j)/sqrt(n)`
  (IDFT with the opposite sign), `DFT(n)*IDFT(n)` and `IDFT(n)*DFT(n)` cancel
  to `Identity(n)` (combine_powers), `.I` swaps the two kinds, `conjugate`
  gives `Adjoint(Transpose(DFT(n)))` (MatrixExpr._eval_conjugate), `Q.unitary`
  holds, subs/xreplace/has see the dimension, latex is `\text{DFT}_{n}`.
- `GenericIdentity()` and `GenericZeroMatrix()`
  (`sympy.matrices.expressions.special`) are shapeless native nodes:
  `shape`, `rows`, `cols`, entries and shape-dependent operations raise
  TypeError "... does not have a specified shape". `MatMul()` and `MatAdd()`
  return them, `MatMul.identity`/`MatAdd.identity` are them, and
  MatMul/MatAdd (constructor and native product/sum) drop them first, so
  `MatMul(G, G)` is 1 and `MatAdd(Z, Z)` raises ValueError as in SymPy.
  `G*x` and `x*G` are `x` (MatMul(G, x).doit()): `2*G` is 2, `G*A` is A,
  as nc_simplify needs. They are instances of Identity/ZeroMatrix. Neither class is exported at
  top level (SymPy does not export them either).

`Matrix == other`: when other is not a matrix the native comparison now
returns NotImplemented (RepMatrix.__eq__), so the other operand's `__eq__`
runs; `smith_normal_form(RawMatrix) == RawMatrix(...)` holds
(test_smith_normal_deprecated). Comparisons with scalars, lists and tuples
still give False.

Held matrix derivatives: `X.diff(A)` with matrix variables that X does not
contain gives SymPy's `ArrayDerivative` (str `Derivative(X, A)`, shape
`A.shape*count + X.shape`, expr, variables, variable_count, latex), and
`Derivative(X, A)` gives the held Derivative whose doit (and
`evaluate=True`) is `ZeroMatrix(*X.shape)`. Both are instances of
Derivative. Derivative(X, A) with A in X still raises NotImplementedError
(SymPy's own doit of it is wrong or raises TypeError).

Matrix probe (sympy/matrices/expressions/tests): same 5305 -> 5331, exc
154 -> 139, notreached 351 -> 337, diff 322 -> 325 (statements now reached),
no statement from same to different.

Intentional divergences (error type or message only):
- `DFT(n).diff(x)` raises NotImplementedError where SymPy raises
  AttributeError ('super' object has no attribute '_eval_derivative').
- `DFT(-1)` raises ValueError with symox's dimension message.
- `MatMul(GenericIdentity(), GenericZeroMatrix())` prints 0 where SymPy's
  printer raises.
- The srepr of a MatrixSymbol inside ArrayDerivative prints its name as a
  string literal (symox-wide srepr of MatrixSymbol).

Remaining (unsupported, never a wrong value):
- Derivative with a symbolic integer order (issue 22210,
  `Derivative(f(x), (x, n))`): not contained. The native node is
  `Kind::Derivative(Expr, Vec<(Expr, u32)>)`; storing the order as an Expr
  touches 90 match sites in 40 files (diff, subs, the printers, ODE, Laplace
  transforms, lambdify, evalf, the matrix DAG serialization, and solve.rs,
  which this wave must not change) plus the order arithmetic that assumes a
  machine count, so symox still raises ValueError where SymPy keeps it held.
- Matrix derivatives by a matrix that need the array algebra: Hadamard
  products/powers (`HadamardProduct(X, A).diff(X)`, `HadamardPower(X, 2)`),
  element-wise functions (`X.applyfunc(sin).diff(X)`) give ArrayDiagonal /
  ArrayAdd forms in SymPy; symox raises TypeError (symbolic shapes in the
  array path). `(x.T*HadamardProduct(X, A)*x).diff(X)` (SymPy:
  `DiagMatrix(x)*A*DiagMatrix(x)`) raises NotImplementedError.
- Matrix element calculus: `X[i, j].diff(X[m, n])` stays a held Derivative
  (SymPy: KroneckerDelta with a range), `(X*A)[i, j]` stays a MatrixElement
  (SymPy: `Sum(X[i, _i_1]*A[_i_1, j], (_i_1, 0, k - 1))`); symox
  KroneckerDelta has no delta_range.
- `Trace(DFT(n))` and `Determinant(DFT(n))` report type Expr instead of
  Trace/Determinant (pre-existing for all matrix expressions).

Tests: `tests/test_mat_rest6_parity.py`.

## limwrong (wave 12)

Wrong limit values after limacc, by root cause:
- Hidden zero leading coefficients (issues 6364, 6366/6843, 21598): the
  series engine kept a leading coefficient such as
  `1 - sin(a)**2 - cos(a)**2` or `log(1/b) + log(b)` (b > 0) as nonzero and
  divided by it. A sum of two or more leading terms whose symbolic
  coefficient vanishes numerically at two sample points is now proven zero
  with expand, the assumption-respecting expand_log, trigsimp and cancel
  (SymPy's Add._eval_as_leading_term uses trigsimp().cancel()); a zero that
  cannot be proven makes the expansion unsupported instead of wrong. The
  Gruntz symbolic-exponent fallback does the same for its sums.
- log branch cut (issue 22220 at -pi and pi+): log_ser now ports STEP 5 of
  log._eval_nseries: a negative leading coefficient followed by a nonreal
  term adds `-2*I*pi` when that term's imaginary part is negative.
- `base -> -oo` with `exponent -> oo` gives zoo (Limit.pow_heuristics,
  issue 18399) instead of oo.
- `AccumBounds**t` with `t -> 0` keeps zero and infinite endpoints
  (`frac(x)**(1/x)` at oo is `AccumBounds(0, 1)`, `(tan(x)**2)**(2/x)` is
  `AccumBounds(0, oo)`) instead of collapsing to 1.
- acoth at a branch point (issue 25847): the Gruntz leading term of
  `acoth(c*w**e)` follows acoth._eval_as_leading_term (`1/arg` for an
  infinite argument, `-sign(c)*I*pi/2` for a vanishing real one) instead of
  evaluating acoth(0) on the cut.
- Issue 22334: symbolic powers of the mrv variable (`(1/w)**k`) are now
  factored out through a positive placeholder when all exponents differ by
  numbers, and `exp(c*log(w) + r)` with a symbolic c folds when the known
  value of `log(w)` cancels exactly. The cancellation in the denominator is
  still not resolved, so the limit stays unevaluated (it was the wrong 0).

Probe (limacc's sympy/series/tests/test_limits.py + test_gruntz.py): same
493 -> 516 of 648, no statement from same to different.

Remaining (unsupported or a different form, never a wrong value):
- Issue 22334: `Limit` instead of `1/(k + 1)` (strict xfail).
- Issue 6364: `1/(-sin(a)**2/2 + cos(a)**2/2)`, equivalent to SymPy's
  `1/(cos(a)**2 - 1/2)`.
- Not addressed in this wave: Bessel functions at infinity, floor/ceiling,
  asech/acsch, LambertW, Mod, asin(exp(x)), `x*Shi(x)*I`, `x**2*Ci(x)`;
  intrest3's dependency (`integrate(atan(a*x)/x**2, (x, 1, oo))`,
  `integrate(DiracDelta(a*x - b), (x, -oo, oo))`); float-vs-exact forms of
  issues 5164 and 5383.

Tests: `tests/test_limwrong_parity.py`.

## simp-tie (wave 12)

simplify returns SymPy's own sequence (`reference_sequence`) whenever it
applies, and count_ops distributes a number over a single sum in a
numerator or denominator, as `fraction` does:
`simplify((x + 5)/(16*x**2 + 28*x + 6))` is `(x + 5)/(2*(8*x**2 + 14*x + 3))`
and issue 17141 keeps the content-primitive square.

Callers that relied on the earlier native forms:
- The period of `(sin(x) + cos(x))**2` is 2*pi, as SymPy (simplify gives
  `2*sin(x + pi/4)**2`).
- `P` of a bounded event of a discrete distribution sums the pmf directly,
  as SymPy's `SingleDiscretePSpace.eval_prob` / `eval_sum_direct`, without a
  final simplify: `P(Poisson('X', 3*t) < 4)` is
  `(9*t**3/2 + 9*t**2/2 + 3*t + 1)*exp(-3*t)`.
- checkinfsol is SymPy's single simplify of the determining residual, so a
  residual can stay unreduced where SymPy's checker leaves it unreduced.
- The native polynomial system solvers (`parametric::nonlinsolve`,
  `systems::solve_poly_system`) run in native-selection mode: SymPy's
  sequence is one more candidate of the native search, as before, since
  SymPy's solvers do not call simplify on their roots and restriction
  factors. Otherwise a quartic restriction factor of
  `nonlinsolve([x**3 - 1/y, x + y - 1])` is undecided by the exact zero test,
  and the x component of `nonlinsolve([x**3 + 2/y - 2, x - y])` becomes a
  form whose inner radicand is a negative real with a signed zero imaginary
  part (mpmath evaluates it on the other sqrt branch).
- The zero-dimensional fallback of nonlinsolve tests a denominator at a
  root on the substituted value before simplifying it.

Tests: `tests/test_simp_tie_parity.py`.

Simplify probe (sympy/simplify/tests): same 1721 -> 1730, diff 126 -> 117,
no statement from same to different. ODE corpus cls/hint/def
301/275/270 -> 303/279/274 of 315; integration corpus SAME 1158 -> 1159,
EQUIV 21 -> 20, DIFF 10, MISSING 21.

## intrest4 (wave 13): integration corpus leftovers

Wrong results fixed (intentional corrections, SymPy is wrong):
- The Meijer G rewrites of a shifted Heaviside (`Heaviside(x - b)`,
  `Heaviside(b - x)`, `Heaviside(2*x + 4)`) hold for `x > 0` only, but
  `meijerint_indefinite` returns them as antiderivatives on the whole real
  line. The reference's `integrate(Heaviside(x - 1/2), (x, -1, 1))` is 2
  (true value 1/2), `integrate(Heaviside(x - a), (x, -1, 1))` with a
  positive is 2 at a = 0.3 (true 0.7), `Heaviside(a - x)` and
  `Heaviside(x + a)` give 0 at a = 0.3 (true 1.3), and the antiderivative
  of `Heaviside(2*x + 4)` has derivative 0 for x > 2. In the term loop of
  `integrate`, a term with a Heaviside of `x` now skips the Meijer G stage
  unless `x` is nonnegative or the definite integral runs over nonnegative
  bounds; it goes to `manualintegrate`, as in the reference's next stage:
  `(x - 1/2)*Heaviside(x - 1/2)`, value 1/2; `(1 - a)*Heaviside(1 - a)` for
  a positive; `(1 - b)*Heaviside(1 - b) - (-b - 1)*Heaviside(-b - 1)`
  without assumptions, where the reference held an integral of the invalid
  rewrite. On (0, 1) the reference's G form is right and is kept (the
  corpus entries stay SAME). The public `meijerint_indefinite` is
  unchanged.

Forms ported:
- `_simpc` of `_laplace_trig_ltex`: the numerator coefficients of a pole
  pair are rewritten in `cos` through `as_real_imag`, so
  `laplace_transform(Heaviside(t - a)*sin(t), t, s)` with a positive is
  `(s*cos(a - pi/2) + cos(a))*exp(-a*s)/(s**2 + 1)` and
  `laplace_transform(sin(t + r), t, s)` with r real is
  `(s*cos(r - pi/2) + cos(r))/(s**2 + 1)` (both EQUIV -> SAME).

Corpus entries reviewed and kept:
- DIFF (10), `laplace_transform(1/(t+1))`: SymPy's `-exp(-s)*Ei(-s)` is
  wrong, symox's `-exp(s)*Ei(-s)` equals e**s*E1(s) and the quadrature at
  s = 1/2, 1, 3. Intentional correction.
- MISSING (17 of 21), `1/(x+a)`, `1/(x-a)`, `1/(1+a*x)` on (0, 1) and `1/x`
  on (a, 1): held on purpose, SymPy's unconditional value is that of a
  divergent integral when the pole reaches the interval (a = -1/2 for
  `1/(x+a)`, `1/(1+a*x)` with a = -2, `1/x` with a = -1, `1/(x-a)` with
  a = 1/2, also for the positive set). The other 4 are the documented
  `1/(x**2+a)` holds on (0, oo).
- MISSING (1) `cos(s*x)` on (0, oo) with s nonzero: depends on
  `limit(sin(s*x)/s, x, oo)` (limits, task limrest).

Remaining (different forms, never a wrong value):
- `mellin_transform(besselj(0, x))` (8 EQUIV): `2**s/2` vs `2**(s - 1)`:
  `simplify(2**(s - 1)*gamma(s/2)/gamma(1 - s/2))` picks the other form
  of a count_ops tie (simplify candidate selection, task simp-tie).
- `hankel_transform(exp(-a*x)/x, x, k, 0)` with real symbols: both libraries
  integrate to `1/(k*sqrt(a**2/k**2 + 1))`; SymPy's `simplify` gives
  `Abs(k)/(k*sqrt(a**2 + k**2))`, symox's keeps it (simplify area).
- `ratint(1/(x**3 + x + 1))`: 82 KB vs 17 KB. `log_to_real` is the same
  algorithm; the forms split at `roots_quadratic` of the `v` polynomial in
  the EX domain, which calls `simplify(-c/a)` on a nested cube-root radical
  and gets a different candidate (simplify area; `poly/roots.rs` also has
  pending user work), and grow from there in `A**2 + B**2`.
- `cos(x)/x` on (1, oo) (5 EQUIV): the hyperexpanded Meijer G terms are
  combined differently (`-2*EulerGamma/sqrt(pi)` vs
  `-(EulerGamma + 2*log(2))/sqrt(pi) - EulerGamma/sqrt(pi)`); not traced
  within the budget.
- `Heaviside(x - a)` on (-1, 1) with an unassumed or real a: symox now
  evaluates (above) where SymPy holds the invalid rewrite.

Tests: `tests/test_intrest4_parity.py`; `tests/test_integrate_methods_corpus.py`
records the corrected `integrate(Heaviside(2*x + 4), x)`.

## limrest (wave 13)

Implemented (ports of the SymPy mechanisms, no case patches):
- Bessel functions at real infinity: the Gruntz tractable rewrite uses the
  Hankel asymptotic expansions (besselj/bessely/besseli/besselk.
  _eval_aseries), with auxiliary asymptotic-series nodes for P, Q and the
  I/K tails. At -oo with an integer order, the reflection formulas
  `J(-z) = (-1)**n*J(z)`, `I(-z) = (-1)**n*I(z)`,
  `Y(-z) = (-1)**n*(Y(z) + 2*I*J(z))`, `K(-z) = (-1)**n*K(z) - I*pi*I(z)`
  are applied (issues 19154, 26040).
- Taylor series at zero of besselj/besseli with an integer order and of
  LambertW. A composed series whose low coefficients all vanish is now
  inexact and truncated, which removes the wrong 0 of issue 21701.
- asec/acsc/asech/acsch series go through acos/asin/acosh/asinh of 1/arg
  (issue 10868 and the asec/acsc/asech/acsch statements of 25582/25847). The
  constant term keeps the original function form at a regular point.
- asin/acos/asinh/acosh at an infinite argument: principal-branch log
  expansions (asin(exp(x)) at oo, issues 25582 and 25847).
- Trigonometric poles whose location carries pi: the constant is cancelled
  before the pole test (issue 14411).
- One-sided limits onto a sqrt/log branch cut when the imaginary part of the
  argument vanishes to first order: the side is taken from the leading term
  of im(arg(z0 + h)) in h (issue 23752).
- `limit(f, x, c*oo)` with a non-real direction substitutes `x -> c*x/|c|`
  and takes the limit at oo, like Limit.doit.
- Gruntz leadterm: a symbolic power `b**p` with b positive near 0 is split
  as `c**p * w**(m*p) * (b/(c*w**m))**p` before retrying the symbolic power
  leading term.

Probe (sympy/series/tests/test_limits.py + test_gruntz.py): same 516 -> 560
of 648, no statement from same to different.

Remaining (unsupported or a different form):
- Issue 22334: still `Limit` instead of `1/(k + 1)`.
- Bessel functions at `I*oo`/`-I*oo` (besselj, bessely, besselk) and
  `besselk(1, x)` at -oo stay unevaluated.
- floor/ceiling limits (200, 201, 236, 237, 958), Mod (25230),
  `acoth(sin(x)/x)` with '+-', Gruntz special-function PoleErrors (digamma,
  loggamma, gamma), `x*Shi(x)*I`, `x**2*Ci(x)`.
- intrest3 dependencies (`integrate(atan(a*x)/x**2, (x, 1, oo))`,
  `integrate(DiracDelta(a*x - b), (x, -oo, oo))`) not addressed.
- Float-vs-exact forms of issues 5164 and 5383 (`4.00000000000000` vs `4`).
- `limit(asec(x), x, 1/2)` gives `asec(1/2)` where SymPy gives
  `I*log(sqrt(3) + 2)` (unchanged from before, same value).
- Deliberate differences: `limit(besseli(1, x)*sqrt(x)/exp(x), x, oo)` is
  `sqrt(2)/(2*sqrt(pi))`; SymPy 1.14 returns `sqrt(2)*sqrt(pi)/2`, which is
  numerically wrong. `limit(besselj(0, x)*sqrt(x), x, oo)` gives
  `sqrt(2)*AccumBounds(-1, 1)/sqrt(pi)` where SymPy stays unevaluated.

Tests: `tests/test_limrest_parity.py`.

## mat-rest7 (wave 13)

Entries of symbolic matrix products and derivatives by a matrix element,
ported from SymPy's `_entry` and `_eval_derivative` code paths:
- `MatMul._entry` with a symbolic inner dimension gives the held
  `Sum(A[i, _i_1]*B[_i_1, j], (_i_1, 0, k - 1))`. Nested entries share one
  dummy counter, as SymPy's `dummy_generator` does. `MatPow._entry` with a
  positive integer exponent expands to the product, and `Identity(k)[i, j]`
  is `KroneckerDelta(i, j, (0, k - 1))`.
- The `MatrixElement(...)` constructor, `xreplace`/`subs` rebuilds and
  `srepr` parsing still keep symbolic products, powers and identities held,
  as SymPy does. Only `parent[i, j]` expands.
- Ranged `KroneckerDelta(i, j, (lo, hi))`: a Tuple argument and the
  `delta_range` property. A delta sum needs no Piecewise when the range is
  inside the summation limits. Sums of an Add and nested Sums with a delta
  are summed term by term or with the summations exchanged.
- `MatrixElement._eval_derivative`: `X[i, j].diff(X[m, n])` is a product of
  ranged deltas. For an `Inverse` parent it is
  `-Sum(M[i, z1]*Y[z1, z2].diff(v)*M[z2, j])`. `Sum.diff` by a matrix element
  differentiates under the sum. `Trace(A).diff(A[m, n])` goes through
  `rewrite(Sum)`. `Trace(X).rewrite(Sum)` works now (the rewrite visitor
  skipped non-Expr arguments).
- In array-to-matrix conversion, a contracted cycle's Trace is evaluated
  (`doit(deep=False)`), so `Trace(Trace(A)*A).diff(A)` is `2*Trace(A)*I`.

Probe (matrices tests, 41 files): same 5331 -> 5348, diff 325 -> 314,
exc 139 -> 133, no statement from same to different.

Remaining (a different form or unsupported, never a wrong value):
- Held matrix `Derivative(X, X)` and `Derivative(X*A, X)` raise
  NotImplementedError. SymPy keeps them held, but its
  `Derivative(X, X).doit()` gives 0, which is wrong.
- Derivatives of HadamardProduct, HadamardPower and applyfunc need
  `ArrayDiagonal` to matrix conversion (`DiagMatrix` forms). These still raise.
- The cookbook quadratic-form element derivative is equal, but SymPy pulls
  `b[n, 0]` out of the first Sum.
- `Trace(Trace(Trace(A)*A)*A).diff(A)` keeps `Trace(Trace(A)*A)`
  unevaluated, because SymPy's MatMul evaluates a scalar Trace coefficient
  when it is built. `(Trace(A)**2*X).diff(A)` is one term where SymPy has two.
- A product entry with mixed finite and symbolic inner dimensions stays a
  held Sum, without the partial expansion that SymPy's `doit` does.
- `MatrixElement(A*B, 0, 0).subs(n, 2)` evaluates to the explicit sum (as
  before). SymPy keeps it held until `doit`.

Tests: `tests/test_mat_rest7_parity.py`.

## limrest2 (wave 14)

Implemented (series engine and Gruntz):
- `polygamma(n, z)`, `gamma(z)` and `loggamma(z)` of an argument whose
  leading term is a constant in `logt` that tends to `+oo` (`digamma(x)`,
  `loggamma(x)`, `log(x)` after the mrv rewrite). They are expanded by Taylor
  around that constant: `polygamma(N, c0 + r) = Sum(polygamma(N + k, c0)*r**k/k!)`.
  This follows the reference `_eval_nseries` at a finite point. The sign of a
  coefficient in `logt` is decided by a Gruntz limit with `logt = -L`, `L -> oo`
  (nesting is limited to depth 4). The Stirling and A&S 6.3.18 branches and
  the "exponentially small" test in `exp` use the same sign.
- The Gruntz `sign` retries with the exponential factors of each product
  merged when the mrv rewrite was undone by the canonical product. This
  covers `exp(x)*(gamma(x + exp(-x)) - gamma(x))` and the gamma/Stirling
  differences.
- A leading coefficient `I*a`, with `a` real and depending on x, gives
  `I*sign(a)*oo`, like the reference `c0*oo` (`x*Shi(x)*I`, `I*x*exp(x)`).
  This also covers `besselj(1, x)` at `+-I*oo` and `besselk(1, x)` at `-oo`.

Probe (`limit`/`gruntz`/`Limit` statements of test_limits.py and
test_gruntz.py): same 560 -> 570 of 648, no statement from same to different.
ODE corpus and integration corpus are unchanged from their baselines.

Intentional divergence:
- `limit(x**2*Ci(x), x, oo)` stays unevaluated. SymPy returns
  `Ci(1/_w)/_w**2`, which is meaningless. The limit does not exist, because
  `x**2*Ci(x) ~ x*sin(x)`.

Remaining (unevaluated, never a wrong value):
- Issue 22334, `(n + 1)**k/((n + 1)**(k + 1) - n**(k + 1))`. After the mrv
  rewrite, the expression has `exp(-k*x + k*log(1 + 1/w))` with `log(w) = -x`.
  A symbolic `k` multiplies `x`, and the symbolic-power leading term does not
  handle that yet. It needs `w**k` factors kept through the expansion.
- `Ei(x - exp(-exp(x))) - Ei(x)` and the `erf` analogue: the leading terms
  cancel symbolically.
- `uppergamma(x, 1)/gamma(x)` and `subfactorial(x)` at oo: the mrv set of
  these functions is not known.
- `log(factorial(x))` Stirling remainder: hidden zero
  `-log(pi)/2 - log(2)/2 + log(2*pi)/2`.
- floor/ceiling/Mod limits at oo (issue 25230). These need an
  `AccumBounds` coefficient in the series engine (reference
  `floor._eval_nseries`). Canonical `AB - AB` cancellation would make that
  unsafe, so it is not done.
- `acoth(sin(x)/x)` with `dir='+-'`, `(-1)**x*2**(-x)` (issue 22181),
  `integrate(cos(s*x), (x, 0, oo))`, `integrate(atan(a*x)/x**2, (x, 1, oo))`,
  and the Float-vs-exact forms of issues 5164 and 5383.

Tests: `tests/test_limrest2_parity.py`. The strict xfail on `x*Shi(x)*I` in
`tests/test_limacc_parity.py` is removed, because that limit now evaluates.

## mat-rest8 (wave 14)

Done:
- Derivatives of HadamardProduct, HadamardPower and `applyfunc` by a
  MatrixSymbol go through SymPy's array-expression path, which is now
  vendored in `symox/_arrayexpr` (`convert_matrix_to_array`, `array_derive`,
  `convert_array_to_matrix`). `ArrayDiagonal` converts to `DiagMatrix`
  forms. `DiagMatrix.doit` pulls scalars out of a vector product and unwraps
  a Transpose, as SymPy does. A HadamardPower with a symbolic shape and a
  scalar exponent gives `hadamard_product(dexp*log.(base) + exp*dlbase, self)`.
- `hadamard_product` sorts its arguments in SymPy's canonical order:
  MatMul, then MatAdd, then the other matrix classes.
- A held matrix `Derivative(X, X)` or `Derivative(X*A, X)` is built and
  stays held, as in SymPy. Its `doit` returns `expr.diff(v)` for each
  variable, or stays held when that diff is unsupported. Here SymPy is
  wrong or raises: `Derivative(X, X).doit()` is 0 and
  `Derivative(X*A, X).doit()` raises TypeError. symox gives the array
  derivative `PermuteDims(ArrayTensorProduct(I, I), (3)(1 2))`, which is
  equal to `X.diff(X)` in both libraries.
- For a product entry with mixed finite and symbolic inner dimensions, the
  leading finite limits are summed (as SymPy's `eval_sum_direct` does). The
  rest stays a Sum over the remaining limits.
- `MatrixElement(A*B, 0, 0).subs(n, 2)` stays held, and `doit` expands it
  (`args[0][args[1], args[2]]`).

Probe (matrices tests, 41 files): same 5350 -> 5366, diff 312 -> 311,
exc 133 -> 126, no statement from same to different.

Remaining (a different form, never a wrong value):
- `hadamard_power(a.T*X*b, S.Half).diff(X)` is
  `a*Lambda(_d, 1/(2*sqrt(_d))).(a.T*X*b)*b.T`. SymPy gives
  `(1/2)*a*(a.T*X*b)**(-1/2)*b.T`, because its ElementwiseApplyFunction
  applies the function directly to a 1x1 argument. In symox a 1x1
  MatrixExpr inside a scalar Lambda is turned into a scalar.
- A nested `applyfunc` derivative is composed only through a `doit` retry.
  For the uncomposed form, SymPy's conversion raises IndexError.
- `Trace(Trace(Trace(A)*A)*A).diff(A)` and the two-term
  `(Trace(A)**2*X).diff(A)` split are still open (see mat-rest7). SymPy's
  `MatrixExpr.__mul__` does a deep `doit` of the product, and symox's `*`
  does not.
- The cookbook quadratic-form element derivative is equal, but SymPy pulls
  `b[n, 0]` out of the first Sum. This was not attempted.

Tests: `tests/test_mat_rest8_parity.py`. Three older tests that expected an
error or the eager entry now check the SymPy values.

## audit14 (wave 14): derivative forms and simp-tie leftovers

Fixed:
- `diff`, `Derivative(..., evaluate=True)` and `Derivative.doit()` now sort
  the differentiation variables as SymPy's `_sort_variable_count` does
  (`diff::sort_variable_count`: `_block` dependencies, ties ordered by
  `(_nodes, default_sort_key)`, then a merge of equal neighbours). A nested
  Derivative is combined as SymPy does. The result goes through the canonical
  rebuild of every inner Derivative.
- Chain rule through an applied function: when the argument is `g(x)` and no
  other argument contains it, the form is now `Derivative(f(g(x)), g(x))`,
  as in SymPy. Before, symox gave a `Subs` here.
- `Derivative._eval_subs` cases:
  - case (0): an invalid new variable that replaces a non-symbol variable
    gives `Subs(self.xreplace({old: xi}), xi, new)`;
  - case (1): a valid new variable that brings another differentiation
    variable into the expression gives `Subs`;
  - case (3): `subs(g, Lambda)` on a variable `g(x)` gives
    `Subs(Derivative(..., xi_i), xi_i, value)`.
- `Subs.doit` of a Derivative follows SymPy's steps:
  1. the substitutions unrelated to the differentiation are done first;
  2. the derivative is taken again with respect to the variables that are
     present;
  3. the remaining substitutions go through `Derivative._eval_subs`.

  The blanket hold of a Subs over a bound derivative is gone. So
  `Subs(Derivative(f(x, y), y), y, z).doit()` is now
  `Derivative(f(x, z), z)`, as in SymPy.
- `Subs._eval_derivative` is ported: a chain-rule term for each point, plus
  the direct term. A non-symbol variable goes through a dummy. Wrong values
  fixed:
  - `diff(Subs(Derivative(f(x, y), y), y, g(x)), g(x))` was 0;
  - `diff(f(x, y, g(x)), y, g(x), x)` missed a cross term;
  - `diff(f(x, g(x)), x, 2)` kept a held `Derivative(Subs(...), x)`.
- `Hyper_Function._is_suitable_origin` is now the full SymPy test (no
  `ai - bj` that is a non-negative integer, no zero `ai`, no non-positive
  integer `bj`); before it was a wrong True. `Hyper_Function.has` now sees
  the parameter Tuples.

Intentional correction:
- `diff(f(g(x), g(x)**2), x)`. SymPy writes the first slot as
  `Derivative(f(g(x), g(x)**2), g(x))`, which also differentiates the
  `g(x)**2` argument. With f -> x*y and g -> sin, SymPy then gives
  `5*sin(x)**2*cos(x)` instead of `3*sin(x)**2*cos(x)`.
  - symox keeps `Subs(Derivative(f(_xi_1, g(x)**2), _xi_1), _xi_1, g(x))`
    whenever another argument contains the variable structurally.
  - Its value agrees with `diff(sin(x)**3, x)`: simplify gives 0, and
    30-digit checks at x = 1/3 and x = 7/5 pass.

simp-tie probe (15 simplify test files): same 1747 -> 1747, diff 117 -> 117
(the Hyper_Function statements changed after the run), no statement from
same to different. One run had `test_prudnikov_1` time out under machine
load. It was rerun with a 300 s timeout, and its 17 statements are the same
as before.

Remaining (documented, not fixed):
- `mellin_transform(besselj(0, x), x, s)` gives `2**s/2` where SymPy gives
  `2**(s - 1)`. SymPy's `transforms._simplify` calls
  `powdenest(expr, polar=True)`, which turns `2*4**(s/2 - 1)` into
  `2**(s - 1)`. symox's `powdenest(polar=True)` raises NotImplementedError,
  so porting that option is the fix.
- `hankel_transform(exp(-a*x)/x, x, k, 0)` now matches SymPy.
- `integrate(cos(x)/x, (x, 1, oo))` is equal to SymPy's result (both are
  `-Ci(1)`), but the unsimplified sum after hyperexpand groups the
  EulerGamma and log terms differently. The responsible step is not traced.
- `ratint(1/(x**3 + x + 1))`: SymPy's result is 16653 characters long,
  symox's is 82704. The size comes from `roots_quadratic` in
  `poly/roots.rs`. That file is in the user's work area and was not touched.
- Dummy naming: `diff(f(g(x)), x).subs(g, Lambda(x, x**2)).doit()` prints
  `_xi_0` where SymPy prints `_xi`. The value is the same.
- simp-tie probe, still different:
  - `test_cse` (75), `test_fu` (13), `test_simplify` (10),
    `test_radsimp` (5), `test_powsimp` (4), `test_trigsimp` (1);
  - `test_hyperexpand`: 4 statements, including
    `test_meijerg_with_Floats`, which raises TypeError, and the
    `test_partial_simp` parameter choice;
  - `test_function`: `Hyper_Function.args` are Python tuples, where SymPy
    has Tuples.

  These were not worked through in this task.

Tests: `tests/test_audit14_parity.py`.

## pw-nc (wave 15): Piecewise and noncommutative leftovers

Implemented:
- The commutativity of a `Piecewise` follows its branch expressions, as
  `ExprCondPair.is_commutative` does in SymPy; conditions are ignored
  (`decompose.rs`). So `Piecewise((A, y < 1), (B, True)).is_commutative` is
  False, `C*p*A` keeps its factor order and
  `piecewise_fold(C*p*A)` is `Piecewise((C*A**2, y < 1), (C*B*A, True))`.
  Before, the Piecewise was commutative and was moved before `C`.
- Partial `evalf` follows `Mul._eval_evalf` and `AssocOp._eval_evalf`:
  - a -1 coefficient is kept: `(-x).evalf()` is `-x`, not `-1.0*x`;
  - `I` is not evaluated inside an Add, a Mul or a Pow: `(x + I).evalf()` is
    `x + I`, and `sin(I*x).evalf()` is `I*sinh(x)`.
  Before, all three were wrong forms with a `1.0` coefficient.
- `adjoint` printing matches the reference. `latex(adjoint(A))` is
  `A^{\dagger}` (no parentheses), and a power is
  `\left(A^{\dagger}\right)^{2}`. `pretty(adjoint(A))` is the base with a
  `†` exponent (`+` in ASCII).
- Matrices with noncommutative entries (`python/symox/_matrix_nc.py`):
  - charpoly and every method built on it (eigenvals, eigenvects,
    is_diagonalizable, diagonalize, jordan_form, exp, log, analytic_func,
    singular_values, pinv, is_nilpotent, left_eigenvects,
    bidiagonal_decomposition) raise
    `PolynomialError('non-commutative expressions are not supported')`, as in
    SymPy. Before, the native kernels treated the entries as commuting and
    returned wrong values.
  - echelon_form, QRdecomposition and QRsolve raise NotImplementedError
    (unsupported, not a result). Before, echelon_form wrote 0 for
    `A*C - C*A`, and QR used conjugate where an operator needs the adjoint.

Intentional corrections (SymPy is wrong):
- `M.inv(method='GE')` and `M.inv(method='ADJ')` with noncommutative
  entries are computed by LU. SymPy applies commutative formulas:
  - for `M = [[A, B], [C, A]]` it gives `A*(A**2 - B*C)**-1`-style entries;
  - with random 2x2 real matrices for A, B, C, the product of SymPy's
    inverse with the 4x4 block matrix differs from the identity by 9 to 13
    (max entry);
  - the LU inverse is off by 1e-15.
  `M.solve(rhs)` (default 'GJ') uses LUsolve for the same reason, for a
  square M only (otherwise NotImplementedError).
- `rref`, `nullspace` and `columnspace` of `[[A, B], [C, A]]`: SymPy's rref
  has two pivots in column 0 and entries with `(A*C - C*A)**-1`, which is
  not a reduced form. symox keeps its generic elimination.
- Piecewise series where SymPy drops branches. SymPy's
  `Piecewise._eval_as_leading_term` picks the default branch when the
  conditions are undecided, so `Order.contains` sees the wrong power:
  - `(1 + Piecewise((x, y > 0), (Piecewise((x**2, y < -1), (x**3, True)), True))).series(x, 0, 3)`
    is `1 + O(x**3)` in SymPy, which is wrong at y = 1, where the
    expression is `1 + x`;
  - `(1 + x*Piecewise((x, y > 0), (x**4, True))).series(x, 0, 3)` is
    `1 + O(x**3)` in SymPy (the true value at y = 1 is `1 + x**2`);
  - `x*Piecewise(...)` above is `O(x**3)` in SymPy;
  - `Piecewise((x, y < 1), (x**5, True)) + O(x**2)` is `O(x**2)` in SymPy,
    and `Piecewise(...)*O(x**2)` is `O(x**7)` (at y = 0 it is `O(x**3)`).
  symox keeps the per-branch terms.
- `(x*y).evalf(subs={x: 1})` is `x*y` in SymPy and `(-x*y).evalf(subs={x: 1})`
  is `-y`: `AssocOp._eval_evalf` ignores `subs`. symox substitutes and
  gives `1.0*y` and `-1.0*y`.

Not a gap:
- The order of the mapping returned by `polarify(expr, subs=True)` comes
  from set iteration over `free_symbols` in SymPy. It is hash-randomized:
  three runs of the same statement gave three orders. The mapping and the
  expression are equal.
- `solveset(A*x - B, x)` is `{B*A**(-1)}` and `solve(A*x - B, x)` is `[]`, as
  in SymPy. The disagreement between the two is the reference's own;
  `solve.rs` is unchanged.
- `adjoint(sin(A))` (`transpose(sin(conjugate(A)))`), the default
  `Matrix.inv` with nc entries (PolynomialError), the per-branch O-term in a
  Piecewise series and the power of a Piecewise already matched.

Remaining (documented, not fixed):
- `sin(Piecewise((x, x < 1), (x**5, True))).series(x, 0, 3)` raises
  NotImplementedError. SymPy gives
  `Piecewise((x, x < 1), (O(x**5), True)) + O(x**3)` through
  `Function._eval_nseries`, which is not ported for a Piecewise argument.
- `Piecewise((x, y < 1), (x**5, True)).series(x, oo, 3)` raises; SymPy
  returns the Piecewise unchanged.
- `exp(Piecewise(...)).series(x, 0, 3)` raises NotImplementedError in both,
  with different messages.
- `O(x**2).subs(x, cos(x))`: SymPy gives `O((x - pi/2)**2, (x, pi/2))`,
  symox gives `O(cos(x)**2)`. The new point comes from solveset of
  `x - cos(d)`, which is a ConditionSet in symox (solve area, out of scope).
- `(1 + Piecewise((x, y > 0), (x**4, True))**2).series(x, 0, 3)` prints
  `Piecewise(...) + 1` where SymPy prints `1 + Piecewise(...)` (same value).
- `adjoint(sqrt(A))` and `adjoint(sqrt(x))`: SymPy keeps
  `adjoint((re(A)**2 + im(A)**2)**(1/4))` unevaluated, symox evaluates it
  to the conjugate (same value, different form).
- `QRdecomposition` of nc matrices: SymPy hits RecursionError for
  `[[A, B], [C, A]]` and gives `A/Abs(A)` for the diagonal case, which is
  wrong for operators; symox raises NotImplementedError.
- The pretty host hook supports one-line forms only; multi-line pretty
  forms of other host classes are not ported.

Probe (statements run in SymPy 1.14 and symox, same/total):

| Probe | Before | After |
|---|---|---|
| Piecewise series, Order, polarify (pw) | 36/44 | 36/44 |
| nested Piecewise series (pw2) | 35/48 | 37/48 |
| noncommutative (nc: limit, nsimplify, evalf, matrices) | 154/165 | 161/165 |
| nc matrices (mat) | 14/33 | 24/33 |
| partial evalf (ev, ev2, ev3) | 20/62 | 61/62 |
| adjoint printing, nc Piecewise (extra) | - | 25/25 |

All remaining differences are listed above.

Tests: `tests/test_pw_nc_parity.py`.

## matnum (wave 15): Float matrix functions and pivot zero tests

Implemented:
- Float matrices: `exp` and `log` use SymPy's Jordan route. `exp`
  computes `jordan_form`, then `_eval_matrix_exp_jblock` on each block,
  then `P*eJ*P.inv()`, and takes `re` when every value is real. `log`
  runs `cancel` before and after. Before this change, the native Hermite
  interpolation raised NotImplementedError except for spectral and small
  antisymmetric cases. A Float with prec > 53 keeps its precision.
  Exact and symbolic matrices still use the native path.
- `inv` of matrices whose entries are Rational, Float or complex with
  Float parts: a port of SymPy's `RR`/`CC` `inv_den`. This is dense
  Gauss-Jordan (`ddm_irref`) with partial pivoting at the largest entry
  precision. Zero entries come back as an exact `0`. Values match
  SymPy's bit for bit, including the last digit of the exp results above.
- `is_diagonalizable` of Float matrices follows SymPy's
  `_is_diagonalizable`: real symmetric, then Hermitian, then the mpmath
  eigenvectors. The native test raised NotImplementedError.
- mpmath eigenvectors: entries are converted exactly from `_mpf_`, and real
  entries stay `mpf`. This is what `mp.matrix(M.evalf(...))` does. The
  decimal-string round trip changed the eigenvector signs for prec-30
  matrices.
- `diag**e` with a non-integer `e` gives exact `0` off the diagonal, as
  SymPy's `a._new(..., lambda i, j: a[i,j]**exp if i == j else 0)` does.
- Pivot zero tests (`rref`, `rank`, `nullspace`, `echelon_form`,
  `columnspace`, `inverse_GE` and the other users of `pivot_simp`): when
  expand, cancel and the Pythagorean rewrites leave a candidate
  undecided, and the candidate is not a rational function of its symbols,
  the native search now asks a host hook. The hook is
  `symox._matrix_pivot.pivot_simplify`, which runs `simplify` and then
  `equals(0)`, as SymPy's `_find_reasonable_pivot` passes 2 and 3 do.
  Example: `Matrix([[cos(2*x) - cos(x)**2 + sin(x)**2, 1], [0, 0]])` used
  that entry as a nonzero pivot, so `rref` and `nullspace` were wrong
  (`-1/(...)`). They now give SymPy's `(Matrix([[0, 1], [0, 0]]), (1,))`
  and `[Matrix([[1], [0]])]`. Rational-function candidates skip the hook,
  because `cancel` is already canonical for them. A thread-local guard
  stops the zero test, which calls `pivot_simp` itself, from re-entering
  the hook.

Probe: 230 statements against SymPy 1.14. Before: same 178, diff 36 (on
the first 214 statements). After: same 222, diff 8. No statement went
from same to different.

Intentional divergences / remaining:
- `Matrix([[0, 1.0], [-1.0, 0]]).log()`: SymPy gets `0.5*pi` because its
  `cancel` expands the generators (`cancel(log(1.0*I))` gives `I*pi/2`).
  symox's `cancel` keeps `log(1.0*I)`. This is a scalar `cancel` gap,
  not a matrix one.
- Symbolic `exp` (native Hermite interpolation, not SymPy's Jordan
  route) gives equal values in other forms. Examples:
  `[[a,1],[0,b]]`, `[[x,1],[1,x]]`, `[[0,x],[-x,0]]`, and
  `[[0,r],[-r,0]]` with real or positive `r`, where SymPy gives
  `cos(r), sin(r)`.
- `singular_values()` and `charpoly()` of entries such as
  `sin(x)**2 + cos(x)**2`: SymPy's `charpoly` runs `simplify` on each
  EXRAW Berkowitz coefficient, which gives `[2, 0]`. symox keeps the
  unsimplified coefficients.
- Bareiss `det` still uses the native zero test without the host
  `simplify`.
- Domain inference and parameter case splits for `inv`, `det` and
  `solve` with symbolic entries (part (3) of the task) were not reached.

- Discrete Markov chains at symbolic times use `T**(a - c)`, as SymPy's
  `trans_probs**(rv.key - given.key)` does, instead of the
  diagonalization kernel. For a Float transition matrix,
  `P(Eq(Y[a], 2), Eq(Y[c], 1))` raised NotImplementedError before. It
  now gives SymPy's Piecewise, including the `1e-16*0.3**(a - c)` term.

Tests: `tests/test_matnum_parity.py`. Two older tests in
`tests/test_matrix_functions.py` and `tests/test_matrix_log.py`, and one
in `tests/test_missing_stats_wave4_parity.py`, asserted the previous
NotImplementedError. They now assert SymPy's values.

## simp-rest6 (wave 15)

The SymPy simplify test probe went from 1797 same / 108 diff / 39 not
reached / 4 exc to 1822 same / 101 diff / 23 not reached / 2 exc. No
statement went from same to different.

Fixed:
- `powdenest(..., polar=True)` is ported. `mellin_transform(besselj(0, x),
  x, s)` gives `2**(s - 1)*gamma(s/2)/gamma(1 - s/2)` as in SymPy.
- `simplify` keeps `2**(s - 1)*y` and `2*3**(x - 1)`. powsimp's
  `multiplicity` subtracts the denominator multiplicity, so `powsimp(2*3**x/3)`
  gives `2*3**(x - 1)`. powsimp returns matrix symbols unchanged and
  rebuilds other matrix expressions from their arguments.
- Noncommutative `expand(c**-1*a**-1)` keeps `c**(-1)*a**(-1)`. Before this
  change it returned `(c*a)**(-1)`, which is the wrong value.
- `nc_simplify(2.0*cos(pi/3)*X)` gives `X` because the coefficient test is
  `equal_valued`.
- MatPow:
  - `MatPow(MatPow(M, 2), -1)` combines to `M**(-2)`.
  - `MatPow(..., -1).doit(inv_expand=False)` returns the Inverse without
    expanding it.
  - MatMul combines a nested product with a following power of its inverse
    before flattening, which is the order of SymPy's canonicalize.

  `nc_simplify(MatMul(B**-1, A**-1, (A*B)**2))` gives `A*B`. Before this
  change it hit a RecursionError.
- `besselsimp(-besselj(1, x) - besselj(-1, x))` gives `0`. Before this
  change it gave the unevaluated form `0/x`.
- Mul with zoo: a number factor that is positive or negative is dropped, and
  so is one whose numeric zero and real facts are known. `zoo*sqrt(5 -
  sqrt(5))` is `zoo`.
- Pow.__new__ negative-base rule: when the base is a number with a negative
  coefficient and the exponent is an integer Symbol of known parity, the
  sign is pulled out. `(-I)**n` is `-I**n` for odd n, `(-2*I)**m` is
  `(2*I)**m` for even m.
- New compatibility modules: `symox.simplify.gammasimp` and
  `symox.polys.domains.realfield`. These modules let the hyperexpand
  lerchphi tests and `test_meijerg_with_Floats` run.

Remaining (documented, not fixed):
- `radsimp(num/den)` for the issue 15119 radical example: SymPy gives
  `nan`, symox gives `zoo*num`. Here den is exactly 0 and `num.is_zero` is
  True in both. The step where SymPy turns 0/0 into nan was not traced.
- `integrate(cos(x)/x, (x, 1, oo))`: the value is equal to SymPy's, but the
  form after hyperexpand groups the EulerGamma and log terms differently.
  The difference is in how the residue/series result is collected. Not
  ported.
- Differences that are form or order only:
  - `test_cse` (75);
  - the dummy numbering of `fu` `hyper_as_trig`;
  - issue 25590, issue 23543, separatevars, `Eq(1, x)` and Boole orderings;
  - reduce_inverses_nc_pow;
  - collect_const L362, fraction L444, issue 13143;
  - powsimp_negated_base `I*(-n)`.
- `Derivative` with a symbolic order (issue 22210) is unsupported because
  core derivative orders are u32.
- `MatrixSymbol.args[0]` is a str, not a Str. `Hyper_Function.args` are
  tuples, not Tuples.
- hyperexpand:
  - plan/ushift crashes;
  - shift_operators timeouts;
  - partial_simp and reduction_operators, whose parameters are random.
- `test_nsimplify` L409, the polymorphism NotImplementedError and
  exptrigsimp_noncommutative are still different.

Tests: `tests/test_simp_rest6_parity.py`.

## audit15 (wave 15)

Replay of lots c1 (stochastic processes, to_meijerg, symbolic fps), c4
(physics, tensor, vector, categories) and the printing "a" lots against
SymPy 1.14. Counts (examples_match / total):
- c1: 27/36 -> 27/36;
- c4: 37/37 -> 37/37;
- a: 127/199 -> 138/199. No object went from match to difference.

Implemented:
- The reference matrix postprocessor after `subs` and `xreplace`: an `Add`,
  `Mul` or `Pow` with matrix operands becomes `MatAdd`, `MatMul` or
  `MatPow`. `Add(x, y).subs({x: A, y: B})` is `A + B`, the sequential
  `(x*y).subs({x: A, y: B})` is `B*A` (the matrix put into the scalar
  coefficient of `MatMul(y, A)` becomes its first factor), the simultaneous
  `(x*y).xreplace({x: A, y: B})` is `A*B`, and mismatched shapes raise
  ShapeError. `x.subs(x, A)` returns the MatrixSymbol itself.
- `ArraySymbol` and `Reshape(...).as_explicit()` elements print as
  `A[0, 1]` and `{{B}_{1, 2}}` in LaTeX.
- `pprint` of `PermutationMatrix` and `MatrixPermute`.
- `pinv_solve` free parameters are Dummies named `_w0_0`, as in SymPy.
- `beta` prints `B(a, b)` / `Β(a, b)` and `\operatorname{B}\left(a, b\right)`
  (stats Beta, BetaPrime, FisherZ and StudentT densities).
- `FiniteSet ∩ ImageSet`: the reference `_handle_finite_sets` with image set
  membership. `{0, ..., 10} ∩ ImageSet(Lambda(x, x**2), Naturals)` is
  `{1, 4, 9}`; undecided symbolic elements stay in an unevaluated
  Intersection, as in SymPy.

Intentional divergences:
- `Subs(x, x, 0).subs(x, y)`: SymPy gives `has(x)` True and `has(y)` False
  through its `@cacheit` cache, because `Subs(x, x, 0) == Subs(y, y, 0)`.
  symox gives False and True, the values of `Subs(y, y, 0)` without the
  cache.
- DDM/SDM display of MPQ entries (already documented).
- The ShapeError message text differs; the exception type is the same.

Remaining (documented, not fixed):
- c1: DiscreteMarkovChain with `Str` states (NotImplementedError, and
  `Eq(Str, ...)` fails in the core), MultivariateNormal with MatrixSymbol
  parameters (NotImplementedError). The other c1 items were already
  documented (apply_finite_diff order, approximants, Expectation order,
  GaussianSymplecticEnsemble, Moyal, WienerProcess, Sum.euler_maclaurin).
- AlgebraicField alias display, Domain type names, UndefinedPredicate,
  set and dict print order, last digits of catalan and of jn/yn at complex
  float points.
- Triangular piecewise branch order (same value), JointEigenDistribution
  repr (`JointDistributionHandmade(...)` in SymPy), level_spacing_distribution
  Dummy name (`_s` versus `_spacing`), solveset ConditionSet dummy `_R` for
  an imaginary symbol.
- `sin(A)` for a MatrixSymbol raises NotImplementedError; SymPy builds
  `sin(A)`. (`sin(x).subs(x, A)` already gives `sin(A)`.)
- wigner_9j (the reference value in the doctest is wrong), Matplotlib
  backend output, devise_plan operator choice.

Tests: `tests/test_audit15_parity.py`.

## matparam (wave 16)

Domain inference and parameter routes for matrices, compared with SymPy 1.14
on a 308-statement probe (symbolic parameters, `QQ<sqrt(2)>`, Gaussian
entries, Floats mixed with symbols, identity entries, errors). Probe:
273 same / 35 different before, 301 same / 7 different after (28 fixed,
0 regressions).

Implemented:
- `to_DM()` without a domain follows `RepMatrix.to_DM`: `construct_domain`
  over the structurally nonzero entries (so `0.0` is kept), EX mapped to
  EXRAW. Domains such as `ZZ[sqrt(a)]`, `ZZ(exp(a))`, `ZZ_I[a]`, `RR[a]`,
  `ZZ[pi]`, `ZZ[sin(x)]` and EXRAW for `[sqrt(2), 1.0]` now match. An
  unsupported option raises `OptionError` as in SymPy.
- `inv()` over an inexact composite domain (`RR[a]`, `CC[a]`): the
  reference route (exact conversion, `inv_den`, conversion back, division by
  the denominator), including `NonInvertibleMatrixError` for a singular
  matrix.
- `charpoly()` over EXRAW applies `simplify` to the Berkowitz coefficients
  (so `sin(x)**2 + cos(x)**2` entries give `_x**2 - 2*_x + 1`); over an
  inexact composite domain the polynomial keeps that domain (`RR[a]`). The
  default charpoly is evaluated under `evaluate(True)`. `eigenvals()` and
  `singular_values()` follow the simplified charpoly when it differs from
  the native one (`{2: 1, 0: 1}`, `[2, 0]`), and Float symbolic matrices
  get their eigenvalues through the roots route.
- `cancel()` expands `log` atoms with Float arguments, as `sring` does:
  `cancel(log(1.0*I))` gives `I*pi/2`, and
  `Matrix([[0, 1.0], [-1.0, 0]]).log()` gives `[[0, 0.5*pi], [-0.5*pi, 0]]`.
- `Matrix.exp()` takes the Jordan route for symbolic matrices too (the
  route was already used for Float matrices), so `[[0, r], [-r, 0]]`,
  `[[x, 1], [1, x]]` and `[[0, 1], [-x**2, 0]]` give the reference
  exponential forms. Exact numeric matrices keep the native Hermite route.
  Declared nonfinite symbols keep the native `ValueError`, where SymPy
  returns `exp(z)`.
- The Bareiss det zero test with hidden-zero entries (`cos(2*x)` forms)
  already matched SymPy.

Remaining (documented, not fixed):
- `eigenvals()` dict order for some algebraic and Gaussian matrices (roots
  order in the core); `CC[a]` eigenvalues are the same values in another
  sqrt normalization.
- `QRdecomposition()` and `pinv()` forms for symbolic entries (`Abs` versus
  `conjugate` products, an unexpanded denominator in the reference).
- `Matrix([[1.0*a, 1], [1, 1.0*a]]).gauss_jordan_solve(...)` form (same
  values), and `Matrix([[1.0, sqrt(2)], [2, 3.0]]).inv()` gives
  `2.0*sqrt(2)` where SymPy gives `2*sqrt(2)`.
- `construct_domain([0.0, a, 1, 1])` gives `RR[a]` (SymPy `ZZ[a]`).
- `cancel(log(4))` stays `log(4)` (SymPy `2*log(2)`); only Float log
  arguments are expanded.
- `exp()` of exact numeric matrices (`[[1, 2], [3, 4]]`,
  `[[sqrt(2), 1], [1, 0]]`, `[[0, 1], [-1, 0]]`) has other forms than the
  reference (same values).
- `eigenvects()` does not follow the simplified charpoly of identity
  entries.
- SymPy raises `TypeError` for `Matrix([[sin(x), 1], [0, log(x)]]).charpoly()`;
  symox returns the polynomial.

Tests: `tests/test_matparam_parity.py`.

## rrcc (wave 16): RR/CC and Gaussian composite domains

Probe of 269 statements against SymPy 1.14 (mixed RR/ZZ[y], several
generators over RR[y], apart/invert/gcdex, explicit CC, RR(x) fields,
Gaussian sqf_list, CC[y] inference, `.rep`, precision with RealField(100)
and 30-digit Floats): 220 identical before, 233 after, plus 2 intentional
corrections. No statement that matched before differs now.

Implemented:
- `invert`, `gcdex` and `half_gcdex` over RR above 53 bits (expressions and
  Poly methods) keep the precision of the inputs, as SymPy:
  `invert(Float('1.5', 30)*x, x**2 + 1)` is
  `-0.666666666666666666666666666667*x` (it gave `-2*x/3`).
- `Poly.trunc` over RR: the floor modulo of every coefficient at the domain
  precision (it raised TypeError).
- Zero Float coefficients make the domain RR (`Poly(0.0, x)`,
  `Poly([0.0, 1], x)`, the largest precision met), and the zero polynomial
  over RR gives `0.0` for `LC`, `TC`, `nth`, `coeffs` and `terms`.
- `numer` and `denom` of RR(x) and CC(x) elements (they returned the whole
  fraction and 1, so printing `x/(x + 1)` in `RR.frac_field(x).field`
  raised CoercionFailed).
- Expression `sqf_list` with float factors builds each component on its
  own domain, as `_symbolic_factor_list`: `sqf_list(2.0*(x - y)**2*(x + y))`
  is `(2.0, [(x + y, 1), (x - y, 2)])` (factors were `1.0*x + 1.0*y`).
- RR or RR[y] against ZZ[y] or QQ(y) (`+ - *`, gcd, compose, rem) unify to
  RR[y]/RR(y) instead of raising TypeError, in one or several generators
  and in both operand orders.
- `content` and `primitive` over RR[y] in several generators: the monic gcd
  over the exact QQ[y], as `dmp_inner_gcd` on inexact grounds
  (`Poly(1.5*y*x*z + 3*y, x, z).content()` is `1.0*y`; it gave `3/2`).

Intentional correction:
- ZZ_I `sqf_list`: SymPy keeps the factors of `dmp_sqf_list` but loses their
  unit, so `Poly((I*x + 1)**2, x).sqf_list()` is `(1, [(-x + I, 2)])`, whose
  product is `-(I*x + 1)**2`. symox keeps the same factors and returns the
  unit that makes the product exact: `(-1, [(-x + I, 2)])`. Same for
  `(I*x + y)**3`, `(I*x + 2*y)**2` and `(2*I*x + 1)**2*(x + I)`; the tests
  check the expanded product.

Remaining (documented, not fixed):
- Explicit `domain=CC` or `ComplexField(prec)` as a Poly domain (and
  `set_domain(CC)`, `unify` with a CC Poly, `factor(..., domain=CC)`, and
  `factor(x**2 - 2, domain=RR)`): NotImplementedError, there is no plain CC
  coefficient backend. The inferred CC Poly (`Poly(x + I/2, x)`) is the
  legacy RR representation with `I` as a domain generator: its gcd,
  invert and gcdex with other CC polys raise TypeError or ValueError.
- gcd of multivariate RR[y] polynomials (ValueError, only ring arithmetic
  is supported in several generators), the lcm scaling over RR[y] in several
  generators (`3.0*x**2*z + ...` in SymPy, `1.0*x**2*z + ...` in symox,
  the same lcm up to a unit), `integrate` over RR[y] (SymPy moves to RR(y)),
  and the sqf_list unit over RR[y] in several generators.
- `RR.frac_field(x)` elements are stored as normalised expressions, so
  `x/2` prints `0.5*x` (SymPy `x/2.0`).
- Zero Floats at high precision print `0.00000000000000000000000000000`
  (SymPy `0.0`): core Float printing.
- `Poly(-(x - I*y)**2, x, y)` builds the legacy `ZZ[I]` domain instead of
  `ZZ_I` after an earlier `sqf_list` or `factor_list` of the same
  expression (state-dependent, also on main; the expanded form is right).
- The TypeError text of `trunc` over QQ names `QQ` (SymPy `PythonMPQ`).

Tests: `tests/test_rrcc_parity.py`.
