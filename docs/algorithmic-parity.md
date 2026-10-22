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
- `is_complex` of a finite non-real sum such as `1 + I` is unknown. The
  finiteness of `exp(f)` for finite `f`, and of `I**(1/3)`, are unknown.
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
