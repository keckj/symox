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
subset is incomplete. These options do not establish all solver-option parity.

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
and a one-equation list. Undecided parameter-only equations retain generic
candidates rather than proving inconsistency. Full parameter strata remain
unproved. A single expression target such as `f(x)`, `exp(x)` or `x**2` can be
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

## Remaining blocking gaps

The versioned per-function ledger retains its own native snapshot. The
downstream candidate does not renew those historical certificates. Current
upstream and regression results are recorded separately in the
[integration report](downstream-integration.md) and its candidate-specific
validation evidence.

The downstream fixes also cover exact Python `Fraction` conversion of native
integers and rationals, two-dimensional NumPy/matrix arithmetic, external
`_sympy_` conversion hooks, and relational constructor identity and instance
checks. Cyclic conversion hooks raise an exception before exhausting the native
stack. These protocols do not establish general Python subclass compatibility:
custom symbolic subclasses, their metadata and container type predicates remain
open obligations. Dense matrix subclasses can be constructed and override Python
methods; preservation of those subclasses across every matrix operation is not
yet certified.

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

Implicit external conversion never invokes a native object's reference exporter.
Conversion chains are traversed iteratively, with retained owners and identity
checks for cycles; reentrant hooks are guarded separately. Matrix `nsimplify`
applies scalar recognition entrywise. These changes do not certify all
recognition options or arbitrary external object protocols. `exp` on a concrete
matrix dispatches to the native matrix function and returns immutable output;
general floating-point matrix eigenvalue support remains incomplete.

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
reported as callable, including after importing differential geometry. This does
not certify the full nominal function or symbolic subclass hierarchy.
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

- The integration engine does not implement the full collection of Risch,
  heurisch and Meijer G methods available in SymPy.
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
- `evalf(strict=True)`, user-controlled `maxn` and the full set of numerical
  evaluation options are not implemented. The former 10,000-decimal-digit cap
  has been removed; native precision still uses a bounded machine counter with
  space reserved for working precision. Binary exponents also use a bounded
  native counter; out-of-range arithmetic must retain an exact expression
  rather than wrap to an unrelated number. Explicit rational output retains
  allocation budgets even though Float literals have compact storage at much
  larger exponents. Removing policy caps does not certify all numerical
  algorithms, counter ranges or precision options.
- Complex algebraic roots do not yet have certified complex isolating boxes;
  complex root ordering and numerical refinement retain the limitations
  described in [the root guide](algebraic-roots-quaternions.md).
- General transcendental solving, all solver options, symbolic parameter case
  distinctions and full set-valued solution-family coverage still need a
  versioned behavioral audit. Supported Lambert families do not establish parity
  for arbitrary mixtures of powers, exponentials and logarithms.
- Wildcard matching does not yet accept the reference `Wild(properties=...)`
  predicates used by ChemPy. Predicate checks must preserve valid alternatives
  during matching; ignoring them or rejecting only the first completed match
  would not establish equivalent behavior.
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
