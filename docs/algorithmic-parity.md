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
order, and `_as_mpf_val` gives `float()`. Remaining gaps:
an instance that is no longer referenced by any Python object is rebuilt
without its instance attributes (constructor results produced inside a native
operation are retained for a bounded window); `Basic.__new__` sympifies its
arguments, whereas the reference stores raw Python values; native assumption
queries, `_eval_subs`, overridden `func`/`free_symbols` and nested Python
`_print_<Class>` printer methods are not consulted by native algorithms; the
reference treats `is_commutative = None` as noncommutative (which changes
`Add`/`Mul` ordering), while native nodes are commutative; `Wild` subclasses
are not hosted; `And`/`Max` are not classes, so `singledispatch` registration
and `LatticeOp.make_args` on them remain unsupported; native operators on a
tree that already contains such a `Number` do not refold its coefficient
(`(w*x)*3`), and native results equal to `oo`, `-oo`, `nan` or `zoo` are not
the `S` singleton objects, so identity tests such as `-oo is S.NegativeInfinity`
fail.

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
see the attributes patched onto the node class of that name.

The explicit `NonElementaryIntegral` constructor builds a held native integral
and retains its marker for reconstruction, pickle and the covered root
structural operations. Constructing this marker does not prove that its
integrand is non-elementary. The integration engine does not automatically use
it to label unsupported computations: the complete Risch decision procedure
remains unimplemented. Preserving that marker inside every composite native
expression remains a separate obligation.

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
inexact input, singular or discontinuous weights, products of deltas and
derivative orders; all established integration methods remain available.
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
now memoizes shared subexpressions within each recursive transformation pass.
Head discovery visits each distinct node once and combines the same head bits;
rewrite rules, acceptance conditions and candidate ranking are retained. Local
caches keep expression owners alive to prevent reuse of identity keys. Exact
identity and unchanged-expression tests cover graphs whose expanded trees have
more than `2**64` occurrences. Public traversal multiplicities are unchanged.
This removes repeated work in these passes, not the inherent cost of expression
growth or pattern matching. The separate public matching/Fu driver and the
complexity of the complete `simplify` pipeline remain separate obligations.

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
