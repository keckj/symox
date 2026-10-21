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
  `Poly` over `ZZ_I[x]`, `QQ_I[x]` or `QQ<sqrt(2)>[x]` raises
  `NotImplementedError`; expression-level `factor` or `gcd` with a
  polynomial-ring `domain=` (e.g. `ZZ[x]`, `ZZ[x][y]`) raise
  `NotImplementedError` (the `Poly` methods work); adding a `ZZ[x][y]` and a
  `ZZ[x]` polynomial raises `TypeError` instead of `GeneratorsError`;
  `ZZ_I[x].is_unit(I)` is `True` (reference `False`);
  `QQ<sqrt(2)>(sqrt(2))` is accepted where the reference rejects it; inferred
  domains with `I` and parameters stay `ZZ[y,I]` where the reference infers
  `ZZ_I[y]` (resultants agree); fraction fields over `ZZ_I`, `QQ_I` or
  algebraic fields are descriptor-only and raise `NotImplementedError` on
  element construction.
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
  Remaining: `RR[pi]` raises `NotImplementedError`; `QQ`, `QQ(x)` and
  `QQ<sqrt(2)>` reject `convert(2.5)` that the reference converts exactly;
  `QQ_I.convert(1j)` fails; sets and dicts raise `AttributeError` or are
  accepted where some reference domains raise `CoercionFailed`;
  `CC.invert` raises `NotImplementedError` (reference `TypeError`);
  `EXRAW` elements are wrapped and do not compare equal to the raw `Expr`.
- `(-2)**(1/6)` with `extension=True` yields the larger field
  `QQ<2**(1/6) + (-1)**(1/6)>` because core power canonicalization splits the
  base; results stay exact but the field is not the reference one.
- `construct_domain` elements of composite domains print with their parent,
  e.g. `ZZ[pi](pi)`, where the reference prints `pi`;
  `construct_domain([1], extension=sqrt(2))` raises `NotImplementedError`
  where the reference returns `ZZ`.
- `EX` root forms and `nroots` digits can differ from the reference (symox
  can be more accurate); exact reduction of root objects (`CRootOf` and
  `RootSum` simplification into an algebraic field) is not ported.

## Remaining blocking gaps

General noncommutative scalar algebra remains unsupported in this snapshot:
`Symbol('a', commutative=False)` raises `ValueError` in the declaration layer.
Dedicated tensor or quantum object protocols do not establish a general
noncommutative scalar engine; neither do user classes with
`is_commutative = None`, whose products keep their written order (see the
paragraph on user subclasses below). Exact arithmetic reuse is restricted to ordinary
commuting symbols and does not change this pre-existing rejection. This gap
must be resolved before certifying general `Symbol` and scalar algebra parity.

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
ordering only when the difference is extended real. Sums, products and powers
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
- Noncommutative symbols are not implemented.
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
- `Heaviside(x)` stores no `H0` argument (the reference args are
  `(x, 1/2)`), and `Heaviside(x, H0)` for another `H0` evaluates to
  `Piecewise((H0, Eq(x, 0)), (Heaviside(x), True))`, so its printed code is a
  nested conditional with the same values. Arithmetic with a string operand
  sympifies the string, where the reference raises `TypeError`.
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
reference raises `IndexError`; `Heaviside(x).args` omits the reference
`1/2` default; `Lambda.nargs` is a raw native set node. SymPy 1.14 builds a
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
otherwise) and `meijerg=False` skips it. `risch`, `heurisch` and `manual`
other than None still raise NotImplementedError: the native engine does not
separate those methods, so it cannot honor them faithfully. An improper
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

`meijerint_indefinite` is ported (`_meijerint_indefinite_1`, the splitting
points, the hyperbolic rewrite and the choice of the simplest result) and
is the last stage of indefinite integration unless `meijerg=False`.
`meijerg=True` uses only it. An error inside the method (a failing
`hyperexpand`) gives no result here, where the reference lets the exception
out of `integrate`.

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
product is not simplified again. Remaining differences:
- the native solver still checks the union against the whole product, so a
  root of one factor where another factor is undefined is rejected:
  `solve(x*log(x), x)` is `[1]` and `solve((x - 1)*log(x - 1), x)` is `[2]`,
  where the reference returns `[0, 1]` and `[1, 2]`; the extra reference
  values are not zeros of the expression;
- `solve(x**5 + a*x + 1, x)` raises `NotImplementedError`, where the
  reference returns the empty partial result of `roots` (its `all_roots`
  does not support symbolic coefficients): no root was proven absent;
- `simplify` lacks the reference step
  `shorter(expand_log(e, deep=True), logcombine(e))`, so the roots of a
  power of `exp` keep `log(sqrt(2))` where the reference has `log(2)/2`
  (`exp(2*x) - 2`, `exp(-2*x) - 2`, `exp(2*x + 1) - 2`, `exp(2*a*x) - 3`,
  `exp(x) - exp(2*x)/3 - 1`, `exp(2*x + a) - 1`);
- `solve(x**3 - 3*x + 1, x, simplify=True)` differs in form: the native
  simplify search does not reach the reference forms of the Cardano roots
  (equal values);
- the native Pow does not fold `sqrt(2**(1/3)/2)` to `2**(2/3)/2`, so
  `exp(-3*x) - 2` gives `log(sqrt(2)*sqrt(2**(1/3))/2)` (equal value);
- the native `log` takes `I*pi` out of a provably negative argument, so one
  root of `exp(3*x) - 3*exp(x) + 1` prints `log(A) + I*pi` where the
  reference keeps `log(-A)` (equal value);
- for `2**(3*x) - 3*2**x + 1` the reference prints one solution with the
  logarithm of an exactly negative number in a form that its evalf puts at
  `-I*pi`; the native value is the principal one, `+I*pi`, so the two
  solutions differ by `2*I*pi/log(2)` (both are solutions; the other
  forms differ by the simplify gap above);
- the reference rewrites a perfect-power integer base before solving
  (`4**x - 8` and `(1/4)**x - 8` give two solutions, through `2**x`); the
  native solver keeps the base and gives one (`3/2`, `-3/2`). `6**x - 8`
  gives `log(8)/log(6)` where the reference has `3*log(2)/log(6)`;
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
- Piecewise results keep `Abs(w**2) < 1` where the reference `simplify`
  rewrites it as `(w > -1) & (w < 1)` (sine and cosine transforms of
  `besselj(0, t)` and `sin(t)/t`, the inverse Fourier transform of
  `sin(k)/k`). The reference rewrite assumes a real `w`. The nested
  Piecewise of the cosine transform of `sin(t)/t` is folded.
- The native definite integral of `exp(-a*x)*sin(b*x)*Heaviside(x)` over the
  real line is not computed. The table gives the reference value with the
  condition `True` instead of the reference half plane condition.
- `fourier_transform(exp(-a*x**2), x, k)` with an unrestricted `a`: the
  reference holds it (the integral of the function is not computed), the
  native table returns the value with the condition `re(a) > 0`.
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
  heurisch and Meijer G methods available in SymPy. Meijer G definite
  integration, its convergence conditions and `hyperexpand` are ported (see
  above), and so are the Meijer G indefinite integrals and the Piecewise
  integration of `Integral.doit`. The `risch`, `heurisch` and `manual`
  method selections and `NonElementaryIntegral` results (which need a Risch
  proof) remain open: those options still raise NotImplementedError, and
  results the reference only reaches through heurisch can differ in form.
  The integration corpus (`tools/integration_corpus.py`) now reports 1033
  identical and 89 equivalent results, against 988 and 93 before the
  Piecewise and Meijer G indefinite ports; 54 cases stay unevaluated where
  the reference answers. Those open cases are oscillatory integrals on
  `(0, oo)` that the reference returns as `AccumBounds`, divergent values
  (`1/sin(x)` on `(0, pi)`, `x**a` on `(1, oo)`), `1/cosh(x)` over the real
  line, and finite intervals with an unassumed parameter in a logarithmic,
  `Heaviside` or `DiracDelta` integrand (`1/(x+a)`, `1/x` on `(a, 1)`,
  `log(a*x)`). For the last group, some reference values such as `-log(a)`
  for `1/x` on `(a, 1)` hold only when the interval avoids the pole, so the
  native result stays held rather than copying them. The remaining
  numerical differences are reference errors: `laplace_transform(1/(t+1))`
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
  Float of lower precision. Remaining divergences: an exhausted sum with no
  correct bit (an exact zero, or `exp(300)+1-exp(300)` with `maxn=10`) stays
  unevaluated instead of becoming SymPy's precision-1 noise Float such as
  `-0.e-124`; evalf without options keeps the native 4096-bit budget and
  resolves cancellations that SymPy's default `maxn=100` reports as noise;
  the accuracy estimate is a rigorous rounding bound rather than SymPy's
  `add_terms` estimate, so achieved precisions can differ by a few bits;
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
  drawn from SymPy's ODE tests, classification output matches for 295, the
  test's hint for 237 and the default solution for 245 (string comparison).
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
  (SymPy does not return for them). It is also much slower
  than SymPy on the large expressions of the Lie group heuristic
  `abaco2_unique_general`, so a default `dsolve` that reaches it (SymPy
  raises after about 80 s) can take far longer. Deliberate
  extensions are the logarithmic Frobenius basis for repeated or resonant
  indicial roots, `ics` that name the generated constants, exact
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
