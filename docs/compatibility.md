# SymPy compatibility

[README](../README.md)

## Implemented features

Implemented:

* expression DAG with cached hashes, canonical `Add`/`Mul`/`Pow` following
  SymPy's automatic simplification rules (like-term collection, power
  merging, exact radicals `sqrt(8) = 2*sqrt(2)`, `I**2 = -1`, numeric
  coefficient distribution, infinities and `nan`, ...)
* exact integers (machine word fast path, transparent big integers),
  rationals, machine floats and arbitrary-precision floating-point expressions
* elementary and special functions with automatic evaluation (special values
  of trig functions, parity, periodicity, inverse compositions, ...)
* `diff`, `expand` (sparse polynomial fast path), `subs`, `xreplace`,
  `evalf`, `free_symbols`, `has`, relationals, `Piecewise`
* `str`/`srepr` printers matching SymPy, a parser (`sympify`) accepting
  SymPy's `str` and `srepr` output
* `lambdify`: compiles expressions into a register machine with common
  subexpression sharing, vectorised over NumPy arrays (GIL released)
* Python bindings (PyO3, abi3 wheels for CPython >= 3.9), conversion to and
  from SymPy objects
* polynomials over integer, rational, finite-field, rational-function and
  supported algebraic domains: `Poly`,
  factorization, GCD, cancellation, partial fractions, roots, division and
  reduced Gröbner bases
* symbolic assumptions, simplification, integration, summation, equation and
  ODE solving, series and limits
* mutable, immutable and sparse matrices, abstract matrix expressions,
  exact linear algebra and supported eigenproblems
* Boolean algebra, normal forms, truth-table minimization, SAT, symbolic sets
  and supported set-based equation solving
* LaTeX and pretty printing, common subexpression elimination, native code
  printers and extensible C++17 generation
* optional PyTorch and JAX code generation in `symox.extra`, standalone source export, automatic
  differentiation and opt-in compilation with `torch.compile` / `jax.jit`
* portable Slang source and optional SlangPy execution in `symox.extra`; standard SYCL C++ source
  and asynchronous kernels for AdaptiveCpp
* isolated differential fuzzing, seeded replay, shrinking and triage

The sections below describe the supported subsets and remaining limitations.
The [public interface guide](public-interfaces.md), [root/quaternion guide](algebraic-roots-quaternions.md)
and [splines/variational guide](splines-variations-integer-equations.md) describe
the latest root exports. Their presence does not establish complete behavioral parity.

Special functions include combinatorial numbers, gamma/zeta families, error
and exponential integrals, Bessel/Airy functions, elliptic integrals, Lambert W,
hypergeometric functions and orthogonal polynomials, with symbolic evaluation,
derivatives and double precision real/complex numerics. `expand_func`,
`besselsimp`, `hyperexpand`, Stirling numbers and multinomial coefficients are
also available. Unsupported symbolic cases remain unevaluated.

`Float(value, digits)`, `expr.evalf(digits)` and `N(expr, digits)` preserve
arbitrary precision in native expression nodes, including arithmetic, partial
symbolic evaluation, matrices, parsing, pickling and SymPy interchange.
Decimal strings are parsed directly; Python floats retain their original
binary value. The native multiprecision evaluator supports elementary
real/complex functions and a subset of special functions. Unsupported functions
remain symbolic; they do not receive invented digits from a double fallback.
The current resource limits are 10,000 decimal digits or 40,000 precision bits.
Exponential and hyperbolic evaluation is limited to a real argument magnitude
of 2,000,000; power results are limited to a decimal exponent magnitude around
1,000,000. Outside supported ranges expressions remain symbolic.
Compiled numerical backends still use their declared machine dtypes.

`evalf(chop=True)` and `N(..., chop=True)` remove absolutely small real or
imaginary numeric components using the requested binary precision. A positive
numeric `chop` value selects a tolerance; `False` or zero disables chopping.
The same option works entrywise for dense and sparse matrices. Partial
symbolic coefficients retain ordinary `evalf` behavior. Threshold decisions use
native arbitrary-precision magnitudes, including tolerances beyond the range
of a machine float. Adaptive error estimates and the `maxn`, `strict=True`,
`verbose=True` and quadrature options remain unsupported.

`special_nstr(expr, dps)` returns an arbitrary precision decimal string for
numeric expressions supported by the multiprecision evaluator, including
factorial, binomial, rising/falling factorials, beta, Catalan, Fibonacci and
Lucas via exact gamma/elementary identities. Other special function families
return `None` when multiprecision evaluation is unsupported; no double
precision fallback is presented as extra significant digits.

## Scalar expressions and calculus

Assumptions distinguish real `nonzero=True` from complex `zero=False`, and
exact rational-power and supported quadratic facts filter solver candidates
without floating-point tolerances. Unquoted Unicode symbol names are accepted.
`Symbol` and `Dummy` retain all 31 standard facts in three-valued form.
`assumptions0` exposes their exact implication closure; contradictions are
rejected, and explicit false facts survive reconstruction, hashing, pickling
and SymPy interchange. Public argument ordering compares the complete inferred
symbol content. Noncommutative symbols and custom facts are explicitly
unsupported. Undefined-function heads also retain all 31 standard facts.
`AlgebraicNumber.is_nonzero` follows the same real-nonzero meaning as its
scalar expression view.
`Q` predicates and `ask(proposition, assumptions=True, context=...)` query
native three-valued facts, propagate Boolean assumptions with the Rust SAT
solver, and handle scalar arithmetic, elementary functions and supported matrix
properties. For example, `ask(Q.real(x), Q.positive(x))` returns `True`, while
`ask(Q.positive(x), Q.real(x))` returns `None`. Contradictory relevant assumptions
raise `ValueError`. `AssumptionsContext`, `global_assumptions` and nested
`assuming(...)` contexts are available; temporary contexts restore their previous
contents even when the body raises an exception. No SymPy import is required.
`refine(expr, assumptions=True)` and `expr.refine(...)` use these facts for
absolute values, real/imaginary parts, arguments, supported powers, signs,
`atan2` quadrants and Boolean predicates. Unsupported rewrites remain symbolic.
Predicate applications preserve structure, substitutions, matrix arguments,
pickles and SymPy interchange. `Predicate(name)` creates a distinct undefined
predicate that can participate in Boolean assumptions. Custom binary relations can register Python handlers. General linear real
arithmetic inference across inequalities remains incomplete; undecidable
queries may return `None`.

Numeric expressions follow Python's hash contract for equal built-in numbers
and `fractions.Fraction` values, so mixed dictionaries and sets work correctly.
Approximate nonzero integral values keep `is_integer=None`; fractional floats
prove `False`, while Float zero is known to be an integer without changing its
type. Proven nonintegers are neither even nor odd. The signs of infinities use
the extended-real predicates, and `nan` keeps its numerical assumptions unknown.
`Dummy` atoms carry a distinct identity index: two `Dummy("x")` objects remain
different through substitutions, differentiation, SymPy conversion and pickling.
`Symbol("_x")` remains an ordinary symbol. Code printers use distinct indexed
identifiers for dummy variables. Tensor, Slang and SYCL printers also accept
explicit symbol-name mappings.
`as_coeff_add`, `as_coeff_mul`, `as_coeff_Mul`, `as_independent` and
`as_content_primitive`
provide structural decompositions, including dependency selection and rational
or supported radical content. Floating-point presence and maximum precision
are cached per node, so repeated inspection preserves shared expression graphs
instead of traversing every occurrence again.
`as_base_exp`, `as_powers_dict`, `args_cnc` and `could_extract_minus_sign`
provide public factor and sign views, preserving held products and exact
exponents. Power dictionaries retain distinct dummy variables, including held
factors whose exponents cancel to zero. Leading coefficient extraction preserves
raw product tails. Negation preserves the factor order of held products, and
held sums are reclassified before regrouping.
`as_coeff_Add(rational=...)` preserves the leading coefficient and raw sum tail.
Scalar and matrix expressions expose `class_key()` and `sort_key()`; sort keys
retain shared children and support the default lexicographic order. Other named
orders are rejected explicitly.
Multiplication normalizes active held-product spines with a shared-node memo,
collecting numerical coefficients and classifying infinities before zero
shortcuts. Powers and function arguments retain their evaluation boundaries.
Early numeric exceptional powers and late complex infinity follow separate
collection phases, preserving held argument trees and coefficient roles.
`.is_number` classifies the held expression independently of finite evaluation:
infinities and NaN are numbers, while undefined functions and held matrix
contractions are not. Integral and Sum use their scoped free variables.
`collect` supports exact or power-based matching, coefficient dictionaries with
`evaluate=False`, recursive coefficient callbacks and order-term distribution.
`exact=None`, Wild patterns and multivariate derivative patterns are explicitly
unsupported.
Scalar equality, ordering, operation counts and pickles also preserve DAG sharing.
`matches(expr, repl_dict=None)` exposes native wildcard matching in the pattern
direction; legacy `old=True` algebraic matching remains unsupported.
`is_hypergeometric(k)` uses the native consecutive-term ratio simplifier and
returns `None` for expressions containing `Piecewise`. Boolean expressions
expose `to_nnf(simplify=True)`. `dir(x, cdir)` computes the first nonzero
derivative coefficient at zero with a search limit of 256 derivatives;
`Limit.dir` still exposes the held limit's direction string.
`as_dummy`, `fromiter`, `aseries`, `fps`, `fourier_series`,
`extract_branch_factor` and `is_meromorphic` complete the measured public
`Expr` method surface. Their supported cases, examples and resource bounds
are described in [series and expression analysis](series.md).
`atoms()` selects scalar nodes by class, function head or implicit atom class;
`NumberSymbol` distinguishes constants such as `pi` from ordinary numbers.
`is_polynomial`, `is_rational_function` and `is_algebraic_expr` inspect held
expressions without expanding powers, including arbitrary scalar generators
and calculus binding scopes. Unknown cases retain SymPy's three-valued result.
These inspection methods currently reject typed matrix parents explicitly.
`count_ops()` counts all occurrences in the mathematical expression using exact
integers, without expanding a shared tree. Assumption queries cache independent
structural and numerical results per node; closed double-precision evaluations
share a temporary cache across the query. Numerical solvers keep their changing
variable environments separate. Default sort keys retain shared tuple children.
Scalar pickles use flat native node records; existing recursive scalar pickles
remain readable.
Python `expand` supports independent `mul`, `multinomial`, `power_exp`,
`power_base`, `log` and `trig` hints, including `deep=False` and `force`.
Scalar and matrix expansion share a work budget and check conservative bounds
before allocating products or expanding powers. The default term bound is
100,000; `max_terms=` changes it. Exceeding a term or work bound raises
`ValueError`. Matrix expansion preserves held transposes, adjoints, inverses
and matrix powers while expanding their descendants when `deep=True`.
`Integral` and `Sum` are native held binding operators with ordered limits,
scope-aware free symbols and substitution, and exact SymPy interchange.
Their `doit()` methods use the native integration and summation engines;
unsupported cases remain held. Differentiation supports the fundamental theorem
and Leibniz rule for integrals and retains discrete endpoint derivatives.
Finite sums use Karr's convention for reversed endpoints.
`Product` is a native held binder with the same scoping, printing and SymPy
interchange; `product()` and `Product.doit()` evaluate constant, polynomial
(via rising factorials), rational, power and exponential-sum terms, multiply
out short definite ranges and use Karr's convention for reversed ranges.
Unsupported terms stay held. `Sum.is_convergent()`,
`Sum.is_absolutely_convergent()` and `Product.is_convergent()` run SymPy's
tests natively on one-limit series (ranges from `-oo` included): exact
rational and hypergeometric terms, then the divergence test, p-series and
`k**p*log(k)**q*log(log(k))**r` comparison, alternating series, bounded
factors times an absolutely convergent series, ratio, Raabe and root tests;
products go through `Sum(log(f))` and `Sum(f - 1)`. Undecided cases (such as
symbolic exponents) raise `NotImplementedError`. Deviations: SymPy's wrong
`False` for `Sum(sin(k)/k)` and `Sum(cos(k)/k)` raises here, and
`Sum(1 - cos(1/k))` is `True` (SymPy says `False`); `Sum.euler_maclaurin`
and `Sum.eval_zeta_function` are not implemented. `index`, `reverse_order`,
`change_index`, `reorder`, `reorder_limit` (raising `ReorderError`),
`has_finite_limits`, `has_reversed_limits`, `has_empty_sequence` and
`Product.term` follow SymPy's `ExprWithIntLimits`; `hypersimp` and
`hypersimilar` are native (the `hypersimp` ratio may differ from SymPy's
in `simplify` shape). `rsolve`, `rsolve_poly`, `rsolve_ratio` and
`rsolve_hyper` solve linear recurrences with polynomial or rational
coefficients and sums of hypergeometric inhomogeneities (Abramov's
polynomial and rational solutions, Petkovsek's hypergeometric kernel), take
initial conditions as a dict or list, return `None` where SymPy does and
raise SymPy's `ValueError`/`TypeError` on malformed input. Differences: an
inconsistent `rsolve_poly` system returns `None` (SymPy returns the bare
ansatz), nonlinear terms such as `y(n)*y(n - 1)` raise `ValueError`, some
results differ in `simplify` shape (`2**n*n/2` for SymPy's `2**(n - 1)*n`),
and with fewer initial conditions than constants the eliminated constant
may differ (SymPy's choice depends on set ordering).
`Lambda` accepts multiple variables, nested tuple signatures and empty signatures.
Its arguments retain the signature tuple; application uses exact simultaneous
replacement, following SymPy's binding behavior. Identity lambdas use the native
`S.IdentityFunction` singleton. Existing scalar and facade Lambda pickles remain
readable.
Undefined functions preserve all 31 standard true or false assumption facts,
canonical head identity and declared arity, including names shared with builtins.
An explicit `commutative=True` distinguishes a head from an implicit declaration;
applied facts retain unknown values and do not imply continuity. Custom facts
and noncommutative heads are rejected explicitly. Arities are normalized to
sorted, unique tuples. Existing function pickles remain readable.
Held `Derivative`, `Subs` and `Limit` expressions support evaluation and round trips; simultaneous symbol substitutions preserve
the original right-hand sides. Limits of unknown functions remain held when
continuity cannot be established.
`Subs` accepts symbol variables and identifies safe renamings of bound variables
when substitution points do not refer to those variables. Cases with such
references retain structural identity to preserve their evaluation semantics.
Simultaneous substitutions also accept algebraic patterns such as `x**2` and
`sin(x)`, masking all matched values before inserting their replacements.
Exact `xreplace` visits public additive terms, numeric coefficients and whole
powers in compact products. Fractional exponent substitutions require an integer
ratio or a provably safe positive/nonnegative real branch.

`tests/compare_sympy.py` runs every case through SymPy and symox and
compares the printed canonical forms and, when they differ, numerical
equivalence. Known differences and remaining gaps:

* Some antiderivatives use a different branch constant from SymPy (for example
  `integrate(exp(-x)/x, x)` returns `Ei(-x)`).
* Simplification: `trigsimp` (methods `matching` and `fu`, `inverse=True`),
  `futrig`, `exptrigsimp`, `fu` with a custom `measure`, and the `TR0` ...
  `TR22`, `TRmorrie`, `TRpower` rules of `sympy.simplify.fu` are native ports.
  Bounded Gröbner reduction and legacy adapters are also available; advanced
  strategy options and noncommutative inputs do not establish full parity.
  See [internal API limitations](api-implementation.md#ode-series-and-simplification).
  `simplify()` implements ratio selection, custom measures, rationalization,
  inverse combination and `doit` controls. Its candidate search remains smaller
  than SymPy's; abstract-matrix custom controls and full strategy equivalence
  remain open obligations.

* Integral transforms (`laplace_transform`, `inverse_laplace_transform`,
  `fourier_transform`, `sine_transform`, `cosine_transform`,
  `mellin_transform`, `hankel_transform`, their inverses and the held
  `*Transform` classes) follow SymPy's `(result, plane/strip, cond)` tuples
  and the `noconds`, `simplify` and `needeval` hints. The Fourier, sine,
  cosine and Hankel transforms and their inverses are computed through the
  definite integral like SymPy, with SymPy's convergence conditions; native
  closed form tables are only a fallback. Laplace and Mellin transforms try
  their rules first and then the integral (`laplace_transform(cos(sqrt(t)),
  t, s)` is evaluated). Some results are equal to SymPy's but printed
  differently, see [the parity notes](algorithmic-parity.md). Legacy Matrix
  results warn with `DeprecationWarning` instead of
  `SymPyDeprecationWarning`.

Absolute-value derivatives respect real and imaginary assumptions. For an
unconstrained complex argument, derivatives of the real and imaginary parts
remain explicit; the real-variable formula is used only when justified.
NumPy code generation preserves complex-sign semantics on NumPy 1.x and 2.x,
including zero, infinities and NaN. NumPy 2 uses its native fast path; older
versions use a standalone compatibility expression.

`Mod` (with SymPy's full canonicalisation, `%`, `//` and `divmod`), `Rem`,
`sech`, `csch`, `asec`, `acsc`, `asech`, `acsch`, `acoth`, `sinc`, `frac`,
`KroneckerDelta`, `LeviCivita` (`Eijk`) and the constants `Catalan`,
`GoldenRatio` and `TribonacciConstant` are native heads with automatic
evaluation, derivatives, double and multiprecision evalf, assumptions, all
printers, parsing, SymPy conversion, pickling and code generation; C code
uses `fmod` for symbolic `Mod` like SymPy. Sums over a `KroneckerDelta`
linear in the index collapse to a `Piecewise`. `real_root`, `Id`, `Atom`,
`AtomicExpr` and the relational classes (`Relational`/`Rel`, `Equality`,
`Unequality`, `GreaterThan`, `LessThan`, `StrictGreaterThan`,
`StrictLessThan`) are available for `isinstance` checks and construction;
the function heads themselves are callables, not classes. `frac(oo)` currently
stays unevaluated. `AccumBounds`, `Dict`,
`FunctionClass` and `WildFunction` are available; their protocols and limits
are described in the [public interface guide](public-interfaces.md) and
[real analysis guide](real-analysis.md).

The trigonometric, hyperbolic and inverse functions, `atan2`, `exp`, `log`
(including `log(x, b)`), `LambertW` and `sinc` evaluate at special values
like SymPy (`cos(pi/5)`, `tan(5*pi/12)`, `asin(sqrt(2)/2)`, `log(8, 2)`,
`exp(7*I*pi/3)`, ...). `Expr.rewrite(*patterns, rule, deep=True)` is native
and follows SymPy's `Basic.rewrite`: the trigonometric and hyperbolic
functions rewrite to `exp`, each other, `sincos`, `Pow`, `besselj` and
`sinc`, `cos(p*pi/q)` rewrites to radicals (`sqrt`/`pow`, including the
Fermat-prime denominators 3, 5 and 17), inverse functions rewrite to `log`
and to each other, and `exp`, `Pow`, `E`, `re`, `im`, `arg`, `sign`, `Abs`,
`Heaviside`, `DiracDelta`, `Max`, `Min`, `floor`, `ceiling`, `frac`, `Mod`
and `Piecewise` (to `KroneckerDelta`) have their SymPy rules. Supported
rewrites and accumulation-bound behavior remain family-specific.
Special-function, Piecewise, branch and numerical coverage must be assessed
against the [versioned behavioral audit](sympy-parity-audit.md); available
function names do not establish complete rewrite parity.

`Expr.as_real_imag(deep=True)` expands complex parts natively, including sums,
products, principal rational powers and common exponential, logarithmic,
trigonometric and hyperbolic functions. Opaque expressions retain `re`/`im`
projections.
Dense and sparse matrices preserve shape, storage and mutability. `ignore` is
supported for scalar dispatch; other expansion hints are rejected explicitly.
`sec` and `csc` split through `1/cos` and `1/sin`.
`expand_complex(expr, deep=True)`, `expand(expr, complex=True, **hints)` and
`Expr.expand(complex=True)` follow SymPy's hint order: the complex pass runs
first (children before parents, `Abs` folded to `sqrt(re**2 + im**2)`), then
the other enabled hints such as `log` and `mul`. Dense and sparse matrices
expand entrywise; sparse matrices only touch stored entries. `adjoint(expr)`
conjugates scalars and conjugate-transposes matrices. `expand_log` rewrites
`log(b**e)` for real `e` with `-1 < e <= 1`, as SymPy does.
Expansion has conservative resource limits, including 20,000 distinct nodes,
4,096 terms, literal exponent numerators through 128 and a million-bit numeric
budget. Oscillatory
values at infinity that require `AccumBounds` are explicitly unsupported.
`expand_complex` raises `NotImplementedError` for `Derivative`, `Piecewise`,
relational, `Integral`/`Sum`/`Limit` and matrix-expression nodes.
Polar branch constructors and conversions are available separately;
see [complex branches and polar conversions](polar-branches.md).
`re`, `im`, `arg`, `Abs`, `sign` and `conjugate` evaluate natively with
SymPy's rules: the `as_coefficient(I)` split of sums, `as_real_imag`
fallbacks, exact `arg` via `atan2` for numbers, `Abs` of products, powers,
`exp` and `conjugate`, the `sign` product rule, and `conjugate` dispatch
through `exp`, trigonometric, hyperbolic, `gamma`, `erf`, `atan2` and powers.
`Abs(r)**k` with real `r` and odd `k` gives `r**(k-1)*Abs(r)`, `arg(z).diff(x)`
follows SymPy, and `re`/`im`/`conjugate`/`Abs` of matrices work entrywise.
Unconstrained `re(x)`, `im(x)` and `Abs(x)` retain unknown finiteness and
realness while preserving their extended-real ranges. `Abs(x+oo)` remains
symbolic; substituting `x=-oo` gives `nan`. Principal half-integer powers of
a strictly negative real argument have a provable imaginary part; for
example, `im(sqrt(n))` retains `-I*sqrt(n)` for negative `n`. These corrections
do not certify every complex rewrite or branch protocol.

## Boolean logic and sets

`And`, `Or`, `Not`, `Xor`, `Implies`, `Equivalent`, `ITE`, CNF/DNF/NNF
conversion, `SOPform`, `POSform`, `simplify_logic`, truth tables and a native
SAT solver are available. Set operations include intervals, finite sets,
unions, intersections, complements, products, image sets, condition sets,
ranges and the standard named sets through `S`. `solveset` handles polynomial
and supported periodic equations; unresolved equations remain `ConditionSet`s.
Periodic solutions exclude denominator poles.

The SAT solver uses propositional atoms. `use_lra_theory=True`, `minimal=True`
and external solver selections raise `NotImplementedError`. Dense truth-table
operations support at most 20 variables; ordinary SAT solving has no such
limit. `simplify_logic(..., dontcare=...)` minimizes over the care set using
atoms from both expressions; `SOPform` and `POSform` also accept don't-care
rows. `to_anf`, `is_anf`, `ANFform` and `anf_coeffs` support algebraic normal
forms over GF(2). Sparse ANF conversion has no variable-count cap; dense
truth-table conversion has the same 20-variable limit. ANF constants survive
pickling and SymPy conversion. Deep conversion also expands `ITE` into ANF,
where SymPy 1.14 may leave it unchanged. Inequality solving
supports one variable. Other unresolved membership queries remain symbolic.

## Additional structural interfaces

`FU` exposes the public Fu rule dictionary backed by native transformations.
`kroneckersimp` removes products whose Kronecker equality constraints are
incompatible, including nested scalar expressions. `EPath` and `epath` select
and transform native expression arguments and nested Python lists/tuples using
root paths, type/attribute filters, indices and slices. Additive and
multiplicative arguments follow printing order. Other iterable container types
and non-root paths are not currently supported.

`MellinTransformStripError` is available from `symox._transforms` as a
`ValueError` subclass. The native transform table engine remains unchanged.
See [structural sets](symbolic-sets.md) for the additional set and ordinal
exports, and [combinatorics](combinatorics.md) for finite group and ranking APIs.
