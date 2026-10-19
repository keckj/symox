# Autonomous parity audit worklog

## Function and numerical protocols checkpoint — 2026-10-19

Production fixes: `eae7922a53332cb064c2cf0fec53e5d000889e40`. The rebuilt wheel and all
native downstream runs are identified in `audits/downstream/candidate.json`.
The complete internal Python suite passed 29,063 tests, with
5 documented skips, no failures and no expected failures.
All 383 optional-backend tests passed with required
JAX/PyTorch CUDA, Slang runtime and AdaptiveCpp OpenMP execution. Rust checks
on the same core source tree
passed 505 unit tests and five doctests, with one ignored test;
workspace compilation and correctness Clippy checks passed. Focused checks
passed 819 cases on Python 3.10 and
764 on Python 3.12. Another 1,047 regression cases cover
restricted evaluation and scalar conversion; the audit configuration checks
passed all 14 cases. Exact commands and compressed logs are indexed in
`audits/downstream/validation.json`.

The changes preserve held derivative constructor identity, undefined-function
instance queries, escaped symbol names and Tuple reconstruction options.
Lambda, predicate, differential-field and zero-tensor callability use specific
native subtypes. Ordinary expressions remain non-callable after geometry
imports. Field evaluation survives pickle reconstruction and identity
operations through semantic atom registration. Native reconstruction reads
the loaded facade without requiring the
caller's import hook. Python float/complex conversion preserves adaptive
precision and single-rounding literal conversion, and numeric evaluation
memoizes shared DAG nodes by precision and refinement context. Numerical
callables preserve complex arrays and foreign operand protocols without
silently coercing them through float64; Fraction/Decimal exact arithmetic and
arbitrary module overrides remain explicit open obligations.

NRPyLaTeX now reproduces all 32 reference tests. ChemPy passes 561 of 573
reference-passing cases; its twelve remaining failures require wildcard
predicate support. The eight passing selections are NRPyLaTeX, Haarpy (1,107),
SfePy (221), unyt (819), pyneqsys (19), pyodesys (45), the Spyder console
selection (one), and the yt unit selection (one). Upstream skips and expected
failures remain outside that evidence. Each successful selection reproduces
every passing reference node ID, with installed reference SymPy execution
blocked by the importer. The unqualified full Devito reference retains its
previous candidate identity; no native run was eligible.

TorchInductor's reference reports 236 collected tests and 146,621 successful
subtests under pytest 9. Its reported call-phase pass count includes both;
native collection still fails on unsupported core classes.

Lcapy completes its entire unchanged 317-test selection: 193 pass and 124 fail.
Both phases use native Python tracebacks, avoiding pytest's recursive local
variable comparisons while formatting symbolic failures. Selected upstream
test hashes are unchanged. No assertion was modified to obtain these results.

The featured replacement gate remains false: 11 native integrations
fail and 3 reference suites remain unqualified. Latexexpr exercises no
symbolic integration. Remaining gaps include wildcard predicates, nominal
symbolic subclasses and metadata, assumptions and Boolean hierarchy, physics
state and units, numerical matrix functions and native runtime bridges.
The versioned global parity snapshots below have not been re-certified for
this candidate. Universal equivalence remains uncertified and the user's
Symox >= SymPy requirement remains unmet.

## Historical scalar and permutation checkpoint — 2026-10-19

Production fixes: `2a0acc18867c1a57d09e95a9e72450f8c037b603`.
The final candidate wheel and every native downstream run are identified by
`audits/downstream/candidate.json`. The full Python suite passed 29,012 tests
with five documented skips and no failures or expected failures. Optional
backends passed all 383 tests, with required JAX/PyTorch CUDA, Slang runtime
and AdaptiveCpp OpenMP execution. Rust core checks passed 504 unit tests and
five doctests, with one ignored test; workspace compilation and correctness
Clippy checks passed. Focused checks passed 800 cases on Python 3.10 and
557 on Python 3.12. Exact commands and compressed logs are indexed in
`audits/downstream/validation.json`.

The final upstream rerun passes Haarpy (1,107 tests), SfePy (221), unyt (819),
pyneqsys (19), pyodesys (45), the Spyder console selection (one), and the yt
unit selection (one). Existing upstream skips and expected failures remain
outside that evidence. Upstream assertions are unchanged, every eligible
native run reproduces its passing SymPy 1.14 reference cases, and the importer
blocks execution of installed reference SymPy code.

The featured replacement gate remains false: twelve native integrations fail,
three reference suites are unqualified, and latexexpr does not exercise SymPy.
The report retains all those outcomes. Remaining work includes nominal
Symbol/Function subclass compatibility, assumptions and Boolean hierarchy,
physics state and units protocols, numerical matrix functions, and native
runtime conversion bridges. The historical global parity snapshots below
have not been re-certified for this candidate. Universal parity is still
uncertified; the user's global Symox >= SymPy requirement remains unmet.

## Historical parity audit — 2026-10-18

User requirement: traverse every function in `docs/sympy-api.md`, including
signatures, and attain at least SymPy's mathematical capability. Continue
autonomously overnight; no global certification until the remaining gaps close.

The current reference is SymPy 1.14.0. `catalog.json.gz` contains the complete
binding inventory, constructor signatures and public class-method inspection.
`source-certificates.json` contains 21 scoped source-equivalence checks.
The generated human report is `docs/sympy-parity-audit.md`.

## CURRENT CHECKPOINT — release29 GREEN

**LATEST SUPERSEDING CHECKPOINT — 06:14UTC/08:14Paris.**
Full29 FINISHED GREEN: **28,808passed,59skipped,three xfailed,15warnings**,
511.95s. `/tmp/symox-parity-workspace-batch29-pytest.log`, session19990 completed.
Release focused952passed/two xfailed,7.90s; audit/provenance17passed,2.21s;
core504passed/oneignored/fivedoctests. No build or fullpytest remains running.
Installed optimized backend29 and production Python wrappers are the final
validated package. Current package fingerprint matches catalog and BOTH current
six-object replay archives. Workspace diff whitespace check passes.

Final changes in this continuation:
- Piecewise public constructor class/type queries/fromiter/aliases, empty error,
  held false-branch/default filtering. Nine constructor-specific tests.
- Native held scalar symbol insertion + ordinary Python evalf delegation fix
  premature rounding; 37precision/signature/component/NoRef tests.
- Public solve relational dispatcher preserves all mixed equation/inequality
  conditions and coherent public aliases; 37full-set/parameter/NoRef tests.
- Canonical-identity replay overlays preserve true per-record provenance and
  unrelated stale observations;17auditor/closed-source checks.
- Catalog6077callables/67753bindings/311modules and all scoped records refreshed;
  report/ledger use full24archives plus exact current29 target overlays.

Universal function certification count is STILL ZERO: all6077 flagsFalse.
The user's global Symox>=SymPy1.14.0 condition remains NOT MET. Major next work:
full Risch/heurisch/MeijerG; maxn/strict and special-adapter numerical substitutions;
certified complex root boxes; general/mixed/rational-power solver families and
parameter strata; EX/Gaussian/inexact Poly backends; inherited Basic protocols.
Current doc-replay followups include infinite numerical Sum inN, nestedBoolean
Piecewise conditions and independent-variable reduce_inequalities systems.
Do not call this completed global parity. All mathematical improvements already
proved over reference must be preserved. No source edits, permissions, external
messages, commits/pushes, subagents or image/viewer opening remain pending.

## PREVIOUS CHECKPOINT — release29 installed; full29 running

**LATEST SUPERSEDING CHECKPOINT — 06:10UTC/08:10Paris.**
RELEASE29 INSTALLED atomically after4m02s optimized build. No Rust build now.
Production Python+backend FROZEN while full29 runs, session19990,
`/tmp/symox-parity-workspace-batch29-pytest.log`, ~53% last seen, no failures.
Release29 focused **952passed/two xfailed**,7.90s, plus core504pass/oneignored/
fivedoctests. New37numerical-substitution and37relational-dispatch tests included.
Full28 previously GREEN28734pass,516.57s. Full27 GREEN28725pass,510.84s.

29 also fixes public solve RELATIONAL dispatch. `_solve_relations.py` routes
any direct inequality/Ne in scalar/sequence inputs to existing native
reduce_inequalities, retaining the ENTIRE original equation+inequality list.
Previously mixed[x+y-3,x>3] silently discarded x>3. Pure algebraic inputs keep
the original native solver. Public loaded aliases synchronized, private
_core.solve retained. Strings/nativeTuple/list/tuple/iterator and inferred/
explicit vars/flags tested. Univariate polynomial/rational/Abs full interval
sets compared; parameter equation conjunction proved exactly equivalent by
reference Boolean simplification, not sampled. Private relational sweep341pass.
Remaining native multi-variable inequality restrictions explicitly recorded.

Current29 fingerprint
`a7240761319c254b4ffaaa16e5d0aee100985b66f9cc7cd60ec5ae8c0bf001a8`.
Catalog/all scoped artifacts refreshed. Radical review now5records (adds
Expr.evalf), integration3, source21, permutation16, group8, lattice1. ALL
universal_function_parity_certified FALSE. Source-hash staleness0 current rows.
Complete targeted29 ordinary/keyword runs archived as
`batch29-final-examples.json.gz` and `batch29-final-keywords.json.gz`.
Six canonical objects: roots AND Expr.evalf examples match, four differ;
solve has25difference observations vs30on28, including remaining expression
printing, implicit flag, mixed targets, rational-power equations etc.
Other remaining docdifferences: infinite numerical Sum inN, nestedBoolean
Piecewise condition forms, integration exception-message wording.

Audit tool now supports --extra-examples/--extra-keyword-examples on ledger AND
report: overlays targeted evidence by canonical reference identity, replacing
older aliases correctly while leaving unrelated evidence stale. Mixed headers
explicitly saymixed, and every record keeps its original package/auditor hashes.
17auditor/closed-source tests pass2.21s including overlay provenance/version
checks. These2newtests were added after full29 collection; focused separately.
CSV6077rows universalFalse; six targeted rows are fresh, older24 fullreplay
others stay stale. Aggregate2916objects:1549match/1365diff/2nosuccess; named
1487match/1427diff/2nosuccess. These are observations, never compliancepercent.
Report and ledger refreshed with full24base + targeted29 overlays. No new long
whole-document replay needed. Await full29; investigate any actual failure,
then update final checkpoint with exact final counts. No source/package changes
now required. Global parity condition remains NOT fulfilled.

## PREVIOUS CHECKPOINT — release28 GREEN; release29 compiling

**LATEST SUPERSEDING CHECKPOINT — 06:02UTC/08:02Paris.**
Full27 GREEN: 28,725passed/59skipped/three xfailed,510.84s.
Full28 GREEN: **28,734passed/59skipped/three xfailed**,15warnings,516.57s.
Logs `/tmp/symox-parity-workspace-batch27-pytest.log` and batch28 equivalent.
28 production Piecewise is now the native-backed CLASS facade, synchronized
aliases/type checks/fromiter/unbound methods/class_key/repr. Empty input rejects.
Held evaluation drops false branches, stops at true, collapses one true branch.
New class tests9passed; production focused568passed,55.42s. Global Function/Expr
nominal inheritance and complete Basic/ExprCondPair protocols remain open.

Targeted ordinary/keyword replays on immutable28 complete6distinct callables
(Expr.integrate is deduplicated): roots matches, other five have differences.
Archives `batch28-final-examples.json.gz` and `batch28-final-keywords.json.gz`.
These expose a real precision loss in ordinary evalf(subs=...), not just printing:
x+y-z at x=z=1e16,y=1 yields0. Both Rust substitution AND Python _public_special
adapter eagerly evaluate arithmetic before refinement. Fixed source29:
- evalf_mp::substitute_for_evalf holds changed scalar arithmetic on sequential
  symbol insertion, with DAG memos retaining source refs. Numeric values remain
  compact, no rationalization/dense exponent allocations.
- Ordinary Python Expr.evalf delegates directly to native before eager .subs.
  Special-function nodes or nodes inserted by substitutions keep legacy adapter.
- Bindings, non-symbol algebraic keys, Boolean/Piecewise/set flow and matrices
  embedded as scalar operands retain legacy substitution; no universal claim.
  Scalar matrix ENTRIES do use the corrected native path.
Private29 wrapper focus **605passed/two xfailed**,9.05s.
Final extended substitution tests **37passed**,0.66s, cover15/50/120digits,
1e16/1e100/1e1000 cancellation, nonlinear sin/square, positional subs, Matrix/N/n,
existing Float leaves, complex component cancellation, genuine1e-100 branch-cut
component, fallback semantics/sequential mapping and NoRef subprocesses.
Reference raw evalf may exhaust its default maxn or return scaled zero at tiny
components: exact cancellation identities and independent400digit mpmath oracle
are used. Old failed oracle logs/snapshots retained; no accuracy weakened.
Core29 final504passed/one ignored/five doctests pass. Later only testdiagnostic
was strengthened; final certified-build core check likewise passes.

RELEASE29 BUILD RUNNING session63193, separatetarget/release, started~06:00UTC.
`/tmp/symox-parity-batch29-release.log`. Do not edit Rust during compilation.
Production Python wrapper29 is applied after full28 finished; installed backend
is still28 until atomic install of29 when build finishes. No fullpytest now.
Private metadata `/tmp/symox-parity/debug29-evalf-wrapper-metadata.json`, snapshot
`/tmp/symox-parity/debug29-evalf-wrapper-569a1uc9`. Backend therein is29debug.
After optimized build: atomic install, release focused + full29 on FROZEN package;
refresh catalog and ALL scoped proofs/CSV/report. Radical review now adds an
Expr.evalf scoped record and runtime cancellation probe; do not generate it on
old backend. Integration review includes root __init__ hash/class scope.
Catalog28 binding counts13708pass/48758instance-required/4503gap/298opaque/486data.
All universal certificates remainFALSE. Keep whole24 documentation runs with
truthful stale fingerprints; targeted29 replay will cover the changed evalf.

No subagents, approvals, commits/pushes or image/viewer openings. User goal NOT
complete: all-function capability/signatures/precision/complexity condition is
still blocked by major recorded families. Continue meaningful fixes honestly.

## OLDER CHECKPOINTS — superseded


**LATEST SUPERSEDING CHECKPOINT — 05:44UTC/07:44Paris.**
Installed RELEASE27. Full27 running session1014,
`/tmp/symox-parity-workspace-batch27-pytest.log`, ~15% last seen, no failure.
Package Python/backend FROZEN until it completes. Release27 focused731passed,
15.14s. Private signed-power regression sweep678passed,18.84s.
Core27:502passed,one ignored,five doctests pass.
Full26 FINISHED GREEN **28,689passed,59skipped,three xfailed,15warnings**,519.62s.
Log `/tmp/symox-parity-workspace-batch26-pytest.log`, session66507 completed.
Release26 focused876passed,17.36s. Full25 green28636passed.

27 generalized gamma moment phase degrees to Expr, no polynomial-degree/i64
conversion: positive/negative rational, provably signed-real symbolic and huge
integer q (2**100). Absolute Jacobian |q|; convergence Re(a)>-1 for q>0,
Re(a)<-1 for q<0; rate Re(beta)>0. Whole-line mirror only if q proved even
and p parity known, AFTER convergence. New moment_degree is local to gamma;
ordinary polynomial degree extraction unchanged. Actual conds modes/zero-input
scalar exception/reverse limits preserved. All35new signed-q tests pass,
including NoRef huge/symbolic even q. No complex/sign-unknown q claim.
No Rust build now running; target/release lib_core.so is27, target/debug is27b.

Audit freshness hardened: source hashes declared in scoped reviews are checked
against working tree even when installed .so fingerprint unchanged; new CSV
review_source_stale column. One real edit-without-rebuild regression added.
Audit+closed-source tests15passed,2.28s. Catalog/reviews refreshed27, ledger/report
regeneration started (latest session), using archived complete24 observations.
Scoped counts remain21source,16permutation,8group,1lattice,4radical,3integration;
all universal flagsfalse. Integration review/domain doc updated signed q.

Complete RELEASE24 full replays FINISHED and archived:
`seventh-examples.json.gz`:2916objects,1548matching,1366different,2no-success,
zero incomplete. `fifth-keyword-examples.json.gz`:1486matching,1428different,
2no-success,zero incomplete. Packages predate Fourier fix and25/26/27 features;
keep truthful stale provenance. Sessions27893/85438 completed (may need final
write_stdin poll to close sessions). No complete doc replays now running.

NEXT STAGED WORK: public Piecewise factory FUNCTION -> native-backed constructor
CLASS, enabling isinstance/class queries, constructor fromiter/unbound methods,
correct class_key/repr, empty-constructor TypeError and synchronized module
aliases. Do NOT touch production Python while full27 is running.
Private preview created from exact installed27 package; only __init__.py patched.
Metadata `/tmp/symox-parity/piecewise-class-final-preview.json` has snapshot path
AND EXACT replacement text. Snapshot
`/tmp/symox-parity/piecewise-class-final-preview-9e2j2wpe`.
Existing broad constructor/calculus/printing/series/nfloat sweep **882passed**,
14.02s (`/tmp/symox-parity-piecewise-class-preview.log`, session67584 complete).
New five dedicated class/type/context/alias/NoRef tests **5passed**
(`/tmp/symox-parity-piecewise-class-specific.log`), staged source
`/tmp/symox-parity/test_piecewise_constructor_class.py`, copied in preview tests.
After full27 GREEN: apply EXACT replacement to production __init__.py (replace
current def Piecewise wrapper, preserve _native_piecewise capture); copy staged
five tests to repo; run relevant focused tests and full28 on FROZEN package.
No Rust rebuild needed (Python-only). Preserve inspect.isclass/isinstance,
p.func identity/roundtrip, aliases, fromiter, unbound subs, has/replace class
queries, repr/class_key, evaluate context and empty-args error.
This remains a factory-class facade like existing Symbol; native Expr cannot
be subclassed(pyclass lacks subclass). Do NOT claim full Expr/Function class
hierarchy, is_Function or every Basic/ExprCondPair protocol. No property flags
changed by preview. Broad/source-positive tests passed as noted.
Refresh catalog/scoped integration review to include __init__.py and constructor
class scope, CSV/report after final28. Optionally targeted doc replay of solve,
roots,integrate,Piecewise,N on final immutable package (whole latest24 logs are
stale by design; don't run more redundant full40min replays unnecessarily).

Older checkpoints below are superseded. Current full27 should finish~05:50UTC;
then full28 about9min. User authorized all-night autonomy; no permissions,
subagents, commits/pushes or image/browser opening. Global per-function capacity
condition still NOT met: major Risch/MeijerG, full solver parameter/strategy/
mixed-target coverage, EX/Gaussian/inexact domains, evalf options and complex
root isolation, inherited Basic protocols remain. Never certify them falsely.


**SUPERSEDING CHECKPOINT (05:19UTC/07:19Paris): INSTALLED release25 GREEN.**
Full25 **28,636passed,59skipped,three xfailed,15warnings**,558.51s.
Log `/tmp/symox-parity-workspace-batch25-pytest.log`, session19379 completed.
Release25 focus823passed,17.22s. Catalog and all scoped reviews refreshed25;
radical record now includes parameter rows/coefficient dispatch/Tuple semantics.
Native Python+binary are FROZEN while full runs execute; no full pytest now.

SOURCE26 adds `crates/symox-core/src/solve/recast.rs`, reexported from solve.rs
and used in Py solve for one non-Symbol scalar/system target. Targets are fresh
hygienic Dummies, using reference ALGEBRAIC ATOM conventions (not analytic range
assumptions). Undefined-function derivative/binder coefficients are masked as
separate exact atoms. Existing algebraic substitution handles literal/power
recasts; integer exp templates handle exp(k*u+delta), incl affine proportional
arguments. A reconstruction certificate restores ALL atoms and requires exact
original equality or ordinary ring-expansion zero BEFORE solving; no numeric
residual or forced power/log identity. Work equations must be rational in the
new atom and actually contain it. Numeric/boolean/noncommutative/matrix targets
and matrix-containing equations skip recasting; baseline fallback retained.
This is one-target scalar/system scope, not arbitrary mixed unknowns or ODEs.
Root formulas/QQ indexed fallbacks and flags run inside the existing lexical
scope; restored results/labels preserve original expressions and metadata.

Private26 initial456passed4fail were a REAL improvement over reference:
solve(exp(4*x+2)-5,exp(2*x+1)) returns ±sqrt(5) while pinned reference returns[].
The differential oracle for this one family was replaced by independent exact
quadratic factorization plus explicit exact preimages for BOTH values. This
proves full value-set completeness, not finite-point sampling. No value count
or accuracy requirement was weakened.
FINAL26 focused **462passed**,4.51s; includes all formal expression targets,
Derivative(f(x),x), exp(x), x², exp shifts/inverses, dictionary/set outputs,
NoRef runs, complete cubic/quintic indexed roots with no Dummy leaks, and prior
solver checks/helpers/parameter/coeff/output/component regressions.
Metadata `/tmp/symox-parity/debug26-atoms-metadata.json`, immutable package
`/tmp/symox-parity/debug26-atoms-_ylabp8l`. Refined tests saved separately under
`tests/test_solve_expression_targets_refined.py` in that package; old erroneous
reference-equality oracle kept as evidence. Repo test file is corrected.
Final core26 **502passed,one ignored,five doctests passed**. Two new core tests
prove exact reconstruction and reject noninteger exponential lifting/fractional
work powers. After initial debug only an extra conservative matrix-equation
guard was added; selected scalar behavior unchanged. Core26 final compiled/tested.

**RELEASE26 BUILD RUNNING**, see latest tool session and
`/tmp/symox-parity-batch26-release.log`, started~05:19UTC. DO NOT EDIT RUST while
it compiles. After completion install atomically, release focused expression
and solver regression sweep + full frozen workspace; refresh catalog/reviews/
CSV/report. Add recast.rs hash/scope to radical review and docs first, generate
artifacts ONLY after installing actual26. No general parity claim.

Long ordinary7/keyword5 full doc replays remain RUNNING on immutable24 packages
BEFORE Fourier fix; snapshot provenance must stay truthful. Current public
CSV/report still catalog24+archived22 observations; refresh to25/26 and keep
stale flags. Scoped integration generator's q scope now explicitly i64 bound.
Remaining q rational/negative/signed-real, full Piecewise class/Basic protocols,
full solve multi-target/parameter/strategy coverage, Risch/MeijerG, domains,
evalf maxn/strict and complex root boxes remain major. No commits/pushes/images/
subagents/permission questions. Older checkpoint below is historical.


Current date2026-10-19, last clock04:59UTC (06:59Paris). User explicitly
requests autonomous overnight work over every documented function, with real
signatures, domains/completeness, precision and complexity. Do NOT claim global
parity or completion. No subagents/commits/pushes/viewers are authorized.
No goal tool active. Preserve all unrelated/preexisting worktree changes.

**INSTALLED binary is release24**, with subsequent Python Fourier repair.
Full24 completed: **28,552passed,2failed,59skipped,3xfailed,15warnings**,527.64s.
Both failures were Fourier traversal seeing new ExprCondPair views rather than
native scalar expressions. Python series/fourier.py now flattens whole pairs
ONLY for its scalar traversal/replacement maps. Exact Fourier targeted sweep
**17passed**,3.55s. Original full24 log `/tmp/symox-parity-workspace-batch24-pytest.log`,
completed session16120. Do NOT call original full24 green. Release24 focused
**1130passed,one skipped**,20.34s (before Fourier repair).
Full23:28456passed/one obsolete reference-refusal assertion. The latter was
replaced with independent cubic minimal-polynomial/Sturm/triple-angle proof;
failed-case plus cubic sweep12passed. Full22 was genuinely green28301passed.

SOURCE25 (Rust solver + Python Fourier repair) is ready/validated privately:
- Unresolved parameter-only constraints do not prove NoSolution. Linear zero
  rows reject only exact is_zero==False; substitution similarly preserves
  generic conditional candidates and original-equation verification. Explicit
  nonzero constants/assumption-proved nonzero parameters still reject.
  solve([x-y,y-3],x) now{x:y} like reference. Complete parameter strata remain
  unproved; retaining an undecidable candidate does not prove all its values.
- Scalar multi-unknown solve tries the existing native Python coefficient
  collector ONLY for rational functions of independent Symbols, accepting
  exactly one nonempty coefficient-solution mapping. The collector calls the
  native SYSTEM solver, so no scalar recursion. Other attempts fall back.
  Canonical coefficient pivot order is separate from user's tuple column order.
  Collection runs outside lexical radical scope; nested native calls scope
  their own flags. No reference fallback. Actual scalar vs [system] distinction
  matches reference: solve(a*x+b-2*x-3,a,b) gives{a:2,b:3}; [expr] remains an
  algebraic family. Complete scalar mappings return a dict by default.
- Native Tuple equations/unknowns/initial points now unwrap as sequences in
  expr_list, affecting solve,nsolve and parametric nonlinsolve consistently.
- Core25:500passed,one ignored,five doctests passed.
- Intermediate25a:244passed/3fail solely reversed coefficient pivot convention;
  sorted collector variables fixed those while preserving output columns.
- 25b private broader solver sweep405passed,4.54s.
- FINAL25c+Fourier repair private sweep **823passed**,33.88s, includes solver
  options/checks/helpers, all43rootcontrols, gamma92, Piecewise protocols,
  all275test_expr_remaining (including Fourier), Tuple solve/nsolve cases,
  NoRef sparse coefficient degree10**6 dispatch0.31s.
  Metadata `/tmp/symox-parity/debug25-complete-metadata.json`, immutable package
  `/tmp/symox-parity/debug25-complete-2udoury9`.

**RELEASE25 BUILD RUNNING**, session41538,
`/tmp/symox-parity-batch25-release.log`, started~04:59UTC. DO NOT EDIT RUST
until compilation finishes. Install release atomically via copy to .pending
then os.replace, NEVER overwrite mapped shared library in-place. Target:
`/tmp/claude-1000/-home-keckj-Documents-symox/7ddb00e6-9f68-4045-a47b-087149c85b90/scratchpad/target/release/lib_core.so`
->`python/symox/_core.abi3.so`. After install: focused solver/ctor/gamma/Fourier/
rootcontrols and full workspace on FROZEN Python+binary. Refresh catalog and
all scoped reviews/CSV/report using final fingerprints. No full pytest running
now. Cargo debug build/test SAME target must run sequentially; release uses
separate target. Source25c last debug built, target/debug/lib_core.so current.

Release24 features (must preserve):
- Exact gamma moments c*x**a*exp(delta-beta*x**q) over(0,infinity), currently
  positive integer literal q accepted by monomial parser(i64 bound), finite
  complex parameters, provably Re(beta)>0. Condition Re(a)>-1. Sparse phase
  grouping handles1+I rate, shifted exponent and q10**6 without dense arrays.
  Complex powers with proved convergence are evaluated. A false condition
  never produces finite analytic continuation. Sums use conjunction of
  sufficient conditions plus ONE fallback ORIGINAL integral (no infinity-
  infinity). Whole-line even-q/integer-power scope unchanged.
- REAL conds piecewise/none/separate, invalid modes ValueError, multiple-limit
  separate rejects. Reversed gamma bounds negate value, keep condition.
  tuple/list/native Tuple limits. Finite/indefinite results stay scalar.
- Piecewise public args are tuple-compatible ExprCondPair views with expr,
  cond,args,func,to_sympy. Native func returns PUBLIC Piecewise constructor;
  p.func(*p.args), replacement, pickle work. Native Piecewise accepts list and
  native Tuple branch pairs. Structural/folding/replacement/legacy adapters
  accept old flat and new paired representation. Full Basic/ExprCondPair
  inheritance is NOT certified. Fourier scalar traversal fixed in25 Python.
- Gamma tests initially failed complex phase parse (fixed sparsely), then exact
  reference exp_polar representation (known log(1+I) normalized exactly), then
  public args/refunc protocol (implemented, NOT bypassed). No tolerances were
  weakened. Sources/tests in test_gamma_integral_conditions.py and
  test_piecewise_branch_protocol.py.

Release23 features: exact reciprocal rational-power product conjugacy,
all3 real cubic roots numerically through200digits/pickle; stable quartic
branches and genuine tiny imaginary parts preserved through10^-1000; actual
roots/solve cubic/quartic options, roots strict/predicate/filter; full QQ indexed
fallback. Multivariate product component unions supplied OR expanded, regular
pivot preference and inferred/ordered/unordered solve outputs. Full parameter
strata and comparative complexity are still open. Exact cubic inequality
x³-3x+1>0 now supported where reference refuses; independent WHOLE interval
proof in test_cubic_conjugate_accuracy.py, not finite sampling.

CATALOG24:311modules,67753bindings,6077distinct callables;13644bindingpass,
48854instance-probe-needed,4471bindinggap,298opaque,486data.
Scoped artifacts24:21source-equivalence,16permutation,8group,1lattice,
4radical/inequality,3integration/Piecewise records. ALL universal flagsfalse.
`tools/sympy_integration_review.py` requires actual gamma/pair feature before
recording, uses imported package fingerprint. Integrate scope is POSITIVE
INTEGER q accepted by current parser; don't imply all q/counters. Need regenerate
on25 after Fourier fix so source hashes match final implementation. Likewise
radical review should add25 conditional-candidate/coefficient dispatch scope
and Python _solver_checks.py hash if extending the solve record. Never certify
all solve/Basic/integration behavior from these scopes.
CSV6077rows/report regenerated24 with newest completed RELEASE22 observations.
They are stale relative24 and explicitly fingerprinted.

Full RELEASE22 replays FINISHED and archived:
`sixth-examples.json.gz`:2916objects,1547match,1367differences,2no-success,
zero incomplete. `fourth-keyword-examples.json.gz`:1485match,1429differences,
2no-success,zero incomplete. NativeSHA
a6cabb3d6be292bfe4cd6594c7d76a667493809d9102305d206b282bbea337ef.
Original logs `/tmp/symox-sixth-examples.log` and
`/tmp/symox-fourth-keyword-examples.log`; sessions29434/4730 completed.
Targeted23 solve/roots/N ordinary+keyword archives also exist. solve still
35different observations of98 examples; printing differences AND actual gaps
(non-Symbol targets, undetermined dispatch before25,etc) are distinguished.

NEW FULL IMMUTABLE RELEASE24 replays RUNNING (snapshot BEFORE Fourier repair):
- seventh ordinary session27893, `/tmp/symox-seventh-examples.log`,
  `/tmp/symox-parity/seventh-examples.json`, snapshot release24-ordinary.
- fifth keywords session85438, `/tmp/symox-fifth-keyword-examples.log`,
  `/tmp/symox-parity/fifth-keyword-examples.json`, snapshot release24-keywords.
Do not mutate their package snapshots or attribute25 fixes to their results.
Archive and refresh reports when complete. Audit workers forbid SymPy imports
on native side and desktop/browser launches on both; NO PNG/image opening.

Next implementation priorities under consideration (no edits yet):
1. Single non-Symbol scalar solve unknowns f(x),exp(x),x²: recast to a hygienic
   Dummy, mask derivative coefficients for undefined f(x), recognize exact
   INTEGER exp/power identities, restore originals. Gate to rational work
   equations, preserve baseline fallback. Algebraic-atom convention is reference
   behavior: exp(x) as target CAN return0; real-valued f(x) as target can return
   ±I. Do NOT copy target range assumptions into Dummy. Numeric targets return[]
   in reference, notValueError; don't add a mistaken validation rule.
   Exact reference examples checked: solve(Derivative(f(x),x)-f(x)-x,f(x))
   ->[-x+Derivative(f(x),x)]; solve(x+exp(2*x),exp(x))->±sqrt(-x);
   solve(x**4-1,x**2)->[-1,1]. Existing Python recast_to_symbols uses subs;
   derivative-dependent cases need careful masking, not naive derivative of a
   replaced constant. Restrict initially to single target/arithmetic equations.
   Fresh native Expr::dummy(name,Assumptions::NONE) exists. xreplace map uses
   FxHashMap<Expr,Expr>. inspection::scalar_nodes is iterative DAG enumeration
   but not topological in reverse for arbitrary shared graphs; use proper
   enter/exit states or safe memoized rewrite. Guard noncommutative/numeric/
   boolean/matrix targets. No recast code written yet.
2. Gamma phase q positive/negative rational or provably signed-real Symbol:
   sparse LOCAL moment-degree view returning Expr q (not dense poly_coeffs
   with10,000degree cap), value Gamma((a+1)/q)/(abs(q)*beta**((a+1)/q)),
   convergence Re(a)>-1 for q>0 or Re(a)<-1 for q<0. Whole-line branch rules
   need their own proof; don't generalize negative-axis fractional powers
   blindly. Current positiveinteger family still scoped. No implementation yet.
3. Public Piecewise is still a factory function rather than full Expr class
   facade; full inherited protocols/constructor signature obligation remains.
Major Risch/MeijerG, EX/Gaussian/inexact Poly, certified complex root boxes,
full evalf options/maxn/strict, full parameter strata and other Basic gaps mean
GLOBAL PARITY IS NOT CERTIFIED. Continue autonomously, don't ask permissions.

Both old full immutable batch17 doc runs FINISHED and archived:
- fifth-examples.json.gz:2916 objects,1547 match,1367 differences,2 no-success,
  zero incomplete.
- third-keyword-examples.json.gz:2916 objects,1485 match,1429 differences,
  two no-success,zero incomplete.
Both native snapshot SHA3bff5181ca52250b9df5ad7d601dc65b1032dcf5b1caa0050a4dc250ba74af41.
They predate later constructor/solver fixes; retain stale provenance. Reports
still use older fourth/second evidence; regenerate with fifth/third after
stabilization. Catalog/reviews were refreshed on release21; CSV6077 rows and
report regenerated21 before latest archives completed, all universal flagsfalse.
21 scoped source checks,16 permutation,8 group,1 solver lattice review exist.

Reference authenticity additionally verified against downloaded official PyPI
1.14.0 wheel:1,561 package files and all1,533 Python sources includingisympy.py
match; zero differences. Wheel SHAe091cc3e99d2141a0ba2847328f5479b05d94a6635cb96148ccb3f34671bd8f5.
`reference-wheel-integrity.json`; native reference unchanged. Source hash same
972084bcdf2aeda02be0f83a39ce771094bd7657efa5a0e81ac4d75cad8737a9.
reference_environment() now attaches wheel attestation when hashes agree.
Viewer/browser launching stays blocked; NO PNG/image opening.

Continue autonomously. Main remaining EX/Gaussian/inexact coefficient domains,
full solver/integration/Basic/root-box and evalf option gaps remain large;
no universal equivalence or completion claim, no commits/pushes. Older sections
below are historical and superseded by this checkpoint.

## Immediate checkpoint — batch 21 debug installed, release compiling

CURRENT INSTALLED backend is **batch21 DEBUG**, for focused correctness checks.
Release21 is compiling, session92960, log `/tmp/symox-parity-batch21-release.log`.
Do not edit Rust during it. Install release atomically when complete and then
run full workspace pytest on frozen Python/backend; no full pytest running now.
Debug21 focused constructor+trig+hyperbolic+checksol+PurePoly:
**409 passed, four skipped,10.21s** (`/tmp/symox-parity-batch21-debug-focus.log`).
Debug21 broader old domains/poly/structure/implicit-domain sweep:
**784 passed in7.38s** (`/tmp/symox-parity-batch21-debug-polynomials.log`).
Core21:495 passed, one ignored, five doctests passed.
Previous release20 regression sweep:787 passed, four skipped (11.57s).

Source21 fixes the identified legacy metadata loss using
Poly::from_coefficient_terms: convert only coefficient expressions, retaining
all formal main-variable exponents. Both RF and algebraic legacy views use it.
A sparse polynomial with generator y+1 and degree10**6 now keeps degree/support
and differentiates without expanding the binomial; NoRef subprocess test passes.
Poly.total_degree(0) is now0 as reference; degree/degree_list remain -infinity.
Raw-data unification20 and new21 metadata keep structured coefficient domains,
formal additive/function/numeric generators, sparse support and flags intact.
Ordinary Poly and PurePoly pickle sparse coefficients plus typed domain; PurePoly
explicit reference interchange likewise uses data, not a destructive as_expr
round trip. New constructor comparison helper checks degrees and flags too.

Source21 also fixes real-domain inversion for complex affine arguments.
A target outside real sine range may still have real x when a*x+b is complex;
range rejection now requires proved-real a and b. `restricted_image` solves
Im(center)+n*Im(period)=0 exactly: all indices, no indices, a proved integral
single candidate, or an unresolved exact intersection. Applied to trig,
hyperbolic and exponential lattices. Verified exp(I*x)=1 over Reals gives the
whole 2*pi*n lattice, exp((1+I)*x)=1 gives{0}, and sinh(I*x+1)=sinh(1) retains
all2*pi*n roots. Tiny nonreal exp targets are rejected by exact proof, not
machine epsilon. Hyperbolic real boundary tests also use certified endpoint
signs. New complex affine coset comparisons cover six heads and three slopes.
Unknown zero slope stays ConditionSet instead of unguarded division by a.

After release21 passes full tests: regenerate catalog, all21 scoped source
certificates, permutation/group/rounding reviews, new solver-review JSON,
per-function CSV and report. Add solver review link/report evidence and retain
all scopes/complexity obligations; no universal flags true. Current long full
doc runs still use immutable batch17 snapshots, before several fixes:
fifth ordinary session68210, keyword third47068. Last progress2062/1841 of2916;
archive on completion, keep stale snapshot provenance. Latest fully green
workspace run remains16; full17's two constructor failures are fixed and pass
regressions. No current full green21 claim yet.

## Immediate checkpoint — source batch 20 compiling, next fix identified

Installed backend is batch19, not18. Batch19 focused suite:
942 passed, four skipped, eight failed. All trig lattice/hyperbolic/checksol
cases now PASS, including the previously unsupported complex reciprocal target
1+I. Eight new structured-domain/composite-generator arithmetic tests fail:
round-tripping a whole polynomial through an expression introduces the symbol
inside a formal generator as a spurious coefficient parameter and can lose
main-variable exponents for y+1. Source20 uses coefficient/exponent data directly
for rational-function and algebraic domain unification, with a matching-domain
fast path. It also allows numeric formal generator labels (valid reference
coefficient-list/dict input). Core20: 495 passed, one ignored, five doctests pass.
Release20 compiling, session59173, log `/tmp/symox-parity-batch20-release.log`.
No Rust edits until release20 completes.

Python PurePoly explicit reference interchange and pickling now use sparse
coefficient data instead of as_expr, preserving generator slots such as y+1.
164 focused checks pass (12 structured cases deselected pending20). Ordinary
Poly Python reducer now likewise uses typed sparse data, preserving the old
native fallback for legacy unsupported domain descriptors. New tests include
numeric formal generators, pickle, and cross-domain/generator unification;
need backend20. Existing poly/domain regression sweep after two Python fixes:
567 passed in4.09s on backend18. Solver defaults now expose actual
`domain=Complexes` and inequality `domain=Reals` (not synthetic signatures).

**NEXT KNOWN FIX after release20:** FractionPoly::legacy and
AlgebraicPoly::legacy still read as_expr and therefore give wrong base metadata
for formal composite generators. Reproduced p=Poly([a,1],y+1,domain=ZZ[a]):
p.as_dict() is degree1 but p.degree()==0. Must preserve formal exponent slots
when building the legacy metadata view. Proposed common Poly helper builds
metadata from coefficient-expression conversion plus original main-variable
exponents, without reading the full expression. Then add degree/monomial/zero
metadata checks to all structured constructor tests. Consider using a debug
backend for focused iteration before the next release to avoid repeated4min
builds. Keep any installed backend frozen for a full pytest run.
No full workspace pytest currently running. Latest full green remains batch16;
latest full17 had two regressions already fixed in Python. Refresh artifacts
and run full suite after metadata/round-trip fixes pass focused tests.

A new tools/sympy_solver_review.py records the scoped analytic proof for full
trigonometric lattices and independent corrections to reference's nan/infinite
pseudo-roots. Generate its JSON only after stabilized backend verification.
Two full doc replays still running on immutable batch17 snapshots (both
verified): ordinary fifth session68210 and keyword third47068. Last progress
about1530/1423 objects; archive on completion, retaining known stale fingerprints.
21 conditional source-equivalence proofs are available; no universal parity.

## Latest checkpoint — batch 19 compiling

CURRENT INSTALLED backend: batch 18. SOURCE batch 19 release compiling:
`/tmp/symox-parity-batch19-release.log`, session 40405. Do not edit Rust
until it finishes. Core19: 494 passed, one ignored, five doctests passed.
Latest fully green workspace run remains batch16: 27,972 passed, 55 skipped,
three xfailed, zero failed, 498.67s. Batch17 full run FINISHED:
**28,089 passed, 55 skipped, three xfailed, two failed**, 569.04s.
Log `/tmp/symox-parity-workspace-batch17-pytest.log`. Both were constructor
regressions: native DomainElement coefficients must use exact parent embedding,
and composite-generator renaming over structured domains needed the old direct
expression fallback. Both Python fixes installed; exact two failures plus
constructor tests: 121 passed (`/tmp/symox-parity-batch17-followup.log`).
No full workspace pytest currently running. No current full green claim.

Batch17 public Poly/PurePoly normalization supports keyword gens, coefficient
container/generator inputs, field=True, sort/wrt, unexpanded factors, option
validation/exclusions and nonmutation. Factory constructors retain options.
General EX/Gaussian/inexact coefficient backends and complete greedy/composite/
series option behavior remain gaps. Do not silently claim those options done.
`_poly_constructor.py` is invoked by the Rust public constructor; private
`_from_normalized` performs backend construction. Never monkeypatch tp_new.

Batch18/19 extend solveset periodic families to numeric complex sine/cosine/
tangent targets and cot/sec/csc reductions, exact empty omitted-value cases,
complex affine slopes and optional exact integer-lattice normalization.
Real range bounds use outward rational intervals, no machine epsilon band.
`real_sign` is a new certified endpoint-sign helper in special/lambert.rs.
Installed18 focused trig/hyperbolic/checksol: **190 passed, four skipped,
three failed** in 5.87s. All three failures were cot/sec/csc target 1+I:
exact is_zero alone could not prove it nonzero. SOURCE19 adds exact real/imag
component sign proof. New test file `test_trigonometric_solution_lattices.py`
compares whole lattices modulo integer shifts, not sampled residuals.
Reference 1.14 returns {nan} for sec/csc=0 and infinite pseudo-roots for
 tan/cot=±I; native exact EmptySet is independently justified and tested,
not an attempted reference match. Preserve this distinction in reporting.
Source19 also removes an unnecessary Symbol-only restriction in structural
rename_gens: exponent slots and coefficient descriptors are preserved exactly.
New structured-domain/composite-generator constructor tests need backend19.

Next: install19 atomically, rerun the whole focused trig/constructor/domain
suite, resolve remaining failures, then freeze Python/backend for full pytest.
Refresh catalog/reviews/CSV/report after stabilization.

21 scoped source-equivalence checks now pass (five manually normalized leaves,
16 additional algorithms). `tools/sympy_stdlib_certificates.py` checks actual
LOAD_GLOBAL bindings, shared standard-library objects, recursively equivalent
helpers, local stdlib imports, defaults and absence of ambient-state/closure
dependencies. Tests prove a matching wrapper does not hide a changed helper
or rebound builtin. Default certificates() still returns the original five;
CLI uses extended=True and writes 21. New certificate package has whole-native
fingerprint. Ledger now marks stale source proofs too; generic report wording
no longer hardcodes five. Source/auditor tests: 14 passed in 2.38s.
Snapshot publication is now atomic and validates source hashes; concurrent
snapshot/corruption tests pass. Existing running replays use frozen prior
AUDITOR_SOURCE, and both snapshots were independently hash-verified.

CURRENT LONG AUDITS (both frozen batch17, before the two constructor fixes):
- ordinary fifth: session 68210, `/tmp/symox-parity/fifth-examples.json`,
  `/tmp/symox-fifth-examples.log`.
- keyword third: session 47068, `/tmp/symox-parity/third-keyword-examples.json`,
  `/tmp/symox-third-keyword-examples.log`.
At last checkpoint about 993/933 of 2,916 objects. Let them finish, archive
with original fingerprints, and never attribute batch18/19 fixes to them.

Completed evidence: fourth ordinary 1,534 matches/1,380 differences,
2 no-success, zero incomplete (batch12 snapshot 675982100eeb…). Second
keyword 1,484 matches/1,430 differences, two no-success, zero incomplete,
2,916 objects (batch15 189ebca9121c8f2f6378c40314a0bf45df78af639018a7d628306b18fd4b84be).
Archive `second-keyword-examples.json.gz`. Report now includes keyword counts,
source integrity and separate snapshot hashes. CSV still needs regeneration.
Catalog is batch17: 311 modules, 67,753 bindings, 6,077 distinct callables;
13,644 binding passes, 48,854 instance probes needed, 4,471 binding gaps,
298 opaque, 486 data. Source reviews/JSON may be stale after backend changes.

Reference: SymPy1.14.0; all 1,533 Python sources verified against RECORD,
zero mismatches/unhashed. Python3.10.12, SymPy GROUND_TYPES=python,
mpmath1.3.0 BACKEND=python. SHA:
972084bcdf2aeda02be0f83a39ce771094bd7657efa5a0e81ac4d75cad8737a9.
Viewer/browser spawning stays blocked; do not open old PNGs or GUI programs.

Group SLPs already work with exact reconstruction witnesses and repeated bases.
Group+SLP tests: 23 passed, all pairs through degree4, random groups, S30.
Precision/huge decimal/mpf fixes are installed and covered by prior focused
503-test suite. Root/integration/solve/Basic domain gaps remain substantial.
No universal flags true, no global equivalence/completion claim, no commits or
pushes. Preserve unrelated preexisting user edits. Older history below is
superseded by this checkpoint.

## Running traversal

The initial complete documentation traversal writes an atomic checkpoint after
every object to `/tmp/symox-parity/nightly.json`, with progress in
`/tmp/symox-parity-nightly.log`. It uses an immutable native package snapshot,
whose fingerprint is stored per result. Corrections made after it started are
not covered by that initial snapshot and need a second pass.

Resume or reproduce with `tools/sympy_parity.py examples`; use `--pending-only`
only when intentionally retaining explicitly stale evidence. `--module` and
`--retry` support targeted rechecks. The native worker blocks SymPy imports,
desktop viewers and browsers. Compiler/TeX subprocesses remain available.
Matching output is finite evidence, not a general mathematical proof; different
output may be a representation difference and requires review.

## First correction batch

- 67 Python functions now accept the reference parameter names as keywords.
  These are real parameter/body changes, not synthetic `__signature__` fields.
- `Integer(i=...)`, `Float(num=...)`, `sqrt(arg=...)`, `cbrt(arg=...)`,
  `log(arg=..., base=...)`, `N(x=...)` and `root(arg,n,k,evaluate)` are supported.
- `Float` retains the inferred precision of long decimal strings and integers,
  and supports automatic exact digit counts with an empty precision string.
- Numeric sets enumerate exact integers/rationals in reference order. Native
  Tuple nodes support iteration, length, indexing and slicing; ordinary scalar
  expressions still reject iteration. Sequence-detection callers were adjusted.
- Public expression constructors expose inherited native methods. Several
  documented import paths, nested finite-set inputs, `Dict.has`, diagonal
  matrix lengths and Dirac delta linearity checks were added.
- Scalar `nsolve` accepts the inferred-variable two-argument form and native
  multiprecision secant iteration with one or two real/complex initial points.
- Logic rebuilding recognizes integer 0/1 introduced by substitutions. Numeric
  truth handling follows the distinct connective conventions.
- DiracDelta/Heaviside constructors use certified sign/real assumptions.
- MP numeric-recognition atoms use a proposed candidate only after checking it
  against the original value in native multiprecision arithmetic.
- `group(..., multiple=False)` now counts one group at a time, preserving the
  reference memory behavior. Prefix/suffix slice indices match reference calls
  even for custom sequence protocols.

The release backend was built and installed atomically. Initial complete pytest
run: 26,932 passed, 55 skipped, three xfailed, two failed. Both failures were
MP `nsimplify` recognition and were subsequently fixed; all 195 related checks
passed. After the native rebuild, 229 focused tests passed. A separate 69-check
batch covers new interfaces and audit safeguards. Four native MP solver tests
passed, including the new secant algorithm. The full suite needs rerunning on
the final frozen package after further corrections.

SymPy 1.14.0 misreads the decimal position in explicit-precision Float strings
containing underscores. Symox parses their mathematical decimal value; the
regression checks compare a normalized reference literal instead of copying
that upstream bug. This is recorded separately from reference equality.

## Second correction batch and regression follow-up

The first documentation replay is complete: 2,916 documented objects and
22,547 examples. Its archived evidence is `initial-examples.json.gz`.
It contains 1,468 matching object replays, 1,440 differing observations, six
incomplete replays and two replays with no successful examples. Representation
differences and cascaded exceptions remain observations requiring review.

A second documentation replay is running in
`/tmp/symox-parity/final-examples.json`, with native snapshot
`b7c0697e7ac4f79b816336f0131ea8ac2cedd79e7602db3b48c40e90174589fd`.
That immutable snapshot predates the sympification/permutation follow-up below;
do not attribute its observations to later fixes. Progress is logged in
`/tmp/symox-parity-final-examples.log`.

Further implemented changes:

- Exact floor/ceiling/truncation no longer pass through an i64 or machine float.
  Outward rational enclosures handle supported constants, powers and exponentials;
  unresolved cases retain existing fallback behavior. Generic numerical
  truncation fallback remains heuristic, not an interval certificate.
- Five independently justified rounding improvements over SymPy 1.14.0 are
  recorded separately in `rounding-improvements.json`.
- Arbitrary-size integer bit operations remain integer operations; symbolic
  Boolean operators enforce the reference operand domains.
- 115 function heads bind named arguments through actual runtime contracts,
  reject invalid arities, and preserve omitted optional arguments. Some native
  named forms extend reference classes that only accept positional calls at
  runtime despite displaying named signatures. This is not evidence that every
  function domain or strategy is implemented.
- Public `sympify` implements locals, xor control, strict mode, rational string
  parsing, evaluation control and converter protocols. The full-suite follow-up
  restored unit/algebraic-wrapper conversion priority, handles native callable
  objects without module names, and preserves foreign symbolic validation errors.
- Textual reconstruction understands normalized matrix slices, matrix name
  atoms and applied binary predicates. Permutations now accept positional cycle
  notation, cycle lists, identity degrees, chained cycles and sequence action.
  Overlapping cycle-list composition uses inverse positions, taking linear time
  in the degree and total number of cycle entries.
- General `simplify` now implements explicit ratio, custom measure,
  rationalization, inverse and doit controls. Its candidate search is still
  smaller than SymPy's and is not certified equivalent. Abstract-matrix custom
  controls remain an explicit gap.
- Circle equation construction, degenerate three-point circles and Curve
  parameter substitution have differential regression coverage.

The second full suite on snapshot b7 completed with 26,732 passes and 455
failures. The failures clustered in the new sympification conversion and textual
reconstruction paths. The corrected regression batch passed 1,571 cases with
71 failures; the remaining permutation/foreign-fact failures were then fixed.
The final focused batch passed all 413 cases. A fresh complete workspace suite
is running in `/tmp/symox-parity-workspace-pytest.log` (backend installed before
the subsequent binary-splitting/expand-func changes).

The refreshed signature catalog has 311 module paths, 67,753 bindings and 6,077
distinct reference callables: 11,119 passing binding probes, 3,650 binding gaps,
52,202 cases requiring instance probes, 296 opaque signatures and 486 data
bindings. These are inspection categories, not mathematical coverage percentages.

Native work after the workspace suite started (not installed yet): exact
arctangent summation for the certified pi enclosure now uses integer binary
splitting and an alternating-tail bound. Three native tests passed, including
exact partial-sum identities and an independent Chudnovsky comparison up to
2,048 bits. `expand(func=True)` is being connected to the existing native
special-function expansion with depth and arithmetic-option controls.

## Next validation

## Overnight follow-up — 2026-10-19

The second complete documentation pass is archived as `second-examples.json.gz`:
2,916 objects, 22,547 examples, 1,491 matching replays, 1,421 differing
observations, two incomplete replays and two without successful examples. Its
b7 snapshot is explicitly older than subsequent corrections.

The first full workspace regression after conversion fixes completed with
27,253 passes, seven failures, 55 skips and three xfails in 561.99 seconds.
The remaining failures were fixed: preserve sparse containers before general
sympification; preserve original keyword arguments through the differential
geometry adapter; reject nonscalar charpoly callback results explicitly; update
the old test that expected the now-implemented simplify ratio option to fail.
All 263 focused follow-up checks passed.

Batch 4 release installed binary-split certified pi enclosures and the func
expansion hint. Native validation: 477 passes, one ignored test, five passing
doctests; 219 focused Python expansion/rounding checks passed. The pi comparison
includes independent Chudnovsky values through 2,048 bits. Timings before/after
for 50/500/2,000 decimal digits are in `/tmp/symox-pi-enclosure-{before,after}.json`;
these few observations are not an asymptotic certificate.

Batch 5 release installed Boolean-versus-arithmetic equality decisions and
preservation of evaluate=False when exporting held relationals. Both new native
Boolean tests and all 138 focused Python constructor/domain checks passed.

56 Python methods/functions in 23 files now use the reference keyword names.
All public aliases were inspected before renaming; conflicting shared methods
were excluded. Constructor classes for immutable/base explicit matrices and
symbolic matrices now expose actual unbound instance methods, without delegating
class factories that would change mutability. Relational constructors accept
lhs/rhs and reference option behavior; Pow accepts b/e; Inverse/Trace/Determinant
accept mat; OneMatrix(m,n,evaluate) no longer swaps dimensions in keyword calls.

Finite permutation algorithms now implement the previously missing own-class
operations: ranks/unranks in three orders, successors, inversions, invariant
queries, Josephus permutations, cycle structure, distances/matrices and group
operations. Factorial-sized ranks remain Python integers. Fenwick order
statistics replace quadratic searches; power evaluation shifts cycles directly.
Josephus uses a deque for small steps and order statistics for larger steps.
Rank/parity/cycle-structure caches preserve repeated-call complexity.

Exhaustive checks cover all 874 permutations of degrees zero through six.
Large tests cover degree-2,000 rank round trips, a 100,000-element inversion
count and exponents with 2,001 decimal digits. All 438 permutation/construction,
matrix-permutation, combinatorics and Galois integration tests passed. Empty
successor handling and out-of-range Trotter-Johnson unranking intentionally fix
reference failure/invalid-permutation cases; these are scoped improvements,
not matching-reference observations.

The latest catalog has 13,216 passing binding probes, 49,498 instance probes,
4,255 binding gaps, 298 opaque signatures and 486 data bindings. More methods
can now be inspected on constructor classes; newly visible gaps do not mean
previously working computations regressed. The module/binding/canonical counts
remain 311 / 67,753 / 6,077.

A fresh full workspace suite is running in
`/tmp/symox-parity-workspace-final-pytest.log` (exec session 9618). Do not edit
native Python sources or replace the backend until it completes.
The third immutable documentation traversal is running in
`/tmp/symox-parity/third-examples.json` and
`/tmp/symox-parity-third-examples.log` (exec session 55363). The package is frozen
independently of later workspace edits. The report and 6,077-row CSV were
refreshed from the newest catalog and explicitly older second-pass evidence.

## Remaining work

## Modular arithmetic and keyword-dispatch follow-up

The next full workspace suite completed with 27,395 passes, three failures,
55 skips and three xfails in 541.34 seconds. All three failures came from
specialized adapters changing a wrapped operation's mandatory first keyword.
`_dispatch_arguments.FirstArgument` now preserves actual binding through vector,
indexed tensor, array and differential-geometry adapters and the special-function
derivative adapter. All six specialized import orders and named scalar/container
calls are tested; 298 integration checks and the 13-check import/control batch
passed.

Derivative is now a native-node constructor class accepting expr/*variables,
with reference instance checks and unbound methods. Calculus constructor
classes expose native operations; Subs accepts the assumptions keywords that
the reference constructor ignores. The 24-check constructor/import batch passed.
Explicit evaluated derivatives with simplify=False remain unsupported.

The modular-combinatorial audit exposed two timeouts, both now fixed in the
installed batch 6 backend. Held factorials use grouped prime valuations and
Wilson complements. Held binomials route to exact modular arithmetic instead
of materializing huge integers. Large prime-power moduli use Granville
interpolation, with extra p-adic precision and the p=2 sign correction, rather
than a loop over p**e. BSD reference provenance remains in the source header.
The modular integer wrapper no longer truncates nonintegral inputs silently.
An unrepresentable integer binomial loop bound no longer returns the false
answer 1; symbolic constructors preserve that enormous output unevaluated.

Native batch 6: 484 tests passed, one ignored, five passing doctests. Tests force
the Granville path against independent small exact unit products, including
exceptional primes/exponents, then compare huge cases. All 639 focused Python
number-theory/constructor/import checks passed. The two formerly timing-out
documentation objects now match every observation: Mod in 0.966 seconds and
binomial in 0.737 seconds. Evidence is `/tmp/symox-parity/modular-examples.json`.
Its native execution blocks reference imports.

The latest full workspace suite is running with the installed batch 6 backend:
`/tmp/symox-parity-workspace-batch6-pytest.log`, exec session 11562. Native Python
sources have been held stable during that run. New test files added after its
collection will be checked separately. Do not replace the backend until done.

Batch 7 is built but NOT INSTALLED: `/tmp/symox-parity-batch7-release.log`, exec
session 37594. It removes the 10,000-decimal/40,000-bit policy caps and permits
one-bit Float/numerical fields and zero decimal digits. A representability
guard reserves doubled precision and guard bits within the u32 native precision
counter. New tests are in `test_large_precision_parity.py`; not run yet.
After installation, update `test_numeric_expression_domains.py` which still
expects precision=1 and 40001 to fail. The remaining architectural counter
limit and other evalf option/algorithm gaps still prohibit global certification.
Two exploratory reference probes hung at invalid Float precision=0/negative
bit precision; only those owned subprocesses were terminated. Do not reproduce
the reference hang or accept zero bit precision.

Sixteen scoped finite-permutation invariant/precision/complexity reviews are
recorded in `permutation-review.json`, generated by
`tools/sympy_permutation_review.py`. The function CSV now retains review scope,
source fingerprints, stale-review flags and stale-example flags. It never turns
these closed-domain reviews into universal Python-protocol certificates. Six
audit guard tests passed. Production edits since the review invalidate its
whole-package fingerprint; regenerate after the next stable snapshot.

The third complete replay remains active in `/tmp/symox-parity/third-examples.json`,
snapshot da0372926e783df8717e61dd0774e9ab5071a1c9a26dc5a03b193b18bc886fbc.
It predates batch 6. Keep its per-result fingerprints; later fixes are not
covered by it. Viewer/browser launching remains blocked.

Next implementation prepared outside production: native group order/membership
currently enumerate at most 100,000 elements. A correct incremental Schreier
chain prototype is `/tmp/symox-parity/schreier_incremental_prototype.py` (imports
raw composition helpers from `schreier_prototype.py`). It uses persistent orbit
transversals and a queue of new Schreier pairs; it only augments levels between
the originating Schreier level+1 and the sift failure, avoiding repeated global
rebuilds. It computes the documented degree-24 Rubik subgroup order 3,674,160
in ~0.002s, and symmetric orders through degree 50 (1.38s at 50), with exact
Python integers. See `/tmp/symox-schreier-incremental-prototype.log`. Not integrated
or fully verified yet. Validate against independent small-group enumeration and
reference membership before changing production. Keep asymptotic certification
pending until the exact implementation bound has been reviewed.

Additional permutation doc gaps to fix after the current suite: integer action
`i ^ p` (__rxor__), rank shifts p+integer and p-integer; support Cycle generators
and empty/dups options in PermutationGroup; make_perm is missing. Permutation
power's direct method keyword is exponent but reference uses n. General inherited
Basic protocols and PurePoly construction/operations are still major gaps.

## 2026-10-19, 01:20 Europe/Paris — batches 7–10

This section supersedes the pending descriptions above. Batch 9 was built and
installed atomically. Numeric fields now support one-bit precision and more
than 10,000 decimal digits. Exact numeric transfer uses binary mantissas and
Python integers, avoiding decimal string limits. Literal complex-part
extraction and field comparison/hash no longer reject high precision values.
The native core passed 485 tests, one ignored, plus five doctests. The focused
Python suite passed 540 tests, including 20,000-digit precision and roundtrips.

The group prototype is integrated in `combinatorics/_schreier.py` and
`perm_groups.py`. Group order and membership use incremental stabilizer chains
instead of enumerating at most 100,000 elements. Lazy generation, orbits,
stabilizers, supplied bases, derived and lower central series are implemented.
Single-generator and named-family fast paths avoid unnecessary chains. All
generator pairs through degree four were checked against independent closure;
random groups through degree nine, S20 and the documented Rubik subgroup were
checked against the reference. The broad algorithmic complexity review and
remaining group APIs are still pending; `slp_dict=True` is unsupported.

The third ordinary documentation replay is complete and archived as
`third-examples.json.gz`: 1,530 matching observations, 1,382 differences, two
incomplete and two without successful examples. Its snapshot predates the
recent fixes. Mod and binomial passed separate replay after their fixes.

Full workspace batch 9 completed: 27,456 passed, 55 skipped, three xfailed,
five failed in 472.76 seconds. Four failures revealed that binary transfer of
Float zero used a reference tuple constructor which canonicalizes to Integer
zero. The fifth test still expected the removed precision cap. Both causes
are being corrected in batch 10; do not report batch 9 as green.

The auditor now supports `examples --keywords`, converting actual reference
documentation calls to accepted named forms, preserving positional-only and
variadic binding. It records per-call coverage and unmatched reference calls.
Eight audit guard tests pass. The full keyword replay runs from frozen native
snapshot `4181ea96259d4168...`, exec session 47438, with output at
`/tmp/symox-parity/keyword-examples.json` and log
`/tmp/symox-keyword-examples.log`. Viewer/browser launching remains blocked.
It exposes genuine runtime keyword gaps behind permissive wrapper signatures.
The snapshot and auditor are frozen, so subsequent fixes are not covered.

Batch 10 is compiling (session 71785, log
`/tmp/symox-parity-batch10-release.log`). It fixes Float-zero interchange,
adds substitution argument binding for arg1/arg2 to native scalar/matrix/poly
types, and accepts row/col keywords in dense matrix deletion. The Python Poly
override is corrected separately, including simultaneous substitution. New
regressions are in `test_substitution_keyword_parity.py`. These changes have
not yet passed their focused suite or been installed. The nfloat precision cap
is removed too. Do not edit Rust while the release build is active.

The current catalog was regenerated after batch 9: 13,561 binding probes pass,
49,002 require instance probes, 4,406 have binding gaps, 298 are opaque, and
486 are data. These are bindings including aliases, not distinct algorithms.
The permutation review was refreshed, but report/CSV need regeneration against
the latest catalog and archived examples. Universal certification remains false.

## 2026-10-19, 01:35 Europe/Paris — batches 10–13

Batch 10 was installed and its expanded focused suite passed 584 tests.
The five full-suite batch-9 failures were fixed: Float-zero binary interchange
now preserves Float and precision, and obsolete precision-cap expectations
were updated against the actual reference behavior. nfloat accepts negative,
zero and more than 10,000 decimal digits. Named/mixed substitution forms,
Poly's simultaneous substitution, sparse iterator materialization and dense
row/col deletion keywords have regression coverage.

Batch 12 is the CURRENT INSTALLED backend and current frozen Python source.
It includes Poly(rep=...), subclassable Poly, domain conversion aliases for
the differing reference real/complex domain names, unify(K1=...,symbols=...),
and corrected EXRAW priority. reduced(f=...,G=...) and
is_zero_dimensional(F=...) accept their reference names. Additional Python
keyword repairs cover det(matexpr), RegularPolygon(c,r), UnevaluatedExpr(arg),
PermutationGroup.contains(g), jacobi_normalized(a,b), Node.attr_params(looking_for).

PurePoly now extends native Poly and retains the class after its algorithms.
Exact-domain arithmetic, unification by generator positions, division/PRS/GCD,
calculus, evaluation, factories, pickle, functional generators and reference-
import-free execution have differential coverage. Existing charpoly data
is retained as `_CharPolyData` for unsupported coefficient backends; do not
claim full PurePoly/domain parity. PurePoly hashing now respects equality
between exact coefficient domains (unlike the reference's domain-dependent
hash), while nontransitive cross-inexact/cross-class reference comparisons
remain an explicit limitation. An old test asserting unequal ZZ/QQ pure
polynomials was corrected against the reference.

CSE now flattens/restores concrete matrices, including immutability. Sparse
inputs operate only on stored entries. A matrix with shape (10**9,10**9) and
two entries is tested without densification. Focused Python batch 12 passed
526 tests; matrix/PurePoly/CSE focused batch passed 202 tests. The complete
workspace suite is RUNNING and production Python/backend must stay frozen:
`/tmp/symox-parity-workspace-batch12-pytest.log`, session 81339 (about 57%).
Do not attribute a final outcome before it completes.

The new complete ordinary replay is RUNNING from immutable snapshot
675982100eeb5967316c934ba7de26cd3527fdd963d090ac3eb0848ebd220479,
`/tmp/symox-parity/fourth-examples.json`, log `/tmp/symox-fourth-examples.log`,
session 24107 (478 objects at last check). The earlier full keyword replay
continues from snapshot 4181ea..., session 47438 (1,743 objects at last check).
Both preserve frozen auditor versions and block viewers and native reference
imports. The CURRENT auditor additionally rewrites imports inside literal
exec/eval snippets; nine auditor guard tests pass. Existing running replays
predate that harness fix and remain explicitly versioned.

The batch-12 catalog has 67,753 bindings / 6,077 distinct reference callables:
13,644 binding probes pass, 48,854 need instance probes, 4,471 binding gaps,
298 opaque signatures, 486 data. Exposing more actual PurePoly methods also
exposes real method signature gaps; do not interpret these counts as a
capability score. `group-review.json` contains eight finite-group invariant
reviews, explicitly retaining pending comparative complexity certification.

Batch 13 RUST SOURCE adds collect Wild-factor matching and exact=None
generator discovery. It is NOT INSTALLED. Its complete native core tests
passed (see `/tmp/symox-parity-batch13-core.log`, session 37061); release build
still runs in `/tmp/symox-parity-batch13-release.log`, session 23045. No Rust
edits during this release compile. New Python regressions, not yet run, are
`tests/test_collect_extended_patterns.py`. After building, test batch 13 via
a separate package copy if workspace batch-12 pytest is still running. Do not
replace its mapped backend. Two old collect tests still expect Wild and
exact=None to fail; change them only after the frozen batch-12 suite ends.

Next: finish both audits; verify batch 13 and fix its real mismatches; refresh
catalog/reviews/report/CSV at a stable installed fingerprint. Report generator
now links scoped permutation/group reviews. Universal parity remains false.

## 2026-10-19, 02:25 Europe/Paris — completed replays, decimal conversion, batch 15

Both long replays completed without timeouts. The fourth ordinary pass is
archived as `fourth-examples.json.gz`: 2,916 objects, 1,534 matching observations,
1,380 differences to review, two without successful examples. Its snapshot is
675982100eeb5967316c934ba7de26cd3527fdd963d090ac3eb0848ebd220479 (batch 12).
The first complete runtime-keyword pass is archived as
`keyword-examples.json.gz`: 1,408 matching observations, 1,506 differences,
two without successful examples. It attempted 69,903 reference named calls,
with 69,785 successful reference calls; native attempted 68,573, with 67,582
successes, and 1,330 reference calls were unmatched. Cascading semantic
differences affect these counts; they are not a signature compliance score.
The keyword snapshot is 4181ea... (batch 9), before many named-call fixes.
Viewer/browser launches stayed blocked. Both frozen auditors predate the
literal exec-import harness fix. Keep individual source/auditor fingerprints.

The full batch-12 workspace suite is GREEN: 27,598 passed, 55 skipped, three
xfailed, 15 warnings in 544.58 seconds. The installed batch-13 focused suite
passed 405 tests; the prior staged test passed 368 with two obsolete unsupported
expectations deselected, then those expectations were replaced by reference
comparisons. collect Wild and exact=None docs now match all 31 observations.
The targeted pass including PurePoly/sympify/CSE completed 134 objects:
111 matching, 23 differences. CSE has one remaining output-order difference.

CURRENT INSTALLED BACKEND is batch 14. Python sources include subsequent
batch-15 fixes, so there is not yet a green full run for this combination.
Batch 14 introduces `bigfloat/decimal.rs`: outward positive dyadic intervals
with i128 temporary exponents, binary power/scaling, adaptive endpoint-round
agreement, and an exact-size path for rounding ties. Binary rounding uses
nearest-even; decimal display uses the reference half-away convention.
No exponent-sized integer is built merely to parse/print a compact Float.
The old one-million decimal exponent constructor policy cap is removed;
the representable native exponent/precision counters remain limitations.
Float parsing+printing at exponent ±10**6 and ±10**12 takes ~0.00003–0.00010s.

Literal eval_mp now rounds finite real/complex numeric literals directly,
without transcendental argument-reduction magnitude budgets. Rational output
still retains the pre-existing numeric engine's allocation range: base10
conversion through exponent magnitude 10**6 and equivalent binary bounds in
domain to_rational. The smaller 100,000-bit/decimal exponent rationalization
caps are removed. Decimal's integer-ratio protocol avoids Python's 4,300-digit
string-to-integer limit. Exact ratio conversion no longer allocates a precision-
sized denominator limit when limit=False. Native integer Float inputs now
avoid int->decimal string conversion; beyond-Decimal-range scientific strings
use the reference default-precision fallback unless exact inference was asked.

Batch-14 focused numerical suite: 467 passed, one failed in 49.83s. The one
failure was the pre-existing fixed/scientific display threshold; source batch
15 now follows reference min_fixed=min(-(dps//3),-5), max_fixed=dps, with strict
inequalities. Four large-integer/exponent Float input tests pass already with
the installed backend and updated Python. Full decimal tests need batch 15.

A real counter-overflow bug was reproduced: with v=Float('1e999999999999999999',20),
v*v fits but v*v*v wrapped the signed exponent and gave a tiny false Float.
SOURCE batch 15 adds checked wide-exponent arithmetic and a cached expression
feature bit. Representable numeric/complex results are computed, and results
outside the counter remain exact held arithmetic. evalf does not evaluate those
held arithmetic nodes through unchecked exponent arithmetic. This preserves
mathematical values but is not full unbounded-exponent Float parity.

Integer Float powers now use adaptive outward dyadic intervals and wide
temporary exponents too. Positive binary rounding ties become exact at enough
working precision; non-dyadic reciprocals cannot be exact binary ties. This
removes the former million-decimal output magnitude cutoff for integer powers.
Independent exact-rational comparisons cover signs, exponents, powers and
precisions. Numeric division now rounds a quotient directly, rather than a
rounded reciprocal times the numerator. A separate performance investigation
found the remaining ~44s in the million-exponent rational test was two ~22s
native Rational constructors; Number::inv now uses a known-coprime raw ratio
for integers/rationals instead of redundant general GCD normalization. This
speed change has not yet been measured in Python (needs the new backend).

Batch-15 full Rust suite: 492 passed, one ignored, five doctests pass.
Release build is RUNNING, session 4321, log
`/tmp/symox-parity-batch15-release.log`. Do not edit Rust during it.
Once complete install atomically, then run
`test_decimal_conversion_parity.py`, `test_float_exponent_arithmetic.py`,
`test_nsimplify_rational_options.py` and the prior numerical focused suite.
The new Python counter-boundary regressions have NOT run against batch 15 yet.
Check the held-expression comparison expectations; do not weaken numeric
correctness to force matches. Then run full workspace tests on stable sources.

Codegen ReplaceOptim callable-query behavior was fixed: it now delegates to
Expr.replace, passing the matched expression and visiting atoms as reference.
Builtin exp2/log1p replacement callbacks were adapted to the proper arity.
Nine focused callback tests pass, including reference-import blocking.
The targeted codegen doc pass had 11 matches and six differences (17 objects);
it predates the atom-query refinement and should be replayed after stabilization.

The CSV generator now accepts --keyword-examples, with separate runtime counts,
evidence identifiers, snapshot fingerprints and stale flags. Nine audit-tool
tests pass. The last generated public CSV/report still use catalog batch 12 and
third ordinary evidence; regenerate with fourth+keyword archives after installing
batch 15 and refreshing catalog and reviews. Source certificates remain five,
permutation invariant reviews sixteen, group reviews eight; all universal flags
must remain false. No global equivalence or completion claim is justified.

Outstanding implementation priorities include full Poly/PurePoly constructors
(keyword gens, dict/list input, option semantics and coefficient backends),
PermutationGroup SLP output, Basic subclass/protocol completeness and the
documented solver/integration/root-isolation gaps. DO NOT monkeypatch native
Poly.__new__: a subprocess experiment showed PyO3 rejects saved tp_new after
replacing that slot ("not safe, use object.__new__"). No production constructor
adapter was implemented. A Rust constructor hook into a Python normalization
helper is a possible future approach. Preserve all unrelated worktree edits.

## Continuing queue

Finish every documentation replay, distinguish semantic gaps from printing
differences, and refresh the catalog/report with native fingerprints. Recheck
the changed surfaces against a final immutable package. Continue fixing open
algorithm/option/domain gaps; `docs/algorithmic-parity.md` lists major solver,
integration, numerical-evaluation and root-isolation limitations. Never turn a
passing name/signature probe or a few examples into a universal certificate.

Existing unrelated/user changes in the worktree must be preserved. Do not
commit, push or publish without authorization.

## 2026-10-19 — optional backends, complex regressions and English documentation

The custom PyTorch, JAX, Slang, SYCL/AdaptiveCpp and CUDA APIs now live in
`python/symox/extra`, with execution tests in `tests/extra`. Genuine compatible
entry points remain available, including root `lambdify(..., modules='jax')`
and `printing.pytorch.TorchPrinter` / `torch_code`. Optional runtime libraries
remain unloaded when importing Symox or generating source. Benchmarks, CI and
code-generation examples use the new paths. Core and extension CI suites run
separately against the pinned SymPy 1.14.0 reference.

The development Python 3.10 environment now contains PyTorch 2.14.1+cu130,
JAX/JAXlib 0.6.2 with CUDA 12 plugins, and SlangPy 0.43.1. NVIDIA RTX 5090
execution works for PyTorch, JAX, native CUDA and Slang CUDA/Vulkan. AdaptiveCpp
commit `31db2e83903c1db9f50b03b620eb0a09224cd0eb` is installed with its OpenMP CPU
profile. `tools/install_adaptivecpp.py` reproduces that pinned compiler build;
it does not establish AdaptiveCpp GPU coverage. `pip check` reports no broken
requirements. Project extras `extra` and `test-extra` declare the Python
execution dependencies.

The first complete core run exposed 33 remaining dynamic test lookups of
old root extension names and one f2py build failure with setuptools 84. Those
lookups now use `symox.extra`. Default f2py builds select Meson when installed
with NumPy >=1.26, retaining explicit flags and the legacy path otherwise.
The compiler subprocess discovers scripts from the active Python environment.
Meson 1.12.1 and Ninja 1.13.2 are installed, declared in the test dependencies,
and documented; PyTorch's setuptools requirement remains satisfied. All 60
focused source-export, dummy-identity and compiled-wrapper tests pass. The
complete core and extension suites were restarted on this final frozen source.

Both remaining complex-evaluation xfails are ordinary passing tests. Projection
range facts no longer imply finiteness for an unconstrained argument:
`Abs(x+oo)` remains held and substitution of `x=-oo` yields `nan`. Negative
real half-integer powers have an exact imaginary-domain proof, and `im`
retains the reference public product form. Infinite numeric terms print as
constant terms. Additional regressions cover assumptions, substitutions and
positive/negative half-integer exponents. The focused check passed 1,332 tests.

Installing the optional dependencies exposed two further issues. Older JAX
requires its JIT device option to preserve CPU placement of differentiated
executables when a GPU is the default backend. Gradients, Hessians, vmap,
constant outputs and transfers between devices now pass. The standalone Torch
source test keeps its import guard around generated execution and restores
normal imports before invoking Torch's assertion helpers, which lazily import
NumPy through DTensor in the installed Torch version.

The optimized extension suite passed **383 tests, zero skips, zero failures**
in 98.95 seconds on the final frozen source. This includes actual AdaptiveCpp compilation/execution and
Slang CUDA/Vulkan execution. Rust passed **504 unit tests and 5 doctests**, with
one existing ignored unit test. The final complete core Python run passed
**28,799 tests, zero failures and zero xfails** in 533.22 seconds. Its five
skips comprise four invalid SymPy nan/infinite pseudo-root oracle cases
(with separate independent checks) and one pytest-style comparison module
collected separately. Combined core/extension validation passed 29,182 tests.

All 60 tracked Markdown files were checked for English prose, current extension
paths and local links. The two remaining French guides are translated. Outdated
capability descriptions were corrected against current implementations, and the
103 Python documentation blocks compile. Historical work-log counts remain
attached to their original snapshots.

The upstream Sphinx inventory incorrectly labels the `Arch` class as a module.
Its importable module is `sympy.physics.continuum_mechanics.arch`, and Symox
exports `Arch` there. The correction is recorded in the inventory metadata;
there are now **310 measured module paths**, 301 with measurable exports and
9 without, with zero native or reference import errors. No algorithmic gap was
resolved by hiding an import failure.

Installed native fingerprint:
`61be54459967a4a22a7e2ed49b804023f02f1d79adbeabc251dac8630ed19c0a`.
The catalog, function ledger, scoped review fingerprints and generated report
were refreshed. The 54 scoped reviews retain their scope limitations; all
6,077 universal function-parity flags remain false. Fifteen complex-function
and Arch documentation objects were replayed in ordinary and keyword modes:
11 match, four retain explicit observation differences (Abs runtime class,
polar expansion, principal-branch simplification and sign numerical evaluation).
Other historical example evidence remains explicitly stale. This checkpoint
is not global algorithmic parity certification.
