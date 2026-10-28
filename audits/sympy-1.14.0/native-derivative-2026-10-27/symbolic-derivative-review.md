# Symbolic Derivative orders: bounded native representation tranche

Private stage: `/home/keckj/.cache/symox-validation/native-algebra-20261026/branches/symbolic-derivative`.
Demeter: `/home/keckj/symox-ci/wt/next-symbolic-derivative`.

This tranche does **not** close issue22210. It adds exact held data and a tested set of native operations. Symbolic-order evaluation formulas and several transformations remain open. Main and the previously delivered stages were not edited.

## Representation and invariants

- The existing `Kind::Derivative(Expr, Vec<(Expr,u32)>)` and literal machine evaluator remain unchanged. The new `Func::SymbolicDerivative` is publicly named `Derivative`; its native arguments are `body, Tuple(variable, Expr count)...`. Public `func`, `args`, `variable_count`, `_wrt_variables`, `derivative_count`, free symbols, printing and DAG pickle retain this shape.
- Counts are never truncated or inferred integral from a Float. Float(0), Float(1), exact noninteger counts, unconstrained/assumed symbols and arbitrarily large integers remain distinct. Integer zero removes a group. Adjacent groups add exact native Numbers/Expr. A specialization to ordinary u32 orders returns to the established machine representation/evaluator.
- Commutation uses the existing dependency graph. Derivative variables rank their exact numeric order, including rational orders, rather than an integer approximation. Counts contribute to free symbols and substitution, but are constant metadata for formal differentiation. Commutativity comes from the body, ignoring a noncommutative count.
- Evaluated construction, count substitution, xreplace, formal differentiation, doit and canonical rebuilding use typed native `construct -> Result`. Negative orders raise; nonreal/nan/Boolean evaluation errors are distinct from accepted held construction. Python pair counts equal to zero are removed before sympification, distinguishing raw False/0.0 from BooleanFalse/Float(0) expressions. The positional count route remains distinct. Integer subclass conversion and pair equality use the reference protocols and preserve their errors.
- Count errors use an RAII stack bounded to the native call. Nested caught calls consume their own frame; panic unwinding and host interruption pop frames. Host exceptions have priority and are consumed before a count error is returned. A NaN cancellation barrier is used **only** within a live frame; an unguarded native substitution retains the entire unevaluated Subs operation instead. Tests cover two invalid terms that would otherwise cancel, repeated errors followed by a healthy call, nested host calls and KeyboardInterrupt.
- Dense Matrix and MatrixExpr diff/xreplace/doit adapters carry the same guard. Matrix subs already uses the checked scalar adapter. The ordinary machine paths are retained.

## Deliberately held / refused transformations

The new operator is opaque to generic map_children, simptools map/rebuild, rewrite and numerical rebuilding. This prevents a generic transformation from partially changing a count and silently turning an error into NaN or zero. Series and complex expansion refuse the new operator where the established machine derivative route already has a refusal/fallback. This fence applies **only to the new head**.

For example, `c=(a+1)**2-a**2-2*a-2` is accepted as a held count by the reference and by this tranche. Expanding c proves -1. Reference `Derivative(sin(x),(x,c)).expand()` raises ValueError; the new operator intentionally remains completely held. Likewise, its body/count are not generically rewritten or numerically approximated. This is an explicit remaining capability gap, not a claim of parity. After substituting a supported machine count, the established machine transformations apply normally.

Trusted generic Func construction is structural only. The approved public construction/substitution routes perform typed validation; scalar return adaptation validates any surviving new Derivative tree. Generic walkers hold the original operator before reaching a changed Func reconstruction.

## Covered behavior

The new regression module covers exact counts and assumptions; 0/integer/negative/noninteger/nonreal/nan/Boolean substitutions; exact Float merging and precision; grouping, commutation and canonical forms; applied undefined functions; elementary valid jet variables; body/count free symbols and commutativity; differentiation inside sums/products/matrices; guarded nonfinite-context substitutions; head/args; str/srepr/latex/pretty for the certified cases; count_ops; CSE; pickle; textual parsing; explicit codegen refusal; and error isolation/priority.

Two old test expectations that declared orders greater than u32 invalid are removed. Positive regressions now require exact held acceptance of 2**32 with the complete count preserved. No tolerance for a mathematically different answer was added to an old test.

## Remaining gaps for TODO / parity contract

1. **Symbol and constant n-th-order formulas:** evaluation/doit of `Derivative(x,(x,n))` and nonzero constants should produce the reference Piecewise formulas; this tranche holds the operator. Their derived forms also remain held when reference specializes them. The 4,536-case differential audit over 6 assumption sets, 7 bodies, 12 orders and 9 operations reports **291 differences, all in x/constant evaluation/doit/diff formulas**. Full records are in `probe-order-differences.json`.
2. **Large literal evaluation formulas:** counts greater than u32 are represented exactly, but no closed formula for e.g. the 2**32-th derivative of sin is added. `variables` cannot materialize these large new counts and explicitly refuses them; derivative_count and variable_count remain exact. Existing machine-order behavior is unchanged.
3. **Generic transformations of the new operator:** expand/cancel/rewrite/evalf and other generic rebuilding currently hold the complete operator or explicitly refuse it, rather than transforming its body/count. The negative polynomial count example above records the intentional difference. Series/complex expansion are not claimed equivalent.
4. **simplify=False symbolic orders:** the preexisting unsimplified-diff adapter explicitly rejects nonmachine counts; no Python derivative algorithm was added. evaluate=False held construction is supported.
5. **Special variable protocols:** new symbolic-order variables requiring host `_diff_wrt` or specialized Symbol identity use NotImplementedError rather than assuming the protocol true. A jet whose body is a host Function also reports NotImplementedError for its unsupported `_diff_wrt` protocol. MatrixElement jets and other unported non-AppliedUndef jet-variable classes are not certified; the ordinary machine route remains unchanged.
6. **Printer error parity outside the certified cases:** str/srepr/pickle preserve held Boolean/nonreal data, but latex/pretty do not yet reproduce every Boolean/nonreal-order failure of SymPy. Do not claim general printing-error equivalence. The certified symbolic/rational/Float forms in the tests match exactly.
7. **Existing matrix coefficient evaluation:** ordinary `D*A` construction may retain a held derivative where reference MatMul immediately calls doit (and can error). Held MatMul(diff/doit) propagation is covered; this tranche does not globally change the established matrix constructor strategy.

## Final private validation

- After the successful core suite, only the Python binding file and new tests changed; all core source hashes stayed unchanged. Focused and broad Python were rerun against the final extension.
- DEBUG cargo build succeeds; full core suite: **643 passed, 1 preexisting ignored**, plus **5 doctests passed**.
- Focused Python: **423 passed**, including **167 new tests**, old differentiation/parser tests and host hooks.
- Broad existing Python: **4217 passed, 3 skipped**, one existing SymPy aesara deprecation warning. Host hooks are included in both focused and broad runs; the union is **4548 unique passing Python tests**.
- Differential audit: **4536 cases**, with the **291 explicitly open Symbol/constant-formula differences** described above. No runtime SymPy fallback was added.
- These are correctness checks in the private copy, not benchmark timings or the final combined-main full suite.

## Integration

Use `other-sources.patch`, then review/apply `lib-rs.patch` and `functions-rs.patch` separately. `symbolic-derivative.patch` is the complete alternative patch; `source-no-lib.patch` includes functions.rs. Do not overwrite the whole main lib.rs: the parent has added ratint registration elsewhere. Do not overwrite functions.rs: the conditional-integrals agent has a separate Piecewise doit hunk. This tranche changes no Piecewise arm. No new module registration is needed. Parent changes in assumptions.rs, simptools/radsimp.rs and integration_methods.rs are outside this patch.

`source-manifest.json` contains per-file base and final SHA256 hashes; base hashes are authoritative because the base main worktree contains prior uncommitted lots. `runtime.json` links the native binary, source digest, private import paths, Python and reference versions. Final validation logs are `cargo-build-final.log`, `cargo-test-final-all.log`, `pytest-final-focused.log`, `pytest-broad.log` and `probe-orders-final.log`. Earlier exploratory logs are not final validation evidence.
