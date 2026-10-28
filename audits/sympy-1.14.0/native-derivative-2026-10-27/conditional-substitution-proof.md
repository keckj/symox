# Authorization is the exact certified rule

A Piecewise followed by an Integral is not by itself a mathematical certificate. The recognizer first rejects host objects and held arithmetic, then requires the exact complete value/domain construction of one of the delivered fixed rational-integral rules. It extracts the single original definite Integral, reconstructs that rule using the existing native constructor, and compares the entire resulting expression, including all formulas, domain conditions, zero cases, and the final Integral. It performs no general integration, heuristic equivalence test or host callback.

Thus a user Piecewise with the same finite-complex/interval guard but value a+1 is rejected. So are a different interval, a different original Integral, an extra branch or a changed formula. A mathematically identical canonical rule reconstructed or unpickled is still recognized; there is no process-local origin registry or private predicate.

# Conservative specialization

The complete certified Piecewise represents its original Integral: on its proven regular domain the first branch equals that Integral by the archived branch/convergence proofs, and elsewhere the original Integral is its own default. Returning to this full Integral therefore preserves the represented computation while avoiding a lost guard.

The fallback is allowed only for a pattern involving free parameters of that original Integral, and only if the normal native substitution actually modifies that Integral. A new symbol can be captured by its bound variable, so a symbolic right-hand side also retains the full path even if the resulting Integral has no free parameter. A numeric partial substitution that leaves another free parameter also returns to the complete Integral: otherwise an assumption on the remaining parameter could erase a needed cut guard.

A complete numeric specialization uses ordinary Piecewise substitution and keeps its normal formula/condition behavior. Identity, irrelevant and bound-variable substitutions do not trigger the rule. Exact whole-node replacements keep their earlier priority. Exact maps that explicitly override a value, domain guard, endpoint constant or fallback node retain ordinary exact-replacement semantics; they are not misclassified as consistent parameter specialization. Bound-variable edits and irrelevant keys do not by themselves trigger the fallback.

Both subs and xreplace apply their existing transformation to the original Integral; binder scope, capture, sequential versus simultaneous mapping, validation, and exceptions are not reimplemented. The single-substitution hook dispatch remains before the new rule. Host-bearing source/replacement objects retain their existing protocol route.

# Simultaneous and later computation

There is no new context, marker, cache key or memo change. The existing simultaneous substitution path can internally mask a numerical replacement with symbolic Dummy objects. That can conservatively retain the complete transformed Integral even when the final values are numeric; its value is obtained with doit, exactly through the established integral path. This is a limit of immediate simplification, not a loss of the capabilities of release 146.

After a symbolic replacement, explicit doit remains free to run the complete integral algorithm, including established Regular results under its new assumptions. The hardening does not turn those old Regular results into a new global mechanism for rejecting assumption-contradicting operations. It prevents the newly conditional formulas from losing their default before that complete path has been chosen.

No global Contains behavior is changed, and unrelated Piecewise retain their previous transformations and host errors. Memo, functions.rs and Python product code are unchanged. All new recognition and substitution selection are Rust.

# Evidence

Targeted tests cover all five families, generic and imaginary symbolic replacements, expressions, partial and complete numeric substitutions, sequential/simultaneous maps, algebraic patterns, capture/bound variables, exact replacement priorities, foreign Piecewise and explicit formula/domain/default overrides. They also check repeated doit, reconstruction/pickle, host KeyboardInterrupt, error propagation and operation without runtime SymPy.

Independent fresh-process probes compare foreign/control outcomes exactly with the immutable predecessor and exercise the former false-finite chains against the complete original Integral. A broader corpus includes integration conditions, Piecewise behavior, host protocols, binders, scalar/matrix substitution, precision-sensitive evalf substitution and existing substitution regressions.
