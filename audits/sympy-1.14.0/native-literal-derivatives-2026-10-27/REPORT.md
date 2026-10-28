# Literal BigInt derivative orders — private follow-up to r1b

Base: immutable `/home/keckj/.cache/symox-validation/native-derivative-r1b-20261027/after` (release extension cb5c0c98572e4b05eadca9f234fb6df3bbabea93a202e4a2ac292f5541b5748e). This package does not include or overwrite the subsequent LRT/Berkowitz/limit/E1/inexact lots. No main checkout or previously delivered package was modified.

## Exact scope and proofs

The optional Rust certificate accepts literal nonnegative Int/BigInt orders without a u32 conversion. Its variables are ordinary commutative Symbols with no host/scalar identity and no proven infinity. Its coefficients and affine intercepts are exact closed finite numeric expressions; floats, unknown closed functions, free coefficient parameters, host values and nonfinite expressions retain the established fallback. `n=0` retains the original expression. Float/Boolean/negative/noninteger/undecidable counts retain the existing parser and held/error behavior.

* For affine u=A*x+B, the chain rule gives D^n(exp(u))=A^n*exp(u); sin/cos cycle modulo4, sinh/cosh modulo2. These identities hold on C with A,B finite, including canonical sinh/cosh images of sin(I*x)/cos(I*x). A=0 is the exact constant case. Slopes ±1, ±I and E keep their constant-size exact ordinary images (parity/modulo4 and exp(n)).
* For u affine and literal m>=0, induction gives D^n(u^m)=(m)_n*A^n*u^(m-n) for n<=m, and zero for n>m. For n=m-1, (m)_(m-1)=m!, so the same factorial is correct at both top orders.
* Structural polynomial degree upper bounds use max for sums, sums for products and multiplication by literal nonnegative power exponents. They prove zero only above that bound. At the bound d, D^d(P)=d!*[x^d]P; cancellation of leading coefficients is harmless. The bound is obtained before coefficient powers, so enormous factored monomials do not expand just to prove zero.
* Unsupported lower derivatives of nontrivial polynomial products, other functions or variable groups decline the certificate entirely. Node/depth budgets decline this optional optimization; they do not claim absence or zero and preserve the old solver.

All derivative decisions and recipes are native Rust. Python changes only invoke native capabilities or adapt reconstruction; there is no added Python differentiation algorithm and no runtime SymPy mathematics.

## Representation and complete arithmetic

Large coefficient powers may be stored as an exact Pow with a local immutable recipe flag. The flag is absent from kind, public args, equality, structural hash and printing, but serialized in a dedicated private matrix-DAG opcode. A separate inherited feature bit locates recipes in parents. There is no TLS, global evaluation switch or new ordinary power bound. The complete `num_pow` function is byte-identical to r1b (SHA256 in direct-power-proof.json); its pre-existing fractional-radical MAX_POWER_BITS policy is unchanged.

The representation estimate applies only when producing/reconstructing these new derivative recipes. It cannot prevent an ordinary direct Pow from materializing. Numeric coefficients small enough for the recipe estimate are evaluated, including factorial/falling-factorial coefficients around4096; huge ones remain exact Factorial/FF/Pow expressions. The public factorial/FF evaluators and their existing limitations are unchanged. A held Factorial/FF can be explicitly reconstructed through its existing public constructor; the pre-existing behavior of generic held Factorial.doit is not a newly implemented capability here.

For powers, complete direct Pow, int(recipe), recipe.doit(deep=True/False), and public func(*args) are tested on 2**524289 (a reasonable65KiB integer). A tiny Rust test independently distinguishes an explicitly tagged2**12 from ordinary4096 and verifies full materialization. Giant powers are never allocated just for validation. Public func(*args) deliberately follows the ordinary explicit constructor; it need not preserve the recipe flag.

Mul/Add/negation, substitutions, xreplace, pickle/copy, native generic reconstruction, expand and CSE retain the local representation. The flag is authenticated on pickle restore: only exact finite numeric integer powers acquire it, with a nonzero base required for negative exponents. Trivial bases0/±1, exponents0/1, Float substitutions and zero poles use ordinary arithmetic. Zero absorption is certified for finite coefficients; polynomial roots remain regular even with held factorial and slope powers.

Public CSE already reconstructs through Python func(*args). Its added adaptation calls the native reconstructing method only for a genuinely marked Pow, under authenticated bootstrap stock constructors/function code/dependencies/descriptors and ordinary scalar operands. Constructor, func, code, metaclass, native-node and matrix/host dispatch customization use the established route. Regression tests compare actual callback sequences with the established forced fallback, rather than assuming a hook is applicable. Ordinary CSE expressions retain their existing route.

Explicit polynomial operations such as cancel/simplify may still request complete coefficient materialization through their pre-existing evaluators. This lot does not make those arbitrary algorithms compact, suppress their protocols or add a new refusal.

## Compact traversal

Several validation-only walks used the machine Derivative.args form, which repeats each variable by its count. They now inspect compact scalar children (body plus variable/count tuples). The presence/domain checks retain their meaning; count_ops instead multiplies the original per-variable contribution by the exact stored count. Included paths: special-function substitution validation, binders, Lambda/Limit validation, matrix scalar validation, Order-variable validation and evalf's has_dirac_delta inspection. Public args/variables APIs retain their established representation. A later symbolic-order specialization to2**32-1 is explicitly tested.

## Still open

Issue22210 is not closed. Symbolic-order x/constant formulae need a publicly reviewed typed order guard: mathematical membership in Naturals0 would wrongly accept substituted Float(0/1), although the Derivative API retains those counts. The earlier proposal remains archived in the preceding private stage. Nonliteral orders and unsupported large-order functions/coefficient domains remain held or explicitly fail according to the existing route. Matrix.diff's separate count parser still rejects orders beyond its accepted machine range; scalar entrywise differentiation is covered. Full factorial materialization beyond the existing special-function evaluator's own limits remains a pre-existing gap, not a claimed gain.

## Validation provenance

Final hashes and exact gates are recorded separately. Intermediate exploratory runs with a global power bound were discarded; no such bound is present in this package. The debug/release distinction matters for the solve_triangulated_extension timing probe; only equal-build comparisons can establish a performance regression. New tests compare entire expressions/factor multisets, exact periodic classes, degree certificates, substitution behavior, protocols and serialization, rather than only sample residuals.
