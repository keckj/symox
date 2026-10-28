# Review and correction: winding symbolic scales in Gruntz

This follow-up is based on the frozen future-limit package, native SHA58b25015 and gruntz source83912328. Main, that package and every previously delivered stage remain unchanged. The new private stage is `branches/limit-symbolic-winding-guard`; Demeter is `/home/keckj/symox-ci/wt/next-limit-symbolic-winding-guard`.

## Concrete mathematical regression

Let delta(n)=(n+1)^(1+i)-n^(1+i), with real n>0. Baseline cb5 holds the limit of sqrt(delta)/sqrt(-delta). The uncorrected future package58 returns sqrt(1+i)/sqrt(-1-i)=i. Its inverse similarly becomes -i. Both finite answers are false.

Set h(t)=((1+t)^(1+i)-1)/t. The positive real base1+t admits an analytic principal logarithm near t=0, and differentiation proves h(0)=1+i. Exactly delta(n)=n^i*h(1/n).

For integer j tending to infinity take n_j=exp(2*pi*j+pi/4) and m_j=exp(2*pi*j+5*pi/4). Their n^i carriers are exactly exp(i*pi/4) and exp(5*i*pi/4). Thus delta(n_j) tends i*sqrt(2), and delta(m_j) tends -i*sqrt(2). Both limiting values and their negatives avoid the principal square-root cut and zero. Continuity gives ratio limits i and -i. No numerical samples establish this result. The new test checks the derivative, exact carriers and exact principal-root values independently with SymPy1.14.

The uncorrected package also preserves sqrt(delta**2)/delta as held; no regression is claimed for that case. It is included as a noninteger-power case whose total formal scaling degree is integer but whose principal root is not rational in the scale.

## Cause and corrected proof

The old transformation introduces a positive P for w^s0. Positivity is justified when s0 is provably real, since w approaches0 through positive reals. With nonreal s0, P can wind. Homogeneity for positive P (including total degree0) does not prove homogeneity for the actual complex scale.

The check must precede the positive substitution: constructing expressions with positive P can already distribute roots or cancel their P factors. The follow-up first substitutes a NEUTRAL dummy when s0 is not certified real. It accepts that trial only if its scale-dependent operations are Add/Mul and exact integer powers. Powers with scale-dependent exponents and other scale-dependent functions are refused. Scale-independent coefficients/functions remain allowed. This is an exact rational-function certificate, with no numerical threshold or Float-to-integer conversion.

Together with the existing homogeneous-degree check, rationality proves g(P,w)=P^m*g(1,w) for every nonzero complex P on its defined domain. The existing reciprocal/nonintegral-base and leading-coefficient domain checks are retained. Pure rational cancellation for k=i or1+i remains valid: (n+1)^k/delta_k tends1/(k+1). A scale-independent sqrt(1+i) also remains supported. Proven real scales retain the original fractional-power route.

This guards the entire positive-P calculation, including its previous-success arm, rather than merely the new P=1 fallback. A refusal returns to the established limit machinery and preserves the original held calculation when it remains undecidable.

## Scope and remaining gaps

Only gruntz.rs and a new Python test are changed relative to the frozen limit package. Three Rust tests cover the certificate and safe-real/unsafe-complex distinction; eleven new Python tests cover principal-root refusals, the exact two-subsequence proof, certified complex rational cases, integer scale powers, an independent root coefficient and a safe real scale. The original43 tests, including all stored Float payload/precision checks, are unchanged and pass.

Two extra finite examples were probed before being used as gates. The first-power ratio of two shifted differences and a root ratio with a symbolic positive exponent remain held in both58 and the follow-up, despite mathematical limits1/2 and -i. They are preexisting capability gaps, not regressions or proven absence. No test now locks those finite examples to a held result, and no mathematical assertion in an existing test was weakened. The squared/cubed difference ratios already computed by58 remain computed correctly.

The original limit package still needs its documented integration edit: remove the strict xfail decorator only from test_limwrong_parity.py::test_issue_22334_value. This follow-up does not edit that existing file. Its body passes under --runxfail.

## Validation and artifacts

- Local private DEBUG build succeeds, native SHA92ca9eebb92c6e3d0d67cea9edf8260a1bc82865cd6a8863f80bf01611932d22.
- Rust limits filter:11 passed (original8 plus new3).
- Local new+original files:54 passed.
- Demeter private same binary:307 passed, comprising54 new/original homogeneous/winding cases and253 existing limits, sequences, Gamma and Heaviside cases; absolute private PYTHONPATH, existing22334 body included with --runxfail.
- Only gruntz.rs differs among the frozen package's328 native build inputs. Final cohort hashd99f8afbf8e20e166b1f0e4af91a3f31ba9cea2360a4b05cfbb30be5b9badf4e; individual base/final hashes in source-manifest.json.

Apply symbolic-power-winding-followup.patch AFTER the original limit patch12086902. Its SHA is cafb4806a762b53935797da44785928e034d2e3c2991cb5cf13f234cd176bdb3. It passes git apply --check against the frozen owner source. The optional symbolic-power-with-winding-combined.patch contains both packages against immutable02fa and has SHA6bf3361c6a15c8a3982aebd239116b1a5802d8da53d566d6b695c11ecba8382d.

No runtime Python limit algorithm was added. No shared build/install, maturin, commit or performance claim. Parent owns integration and a new combined full-suite gate; the gold r1b suite predates this future limit package.
