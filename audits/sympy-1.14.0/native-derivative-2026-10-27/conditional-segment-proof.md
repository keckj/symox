# Certified scope

Only default-method, `conds='piecewise'`, exactly ordered bounds 0..1 and an `Analysis::Unknown` rational-pole result select these rules. Established regular or divergent results never pass through this extension. The integrand must have numerator one and denominator either x+c or 1+a*x, with commutative parameters independent of x. Other bounds, nonlinear denominators, numerators and domains keep the complete original route. For the affine family, a parametric slope whose zero test is already False also retains the original Integral: otherwise normal Piecewise simplification would discard its zero branch, and a later substitution contradicting those assumptions could select an undefined continuation. Constants have no such parameter substitution and remain admissible.

Explicit method switches, including meijerg=False/True, and conds='none'/'separate' retain their full prior path. A differential probe on immutable release 146 and the private extension compares exact recursive expression trees and errors for all 30 combinations of the 3 families with 10 such options.

# Exact domain and branch proof

For denominator x+c, the unique zero is x=-c, so the denominator is nonzero on the integration interval exactly when -c is outside the closed real interval [0,1]. The condition also requires c to be a finite complex number, represented by Contains(c,Complexes). Symbolic membership is retained even when assumptions imply finiteness, to keep nonfinite substitutions guarded.

If Im(c) is nonzero, x+c stays strictly in one open half-plane for real x in [0,1]. It crosses neither zero nor the principal logarithm's negative-real-axis cut. Consequently, the principal logarithm can be differentiated continuously along the entire path, giving derivative 1/(x+c), and the integral is log(1+c)-log(c).

If c is real and -c is outside [0,1], either c>0 or c<-1. In the positive case the real logarithm gives the same formula. In the negative case both endpoint principal logs have the same +i*pi component, while their real components are log(abs(1+c)) and log(abs(c)); that imaginary constant cancels. Differentiation along this real segment still gives1/(x+c). These cases exhaust the certified condition. No principal-value interpretation is used.

For denominator 1+a*x, a=0 gives the identically-one integrand and integral 1. For finite complex nonzero a, the exact identity 1+a*x=a*(x+1/a) gives

    [log(1+1/a)-log(1/a)]/a

under the condition -1/a outside[0,1]. The same horizontal-segment proof applies to c=1/a. The two logs remain separate; no identity involving log(1+a) is assumed.

The default branch is always the complete original Integral, including endpoint and interior poles. At NaN or infinities the finite-complex condition fails; the substituted Integral follows its original evaluation behavior. Substitution of NaN into a held Integral can retain Integral(nan) until doit, whereas constructing Integral(nan) directly immediately gives nan; the tests compare the original substituted held object and then its evaluated result.

# Piecewise doit stability

A deep doit still evaluates every branch and condition. If evaluating one branch produces exactly the entire parent Piecewise again, the original branch is retained at that fixed point. This prevents repeated embedding of the conditional value into its own fallback Integral. No host-error state is cleared. If all children retain their public presentation, the original Piecewise is returned, as in the previous map_children path; thus raw evaluate=False presentation is preserved. Otherwise the existing Piecewise constructor performs its usual rebuilding.

# Evidence limits

Symbolic certificates verify the unique linear pole, exact denominator factorization, primitive derivative and exact output formulas. Exact endpoint/interior substitutions exercise exclusion; regular real and complex substitutions compare principal endpoint logarithms. These checks accompany the branch proof rather than replace it with sampling.

The two historical tests requiring bare Integral for these families are updated to require Piecewise with the exact original Integral default. Unported families still require bare Integral. Existing improper-integral differences from SymPy, including independent one-sided infinities rather than Cauchy principal values, are preserved.

All new mathematical work and the fixed-point guard are Rust. No runtime SymPy and no Python algorithm are introduced. No release performance claim is made from the private DEBUG build.
