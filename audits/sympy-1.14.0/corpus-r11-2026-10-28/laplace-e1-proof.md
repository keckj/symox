# Principal E1 in the existing affine Laplace rules

Base: immutable r1b gold, native cb5c0c98572e4b05eadca9f234fb6df3bbabea93a202e4a2ac292f5541b5748e,
source manifest identity db859db2e5949ea295c362ae9b06b1b7819f0f28da0753790edda29b45445d2f.
Only the two existing Ei-producing affine rules change their representation.
Their eligibility guards, abscissa zero, conditions, program rules and all
fallbacks remain unchanged. No new Python mathematical logic or runtime SymPy.

Let beta=c/a>0 be the existing reciprocal-rule guard. Then a!=0 and the positive
real time ray has no denominator zero. For Re(s)>0,

    F(s)=Integral(exp(-s*t)/(a*t+c),(t,0,oo))
        = exp(beta*s)/a * Integral(exp(-s*u)/u,(u,beta,oo)).

For positive real s, the substitution v=s*u identifies the last integral as
E1(beta*s). Both sides are analytic on Re(s)>0: uniform absolute convergence on
compact subsets follows from exp(-delta*t)/(|a|*(t+beta)), delta>0, and principal
E1 is analytic there because beta maps this half-plane to itself. The identity
theorem therefore proves F(s)=exp(beta*s)*E1(beta*s)/a on the entire domain.
The argument proves the same formula with a complex or negative scalar a when
beta is positive. It does not assert a value for unguarded negative offsets.

For the existing logarithm rule, a,c>0. Its real-time integrand log(a*t+c) is
unambiguous. Integration by parts for Re(s)>0 has vanishing upper endpoint
exp(-s*t)*log(a*t+c) and lower endpoint log(c), giving

    G(s)=Integral(exp(-s*t)*log(a*t+c),(t,0,oo))
        = log(c)/s + a*F(s)/s
        = (log(c)+exp(beta*s)*E1(beta*s))/s.

Compact uniform bounds on logarithmic growth establish analyticity and justify
the same frequency differentiation and shifts already used by the engine.
Multiplication by exp(-lambda*t) replaces s by s+lambda and the abscissa by
-Re(lambda); on that half-plane Re(beta*(s+lambda))>0, so the same principal E1
branch remains valid. Time moments use (-d/ds)^n with no branch change.

The old -Ei(-z) shortcut agrees with E1(z) only at positive real z using the
real cut-average Ei value. For Re(z)>0, Im(z)>0, Log(-z)=Log(z)-i*pi and the
convergent series imply Ei(-z)=-E1(z)-i*pi; the opposite half-plane gives +i*pi.
Consequently that shortcut added +/-i*pi*exp(z) (and the appropriate scalar or
1/s factor). Replacing it with principal E1 removes this discontinuity. No
combination or forced rewriting of principal logarithms is used.

E1 is constructed as expint(1,z) by the existing native special-function kernel.
Its evaluation does not rewrite it into Ei based on variable assumptions, so a
later complex specialization of an initially positive/real frequency retains
the analytic expression. No branch-changing simplify operation is introduced.
The existing arbitrary-precision E1 path handles numeric evaluation.

Validation is prepared around exact symbolic formula/ROC equality, the ODE

    a*F'(s)-c*F(s)+1/s = 0,

and the integration-by-parts certificate s*G(s)-log(c)-a*F(s)=0. The analytic
integral representation above fixes the branch and normalization; derivative
identities alone are not the proof. Independent high-precision E1 evaluation
controls check requested precision in both complex half-planes. Such controls
are not used to infer a branch, convergence or equality from samples.

The known incorrect SymPy1.14 reciprocal Laplace table is not an oracle for the
changed result; tests use its independently defined principal expint function
and exact differential identities. General Laplace branches and unsupported
coefficient domains are outside this bounded correction and retain their paths.
