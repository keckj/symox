# Rational parameter branches: exact review plan for the r1 integration replay

This preparation concerns precisely the 21 retained-parameter entries identified
in the immutable 570f triage. It is not a new corpus result. The final frozen r1
outputs and conditions must be recorded before assigning the proposed classes.
Expected scope: 20 certified conditional results, one retained affine nonzero
case. The complete 1210-entry replay also retains every other input and result.
No sampled values or quadrature certify convergence, a branch or equivalence.

## Unit interval with one shifted pole

For finite complex c, x+c has its unique zero at x=-c. If -c is outside the
closed real interval [0,1], its segment avoids zero. For Im(c)!=0, its imaginary
part keeps a fixed nonzero sign, so the principal logarithm is analytic along
that segment. For real c>0 it lies on the positive axis. For real c<-1 it lies
on the negative axis; both endpoint logarithms have the same +i*pi term,
which cancels in their difference. Therefore

    Integral(1/(x+c),(x,0,1)) = log(1+c)-log(c)

on Contains(c,Complexes) & ~Contains(-c,Interval(0,1)). Choosing c=a or c=-a
proves both shifted families without logarithm combination. The omitted poles
include both endpoints and the interior. Interior poles are never interpreted
as a Cauchy principal value. Outside the condition the original Integral is the
default, preserving its existing divergence/error behavior.

## Affine slope

At a=0 the integrand is identically one. For finite nonzero a,

    1+a*x = a*(x+1/a),

so the previous segment proof gives

    (log(1+1/a)-log(1/a))/a

under -1/a outside [0,1]. No multiplication identity of principal logarithms is
assumed, and the expression is not replaced by log(1+a)/a. The zero branch must
precede the nonzero branch. For a declared nonzero free symbol Eq(a,0) already
becomes False during ordinary Piecewise construction. Adding this formula would
lose its zero branch under later substitutions, so the existing retained
Integral remains unchanged. This is exactly the one residual capability gap
`int[nz] 1/(1+a*x) (x,0,1)`; no resolution is claimed.

## Segment from a to one

For finite a off the closed negative-real ray (-infinity,0], the straight segment
from a to 1 has no zero and stays in the principal-log slit domain. For nonreal
a its imaginary part preserves its sign except at endpoint 1, which is regular;
for positive real a it stays positive, including the reversed orientation a>1.
The primitive log(x) then gives -log(a). The original Integral remains the
fallback on the closed cut or at nonfinite parameters.

## Quadratic ray

For finite a off (-infinity,0], s=sqrt(a) on the principal branch has Re(s)>0.
The roots +/-i*s are not on the positive real ray. Exact partial fractions give

    1/(x*x+s*s) = (1/(x-i*s)-1/(x+i*s))/(2*i*s).

For x>=0, x-i*s remains in the lower half-plane and x+i*s in the upper one, so
both principal logarithms are continuous and no branch jump occurs. If arg(s)
is theta in (-pi/2,pi/2), arg(i*s)=theta+pi/2 and arg(-i*s)=theta-pi/2.
Consequently log(i*s)-log(-i*s)=i*pi, with equal real logarithmic moduli. At a
large positive endpoint R the difference log(R-i*s)-log(R+i*s) tends to zero.
On compact intervals the denominator has no zero; once x*x>=2*abs(a),
|x*x+a|>=x*x/2, so the tail is absolutely bounded by 2/(x*x). This proves both
convergence and the exact value pi/(2*sqrt(a)). Negative real a has an interior
positive pole, and a=0 an endpoint pole, so these remain in the original path.

## Substitution and branch-construction obligations

The Complexes membership remains symbolic even when assumptions establish
finiteness, keeping NaN and infinities out after substitution. A new rule is
withheld when its cut-exclusion condition is already True for a still-free
parameter, because that guard could disappear before a later contradictory
real specialization. Existing Regular paths predating these rules are not
changed; global correctness of substitutions contradicting symbol assumptions
is not claimed by this scoped review.

Only Piecewise objects whose entire formula, domain and original Integral equal
the canonical certified construction are subject to the hardening. If a
substitution actually changes the original Integral and leaves free parameters,
it returns that transformed Integral, including symbolic chains and partial
numeric specializations. Fully numeric direct substitutions follow the ordinary
conditions. A symbolic Dummy masking a final simultaneous numeric replacement
may conservatively leave the Integral until doit. Foreign Piecewise objects and
explicit whole-node/value/condition replacements retain their own semantics.
Repeated doit uses the exact fixed-point guard and must not self-embed the
fallback or swallow errors in another branch.

The final archive will check the actual 21 output trees, retained domains and
original fallbacks, then record targeted substitution/idempotence controls
separately from the exactly-once complete replay. These controls exercise the
proven rules; they are not sampled proofs of an analytic identity or convergence.
The established fallback computation can retain its pre-existing limitations.
