# Additional inherited scope error found while reviewing the old triage

This is a validation finding; no product code or frozen snapshot is modified.
The old triage calls all ten 1/(t+1) Laplace differences intentional corrections
on Re(s)>0. The exp(s) prefactor is correct, but extending the displayed
-Ei(-s) principal-branch expression to that entire half-plane needs care.

For Re(s)>0, the integral converges absolutely and defines an analytic F(s).
Writing u=t+1 gives F(s)=exp(s)*E1(s), initially by a positive-real substitution;
E1 is analytic on the right half-plane, so analytic continuation establishes
this identity there. It does not establish that -Ei(-s), with its principal
branch and its cut-average real values, is analytic on the same half-plane.

The exact convergent expansions on their slit domains are

    E1(s) = -gamma-Log(s)-sum_{k>=1} (-s)^k/(k*k!),
    Ei(z) = gamma+Log(z)+sum_{k>=1} z^k/(k*k!).

For Re(s)>0 and Im(s)>0, Arg(-s)=Arg(s)-pi, so
Log(-s)=Log(s)-i*pi. Hence Ei(-s)=-E1(s)-i*pi, and the recorded native result
-exp(s)*Ei(-s) differs from F(s) by +i*pi*exp(s). For Im(s)<0 the difference is
-i*pi*exp(s). At real s>0 the real Ei value on the negative cut equals -E1(s),
and the native formula is correct. These statements follow from exact branch
arguments and power series, not numerical comparisons or quadrature.

Thus six entries (positive, real and nonzero-real symbol setups, each with and
without transform metadata) retain the established real-domain prefactor
correction. Four entries (unconstrained and complex setups) expose a native
complex-branch gap as well as the reference's erroneous exp(-s) prefactor.
The correct analytic expression on Re(s)>0 is exp(s)*expint(1,s). The old archive
is left untouched; this narrower classification is recorded only in the new
27 October triage. It does not claim a general audit of transform branches.
