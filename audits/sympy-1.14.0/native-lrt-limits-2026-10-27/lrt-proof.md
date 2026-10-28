# Stock QQ/ZZ Lazard–Rioboo–Trager reduction

The native selection accepts only canonical univariate Poly objects over ZZ/QQ,
a canonical commutative Symbol/Dummy t distinct from the main generator, proper
nonzero numerator, positive-degree denominator, degree at most 32 and rational
coefficient numerator/denominator sizes at most 4096 bits. These are dispatch
bounds: every refusal runs the existing complete rationaltools algorithm. No
custom object is converted during eligibility. Initial public Poly construction
and `t = t or Dummy('t')` occur before selection exactly as before.

The kernel forms b=f-t*g' in QQ[x,t] and calls the existing checked Brown–Traub
subresultant PRS, retaining its actual scalar resultant and all representatives.
For m=deg(g), the leading coefficient in t of R=Res_x(g,b) is
(-1)^m Res(g,g'). Thus degree_t(R)=m certifies squarefreeness in characteristic
zero. The constant coefficient is lc(g)^(m-1-deg(f))*Res(g,f), so its nonvanishing
certifies coprimality. This reuses the calculated resultant and avoids two extra
rational Euclidean chains. Invalid inputs are declined; their previous results
and exceptions remain available, including the established empty result for
ratint_logpart(1,x**2,x), whose inputs violate the documented squarefree premise.
No claim of mathematical correctness is made for such nonpreconditioned results.

The result is the same LRT construction as SymPy 1.14 and the old Symox port:

* infer the scalar resultant coefficient domain from its exact coefficients;
* use the checked squarefree view with its ZZ versus QQ monic conventions;
* include the complete negative scalar unit in the first group, then take each
  primitive part; preserve group order by multiplicity;
* use the same last PRS entry of the requested degree;
* divide the same powers of monic QQ gcds of LC factors with q;
* compute the same reduced Euclidean inverse of the original LC modulo q;
* keep the leading x coefficient equal to one and reduce the other, nonzero
  coefficients modulo q in the same descending x-degree order;
* infer the final coefficient ring QQ[t]/ZZ[t] (or plain QQ/ZZ when t disappears).
  If the group multiplicity is m, retain the original g and its original domain.

The LC squarefree scalar sign need not be explicitly attached to the first LC
factor: its only use is a monic gcd over QQ, which is invariant under a nonzero
rational unit. All other signed group normalizations are preserved explicitly.
The returned h,q representatives—not just their equivalence classes or the
sum's derivative—are compared to both implementations. Consequently subsequent
root identification uses exactly the same log arguments, with no combination of
logarithms or arctangents and no new choice of branch constants.

The independent derivative certificate works in K=QQ(x). For each squarefree q,
compute r=t*h_x/h modulo q in K[t]/q. The trace of multiplication by r, obtained
as the sum of coefficients of t^j in r*t^j mod q, equals the sum over all roots
of q of t*h_x/h. The sum of these traces is checked by exact cancellation to
f/g. This uses no sampled roots, floating tolerance or quadrature. Curated
repeated residues and mixed multiplicity groups supplement random coprime
squarefree examples, nonmonic inputs, rational coefficients and 80/200-bit data.

Protocol guards preserve the full old computation for replaced public Poly or
resultant constructors, constructor normalization, changed public sympify,
changed Expr.as_poly, Poly arithmetic/coefficient/factorization/Euclid methods,
custom attribute access, other coefficient domains, aliases and held evaluation.
Native Poly descriptor identities are captured before Python bootstrap adapters;
stock LC, coeffs, monoms, terms and as_expr adapters are recognized with their
underlying native descriptors. Constructor normalization has a canonical module
reference independent of rationaltools import order. Expr.as_poly is recognized
by its exact stock code object and globals, rejecting a user @wraps wrapper.

Signal checks surround checked PRS, factorization and exact division steps,
run between groups/coefficients and precede the completed result. A callback
error propagates without returning partial pairs. Budget or unsupported polynomial errors decline
the whole route and restart the complete existing implementation. No new Python
mathematical loop, runtime SymPy dependency, global domain change or precision
conversion was introduced. Tests use SymPy solely as an independent oracle.

This is not a full native port of general rationaltools: multivariate inputs,
parameters, algebraic/finite/inexact/expression/custom domains, unusual t,
nonpreconditioned inputs and larger dispatches retain their complete old path.
No DEBUG performance claim is made. A release-only runner records package SHA
and verifies full pairs/domains against SymPy outside timed loops; run it in
fresh alternating immutable baseline/after processes when the host is quiet.
