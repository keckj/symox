# Scope and preservation

This patch extends the delivered conditional-rational package in a new private copy. Only the default methods with conds='piecewise' enrich results that were retained because rational-pole analysis returned Unknown. Already established regular or divergent results keep their old paths. Explicit method switches and conds='none'/'separate' are unchanged. No new Python algorithm, numeric threshold or quadrature is used.

These dispatch limits concern the new Unknown routes; the already Regular cases of release 146, including the univariate imaginary-parameter case, remain unchanged.

The selected shapes are exactly 1/x with limits (x,a,1), and 1/(x**2+a) with limits (x,0,oo), with a commutative parameter independent of the bound variable. The first result is -log(a); the second is pi/(2*sqrt(a)). Both require

    Contains(a, Complexes) & ~Contains(a, Interval(-oo,0))

where the right endpoint is closed. The default branch is the complete original Integral. Complexes membership remains symbolic for parameters, including real/complex/nonzero symbols, so later nonfinite substitutions remain guarded. The predecessor's conservative affine nonzero-parameter fallback remains unchanged. A new slit-domain case also declines when cut exclusion is already True for a free parameter (for example imaginary=True); otherwise ordinary condition simplification would discard that guard before a later contradictory substitution. Nonzero constants remain admissible, and established Regular paths are untouched.

# Straight logarithmic segment

Let a be finite complex and outside the closed negative real ray. Parametrize the straight segment a->1 by z(t)=(1-t)*a+t, 0<=t<=1. If Im(a) is nonzero, Im(z(t))=(1-t)*Im(a) has its fixed nonzero sign for t<1 and the endpoint is the positive real number 1. Thus the segment never reaches zero or the principal log cut. If Im(a)=0, the domain condition implies a>0, and the whole segment is positive real.

The principal logarithm is consequently continuous and differentiable along the segment, including its endpoints, and d(log(z(t)))/dt=z'(t)/z(t). Its integral is log(1)-log(a)=-log(a). The compact segment stays a positive distance from zero, so the integral is finite. The domain proof applies to complex a, including negative real parts, without replacing log by a branch-insensitive primitive.

If a=0 the starting endpoint is a pole. If a is negative real, the straight real segment crosses the pole at zero. Neither case is assigned an analytic continuation or a principal value; the old full Integral evaluates or remains held exactly as before.

# Quadratic principal-root branch and absolute convergence

For finite a off the closed negative real ray, write a=r*exp(i*theta), r>0, -pi<theta<pi. The principal root s=sqrt(a)=sqrt(r)*exp(i*theta/2) has Re(s)>0. Thus the roots of x**2+a, namely +/-i*s, have nonzero imaginary parts +/-Re(s), and there is no pole on the real integration ray.

The exact partial fraction identity is

    1/(x**2+s**2) = [1/(x-i*s)-1/(x+i*s)]/(2*i*s).

For x>=0, x-i*s stays in the lower open half-plane and x+i*s in the upper open half-plane. Their principal logarithms therefore have no jumps along these horizontal paths. Integration from 0 to R gives

    [log(R-i*s)-log(R+i*s)-log(-i*s)+log(i*s)]/(2*i*s).

As R tends to positive infinity, the first two logs have equal log(R) contributions and their difference tends to zero. If Arg(s)=phi, then -pi/2<phi<pi/2; hence Arg(i*s)=phi+pi/2 is in (0,pi), while Arg(-i*s)=phi-pi/2 is in (-pi,0). Both numbers have the same modulus, so log(i*s)-log(-i*s)=i*pi exactly. The integral is therefore pi/(2*s), with the principal root, and no unknown multiple of 2*pi.

Absolute convergence follows independently: for x**2>=2*abs(a),

    abs(x**2+a) >= x**2-abs(a) >= x**2/2,

so the absolute integrand is at most 2/x**2, an integrable tail. Near zero, when x**2<=abs(a)/2, the denominator has modulus at least abs(a)/2. On the intervening compact interval it is continuous and nonzero. This proves convergence without numerical integration.

For a<0 real a positive real pole exists; a=0 gives the nonintegrable endpoint 1/x**2. Those values remain in the complete old Integral branch, preserving independent one-sided divergence rather than a Cauchy principal value.

# Validation and retained paths

Tests certify the unique linear pole, quadratic roots, exact partial fraction identity and primitive derivatives symbolically. They compare the exact requested output formulas and complete Integral defaults, regular real/complex substitutions, both sides within 10**-100 of the cut, cut/zero/nonfinite values, and repeated doit. These examples accompany the proof above rather than substitute sampling for it.

The fixed-point Piecewise.doit guard is inherited without modification from the previous delivered package. No functions.rs delta is part of this patch. Historical assertions for the newly supported 1/x lower-bound family now require the exact Piecewise/Integral default instead of a bare Integral; unported cases remain strict.

Separate differential probes preserve the full result trees and exception types/messages of the predecessor binary for all 20 explicit-option combinations and all 112 nonfinite substitutions (generic/real/complex/nonzero/positive/imaginary/nonreal parameters, both families, before and after doit). Existing positive-parameter results remain unchanged. No DEBUG performance comparison is claimed.

# Remaining substitution-system limitation

These guards cover the original parameter assumptions and the tested direct numeric substitutions before/after doit. They are not a general port of assumption-contradicting symbolic substitution chains. For example, replacing the generic parameter by a fresh b declared imaginary can make ordinary Contains evaluation discard the cut guard; subsequently replacing that b by -1 contradicts its newly introduced assumption and can give -i*pi for the lower-bound family. The concrete witness is archived in contradictory-symbolic-chain.log. Preserving such guards through arbitrary symbolic replacements requires separate substitution/held-condition work; it is not claimed resolved by this bounded mathematical extension.

The analogous initial imaginary-parameter unit-interval example is already Regular in release 146, with identical expression and negative-half specialization before and after the delivered unit package. The two unit-imaginary JSON records document that inherited behavior; no modification of the frozen previous package is included here.
