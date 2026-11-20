# Gamma, hypergeometric, Bessel and Mathieu functions

These families use native Rust expression heads, scalar operations and
numerical kernels. SymPy is used by differential tests and explicit conversions;
it is not a runtime computation backend. Compatible imports are available under
`symox.functions.special.gamma_functions`, `.hyper`, `.bessel` and
`.mathieu_functions`, as well as at the package root.

The coverage table counts exported names, including abstract bases and helpers.
It does not certify complete signatures, algorithms, precision or branch parity.

## Gamma completion

`multigamma(x, p)` evaluates the defining finite gamma product for a positive
integer dimension, including the power of pi. Explicit dimensions are limited
to 1024. Symbolic dimensions are represented by a native `Product`; dimensions
known to be nonpositive or nonintegral raise `ValueError`. The resulting gamma
expression supports native differentiation and evaluation.

`intlike` follows SymPy's helper convention for integer-valued inputs. It tests
integer conversion followed by the residual, so Python integer-valued floats
are accepted. Numeric strings are not numeric inputs to the residual operation.

## Hypergeometric functions

The existing `hyper` function is joined by `appellf1` and `meijerg`.

`appellf1(a, b1, b2, c, x, y)` canonicalizes the interchangeable parameter/argument
pairs and provides native partial derivatives in `x` and `y`. Both zero arguments
reduce to one. Numerical evaluation uses a diagonal convolution of the double
series: nonterminating series require `abs(x), abs(y) < 0.98`; terminating series
can use larger arguments. Evaluation is limited to 2048 diagonals. Analytic
continuation outside these domains and parameter derivatives remain unsupported.

`meijerg(an, aother, bm, bother, z)` also accepts the three-argument form
`meijerg((an, aother), (bm, bother), z)`. Parameters are sorted and checked for
finite values and overlapping upper/lower pole families. The native expression
preserves both parameter partitions and supports argument differentiation,
substitution, pickling and explicit conversion to SymPy.

Meijer G properties include `an`, `aother`, `ap`, `bm`, `bother`, `bq`, `argument`,
`nu` and `delta`. `integrand(s)` gives the defining Mellin-Barnes integrand.
`get_period()` computes the rational-parameter monodromy period where the simple
pole-family criterion applies; unresolved or integer-separated families return
infinity. General parameter differentiation is held.

`hyperexpand(f, allow_hyper=False, rewrite='default', place=None)` and
`expand_func` are a Rust port of SymPy's `sympy.simplify.hyperexpand`
(`crates/symox-core/src/special/hyperexpand/`): order reduction, the table of
known origins, Lerch phi and polylog bases, terminating and shifted sums, and
plans of differential shift operators from an origin to the target. Meijer G
functions use a table formula when one matches and Slater's theorem otherwise,
in `z` and in `1/z`, with SymPy's convergence conditions, `Piecewise`
combination and `place` keyword. For example
`hyperexpand(meijerg([1, 1], [], [1], [0], z))` is `log(z + 1)` and
`hyperexpand(hyper([1, 1, 1], [2, 2], z))` is `polylog(2, z)/z`. As in SymPy,
Slater's theorem is skipped when a pole family has more repeated or
integer-separated lower parameters than matching upper ones; without a table
formula the function stays held (`meijerg([], [], [0, 0], [], z)`). The
compatible `symox.simplify.hyperexpand` module exposes the operator classes,
`devise_plan` and the formula collections over the native code. Remaining
form differences and timeouts are tracked in
[the parity document](algorithmic-parity.md). Numerical Meijer G
evaluation at any precision follows SymPy's `meijerg._eval_evalf`; see
[Multiprecision evaluation](#multiprecision-evaluation).

`TupleArg` supplies tuple-wise leading terms, limits and substitution.
`TupleParametersBase` provides the shared tuple-parameter function protocol.
The `HyperRep_*` helpers preserve native symbolic heads and implement the
small-argument and winding-dependent analytic continuation formulas used by
hypergeometric expansion. `rewrite('nonrepsmall')` selects the local formula;
`rewrite('nonrep')` retains the large/small-argument distinction. Continuation
helpers require a concrete integer winding number. They do not add a general
polar-function or analytic-continuation engine to symox.

## Bessel completion

`hn1(n, z)` and `hn2(n, z)` are native spherical Hankel functions. They support
spherical derivative recurrences, conjugation, integer-order elementary
expansion, and rewrites to spherical and ordinary Bessel/Hankel functions.
The ordinary Hankel relation uses order `n + 1/2`, as required by the defining
identity. SymPy 1.14's direct `hankel1`/`hankel2` rewrite uses `n` instead; the
native implementation uses the defining mathematical relation.

`BesselBase`, `SphericalBesselBase`, `SphericalHankelBase` and `AiryBase` expose
family checks and shared argument/derivative protocols on native expressions.
`assume_integer_order` gates a decorated rewrite on a provably integral order.

`jn_zeros(n, k, method='sympy', dps=15)` returns the first `k` positive spherical
Bessel zeros, using native Bessel evaluation, bracketing and bisection. Orders
must lie between zero and 128, with at most 10,000 roots and one through 15
requested decimal digits. The compatibility method names `sympy` and `scipy`
select the same native algorithm; neither library is imported. Higher precision
and other method names raise `NotImplementedError`.

`marcumq(m, a, b)` supports the zero-argument and coincident-argument reductions,
argument derivatives, and rewrites as an integral, infinite Bessel sum or
applicable finite Bessel sum. Numerical evaluation uses a Poisson mixture of
regularized incomplete gamma values for real, nonnegative parameters, with
`a*a/2 <= 4096`. Complex continuation and general order derivatives remain held.

`marcumq` numerical evaluation is double precision; SymPy 1.14 has no numerical
kernel for it either. The other families evaluate at any precision as described
in [Multiprecision evaluation](#multiprecision-evaluation).

## Multiprecision evaluation

Special functions that SymPy evaluates through mpmath evaluate through the same
mpmath function in the Python package: Bessel, Hankel, spherical Bessel and
Airy functions with their derivatives, gamma, loggamma, polygamma, incomplete
gamma, error functions, exponential, logarithmic, sine and cosine integrals,
Fresnel integrals, zeta, Dirichlet eta, polylog, Lerch phi, beta, elliptic
integrals, Lambert W, `hyper`, `appellf1`, `meijerg`, orthogonal polynomials,
spherical harmonics, and the Bernoulli, harmonic, Euler and Andre extensions to
non-integer arguments. Real erf, erfc and gamma, and the functions with native
multiprecision kernels, stay native.

The native evaluator passes exact binary arguments and certifies the result:
it evaluates at the requested precision plus 10 bits, then at a larger guard
(10, 30, 70, ... bits), until two values agree to the requested precision,
with at most `2*prec + 512` bits. Arguments are rounded at the working
precision plus 12 bits, so ill-conditioned arguments such as
`elliptic_k(1 - 10**-30)` or `zeta(1 + 10**-20)` keep every requested digit
where SymPy's single evaluation loses some. If the values never agree the
expression stays unevaluated; nothing is completed in double precision. When
an inexact argument next to a pole, branch point or domain end (such as
`gamma(-3 + 10**-20)`) was rounded onto it, the evaluation is retried with
more bits before the error is reported.

`elliptic_pi` uses `symox._ellippi`, which follows mpmath's `ellippi` step by
step but moves the initial Carlson `R_J` integral away from the integrand's
singular points, where mpmath's quadrature was slow or inaccurate.

Errors raised by mpmath propagate as in SymPy: `erfinv(2.7)`, `zeta(1.0)`,
`loggamma(-2.0)` and `lerchphi(2, 1, 0).evalf()` raise `ValueError`.
Assumption queries and internal numerical probes do not raise. A non-finite
mpmath result leaves the expression unevaluated.

A special function with Float arguments evaluates on construction under
SymPy's `_should_evalf` rule: every argument is a Float or a pure complex
number `a + b*I` with Float parts, and the evaluation uses the largest Float
precision. `erf(2.0*I)` is a product, not a pure complex number, so it stays
held as in SymPy.

`meijerg(..., z*exp_polar(theta*I))` follows SymPy's branch handling: a
principal `theta` in `(-pi, pi]` evaluates on the principal sheet, and
another branch is carried by mpmath's `r` argument with
`n = ceiling(abs(theta/pi)) + 1`. `exp_polar` itself evaluates only on the
principal sheet, and `exp_polar(I*pi)` evaluates to `-1.0`.

Without the Python package (the bare Rust library), special functions keep
the native kernels, which are double precision for the families above.
`hyper([2, 3, 4], [1], 1/2)` and other cases where mpmath converges slowly
cost a few times SymPy's single evaluation because of the certification pass.

## Mathieu functions

`mathieuc`, `mathieus`, `mathieucprime` and `mathieusprime` are native functions
with the shared `MathieuBase` protocol. Their derivatives with respect to `z`
satisfy the defining equation exactly:

```python
from symox import symbols, mathieuc, cos, expand

a, q, z = symbols('a q z')
f = mathieuc(a, q, z)
assert expand(f.diff(z, 2) + (a - 2*q*cos(2*z))*f) == 0
```

They support parity under `z -> -z`, conjugation and the elementary `q = 0`
reductions. Parameter derivatives and general numerical evaluation at nonzero
`q` stay held, following the symbolic scope of SymPy's Mathieu interface.
Held Mathieu construction with `evaluate=False` is currently unsupported.
Characteristic values, periodic normalization and a Mathieu spectral solver
are not implemented by these exports.

## Additional root exports

`betainc` and `betainc_regularized` support exact positive-integer polynomial
reductions, endpoint identities, endpoint/parameter partial derivatives, symbolic
substitution and arbitrary-precision real/complex mpmath evaluation. Composed
expressions use the same differentiation and evaluation adapters. `riemann_xi`
and `stieltjes` follow SymPy's `eval`: `riemann_xi` is `1/2` at zero and one
and evaluates whenever `zeta(s)` does (`riemann_xi(2)` is `pi/6`), and
`stieltjes(0)` is `EulerGamma` while negative or non-integer numeric indices
give `zoo`. Both evaluate numerically with mpmath.
These constructors are Python adapters over native scalar nodes
(`python/symox/_public_special.py`); they never delegate calculations to SymPy.
The general adaptive error estimates of SymPy's evalf machinery are not provided.

`erf2inv(x,y)` applies SymPy's special values and otherwise returns
`erfinv(y+erf(x))`, where SymPy keeps `erf2inv` held. `Ynm_c` conjugates
spherical harmonics; `Znm` forms real spherical harmonics for orders of known
sign. `chebyshevt_root` and `chebyshevu_root` implement the exact
trigonometric root formulas with validated integral degree/index.
`hermite_prob` follows SymPy's `eval`: explicit polynomials for nonnegative
integral degrees, parity in the argument and the values at zero and infinity
for symbolic degrees; negative numeric degrees stay held.
`jacobi_normalized` divides the native Jacobi polynomial by its gamma-function
normalization. Singular parameter cases inherit the underlying gamma/Jacobi
limits; the adapter does not resolve every removable parameter singularity.

The [splines and variational calculus guide](splines-variations-integer-equations.md)
describes B-splines, interpolation and singularity functions.
