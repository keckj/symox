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
`nu` and `delta`. `integrand(s)` gives the defining Mellin–Barnes integrand.
`get_period()` computes the rational-parameter monodromy period where the simple
pole-family criterion applies; unresolved or integer-separated families return
infinity. General parameter differentiation is held.

`hyperexpand` and `expand_func` use a bounded Slater expansion into native
hypergeometric functions. The expansion supports simple pole families with at
most 16 denominator parameters. Integer-separated or repeated pole families,
which require logarithmic residues, remain held. The `p > q` case uses inversion;
the `p == q` case retains a `Piecewise` distinction between `abs(z) < 1` and
`abs(z) > 1`, leaving the unit-circle boundary held. Double precision evaluation
uses the same certified regions and supported hypergeometric numerical kernels.
Arbitrary-precision general Meijer G evaluation is not implemented.

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

Numerical additions in these families are double precision. Unsupported high
precision requests do not gain arbitrary-precision special-function kernels;
use exact reductions where available and check whether the result stayed held.

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
