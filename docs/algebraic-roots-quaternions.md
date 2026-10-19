# Algebraic roots and quaternions

[README](../README.md)

These interfaces operate on native symbolic scalars. They do not call SymPy.
The numerical root refinement uses the declared mpmath dependency.

## Polynomial roots

`RootOf`, `rootof`, `CRootOf` and `ComplexRootOf` accept a univariate rational
polynomial, an optional generator and a root index. Negative indices count from
the end. Roots are ordered with increasing real roots first, followed by complex
roots ordered by real and imaginary part. Multiplicities are retained.

```python
from symox import symbols, RootOf, CRootOf
x = symbols('x')
RootOf(x**2 - 2, 0)                # -sqrt(2)
r = CRootOf(x**5 - x - 1, 0)
r.evalf(40)                        # 1.167303978261418684256045899854842180721
r.index                            # 0
r.free_symbols                     # set()
```

`RootOf` requests native radical solutions when the complete radical list is
available. `CRootOf` retains nonrational roots as exact symbolic objects.
Rational roots are returned directly. A root's coefficients and index live in
the native expression, so the root is independent of the polynomial generator.
The native `solve` pipeline uses these objects when a rational-coefficient
polynomial cannot be solved completely in radicals. For example,
`solve(x**5-x-1, x)` returns all five indexed roots. Repeated factors are
deduplicated, and real/sign/rational assumptions are checked using exact root
isolation facts. The numerical limitations below still apply.
The `poly` and `index` properties recover its defining data. The printed
coefficient-tuple representation differs from SymPy's polynomial representation.

Real numerical evaluation refines rational Sturm isolating intervals to the
requested decimal precision. Complex evaluation solves square-free rational
factors with mpmath and can fall back to refinement from native numerical seeds.
Complex approximations do not have certified complex isolating boxes; ordering
nearly coincident complex roots is numerical. This is not a full replacement
for SymPy's complex isolation, algebraic-field conversion or root cache APIs.
Polynomials with symbolic, floating or nonrational coefficients raise explicitly.
Isolation also inherits the native polynomial work/resource limits.

`RootSum(p, Lambda(t, r(t)))` evaluates rational functions of the roots exactly
using the trace of the rational function of a companion matrix. For example,
`RootSum(x**2-2, Lambda(x, 1/(y-x)))` gives `2*y/(y**2-2)`.
For other functions, or with `auto=False`, the sum remains symbolic and `.doit()`
expands it into individual roots. Singular denominators raise through exact
matrix inversion. Algebraic reduction of arbitrary expressions involving root
objects and complex isolation certification remain future work.

## Hamilton quaternions

`Quaternion(a,b,c,d)` stores four native scalars and implements Hamilton
multiplication: `i*j=k`, `j*i=-k` and `i**2=j**2=k**2=-1`. Addition, subtraction,
scalar operations, integral powers, conjugation, norm, inverse and normalization
are available. Zero inversion and normalization raise `ValueError`.

`subs`, `expand`, `simplify`, `diff`, `integrate` and `evalf` act componentwise.
The exponential uses the scalar/vector decomposition. The logarithm and `pow`
use the principal axis formula; a negative real quaternion with zero vector
part requires an axis choice and raises rather than inventing one.

Rotations support `from_axis_angle`, `to_axis_angle`, `to_rotation_matrix`,
`from_rotation_matrix` and `rotate_point`. A rotation center passed to
`to_rotation_matrix` produces a 4 by 4 homogeneous transformation. The
matrix-to-quaternion conversion includes rotations through pi. Its input is
expected to be a proper 3 by 3 rotation matrix; it does not certify orthogonality
for arbitrary symbolic matrices. General Euler-angle conventions and the full
SymPy Quaternion method inventory are not implemented.
