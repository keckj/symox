# Algebraic roots and quaternions

[README](../README.md)

These interfaces operate on native symbolic scalars. They do not call SymPy.
Root isolation runs in Rust; the `CRootOf`, `RootSum` and `Quaternion` classes
are ports of SymPy 1.14's reference classes on top of the native kernels.
`eval_approx` uses the declared mpmath dependency, as in SymPy.

## Polynomial roots

`RootOf`, `rootof`, `CRootOf` and `ComplexRootOf` accept a univariate
polynomial, an optional generator and a root index. Negative indices count from
the end. Roots are ordered with increasing real roots first, followed by complex
roots in SymPy's order. Multiplicities are retained. Floating coefficients are
converted to exact rationals first; after SymPy's root preprocessing, domains
other than `ZZ` (algebraic, symbolic) raise `NotImplementedError`, as in SymPy.

```python
from symox import symbols, RootOf, CRootOf
x = symbols('x')
RootOf(x**2 - 2, 0)                # -sqrt(2)
r = CRootOf(x**5 - x - 1, 0)
r.evalf(40)                        # 1.167303978261418684256045899854842180721
r.index                            # 0
r.free_symbols                     # set()
```

`RootOf` requests radical solutions when the complete radical list is
available. `CRootOf` retains nonrational roots as exact symbolic objects.
Rational roots are returned directly. A root is a native expression node whose
arguments are the defining polynomial expression and the index, printed and
serialized as in SymPy (`CRootOf(x**5 - x - 1, 0)`). Its generator is bound for
dependence queries but visible to `has()`.
The native `solve` pipeline uses these objects when a rational-coefficient
polynomial cannot be solved completely in radicals. For example,
`solve(x**5-x-1, x)` returns all five indexed roots. Repeated factors are
deduplicated, and real/sign/rational assumptions are checked using exact root
isolation facts. The `poly` and `index` properties recover its defining data.

Real roots are isolated by continued fractions over Mobius transforms and
complex roots by certified quadrant counting and rectangle bisection, as in
SymPy's `rootisolation`. `RealInterval` and `ComplexInterval` wrap the native
isolation states and refine them exactly; the SymPy root caches
(`clear_cache`, `_get_interval`, `_set_interval`, `_reset`) are available.
`evalf` refines the isolating interval or rectangle with exact rational
arithmetic (`eval_rational`) and rounds once at the requested precision;
`eval_approx` follows SymPy's mpmath `findroot` route. Isolation inherits the
native polynomial work/resource limits.
CRootOf generators of algebraic fields are described in
[coefficient domains](domains.md#algebraic-field-constructors).

`RootSum(p, Lambda(t, r(t)))` follows SymPy's constructor: linear factors are
evaluated directly, and rational functions of the roots are reduced exactly
through symmetric functions of the roots (`symmetrize`/`viete`). For example,
`RootSum(x**2-2, Lambda(x, 1/(y-x)))` gives `2*y/(y**2-2)`.
For other functions, or with `auto=False`, the sum remains symbolic and `.doit()`
expands it into individual roots. `evalf`, `N`, `latex`, `subs` and
differentiation are supported.

## Hamilton quaternions

`symox.algebras.Quaternion(a, b, c, d, real_field=True, norm=None)` is a port
of SymPy 1.14's `Quaternion` as a host expression class over native scalars. It
implements Hamilton multiplication: `i*j=k`, `j*i=-k` and
`i**2=j**2=k**2=-1`. Addition, subtraction, scalar operations, integer powers,
conjugation, norm (including a declared `norm`/`set_norm`), inverse and
normalization follow SymPy. Inverting a zero quaternion raises `ValueError`;
as in SymPy, normalizing it or taking the logarithm of a real quaternion with
zero vector part yields `nan` components.

`subs`, `expand`, `simplify`, `diff`, `integrate` and `evalf` act componentwise.
`exp`, `log`, `pow_cos_sin`, `scalar_part`, `vector_part`, `axis`, `angle`,
`mensor`, `index_vector`, the purity predicates (`is_pure`, `orthogonal`,
`parallel`, `vector_coplanar`, `arc_coplanar`) and the left/right product
matrices are available. Noninteger powers return `NotImplemented`, as in SymPy.

Rotations support `from_axis_angle`, `to_axis_angle`, `to_rotation_matrix`,
`from_rotation_matrix`, `rotate_point` and Euler angles through
`from_euler`/`to_euler` with SymPy's sequence conventions. A rotation center
passed to `to_rotation_matrix` produces a 4 by 4 homogeneous transformation.
The matrix-to-quaternion conversion includes rotations through pi. Its input is
expected to be a proper 3 by 3 rotation matrix; it does not certify orthogonality
for arbitrary symbolic matrices.
