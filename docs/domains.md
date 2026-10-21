# Coefficient domains

[README](../README.md)

The native typed domains are `ZZ`, `QQ`, `GF(p)`, `ZZ_I`, `QQ_I`, real and complex
numerical fields, `EX`, `EXRAW`, supported algebraic extensions of `QQ`,
polynomial rings over `ZZ`/`QQ`, and their rational-function fields.
Elements retain their parent domain through arithmetic, conversion and pickling.
All operations below run without importing SymPy.

`construct_domain` accepts a scalar, iterable or dictionary of coefficients
and returns their inferred native domain and converted elements, with the
reference options (`field`, `extension`, `composite`, and the other `Poly`
options, which are validated and otherwise ignored); dictionaries retain
their keys and other iterables return lists. Unlike the raw expression
fallback of matrix inference, the general constructor uses `EX`. Only
`extension=True` enables algebraic inference, as in the reference. Composite
generators may be Symbols or other non-algebraic expressions (`ZZ[pi]`,
`ZZ[y,sin(x)]`, `ZZ(x,pi)`), and their elements print like the reference
`PolyElement`. Floats with a generator give `EX` where the reference builds
`RR[x]` or `RR[pi]`: composite domains over `RR` and `CC` are not
implemented.

Constructor facades include `IntegerRing`, `RationalField`, `FiniteField`,
`RealField`, `ComplexField`, Gaussian and expression domains, `PolynomialRing`,
`FractionField` and `AlgebraicField`. They return native `Domain` descriptors;
their elements use the common native `DomainElement` protocol. Backend-specific
storage engines from SymPy are not reproduced.

### Legacy backend constructor names

`PythonIntegerRing`, `GMPYIntegerRing`, `PythonRationalField`,
`GMPYRationalField`, `PythonFiniteField` and `GMPYFiniteField` construct the
corresponding native domains. `ZZ_python`, `ZZ_gmpy`, `QQ_python`, `QQ_gmpy`,
`FF_python` and `FF_gmpy` are their compatible constructor aliases. The usual
`pythonintegerring`, `pythonrationalfield`, `pythonfinitefield` and corresponding
`gmpy...` module paths are available. These spellings do not select Python or
GMPY storage; coefficients use symox's exact native kernels, with no GMPY
dependency. Prime-modulus and symmetric-display restrictions of native `GF`
still apply; `symmetric=False` raises `NotImplementedError`.

`PythonRational(numerator, denominator=None)` is an immutable exact rational
adapter backed by Python integer fractions. It normalizes signs and common
factors, supports integer, string, decimal, float and fraction construction,
and preserves its type through rational arithmetic and integer powers.
Comparisons and hashing follow exact fraction semantics. Its native conversion
hook allows `sympify` and native domain constructors to consume it directly.
Floor division and remainder are deliberately unsupported, following the
rational coefficient API. This adapter is available from `pythonrational`;
explicit `to_sympy()` is optional and is the only operation that imports SymPy.

## Exact arithmetic

`abs`, sign predicates, `canonical_unit`, `pos`, `almosteq`, `cofactors`,
`half_gcdex`, `revert`, `sum`, recursive `map`, `new`, `normal` and `tp`
extend the native arithmetic protocols. Exact domains compare exactly even
when `almosteq` receives a tolerance. Polynomial signs follow the leading
coefficient in the domain's monomial order; finite-field signs follow SymPy's
formal convention. `real` and `imag` follow the coefficient-domain convention,
including for algebraic fields; they are not complex conjugation projections.

`exsqrt` returns an exact domain element or `None`: rational roots use integer
square roots of numerator and denominator, finite fields use modular square
roots, and polynomial/function domains factor numerator and denominator.
Algebraic square roots use the supported native embedding conversion and verify
the result by squaring. `is_square` tests this result. `ZZ.sqrt` returns the
integer floor square root; `sqrt` outside `ZZ` has no inexact convention and
raises `NotImplementedError`, as on SymPy's rational domain. `ZZ.log` computes
an integer floor logarithm exactly, including very large integers.

`convert_from` validates the native source domain and raises `CoercionFailed`
when conversion fails. Named conversion hooks for implemented integer,
rational, finite, algebraic and polynomial domains return `None` when there is
no embedding. Finite fields have no canonical embedding into `QQ`. A constant
polynomial can convert to a ground domain; a rational function is not
specialized to a value by a named fraction-field hook.

`drop`, `inject`, `unify_composite` and `unify_with_symbols` preserve native
generator identities and monomial orders. `old_poly_ring` and `old_frac_field`
construct the same native rings as their current counterparts; they do not
select a second legacy representation.

## Modules, ideals and quotients

```python
from symox import QQ, ZZ, symbols

x, y = symbols('x y')
module = QQ.free_module(3)
subspace = module.submodule([1, 2, 0], [0, 1, 3])
assert [2, 5, 3] in subspace

ring = QQ.old_poly_ring(x, y)
ideal = ring.ideal(x**2-y, y**2-1)
assert ideal.contains(x**4-1)
quotient = ring / ideal
assert quotient(x)**4 == quotient.one

integers_mod_15 = ZZ / [15]
assert integers_mod_15(2)**-1 == integers_mod_15(8)
```

Free module coordinates are native elements. Vector arithmetic and the standard
basis work over supported domains. Submodule membership and reduction use typed
Gaussian elimination over fields, including finite and supported algebraic
fields. Submodule elimination over rings is not implemented.

Ideals support membership, reduction, equality, containment, sums and products.
Integer ideals reduce modulo their positive gcd. Field ideals are zero or the
whole field. Polynomial ideals over a field use the native reduced Gröbner
basis. Ideals over `ZZ[x]` are explicitly unsupported; they are not computed by
silently extending coefficients to `QQ`.

Quotient elements support canonical reduction, addition, subtraction,
multiplication, integer powers and compatible-ring comparisons. Inversion
supports integer residue rings, zero ideals in fields, and principal
univariate polynomial quotients. General multivariate quotient inversion is
not implemented. Quotient types are Python protocols over native coefficients;
they are not yet accepted as coefficient backends by `Poly` or `DomainMatrix`.

## Algebraic field constructors

`QQ.alg_field_from_poly(poly, alias=None, root_index=-1)` selects a root using
the native `all_roots` implementation, then constructs the field generated by
that root. Arbitrary irreducible polynomials are accepted, as in the reference:
`QQ.alg_field_from_poly(Poly(x**5 - x - 1))` is `QQ<CRootOf(x**5 - x - 1, 4)>`.
`QQ.cyclotomic_field` applies this to the cyclotomic polynomial of any order
and supports alias subscripts; orders 1 and 2 yield `QQ`.

`QQ.algebraic_field` and `extension=` accept any finite list of algebraic
generators (radicals of any index, nested radicals, `CRootOf`). The primitive
element, minimal polynomials and generator embeddings are exact; when several
candidate factors remain, the factor is chosen with a certified enclosure.
Membership and conversion into a field factor the minimal polynomial over the
field and identify the matching factor in the same certified way.

## Numerical fields

`RR` and `CC` use 53-bit native arbitrary precision arithmetic by default.
`RealField(200)` and `ComplexField(dps=60)` select a precision. Supported limits
are 2-40,000 bits or 1-10,000 decimal digits. Custom tolerance constructor
options are explicitly unsupported.

Decimal strings, including decimal literals in compound numerical expressions,
are parsed directly at the field precision. Python floats retain their actual
binary value. Numerical conversion rejects free symbols, nonfinite inputs and
unsupported numerical functions. `CC(real, imag)` constructs a complex value;
real/imaginary element properties return real-field elements at the same
precision. Numerical equality and hashing compare stored values, including
compatible Python integers, floats, fractions and complex values.

`almosteq` uses an absolute/relative tolerance derived from field precision;
callers can supply a nonnegative real tolerance. `exsqrt` computes numerical
principal roots; negative real arguments have no real root. As with SymPy,
the separate inexact `sqrt` domain convention remains defined only for `ZZ`.

`to_rational(value, limit=False)` returns the exact binary value as an integer
pair. The default limits the denominator using the real-field convention and
preserves small nonzero values when a bounded denominator would erase them.
Conversions to `QQ` and `QQ_I` use this convention. Associated exact fields are
`QQ` for real fields and `QQ_I` for complex fields.

## Gaussian and expression domains

`ZZ_I` and `QQ_I` store exact real/imaginary rational coordinates in the native
quadratic-field kernel. `domain(real, imag)` constructs a value. Integer-ring
division computes a nearest Gaussian quotient and a remainder of smaller norm;
GCD and extended GCD provide exact Bézout identities. Canonical units follow
the Gaussian quadrant convention. LCM uses the normalized GCD and preserves
the product's remaining unit, matching SymPy's convention. Rational-field
numerators and denominators belong to `ZZ_I`.

Exact Gaussian square roots test rational norms and coordinates. Ring division
stays in its parent and requires an exact quotient; use `get_field()` for
fractions. Gaussian integer ideals and modular inversion use Euclidean GCD.
As with other exact typed coefficients, an absolute value outside the selected
field cannot be converted back into that field.

`EX` normalizes scalar expressions with the native simplifier; `EXRAW` retains
their native presentation. Both provide closed expression arithmetic and an
`ex` scalar view. Sign predicates follow the expression-domain coefficient
conventions. `ZZ_I`, `QQ_I` and `EX` are polynomial coefficient backends
with the reference ring algorithms (division, GCD, square-free and complete
factorization, roots). Custom numerical domain objects are not; unsupported
selections are rejected explicitly. Ring factories accept `ZZ_I`, `QQ_I`, algebraic
fields and polynomial rings as ground domains (`ZZ_I[x]`, `QQ<sqrt(2)>[x]`,
`ZZ[x][y]`), with reference printing, conversion and `unify`; `Poly` works
over nested rational rings such as `ZZ[x][y]` and over the Gaussian
composite grounds `ZZ_I[y]`, `QQ_I[y]`, `ZZ_I(y)` and `QQ_I(y)`, which
`Poly` inference selects as the reference does (`Poly(I*x*y, x)` is over
`ZZ_I[y]`, `Poly(-2*I*pi*n/N, n)` over `ZZ_I(N,pi)`). These polynomials
print, divide, factor and compute `gcd`, `sqf` and `cancel` like the
reference; `ZZ_I[x].is_unit(I)` is `False`. `Poly` over `QQ<sqrt(2)>[x]`
still raises `NotImplementedError`. Fraction fields over Gaussian or
algebraic grounds have no domain elements and raise `NotImplementedError`
when elements are built, so `construct_domain(..., field=True)` falls back
to `EX` for them. See
[the parity gaps](algorithmic-parity.md) for the remaining differences.

## Monogenic finite extensions

`symox.polys.agca.extensions.MonogenicFiniteExtension` constructs a quotient
of a native univariate polynomial ring over `ZZ` or `QQ` by a monic polynomial.
Exact leading-coefficient division normalizes the modulus without promoting
an integer ground ring. Extension rank is limited to 4096.

Elements support exact addition, multiplication, integer powers, reduction,
unit inversion, basis coordinates and pickling. Reducible moduli retain their
zero divisors and are not marked as fields. Integer-ground inversion requires
integral inverse coefficients. These are abstract quotients: a nonconstant
element has no chosen scalar algebraic embedding. Named conversion hooks and
`convert_from` convert ground constants and reject a missing embedding.
