# Algebraic numbers

[README](../README.md) · [Coefficient domains](domains.md)

`AlgebraicNumber` retains a native exact field embedding, reduced rational
power-basis coordinates and an optional alias. Structural expression methods
inspect this representation, while mathematical calculations use its exact
scalar value. No SymPy import is needed for calculations.

## Structure

`atoms`, `find`, `count`, `has`, free-argument queries, class/sort keys,
comparison, matching, copying, pickling and `fromiter` preserve the wrapper.
`args` contains the field generator, the coordinate tuple and optional alias.
The mathematical `free_symbols` set is empty; structural queries can still find
the alias or numbers in the defining representation.

Coefficient, content, factor and numerator/denominator views treat the number
as a single commutative constant. They do not replace it by a displayed alias
polynomial. `as_terms` follows the public numerical-coefficient convention;
`count_ops` counts operations in the defining arguments. `as_poly` retains the
existing power-basis polynomial convention.

`subs`, `xreplace` and `replace` rebuild the field generator, coordinates and
alias when their defining arguments change. Replacing the entire number can
return another expression. Substitution mappings process larger patterns first;
sequences preserve their specified order. Simultaneous algebraic substitutions
mask their right-hand sides. A resulting generator must still belong to the
supported native field models.

## Mathematical operations

Differentiation, integration, limits, inversion, rounding and simplification
use native scalar kernels. Unchanged values retain their field and alias through
simplification. For example, simplifying an aliased element with coordinates
`[2, 3]` preserves its value and coordinates, including cases where SymPy 1.14's
root simplifier drops the coordinates.

Constant series, leading terms, Taylor terms, formal/Fourier series and local
regularity queries preserve the algebraic constant. Complex conjugation
constructs the conjugate exact embedding with the same rational coordinates.
Real/imaginary parts use native scalar projections; real elements retain the
original wrapper as their real part.

`to_algebraic_integer` follows the field-generator convention: a generator with
nonmonic primitive integer minimal polynomial is scaled by its leading
coefficient. Power-basis coordinates remain unchanged. This changes the field
generator and can change the represented element; it is not an assertion that
arbitrary rational power-basis coordinates are integral.

Explicit `to_sympy()` / `_sympy_()` conversion preserves `AlgebraicNumber`, its
coordinates and alias. Scalar expressions produced by arithmetic or calculus
use the existing native scalar embedding; they do not retain an opaque
AlgebraicNumber node inside the scalar DAG. Thus name presence does not imply
full SymPy tree parity.

Structural matching accepts `old=True` as well as the default mode. Alias
wildcards bind native symbols, and a wildcard matching the whole number returns
the original wrapper. Exclusions inspect the defining root, coordinates and
alias. General scalar patterns containing an algebraic-number literal remain
limited by the absence of opaque algebraic-number nodes in the scalar DAG.

## Field models

The root may be any generator accepted by the general algebraic field
constructors (radicals of any index, nested radicals, `CRootOf`); see
[algebraic field constructors](domains.md#algebraic-field-constructors).
Minimal polynomials, primitive elements and embeddings are exact.
`to_root`, `is_aliased` and `to_primitive_element` follow SymPy, without
rewriting the root into radicals. A purely symbolic generator such as
`AlgebraicNumber((T, theta))` with a `Symbol` `theta` raises `ValueError`.
