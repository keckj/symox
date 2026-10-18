# Polynomial ideals and linear systems

[README](../README.md)

`ZZ`, `QQ` and `GF(prime)` are native domain objects, also available from
`symox.polys.domains`. `QQ.algebraic_field(...)`, `QQ.frac_field(...)` and
`QQ.poly_ring(...)` construct supported exact extensions. Domain elements retain
their domain under arithmetic; conversions check membership instead of merely
changing a label.

```python
import symox as sm

a, b, x = sm.symbols("a b x")
K = sm.QQ.frac_field(a, b, order="grevlex")
u = K.convert((a + 1)/(b - 1))
assert u * u**-1 == K.one
p = sm.Poly(x**2 + x/a + b, x, domain=K)
assert p.domain == K and p.domain.is_Field
```

For supported coefficient domains, `Poly.domain` and `GroebnerBasis.domain`
return domain objects. Use `str(domain)`
for their display names; a domain does not compare equal to a string.
Typed domains preserve generator assumptions and order, unused parameters,
algebraic embeddings and aliases through polynomial operations and pickling.
`set_domain`, `to_field` and `to_ring` validate the resulting coefficients.
`monic(auto=False)` keeps the coefficient ring and requires exact division by
the leading coefficient; the default promotes rings to fields, including for
the zero polynomial. `clear_denoms(convert=True)` places the numerator in the
associated ring, while the default retains its coefficient domain. Denominator
units respect the ground domain and monomial order: `QQ(a)` clears into `QQ[a]`
without clearing rational constants, whereas `ZZ(a)` clears them too.
`Poly.diff` accepts repeated `(generator, order)` specifications and signed
generator indices; invalid indices are checked before native array access.
`eval` accepts a value, a generator/value pair, a list or tuple of leading
generator values, or a generator mapping. Partial evaluation returns a `Poly`
with the remaining generators and the selected domain, including unused
parameters, aliases and monomial orders. Finite-field results reduce in their
coefficient field. `auto=False` retains the coefficient domain and rejects
points outside it; this also applies within lists and dictionaries. The default
promotes to an implemented common coefficient domain. General scalar results
can use expression or Gaussian domains; partial results still require a
supported `Poly` coefficient backend. Legacy coefficient representations without
a typed domain keep their native evaluation fallback and reject `auto=False`.
Calling a polynomial with several positional values uses the same partial
evaluation protocol and checks arity.

`as_expr(*values)` replaces all generator positions. `as_expr(mapping)` replaces
specified original positions simultaneously, including swaps. SymPy 1.14's
mapping implementation can accidentally cancel a swap by searching its already
modified generator list; symox uses the original positions. Coefficient-domain
parameters are not keys of this mapping.
`Poly.integrate` promotes coefficient rings to their fraction fields by default;
`auto=False` retains the ring and raises `ExactQuotientFailed` when division
is not exact, including multivariate cases where SymPy 1.14 silently truncates
integer coefficients. Finite-field integration uses modular inverses and raises
`NotInvertible` when the required denominator vanishes in the field.
Finite fields currently require a prime characteristic below `2**32`. Domain
elements share a native `DomainElement` representation; `K.dtype` constructs and
recognizes elements of `K`. Converting real-field elements into `QQ` uses
bounded-denominator rational reconstruction; requesting the exact stored binary
ratio is available through the field's `to_rational(..., limit=False)`.
Rational-function and
polynomial coefficient domains use `ZZ` or `QQ` ground coefficients; additional
towers require separate support. Finite-field elements compare within their
own domain, avoiding SymPy's inconsistent equality/hash behavior with integers.
Polynomial coefficient rings support division with remainder and exact division, with
coefficient content preserved in GCDs, LCMs and cofactors. `Poly.div` and
`Poly.exquo` honor `auto=False` when retaining the coefficient ring.
The global `div`, `rem`, `quo` and `exquo` functions honor explicit generators,
`auto`, `polys` and supported coefficient-domain options. They preserve typed
domains when returning polynomials and reject unknown or incompatible options.
`Poly.pdiv`, `prem`, `pquo` and `pexquo`, and their global counterparts, perform
pseudo-division in the selected coefficient ring. The quotient and remainder
satisfy `LC(g)**(degree(f)-degree(g)+1)*f == q*g+r` when `degree(f) >= degree(g)`.
The implementation checks exponent, coefficient, term and work bounds before
allocating large intermediates. The global pseudo-division functions reject the
`auto` option with `FlagError`, following SymPy's option contract.
Fraction fields over finite or algebraic ground domains require further support.
Existing polynomial calculations with general expression coefficients remain
available; requesting a typed domain for unsupported expression or inexact
coefficient rings raises `NotImplementedError`.

`groebner` computes reduced Gröbner bases with native Buchberger and F5B
algorithms over `QQ`, `ZZ`, prime finite fields and rational-function coefficient fields
such as `QQ(a,b)`. Polynomial coefficient rings `QQ[a,b]` and `ZZ[a,b]`
are also supported, with deterministic parameter inference and coefficient
normalization. It supports deterministic `lex`,
`grlex` and `grevlex` orders:

```python
import symox as sm

x, y = sm.symbols("x y")
g = sm.groebner([x*y - 1, y**2 - 1], x, y, method="f5b")
g.reduce(x**2 - 1)               # quotients and remainder
g.contains(x**2 - 1)             # True
```

Bases support iteration, indexing, `exprs`, `polys`, `gens`, `domain`,
`order` and `is_zero_dimensional`. `g.fglm(order="lex")` converts a
zero-dimensional basis using exact quotient-algebra multiplication matrices
and linear dependence, over the supported rational, prime finite,
rational-function and supported algebraic fields.
Changing the monomial order requires a zero-dimensional ideal; requesting the
existing order returns the basis unchanged. Polynomial coefficient rings use
their fraction fields internally and restore integral coefficient normalization.
`method="buchberger"` remains the default; `method="f5b"` selects signature-based
reduction with rewrite and syzygy criteria. Both return reduced bases over all
supported coefficient fields; their relative speed depends on the system.

```python
a, b, x, y = sm.symbols("a b x y")
g = sm.groebner([a*x + y, x*y - b], x, y, domain="QQ(a,b)")
g.reduce(x*y - b)                 # zero remainder
p = sm.Poly(x**2 + x/a + b, x, domain="QQ(a,b)")
p.diff()                         # native exact coefficient arithmetic
```

`Poly.subresultants`, `resultant(..., includePRS=True)` and `discriminant`
use a checked native Brown–Traub PRS over rational, finite-field, supported
algebraic and rational-function coefficients, including polynomial coefficient
rings. Multivariate elimination removes only the first main generator and
preserves the remaining generators and coefficient-domain metadata. PRS terms
retain their exact normalization, including abnormal degree drops; the last
constant PRS term need not equal the resultant. Resultants preserve the
mathematical sign for the original argument order.

Algebraic coefficient fields execute exact arithmetic in
`QQ[T]/(m(T))`, with the chosen radical embedding retained. `extension=sqrt(2)`,
`extension=I`, `extension=True` and domains such as `"QQ<sqrt(2)>"` select
these fields. `Poly`, Gröbner bases, reduction and order conversion use native
algebraic coefficients. Univariate polynomial operations include exact GCD,
factorization, square-free decomposition, resultants and discriminants.
`Poly.decompose()` also uses exact field coefficients for univariate functional
composition. Rational-function fields such as `QQ(t)`, polynomial coefficient
rings such as `QQ[t]` and `ZZ[t]`, and the supported algebraic number fields
retain their domain and generators in every returned component. Normalized
components must stay in the selected coefficient ring; a `ZZ` polynomial is
not silently decomposed into polynomials over `QQ`.

```python
from symox import Poly, symbols

x, t = symbols('x t')
polynomial = Poly((x*x + x/t)**3 + x*x + x/t, x, domain='QQ(t)')
outer, inner = polynomial.decompose()
assert outer.compose(inner) == polynomial
```

The global `decompose` follows the same domain, extension and modulus options,
and accepts `polys=True`. Prime fields use modular arithmetic; a required
division by zero in the field raises `NotInvertible`. Multivariate functional
decomposition remains unsupported. Degree above 4096 is rejected before dense
coefficient allocation.

The supported constructions also include positive rational cube roots,
rational transformations of generators, composita of up to three independent
quadratic fields, and a positive rational cube root combined with a quadratic
field. Examples include `sqrt(2) + sqrt(3) + sqrt(5)` (degree 8) and
`2**Rational(1, 3) + sqrt(3)` (degree 6). Minimal polynomials and constituent
embeddings are certified with exact rational arithmetic. `minimal_polynomial`
(`minpoly`) and `primitive_element` expose these constructions; `Poly.norm()`
uses a resultant over the field's minimal polynomial. `AlgebraicNumber` retains
the chosen root, exact rational power-basis coefficients and optional alias;
`to_number_field` constructs primitive generators or tests exact membership in
a specified supported field. Arbitrary algebraic
towers and mixtures with transcendental parameters remain unsupported.
`ground_roots` extracts ground-field roots with exact multiplicities, including
roots in the associated field of a coefficient ring. `Poly.factor_list_include()`
and `sqf_list_include()` preserve typed domains and certify reconstruction.
Units are absorbed only into factors of multiplicity one; a separate constant
factor avoids SymPy 1.14's repeated-unit reconstruction bug. These checked
operations bound input and intermediate degrees, coefficient sizes, storage and
work, including finite-field splitting, Hensel lifting and algebraic norm
shifts. Unsupported domains or exhausted limits raise explicit errors.
Typed `factor_list`, `sqf_list` and `sqf_part`, both global and on `Poly`, apply
the selected domain and validate options, including `polys=` result types.
They preserve field versus coefficient-ring normalization and correctly remove
inseparable multiplicities in finite fields. A preflight bounds expression
expansion before construction; exact certificates check the resulting factors.
The natural global expression form keeps its rational scalar content outside the
polynomial factors. Global `factor_list(..., frac=True)` and
`sqf_list(..., frac=True)` return `(coefficient, numerator_factors,
denominator_factors)`, with positive multiplicities on both sides. A bounded
DAG splitter uses structural common denominators without polynomial
cancellation: common numerator/denominator factors remain visible. Domain
selection preserves typed descriptors; splitting, both factorizations and their
exact reconstruction certificates share one native work allowance. Python
coercion and rendering run outside the native cancellation scope. A denominator
that becomes zero after coefficient coercion raises an error; finite-field
units are reduced in their field.
In particular, an unreduced rational unit from a SymPy finite-field list is
returned as its canonical field residue instead.
In fraction mode, numerical algebraic units such as `I` and `sqrt(2)` are
extracted outside natural expression factors. A `Poly` input remains a formal
numerator.
Natural fractional lists with noninteger or symbolic power multiplicities
remain explicitly unsupported; formal polynomial generators can be carried
by a `Poly` input.
For already factored expressions, `sqf_list` computes multiplicities over the
whole polynomial, so repeated factors cannot remain inside a multiplicity-one
component.
This also avoids SymPy 1.14's incorrect content scaling for some square-free
lists over rational polynomial coefficient rings.
`gff_list` and `Poly.gff_list()` compute monic greatest factorial factorizations
over supported characteristic-zero univariate domains. A pair `(h, k)` represents
`h(x)*h(x-1)*...*h(x-k+1)`; the product reconstructs the input divided by its
leading coefficient. Exact reconstruction, maximality and preservation of the
original typed domain are checked within the same resource budget. Zero inputs,
finite fields and multivariate inputs are rejected explicitly.
Norms use the monic defining polynomial, so they equal the determinant of
multiplication over `QQ`. SymPy's `Poly.norm()` can differ by a rational factor
when its integral defining polynomial is not monic.
Univariate `cancel`, `apart`, `invert`, `gcdex` and `half_gcdex` also operate over
these fields. Partial fractions handle repeated and irreducible factors using
exact CRT arithmetic; `full=True` and symbolic parameters in algebraic
coefficients remain unsupported. `Poly` exposes the corresponding Bézout,
inverse and cancellation methods over `QQ` and supported algebraic fields,
plus `coeff_monomial` for the existing coefficient domains. A noninvertible
residue class raises the native `NotInvertible` exception.

`Poly.deflate()`, `terms_gcd()`, `exclude()` and `reorder()` transform sparse
exponent vectors without rebuilding coefficients. Inverse exponent certificates
preserve the polynomial, generator identities and typed coefficient domains,
including unused parameters and algebraic aliases. `reorder()` accepts explicit
permutations or the `wrt`/`sort` priorities; duplicate generators are rejected
instead of reproducing SymPy 1.14's expression-changing result.
`homogeneous_order()` and `has_only_gens()` inspect main-generator support only.
Univariate `all_monoms()`, `all_terms()` and `all_coeffs()` share a linear dense
view capped at 100,000 entries. Sparse transforms allow large exponents but
bound exponent slots and copied coefficient storage. Default generator sorting
also bounds label work; matrix-backed atoms require an explicit permutation.

Poly also provides SymPy's structural methods (`EC`/`EM`/`ET`, ordered
`LC`/`LM`/`LT`, ground arithmetic, `inject`/`eject`, `lift`, `retract`,
`set_modulus`, `rat_clear_denoms`, `sqf_norm`, `revert`, `dispersion`,
`nth_power_roots_poly`, `real_roots`/`all_roots`, `homogenize`, `slice`,
`ltrim`, ...) and the Basic protocol (`atoms`, `has`, `subs`, `xreplace`,
`match`, `compare`, `sort_key`, ...). The top-level `poly_from_expr`,
`parallel_poly_from_expr`, `degree_list`, `total_degree`, `trunc`, `monic`,
`content`, `sqf_norm`, `all_roots` and `nth_power_roots_poly`, the polyfuncs
`horner`, `interpolate`, `rational_interpolate`, `symmetrize` and `viete`, the
special and orthogonal polynomials (`cyclotomic_poly`, `swinnerton_dyer_poly`,
`symmetric_poly`, `interpolating_poly`, `random_poly`, Chebyshev, Hermite,
Jacobi, Laguerre, Legendre), `itermonomials`, `Monomial`, the monomial orders
and `apart_list`/`assemble_partfrac_list` follow SymPy 1.14, together with the
full `polyerrors` hierarchy. The orthogonal recurrences (including
`spherical_bessel_fn`), the `apart_list` full decomposition over `ZZ`/`QQ`,
the dispersion set, `is_cyclotomic` and SymPy's `Poly.cancel` sign/content
normalization run natively. Poly also carries `which_real_roots`,
`which_all_roots`, `make_monic_over_integers_by_scaling_roots`, `refine`,
`count`, `has_free`, `is_same`, `matches` and the rest of the Basic protocol.
`all_roots` lists real roots first.
`galois_group(f, *gens)` and `Poly.galois_group()` classify irreducible
univariate polynomials over `ZZ` and `QQ` through degree four using exact
discriminants, cubic resolvents and native quadratic-field factorization.
They return `(PermutationGroup, alternating)` or a transitive-subgroup enum
with `by_name=True`. The group results expose generators, exact order,
membership, transitivity, commutativity and solvability. Permutations support
array construction, composition, inversion and integer powers. This is a
bounded small-group API, not the full SymPy combinatorics module. Degrees five
and six raise `NotImplementedError`; they are supported by
[SymPy's Galois algorithms](https://docs.sympy.org/latest/modules/polys/numberfields.html).
Other limitations:
`set_modulus` accepts primes only; `apart_list` rejects parametric
(`ZZ[a]`) and floating coefficient domains; `assemble_partfrac_list` expands
only linear root factors (no `RootSum`); `jacobi_poly` with algebraic
parameters returns an equal but differently arranged `EX` expression;
`which_real_roots` keeps exact non-real candidates where SymPy 1.14
overflows; `gff` raises `NotImplementedError` like SymPy; `homogenize` with a
generator already present keeps every term instead of reproducing SymPy
1.14's overwritten coefficient.

```python
theta = sm.sqrt(2) + sm.sqrt(3)
a = sm.to_number_field(sm.sqrt(2), theta)
a.coeffs()                      # [1/2, 0, -9/2, 0] in descending powers of theta
a.as_expr()                     # represents sqrt(2) exactly through powers of theta
a.minpoly_of_element()          # Poly(_x**2 - 2, _x, domain='QQ')
```

`sturm`, `count_roots`, `intervals` and `refine_root` provide certified real-root
analysis for univariate polynomials over `QQ` and `ZZ`, including their `Poly`
methods. Counts use distinct roots and closed rational bounds; isolating
intervals retain multiplicities and exact rational roots. Refinement preserves
the isolated root and applies `eps`, `steps`, `fast` and `check_sqf`.
The arithmetic stays exact and sparse, with a degree bound of 4096 and shared
work, storage and rational-size limits. Sturm sequences promote `ZZ` to `QQ`
unless `auto=False`. Isolating fractions can differ from SymPy's continued
fraction algorithm. Complex rectangles, other coefficient domains,
`intervals(all=True)`, `strict=True` and simultaneous polynomial lists are
explicitly unsupported. Reversed bounds and singleton intervals that are not
roots raise errors.

`linear_eq_to_matrix` extracts exact linear coefficients; `linsolve` accepts
equations with explicit unknowns, an augmented matrix, or an `(A, b)` matrix
pair. Underdetermined solutions retain free variables or receive generated
parameters when none are supplied. `solve_linear_system` provides the
dictionary interface. Nonlinear terms raise `NonlinearError`. Symbolic row
reduction uses generic nonzero pivots; it does not enumerate exceptional
parameter values.

`solve_poly_system` uses Gröbner elimination and exact back substitution for
supported finite polynomial systems. `nonlinsolve` also decomposes supported
factorable and triangular rational systems with free parameters. Elimination
retains both zero and nonzero coefficient branches; pole exclusions are
represented with `Complement(Complexes, FiniteSet(...))`. All branches must be
resolved before a result is returned. External parameter conditions, unresolved
components and roots outside the exact algorithms raise `NotImplementedError`;
no partial list is presented as a complete solution.
Finite radical roots and residuals require exact zero certificates. Supported
closed radical identities are proved by polynomial ideal membership, including
nonvanishing proofs for their denominators; unresolved certificates stay
explicitly unsupported.

`solve(..., set=True)` infers omitted unknowns and returns the symbols and
solution tuples. Candidate checking respects known symbol assumptions;
`check=False` retains candidates before pole, residual and assumption filtering.
`nsolve`
supports real and complex Newton iteration, arbitrary precision through `prec`,
real bracketed bisection, `tol`, `maxsteps`, `verify` and dictionary output.
Overdetermined systems use complex QR least squares in the arbitrary-precision
solver. It verifies the final residual by default and rejects failed convergence.
These solvers execute natively without SymPy or mpmath. Unsupported numerical
solver selections are explicit errors.
Unsupported `solve` flags, ODE hints and integration options are rejected
instead of being silently ignored.

`residue(expr, x, point)` extracts the coefficient of the inverse first power
from the native Laurent expansion, including higher-order poles. Unsupported
branch points and essential singularities raise `NotImplementedError`.
Rational definite integrals cancel removable factors and analyze one-sided
divergences at real poles and infinite tails. Exact Sturm root counts certify
pole coverage; unresolved cases remain unevaluated. Opposing infinities produce
`nan`, rather than a finite Cauchy principal value.
