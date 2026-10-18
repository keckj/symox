# Polynomial rings, rational function fields and number fields

The constructors and arithmetic below use native symox kernels. SymPy is only
used by differential tests and explicit `to_sympy()` conversions.

## Sparse ring constructors

```python
from symox import ring, field, QQ, ZZ

R, x, y = ring('x,y', QQ, order='grevlex')
f = (x+y)*(x-y)
f.exquo(x-y)                 # x+y
f.gcd((x+y)*(y+1))           # x+y
q, r = f.div([x, y])
assert q[0]*x + q[1]*y + r == f

F, u, v = field('u,v', ZZ)
g = (u*u-v*v)/(u-v)          # u+v, with exact cancellation
assert g.denom == 1
```

`ring` and `field` return the parent followed by its generators. `xring` and
`xfield` instead return `(parent, generators_tuple)`. `vring` and `vfield`
inject generators into the calling module's globals and return the parent.
This namespace mutation occurs only when these explicit constructors are used.

`sring(expressions, *generators, **options)` and `sfield` infer generators and
coefficient domains. A scalar produces one element; a sequence produces a list.
Supported options are `domain`, `order`, `field` and `extension`, subject to the
underlying domain restrictions. Empty generator lists represent constant rings.

Parents expose `symbols`, `gens`, `ngens`, `domain`, `order`, `zero`, `one`,
`from_expr`, `is_element` and conversions between ring and field parents.
Polynomial elements expose sparse coefficient dictionaries, ordered terms,
degrees, leading coefficients/monomials, exact arithmetic, division with
remainder, GCD, LCM, derivatives, evaluation and generator dropping.
Fraction elements expose reduced `numer` and `denom`, arithmetic including
negative integer powers, derivatives and evaluation. Pickling preserves parents.

The current native composite-domain backend supports ground domains `ZZ` and
`QQ`, with `lex`, `grlex` and `grevlex`. Finite, algebraic, numerical and nested
ground domains, local orders and arbitrary callable orders are not implemented
by these constructors. They fail explicitly rather than converting coefficients
to a different domain. This restriction does not remove the broader coefficient
domains supported by `Poly`; see [polynomials](polynomials.md).

## Embeddings and real isolation

`field_isomorphism(a, b, fast=True)` returns rational coefficients in falling
powers of the target generator, or `None` when there is no embedding in the
supported native field model. AlgebraicNumber inputs are converted to primitive
elements. The `fast=False` option uses the same exact native embedding engine.
Existing limits on radical towers and algebraic field models still apply.

`isolate(a, eps=None, fast=False)` returns an interval with rational endpoints
containing exactly the selected real algebraic root. Native Sturm isolation
supplies the intervals; exact sign queries certify the selected embedding.
Rational inputs return `(a, a)`. `eps` must be a positive rational. `fast` is
accepted without changing the exact isolation algorithm. Complex inputs and
embeddings whose signs cannot be certified raise `NotImplementedError`.

## Maximal orders and prime ideals

```python
from symox import symbols, Poly, round_two, prime_decomp, prime_valuation

x = symbols('x')
T = Poly(x*x-5, x)
ZK, discriminant = round_two(T)
assert discriminant == 5
for P in prime_decomp(11, ZK=ZK):
    assert (P.e, P.f) == (1, 1)
    assert prime_valuation(P**3, P) == 3
```

`round_two` provides the maximal-order result through an exact bounded prime
saturation algorithm. It enumerates projective residue classes in `(1/p)L/L`;
every enlargement must have an integral characteristic polynomial. The final
lattice discriminant is checked independently using the trace pairing.

Inputs must be irreducible over `ZZ` or `QQ`, of degree one through eight, and
have integral coefficients after making the polynomial monic. Enumeration is
limited to 100,000 projective classes per prime. This is suitable for small
fields; it is not an implementation of the efficient general Round Two
algorithm. The optional `radicals` cache container is accepted and left unchanged.

Orders expose `parent`, `matrix` (an integer lattice matrix), `denom`,
`QQ_matrix`, `basis_elements` and exact membership/coordinate calculations.
Matrix views provide `to_Matrix()`. Power basis elements expose ascending
integer `coeffs`, a common `denom`, exact products and polynomial conversion.
Full fractional ideal lattices support sum, product, positive powers and rational
scaling, normalized by integer Hermite column operations. These small matrix
views do not implement the full `DomainMatrix` or number-field module APIs.

`prime_decomp` factors all rational primes in quadratic fields, including primes
dividing the power-basis index. In higher degrees it currently requires a prime
coprime to that index. It checks prime norms, `sum(e*f) == degree`, and the ideal
product against `(p)`. A supplied maximal order and discriminant must agree with
the certified values. Precomputed prime radicals are not accepted yet.

Prime objects expose `p`, `e`, `f`, `alpha`, `ZK`, `as_submodule`, multiplication,
nonnegative powers and `valuation`. `prime_valuation` checks fractional ideal
closure and computes the exponent by exact lattice containment, including
negative exponents after rational scaling. Zero-rank ideals, inverse ideal
powers and valuations requiring more than 1024 containment steps are unsupported.
