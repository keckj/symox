# Polynomial rings, rational function fields and number fields

The ring and field constructors and their arithmetic use native symox kernels.
The number-field module layer is a Python port of SymPy's algorithms over
native domains, with its Hermite normal forms in Rust. SymPy is only used by
differential tests and explicit `to_sympy()` conversions.

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

The native composite-domain backend accepts the ground domains `ZZ`, `QQ`,
`ZZ_I`, `QQ_I`, `RR`/`CC` (including other precisions), polynomial rings and
rational-function fields (`ZZ[y]`, `QQ(y)`), and, for rings only, supported
algebraic fields such as `QQ<sqrt(2)>`, with `lex`, `grlex` and `grevlex`.
Finite fields and `EX` grounds, fraction fields over algebraic grounds, local
orders and arbitrary callable orders are not implemented by these
constructors. They fail explicitly with `NotImplementedError` rather than
converting coefficients to a different domain. This restriction does not
remove the broader coefficient domains supported by `Poly`; see
[polynomials](polynomials.md) and [coefficient domains](domains.md).

## Embeddings and real isolation

`field_isomorphism(a, b, fast=True)` returns rational coefficients in falling
powers of the target generator, or `None` when there is no embedding in the
native field model. AlgebraicNumber inputs are converted to primitive
elements. The `fast=False` option uses the same exact native embedding engine
(`to_number_field`), which handles the general algebraic fields described in
[coefficient domains](domains.md#algebraic-field-constructors).

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

`symox.polys.numberfields.modules`, `basis` (`round_two`), `primes`
(`prime_decomp`, `prime_valuation`, `PrimeIdeal`), `utilities` and `exceptions`
are ports of SymPy's pure-Python modules running on native symox domains and
`DomainMatrix`. `round_two` is the reference Round Two algorithm (Dedekind
criterion and `p`-radical enlargement); the optional `radicals` dictionary
receives the computed nilradicals, and `prime_decomp(p, T, ZK, dK, radical)`
accepts them back. Integer Hermite normal forms (ordinary and modulo `D`), used
by `round_two`, `prime_decomp` and the modules, run on the native Rust HNF
kernels. As in SymPy, `round_two` can raise `ClosureFailure` on some inputs.

The reference module API is available: `PowerBasis`, `Submodule`, module
elements built from columns, `element_from_rational`, `represent`,
`find_min_poly`, and module homomorphisms with `kernel(modulus)`. Orders expose `matrix` as a
`DomainMatrix` over `ZZ`, `denom`, `QQ_matrix` and `basis_elements`.
`AlgebraicField` also provides `maximal_order`, `integral_basis(fmt)`,
`discriminant` and `primes_above` through these modules.

`prime_valuation` follows the reference loop (Cohen 4.8.17). Where SymPy 1.14
reconverts a non-integral intermediate matrix to `ZZ` and raises
`CoercionFailed` (for example `2*ZK` in `QQ(sqrt(2))`), symox instead returns
the valuation given by exact ideal containment, bounded by the norm. A
symbolic generator in `QQ.algebraic_field((T, theta))` is accepted for
`prime_decomp` and `primes_above`; a reducible `T` raises `NotImplementedError`.
