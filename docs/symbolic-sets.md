# Structural symbolic sets

Native scalar sets (`Interval`, `FiniteSet`, `Union`, `Intersection`,
`Complement`, `ProductSet`, `Range`, `ImageSet`, `ConditionSet` and named sets)
are described in [compatibility](compatibility.md). The additional structures
below evaluate membership using native scalar predicates.

## Images, power sets and disjoint unions

`imageset` accepts a native `Lambda`, a callable, or variables followed by an
expression and its base sets. It follows SymPy's `imageset` and
`_set_function` handlers: images of `FiniteSet`, `Interval` and `Reals`
(from singularities, limits and critical points), `Range`, `Integers`,
`Naturals`, `Union` and invertible `Intersection` bases are evaluated, and
nested image sets are composed. Several base sets give an unevaluated
`ImageSet`, and other images are represented by native `ImageSet` nodes.

`PowerSet(A)` preserves a symbolic power-set object. Membership tests whether
the proposed element is a subset of `A`. Finite bases support iteration,
cardinality and `rewrite(FiniteSet)`. Enumeration is limited to 16 base elements
(65,536 subsets). Infinite power-set iteration and undecidable subset queries
raise `NotImplementedError`.

`DisjointUnion(A, B, ...)` tags each element with its component index: an element
is the tuple `(value, index)`. It provides membership, finite cardinality,
round-robin iteration across iterable components, and `rewrite(Union)` into
Cartesian products. Equal values in different components remain distinct.

Power sets, disjoint unions and complex regions are structural Python sets,
not native scalar AST nodes. Their membership methods, stated rewrites,
substitution and explicit `to_sympy()` conversion are supported. General set
algebra mixing these structures with native sets, native scalar coercion,
and the wider SymPy Set method surface require further implementation.

## Complex regions

`ComplexRegion(ProductSet(A, B), polar=False)` interprets pairs as `a + I*b`;
as in SymPy, the constructor returns a `CartesianComplexRegion` or a
`PolarComplexRegion` (`Reals*Reals` gives `Complexes`, a product of two
finite sets gives a `FiniteSet`). Unions of two-factor products are also
accepted. `a_interval`, `b_interval`, `psets`, `from_real`, membership,
finite enumeration and rectangular area are available. `union` and
`intersect` follow SymPy's rules for two regions of the same kind and for
real sets. Area uses inclusion-exclusion, so overlapping rectangles are
counted once, with at most 12 component rectangles.

With `polar=True`, pairs represent `r*exp(I*theta)`. Angles given as
intervals or finite sets are normalized into `[0, 2*pi)` by SymPy's
`normalize_theta_set` (`symox.sets.fancysets`); radii must be nonnegative.
Other angle sets and symbolic angles whose sign cannot be certified raise
`NotImplementedError`. Polar interval regions also support area using the
radial Jacobian. General region topology and implicit complex regions are
not implemented.

## Ordinals

`Ordinal`, `OmegaPower`, `ord0` and `symox.sets.ordinals.omega` represent finite
Cantor normal forms below epsilon-zero. For example:

```python
from symox.sets.ordinals import omega as w

assert 3+w == w
assert (w+1)*w == w**2
assert w*(w+1) == w**2+w
```

Terms have strictly decreasing ordinal exponents and positive integer
multiplicities. Addition and multiplication respect ordinal order; comparison,
finite ordinal conversion, limit/successor predicates, finite powers and
`omega**ordinal` are supported. These are structural ordinal values, separate
from native scalar expressions. General transfinite exponentiation with bases
other than `omega` is not implemented.
