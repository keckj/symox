# Indexed tensor algebra

[README](../README.md) · [Concrete arrays](arrays.md)

`symox.tensor.tensor` provides typed abstract indices, tensor heads, Einstein
products and sums, signed slot symmetries, metric contraction, component
selection and wildcard matching. Scalar coefficients, matrix inverses and
array operations use native symox algorithms. Calculations do not import
SymPy; explicit `to_sympy()` conversion is optional.

## Indices and expressions

```python
from symox import Array, diag
from symox.tensor.tensor import (
    TensorIndexType, TensorIndex, TensorHead, TensorSymmetry,
    tensor_indices, tensor_heads, canon_bp, contract_metric,
)

E = TensorIndexType('E', dim=3)
i, j, k, l = tensor_indices('i j k l', E)
v, w = tensor_heads('v w', [E])
S = TensorHead('S', [E, E], TensorSymmetry.fully_symmetric(2))
F = TensorHead('F', [E, E], TensorSymmetry.fully_symmetric(-2))

assert canon_bp(S(i, j)*F(-i, -j)) == 0
assert contract_metric(E.metric(i, j)*v(-j), E.metric) == v(i)
```

Indices with `is_up=True` are contravariant; negating an index changes its
variance. An Einstein contraction pairs two occurrences of the same typed
index with opposite variance. Repeated indices of equal variance and more
than two occurrences in one product are rejected. Different index types cannot
be substituted for one another, even if their dimensions agree.

Tensor application creates a `Tensor`; products and sums create `TensMul` and
`TensAdd`. The abstract `TensExpr` base requires concrete subclasses for
evaluation. Sums require identical free index sets. Products distribute over
sums eagerly, and scalar coefficients combine natively. Ordinary scalar
arithmetic and `Add`, `Mul`, `simplify` and `expand` dispatch to the tensor
layer when indexed expressions are present.

`get_indices()` reports slot order, `get_free_indices()` preserves the order
of uncontracted slots and `free_args` sorts free indices by type and name.
`rank` counts free indices; `ext_rank` counts all slots. `free` and `dum`
describe a monomial's free and contracted positions. `coeff`, `nocoeff`,
`components`, `index_types`, `types`, structural replacement, hashing and
pickling are available on the corresponding concrete classes.

Bound dummy names are generated from the index type's `dummy_name`, avoiding
free names. Multiplying expressions that already contain contractions first
renames their bound indices, so scalar contractions cannot capture another
factor's free indices. `substitute_indices(expr, (old, new), ...)` makes typed,
simultaneous substitutions and avoids capture by existing dummies.

## Symmetries and canonical forms

`TensorSymmetry` accepts a base and signed native `Permutation` generators.
The last two permutation points encode the sign. Generators can exchange only
slots with the same index type. Convenience constructors are:

- `fully_symmetric(n)`: symmetric slots, or antisymmetric slots when `n < 0`.
- `direct_product(*blocks)`: independent symmetric/antisymmetric slot blocks.
- `no_symmetry(n)`: an unconstrained rank.
- `riemann()`: the rank-four monoterm Riemann symmetries.

`canon_bp(expr)` and `expr.canon_bp()` use an exact signed orbit of slot
permutations, permitted component commutations and metric-compatible dummy
variance changes. Opposite signs on one orbit state prove zero. Sums combine
their canonical monomials. The compatible function name does not imply the
large-group Butler-Portugal algorithm: this implementation enumerates a
bounded orbit. Representatives can differ structurally from SymPy's; their
signed tensor identities are preserved.

Commutation groups are `0` (commuting), `1` (anticommuting) and `2`
(noncommuting). Group zero commutes with all groups. Arbitrary user-defined
commutation tables are not implemented.

`riemann_cyclic_replace(R(i,j,k,l))` constructs the standard three-term cyclic
projection. `riemann_cyclic(expr)` applies it to Riemann-symmetric components,
expands products and canonicalizes the result. Monoterm symmetry alone does
not impose the Bianchi identity; the cyclic projection is explicit.

## Metrics and explicit components

`TensorIndexType` supports symbolic or finite dimensions, custom dummy and
metric names, and `metric_symmetry=1`, `-1`, `0` or `None`. `metric`, `delta`
and, for concrete dimensions, `epsilon` create the corresponding heads.
`set_metric` can replace a metric with the same declared slot symmetry.

`contract_metric(expr, g)` eliminates symmetric and antisymmetric metrics,
including traces equal to the dimension with the appropriate variance sign.
Signs are independent of the position of a commuting metric factor.
`expr.contract_delta(delta)` applies the same elimination to a Kronecker head.
Automatic contraction with a nonsymmetric metric is not implemented.

```python
assert (v(i)*w(-i)).replace_with_arrays({
    v(i): [1, 2, 3], w(-i): [4, 5, 6],
}) == 32

assert v(-i).replace_with_arrays({
    v(i): [1, 2, 3], E: diag(1, 2, 3),
}).tolist() == [1, 4, 9]
```

`replace_with_arrays(mapping, indices=None)` consumes full component arrays
keyed by tensor applications. Their index variances specify the supplied
components. A type key supplies its covariant metric; symmetric metrics and
their exact inverses lower and raise indices. Metric and Kronecker tensor
components can be derived from the type's metric and dimension. Explicit
tensor-head keys are also accepted, with contravariant supplied axes.

Contractions use native array tensor products and axis contraction. Optional
`indices` specifies the output axis order and can request opposite variances
when metric data are supplied. Rank and known dimensions are validated;
declared slot symmetries are checked on supplied component axes with the
same variance. Mixed-variance component data remain the caller's tensor
assignment. Anticommuting components require a Grassmann coefficient algebra
and are rejected by the ordinary scalar array evaluator. Automatic numeric
raising/lowering with an antisymmetric or nonsymmetric metric is unsupported;
explicit components already carrying the requested variances can still be used.

`TensorElement(expr, {index: position, ...})` selects fixed components, removes
their free axes and evaluates through `replace_with_arrays`. It accepts typed
indices or their symbol names and validates concrete positions. General
symbolic products of component-selection objects are not implemented.

`PartialDerivative(expr, *variables)` from `symox.tensor.toperators` records
free and contracted indices and evaluates through `replace_with_arrays` from
supplied components; symbolic product algebra with these objects is not fully
adapted. A contracted derivative index is summed with the matching axis of
the components; SymPy 1.14 returns a result that disagrees with its own manual
contraction there (an intentional divergence recorded in
[the parity contract](algorithmic-parity.md)). `Tensor.split()` and
`TensMul.split()` return the factor list, and string printing of tensor
objects follows SymPy.

## Matching and compatibility helpers

`WildTensorIndex` matches an index type and variance, with
`ignore_updown=True` allowing either variance. Existing wildcard bindings and
their negations must remain consistent. `WildTensorHead` creates `WildTensor`
patterns; a call with no indices captures a whole tensor expression.
`unordered_indices=True` permits index-order matching. Captured expressions
can be substituted back into patterns, preserving their index alignment.

Products support typed fixed factors, scalar coefficient patterns and at most
one unrestricted tensor wildcard. Commuting factors use bounded backtracking;
noncommuting factors retain their order, including a prefix or suffix captured
by the wildcard. This is not the complete SymPy tensor pattern engine.

`get_indices`, `get_free_indices`, `get_dummy_indices`, `get_coeff` and
`get_index_structure` expose expression metadata; a sum has no single
monomial index structure. `components_canon_args` and `perm2tensor` provide
signed permutation representations. `get_lines` recognizes directed matrix
lines and traces for an index type; transposed-line flags are unsupported.
`get_postprocessor` dispatches `Add` and `Mul` nodes containing tensors.

Legacy `TensorType`, `tensorhead`, `tensorsymmetry`, expression calls and
`fun_eval` emit Python deprecation warnings and delegate to the supported
constructors or substitutions. The exported `deprecate_call`, `deprecate_data`
and `deprecate_fun_eval` warning functions are also available. Legacy `.data`
storage is not implemented; use explicit component mappings.

## Computational bounds

Expressions have at most 32 slots per monomial, sums and expansions at most
10,000 terms, and canonical symmetry orbits at most 100,000 states. Exceeding
an orbit bound raises `NotImplementedError`, rather than returning an
unverified canonical form. Wildcard factor/sum matching is limited to eight
factors/terms, and unordered matching to eight free indices. Concrete arrays
retain their 100,000-position materialization bounds. These Python structural
objects do not implement SymPy's complete `Basic` hierarchy, arbitrary
Young-projector algebra or a general large-group tensor canonicalizer.
