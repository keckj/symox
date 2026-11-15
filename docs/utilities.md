# Utilities

[README](../README.md)

The following helpers are available from `symox` and `symox.utilities`, without
importing SymPy.

## Containers and iterators

`flatten`, `unflatten` and `reshape` handle nested containers. `flatten`
supports a depth limit or an expression class such as `Add`; ordinary scalar
expressions remain atomic unless a class is supplied. Arrays flatten through
their entries. `reshape` uses a nested template of counts and container types.

`sift` groups values by a key and returns a `defaultdict(list)`. Its binary mode
returns `(true_values, false_values)` and rejects nonbinary keys. `group`
collects adjacent equal values, optionally returning `(value, count)` pairs.

`subsets`, `variations`, `cartes`, `prefixes`, `postfixes` and `rotations` return
iterators. `subsets` and `variations` support repetition. `has_dups` accepts
unhashable values, `has_variety` tests whether values differ, and `take`
consumes at most the requested number of items. `prod` retains ordinary Python
integer arithmetic for integer-only inputs. `dict_merge` applies mappings in
order, with later values taking precedence.

## Ordering and traversal

`ordered(seq, keys=None, default=True, warn=False)` applies successive sort keys
only to groups still tied by earlier keys. Its default keys compare node counts
and then `default_sort_key`. Node counts follow SymPy's expanded-tree count
(a `Float` counts one half); a shared subexpression is counted once per call
and its count reused. `warn=True` rejects unresolved ties between distinct
values when no more keys are available. Generic sort keys support strings and
nested lists, tuples, dictionaries and sets as well as native scalar expressions.
Explicit monomial orders propagate to native expression children, including
expressions nested in containers. Native symbolic sets retain their symbolic
class key rather than being treated as ordinary iterable containers.

`preorder_traversal` visits parents before children. Calling its `skip()` method
omits the last yielded node's descendants. `postorder_traversal` visits children
before parents. Both support deterministic child ordering via `keys=True` or
custom ordering keys, and use an explicit stack rather than Python recursion.
They traverse the native public arguments, including Python containers used
for constructor metadata. Repeated expression occurrences remain repeated.

`bottom_up` transforms expression children before their parents; `atoms=True`
also transforms leaves. `nonbasic=True` attempts the callback on ordinary Python
objects. `use(expr, func, level=n)` applies a transformation only at depth `n`.

```python
from symox import symbols, ordered, preorder_traversal, use, expand

x, y = symbols('x y')
assert list(ordered([y, x])) == [x, y]
expression = (x+y)**2*x+1
assert use(expression, expand, level=2) == x*(x*x+2*x*y+y*y)+1
```

`topological_sort((vertices, edges), key=None)` orders a finite directed acyclic
graph using a priority queue. Edges point from a predecessor to its successor.
The key orders currently available vertices; cycles and undeclared edge endpoints
raise `ValueError`.

## Memoization

`cacheit` applies a typed LRU cache with at most 1,000 entries per decorated
function. Unhashable inputs execute without caching. Function exceptions are
propagated once and are not cached. Decorated functions expose `cache_info()`
and `cache_clear()`; `symox.utilities.cache.clear_cache()` clears registered
caches. Callers should cache immutable results.

`memoize_property` computes a property once per instance and stores its value in
`_<property_name>`. This also supports slotted classes that declare that slot.

## Combinatorial enumeration and infinite iterators

The complete measured iterable surface is exposed by
`symox.utilities.iterables` and `symox.utilities`. Multiset combinations,
permutations, partitions and derangements enumerate distinct results in a
canonical order. Integer `partitions` yields multiplicity dictionaries;
`ordered_partitions` yields nondecreasing lists, not compositions.
`binary_partitions` restricts parts to powers of two. Several enumeration
algorithms reuse their yielded lists or dictionaries: copy a result before
retaining it while advancing the iterator.

`generate_bell` enumerates permutations in adjacent-transposition order;
`generate_involutions` selects self-inverse permutations. `necklaces`,
`bracelets`, `generate_oriented_forest`, `kbins`, `permute_signs` and
`signed_permutations` provide finite combinatorial generators. These outputs
can grow exponentially or factorially; the functions do not impose a size cap.

`random_derangement` uses a bounded randomized group rotation, preserving the
input multiset. If a derangement is impossible it raises `ValueError`, or with
`strict=False` minimizes the number of fixed points. Sampling is not uniform
over all possible derangements. A custom `choice` function controls its
randomness without introducing a rejection loop.

`iproduct` fairly enumerates Cartesian products of finite or infinite iterators.
It caches prefixes of its inputs, so infinite inputs can require unbounded
memory. An empty factor immediately gives an empty product.
`roundrobin` interleaves inputs until all are exhausted.

`strongly_connected_components` returns directed graph components in reverse
topological order; `connected_components` treats edges as undirected. Both use
explicit stacks rather than recursive graph traversal.

`ibin`, `multiset`, `common_prefix`, `common_suffix`, rotations, `least_rotation`,
`minlex`, `runs`, `sequence_partitions`, `sequence_partitions_empty`, `uniq`,
`iterable`, `is_sequence` and `is_palindromic` handle ordinary containers.
`NotIterable` marks objects that should be excluded by `iterable`.
`numbered_symbols` generates native symbols, honoring assumptions and excluded
symbols; `filter_symbols` removes excluded or already seen values. `capture`
collects standard output and restores it if the callable raises.
`interactive_traversal` lets a terminal user select expression subtrees.

## Compatibility submodules

`symox.utilities` also exports `filldedent`, `lambdify`, `public`, `threaded`,
`xthreaded` and `timed`. The SymPy import paths `symox.utilities.misc`,
`decorator`, `exceptions`, `timeutils`, `source`, `magic` and `enumerative`
provide the corresponding helpers, such as `ordinal`, `deprecated`,
`SymPyDeprecationWarning`, `timethis` and `multiset_partitions_taocp`.
`symox.utilities.memoization` provides `recurrence_memo` and
`assoc_recurrence_memo`. As in SymPy, an index below the cache length,
including a negative one, reads the cache directly, so negative indices count
from the end of the computed cache, and an index outside it raises
`IndexError`.
`symox.utilities.lambdify`, `codegen` and `autowrap` are described with
[code generation](codegen/overview.md) and the
[API implementation notes](api-implementation.md).
