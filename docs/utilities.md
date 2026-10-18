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
and then `default_sort_key`. `warn=True` rejects unresolved ties between distinct
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
