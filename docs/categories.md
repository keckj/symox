# Categories and finite diagrams

[README](../README.md) · [Printing](printing.md) · [Lie algebras](lie-algebras.md)

`symox.categories` provides all eleven exports of SymPy 1.14's corresponding
package. The supported structural operations are described below; this is not
a claim of full symbolic-object or layout compatibility. The public interfaces
follow the [SymPy category documentation](https://docs.sympy.org/latest/modules/categories.html).
Objects and morphisms are immutable, hashable Python values with the
reference sort keys. Diagrams use Python frozensets and read-only mappings
rather than native scalar `FiniteSet` and `Dict` containers; these sets print
like SymPy's `FiniteSet`/`EmptySet` and are ordered like the reference `Dict`.
No SymPy runtime dependency is required.

## Objects, morphisms and composition

```python
from symox.categories import Object, NamedMorphism, IdentityMorphism
from symox.categories import CompositeMorphism, Diagram, Category

A, B, C = map(Object, 'ABC')
f = NamedMorphism(A, B, 'f')
g = NamedMorphism(B, C, 'g')
assert (g*f).domain == A
assert (g*f).codomain == C
assert CompositeMorphism(f, g) == g*f
assert IdentityMorphism(B)*f == f
```

`NamedMorphism(domain, codomain, name)` is a labelled arrow. As in the
reference, the base class `Morphism(domain, codomain)` cannot be instantiated
and raises `NotImplementedError`. Endpoints must be category `Object` instances. `compose` and `*` compose in categorical order:
`g*f` applies `f` first. `CompositeMorphism` constructor arguments instead use
path order `(f, g)`. Nested compositions flatten, compatible identity arrows
cancel, and incompatible endpoints raise an error. `flatten(name)` replaces
a composite with a named arrow having the same endpoints. Morphism equality
is structural; a composite and a separately named arrow are not automatically
identified by an asserted commutative diagram.

`Category(name, objects=(), commutative_diagrams=())` records an abstract named
category, its supplied object set and asserted diagrams. As in the abstract
reference interface, `hom` and `all_morphisms` raise `NotImplementedError`:
an arbitrary abstract category does not specify enumerable hom sets. Functors,
natural transformations, limits and categorical equivalence inference are not
implemented.

## Finite diagrams and properties

`Diagram(premises=(), conclusions=())` accepts lists of morphisms or mappings
from morphisms to properties. A property string is one property, not a sequence
of characters. Diagram premises include endpoint identities, components of
explicit composites and composable simple paths. Path properties are the
intersection of their component properties. Identity arrows cannot carry
additional properties. Conclusions on objects outside the premise object set
are excluded, following the reference convention.

```python
d = Diagram({f: ['mono', 'epi'], g: ['mono']}, {g*f: 'unique'})
assert d.premises[g*f] == frozenset({'mono'})
assert d.conclusions[g*f] == frozenset({'unique'})
assert g*f in d.hom(A, C)[0]
K = Category('K', [A, B, C], [d])
assert d in K.commutative_diagrams
```

Diagrams remain finite in the presence of cycles: explicitly supplied cyclic
words and loop arrows are retained, but their unbounded powers are not
enumerated. Closure generates simple paths without revisiting objects; it
does not construct an infinite free path category. A 10,000 generated-path
limit raises an error rather than silently truncating. Premise order does not
change the acyclic path set.

`hom(A, B)` returns premise and conclusion morphism sets with those endpoints.
`d.is_subdiagram(other)` follows the reference's direction: it asks whether
`other` is contained in `d`, with matching properties. `subdiagram_from_objects`
retains morphisms whose endpoints belong to the selected subset, including
composites whose internal factor objects need not belong to that subset.

## Layout and rendering

`DiagramGrid`, `ArrowStringDescription`, `XypicDiagramDrawer` and
`xypic_draw_diagram` are a port of SymPy 1.14's `categories/diagram_drawing.py`,
and the 153 reference doctests of that module pass. `DiagramGrid(diagram,
groups=None, **hints)` uses the reference generic layout (skeleton triangles
welded onto a growing grid), lays out disconnected components side by side,
and accepts the `layout='sequential'` and `transpose=True` hints. Groups lay
out nested subdiagrams as in the reference. The reference errors are kept:
`IndexError` for an empty diagram and `KeyError` for groups that miss an
object. Grid dimensions, indexed cells, the displayed morphism mapping and
pretty printing of the grid are available.

`XypicDiagramDrawer.draw` and `xypic_draw_diagram` produce executable
`\xymatrix` markup, with native LaTeX object/arrow labels, composite labels,
loops, parallel-arrow curves, masks and a supplied matrix format string.
Property callbacks in `drawer.arrow_formatters` and a `default_arrow_formatter`
can change the arrow description's label, position, curve and style.

```python
from symox.categories import xypic_draw_diagram

source = xypic_draw_diagram(d, diagram_format='@R=1cm')
assert source.startswith(r'\xymatrix@R=1cm')
assert r'\circ' in source
```

`preview_diagram` passes this markup to the native printing preview pipeline,
loading the Xy-pic extensions needed for matrices, curves and loops
(`\xyoption{all}`, which SymPy omits, so its LaTeX run fails on curved arrows
and loops). It accepts the output/viewer options and forwards `filename`,
`outputbuffer` and `dvioptions` to the preview, where SymPy hands them to
`DiagramGrid` as ignored hints. For example, `viewer='file',
filename='diagram.png'` writes a PNG when TeX, Xy-pic and the conversion
executable are installed. Rendering tools are optional system dependencies;
no GUI is required for file output.

Tests cover associativity and identities, endpoint validation, diagram/path
properties, selected-object subdiagrams, finite loops, masks and arrow
formatters, and compare grids and Xy-pic markup with SymPy 1.14 output. When
TeX is available, a diagram with parallel arrows and a loop is actually
compiled to PNG. A separate subprocess renders diagram markup with all SymPy
imports blocked.
