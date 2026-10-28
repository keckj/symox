# CSE attribute-protocol follow-up (Python only)

Base is the frozen21-file BigInt derivative package at `../symbolic-derivative-bigint-r1b`, not main. No Rust file or binary changed. Native release remains b83cdadd3d15aefd1bd3c290f32642a59bf99e72466844a9e207ed3d3e11acd7.

An independent post-freeze probe showed the stock-constructor guard still accepted an Integer base with a customized `__getattribute__`: the native reconstruction observed0 callbacks, whereas the established CSE reconstruction observed2 (`_symox_matrixexpr`, `_symox_sparse`). The same concern applies to the Pow's `func` attribute and absent base attributes handled by `__getattr__`.

The guard now statically checks the MRO dictionaries of both operand types. It declines if their effective `__getattribute__` differs from object.__getattribute__ or if any MRO class defines `__getattr__`. This inspection does not call user descriptors. Declining runs the exact established func/constructor route; no mathematical algorithm or domain was changed.

Three regression cases customize the root's getattribute, the base's getattribute and the base's getattr. Each compares the complete returned CSE result and observed callbacks against the established forced fallback. All1171 focused release tests pass locally (16.43s) and on Demeter (18.96s), including the original673 new BigInt/recipe cases, existing derivative/parser/host/query tests and the entire Holo regression file. The preceding frozen package's broader10539-test release gate,689 Rust tests+5doc tests and clippy gates remain separate evidence.

Exact bases/finals and patch digest are in package.json. The patch applies cleanly to the delivered original stage. Both this follow-up and the original package are now source-frozen; proof logs may be copied independently.
