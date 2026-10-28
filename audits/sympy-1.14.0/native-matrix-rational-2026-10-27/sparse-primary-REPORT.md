# Sparse exact SDM/Berkowitz port (private frozen sources)

Base: immutable R3b `native-inexact-r1-20261027/after`, native SHA256
`8a8a8e126be486e31b83b7cae2c306e54ab21e9d569f01bbe46adeb8bb661f3a`.
The five source paths and all base/final hashes are fixed in
`source-manifest.json`; the patch passes nonmutating `git apply --check` on R3b.
No main edits, shared installation, maturin, commits, runtime SymPy or benchmark.

## Mathematical scope and proof

For stock native ZZ/QQ entries the Rust kernel computes the entire sparse map
of coefficients of det(lambda I - M), indexed from its monic leading term.
Its coefficient arithmetic is exact native Int/Rat, not expression or Python
arithmetic. Other domains and customized protocols retain the existing general
Berkowitz implementation. A declined certificate is None, never an empty or
approximate characteristic polynomial.

The vertices incident to a nonzero entry form an invariant active block after
simultaneously permuting the rows and columns. Every omitted vertex has both
zero row and zero column, so the permutation gives M_active direct-summed
with a zero matrix. Thus det(lambda I-M)=lambda**(n-r)*det(lambda I-M_active).
Leading coefficient indices do not change. A leading zero row alone or column
alone is NOT omitted unless the vertex has neither incoming nor outgoing
nonzeros. The complete coordinate validation precedes ignoring explicit zeros.
The kernel stores O(active vertices + stored entries) matrix/column views,
never an n*n matrix or a declared-n vector. The active block can fill during
arithmetic as mathematically necessary. Output can have n+1 coefficients.

On each suffix M=[[a,R],[C,A]], the Berkowitz Toeplitz first column is
(1,-a,-RC,-RAC,...,-RA**(m-2)C). Its convolution with the complete monic
coefficient map of A, truncated to the current degree, gives the complete map
for M. Induction from the empty suffix proves completeness. Sparse exact
cancellation and stopping at A**j*C=0 do not discard any later term. The
convolution checks i<=degree and j<=degree-i before adding indices, including
usize::MAX. Independently, 96 complete rational matrices match the existing
native dense characteristic polynomial in Rust, and 72 seeded ZZ/QQ matrices
match SymPy1.14, dense exact characteristic polynomials and Cayley-Hamilton.

## Cold public route and callbacks

The selector only inspects authenticated stock protocol metadata. In particular
it does not populate Basic's ABC cache, call a domain conversion or make a
classification/import callback to prove stock behavior. A current exact native
DomainElement negative Basic cache is required, since the old K.sum conversion
would otherwise execute observable metaclass/import work. Callback code,
delegates, closures, imports, aliases and domain/element/SDM protocols are
checked conservatively. Unknown instance state, subclasses or nonstock entries
keep the old computation. The raw binding validates types and coordinates;
errors and interrupts are distinct from declined optimization.

The first public call can start with this cache cold. It then computes its real
old prefix until an original conversion establishes the cache. A recursive
suffix can enter the native kernel, and the separate native convolution can
finish the parent from the ACTUAL legacy q and Tvals. It never recomputes a
matrix snapshot, repeats a callback or replaces a callback-produced prefix.
A temporary sum hook that restores itself and changes either the source matrix
or the actual prefix result has the same complete value and callback trace as
the forced old path.

A fresh process with a normally constructed QQ SDM sparse seven-cycle reports:
`initial=false`, then first SDM charpoly routes `matrix=true, combine=true`,
then first DomainMatrix charpoly `matrix=true` (first-public-route.json).
No warm-up or ABC invalidation was used there. Dedicated dense/cycle/diagonal
cold tests also invalidate only an unrelated registry, preserving Basic rules.
For diagonal input the prefix may descend to the original conversion and native
convolution resumes the return path. A deeply uncoupled cold input can still
reach historical Python recursion limits before its first conversion. This
residual is explicitly NOT a claim of a fully native cold path or a cold speedup.
Warm stock calls bypass all Python dot/matvec/domain-sum/convert arithmetic.

## Tests and preserved behavior

185 focused Python tests pass locally and on Demeter. Coverage includes exact
ZZ/QQ full maps, large numerators/coprime denominators, source alias/read-only
behavior, real SDM/DomainMatrix APIs without dense conversion, NoRef import
blocking, SIGINT, user exceptions/KeyboardInterrupt and full callback traces.
The old authority is forced by disabling the selector, not just the binding.
Leading isolated vertices, noncontiguous labels, stored zeros and nilpotence
are checked on both the full public native route and public prefix route.
Negative/out-of-range rows and columns, including explicit zero entries, are
refused by the binding and preserve the old full result or explicit error on
both routes. Prefix coefficient indices reject bool/nonint/negative/outside
indices. An empty q declines to retain min(q)'s original error; a nonempty
all-zero q correctly convolves to an empty map.

Four Rust unit tests pass (including 96 independent matrix comparisons,
usize::MAX dimension and convolution overflow, cancellation and interruption).
The broad 55-file matrix/domain gate and Demeter full Rust/clippy are in
progress at source freeze. Final results will be archived separately without
altering these five frozen sources or their manifest. No performance claim.
