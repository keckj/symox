# LRT/vector scalar dispatch correction — 2026-10-27

Private source copied from immutable R3 source identity 1ad26aa2fbe980b159a50f38583dce6fc0a301045b92194803207c059ef07875. Native library remains 8a8a8e126be486e31b83b7cae2c306e54ab21e9d569f01bbe46adeb8bb661f3a; no Rust edit or build, main and earlier deliveries untouched.

## Cause and bounded correction

Collecting the fourteen files of R3 shard 29 installs symox.vector's official arithmetic wrapper outside the already-captured Tensor/public-core delegate (and possibly Ore). The old selector accepted only one direct Ore wrapper and refused this legitimate vector layer. The existing native-route assertion then executed the complete Python LRT recurrence. Its mathematical result was unchanged.

The vector installer now records its four scalar arithmetic wrappers, each exact delegate, function code/closure/default snapshots, module namespace, _BasisDependent class/MRO, and the builtin isinstance used by the wrapper. This records protocol metadata only: the vector dispatch branches and their order are unchanged. A static code template and installer anchor are fixed before installation.

rationaltools walks only registry members whose actual function identity matches the current operator. Every layer must carry its matching certificate, exact code, original closure cells and delegates. Vector layers additionally require the captured plain-type _BasisDependent/MRO and unchanged isinstance lookup, and prove that every admitted native scalar class has no such basis class in its MRO. Thus their isinstance branch is false on all native scalar operands and the wrapper forwards unchanged. Ore layers retain all existing type, namespace, marker, priority, attribute and closure checks. The chain must terminate at the bootstrap-captured operator; arbitrary __wrapped__ links are never followed, and cycles decline. Name is part of the per-call verification cache key. Failed, ambiguous, modified or reloaded certificates retain the full existing algorithm.

Existing scalar priority/MRO/attribute checks, LRT bootstrap/resultant/as_poly/constructor checks, Symbol/Dummy restrictions, native-domain eligibility and host fallback remain in place. No mathematical formula or domain is changed, and no host exceptions or interrupts are caught by the dispatch adaptation.

## Tests and limits

New tests check native selection for both LRT and Atan across vector-only, vector→Ore, Ore→vector and Tensor→vector→Ore imports, three rationaltools reloads, exact reference LRT results (existing cases), real vector multiplication and the established division error. They compare observable KeyboardInterrupt callbacks against a forced full fallback for mutated vector delegates, module-global isinstance and a user @wraps wrapper; a cyclic delegate conservatively declines without recursion.

The existing Atan mutation tests assumed Ore was the outer wrapper and raised StopIteration before their assertions when combined with vector collection. Their helper now locates the same Ore record through recorded vector delegates; every mutation, fallback and exception assertion is unchanged. No failing route assertion is removed or relaxed.

Local: existing shard-route reproduction 1 passed / 943 deselected; LRT+forwarding 187 passed; combined Atan/vector 179 passed. Demeter first pass: complete shard 29 944 passed, LRT+forwarding 187 passed. Final exact-source combined Demeter run: **1153 passed in 56.93 s**, all fourteen shard files followed by LRT forwarding and the five Atan/vector files, in one process. Remote four-file and native binary hashes match source-manifest.json. Result is in combined-final-demeter.log.

This is a bounded protocol-selection correction for the official vector/Ore scalar chains, not a global mathematical parity or speed claim. Unsupported/custom protocols continue through the established full implementation. R3's previously archived full-suite failure remains historical evidence; a complete corrected-snapshot suite is the integration gate owned by root.
