# SymPy 1.14 documentation-example triage, wave-8 replay

Replay of every documented object of the catalog (2916 objects, ordinary examples and keyword examples) with `tools/sympy_parity.py examples [--keywords] --timeout 30`, an immutable native snapshot and a memory cap, in shards. A first full replay ran with the native build of integration commit a2a0038 (waves 5-7); it found the regressions below. A power cut then destroyed its raw outputs, so both replays were run again with the build carrying the regression fixes (native fingerprint `27ba41943237554c`). The objects whose examples call `expand` were re-run with the final build (expand fix). Still-differing objects of `triage-wave4-final.json` keep their verdict (`review: carried-*`; their differing outputs were byte-identical in the a2a0038 replay); the 9 objects whose differences changed were reviewed manually (`review: manual`). Intentional corrections (h) are kept apart from matches.

## Counts per category

| category | meaning | wave-4 final | wave 8 |
|---|---|---:|---:|
| g | wrong result | 0 | 0 |
| e | insufficient precision | 0 | 1 |
| d | missing solution / incomplete result | 0 | 0 |
| c | missing domain/feature | 30 | 13 |
| f | different error type/message | 15 | 5 |
| h | reference wrong, symox right (correction candidate) | 21 | 20 |
| b | equivalent result in another form | 68 | 66 |
| i | nondeterministic/environment | 33 | 34 |
| a | printing/form only | 89 | 89 |
| total | | 256 | 228 |

Replay statuses (ordinary examples): examples_match: 2683, incomplete: 2, no_successful_examples: 2, observation_difference: 229.
Comparison with the baseline: fixed: 28, match: 2659, reference_timeout: 1, still_differs: 228.

`reference_timeout`: SymPy itself exceeded the 30 s budget on the loaded machine while symox completed every example (not a symox difference):

- `sympy.physics.mechanics.ReferenceFrame.orient_space_fixed`: SymPy timed out under load in keyword mode; symox completed every example, and the ordinary examples match
- `sympy.solvers.ode.ode.ode_sol_simplicity`: SymPy 1.14 needs 53 s for these examples on the loaded machine (load ~20); a direct SymPy run gives -2, -1, oo, [28, 35], Eq(f(x)/tan(f(x)/(2*x)), C1), identical to the symox outputs

## Fixed since the last triage: 28 objects

By wave-4 final category: c: 16, f: 8, b: 3, a: 1.

`sympy.Abs.nseries`, `sympy.Ge`, `sympy.Gt`, `sympy.Le`, `sympy.Lt`, `sympy.Rational`, `sympy.banded`, `sympy.categories.Category`, `sympy.categories.Diagram.subdiagram_from_objects`, `sympy.diffgeom.CoordinateSymbol`, `sympy.has_dups`, `sympy.loggamma`, `sympy.matrices.MatrixKind`, `sympy.matrices.matrixbase.classof`, `sympy.nsimplify`, `sympy.physics.control.PIDController.eval_frequency`, `sympy.physics.control.phase_margin`, `sympy.physics.hep.gamma_matrices.simplify_gpgp`, `sympy.physics.quantum.spin.couple`, `sympy.physics.quantum.spin.uncouple`, `sympy.prime_decomp`, `sympy.sympify`, `sympy.tensor.array.ArrayKind`, `sympy.tensor.array.expressions.ArrayContraction.sort_args_by_name`, `sympy.tensor.array.expressions.PermuteDims`, `sympy.tensor.array.expressions.convert_array_to_matrix`, `sympy.tensor.array.expressions.convert_indexed_to_array`, `sympy.vector.ImplicitRegion.rational_parametrization`.

## Regressions found with the a2a0038 build: 3 (ordinary examples)

| object | baseline | a2a0038 replay | cause | fix |
|---|---|---|---|---|
| `sympy.FunctionMatrix` | triage-wave4-final (resolved) | incomplete (native timeout) | (Y**2)[10, 10] of a 1000x1000 FunctionMatrix materialised the whole explicit power (O(d**3)): 52 s, 21 s in wave 4 | MatrixExpr::element of a square of a lazy finite matrix multiplies one row by one column (SymPy MatMul._entry), then dotprodsimp as the explicit power; 0.5 s (`crates/symox-core/src/matrixexpr.rs`) |
| `sympy.assumptions.BinaryRelation` | wave4-examples (match) | observation_difference (AttributeError on str of Q.R(2, 1)) | Predicate.__str__ uses _symox_predicate, which a BinaryRelation subclass does not set | BinaryRelation.__str__ = __repr__ (`python/symox/assumptions/relation/binrel.py`) |
| `sympy.physics.wigner.gaunt` | wave4-examples (match) | incomplete (native timeout) | gaunt(1000, 1000, 1200, 9, 3, -12): 93 s (SymPy 2 s, wave 4 6.5 s). Integer._eval_power (int_rational_pow, waves 5-7) calls perfect_power on huge factorial products, which tried an integer root for every prime exponent up to bits(n) (157 s for one sqrt) | perfect_power follows SymPy's perfect_power(factor=True): trial division by small primes, a multiplicity 1 rules out every power, otherwise only divisors of the multiplicity are tried (`crates/symox-core/src/ntheory/factor.rs`) |

All three match SymPy again in the wave-8 replay.

## Keyword examples

Against `fifth-keyword-examples.json.gz` overlaid by the later keyword batches: fixed: 1116, match: 1488, not_run: 16, reference_timeout: 2, regression: 1, still_differs: 309 (`not_run`: alias names absent from the current catalog). Keyword regressions:

- `sympy.assumptions.BinaryRelation`: fixed. same cause as the example regression
- `sympy.physics.wigner.gaunt`: fixed. same cause as the example regression
- `sympy.physics.control.PIDController.eval_frequency`: fixed. expand() rebuilt the held 1/(0.99 + 0.2*I) of TransferFunction.to_expr and evaluated it to 0.980296049406921*(0.99 - 0.2*I); SymPy's multinomial hint leaves base**-1 alone and no other hint applies, so the node is now returned unchanged (crates/symox-core/src/expand.rs)
- `sympy.stats.Zeta`: b (documented). variance(X) is returned over a common denominator, (-pi**8/8100 + zeta(3)*zeta(5))/zeta(5)**2, where SymPy gives -pi**8/(8100*zeta(5)**2) + zeta(3)/zeta(5): the polynomial-moment path of symox.stats simplifies its sum, SymPy's moment-generating-function path does not. Same value; already b in triage-wave4-final for the ordinary examples

Category changes among still-differing objects: c->e: 1, c->i: 1, f->a: 1, f->b: 1, h->c: 1.

## Re-check of the wave-4 final g, e and h objects

| object | wave-4 final | wave 8 | note |
|---|---|---|---|
| `sympy.DiagonalOf` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.combinatorics.Collector.induced_pcgs` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.hn1` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.hn2` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.physics.continuum_mechanics.Cable.draw` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.physics.continuum_mechanics.Cable.plot_tension` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.physics.continuum_mechanics.Cable.solve` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.physics.optics.gaussian_conj` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.BetaNoncentral` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.BoundedPareto` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.Dirichlet` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.GaussianUnitaryEnsemble` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.Levy` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.LogitNormal` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.PoissonProcess` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.PowerFunction` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.Reciprocal` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.RobustSoliton` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.stats.WienerProcess` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.utilities.autowrap.ufuncify` | h | still differs, now h | carried: differing outputs unchanged; symox keeps the correct value or behaviour |
| `sympy.vector.ImplicitRegion.regular_point` | h | still differs, now c | The circle example now gives SymPy's point. Other regions still give another point or NotImplementedError from the regular_point search (documented gap). |

## Wrong results (g): 0

None.


## Precision issues (e): 1

- `sympy.solve` (sympy.solvers, previously c): New (previously masked by a missing solution): the radical equation x**(1/3) - x**(1/5) + 1/7 now solves as SymPy, to CRootOf(7*x**5 - 7*x**3 + 1, k)**15. Evaluating the exactly zero residual at a root reports spurious digits: [abs(eq.subs(x, s).n(2)) for s in sol] gives 6.6e-24 and 0 where SymPy gives its scaled zero 0.e-110 for the two exact roots (the values are correct; only the accuracy claimed for a cancelled zero differs). The other solve differences are a/b (set=True tuple order, list order, -oo < x form).
  - repro: `r = CRootOf(7*x**5 - 7*x**3 + 1, 1); (root(r**15, 3) - root(r**15, 5) + Rational(1, 7)).n(2) # SymPy: 0.e-110`

## Incomplete results (d): 0

None.


## Reference wrong, symox right (h): 20

`sympy.DiagonalOf`, `sympy.combinatorics.Collector.induced_pcgs`, `sympy.hn1`, `sympy.hn2`, `sympy.physics.continuum_mechanics.Cable.draw`, `sympy.physics.continuum_mechanics.Cable.plot_tension`, `sympy.physics.continuum_mechanics.Cable.solve`, `sympy.physics.optics.gaussian_conj`, `sympy.stats.BetaNoncentral`, `sympy.stats.BoundedPareto`, `sympy.stats.Dirichlet`, `sympy.stats.GaussianUnitaryEnsemble`, `sympy.stats.Levy`, `sympy.stats.LogitNormal`, `sympy.stats.PoissonProcess`, `sympy.stats.PowerFunction`, `sympy.stats.Reciprocal`, `sympy.stats.RobustSoliton`, `sympy.stats.WienerProcess`, `sympy.utilities.autowrap.ufuncify`.

## Missing features (c): 13, by module

| module | count | objects |
|---|---:|---|
| sympy.functions | 2 | `assoc_laguerre`, `harmonic` |
| sympy.printing | 2 | `julia_code`, `octave_code` |
| sympy.stats | 2 | `DiscreteMarkovChain`, `MultivariateNormal` |
| sympy.concrete | 1 | `euler_maclaurin` |
| sympy.matrices | 1 | `combine_kronecker` |
| sympy.polys.domains | 1 | `AlgebraicField` |
| sympy.series | 1 | `approximants` |
| sympy.tensor | 1 | `shape` |
| sympy.utilities | 1 | `lambdify` |
| sympy.vector | 1 | `regular_point` |

## Other categories

- f (different error type/message): 5: `sympy.ImmutableDenseMatrix`, `sympy.ImmutableSparseMatrix`, `sympy.nsolve`, `sympy.printing.aesaracode.aesara_function`, `sympy.testing.pytest.warns_deprecated_sympy`.
- i (nondeterministic/environment): 34; b (equivalent form): 66; a (printing/form only): 89.
