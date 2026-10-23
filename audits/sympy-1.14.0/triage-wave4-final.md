# SymPy 1.14 documentation-example triage, final wave-4 round

Replay of all 1167 objects of `triage-wave4-after.json` (its 580 still-differing records and its 587 resolved objects) with the final wave-4 build after waves 4b/4c (integration commit 29524a8, native fingerprint `2970e9f7a7e23507`), same harness as the after round (`tools/sympy_parity.py examples --timeout 30`, immutable native snapshot, two halves under a 4 GB memory cap). Objects whose examples now match are listed under `resolved` in `triage-wave4-final.json`. Still-differing objects whose differing outputs (native and reference) are byte-identical to the after round keep their verdict (`review: carried-*`). Changed objects were re-classified heuristically, and every changed object that was or became a g/e/d/h candidate, or was uncertain (39 objects), was reviewed manually by running it in both libraries (`review: manual`).

## Counts per category

| category | meaning | after round | final |
|---|---|---:|---:|
| g | wrong result | 7 | 0 |
| e | insufficient precision | 1 | 0 |
| d | missing solution / incomplete result | 0 | 0 |
| c | missing domain/feature | 155 | 30 |
| f | different error type/message | 41 | 15 |
| h | reference wrong, symox right (correction candidate) | 19 | 21 |
| b | equivalent result in another form | 128 | 68 |
| i | nondeterministic/environment | 30 | 33 |
| a | printing/form only | 199 | 89 |
| total | | 580 | 256 |

Resolved in total: 911 objects (587 already resolved in the after round, all still matching or without runnable example; 324 newly resolved, by after-round category: g: 5, c: 110, f: 23, b: 61, a: 125). Regressions (resolved before, differing now): 0.

Category changes among still-differing objects: a->b: 1, b->a: 7, c->a: 6, c->b: 6, c->f: 4, c->h: 2, e->c: 1, f->a: 2, f->b: 1, f->c: 1, f->i: 3, g->a: 1, g->c: 1.

## Re-check of the after-round g, e and h objects

| object | after round | final | note |
|---|---|---|---|
| `sympy.DiagonalOf` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.Integral.doit` | g | resolved (examples match) | fixed: 1/3 with no free symbols, as SymPy |
| `sympy.Range` | g | still differs, now a | fixed: inf/sup of Range(ii, jj, kp) raise ValueError like SymPy, Range(0, 10, 3) gives 0 and 9; only pprint of a symbolic Range differs (a) |
| `sympy.combinatorics.Collector.induced_pcgs` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.hn1` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.hn2` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.lambdify` | e | still differs, now c | fixed: lambdify(x, x + 1)(2**60) == 1152921504606846977; remaining differences are c/b/h (see the record) |
| `sympy.physics.continuum_mechanics.Cable.draw` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.physics.continuum_mechanics.Cable.plot_tension` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.physics.continuum_mechanics.Cable.solve` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.physics.optics.gaussian_conj` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.physics.quantum.TensorProduct` | g | resolved (examples match) | fixed: Dagger(A) and Dagger(B)*Dagger(A) |
| `sympy.solve` | g | still differs, now c | fixed: solve(eqs, eqs.atoms(Indexed)) gives {A[1]: 1, A[2]: 2}; remaining differences are c (implicit=True, radical equations, multivariate reduce_inequalities) |
| `sympy.stats.BetaNoncentral` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.BoundedPareto` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.CrossCovarianceMatrix` | g | resolved (examples match) | fixed: (A*X).args == (A, X), shape (k, 1); CCM(A*X, Y).shape == (k, k) |
| `sympy.stats.Dirichlet` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.ExpectationMatrix` | g | resolved (examples match) | fixed: A*X keeps its order and shape; expand gives A*ExpectationMatrix(X) |
| `sympy.stats.GaussianUnitaryEnsemble` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.Levy` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.LogitNormal` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.PoissonProcess` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.PowerFunction` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.Reciprocal` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.RobustSoliton` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |
| `sympy.stats.VarianceMatrix` | g | resolved (examples match) | fixed: A*X keeps its order and shape; expand gives A*VarianceMatrix(X)*A.T |
| `sympy.vector.ImplicitRegion.regular_point` | h | still differs, now h | still correct in symox: repro re-run, symox keeps the correct value or behaviour, SymPy keeps the wrong one |

## Remaining wrong results (g): 0

None.


## Remaining precision issues (e): 0

None.


## Remaining incomplete results (d): 0

None.


## Reference wrong, symox right (h): 21

- `sympy.DiagonalOf`: x is n x 3 with n only integer: true diagonal length is Min(3, n); SymPy assumes the symbolic side is larger (3) and returns None for n x m; native Min forms are exact
- `sympy.combinatorics.Collector.induced_pcgs`: Sylow 3-subgroup of S8: same gens (0 1 2),(3 4 5) in both, generating order 9; induced pcgs must have relative orders multiplying to 9: symox [3,3] right, sympy [3]
- `sympy.hn1`: SymPy rewrite(hankel1) uses hankel1(nu, z), wrong: h1_n = sqrt(pi/(2z)) H1_{n+1/2}; at nu=1,z=2 SymPy gives 0.5111-0.0949i vs hn1=0.4354-0.3506i; symox matches. jn/yn rewrites equivalent (b).
- `sympy.hn2`: SymPy rewrite(hankel2) gives sqrt(pi/(2z))*hankel2(nu,z), wrong (correct order nu+1/2, mpmath-checked); symox right. SymPy rewrite(jn/yn) valid only for integer nu; symox jn - I*yn valid.
- `sympy.physics.continuum_mechanics.Cable.draw`: c.solve() for overdetermined fixed-geometry cable: no tension set satisfies node equilibrium (H=4.05 at Z vs 6.17 at X); SymPy returns inconsistent tensions. Plot diffs are i.
- `sympy.physics.continuum_mechanics.Cable.plot_tension`: First solve(): loads at fixed points Z(2,7.26),X(4,6) give inconsistent H (node equilibrium needs H=4.05 at Z, 6.17 at X); native rightly reports no solution. Plot class diff is env.
- `sympy.physics.continuum_mechanics.Cable.solve`: Ex1 overdetermined: SymPy tensions give H=5.255 on AZ vs 4.054 on ZX (no equilibrium); native says no solution. Ex2: SymPy R_B has |R_B|=T_A=61717, native R_B matches tension_at(100)=50687.
- `sympy.physics.optics.gaussian_conj`: Unchanged verdict, new printed form: symox z_r_out = f**2*z_r/(z_r**2 + (f - s)**2) is the (1 - s/f)**2 thin-lens formula; SymPy uses 1 - (s/f)**2
- `sympy.stats.BetaNoncentral`: Unchanged: cdf(BetaNoncentral(1,1,1))(2) is 1 (SymPy 2*exp(1/2) > 1); Sum index j vs k in pprint (a)
- `sympy.stats.BoundedPareto`: Same formula on [L,H] (checked algebraically); SymPy cdf lacks the x>=H branch so exceeds 1 beyond H; symox gives 1 there
- `sympy.stats.Dirichlet`: Marginal of 2-component Dirichlet is Beta(a1,a2): symox x**(a1-1)*(1-x)**(a2-1)/beta(a1,a2) is correct; SymPy integrates x2 over [0,1] ignoring simplex, result does not integrate to 1
- `sympy.stats.GaussianUnitaryEnsemble`: n=2 weight exp(-Tr H^2) over (h11,h22,Re h12,Im h12) integrates to pi**2/2, so 2/pi**2 (symox) normalizes; SymPy 1/(2*pi**2) uses Z for exp(-Tr/2), integrates to 1/4
- `sympy.stats.Levy`: SymPy cdf erfc(sqrt(c)*sqrt(1/(2z-2mu))) is complex for z<mu; support is [mu,oo), symox adds 0 for z<mu
- `sympy.stats.LogitNormal`: SymPy cdf uses log(z/(1-z)) for all z (complex outside (0,1)); symox gives 0 for z<0 and 1 for z>=1, correct
- `sympy.stats.PoissonProcess`: symox correct: P(X(2)<=3|X(1)>1)=45e^-3/(2(e^3-4))~0.0696 (sum over X(1)=2,3 of independent increments); sympy 5e^-3; P(X!=4) for Poi(6)=1-54e^-6, sympy 1-36e^-6
- `sympy.stats.PowerFunction`: Same (z-a)**2/(b-a)**2 on [a,b]; SymPy cdf has no z>=b branch (exceeds 1), symox returns 1 there
- `sympy.stats.Reciprocal`: Same log(x/a)/log(b/a) on [a,b]; SymPy cdf has no x>=b branch (exceeds 1 beyond b), symox gives 1
- `sympy.stats.RobustSoliton`: sympy Z sums tau over i up to k/R > k, so its pmf over support 1..k sums to 0.9941; symox normalizes over 1..k (sums to 1). E/P differences follow
- `sympy.stats.WienerProcess` (new, previously c): P((X(t1) > 2) | (X(t1) < 4)) is the probability of a sure event: symox 1 is right, SymPy -erf(1)/2 + erf(2)/2 + 1 > 1 is wrong
- `sympy.utilities.autowrap.ufuncify` (new, previously c): The helper printing gap is fixed; with helpers, symox gives the correct array([4., 36.]) (x**2 + y*x**3) where SymPy raises TypeError: return arrays must be of ArrayType; the f2py/Cython differences are environment (i)
- `sympy.vector.ImplicitRegion.regular_point`: Unchanged: no regular point of (x+y+z)**4 exists, SymPy returns a singular point; circle point (-6, 3) lies on (x+2)**2 + (y-3)**2 = 16 (b)

## Missing features (c): 30, by module

Documented: the gap is described as unsupported or different in `docs/algorithmic-parity.md`.

| module | count | objects (documented: yes/no) |
|---|---:|---|
| sympy.tensor | 6 | `shape` (yes), `ArrayKind` (yes), `sort_args_by_name` (yes), `PermuteDims` (yes), `convert_array_to_matrix` (yes), `convert_indexed_to_array` (yes) |
| sympy.functions | 3 | `assoc_laguerre` (yes), `harmonic` (no), `loggamma` (no) |
| sympy.utilities | 3 | `has_dups` (yes), `lambdify` (no), `autowrap` (no) |
| sympy.categories | 2 | `Category` (yes), `subdiagram_from_objects` (yes) |
| sympy.matrices | 2 | `MatrixKind` (yes), `combine_kronecker` (yes) |
| sympy.physics.quantum | 2 | `couple` (yes), `uncouple` (yes) |
| sympy.printing | 2 | `julia_code` (yes), `octave_code` (yes) |
| sympy.stats | 2 | `DiscreteMarkovChain` (yes), `MultivariateNormal` (yes) |
| sympy.concrete | 1 | `euler_maclaurin` (yes) |
| sympy.diffgeom | 1 | `CoordinateSymbol` (yes) |
| sympy.physics.hep | 1 | `simplify_gpgp` (yes) |
| sympy.polys.domains | 1 | `AlgebraicField` (yes) |
| sympy.polys.numberfields | 1 | `prime_decomp` (yes) |
| sympy.series | 1 | `approximants` (yes) |
| sympy.solvers | 1 | `solve` (no) |
| sympy.vector | 1 | `rational_parametrization` (yes) |

## Undocumented c objects: 5

These c objects have at least one missing feature that `docs/algorithmic-parity.md` does not list as unsupported. All of them fail loudly (NotImplementedError, PoleError, ValueError) or leave the expression unevaluated; none returns a wrong value.

- sympy.functions:
  - `sympy.harmonic`: harmonic at rational arguments (harmonic(1/3), harmonic(25/7) with doit/expand_func) is not evaluated in closed form; the limit and expand_func(harmonic(n - k)) gaps are documented
  - `sympy.loggamma`: series(loggamma(x), x, 0, 4) raises PoleError (SymPy gives -log(x) - EulerGamma*x + pi**2*x**2/12 + ...), and expand_func(loggamma(16/3)) is not expanded
- sympy.solvers:
  - `sympy.solve`: solve(..., implicit=True), the radical equations x**(1/3) - x**(1/5) + 1/7 and (x**3 - 3*x**2)**(1/3) + 1 - x (NotImplementedError), and reduce_inequalities over several symbols ([0 <= x - 1, y < 3], [x, y])
- sympy.utilities:
  - `sympy.lambdify`: symox.testing.pytest.ignore_warnings is missing and CompiledFunction has no __globals__ (the missing inspect.getsource is documented)
  - `sympy.utilities.autowrap.autowrap`: a single helper tuple helpers=("helper_func", 4*x, [x]) raises ValueError; SymPy accepts it (a list of tuples works)

## Other categories

- f (different error type/message): 15: `sympy.Abs.nseries`, `sympy.Derivative`, `sympy.Ge`, `sympy.Gt`, `sympy.ImmutableDenseMatrix`, `sympy.ImmutableSparseMatrix`, `sympy.Le`, `sympy.Lt`, `sympy.Predicate`, `sympy.Rational`, `sympy.banded`, `sympy.nsolve`, `sympy.printing.aesaracode.aesara_function`, `sympy.sympify`, `sympy.testing.pytest.warns_deprecated_sympy`.
- i (nondeterministic/environment): 33; b (equivalent form): 68; a (printing/form only): 89.
