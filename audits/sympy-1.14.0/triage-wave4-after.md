# SymPy 1.14 documentation-example triage after wave 4

Replay of every object listed in `triage-wave4.json` (1167 objects, all categories) with the wave-4 build (native fingerprint `622001bd260ea1c8`), same harness (`tools/sympy_parity.py examples --timeout 30`, immutable native snapshot). Objects whose examples now match are listed under `resolved` in `triage-wave4-after.json`. Still-differing objects whose differing outputs are byte-identical to the previous round keep their previous verdict (`review: carried-*`); changed ones were re-classified heuristically, and every changed object that was or became a g/e/d/h candidate, or was uncertain, was reviewed manually (`review: manual`).

## Counts per category

| category | meaning | before | after |
|---|---|---:|---:|
| g | wrong result | 71 | 7 |
| e | insufficient precision | 4 | 1 |
| d | missing solution / incomplete result | 3 | 0 |
| c | missing domain/feature | 557 | 155 |
| f | different error type/message | 45 | 41 |
| h | reference wrong, symox right (correction candidate) | 15 | 19 |
| b | equivalent result in another form | 137 | 128 |
| i | nondeterministic/environment | 23 | 30 |
| a | printing/form only | 312 | 199 |
| total | | 1167 | 580 |

Resolved (all examples now match, or no example is runnable): 587 objects; by previous category: g: 53, e: 3, d: 1, c: 354, f: 5, b: 21, i: 9, a: 141.

Category changes among still-differing objects: a->b: 1, c->a: 26, c->b: 8, c->f: 1, c->g: 4, c->h: 3, c->i: 16, d->a: 1, d->b: 1, f->a: 1, g->a: 1, g->b: 2, g->c: 10, g->f: 1, g->h: 1.

## Remaining wrong results (g): 7

- `sympy.Integral.doit` (sympy.integrals, previously g): Native result keeps the integration variable t free in the condition x/t < 0 (free_symbols {t}) outside the t-integral; true value 1/3 is constant.
  - repro: `from symox import *; from symox.abc import x,t; r=(x**2+Piecewise((0,x/t<0),(1,True))).integrate((t,S(4)/5,1),(x,-1,1)); print(r, r.free_symbols) # sympy: 1/3 set(), symox: Integral(Piecewise(...x/t < 0...)) {t}`
- `sympy.Range` (sympy.sets, previously c): Range(i,j,k).inf returns j - k where sympy raises; for positive step inf/sup are swapped: Range(ii,jj,kp).inf=jj-kp, .sup=ii; subs to (0,10,3) gives 7 (true inf 0)
  - repro: `from symox import *; ii,jj=symbols('ii jj',integer=True); kp=Symbol('kp',positive=True,integer=True); print(Range(ii,jj,kp).inf.subs({ii:0,jj:10,kp:3}), Range(ii,jj,kp).sup) # sympy: ValueError (cannot tell if Range is null), symox: 7 ii (inf of Range(0,10,3) is 0)`
- `sympy.physics.quantum.TensorProduct` (sympy.physics.quantum, previously g): Dagger of noncommutative symbol returns conjugate(A), not Dagger(A); Dagger(A*B) gives conjugate(A)*conjugate(B) (order not reversed). expand(tensorproduct=True) also not distributed.
  - repro: `from symox import *; from symox.physics.quantum import Dagger; A,B=symbols('A B',commutative=False); print(Dagger(A), Dagger(A*B)) # sympy: Dagger(A) Dagger(B)*Dagger(A), symox: conjugate(A) conjugate(A)*conjugate(B)`
- `sympy.solve` (sympy.solvers, previously g): solve(eqs, eqs.atoms(Indexed)) returns [] (claims no solution) because atoms(Indexed) is empty (Indexed class mismatch); SymPy gives {A[1]:1, A[2]:2}. Also c: implicit, radical eqs, multivar ineq.
  - repro: `from symox import *; A=IndexedBase('A'); eqs=Tuple(A[1]+A[2]-3, A[1]-A[2]+1); print(solve(eqs, eqs.atoms(Indexed)), eqs.atoms(Indexed)) # sympy: {A[1]: 1, A[2]: 2} {A[1], A[2]}, symox: [] set()`
- `sympy.stats.CrossCovarianceMatrix` (sympy.stats, previously c): A*X with RandomMatrixSymbol X is built as X*A (shape (k,k) instead of (k,1)); CCM(A*X,Y).shape (k,1) vs (k,k). Also expand missing/raises
  - repro: `from symox import *; from symox.stats import *; from symox.stats.rv import RandomMatrixSymbol; k=Symbol('k'); A=MatrixSymbol('A',k,k); X=RandomMatrixSymbol('X',k,1); print((A*X).args, (A*X).shape) # sympy: (A, X) (k, 1), symox: (X, A) (k, k); print(CrossCovarianceMatrix(A*X, RandomMatrixSymbol('Y',k,1)).shape) # sympy (k, k), symox (k, 1)`
- `sympy.stats.ExpectationMatrix` (sympy.stats, previously c): A*X with RandomMatrixSymbol X is stored as X*A (wrong order/shape (k,k)); also no .shape, expand not linear/raises
  - repro: `from symox import *; from symox.stats import *; from symox.stats.rv import RandomMatrixSymbol; k=Symbol('k'); A=MatrixSymbol('A',k,k); X=RandomMatrixSymbol('X',k,1); print((A*X).args, (A*X).shape) # sympy: (A, X) (k, 1), symox: (X, A) (k, k)`
- `sympy.stats.VarianceMatrix` (sympy.stats, previously c): A*X with RandomMatrixSymbol X is stored as X*A (wrong order, shape (k,k)); expand also not done
  - repro: `from symox import *; from symox.stats import *; from symox.stats.rv import RandomMatrixSymbol; k=Symbol('k'); A=MatrixSymbol('A',k,k); X=RandomMatrixSymbol('X',k,1); print((A*X).args, (A*X).shape) # sympy: (A, X) (k, 1), symox: (X, A) (k, k)`

## Remaining precision issues (e): 1

- `sympy.lambdify` (sympy.utilities, previously e): default (numpy) lambdify converts int args to float: lambdify(x, x + 1)(2**60) gives 1.152921504606847e+18 vs exact 1152921504606846977; also Eq gives 1.0 not bool, modules=sympy Matrix -> list
  - repro: `from symox import *; from symox.abc import x; print(lambdify(x, x + 1)(2**60)) # sympy: 1152921504606846977, symox: 1.152921504606847e+18`

## Remaining incomplete results (d): 0

None.


## Reference wrong, symox right (h): 19

- `sympy.DiagonalOf`: x is n x 3 with n only integer: true diagonal length is Min(3, n); SymPy assumes the symbolic side is larger (3) and returns None for n x m; native Min forms are exact
- `sympy.combinatorics.Collector.induced_pcgs`: Sylow 3-subgroup of S8: same gens (0 1 2),(3 4 5) in both, generating order 9; induced pcgs must have relative orders multiplying to 9: symox [3,3] right, sympy [3]
- `sympy.hn1`: SymPy rewrite(hankel1) uses hankel1(nu, z), wrong: h1_n = sqrt(pi/(2z)) H1_{n+1/2}; at nu=1,z=2 SymPy gives 0.5111-0.0949i vs hn1=0.4354-0.3506i; symox matches. jn/yn rewrites equivalent (b).
- `sympy.hn2`: SymPy rewrite(hankel2) gives sqrt(pi/(2z))*hankel2(nu,z), wrong (correct order nu+1/2, mpmath-checked); symox right. SymPy rewrite(jn/yn) valid only for integer nu; symox jn - I*yn valid.
- `sympy.physics.continuum_mechanics.Cable.draw`: c.solve() for overdetermined fixed-geometry cable: no tension set satisfies node equilibrium (H=4.05 at Z vs 6.17 at X); SymPy returns inconsistent tensions. Plot diffs are i.
- `sympy.physics.continuum_mechanics.Cable.plot_tension`: First solve(): loads at fixed points Z(2,7.26),X(4,6) give inconsistent H (node equilibrium needs H=4.05 at Z, 6.17 at X); native rightly reports no solution. Plot class diff is env.
- `sympy.physics.continuum_mechanics.Cable.solve`: Ex1 overdetermined: SymPy tensions give H=5.255 on AZ vs 4.054 on ZX (no equilibrium); native says no solution. Ex2: SymPy R_B has |R_B|=T_A=61717, native R_B matches tension_at(100)=50687.
- `sympy.physics.optics.gaussian_conj`: SymPy uses (1-(s/f)**2) instead of (1-s/f)**2 in m and z_r_out (Self 1983: m=1/sqrt((1-s/f)^2+(z_r/f)^2)); symox formula correct; s_out equal (s=3,z=2,f=5: 2.5 vs 6.25)
- `sympy.stats.BetaNoncentral`: cdf(BNC(1,1,1))(2): SymPy 2*exp(1/2)>1 integrates density beyond support [0,1]; symox 1 is the correct CDF. pprint of Sum falls back to str (a)
- `sympy.stats.BoundedPareto`: Same formula on [L,H] (checked algebraically); SymPy cdf lacks the x>=H branch so exceeds 1 beyond H; symox gives 1 there
- `sympy.stats.Dirichlet`: Marginal of 2-component Dirichlet is Beta(a1,a2): symox x**(a1-1)*(1-x)**(a2-1)/beta(a1,a2) is correct; SymPy integrates x2 over [0,1] ignoring simplex, result does not integrate to 1
- `sympy.stats.GaussianUnitaryEnsemble`: n=2 weight exp(-Tr H^2) over (h11,h22,Re h12,Im h12) integrates to pi**2/2, so 2/pi**2 (symox) normalizes; SymPy 1/(2*pi**2) uses Z for exp(-Tr/2), integrates to 1/4
- `sympy.stats.Levy`: SymPy cdf erfc(sqrt(c)*sqrt(1/(2z-2mu))) is complex for z<mu; support is [mu,oo), symox adds 0 for z<mu
- `sympy.stats.LogitNormal`: SymPy cdf uses log(z/(1-z)) for all z (complex outside (0,1)); symox gives 0 for z<0 and 1 for z>=1, correct
- `sympy.stats.PoissonProcess`: symox correct: P(X(2)<=3|X(1)>1)=45e^-3/(2(e^3-4))~0.0696 (sum over X(1)=2,3 of independent increments); sympy 5e^-3; P(X!=4) for Poi(6)=1-54e^-6, sympy 1-36e^-6
- `sympy.stats.PowerFunction`: Same (z-a)**2/(b-a)**2 on [a,b]; SymPy cdf has no z>=b branch (exceeds 1), symox returns 1 there
- `sympy.stats.Reciprocal`: Same log(x/a)/log(b/a) on [a,b]; SymPy cdf has no x>=b branch (exceeds 1 beyond b), symox gives 1
- `sympy.stats.RobustSoliton`: sympy Z sums tau over i up to k/R > k, so its pmf over support 1..k sums to 0.9941; symox normalizes over 1..k (sums to 1). E/P differences follow
- `sympy.vector.ImplicitRegion.regular_point`: (x+y+z)**4: every point of the surface has zero gradient, so no regular point exists; SymPy returns singular (-10,-10,20), symox raises. Circle point (0,3-2sqrt3) is valid

## Missing features (c): 155, by module

| module | count | objects |
|---|---:|---|
| sympy.stats | 20 | `BernoulliProcess`, `Binomial`, `ContinuousMarkovChain`, `Die`, `DiscreteMarkovChain`, `canonical_form`, `communication_classes`, `decompose`, `stationary_distribution`, `DiscreteRV`, `Expectation`, `GammaProcess` ... |
| sympy.functions | 18 | `FallingFactorial`, `airyai`, `airyaiprime`, `airybi`, `airybiprime`, `assoc_laguerre`, `elliptic_e`, `elliptic_f`, `elliptic_k`, `elliptic_pi`, `expint`, `harmonic` ... |
| sympy.matrices | 11 | `LUdecomposition`, `schur`, `diag`, `has`, `has_matching_shape`, `structurally_equal`, `copyin_matrix`, `MutableSparseMatrix`, `MatrixKind`, `combine_kronecker`, `SparseRepMatrix` |
| sympy.printing | 10 | `ccode`, `fcode`, `glsl_code`, `jscode`, `julia_code`, `latex`, `octave_code`, `rcode`, `rust_code`, `smtlib_code` |
| sympy.physics.quantum | 8 | `Wavefunction`, `integrate_result`, `IntQubit`, `represent`, `JzKet`, `JzKetCoupled`, `couple`, `uncouple` |
| sympy.physics.control | 7 | `MIMOFeedback`, `doit`, `MIMOParallel`, `MIMOSeries`, `TransferFunctionMatrix`, `from_Matrix`, `phase_margin` |
| sympy.simplify | 6 | `nsimplify`, `opt_cse`, `TR1`, `devise_plan_meijer`, `sqrt_biquadratic_denest`, `trigsimp_old` |
| sympy.tensor | 6 | `shape`, `ArrayKind`, `sort_args_by_name`, `PermuteDims`, `convert_array_to_matrix`, `convert_indexed_to_array` |
| sympy.vector | 6 | `orient_new`, `orient_new_quaternion`, `gradient`, `multiplicity`, `rational_parametrization`, `VectorKind` |
| sympy.categories | 5 | `Category`, `Diagram`, `hom`, `subdiagram_from_objects`, `DiagramGrid` |
| sympy.core | 5 | `nseries`, `Derivative`, `Rational`, `expand`, `symbols` |
| sympy.solvers | 5 | `diophantine`, `factor_system`, `nonlinsolve`, `solve_triangulated`, `classify_sysode` |
| sympy.utilities | 5 | `has_dups`, `public`, `ufuncify`, `get_datatype`, `codegen` |
| sympy.codegen | 4 | `newtons_method`, `SeriesApprox`, `MatrixSolve`, `optimize` |
| sympy.concrete | 4 | `Product`, `Sum`, `euler_maclaurin`, `summation` |
| sympy.holonomic | 4 | `HolonomicFunction`, `to_expr`, `to_meijerg`, `to_sequence` |
| sympy.integrals | 3 | `as_sum`, `principal_value`, `transform` |
| sympy.polys.partfrac | 3 | `apart`, `apart_list`, `assemble_partfrac_list` |
| sympy.polys.polytools | 3 | `Poly`, `per`, `which_all_roots` |
| sympy.series | 3 | `O`, `approximants`, `fps` |
| sympy.sets | 3 | `powerset`, `ProductSet`, `SetKind` |
| sympy.assumptions | 2 | `Predicate`, `refine_matrixelement` |
| sympy.calculus | 2 | `apply_finite_diff`, `periodicity` |
| sympy.physics.optics | 2 | `RayTransferMatrix`, `jones_2_stokes` |
| sympy.diffgeom | 1 | `CoordinateSymbol` |
| sympy.geometry | 1 | `Line` |
| sympy.logic | 1 | `entails` |
| sympy.ntheory | 1 | `jacobi_symbol` |
| sympy.physics.hep | 1 | `simplify_gpgp` |
| sympy.plotting | 1 | `plot_implicit` |
| sympy.polys.domains | 1 | `AlgebraicField` |
| sympy.polys.fields | 1 | `sfield` |
| sympy.polys.numberfields | 1 | `prime_decomp` |
| sympy.polys.rootoftools | 1 | `CRootOf` |
