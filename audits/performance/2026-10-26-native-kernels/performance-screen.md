# Public algorithm performance screening

Reference: SymPy 1.14.0. Symox distribution: 0.1.1.
Package fingerprint: `d564f63a251e256d46650aa3d354ecc50e502d4200ad3f69aed00b8e80c5fc52`.

137 workloads; 132 verified and comparable; 5 excluded.
31 comparable workloads have a first-call or repeated-call speedup ≤ 1.5.

Speedup means reference time / native time; values below 1 mean Symox is slower.
Fresh processes alternate backend order on one CPU. Inputs and result verification are outside timing.
Repeated calls retain both libraries' caches. These measurements are not a universal parity or complexity certificate.
Three-process screenings are provisional; bootstrap intervals describe observed variation, not host isolation.

## Weak or marginal cases

| Workload | First speedup | Repeated speedup | Native first, ms | Native repeated, ms | Processes |
|---|---:|---:|---:|---:|---:|
| `lambdify.evaluate.20` | 0.113 | 0.489 | 18.0428 | 3.83603 | 1 |
| `solvers.polynomial.7` | 1.368 | 0.254 | 3.0887 | 6.78371 | 1 |
| `number_theory.perfect_power` | 0.359 | 16.018 | 0.248925 | 0.00253616 | 1 |
| `combinatorics.partition.200` | 5.615 | 0.490 | 0.084624 | 0.000548731 | 1 |
| `number_theory.primepi.100000` | 10.020 | 0.491 | 0.068998 | 0.0005473 | 1 |
| `combinatorics.partition.50` | 3.550 | 0.491 | 0.053428 | 0.000545113 | 1 |
| `number_theory.primepi.1000000` | 12.670 | 0.493 | 0.139828 | 0.000548498 | 1 |
| `numerics.large_argument.50` | 0.537 | 32.313 | 0.214197 | 0.000734509 | 1 |
| `numerics.large_argument.100` | 0.586 | 44.323 | 0.221578 | 0.000712515 | 1 |
| `lambdify.evaluate.100` | 0.671 | 0.671 | 26.3342 | 26.0988 | 1 |
| `matrices.charpoly_integer.6` | 0.691 | 2.879 | 0.820591 | 0.119457 | 1 |
| `numerics.tiny_sine.100` | 0.699 | 58.839 | 0.098044 | 0.000751145 | 1 |
| `substitution.subs.12` | 40.706 | 0.716 | 4.36342 | 0.111858 | 1 |
| `parsing.sympify.50` | 16.885 | 0.722 | 5.71361 | 6.50067 | 1 |
| `numerics.tiny_sine.50` | 0.733 | 60.686 | 0.092648 | 0.000712433 | 1 |
| `construction.incremental_sum.100` | 21.552 | 0.756 | 2.96654 | 0.470985 | 1 |
| `limits.removable` | 65.933 | 0.764 | 2.26383 | 0.241012 | 1 |
| `numerics.constant.300` | 0.775 | 284.445 | 1.26884 | 0.000711255 | 1 |
| `polynomials.resultant.3` | 0.890 | 1.658 | 3.62649 | 0.39805 | 1 |
| `numerics.constant.100` | 0.902 | 162.826 | 0.509337 | 0.000718178 | 1 |
| `combinatorics.symmetric_group.8` | 4.238 | 1.010 | 0.047811 | 0.000115053 | 1 |
| `combinatorics.symmetric_group.12` | 4.276 | 1.011 | 0.047365 | 0.00011419 | 1 |
| `solvers.polynomial.4` | 10.523 | 1.096 | 2.17519 | 2.82684 | 1 |
| `matrices.charpoly_integer.12` | 1.149 | 4.433 | 1.24651 | 0.277759 | 1 |
| `expansion.multinomial.12` | 13.447 | 1.228 | 5.87753 | 0.0011417 | 1 |
| `numerics.constant.50` | 1.259 | 135.004 | 0.289437 | 0.000717976 | 1 |
| `lambdify.compile_and_evaluate.20` | 1.272 | 6.176 | 100.68 | 4.19385 | 1 |
| `assumptions.tiny_positive` | 5.729 | 1.343 | 0.064705 | 9.64656e-05 | 1 |
| `solvers.nsolve` | 54.475 | 1.353 | 0.782398 | 0.611768 | 1 |
| `matrices.charpoly_generic.5` | 1.403 | 1.489 | 5.93024 | 4.06382 | 1 |
| `lambdify.compile_and_evaluate.100` | 1.486 | 3.811 | 152.743 | 27.3724 | 1 |

## Excluded results

These statuses require diagnosis; a validation timeout or alternate formatting is not evidence of an incorrect native algorithm.

| Workload | Status |
|---|---|
| `assumptions.trig_zero` | undecided_assumption |
| `differentiation.composite.8` | validation_timeout |
| `integration.rational.8` | not_comparable |
| `series.transcendental.20` | not_comparable |
| `series.transcendental.40` | not_comparable |

## All comparable cases

| Workload | First speedup | Repeated speedup |
|---|---:|---:|
| `assumptions.tiny_positive` | 5.729 | 1.343 |
| `combinatorics.partition.200` | 5.615 | 0.490 |
| `combinatorics.partition.50` | 3.550 | 0.491 |
| `combinatorics.permutation.100` | 1.556 | 2.160 |
| `combinatorics.permutation.1000` | 1.750 | 1.974 |
| `combinatorics.symmetric_group.12` | 4.276 | 1.011 |
| `combinatorics.symmetric_group.8` | 4.238 | 1.010 |
| `construction.bulk_add.100` | 43.553 | 1.874 |
| `construction.bulk_add.1000` | 203.499 | 5.456 |
| `construction.incremental_sum.100` | 21.552 | 0.756 |
| `construction.incremental_sum.500` | 86.126 | 116.558 |
| `differentiation.composite.1` | 102.158 | 2.452 |
| `differentiation.composite.4` | 24.267 | 7.703 |
| `expansion.multinomial.12` | 13.447 | 1.228 |
| `expansion.multinomial.16` | 55.257 | 6.087 |
| `expansion.multinomial.8` | 48.031 | 1.572 |
| `expansion.univariate.100` | 42.594 | 1.579 |
| `expansion.univariate.30` | 67.098 | 1.553 |
| `integration.dirac_affine` | 87.433 | 65.869 |
| `integration.exponential_trig` | 412.814 | 124.165 |
| `integration.gamma_moment.10` | 326.357 | 4.444 |
| `integration.gamma_moment.50` | 329.325 | 7.271 |
| `integration.polynomial.100` | 40.975 | 1.536 |
| `integration.polynomial.30` | 49.611 | 1.754 |
| `integration.rational.4` | 202.186 | 172.019 |
| `lambdify.compile_and_evaluate.100` | 1.486 | 3.811 |
| `lambdify.compile_and_evaluate.20` | 1.272 | 6.176 |
| `lambdify.evaluate.100` | 0.671 | 0.671 |
| `lambdify.evaluate.20` | 0.113 | 0.489 |
| `limits.at_infinity` | 161.677 | 164.817 |
| `limits.removable` | 65.933 | 0.764 |
| `matrices.charpoly_generic.3` | 1.516 | 1.627 |
| `matrices.charpoly_generic.4` | 1.838 | 1.540 |
| `matrices.charpoly_generic.5` | 1.403 | 1.489 |
| `matrices.charpoly_generic.6` | 1.872 | 1.955 |
| `matrices.charpoly_integer.12` | 1.149 | 4.433 |
| `matrices.charpoly_integer.6` | 0.691 | 2.879 |
| `matrices.charpoly_shared.4` | 2.343 | 2.713 |
| `matrices.det_generic.3` | 249.802 | 6.205 |
| `matrices.det_generic.4` | 641.147 | 625.013 |
| `matrices.det_integer.20` | 366.514 | 668.553 |
| `matrices.det_integer.8` | 125.668 | 432.887 |
| `matrices.eigenvalues` | 74.046 | 8.122 |
| `matrices.inverse_generic.3` | 56.786 | 19.317 |
| `matrices.inverse_rational.10` | 9.481 | 9.194 |
| `matrices.inverse_rational.20` | 4.736 | 4.555 |
| `matrices.inverse_rational.4` | 12.381 | 58.523 |
| `matrices.jacobian.16` | 465.524 | 16.060 |
| `matrices.jacobian.8` | 383.401 | 27.784 |
| `matrices.linear_solve.20` | 14.034 | 9.600 |
| `matrices.linear_solve.8` | 95.177 | 16.563 |
| `matrices.multiply_symbolic.12` | 159.696 | 494.334 |
| `matrices.multiply_symbolic.6` | 226.413 | 10.146 |
| `matrices.nullspace.12` | 13.958 | 20.758 |
| `matrices.power_symbolic.4` | 134.197 | 3.592 |
| `matrices.rref_deficient.12` | 3.284 | 3.033 |
| `matrices.rref_integer.10` | 6.005 | 17.495 |
| `matrices.rref_integer.30` | 9.539 | 12.272 |
| `number_theory.factorint` | 2.663 | 46.637 |
| `number_theory.perfect_power` | 0.359 | 16.018 |
| `number_theory.primepi.100000` | 10.020 | 0.491 |
| `number_theory.primepi.1000000` | 12.670 | 0.493 |
| `numerics.complex_log.100` | 12.730 | 125.509 |
| `numerics.complex_log.50` | 17.773 | 115.878 |
| `numerics.constant.100` | 0.902 | 162.826 |
| `numerics.constant.300` | 0.775 | 284.445 |
| `numerics.constant.50` | 1.259 | 135.004 |
| `numerics.large_argument.100` | 0.586 | 44.323 |
| `numerics.large_argument.50` | 0.537 | 32.313 |
| `numerics.tiny_sine.100` | 0.699 | 58.839 |
| `numerics.tiny_sine.50` | 0.733 | 60.686 |
| `ode.harmonic` | 2.106 | 1.507 |
| `ode.linear_first` | 1.705 | 1.682 |
| `parsing.sympify.300` | 23.253 | 16.111 |
| `parsing.sympify.50` | 16.885 | 0.722 |
| `polynomials.apart.4` | 38.946 | 44.800 |
| `polynomials.apart.8` | 82.503 | 79.774 |
| `polynomials.cancel_sum.15` | 68.560 | 3.309 |
| `polynomials.cancel_sum.5` | 14.938 | 5.419 |
| `polynomials.factor_cyclotomic.100` | 27.107 | 48.924 |
| `polynomials.factor_cyclotomic.20` | 12.814 | 21.485 |
| `polynomials.factor_cyclotomic.50` | 14.250 | 34.561 |
| `polynomials.factor_modular.20` | 12.780 | 20.560 |
| `polynomials.factor_modular.50` | 9.357 | 10.530 |
| `polynomials.factor_multivariate.2` | 31.033 | 32.483 |
| `polynomials.factor_multivariate.4` | 36.793 | 39.878 |
| `polynomials.factor_wilkinson.10` | 15.979 | 18.682 |
| `polynomials.factor_wilkinson.20` | 19.797 | 20.448 |
| `polynomials.gcd.10` | 108.030 | 9.499 |
| `polynomials.gcd.30` | 85.112 | 6.670 |
| `polynomials.groebner.3` | 15.454 | 34.644 |
| `polynomials.groebner.4` | 19.954 | 39.772 |
| `polynomials.nroots.10` | 8.531 | 15.454 |
| `polynomials.nroots.30` | 39.716 | 41.209 |
| `polynomials.radical_division` | 73.515 | 32.310 |
| `polynomials.resultant.3` | 0.890 | 1.658 |
| `polynomials.resultant.6` | 5.180 | 2.638 |
| `polynomials.roots` | 19.484 | 17.020 |
| `printing.latex.5` | 79.219 | 48.952 |
| `printing.latex.8` | 18.817 | 62.599 |
| `printing.str.5` | 31.839 | 20.981 |
| `printing.str.8` | 101.668 | 24.142 |
| `series.rational.10` | 759.019 | 1196.932 |
| `series.rational.30` | 132.549 | 938.441 |
| `series.transcendental.10` | 3920.362 | 6849.222 |
| `simplification.powsimp.100` | 10.934 | 7.909 |
| `simplification.powsimp.20` | 10.285 | 9.913 |
| `simplification.radsimp.3` | 37.840 | 36.955 |
| `simplification.radsimp.5` | 70.684 | 240.135 |
| `simplification.rational.12` | 80.749 | 8.000 |
| `simplification.rational.5` | 79.169 | 2.589 |
| `simplification.tangent` | 100.884 | 4.637 |
| `simplification.trig_identity.4` | 97.490 | 13.173 |
| `simplification.trig_identity.8` | 48.935 | 12.884 |
| `simplification.trigsimp.2` | 45.918 | 33.440 |
| `simplification.trigsimp.4` | 2.097 | 20.349 |
| `simplification.trigsimp.8` | 59.702 | 20.054 |
| `solvers.linear_system.3` | 5.681 | 9.857 |
| `solvers.linear_system.8` | 3.554 | 5.542 |
| `solvers.nonlinear_system` | 222.750 | 7.030 |
| `solvers.nsolve` | 54.475 | 1.353 |
| `solvers.polynomial.4` | 10.523 | 1.096 |
| `solvers.polynomial.7` | 1.368 | 0.254 |
| `solvers.solveset_real` | 7.745 | 13.294 |
| `statistics.binomial_variance` | 17.969 | 5.135 |
| `statistics.normal_moment` | 117.370 | 96.273 |
| `substitution.subs.12` | 40.706 | 0.716 |
| `substitution.subs.8` | 25.100 | 1.544 |
| `substitution.xreplace.100` | 77.252 | 4.199 |
| `substitution.xreplace.500` | 28.779 | 23.378 |
| `summation.geometric.100` | 106.466 | 16.041 |
| `summation.geometric.10000` | 182.746 | 26.299 |
