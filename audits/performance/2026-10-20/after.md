# Public algorithm performance screening

Reference: SymPy 1.14.0. Symox distribution: 0.1.1.
Package fingerprint: `d1c023c144e08b3f3bb67d840f975ee733069263ba0ce8272507870019a43fa4`.

137 workloads; 132 verified and comparable; 5 excluded.
24 comparable workloads have a first-call or repeated-call speedup ≤ 1.5.

Speedup means reference time / native time; values below 1 mean Symox is slower.
Fresh processes alternate backend order on one CPU. Inputs and result verification are outside timing.
Repeated calls retain both libraries' caches. These measurements are not a universal parity or complexity certificate.
Three-process screenings are provisional; bootstrap intervals describe observed variation, not host isolation.

## Weak or marginal cases

| Workload | First speedup | Repeated speedup | Native first, ms | Native repeated, ms | Processes |
|---|---:|---:|---:|---:|---:|
| `construction.incremental_sum.100` | 87.648 | 0.554 | 0.741564 | 0.653203 | 7 |
| `expansion.univariate.30` | 59.004 | 0.594 | 0.157299 | 0.00231797 | 7 |
| `expansion.multinomial.12` | 76.099 | 0.607 | 1.04411 | 0.00230397 | 7 |
| `expansion.multinomial.8` | 68.059 | 0.611 | 0.404479 | 0.00228983 | 7 |
| `expansion.univariate.100` | 43.738 | 0.622 | 0.719692 | 0.00225096 | 7 |
| `parsing.sympify.50` | 7.531 | 0.670 | 7.18542 | 6.70182 | 7 |
| `number_theory.primepi.100000` | 10.421 | 0.707 | 0.074846 | 0.000386878 | 7 |
| `number_theory.primepi.1000000` | 12.780 | 0.707 | 0.145618 | 0.000387351 | 7 |
| `combinatorics.partition.50` | 3.670 | 0.710 | 0.060629 | 0.000383314 | 7 |
| `combinatorics.partition.200` | 5.988 | 0.713 | 0.09257 | 0.000383649 | 7 |
| `construction.bulk_add.100` | 104.734 | 0.810 | 0.198886 | 0.0097331 | 7 |
| `numerics.tiny_sine.50` | 0.910 | 88.617 | 0.077489 | 0.000497224 | 7 |
| `numerics.tiny_sine.100` | 0.919 | 88.783 | 0.075635 | 0.000503905 | 7 |
| `numerics.large_argument.100` | 0.932 | 64.614 | 0.145147 | 0.000487513 | 7 |
| `lambdify.evaluate.20` | 0.965 | 0.981 | 2.07569 | 1.93901 | 7 |
| `combinatorics.symmetric_group.12` | 93.313 | 0.984 | 0.001588 | 0.000118944 | 7 |
| `combinatorics.symmetric_group.8` | 77.950 | 0.985 | 0.001969 | 0.000119457 | 7 |
| `lambdify.evaluate.100` | 0.985 | 0.991 | 9.42731 | 9.49306 | 7 |
| `numerics.large_argument.50` | 1.004 | 47.538 | 0.123896 | 0.000511737 | 7 |
| `construction.bulk_add.1000` | 91.296 | 1.136 | 0.596155 | 0.0696089 | 7 |
| `statistics.binomial_variance` | 5.068 | 1.159 | 4.80478 | 4.3559 | 7 |
| `combinatorics.permutation.100` | 1.183 | 1.279 | 0.030988 | 0.0182479 | 7 |
| `combinatorics.permutation.1000` | 1.287 | 1.279 | 0.15824 | 0.146914 | 7 |
| `matrices.charpoly_integer.6` | 1.308 | 4.520 | 0.392191 | 0.0737768 | 7 |

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
| `assumptions.tiny_positive` | 5.759 | 1.528 |
| `combinatorics.partition.200` | 5.988 | 0.713 |
| `combinatorics.partition.50` | 3.670 | 0.710 |
| `combinatorics.permutation.100` | 1.183 | 1.279 |
| `combinatorics.permutation.1000` | 1.287 | 1.279 |
| `combinatorics.symmetric_group.12` | 93.313 | 0.984 |
| `combinatorics.symmetric_group.8` | 77.950 | 0.985 |
| `construction.bulk_add.100` | 104.734 | 0.810 |
| `construction.bulk_add.1000` | 91.296 | 1.136 |
| `construction.incremental_sum.100` | 87.648 | 0.554 |
| `construction.incremental_sum.500` | 94.224 | 85.598 |
| `differentiation.composite.1` | 144.589 | 8.540 |
| `differentiation.composite.4` | 285.941 | 115.750 |
| `expansion.multinomial.12` | 76.099 | 0.607 |
| `expansion.multinomial.16` | 76.404 | 2.505 |
| `expansion.multinomial.8` | 68.059 | 0.611 |
| `expansion.univariate.100` | 43.738 | 0.622 |
| `expansion.univariate.30` | 59.004 | 0.594 |
| `integration.dirac_affine` | 242.096 | 552.361 |
| `integration.exponential_trig` | 1261.355 | 811.341 |
| `integration.gamma_moment.10` | 589.490 | 97.278 |
| `integration.gamma_moment.50` | 651.766 | 82.948 |
| `integration.polynomial.100` | 336.236 | 11.922 |
| `integration.polynomial.30` | 461.322 | 20.408 |
| `integration.rational.4` | 447.747 | 422.054 |
| `lambdify.compile_and_evaluate.100` | 5.711 | 5.634 |
| `lambdify.compile_and_evaluate.20` | 13.630 | 4.729 |
| `lambdify.evaluate.100` | 0.985 | 0.991 |
| `lambdify.evaluate.20` | 0.965 | 0.981 |
| `limits.at_infinity` | 700.210 | 1312.782 |
| `limits.removable` | 545.321 | 8.499 |
| `matrices.charpoly_generic.3` | 2.349 | 2.271 |
| `matrices.charpoly_generic.4` | 2.676 | 1.821 |
| `matrices.charpoly_generic.5` | 1.869 | 1.531 |
| `matrices.charpoly_generic.6` | 2.116 | 2.419 |
| `matrices.charpoly_integer.12` | 2.836 | 6.482 |
| `matrices.charpoly_integer.6` | 1.308 | 4.520 |
| `matrices.charpoly_shared.4` | 3.906 | 3.636 |
| `matrices.det_generic.3` | 248.351 | 6.472 |
| `matrices.det_generic.4` | 595.783 | 619.902 |
| `matrices.det_integer.20` | 338.303 | 658.568 |
| `matrices.det_integer.8` | 106.627 | 447.548 |
| `matrices.eigenvalues` | 315.489 | 236.617 |
| `matrices.inverse_generic.3` | 139.035 | 55.634 |
| `matrices.inverse_rational.10` | 8.671 | 8.619 |
| `matrices.inverse_rational.20` | 4.662 | 4.356 |
| `matrices.inverse_rational.4` | 13.898 | 73.954 |
| `matrices.jacobian.16` | 564.791 | 74.033 |
| `matrices.jacobian.8` | 603.212 | 53.429 |
| `matrices.linear_solve.20` | 13.818 | 9.838 |
| `matrices.linear_solve.8` | 95.445 | 17.321 |
| `matrices.multiply_symbolic.12` | 114.229 | 70.880 |
| `matrices.multiply_symbolic.6` | 139.654 | 2.754 |
| `matrices.nullspace.12` | 14.370 | 20.997 |
| `matrices.power_symbolic.4` | 69.552 | 3.926 |
| `matrices.rref_deficient.12` | 3.463 | 3.279 |
| `matrices.rref_integer.10` | 7.215 | 19.026 |
| `matrices.rref_integer.30` | 9.599 | 13.051 |
| `number_theory.factorint` | 2.662 | 47.650 |
| `number_theory.perfect_power` | 1.652 | 16.650 |
| `number_theory.primepi.100000` | 10.421 | 0.707 |
| `number_theory.primepi.1000000` | 12.780 | 0.707 |
| `numerics.complex_log.100` | 21.089 | 174.834 |
| `numerics.complex_log.50` | 33.569 | 164.480 |
| `numerics.constant.100` | 2.313 | 227.708 |
| `numerics.constant.300` | 1.529 | 415.623 |
| `numerics.constant.50` | 2.361 | 190.355 |
| `numerics.large_argument.100` | 0.932 | 64.614 |
| `numerics.large_argument.50` | 1.004 | 47.538 |
| `numerics.tiny_sine.100` | 0.919 | 88.783 |
| `numerics.tiny_sine.50` | 0.910 | 88.617 |
| `ode.harmonic` | 122.297 | 201.653 |
| `ode.linear_first` | 108.541 | 225.587 |
| `parsing.sympify.300` | 10.376 | 10.917 |
| `parsing.sympify.50` | 7.531 | 0.670 |
| `polynomials.apart.4` | 86.168 | 126.113 |
| `polynomials.apart.8` | 182.808 | 185.123 |
| `polynomials.cancel_sum.15` | 222.607 | 12.614 |
| `polynomials.cancel_sum.5` | 74.992 | 36.356 |
| `polynomials.factor_cyclotomic.100` | 32.592 | 50.325 |
| `polynomials.factor_cyclotomic.20` | 16.577 | 27.800 |
| `polynomials.factor_cyclotomic.50` | 22.009 | 40.430 |
| `polynomials.factor_modular.20` | 14.831 | 24.285 |
| `polynomials.factor_modular.50` | 10.953 | 11.872 |
| `polynomials.factor_multivariate.2` | 35.875 | 35.951 |
| `polynomials.factor_multivariate.4` | 40.779 | 42.971 |
| `polynomials.factor_wilkinson.10` | 19.328 | 21.061 |
| `polynomials.factor_wilkinson.20` | 34.263 | 33.724 |
| `polynomials.gcd.10` | 146.027 | 15.213 |
| `polynomials.gcd.30` | 109.340 | 8.531 |
| `polynomials.groebner.3` | 20.019 | 69.903 |
| `polynomials.groebner.4` | 30.678 | 54.451 |
| `polynomials.nroots.10` | 19.113 | 18.126 |
| `polynomials.nroots.30` | 40.790 | 39.861 |
| `polynomials.radical_division` | 170.894 | 89.339 |
| `polynomials.resultant.3` | 4.870 | 2.322 |
| `polynomials.resultant.6` | 5.944 | 3.374 |
| `polynomials.roots` | 34.866 | 30.023 |
| `printing.latex.5` | 56.996 | 76.734 |
| `printing.latex.8` | 76.309 | 69.065 |
| `printing.str.5` | 23.524 | 24.068 |
| `printing.str.8` | 25.049 | 21.757 |
| `series.rational.10` | 1278.569 | 3540.528 |
| `series.rational.30` | 1115.592 | 1923.181 |
| `series.transcendental.10` | 4735.019 | 10336.077 |
| `simplification.powsimp.100` | 12.127 | 8.835 |
| `simplification.powsimp.20` | 10.502 | 10.673 |
| `simplification.radsimp.3` | 49.719 | 47.029 |
| `simplification.radsimp.5` | 103.430 | 351.481 |
| `simplification.rational.12` | 263.392 | 57.362 |
| `simplification.rational.5` | 261.481 | 46.163 |
| `simplification.tangent` | 122.625 | 11.280 |
| `simplification.trig_identity.4` | 119.863 | 31.494 |
| `simplification.trig_identity.8` | 59.853 | 15.931 |
| `simplification.trigsimp.2` | 31.195 | 29.958 |
| `simplification.trigsimp.4` | 31.105 | 30.364 |
| `simplification.trigsimp.8` | 30.161 | 30.561 |
| `solvers.linear_system.3` | 28.777 | 24.507 |
| `solvers.linear_system.8` | 13.939 | 8.099 |
| `solvers.nonlinear_system` | 388.102 | 37.082 |
| `solvers.nsolve` | 67.974 | 3.627 |
| `solvers.polynomial.4` | 17.060 | 29.561 |
| `solvers.polynomial.7` | 11.723 | 8.615 |
| `solvers.solveset_real` | 76.374 | 38.202 |
| `statistics.binomial_variance` | 5.068 | 1.159 |
| `statistics.normal_moment` | 91.712 | 134.604 |
| `substitution.subs.12` | 22.258 | 6.722 |
| `substitution.subs.8` | 22.961 | 6.724 |
| `substitution.xreplace.100` | 35.937 | 2.088 |
| `substitution.xreplace.500` | 39.879 | 29.345 |
| `summation.geometric.100` | 161.695 | 23.216 |
| `summation.geometric.10000` | 230.580 | 128.591 |
