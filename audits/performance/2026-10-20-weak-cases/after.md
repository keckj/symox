# Public algorithm performance screening

Reference: SymPy 1.14.0. Symox distribution: 0.1.1.
Package fingerprint: `6f23f3f8e68e201204afa14900aeedbb0b6a4ca554d61b9839c748763040550b`.

137 workloads; 132 verified and comparable; 5 excluded.
17 comparable workloads have a first-call or repeated-call speedup ≤ 1.5.

Speedup means reference time / native time; values below 1 mean Symox is slower.
Fresh processes alternate backend order on one CPU. Inputs and result verification are outside timing.
Repeated calls retain both libraries' caches. These measurements are not a universal parity or complexity certificate.
Three-process screenings are provisional; bootstrap intervals describe observed variation, not host isolation.

## Weak or marginal cases

| Workload | First speedup | Repeated speedup | Native first, ms | Native repeated, ms | Processes |
|---|---:|---:|---:|---:|---:|
| `number_theory.primepi.1000000` | 12.703 | 0.701 | 0.138766 | 0.00037742 | 7 |
| `number_theory.primepi.100000` | 10.227 | 0.705 | 0.068279 | 0.000377567 | 7 |
| `combinatorics.partition.200` | 5.476 | 0.707 | 0.088988 | 0.000381459 | 7 |
| `combinatorics.partition.50` | 4.420 | 0.710 | 0.056563 | 0.000378797 | 7 |
| `numerics.tiny_sine.50` | 0.896 | 90.663 | 0.071536 | 0.0004701 | 7 |
| `numerics.tiny_sine.100` | 0.908 | 90.615 | 0.075172 | 0.000486971 | 7 |
| `construction.incremental_sum.100` | 148.415 | 0.936 | 0.431213 | 0.37258 | 7 |
| `lambdify.evaluate.20` | 0.958 | 0.970 | 2.01358 | 1.89315 | 7 |
| `combinatorics.symmetric_group.8` | 92.924 | 0.992 | 0.001508 | 0.000117344 | 7 |
| `lambdify.evaluate.100` | 0.994 | 0.994 | 9.12111 | 9.1608 | 7 |
| `combinatorics.symmetric_group.12` | 124.845 | 0.995 | 0.001461 | 0.000113414 | 7 |
| `numerics.large_argument.50` | 0.998 | 49.545 | 0.120043 | 0.000475606 | 7 |
| `numerics.large_argument.100` | 1.108 | 62.578 | 0.122818 | 0.000486701 | 7 |
| `statistics.binomial_variance` | 5.579 | 1.242 | 4.22553 | 3.79645 | 7 |
| `matrices.charpoly_integer.6` | 1.340 | 4.501 | 0.358006 | 0.0723693 | 7 |
| `parsing.sympify.50` | 10.747 | 1.442 | 4.88047 | 2.94195 | 7 |
| `number_theory.perfect_power` | 1.497 | 15.813 | 0.056875 | 0.00255881 | 7 |

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
| `assumptions.tiny_positive` | 6.161 | 1.536 |
| `combinatorics.partition.200` | 5.476 | 0.707 |
| `combinatorics.partition.50` | 4.420 | 0.710 |
| `combinatorics.permutation.100` | 1.638 | 2.146 |
| `combinatorics.permutation.1000` | 1.682 | 1.952 |
| `combinatorics.symmetric_group.12` | 124.845 | 0.995 |
| `combinatorics.symmetric_group.8` | 92.924 | 0.992 |
| `construction.bulk_add.100` | 109.865 | 2.063 |
| `construction.bulk_add.1000` | 104.593 | 4.828 |
| `construction.incremental_sum.100` | 148.415 | 0.936 |
| `construction.incremental_sum.500` | 128.023 | 116.011 |
| `differentiation.composite.1` | 173.558 | 10.143 |
| `differentiation.composite.4` | 297.028 | 119.330 |
| `expansion.multinomial.12` | 75.329 | 2.387 |
| `expansion.multinomial.16` | 77.239 | 9.683 |
| `expansion.multinomial.8` | 69.410 | 2.419 |
| `expansion.univariate.100` | 45.246 | 2.442 |
| `expansion.univariate.30` | 67.300 | 2.471 |
| `integration.dirac_affine` | 248.498 | 531.374 |
| `integration.exponential_trig` | 1183.238 | 872.875 |
| `integration.gamma_moment.10` | 652.170 | 105.864 |
| `integration.gamma_moment.50` | 575.518 | 92.043 |
| `integration.polynomial.100` | 326.761 | 12.192 |
| `integration.polynomial.30` | 453.589 | 20.587 |
| `integration.rational.4` | 441.431 | 416.647 |
| `lambdify.compile_and_evaluate.100` | 5.957 | 5.400 |
| `lambdify.compile_and_evaluate.20` | 14.235 | 4.758 |
| `lambdify.evaluate.100` | 0.994 | 0.994 |
| `lambdify.evaluate.20` | 0.958 | 0.970 |
| `limits.at_infinity` | 656.850 | 1584.746 |
| `limits.removable` | 579.551 | 10.852 |
| `matrices.charpoly_generic.3` | 2.199 | 2.288 |
| `matrices.charpoly_generic.4` | 2.542 | 1.751 |
| `matrices.charpoly_generic.5` | 1.875 | 1.583 |
| `matrices.charpoly_generic.6` | 2.437 | 2.659 |
| `matrices.charpoly_integer.12` | 2.834 | 6.522 |
| `matrices.charpoly_integer.6` | 1.340 | 4.501 |
| `matrices.charpoly_shared.4` | 4.066 | 3.744 |
| `matrices.det_generic.3` | 257.195 | 6.945 |
| `matrices.det_generic.4` | 654.563 | 633.452 |
| `matrices.det_integer.20` | 336.369 | 640.164 |
| `matrices.det_integer.8` | 114.032 | 433.823 |
| `matrices.eigenvalues` | 331.258 | 232.943 |
| `matrices.inverse_generic.3` | 135.827 | 54.629 |
| `matrices.inverse_rational.10` | 8.033 | 8.302 |
| `matrices.inverse_rational.20` | 4.665 | 4.333 |
| `matrices.inverse_rational.4` | 13.164 | 73.924 |
| `matrices.jacobian.16` | 632.524 | 74.057 |
| `matrices.jacobian.8` | 562.363 | 52.443 |
| `matrices.linear_solve.20` | 14.182 | 10.060 |
| `matrices.linear_solve.8` | 92.057 | 16.972 |
| `matrices.multiply_symbolic.12` | 113.985 | 70.155 |
| `matrices.multiply_symbolic.6` | 136.345 | 2.820 |
| `matrices.nullspace.12` | 14.614 | 21.133 |
| `matrices.power_symbolic.4` | 70.678 | 3.924 |
| `matrices.rref_deficient.12` | 3.192 | 3.235 |
| `matrices.rref_integer.10` | 6.888 | 18.761 |
| `matrices.rref_integer.30` | 9.192 | 12.842 |
| `number_theory.factorint` | 2.664 | 46.933 |
| `number_theory.perfect_power` | 1.497 | 15.813 |
| `number_theory.primepi.100000` | 10.227 | 0.705 |
| `number_theory.primepi.1000000` | 12.703 | 0.701 |
| `numerics.complex_log.100` | 21.037 | 182.923 |
| `numerics.complex_log.50` | 33.747 | 165.753 |
| `numerics.constant.100` | 2.204 | 237.625 |
| `numerics.constant.300` | 1.579 | 423.476 |
| `numerics.constant.50` | 2.276 | 202.453 |
| `numerics.large_argument.100` | 1.108 | 62.578 |
| `numerics.large_argument.50` | 0.998 | 49.545 |
| `numerics.tiny_sine.100` | 0.908 | 90.615 |
| `numerics.tiny_sine.50` | 0.896 | 90.663 |
| `ode.harmonic` | 122.504 | 204.960 |
| `ode.linear_first` | 125.696 | 232.865 |
| `parsing.sympify.300` | 15.064 | 16.238 |
| `parsing.sympify.50` | 10.747 | 1.442 |
| `polynomials.apart.4` | 92.546 | 130.651 |
| `polynomials.apart.8` | 179.897 | 187.170 |
| `polynomials.cancel_sum.15` | 229.464 | 13.326 |
| `polynomials.cancel_sum.5` | 78.395 | 39.714 |
| `polynomials.factor_cyclotomic.100` | 34.423 | 50.751 |
| `polynomials.factor_cyclotomic.20` | 16.432 | 28.947 |
| `polynomials.factor_cyclotomic.50` | 23.515 | 42.189 |
| `polynomials.factor_modular.20` | 15.427 | 27.128 |
| `polynomials.factor_modular.50` | 11.800 | 12.427 |
| `polynomials.factor_multivariate.2` | 34.799 | 37.187 |
| `polynomials.factor_multivariate.4` | 41.001 | 43.040 |
| `polynomials.factor_wilkinson.10` | 19.503 | 21.628 |
| `polynomials.factor_wilkinson.20` | 33.256 | 33.823 |
| `polynomials.gcd.10` | 140.545 | 16.116 |
| `polynomials.gcd.30` | 110.950 | 8.506 |
| `polynomials.groebner.3` | 20.306 | 70.133 |
| `polynomials.groebner.4` | 32.375 | 53.407 |
| `polynomials.nroots.10` | 19.371 | 17.853 |
| `polynomials.nroots.30` | 41.204 | 41.095 |
| `polynomials.radical_division` | 167.604 | 85.235 |
| `polynomials.resultant.3` | 4.831 | 2.291 |
| `polynomials.resultant.6` | 6.227 | 3.295 |
| `polynomials.roots` | 33.465 | 28.680 |
| `printing.latex.5` | 58.180 | 76.471 |
| `printing.latex.8` | 78.928 | 70.481 |
| `printing.str.5` | 23.364 | 24.245 |
| `printing.str.8` | 25.546 | 22.005 |
| `series.rational.10` | 1138.103 | 3451.545 |
| `series.rational.30` | 1147.256 | 1969.085 |
| `series.transcendental.10` | 4765.417 | 10924.239 |
| `simplification.powsimp.100` | 12.921 | 8.745 |
| `simplification.powsimp.20` | 11.148 | 10.593 |
| `simplification.radsimp.3` | 52.855 | 47.007 |
| `simplification.radsimp.5` | 104.331 | 357.772 |
| `simplification.rational.12` | 243.705 | 59.672 |
| `simplification.rational.5` | 273.445 | 50.407 |
| `simplification.tangent` | 117.278 | 11.303 |
| `simplification.trig_identity.4` | 120.315 | 32.438 |
| `simplification.trig_identity.8` | 59.583 | 15.670 |
| `simplification.trigsimp.2` | 31.525 | 31.125 |
| `simplification.trigsimp.4` | 30.060 | 30.607 |
| `simplification.trigsimp.8` | 30.234 | 30.690 |
| `solvers.linear_system.3` | 26.990 | 25.092 |
| `solvers.linear_system.8` | 13.941 | 8.318 |
| `solvers.nonlinear_system` | 414.638 | 37.486 |
| `solvers.nsolve` | 66.050 | 3.624 |
| `solvers.polynomial.4` | 17.883 | 29.466 |
| `solvers.polynomial.7` | 11.991 | 8.511 |
| `solvers.solveset_real` | 75.288 | 39.331 |
| `statistics.binomial_variance` | 5.579 | 1.242 |
| `statistics.normal_moment` | 99.585 | 142.843 |
| `substitution.subs.12` | 28.405 | 8.435 |
| `substitution.subs.8` | 29.166 | 8.285 |
| `substitution.xreplace.100` | 49.247 | 2.135 |
| `substitution.xreplace.500` | 50.484 | 26.109 |
| `summation.geometric.100` | 157.802 | 24.478 |
| `summation.geometric.10000` | 244.232 | 149.893 |
