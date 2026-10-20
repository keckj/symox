# Public algorithm performance screening

Reference: SymPy 1.14.0. Symox distribution: 0.1.1.
Package fingerprint: `3320a259cd9b79511eb133750086d0ca9925fd658fb3c9275617db257927b1ab`.

137 workloads; 132 verified and comparable; 5 excluded.
16 comparable workloads have a first-call or repeated-call speedup ≤ 1.5.

Speedup means reference time / native time; values below 1 mean Symox is slower.
Fresh processes alternate backend order on one CPU. Inputs and result verification are outside timing.
Repeated calls retain both libraries' caches. These measurements are not a universal parity or complexity certificate.
Three-process screenings are provisional; bootstrap intervals describe observed variation, not host isolation.

## Weak or marginal cases

| Workload | First speedup | Repeated speedup | Native first, ms | Native repeated, ms | Processes |
|---|---:|---:|---:|---:|---:|
| `number_theory.primepi.1000000` | 13.183 | 0.706 | 0.142947 | 0.000388784 | 7 |
| `combinatorics.partition.200` | 6.247 | 0.708 | 0.098659 | 0.000386951 | 7 |
| `combinatorics.partition.50` | 2.959 | 0.708 | 0.067263 | 0.000387027 | 7 |
| `number_theory.primepi.100000` | 11.019 | 0.710 | 0.072011 | 0.000389611 | 7 |
| `numerics.tiny_sine.50` | 0.964 | 85.955 | 0.075741 | 0.000514589 | 7 |
| `lambdify.evaluate.20` | 0.975 | 0.985 | 2.12717 | 1.98428 | 7 |
| `combinatorics.symmetric_group.8` | 89.682 | 0.980 | 0.001762 | 0.000120253 | 7 |
| `lambdify.evaluate.100` | 0.981 | 0.998 | 9.68286 | 9.54001 | 7 |
| `combinatorics.symmetric_group.12` | 134.151 | 0.983 | 0.00164 | 0.000120214 | 7 |
| `numerics.tiny_sine.100` | 1.003 | 91.892 | 0.072891 | 0.000496739 | 7 |
| `numerics.large_argument.100` | 1.022 | 62.768 | 0.131826 | 0.000503497 | 7 |
| `numerics.large_argument.50` | 1.023 | 48.912 | 0.118498 | 0.000525527 | 7 |
| `statistics.binomial_variance` | 5.603 | 1.247 | 4.44811 | 3.99317 | 7 |
| `matrices.charpoly_integer.6` | 1.332 | 4.572 | 0.389784 | 0.0741227 | 7 |
| `parsing.sympify.50` | 10.778 | 1.471 | 4.9908 | 3.00919 | 7 |
| `numerics.constant.300` | 1.484 | 408.027 | 0.638932 | 0.000504162 | 7 |

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
| `assumptions.tiny_positive` | 7.259 | 1.519 |
| `combinatorics.partition.200` | 6.247 | 0.708 |
| `combinatorics.partition.50` | 2.959 | 0.708 |
| `combinatorics.permutation.100` | 1.559 | 2.162 |
| `combinatorics.permutation.1000` | 1.710 | 1.952 |
| `combinatorics.symmetric_group.12` | 134.151 | 0.983 |
| `combinatorics.symmetric_group.8` | 89.682 | 0.980 |
| `construction.bulk_add.100` | 110.276 | 2.053 |
| `construction.bulk_add.1000` | 103.717 | 4.765 |
| `construction.incremental_sum.100` | 102.832 | 1.830 |
| `construction.incremental_sum.500` | 102.255 | 140.150 |
| `differentiation.composite.1` | 214.862 | 10.358 |
| `differentiation.composite.4` | 304.442 | 114.513 |
| `expansion.multinomial.12` | 78.493 | 2.370 |
| `expansion.multinomial.16` | 75.432 | 10.585 |
| `expansion.multinomial.8` | 73.301 | 2.385 |
| `expansion.univariate.100` | 46.422 | 2.342 |
| `expansion.univariate.30` | 73.826 | 2.485 |
| `integration.dirac_affine` | 262.816 | 459.240 |
| `integration.exponential_trig` | 1324.285 | 784.637 |
| `integration.gamma_moment.10` | 658.879 | 98.492 |
| `integration.gamma_moment.50` | 622.148 | 85.941 |
| `integration.polynomial.100` | 308.435 | 11.900 |
| `integration.polynomial.30` | 443.183 | 21.981 |
| `integration.rational.4` | 443.795 | 399.304 |
| `lambdify.compile_and_evaluate.100` | 6.275 | 5.677 |
| `lambdify.compile_and_evaluate.20` | 15.215 | 4.822 |
| `lambdify.evaluate.100` | 0.981 | 0.998 |
| `lambdify.evaluate.20` | 0.975 | 0.985 |
| `limits.at_infinity` | 694.804 | 1825.568 |
| `limits.removable` | 612.707 | 10.295 |
| `matrices.charpoly_generic.3` | 2.307 | 2.238 |
| `matrices.charpoly_generic.4` | 2.569 | 1.674 |
| `matrices.charpoly_generic.5` | 1.801 | 1.522 |
| `matrices.charpoly_generic.6` | 2.242 | 2.395 |
| `matrices.charpoly_integer.12` | 2.903 | 6.587 |
| `matrices.charpoly_integer.6` | 1.332 | 4.572 |
| `matrices.charpoly_shared.4` | 3.826 | 3.815 |
| `matrices.det_generic.3` | 301.425 | 6.574 |
| `matrices.det_generic.4` | 676.420 | 627.379 |
| `matrices.det_integer.20` | 348.855 | 681.043 |
| `matrices.det_integer.8` | 101.524 | 439.485 |
| `matrices.eigenvalues` | 334.260 | 230.874 |
| `matrices.inverse_generic.3` | 146.521 | 54.605 |
| `matrices.inverse_rational.10` | 8.555 | 8.596 |
| `matrices.inverse_rational.20` | 4.731 | 4.348 |
| `matrices.inverse_rational.4` | 12.557 | 72.063 |
| `matrices.jacobian.16` | 477.238 | 70.723 |
| `matrices.jacobian.8` | 569.648 | 51.449 |
| `matrices.linear_solve.20` | 13.260 | 9.753 |
| `matrices.linear_solve.8` | 92.573 | 16.791 |
| `matrices.multiply_symbolic.12` | 92.042 | 234.297 |
| `matrices.multiply_symbolic.6` | 125.354 | 8.313 |
| `matrices.nullspace.12` | 16.175 | 21.023 |
| `matrices.power_symbolic.4` | 67.586 | 3.804 |
| `matrices.rref_deficient.12` | 3.146 | 3.227 |
| `matrices.rref_integer.10` | 7.045 | 20.321 |
| `matrices.rref_integer.30` | 9.959 | 12.964 |
| `number_theory.factorint` | 2.683 | 49.865 |
| `number_theory.perfect_power` | 1.527 | 16.295 |
| `number_theory.primepi.100000` | 11.019 | 0.710 |
| `number_theory.primepi.1000000` | 13.183 | 0.706 |
| `numerics.complex_log.100` | 20.275 | 173.091 |
| `numerics.complex_log.50` | 31.104 | 175.039 |
| `numerics.constant.100` | 2.207 | 234.464 |
| `numerics.constant.300` | 1.484 | 408.027 |
| `numerics.constant.50` | 2.214 | 198.755 |
| `numerics.large_argument.100` | 1.022 | 62.768 |
| `numerics.large_argument.50` | 1.023 | 48.912 |
| `numerics.tiny_sine.100` | 1.003 | 91.892 |
| `numerics.tiny_sine.50` | 0.964 | 85.955 |
| `ode.harmonic` | 118.637 | 176.868 |
| `ode.linear_first` | 128.037 | 209.293 |
| `parsing.sympify.300` | 14.757 | 16.618 |
| `parsing.sympify.50` | 10.778 | 1.471 |
| `polynomials.apart.4` | 86.398 | 122.951 |
| `polynomials.apart.8` | 183.329 | 179.759 |
| `polynomials.cancel_sum.15` | 230.480 | 12.917 |
| `polynomials.cancel_sum.5` | 80.829 | 38.457 |
| `polynomials.factor_cyclotomic.100` | 36.053 | 54.255 |
| `polynomials.factor_cyclotomic.20` | 14.583 | 29.999 |
| `polynomials.factor_cyclotomic.50` | 23.427 | 44.432 |
| `polynomials.factor_modular.20` | 15.742 | 26.868 |
| `polynomials.factor_modular.50` | 12.261 | 12.890 |
| `polynomials.factor_multivariate.2` | 36.281 | 35.591 |
| `polynomials.factor_multivariate.4` | 39.826 | 44.924 |
| `polynomials.factor_wilkinson.10` | 19.136 | 21.606 |
| `polynomials.factor_wilkinson.20` | 34.791 | 34.741 |
| `polynomials.gcd.10` | 145.851 | 15.914 |
| `polynomials.gcd.30` | 108.192 | 8.832 |
| `polynomials.groebner.3` | 21.822 | 71.717 |
| `polynomials.groebner.4` | 33.803 | 54.069 |
| `polynomials.nroots.10` | 19.847 | 17.748 |
| `polynomials.nroots.30` | 41.552 | 40.866 |
| `polynomials.radical_division` | 161.663 | 82.169 |
| `polynomials.resultant.3` | 4.946 | 2.345 |
| `polynomials.resultant.6` | 6.143 | 3.369 |
| `polynomials.roots` | 35.948 | 28.484 |
| `printing.latex.5` | 54.034 | 77.654 |
| `printing.latex.8` | 75.022 | 69.503 |
| `printing.str.5` | 22.819 | 23.627 |
| `printing.str.8` | 21.441 | 20.443 |
| `series.rational.10` | 1063.737 | 3318.196 |
| `series.rational.30` | 1106.008 | 1982.556 |
| `series.transcendental.10` | 4088.412 | 10359.249 |
| `simplification.powsimp.100` | 11.907 | 8.658 |
| `simplification.powsimp.20` | 10.947 | 10.395 |
| `simplification.radsimp.3` | 49.765 | 46.181 |
| `simplification.radsimp.5` | 111.297 | 368.206 |
| `simplification.rational.12` | 270.465 | 63.513 |
| `simplification.rational.5` | 282.552 | 65.461 |
| `simplification.tangent` | 117.154 | 11.580 |
| `simplification.trig_identity.4` | 120.857 | 32.494 |
| `simplification.trig_identity.8` | 60.103 | 15.705 |
| `simplification.trigsimp.2` | 30.751 | 31.263 |
| `simplification.trigsimp.4` | 31.691 | 32.317 |
| `simplification.trigsimp.8` | 30.533 | 30.626 |
| `solvers.linear_system.3` | 25.912 | 22.444 |
| `solvers.linear_system.8` | 13.146 | 7.522 |
| `solvers.nonlinear_system` | 296.957 | 40.375 |
| `solvers.nsolve` | 66.968 | 3.735 |
| `solvers.polynomial.4` | 19.053 | 28.769 |
| `solvers.polynomial.7` | 12.680 | 8.384 |
| `solvers.solveset_real` | 78.597 | 38.910 |
| `statistics.binomial_variance` | 5.603 | 1.247 |
| `statistics.normal_moment` | 99.000 | 155.685 |
| `substitution.subs.12` | 27.686 | 9.026 |
| `substitution.subs.8` | 28.598 | 9.287 |
| `substitution.xreplace.100` | 47.581 | 2.042 |
| `substitution.xreplace.500` | 47.010 | 22.839 |
| `summation.geometric.100` | 162.213 | 24.704 |
| `summation.geometric.10000` | 242.080 | 132.957 |
