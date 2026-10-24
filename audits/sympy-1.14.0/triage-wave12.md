# Triage wave 12 (audit12): c lots c1, c2 and c4

Replay (`tools/sympy_parity.py examples --timeout 30`) of every object of the
c lots of `triage-wave4-after.json`: c1 (stats, holonomic, sums, series,
integrals, calculus: 36 objects), c2 (functions, core, sets, solvers,
simplify, logic, assumptions: 40) and c4 (physics, tensor, vector,
categories, diffgeom, geometry: 37).

| lot | before (match/total) | after |
|---|---:|---:|
| c1 | 27/36 | 27/36 |
| c2 | 30/40 | 31/40 |
| c4 | 37/37 | 37/37 |

The gaps named in TODO.md were already closed by earlier waves when this
replay ran: `to_meijerg`, `fps` with a symbolic power prefactor, the
stochastic-process doctests (BernoulliProcess, ContinuousMarkovChain,
GammaProcess and the `dependent`/`independent` examples), multiple
derivatives ended with `factor_terms(signsimp(...))`, `nseries` of `x**y`,
`trigsimp(old=True)`, `devise_plan_meijer`, `solve_triangulated` with
`extension=True`, and every c4 doctest (the Python-level pretty printing).

Still-differing objects. a = printing/order only, b = equivalent form,
c = missing feature, h = SymPy is wrong and symox keeps the correct result,
i = environment (Python set or dict order).

| object | before | after | note |
|---|---|---|---|
| expand | b | match | power_base is the reference `_eval_expand_power_base` (nonnegative factors split for any exponent), applied after mul, factors of a product visited as powers. |
| Derivative | b | b | [6] variables are not sorted canonically (`_sort_variable_count`): `Derivative(f(x, y), x, y)` vs `(y, x)`. [8] the chain rule keeps `Subs(Derivative(f(_xi_1), _xi_1), _xi_1, g(x))` where 1.14 writes `Derivative(f(g(x)), g(x))`. |
| Predicate | a | a | `type(Q.P)` is not `UndefinedPredicate`; module name in the error message. |
| symbols | i + a | i + a | set order; class of `symbols('f', cls=Function)` items. |
| ProductSet | i | i | Python set order. |
| classify_sysode | i | i | dict key order of `order`. |
| diophantine | i + a | i + a | set order; the `x = 0` factor gives `(0, t_1, t_2)` where SymPy writes `(0, n1, n2)`. |
| nonlinsolve | a | a | free parameter shown as `Complexes \ {0}` where SymPy writes `{d} \ {0}`. |
| harmonic | b + c | b + c | unchanged (see audit9). |
| assoc_laguerre | h | h | unchanged (see audit9). |
| Sum.euler_maclaurin | h | h | unchanged. |
| approximants | h | h | unchanged. |
| apply_finite_diff | a | a | term order. |
| stats.Expectation | a | a | factor order `X*a`, nested argument order. |
| stats.Moyal | b | b | `erfc(...)` for `1 - erf(...)`. |
| stats.GaussianSymplecticEnsemble | b | b | closed normalization `32/pi**3` instead of the held `Integral`. |
| stats.WienerProcess | h | h | `P((X(t1) > 2) | (X(t1) < 4), ...)` is 1; SymPy gives `1 - erf(1)/2 + erf(2)/2 > 1`. |
| stats.DiscreteMarkovChain | c | c | `Str` states (documented). |
| stats.MultivariateNormal | c | c | `MatrixSymbol` parameters (documented). |
