# Ordinary and partial differential equations

[README](../README.md) · [Solution checks](solver-checks.md)

The differential equation interfaces do not import SymPy at runtime. Their
Python layer is adapted from SymPy 1.14 (see `LICENSE-SYMPY`) and runs on
the native expression, polynomial, integration and solving algorithms.

## Scalar ordinary differential equations

`dsolve`, `classify_ode`, `checkodesol`, `constantsimp`, `constant_renumber`,
`homogeneous_order`, `infinitesimals`, `checkinfsol` and `solve_ics` follow
SymPy's `sympy.solvers.ode` modules: `ode.py`, `single.py`, `riccati.py`,
`nonhomogeneous.py`, `hypergeometric.py`, `lie_group.py`, `subscheck.py`
and `systems.py`. Classification returns SymPy's hint tuple in SymPy's
order, `dict=True` returns the same matching dictionaries, and
`hint='all'`, `'all_Integral'`, `'best'`, `simplify=False`, `ics`, `xi`,
`eta`, `n` and `x0` have SymPy's meaning. Every hint of
`symox.solvers.ode.allhints` is available, including the Airy and Bessel
families (`2nd_linear_airy`, `2nd_linear_bessel`), the hypergeometric
reduction, Riccati special forms, the Liouville form, the nth-order
reducible form, undetermined coefficients, variation of parameters, the
Euler families, the power series and the Lie group method.

```python
from symox import Function, Symbol, airyai, cos, classify_ode, dsolve, checkodesol

x = Symbol('x')
f = Function('f')(x)
equation = (x + cos(f))*f.diff(x) + f + 2*x
assert '1st_exact' in classify_ode(equation, f)
solution = dsolve(equation, f, hint='1st_exact')
assert checkodesol(equation, solution, f, solve_for_func=False) == (True, 0)
assert dsolve(f.diff(x, 2) - x*f, f).rhs.has(airyai)
```

The result forms (constant names and numbering, sign and factor choices,
implicit or explicit solutions) are those of SymPy whenever the underlying
`solve`, `simplify`, `integrate` and `factor` produce SymPy's forms. On a
corpus of 315 equations drawn from SymPy's ODE test suite, the last
recorded run (2026-10-30) matched SymPy's output string for 303
classifications, 279 results of the test suite's hint and 274 default
`dsolve` results; the corpus has not been replayed on the later
integration waves. The classification mismatches of that run are SymPy
reference timeouts (run alone, SymPy gives the symox lists). The other
differences come from the underlying algorithms: the native `simplify` is
a candidate search, not SymPy's sequence of rewrites, so it can stop at
another equivalent form (this can also decide whether the `factorable`
hint is listed, since its match runs `solve` on the simplified form);
integration results and integrals left unevaluated; and root ordering and
radical forms of the characteristic roots. The corpus compares strings
only; it does not prove that each differing result is an equivalent
solution.

### Native kernels

The hint matchers and the linear solution steps run in Rust
(`crates/symox-core/src/ode/`), called from the reference Python bodies
through admission hooks (`_classify_native.py`, `_odesolve_native.py`):

- classification: the first-order matches (`separable`, `1st_exact`,
  `1st_linear`, `Bernoulli`, `almost_linear`, the homogeneous-coefficient
  and power series forms, `Riccati_special_minus2`, the `match_riccati`
  coefficients of `1st_rational_riccati`), the `nth_linear_*` and Euler
  detections with `get_linear_coefficients`, `match_2nd_hypergeometric`,
  the `2nd_linear_airy` match and `2nd_linear_bessel` rejection, the
  order-two power series branch, `nth_order_reducible`, `Liouville`,
  `2nd_nonlinear_autonomous_conserved` and the `allhints` ordering;
- solutions: the characteristic and indicial polynomials, their roots in
  the order of `roots(multiple=True)` (`rootof` when radicals do not give
  every root), the fundamental solutions and the collect/powsimp step of
  the constant-coefficient and Euler hints; the trial set and exact
  coefficient solve of undetermined coefficients; the Wronskians, held
  integrals and particular solution of variation of parameters.

Each kernel admits a documented exact domain (rational or symbolic
constant coefficients, forcing terms `q*x**k*exp(a*x)*sin|cos(b*x)` with
rational parameters, no Floats, `I` or unmodelled content) and declines
otherwise; a decline runs the reference body unchanged, never a no-match.
Still in Python: `factorable`, `nth_algebraic`, `lie_group`,
`equivalence_hypergeometric`/`match_2nd_2F1`, the positive Bessel match,
`eq_preprocessed`, the `simplify`/`trigsimp` of the particular solutions,
the evaluation of their integrals, `get_numbered_constants` and the
`odesimp`/`constantsimp` steps of `dsolve`. The domains and remaining
gaps are listed in [the parity notes](algorithmic-parity.md#odeclass-wave-35-native-classify_ode-matches-2026-10-30).

### Native fallback

When every reference method raises `NotImplementedError` for
`hint='default'`, `dsolve` tries the native Rust solvers (linear constant
coefficient and Cauchy-Euler equations, first-order linear, separable,
Bernoulli and exact equations, Lie point symmetries). If they also fail,
the reference `NotImplementedError` is raised again: a failure is never
reported as an absence of solutions. The bounds of the native solvers (at
most 64 series coefficients, matrix extraction of at most 64 functions,
matrix exponentials of at most eight functions, differential order 256)
apply only to this fallback; the reference methods have no such bounds,
so `n=100` series and large systems behave as in SymPy. `classify_ode`
raises `NotImplementedError` for differential orders above `10**6`, where
SymPy differentiates once per order and does not return.

### Deliberate extensions

These deviations return a correct result where SymPy returns an
incomplete one or fails:

- Initial conditions may use names that coincide with the generated
  constants (`ics={f(0): C1}` gives `C1*exp(x)` for `f' = f`). SymPy
  substitutes its own constants into the given values.
- `2nd_power_series_regular` returns a complete basis when the indicial
  roots coincide or differ by an integer with an inconsistent recurrence.
  The second solution is the Frobenius derivative
  `d/dr[(x - x0)**r * sum(b_j(r)*(x - x0)**j)]` at the smaller root, with
  `b_j = (r - m)**k * a_j` (`k = 0` for a repeated root, `k = 1` at a
  resonance), so it carries a `log(x - x0)` term. SymPy 1.14 returns only
  the first series in these cases. The reference expression is simplified
  exactly as in SymPy and `C2` times the logarithmic solution is appended
  afterwards, so `sol.subs(C2, 0)` is the SymPy output.
- Non-polynomial coefficients `p`, `q` are expanded in `x - x0`. SymPy
  1.14 calls `series(ser, x=x0, n=n+1)`, which does not expand them, and
  returns series that do not solve the equation (`2F1(1, 2; 1/2; x)`
  starts `1 - 4*x` there instead of `1 + 4*x`).
- Irrational indicial roots return the reference's raw series, where
  SymPy raises `NotImplementedError` from `Order.contains` during
  constant simplification.
- Characteristic roots in the native fallback take their multiplicities
  from an exact square-free decomposition, not a numerical derivative
  test that could merge close simple roots.

```python
from symox import log

eq = x**2*f.diff(x, 2) + x*f.diff(x) + x**2*f       # Bessel of order 0
solution = dsolve(eq, f, hint='2nd_power_series_regular', n=6)
assert solution.rhs.has(log)
```

When `factorable` is the default hint and its sub-equation reaches the
Lie group heuristics, SymPy spends about 85 s before it raises
`NotImplementedError` for `(x + cos(f))*f' + f + 2*x`; the `simplify`
calls of `abaco2_unique_general` take about 78 s in a release build of
symox (before: 210 s), most of it in the native trigonometric candidates
of the third call. Name the method (`hint='1st_exact'`) for such
equations.

### Arbitrary constants and Lie point symmetries

`constantsimp(expr, constants)` and `constant_renumber` follow SymPy,
including the `terms_gcd(clear=False, deep=True, expand=False)` factoring
of terms that share constants. Only symbols listed in `constants` are
renamed or absorbed.

`infinitesimals(eq, func, hint='default')` runs SymPy's heuristics
(`abaco1_simple`, `abaco1_product`, `abaco2_similar`,
`abaco2_unique_unknown`, `abaco2_unique_general`, `linear`,
`function_sum`, `bivariate`, `chi`); `hint='all'` combines their results.
`checkinfsol(eq, generators, func)` substitutes the generators and returns
`(is_valid, residual)` pairs. `dsolve(..., hint='lie_group')` accepts
`xi=` and `eta=` generators. Like SymPy, these are first-order searches,
not a complete classification of Lie symmetries.

## Systems

`dsolve` on a list of equations follows SymPy's `dsolve_system`: linear
first-order systems of types 1 to 6 (constant coefficients, commuting
antiderivatives, and the reducible non-constant forms), higher-order
linear systems reduced to first order, decomposition into weakly and
strongly connected components, and SymPy's nonlinear first-order forms
for two and three equations. `linear_ode_to_matrix`, `canonical_odes`,
`linodesolve` and `classify_sysode` follow `systems.py` and `ode.py`.
Systems that SymPy cannot solve raise `NotImplementedError`; then the
native fallback tries its exact spectral and triangular methods.

```python
from symox import symbols, Function, dsolve

t = symbols('t')
f, g = Function('f')(t), Function('g')(t)
solutions = dsolve([f.diff(t) - t*g, g.diff(t) - t*f], [f, g])
```

## Partial differential equations

`pdsolve`, `classify_pde`, `checkpdesol`, `pde_separate`,
`pde_separate_add` and `pde_separate_mul` follow SymPy's `pde.py`, with
the hints `1st_linear_constant_coeff_homogeneous`,
`1st_linear_constant_coeff`, `1st_linear_constant_coeff_Integral` and
`1st_linear_variable_coeff` for first-order linear equations in two
variables.

```python
from symox import symbols, Function, pdsolve, checkpdesol

x, y = symbols('x y')
u = Function('u')(x, y)
equation = 2*u.diff(x) + 3*u.diff(y) + 4*u
solution = pdsolve(equation, u)
assert checkpdesol(equation, solution, u) == (True, 0)
```

Deviations from SymPy:

- When the reference methods fail or return `False`, `pdsolve` tries the
  native characteristic solver. It covers isotropic complex directions
  (`u_x + I*u_y`), mixed fields such as `(x + y)*u_x + u_y` and variable
  zeroth-order terms such as `x*u_x - y*u_y + y**2*u - g(x, y)`.
  Otherwise the reference error is raised. The test suite checks these
  results with `checkpdesol`; `pdsolve` does not check them at runtime.
- The characteristic integration variable is renamed when it would
  capture a symbol of the equation, the characteristic coordinates of
  the constant coefficient method are fresh when the equation contains
  symbols named `xi` or `eta`, and the arbitrary function is renamed
  (`F_`) when the equation contains a function `F`. SymPy returns wrong
  or incomplete solutions in these three cases.
- `checkpdesol` combines integrals with identical limits by linearity when
  the coefficients do not depend on the bound variables, so a held
  characteristic integral can be checked.

Second-order PDE solvers, general characteristic systems and boundary
conditions are not implemented, as in SymPy. A failed separation returns
`None`; it does not claim that no separated solution exists.

## Remaining gaps

- Output forms differ from SymPy where `simplify`, `integrate` or
  `factor` choose another equivalent candidate (see the corpus figures
  above); the corpus compares strings, so this is not proved case by case.
- Known linear-hint gaps: `collect` does not match partial exponents
  (`x*exp(2*x)` on the pattern `x*exp(x)`), casus irreducibilis cubic
  roots keep radicals where SymPy uses `cos`/`sin(atan(...)/3)` forms (and
  then sometimes raises), and some variation-of-parameters integrals
  (`var_of_parameters_15`, `f'' + 6*f' + 16*f = log(x)` with
  `simplify=False`) exceed 6 GB inside integration.
- The Frobenius logarithmic solution and the PDE fallback are extensions;
  SymPy's own results in those cases are a subset of these.
- First-order PDEs beyond SymPy's hints and the native characteristic
  solver (nonlinear first-order PDEs, Charpit's method) are not
  implemented.
