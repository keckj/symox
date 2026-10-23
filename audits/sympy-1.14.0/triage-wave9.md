# Triage wave 9 (audit9): remaining c and f objects

Replay of the 16 objects (examples and keyword modes give the same result).
Before: 16 observation_difference. After: 5 examples_match, 11
observation_difference. The remaining differences are documented below:
a = cosmetic/representation, b = equivalent form, c = missing feature,
h = SymPy is wrong and symox keeps the correct result, i = environment.

| object | before | after | note |
|---|---|---|---|
| AlgebraicField | c | match | Degree 5 and 6 Galois groups ported (resolvent lookup, Tschirnhausen). |
| assoc_laguerre | c | h | Alpha derivative implemented as Sum(L_k^a(x)/(n - k)). SymPy divides by n - a, which is wrong (0.319 off at n=2, a=1/3, x=2/7). |
| harmonic | c | b + c | harmonic(oo, m) evaluates for symbolic m. [25] equivalent form of the polygamma rewrite (b). [34] limit(harmonic(n, m+1), n, oo) stays unevaluated (c, needs Hurwitz zeta limits in gruntz). |
| julia_code | c | match | User functions accept matrix arguments. |
| octave_code | c | match | Same fix as julia_code. |
| shape | c | match | Expr subclasses keep matrix arguments. |
| lambdify | c | a + h + c + i | Missing-module ImportError fixed [55]; [62, 65] tensorflow-only (i). [17, 34, 36, 38, 98] numpy scalar repr (a). [72] SymPy's ZeroDivisionError comes from evaluating every Piecewise branch; symox returns the selected branch (h). [75] 0.0 vs 0 (a). [92] inspect.getsource on CompiledFunction (c). |
| DiscreteMarkovChain | c | c | Str states are not expression atoms in symox. |
| MultivariateNormal | c | c | MatrixSymbol parameters (symbolic dimension) not supported. |
| Sum.euler_maclaurin | h | h | Already documented. |
| approximants | h | h | Already documented. |
| combine_kronecker | h | h | SymPy drops the cross terms; 1x1 proof in tests (12 vs 25). |
| ImplicitRegion.regular_point | h | h | SymPy returns a singular point. Already documented. |
| nsolve | f | f (reduced) | [25] fixed (mpmath-style secant/bisect, step limits, messages). [17] real part about 1e-61 instead of 0; [20] residual of a diverging run differs; [27] residual printed as 10000.0000000000000000 instead of 10000.0. |
| aesara_function | f | match | Output identical to SymPy. |
| warns_deprecated_sympy | f | a + i | Warning class and warns ported. [2] doctest stacklevel (i) now passes; [3] module name symox vs sympy in the message (a). |
