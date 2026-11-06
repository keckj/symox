"""Corpus for the transcendental fallbacks of SymPy 1.14's ``_solve``.

Each case is ``(id, expression, symbols, flags)``; ``expression`` and
``symbols`` are evaluated in a namespace holding the public API of the
backend (``symox`` or ``sympy``) and the symbols below.

Families:

- ``mix_*``: radicals of the unknown mixed with other functions of it
  (``sqrt(x)*exp(x)``, ``cosh(x)/x - sqrt(2)``), where the reference takes
  its unrad fallback before ``_tsolve``;
- ``sym_*``: symbolic exponents;
- ``ts_*``: the ``solve`` calls of ``test_solve_transcendental``,
  ``test_uselogcombine``, ``test_issue_15731``, ``test_issue_17452`` and
  neighbours in SymPy's ``test_solvers.py``;
- ``lw_*``: the ``solve`` calls of ``test_lambert_multivariate`` and
  ``test_lambert_bivariate`` (including issue 23253);
- ``flag_*``: the same families with ``check=False``, ``dict``,
  ``simplify=False``, ``_unrad=False`` and ``force=False``.

The radical corpus of the unrad lot (``2026-11-06-unrad/scripts/cases.py``)
is replayed alongside (``run_cases.py --corpus unrad``) to check that it
does not regress.
"""

SYMBOLS = "x y z a b t C1 F"
POSITIVE = "p q"

CASES = [
    # radicals mixed with functions of the unknown
    ("mix_sqrt_exp", "sqrt(x)*exp(x) - 1", "x", {}),
    ("mix_sqrt_exp2", "sqrt(x)*exp(x) - 2", "x", {}),
    ("mix_sqrt_minus_exp", "sqrt(x) - exp(x)", "x", {}),
    ("mix_sqrt_plus_exp", "sqrt(x) + exp(x)", "x", {}),
    ("mix_cbrt_exp", "x**Rational(1, 3)*exp(x) - 1", "x", {}),
    ("mix_exp_sqrt_arg", "exp(sqrt(x)) - 3", "x", {}),
    ("mix_x_exp_sqrt", "x*exp(sqrt(x)) - 1", "x", {}),
    ("mix_exp_pm_sqrt", "exp(sqrt(x)) + exp(-sqrt(x)) - 3", "x", {}),
    ("mix_sqrt_log", "sqrt(x) + log(x)", "x", {}),
    ("mix_sqrt_log_prod", "sqrt(x)*log(x) - 1", "x", {}),
    ("mix_sqrt_of_exp", "sqrt(exp(x)) - 2", "x", {}),
    ("mix_sqrt_exp_plus1", "sqrt(exp(x) + 1) - 3", "x", {}),
    ("mix_sqrt_of_log", "sqrt(log(x)) - 2", "x", {}),
    ("mix_sqrt_log_two", "sqrt(log(x) + 1) - sqrt(log(x)) - 1", "x", {}),
    ("mix_sqrt_2pow", "sqrt(x)*2**x - 1", "x", {}),
    ("mix_sqrt1p_exp", "sqrt(1 + exp(x)) - exp(x)", "x", {}),
    ("mix_log_sqrt_arg", "log(sqrt(x) + 2) - 1", "x", {}),
    ("mix_cosh_over_x", "-sqrt(2) + cosh(x)/x", "x", {}),
    ("mix_sqrt_sin", "sqrt(sin(x)) - Rational(1, 2)", "x", {}),
    ("mix_sqrt_x1_exp", "sqrt(x + 1)*exp(x) - 1", "x", {}),
    ("mix_sqrt_exp_param", "sqrt(x)*exp(x) - a", "x", {}),
    ("mix_sqrt_exp_minus", "sqrt(exp(x)) + sqrt(exp(-x)) - 4", "x", {}),
    ("mix_sqrt_x_2x", "sqrt(x) - 2**x", "x", {}),
    ("mix_x_sqrt_log", "x - sqrt(log(x))", "x", {}),
    ("mix_exp_sqrt_23707", "sqrt(x - y)*exp(t*sqrt(x - y)) - exp(t*sqrt(x - y))", "y", {}),
    ("mix_sqrt_lambert", "sqrt(x)*exp(sqrt(x)) - 1", "x", {}),
    ("mix_sqrt_pos", "sqrt(p)*exp(p) - 1", "p", {}),
    ("mix_two_sqrt_exp", "sqrt(x) + sqrt(x + 1) - exp(x)", "x", {}),
    # symbolic exponents
    ("sym_xa", "x**a - 2", "x", {}),
    ("sym_x2a", "x**(2*y) - 13*x**y + 36", "x", {}),
    ("sym_xy_x2y", "x**y + x**(2*y) - 1", "x", {}),
    ("sym_xy1", "x**y - 1", "x", {}),
    ("sym_xa_exp", "x**a*exp(x) - 1", "x", {}),
    ("sym_sqrt_pow_a", "sqrt(x)**a - 2", "x", {}),
    ("sym_inv_a", "x**(1/a) - b", "x", {}),
    ("sym_x1a", "(x + 1)**a - 2", "x", {}),
    ("sym_half_a", "x**(a/2) - 3", "x", {}),
    ("sym_a_2a", "x**a + x**(2*a) - 2", "x", {}),
    ("sym_axb", "y - a*x**b", "x", {}),
    ("sym_zx", "z**x - y", "x", {}),
    ("sym_y1x", "y**(1/x) - z", "x", {}),
    ("sym_xx2", "x**x - 2", "x", {}),
    ("sym_xx_x", "x**x - x", "x", {}),
    ("sym_x3x2", "(x**3)**(x/2) + pi/2", "x", {}),
    ("sym_x_x4", "x**(x + 4) - 4", "x", {}),
    ("sym_mx", "(-x)**(-x + 4) - 4", "x", {}),
    ("sym_x2x", "x**(2/x) - 2", "x", {}),
    ("sym_x_half", "x**(x + S.Half) - 4*sqrt(2)", "x", {}),
    ("sym_x2p1", "(x**2 + 1)**x - 25", "x", {}),
    ("sym_xo2", "(x/2)**(2/x) - sqrt(2)", "x", {}),
    ("sym_7xx", "(7**x)**x + pi", "x", {}),
    ("sym_x11", "x**(x/11) + pi/11", "x", {}),
    ("sym_bx", "b**x - b**2", "x", {}),
    ("sym_neg_sqrt", "(-sqrt(sqrt(2)))**x - 2", "x", {}),
    ("sym_sqrt2x", "(sqrt(2))**x - sqrt(sqrt(2))", "x", {}),
    ("sym_x2y_param", "x**(y/2) - x**(y/3)", "x", {}),
    # test_solve_transcendental and neighbours
    ("ts_exp3", "exp(x) - 3", "x", {}),
    ("ts_prod", "(a*x + b)*(exp(x) - 3)", "x", {}),
    ("ts_cosh_y", "exp(x) + exp(-x) - y", "x", {}),
    ("ts_log3", "log(x) - 3", "x", {}),
    ("ts_sqrt3x", "sqrt(3*x) - 4", "x", {}),
    ("ts_3x2", "3**(x + 2)", "x", {}),
    ("ts_x2x", "x + 2**x", "x", {}),
    ("ts_log3x", "2*x + 5 + log(3*x - 2)", "x", {}),
    ("ts_log4x", "3*x + log(4*x)", "x", {}),
    ("ts_prod8", "(2*x + 8)*(8 + exp(x))", "x", {}),
    ("ts_exp34", "2*exp(3*x + 4) - 3", "x", {}),
    ("ts_log34", "2*log(3*x + 4) - 3", "x", {}),
    ("ts_exp1", "exp(x) + 1", "x", {}),
    ("ts_poly7", "2*(3*x + 4)**5 - 6*7**(3*x + 9)", "x", {}),
    ("ts_bexp", "y - b*exp(a/x)", "x", {}),
    ("ts_2x10", "2**x - 10", "x", {}),
    ("ts_5x2x", "exp(log(5)*x) - 2**x", "x", {}),
    ("ts_loop", "(x + 2)**y*x - 3", "x", {}),
    ("ts_sin_sqrt", "sin(sqrt(x))", "x", {}),
    ("ts_sinh_ni", "sinh(x)*sinh(sinh(x)) + cosh(x)*cosh(sinh(x)) - 3", "x", {}),
    ("ts_logcomb", "z - log(x) + log(y/(x*(-1 + y**2/x**2)))", "x", {"force": True}),
    ("ts_logcomb2", "log(x + 3) + log(1 + 3/x) - 3", "x", {}),
    ("ts_logcomb3", "log(exp(2*x) + 1) + log(-tanh(x) + 1) - log(2)", "x", {}),
    ("ts_log_sum", "log(x - 3) + log(x + 3)", "x", {}),
    ("ts_erf", "erf(x) - y", "x", {}),
    # Lambert W
    ("lw_sub1", "(x**2 - 2*x + 1).subs(x, log(x) + 3*x)", "x", {}),
    ("lw_sub2", "(x**2 - 2*x + 1).subs(x, (log(x) + 3*x)**2 - 1)", "x", {}),
    ("lw_sub3", "(x**2 - 2*x - 2).subs(x, log(x) + 3*x)", "x", {}),
    ("lw_nested", "(x*exp(x) - 3).subs(x, x*exp(x))", "x", {}),
    ("lw_sin_ni", "x - sin(x)*log(y - x)", "x", {}),
    ("lw_x3_3x", "x**3 - 3**x", "x", {}),
    ("lw_3log", "3*log(x) - x*log(3)", "x", {}),
    ("lw_W", "LambertW(2*x) - y", "x", {}),
    ("lw_biv1", "(x**2 + x)*exp(x**2 + x) - 1", "x", {}),
    ("lw_ax", "a/x + exp(x/2)", "x", {}),
    ("lw_ax_diff", "(a/x + exp(x/2)).diff(x)", "x", {}),
    ("lw_xlog", "x*log(x) + 3*x + 1", "x", {}),
    ("lw_2x_x2", "-x**2 + 2**x", "x", {}),
    ("lw_3x5", "3*x + 5 + 2**(-5*x + 3)", "x", {}),
    ("lw_5x", "5*x - 1 + 3*exp(2 - 7*x)", "x", {}),
    ("lw_logsub", "(log(x) + x).subs(x, x**2 + 1)", "x", {}),
    ("lw_cos", "3**cos(x) - cos(x)**3", "x", {}),
    ("lw_uniq", "2*log(x) - 2*log(z) + log(z + log(x) + log(z))", "x", {}),
    ("lw_ax_diff2", "(a/x + exp(x/2)).diff(x, 2)", "x", {}),
    ("lw_23253a", "1/log(sqrt(x) + 2)**2 - 1/x", "x", {}),
    ("lw_23253b", "1/log(1/sqrt(x) + 2)**2 - x", "x", {}),
    ("lw_23253c", "1/log(x**2 + 2)**2 - x**-4", "x", {}),
    ("lw_xexp", "x*exp(x) - 1", "x", {}),
    ("lw_xexp2", "x*exp(x**2) - 1", "x", {}),
    ("lw_logx2", "log(x) + x**2", "x", {}),
    # flags
    ("flag_mix_nocheck", "sqrt(x)*exp(x) - 1", "x", {"check": False}),
    ("flag_mix_dict", "sqrt(x)*exp(x) - 1", "x", {"dict": True}),
    ("flag_mix_nosimp", "sqrt(x)*exp(x) - 1", "x", {"simplify": False}),
    ("flag_mix_nounrad", "sqrt(x)*exp(x) - 1", "x", {"_unrad": False}),
    ("flag_mix_eq", "Eq(sqrt(x)*exp(x), 1)", "x", {}),
    ("flag_23253_nocheck", "1/log(sqrt(x) + 2)**2 - 1/x", "x", {"check": False}),
    ("flag_sym_noforce", "x**(2*y) - 13*x**y + 36", "x", {"force": False}),
    ("flag_lw_dict", "x*exp(x) - 1", "x", {"dict": True}),
    ("flag_sqrt_log_nocheck", "sqrt(x) + log(x)", "x", {"check": False}),
]
