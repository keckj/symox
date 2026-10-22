"""SymPy's rewrite rules for the exponential/trigonometric integrals, the
Airy, Bessel, zeta, elliptic and incomplete beta families and Catalan
numbers.

The formulas and their applicability conditions follow the
``_eval_rewrite_as_<rule>`` methods of SymPy 1.14 (BSD-3-Clause; see
LICENSE-SYMPY) in ``functions/special/error_functions.py``, ``bessel.py``,
``beta_functions.py``, ``zeta_functions.py``, ``elliptic_integrals.py`` and
``functions/combinatorial/numbers.py``. No SymPy code is imported.

Each body receives the node and its rewritten arguments; returning None
keeps the rebuilt node, as SymPy does when a method returns None.
"""


def _s():
    import symox
    return symox


def _t(*args):
    """SymPy's ``Dummy(uniquely_named_symbol('t', args).name)``."""
    s = _s()
    names = set()
    for a in args:
        for f in getattr(a, 'free_symbols', ()):
            names.add(f.name)
    name, i = 't', 0
    while name in names:
        i += 1
        name = 't%d' % i
    return s.Dummy(name)


def _E1(z):
    return _s().expint(1, z)


# --- exponential and logarithmic integrals -----------------------------------

def _Ei_expint(e, z, **hints):
    s = _s()
    return -s.expint(1, s.polar_lift(-1)*z) - s.I*s.pi


def _Ei_uppergamma(e, z, **hints):
    s = _s()
    return -s.uppergamma(0, s.polar_lift(-1)*z) - s.I*s.pi


def _Ei_li(e, z, **hints):
    s = _s()
    if getattr(z, 'class_name', None) == 'log':
        return s.li(z.args[0])
    return s.li(s.exp(z))


def _Ei_Si(e, z, **hints):
    s = _s()
    if z.is_negative:
        return s.Shi(z) + s.Chi(z) - s.I*s.pi
    return s.Shi(z) + s.Chi(z)


def _Ei_Integral(e, z, **hints):
    s = _s()
    t = _t(z)
    return s.Integral(s.E**t/t, (t, s.S.NegativeInfinity, z))


def _expint_uppergamma(e, nu, z, **hints):
    s = _s()
    return z**(nu - 1)*s.uppergamma(1 - nu, z)


def _expint_Ei(e, nu, z, **hints):
    s = _s()
    if nu == 1:
        return -s.Ei(z*s.exp_polar(-s.I*s.pi)) - s.I*s.pi
    if nu.is_Integer and nu > 1:
        # DLMF 8.19.7
        x = -s.unpolarify(z)
        k1 = int(nu) - 1
        return (x**(nu - 1)/s.factorial(nu - 1)*_expint_Ei(None, s.Integer(1), z)
                + s.exp(x)/s.factorial(nu - 1)
                * s.Add(*[s.factorial(nu - k - 2)*x**k for k in range(k1)]))
    return e


def _expint_Si(e, nu, z, **hints):
    s = _s()
    if nu != 1:
        return e
    return s.Shi(z) - s.Chi(z)


def _expint_Integral(e, nu, z, **hints):
    s = _s()
    t = _t(nu, z)
    return s.Integral(t**-nu*s.exp(-t*z), (t, 1, s.oo))


def _li_Li(e, z, **hints):
    s = _s()
    return s.Li(z) + s.li(2)


def _li_Ei(e, z, **hints):
    s = _s()
    return s.Ei(s.log(z))


def _li_uppergamma(e, z, **hints):
    s = _s()
    L = s.log(z)
    return (-s.uppergamma(0, -L) + s.S.Half*(s.log(L) - s.log(1/L)) - s.log(-L))


def _li_Si(e, z, **hints):
    s = _s()
    L = s.log(z)
    return (s.Ci(s.I*L) - s.I*s.Si(s.I*L) - s.S.Half*(s.log(1/L) - s.log(L)) - s.log(s.I*L))


def _li_Shi(e, z, **hints):
    s = _s()
    L = s.log(z)
    return s.Chi(L) - s.Shi(L) - s.S.Half*(s.log(1/L) - s.log(L))


def _li_hyper(e, z, **hints):
    s = _s()
    L = s.log(z)
    return L*s.hyper((1, 1), (2, 2), L) + s.S.Half*(s.log(L) - s.log(1/L)) + s.EulerGamma


def _li_meijerg(e, z, **hints):
    s = _s()
    L = s.log(z)
    return (-s.log(-L) - s.S.Half*(s.log(1/L) - s.log(L))
            - s.meijerg(((), (1,)), ((0, 0), ()), -L))


def _Li_li(e, z, **hints):
    s = _s()
    return s.li(z) - s.li(2)


# --- trigonometric integrals -------------------------------------------------

def _Si_expint(e, z, **hints):
    s = _s()
    return s.pi/2 + (_E1(s.polar_lift(s.I)*z) - _E1(s.polar_lift(-s.I)*z))/2/s.I


def _Ci_expint(e, z, **hints):
    s = _s()
    return -(_E1(s.polar_lift(s.I)*z) + _E1(s.polar_lift(-s.I)*z))/2


def _Shi_expint(e, z, **hints):
    s = _s()
    return (_E1(z) - _E1(s.exp_polar(s.I*s.pi)*z))/2 - s.I*s.pi/2


def _Chi_expint(e, z, **hints):
    s = _s()
    return -s.I*s.pi/2 - (_E1(z) + _E1(s.exp_polar(s.I*s.pi)*z))/2


_TRIG_EXPINT = {'Si': _Si_expint, 'Ci': _Ci_expint, 'Shi': _Shi_expint, 'Chi': _Chi_expint}


def _trig_Ei(e, z, **hints):
    s = _s()
    return _TRIG_EXPINT[e.class_name](e, z).rewrite(s.Ei)


def _trig_uppergamma(e, z, **hints):
    s = _s()
    return _TRIG_EXPINT[e.class_name](e, z).rewrite(s.uppergamma)


def _Si_Integral(e, z, **hints):
    s = _s()
    t = _t(z)
    return s.Integral(s.sinc(t), (t, 0, z))


# --- Bessel functions --------------------------------------------------------

def _besselj_besseli(e, nu, z, **hints):
    s = _s()
    return s.exp(s.I*s.pi*nu/2)*s.besseli(nu, s.polar_lift(-s.I)*z)


def _besselj_bessely(e, nu, z, **hints):
    s = _s()
    if nu.is_integer is False:
        return s.csc(s.pi*nu)*s.bessely(-nu, z) - s.cot(s.pi*nu)*s.bessely(nu, z)
    return None


def _besselj_jn(e, nu, z, **hints):
    s = _s()
    return s.sqrt(2*z/s.pi)*s.jn(nu - s.S.Half, z)


def _bessely_besselj(e, nu, z, **hints):
    s = _s()
    if nu.is_integer is False:
        return s.csc(s.pi*nu)*(s.cos(s.pi*nu)*s.besselj(nu, z) - s.besselj(-nu, z))
    return None


def _bessely_besseli(e, nu, z, **hints):
    s = _s()
    aj = _bessely_besselj(e, nu, z)
    if aj is not None:
        return aj.rewrite(s.besseli)
    return None


def _bessely_yn(e, nu, z, **hints):
    s = _s()
    return s.sqrt(2*z/s.pi)*s.yn(nu - s.S.Half, z)


def _besseli_besselj(e, nu, z, **hints):
    s = _s()
    return s.exp(-s.I*s.pi*nu/2)*s.besselj(nu, s.polar_lift(s.I)*z)


def _besseli_bessely(e, nu, z, **hints):
    s = _s()
    return _besseli_besselj(e, nu, z).rewrite(s.bessely)


def _besseli_jn(e, nu, z, **hints):
    s = _s()
    return _besseli_besselj(e, nu, z).rewrite(s.jn)


def _besselk_besseli(e, nu, z, **hints):
    s = _s()
    if nu.is_integer is False:
        return s.pi*s.csc(s.pi*nu)*(s.besseli(-nu, z) - s.besseli(nu, z))/2
    return None


def _besselk_besselj(e, nu, z, **hints):
    s = _s()
    ai = _besselk_besseli(e, nu, z)
    return ai.rewrite(s.besselj) if ai is not None else None


def _besselk_bessely(e, nu, z, **hints):
    s = _s()
    aj = _besselk_besselj(e, nu, z)
    return aj.rewrite(s.bessely) if aj is not None else None


def _besselk_yn(e, nu, z, **hints):
    s = _s()
    ay = _besselk_bessely(e, nu, z)
    return ay.rewrite(s.yn) if ay is not None else None


def _jn_besselj(e, nu, z, **hints):
    s = _s()
    return s.sqrt(s.pi/(2*z))*s.besselj(nu + s.S.Half, z)


def _jn_bessely(e, nu, z, **hints):
    s = _s()
    return s.S.NegativeOne**nu*s.sqrt(s.pi/(2*z))*s.bessely(-nu - s.S.Half, z)


def _jn_yn(e, nu, z, **hints):
    s = _s()
    return s.S.NegativeOne**nu*s.yn(-nu - 1, z)


def _yn_besselj(e, nu, z, **hints):
    s = _s()
    if not nu.is_integer:
        return None
    return s.S.NegativeOne**(nu + 1)*s.sqrt(s.pi/(2*z))*s.besselj(-nu - s.S.Half, z)


def _yn_bessely(e, nu, z, **hints):
    s = _s()
    if not nu.is_integer:
        return None
    return s.sqrt(s.pi/(2*z))*s.bessely(nu + s.S.Half, z)


def _yn_jn(e, nu, z, **hints):
    s = _s()
    return s.S.NegativeOne**(nu + 1)*s.jn(-nu - 1, z)


# --- Airy functions ----------------------------------------------------------

def _airyai_besselj(e, z, **hints):
    s = _s()
    ot, tt = s.Rational(1, 3), s.Rational(2, 3)
    a = s.Pow(-z, s.Rational(3, 2))
    if s.re(z).is_negative:
        return ot*s.sqrt(-z)*(s.besselj(-ot, tt*a) + s.besselj(ot, tt*a))
    return None


def _airyai_besseli(e, z, **hints):
    s = _s()
    ot, tt = s.Rational(1, 3), s.Rational(2, 3)
    a = s.Pow(z, s.Rational(3, 2))
    if s.re(z).is_positive:
        return ot*s.sqrt(z)*(s.besseli(-ot, tt*a) - s.besseli(ot, tt*a))
    return ot*(s.Pow(a, ot)*s.besseli(-ot, tt*a) - z*s.Pow(a, -ot)*s.besseli(ot, tt*a))


def _airyai_hyper(e, z, **hints):
    s = _s()
    R = s.Rational
    pf1 = s.S.One/(3**R(2, 3)*s.gamma(R(2, 3)))
    pf2 = z/(s.root(3, 3)*s.gamma(R(1, 3)))
    return pf1*s.hyper([], [R(2, 3)], z**3/9) - pf2*s.hyper([], [R(4, 3)], z**3/9)


def _airybi_besselj(e, z, **hints):
    s = _s()
    ot, tt = s.Rational(1, 3), s.Rational(2, 3)
    a = s.Pow(-z, s.Rational(3, 2))
    if s.re(z).is_negative:
        return s.sqrt(-z/3)*(s.besselj(-ot, tt*a) - s.besselj(ot, tt*a))
    return None


def _airybi_besseli(e, z, **hints):
    s = _s()
    ot, tt = s.Rational(1, 3), s.Rational(2, 3)
    a = s.Pow(z, s.Rational(3, 2))
    if s.re(z).is_positive:
        return s.sqrt(z)/s.sqrt(3)*(s.besseli(-ot, tt*a) + s.besseli(ot, tt*a))
    b, c = s.Pow(a, ot), s.Pow(a, -ot)
    return s.sqrt(ot)*(b*s.besseli(-ot, tt*a) + z*c*s.besseli(ot, tt*a))


def _airybi_hyper(e, z, **hints):
    s = _s()
    R = s.Rational
    pf1 = s.S.One/(s.root(3, 6)*s.gamma(R(2, 3)))
    pf2 = z*s.root(3, 6)/s.gamma(R(1, 3))
    return pf1*s.hyper([], [R(2, 3)], z**3/9) + pf2*s.hyper([], [R(4, 3)], z**3/9)


def _airyaiprime_besselj(e, z, **hints):
    s = _s()
    tt = s.Rational(2, 3)
    a = s.Pow(-z, s.Rational(3, 2))
    if s.re(z).is_negative:
        return z/3*(s.besselj(-tt, tt*a) - s.besselj(tt, tt*a))
    return None


def _airyaiprime_besseli(e, z, **hints):
    s = _s()
    ot, tt = s.Rational(1, 3), s.Rational(2, 3)
    a = tt*s.Pow(z, s.Rational(3, 2))
    if s.re(z).is_positive:
        return z/3*(s.besseli(tt, a) - s.besseli(-tt, a))
    a = s.Pow(z, s.Rational(3, 2))
    b, c = s.Pow(a, tt), s.Pow(a, -tt)
    return ot*(z**2*c*s.besseli(tt, tt*a) - b*s.besseli(-ot, tt*a))


def _airyaiprime_hyper(e, z, **hints):
    s = _s()
    R = s.Rational
    pf1 = z**2/(2*3**R(2, 3)*s.gamma(R(2, 3)))
    pf2 = 1/(s.root(3, 3)*s.gamma(R(1, 3)))
    return pf1*s.hyper([], [R(5, 3)], z**3/9) - pf2*s.hyper([], [R(1, 3)], z**3/9)


def _airybiprime_besselj(e, z, **hints):
    s = _s()
    tt = s.Rational(2, 3)
    a = tt*s.Pow(-z, s.Rational(3, 2))
    if s.re(z).is_negative:
        return -z/s.sqrt(3)*(s.besselj(-tt, a) + s.besselj(tt, a))
    return None


def _airybiprime_besseli(e, z, **hints):
    s = _s()
    ot, tt = s.Rational(1, 3), s.Rational(2, 3)
    a = tt*s.Pow(z, s.Rational(3, 2))
    if s.re(z).is_positive:
        return z/s.sqrt(3)*(s.besseli(-tt, a) + s.besseli(tt, a))
    a = s.Pow(z, s.Rational(3, 2))
    b, c = s.Pow(a, tt), s.Pow(a, -tt)
    return s.sqrt(ot)*(b*s.besseli(-tt, tt*a) + z**2*c*s.besseli(tt, tt*a))


def _airybiprime_hyper(e, z, **hints):
    s = _s()
    R = s.Rational
    pf1 = z**2/(2*s.root(3, 6)*s.gamma(R(2, 3)))
    pf2 = s.root(3, 6)/s.gamma(R(1, 3))
    return pf1*s.hyper([], [R(5, 3)], z**3/9) + pf2*s.hyper([], [R(1, 3)], z**3/9)


# --- incomplete beta ---------------------------------------------------------

def _betainc_hyper(e, a, b, x1, x2, **hints):
    s = _s()
    return (x2**a*s.hyper((a, 1 - b), (a + 1,), x2) - x1**a*s.hyper((a, 1 - b), (a + 1,), x1))/a


def _betainc_Integral(e, a, b, x1, x2, **hints):
    s = _s()
    t = _t(a, b, x1, x2)
    return s.Integral(t**(a - 1)*(1 - t)**(b - 1), (t, x1, x2))


def _betainc_regularized_hyper(e, a, b, x1, x2, **hints):
    s = _s()
    return _betainc_hyper(e, a, b, x1, x2)/s.beta(a, b)


def _betainc_regularized_Integral(e, a, b, x1, x2, **hints):
    s = _s()
    t = _t(a, b, x1, x2)
    integrand = t**(a - 1)*(1 - t)**(b - 1)
    return s.Integral(integrand, (t, x1, x2))/s.Integral(integrand, (t, 0, 1))


# --- zeta family ---------------------------------------------------------------

def _polylog_lerchphi(e, s_, z, **hints):
    s = _s()
    return z*s.lerchphi(z, s_, 1)


def _zeta_dirichlet_eta(e, s_, a=None, **hints):
    s = _s()
    if a is not None and a != 1:
        return e
    return s.dirichlet_eta(s_)/(1 - 2**(1 - s_))


def _zeta_lerchphi(e, s_, a=None, **hints):
    s = _s()
    return s.lerchphi(1, s_, 1 if a is None else a)


def _dirichlet_eta_zeta(e, s_, a=None, **hints):
    s = _s()
    if a is None or a == 1:
        return s.Piecewise((s.log(2), s.Eq(s_, 1)), ((1 - 2**(1 - s_))*s.zeta(s_), True))
    return s.Piecewise((s.log(2) - s.digamma(a) + s.digamma((a + 1)/2), s.Eq(s_, 1)),
                       (s.zeta(s_, a) - 2**(1 - s_)*s.zeta(s_, (a + 1)/2), True))


# --- elliptic integrals ------------------------------------------------------

def _elliptic_k_hyper(e, m, **hints):
    s = _s()
    return s.pi*s.S.Half*s.hyper((s.S.Half, s.S.Half), (s.S.One,), m)


def _elliptic_k_meijerg(e, m, **hints):
    s = _s()
    return s.meijerg(((s.S.Half, s.S.Half), []), ((s.S.Zero,), (s.S.Zero,)), -m)/2


def _elliptic_k_Integral(e, m, **hints):
    s = _s()
    t = _t(m)
    return s.Integral(1/s.sqrt(1 - m*s.sin(t)**2), (t, 0, s.pi/2))


def _elliptic_f_Integral(e, z, m, **hints):
    s = _s()
    t = _t(z, m)
    return s.Integral(1/(s.sqrt(1 - m*s.sin(t)**2)), (t, 0, z))


def _elliptic_e_hyper(e, *args, **hints):
    s = _s()
    if len(args) == 1:
        return (s.pi/2)*s.hyper((s.Rational(-1, 2), s.S.Half), (s.S.One,), args[0])
    return None


def _elliptic_e_meijerg(e, *args, **hints):
    s = _s()
    if len(args) == 1:
        return -s.meijerg(((s.S.Half, s.Rational(3, 2)), []), ((s.S.Zero,), (s.S.Zero,)), -args[0])/4
    return None


def _elliptic_e_Integral(e, *args, **hints):
    s = _s()
    z, m = (s.pi/2, args[0]) if len(args) == 1 else args
    t = _t(*args)
    return s.Integral(s.sqrt(1 - m*s.sin(t)**2), (t, 0, z))


def _elliptic_pi_Integral(e, *args, **hints):
    s = _s()
    if len(args) == 2:
        n, m, z = args[0], args[1], s.pi/2
    else:
        n, z, m = args
    t = _t(*args)
    return s.Integral(1/((1 - n*s.sin(t)**2)*s.sqrt(1 - m*s.sin(t)**2)), (t, 0, z))


# --- Catalan numbers ---------------------------------------------------------

def _catalan_binomial(e, n, **hints):
    s = _s()
    return s.binomial(2*n, n)/(n + 1)


def _catalan_factorial(e, n, **hints):
    s = _s()
    return s.factorial(2*n)/(s.factorial(n + 1)*s.factorial(n))


def _catalan_gamma(e, n, piecewise=True, **hints):
    s = _s()
    return 4**n*s.gamma(n + s.S.Half)/(s.gamma(s.S.Half)*s.gamma(n + 2))


def _catalan_hyper(e, n, **hints):
    s = _s()
    return s.hyper([1 - n, -n], [2], 1)


def _catalan_Product(e, n, **hints):
    s = _s()
    if not (n.is_integer and n.is_nonnegative):
        return e
    k = s.Dummy('k', integer=True, positive=True)
    return s.Product((n + k)/k, (k, 2, n))


# --- erfc and singularity functions ------------------------------------------

def _erfc_erf(e, z, **hints):
    return 1 - _s().erf(z)


def _sf_Piecewise(e, x, a, n, **hints):
    s = _s()
    if n in (-1, -2, -3, -4):
        return s.Piecewise((s.oo, s.Eq(x - a, 0)), (0, True))
    if n.is_nonnegative:
        return s.Piecewise(((x - a)**n, x - a >= 0), (0, True))
    return None


def _sf_Heaviside(e, x, a, n, **hints):
    s = _s()
    if n in (-1, -2, -3, -4):
        return s.diff(s.Heaviside(x - a), x.free_symbols.pop(), -int(n))
    if n.is_nonnegative:
        return (x - a)**n*s.Heaviside(x - a, 1)
    return None


RULES = {
    ('Ei', 'expint'): _Ei_expint,
    ('Ei', 'uppergamma'): _Ei_uppergamma,
    ('Ei', 'li'): _Ei_li,
    ('Ei', 'Si'): _Ei_Si,
    ('Ei', 'Ci'): _Ei_Si,
    ('Ei', 'Shi'): _Ei_Si,
    ('Ei', 'Chi'): _Ei_Si,
    ('Ei', 'Integral'): _Ei_Integral,
    ('expint', 'uppergamma'): _expint_uppergamma,
    ('expint', 'Ei'): _expint_Ei,
    ('expint', 'Si'): _expint_Si,
    ('expint', 'Ci'): _expint_Si,
    ('expint', 'Shi'): _expint_Si,
    ('expint', 'Chi'): _expint_Si,
    ('expint', 'Integral'): _expint_Integral,
    ('li', 'Li'): _li_Li,
    ('li', 'Ei'): _li_Ei,
    ('li', 'uppergamma'): _li_uppergamma,
    ('li', 'Si'): _li_Si,
    ('li', 'Ci'): _li_Si,
    ('li', 'Shi'): _li_Shi,
    ('li', 'Chi'): _li_Shi,
    ('li', 'hyper'): _li_hyper,
    ('li', 'meijerg'): _li_meijerg,
    ('Li', 'li'): _Li_li,
    ('Si', 'expint'): _Si_expint,
    ('Ci', 'expint'): _Ci_expint,
    ('Shi', 'expint'): _Shi_expint,
    ('Chi', 'expint'): _Chi_expint,
    ('Si', 'Ei'): _trig_Ei,
    ('Ci', 'Ei'): _trig_Ei,
    ('Shi', 'Ei'): _trig_Ei,
    ('Chi', 'Ei'): _trig_Ei,
    ('Si', 'uppergamma'): _trig_uppergamma,
    ('Ci', 'uppergamma'): _trig_uppergamma,
    ('Shi', 'uppergamma'): _trig_uppergamma,
    ('Chi', 'uppergamma'): _trig_uppergamma,
    ('Si', 'sinc'): _Si_Integral,
    ('besselj', 'besseli'): _besselj_besseli,
    ('besselj', 'bessely'): _besselj_bessely,
    ('besselj', 'jn'): _besselj_jn,
    ('bessely', 'besselj'): _bessely_besselj,
    ('bessely', 'besseli'): _bessely_besseli,
    ('bessely', 'yn'): _bessely_yn,
    ('besseli', 'besselj'): _besseli_besselj,
    ('besseli', 'bessely'): _besseli_bessely,
    ('besseli', 'jn'): _besseli_jn,
    ('besselk', 'besseli'): _besselk_besseli,
    ('besselk', 'besselj'): _besselk_besselj,
    ('besselk', 'bessely'): _besselk_bessely,
    ('besselk', 'yn'): _besselk_yn,
    ('jn', 'besselj'): _jn_besselj,
    ('jn', 'bessely'): _jn_bessely,
    ('jn', 'yn'): _jn_yn,
    ('yn', 'besselj'): _yn_besselj,
    ('yn', 'bessely'): _yn_bessely,
    ('yn', 'jn'): _yn_jn,
    ('airyai', 'besselj'): _airyai_besselj,
    ('airyai', 'besseli'): _airyai_besseli,
    ('airyai', 'hyper'): _airyai_hyper,
    ('airybi', 'besselj'): _airybi_besselj,
    ('airybi', 'besseli'): _airybi_besseli,
    ('airybi', 'hyper'): _airybi_hyper,
    ('airyaiprime', 'besselj'): _airyaiprime_besselj,
    ('airyaiprime', 'besseli'): _airyaiprime_besseli,
    ('airyaiprime', 'hyper'): _airyaiprime_hyper,
    ('airybiprime', 'besselj'): _airybiprime_besselj,
    ('airybiprime', 'besseli'): _airybiprime_besseli,
    ('airybiprime', 'hyper'): _airybiprime_hyper,
    ('betainc', 'hyper'): _betainc_hyper,
    ('betainc', 'Integral'): _betainc_Integral,
    ('betainc_regularized', 'hyper'): _betainc_regularized_hyper,
    ('betainc_regularized', 'Integral'): _betainc_regularized_Integral,
    ('polylog', 'lerchphi'): _polylog_lerchphi,
    ('zeta', 'dirichlet_eta'): _zeta_dirichlet_eta,
    ('zeta', 'lerchphi'): _zeta_lerchphi,
    ('dirichlet_eta', 'zeta'): _dirichlet_eta_zeta,
    ('elliptic_k', 'hyper'): _elliptic_k_hyper,
    ('elliptic_k', 'meijerg'): _elliptic_k_meijerg,
    ('elliptic_k', 'Integral'): _elliptic_k_Integral,
    ('elliptic_f', 'Integral'): _elliptic_f_Integral,
    ('elliptic_e', 'hyper'): _elliptic_e_hyper,
    ('elliptic_e', 'meijerg'): _elliptic_e_meijerg,
    ('elliptic_e', 'Integral'): _elliptic_e_Integral,
    ('elliptic_pi', 'Integral'): _elliptic_pi_Integral,
    ('catalan', 'binomial'): _catalan_binomial,
    ('catalan', 'factorial'): _catalan_factorial,
    ('catalan', 'gamma'): _catalan_gamma,
    ('catalan', 'hyper'): _catalan_hyper,
    ('catalan', 'Product'): _catalan_Product,
    ('erfc', 'erf'): _erfc_erf,
    ('SingularityFunction', 'Piecewise'): _sf_Piecewise,
    ('SingularityFunction', 'Heaviside'): _sf_Heaviside,
    ('SingularityFunction', 'DiracDelta'): _sf_Heaviside,
    ('SingularityFunction', 'HeavisideDiracDelta'): _sf_Heaviside,
}
