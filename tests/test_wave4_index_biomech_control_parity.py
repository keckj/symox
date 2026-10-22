"""get_contraction_structure nesting, musculotendon constants and exact gbt, as in SymPy 1.14."""
import pytest

import symox as S
from symox.tensor import IndexedBase, Idx
from symox.tensor.index_methods import get_contraction_structure, get_indices


def _norm(value):
    # Order-free text form: keys and set members as strings, lists sorted.
    if isinstance(value, dict):
        keys = []
        for key, entry in value.items():
            text = str(tuple(sorted(map(str, key)))) if isinstance(key, tuple) else str(key)
            keys.append((text, _norm(entry)))
        return '{' + ', '.join('%s: %s' % item for item in sorted(keys)) + '}'
    if isinstance(value, list):
        return '[' + ', '.join(sorted(_norm(entry) for entry in value)) + ']'
    if isinstance(value, set):
        return '{' + ', '.join(sorted(map(str, value))) + '}'
    return str(value)


x, y, A, B, C = map(IndexedBase, 'xyABC')
i, j, k, l = map(Idx, 'ijkl')


def test_contraction_structure_of_products_of_contracted_factors():
    expr = x[i, i]*y[j, j]
    d = get_contraction_structure(expr)
    assert set(map(str, d)) == {'None', str(expr)}
    assert d[None] == {expr}
    assert _norm(d[expr]) == "[{('i',): {x[i, i]}}, {('j',): {y[j, j]}}]"
    assert _norm(get_contraction_structure(A[j, j]**A[i, i])) == (
        "{A[j, j]**A[i, i]: [{('i',): {A[i, i]}}, {('j',): {A[j, j]}}], None: {A[j, j]**A[i, i]}}")


@pytest.mark.parametrize('build, expected', [
    (lambda: x[i]*y[i], "{('i',): {x[i]*y[i]}}"),
    (lambda: x[i]*y[i] + A[j, j], "{('i',): {x[i]*y[i]}, ('j',): {A[j, j]}}"),
    (lambda: A[i, j]*x[j] + B[i, j]*y[j], "{('j',): {A[i, j]*x[j], B[i, j]*y[j]}}"),
    (lambda: A[i, j]*B[j, k]*C[k, l], "{('j', 'k'): {A[i, j]*B[j, k]*C[k, l]}}"),
    (lambda: x[i] + y[i], "{None: {x[i], y[i]}}"),
    (lambda: S.sin(x[i]*y[i]), "{None: {sin(x[i]*y[i])}, sin(x[i]*y[i]): [{('i',): {x[i]*y[i]}}]}"),
    (lambda: S.exp(x[i]*y[i]), "{None: {exp(x[i]*y[i])}, exp(x[i]*y[i]): [{('i',): {x[i]*y[i]}}]}"),
    (lambda: (x[i]*y[i])**2, "{('i',): {x[i]**2*y[i]**2}}"),
    (lambda: A[i, i], "{('i',): {A[i, i]}}"),
])
def test_contraction_structure_matches_reference(build, expected):
    assert _norm(get_contraction_structure(build())) == expected


def test_contraction_structure_nested_function_and_sum():
    f = S.Function('f')
    expr = f(x[j, k, k])*y[j]
    d = get_contraction_structure(expr)
    assert d[(j,)] == {expr}
    assert _norm(d[expr]) == "[{None: {f(x[j, k, k])}, f(x[j, k, k]): [{('k',): {x[j, k, k]}}]}]"
    expr = x[i]*y[i]*(A[j, j] + B[k, k]*C[l, l])
    d = get_contraction_structure(expr)
    assert d[(i,)] == {expr}
    assert _norm(d[expr]) == ("[{('j',): {A[j, j]}, B[k, k]*C[l, l]: [{('k',): {B[k, k]}}, "
                              "{('l',): {C[l, l]}}], None: {B[k, k]*C[l, l]}}]")


def test_get_indices_matches_reference():
    assert get_indices(A[i, j]*x[j]) == ({i}, {})
    assert get_indices(A[i, j]*B[j, k]*C[k, l]) == ({i, l}, {})
    assert get_indices(x[i, i, i]) == (set(), {})
    assert get_indices(x[i]*y[i]*A[i]) == (set(), {})
    assert get_indices(S.exp(x[i]*y[j])) == ({i, j}, {})
    with pytest.raises(ValueError):
        get_indices(x[i] + y[j])


def _muscle(**kwargs):
    from symox.physics.mechanics import LinearPathway, Point, ReferenceFrame, dynamicsymbols
    from symox.physics.biomechanics import FirstOrderActivationDeGroote2016, MusculotendonDeGroote2016
    N = ReferenceFrame('N')
    O, I = S.symbols('O, I', cls=Point)
    q, u = dynamicsymbols('q u')
    I.set_pos(O, q*N.x)
    O.set_vel(N, 0)
    I.set_vel(N, u*N.x)
    activation = FirstOrderActivationDeGroote2016('muscle')
    defaults = kwargs.pop('defaults', False)
    maker = MusculotendonDeGroote2016.with_defaults if defaults else MusculotendonDeGroote2016
    return maker('muscle', LinearPathway(O, I), activation, **kwargs), q


_P = (['l_T_slack', 'F_M_max', 'l_M_opt', 'v_M_max', 'alpha_opt', 'beta', 'tau_a', 'tau_d', 'b']
      + ['c_%d_fl_T' % n for n in range(4)] + ['c_%d_fl_M_pas' % n for n in range(2)]
      + ['c_%d_fl_M_act' % n for n in range(12)] + ['c_%d_fv_M' % n for n in range(4)])


@pytest.mark.parametrize('form', [0, 1, 2])
def test_musculotendon_constants_in_reference_order(form):
    muscle, _ = _muscle(musculotendon_dynamics=form)
    assert [str(c) for c in muscle.p] == [name + '_muscle' for name in _P]
    assert muscle.p[0] == S.Symbol('l_T_slack_muscle')
    assert muscle.p[6] == S.Symbol('tau_a_muscle')
    assert muscle.x.rows == (1 if form == 0 else 2)
    assert muscle.r.rows == 1


@pytest.mark.parametrize('form, force', [(0, -0.1350588818), (1, 0.08187813221), (2, -0.342)])
def test_musculotendon_with_defaults_classmethod(form, force):
    from symox.physics.mechanics import dynamicsymbols
    muscle, q = _muscle(musculotendon_dynamics=form, defaults=True)
    assert [str(c) for c in muscle.p] == ['l_T_slack_muscle', 'F_M_max_muscle', 'l_M_opt_muscle',
                                         'tau_a_muscle', 'tau_d_muscle', 'b_muscle']
    assert muscle.v_M_max == S.Float('10.0') and muscle.beta == S.Float('0.1')
    values = {c: 0.37 + 0.01*n for n, c in enumerate(muscle.p)}
    values.update({q: 0.7, muscle.a: 0.4, muscle.e: 0.6})
    for state in muscle.x:
        values.setdefault(state, 0.9)
    value = S.sympify(muscle.force).doit().subs(q.diff(dynamicsymbols._t), 0.2).subs(values)
    assert float(S.N(value)) == pytest.approx(force, rel=1e-8)


def test_explicit_musculotendon_with_symbolic_damping():
    # SymPy's explicit formulations neglect fiber damping instead of raising.
    muscle, _ = _muscle(musculotendon_dynamics=1)
    assert muscle.rhs().rows == 2


def test_gbt_keeps_exact_alpha():
    from symox.physics.control import TransferFunction, gbt, bilinear, forward_diff, backward_diff
    L, R, T, s = S.symbols('L R T s')
    tf = TransferFunction(1, s*L + R, s)
    num, den = gbt(tf, T, 0.5)
    half = S.Mul(S.Rational(1, 2), T, S.Pow(L + R*T/2, -1))
    assert num == [half, half]
    assert den == [1, S.Mul(-L + R*T/2, S.Pow(L + R*T/2, -1))]
    assert not any(v.has(S.Float) for v in num + den)
    assert gbt(tf, T, 0.5) == bilinear(tf, T)
    assert forward_diff(tf, T) == ([T/L], [1, (-L + R*T)/L])
    assert backward_diff(tf, T) == ([T/(L + R*T), 0], [1, -L/(L + R*T)])
    num, den = gbt(tf, T, 0.3)
    assert not any(v.has(S.Float) for v in num + den)
    assert bilinear(TransferFunction(s + 2, s**2 + 3*s + 5, s), S.Rational(1, 10)) == (
        [S.Rational(22, 465), S.Rational(4, 465), S.Rational(-6, 155)],
        [1, S.Rational(-158, 93), S.Rational(23, 31)])
