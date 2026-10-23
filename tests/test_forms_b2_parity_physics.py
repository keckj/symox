"""Equivalent-form parity for the physics records of audit batch b2.

Expected strings come from SymPy 1.14.
"""
import pytest

from symox import N, Rational, Symbol, factor, oo, symbols
from symox.physics.hydrogen import E_nl_dirac
from symox.physics.sho import R_nl
from symox.physics.wigner import wigner_d_small
from symox.physics.optics import Medium, TWave, conjugate_gauss_beams, geometric_conj_ab, gaussian_conj
from symox.physics.biomechanics import (FiberForceLengthPassiveDeGroote2016, FirstOrderActivationDeGroote2016, MusculotendonDeGroote2016, TendonForceLengthDeGroote2016, TendonForceLengthInverseDeGroote2016)
from symox.physics.mechanics import dynamicsymbols, ReferenceFrame, Point, LinearPathway
from symox.physics.units import farad, meter, henry

r, l, nu, beta = symbols('r l nu beta')
a, b, f, wi, wo, lam = symbols('a b f w_i w_o l')
A1, phi1, A2, phi2 = symbols('A1 phi1 A2 phi2')
eps, mu = symbols('epsilon mu')
cases = {
 'dirac_30': lambda: E_nl_dirac(3, 0),
 'dirac_21_down': lambda: E_nl_dirac(2, 1, False),
 'dirac_sym': lambda: E_nl_dirac(2, 0, Z=3),
 'rnl_0l': lambda: R_nl(0, l, 1, r),
 'rnl_10': lambda: R_nl(1, 0, 1, r),
 'rnl_00nu': lambda: R_nl(0, 0, nu, r),
 'wigner_1': lambda: wigner_d_small(1, beta),
 'wigner_half': lambda: wigner_d_small(Rational(1, 2), beta),
 'conj_ab': lambda: geometric_conj_ab(a, b),
 'conj_ab_oo': lambda: (geometric_conj_ab(oo, b), geometric_conj_ab(a, oo)),
 'cgb_0': lambda: conjugate_gauss_beams(lam, wi, wo, f=f)[0],
 'cgb_1': lambda: factor(conjugate_gauss_beams(lam, wi, wo, f=f)[1]),
 'gconj_sout': lambda: gaussian_conj(3, 2, 5)[0],
 'medium_imp': lambda: Medium('m1', permittivity=10*farad/meter, permeability=1*henry/meter).intrinsic_impedance,
 'medium_n': lambda: Medium('m2', eps, mu).refractive_index,
 'medium_repr': lambda: Medium('m', n=2),
 'twave_add': lambda: TWave(A1, f, phi1) + TWave(A2, f, phi2),
 'twave_amp': lambda: (TWave(A1, f, phi1) + TWave(A2, f, phi2)).amplitude,
 'twave_neg': lambda: -TWave(A1, f, phi1),
}
def mt():
    q = dynamicsymbols('q'); N = ReferenceFrame('N'); O, P = Point('pO'), Point('pI')
    P.set_pos(O, q*N.x); pw = LinearPathway(O, P)
    a_ = FirstOrderActivationDeGroote2016.with_defaults('muscle')
    m = MusculotendonDeGroote2016.with_defaults('muscle', pw, a_, tendon_slack_length=Symbol('l_T_slack'), peak_isometric_force=Symbol('F_M_max'), optimal_fiber_length=Symbol('l_M_opt'), maximal_fiber_velocity=Symbol('v_M_max'), optimal_pennation_angle=Symbol('alpha_opt'), fiber_damping_coefficient=Symbol('beta'))
    return (m.x, m.r, m.p, m.constants)
cases['musculotendon'] = mt
lt = Symbol('l_T_tilde'); lm = Symbol('l_M_tilde'); c = symbols('c0:12')
cases['fl_T_diff'] = lambda: TendonForceLengthDeGroote2016(lt, *c[:4]).diff(lt)
cases['fl_M_pas_diff'] = lambda: FiberForceLengthPassiveDeGroote2016(lm, *c[:2]).diff(lm)
cases['fl_M_pas_doit'] = lambda: FiberForceLengthPassiveDeGroote2016.with_defaults(lm).doit()
cases['fl_T_inv'] = lambda: TendonForceLengthInverseDeGroote2016.with_defaults(Symbol('fl_T')).doit()

EXPECTED = {
    'dirac_30': '-0.0555562951740285',
    'dirac_21_down': '-0.125002080189006',
    'dirac_sym': '-1.12516853302077',
    'rnl_0l': 'r**l*sqrt(2**(l + 3/2)*2**(l + 2)/factorial2(2*l + 1))*exp(-r**2)/pi**(1/4)',
    'rnl_10': '4*2**(1/4)*sqrt(3)*(3/2 - 2*r**2)*exp(-r**2)/(3*pi**(1/4))',
    'rnl_00nu': '2*2**(3/4)*sqrt(nu**(3/2))*exp(-nu*r**2)/pi**(1/4)',
    'wigner_1': 'Matrix([[cos(beta/2)**2, sqrt(2)*sin(beta/2)*cos(beta/2), sin(beta/2)**2], [-sqrt(2)*sin(beta/2)*cos(beta/2), -sin(beta/2)**2 + cos(beta/2)**2, sqrt(2)*sin(beta/2)*cos(beta/2)], [sin(beta/2)**2, -sqrt(2)*sin(beta/2)*cos(beta/2), cos(beta/2)**2]])',
    'wigner_half': 'Matrix([[cos(beta/2), sin(beta/2)], [-sin(beta/2), cos(beta/2)]])',
    'conj_ab': 'a*b/(a + b)',
    'conj_ab_oo': '(b, a)',
    'cgb_0': 'f*(1 - sqrt(w_i**2/w_o**2 - pi**2*w_i**4/(f**2*l**2)))',
    'cgb_1': 'f*w_o**2*(w_i**2/w_o**2 - sqrt(w_i**2/w_o**2 - pi**2*w_i**4/(f**2*l**2)))/w_i**2',
    'gconj_sout': '-5/4',
    'medium_imp': 'sqrt(10)*sqrt(henry)/(10*sqrt(farad))',
    'medium_n': '299792458*meter*sqrt(epsilon*mu)/second',
    'medium_repr': 'MediumN(m, 2)',
    'twave_add': 'TWave(sqrt(A1**2 + 2*A1*A2*cos(phi1 - phi2) + A2**2), f, atan2(A1*sin(phi1) + A2*sin(phi2), A1*cos(phi1) + A2*cos(phi2)), 1/f, n)',
    'twave_amp': 'sqrt(A1**2 + 2*A1*A2*cos(phi1 - phi2) + A2**2)',
    'twave_neg': 'TWave(-A1, f, phi1, 1/f, n)',
    'musculotendon': '(Matrix([[a_muscle(t)]]), Matrix([[e_muscle(t)]]), Matrix([\n[l_T_slack],\n[  F_M_max],\n[  l_M_opt],\n[  v_M_max],\n[alpha_opt],\n[     beta]]), Matrix([\n[l_T_slack],\n[  F_M_max],\n[  l_M_opt],\n[  v_M_max],\n[alpha_opt],\n[     beta]]))',
    'fl_T_diff': 'c0*c3*exp(c3*(-c1 + l_T_tilde))',
    'fl_M_pas_diff': 'c1*exp(c1*(l_M_tilde - 1)/c0)/(c0*(exp(c1) - 1))',
    'fl_M_pas_doit': '2.37439874427164e-5*exp(6.66666666666667*l_M_tilde) - 0.0186573603637741',
    'fl_T_inv': '0.02946663003430684*log(5.0*fl_T + 1.25) + 0.995',
}


@pytest.mark.parametrize('name', sorted(EXPECTED))
def test_matches_sympy_form(name):
    assert str(cases[name]()) == EXPECTED[name]


def test_conjugate_gauss_beams_arguments():
    with pytest.raises(ValueError):
        conjugate_gauss_beams(1, 2, 3)
    with pytest.raises(ValueError):
        conjugate_gauss_beams(1, 2, 3, f=1, s_in=2)
    with pytest.raises(NotImplementedError):
        conjugate_gauss_beams(1, 2, 3, dist=4)
    # Symox extension: s_in recovers the focal length.
    s_in = conjugate_gauss_beams(1, 2, 3, f=10)[0]
    assert abs(N(conjugate_gauss_beams(1, 2, 3, s_in=s_in)[2]) - 10) < 1e-12


@pytest.mark.parametrize('form', [3, 4])
def test_implicit_musculotendon_formulations_are_kept(form):
    q = dynamicsymbols('q'); N_ = ReferenceFrame('N'); O, P = Point('pO'), Point('pI')
    P.set_pos(O, q*N_.x)
    act = FirstOrderActivationDeGroote2016.with_defaults('muscle')
    m = MusculotendonDeGroote2016('muscle', LinearPathway(O, P), act, musculotendon_dynamics=form, tendon_slack_length=1, peak_isometric_force=1, optimal_fiber_length=1, maximal_fiber_velocity=10, optimal_pennation_angle=0, fiber_damping_coefficient=Rational(1, 10))
    assert m.algebraic_equations.shape == (1, 1)
    assert m.r.rows == 2
