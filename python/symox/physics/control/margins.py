"""Gain and phase margins from positive real crossover frequencies."""
import symox as S
from .lti import _tf


def _crossovers(expression,omega):
    numerator=S.cancel(expression).as_numer_denom()[0]
    if numerator==0:return []
    roots=S.solve(numerator,omega)
    return sorted([r for r in roots if r.is_real is not False and r.is_positive is True],key=lambda r:float(r.evalf()))


def phase_margin(system):
    system=_tf(system);omega=S.Symbol('omega',positive=True);response=system.eval_frequency(S.I*omega)
    crossovers=_crossovers(S.simplify(response*S.conjugate(response)-1),omega)
    if not crossovers:return S.oo
    return S.simplify(180+180*S.arg(response.subs(omega,crossovers[0]))/S.pi)


def gain_margin(system):
    system=_tf(system);omega=S.Symbol('omega',positive=True);response=system.eval_frequency(S.I*omega)
    crossovers=_crossovers(S.im(response),omega)
    crossovers=[w for w in crossovers if S.re(response.subs(omega,w)).is_negative is True]
    if not crossovers:return S.oo
    return S.simplify(-20*S.log(S.Abs(response.subs(omega,crossovers[0])),10))


__all__=['phase_margin','gain_margin']
