"""Gain and phase margins of SISO systems (reference algorithm over positive frequencies)."""
import symox as S
from .lti import TransferFunction


def _frequency_expression(system):
    siso = system.doit() if hasattr(system,'doit') else system
    if not isinstance(siso,TransferFunction):raise ValueError("Margins are only applicable for SISO LTI systems.")
    _w = S.Dummy("w",real=True)
    expr = system.to_expr()
    if expr.has(S.exp):raise NotImplementedError("Margins for systems with Time delay terms are not supported.")
    if len(expr.free_symbols) > 1:
        raise ValueError("Extra degree of freedom found. Make sure"
            " that there are no free symbols in the dynamical system other"
            " than the variable of Laplace transform.")
    return _w,expr.subs({siso.var:S.I*_w})


def _decided(value,positive):
    """Decide a candidate root numerically at 50 digits, refusing undecidable ties."""
    imag,real = S.im(value).evalf(50),S.re(value).evalf(50)
    if abs(imag) > S.Float('1e-40'):return False
    if S.simplify(S.im(value)) != 0:raise NotImplementedError('cannot decide whether %s is real' % value)
    if abs(real) <= S.Float('1e-40'):
        if S.simplify(S.re(value)) == 0:return False
        raise NotImplementedError('cannot decide the sign of %s' % value)
    return real > 0 if positive else True


def _members(solutions):
    """Elements of a solveset result over (0, oo) whose interval membership was left unevaluated."""
    if solutions == S.EmptySet:return []
    if isinstance(solutions,S.Interval):raise TypeError("'Interval' object is not iterable")
    if isinstance(solutions,S.FiniteSet):return list(solutions)
    if isinstance(solutions,S.Union):return [v for a in solutions.args for v in _members(a)]
    if isinstance(solutions,S.Intersection):
        finite = [a for a in solutions.args if isinstance(a,S.FiniteSet)]
        others = [a for a in solutions.args if not isinstance(a,S.FiniteSet)]
        if len(finite) == 1 and others == [S.Interval(0,S.oo,left_open=True)]:
            return [v for v in finite[0] if _decided(v,True)]
    raise NotImplementedError('cannot enumerate the positive solutions %s' % solutions)


def _positive_roots(equation,_w,polynomial):
    domain = S.Interval(0,S.oo,left_open=True)
    solutions = S.solveset(equation,_w,domain)
    if isinstance(solutions,S.ConditionSet):
        # 20*log10|G(i*w)| = 0 exactly when |G(i*w)|**2 = G(i*w)*conjugate(G(i*w)) = 1 for real w.
        # After cancel the numerator is coprime with the denominator, so its roots keep G finite.
        solutions = S.solveset(polynomial(),_w,domain)
    if isinstance(solutions,S.ConditionSet):
        raise NotImplementedError('cannot solve %s = 0 over positive frequencies' % equation)
    return list(S.FiniteSet(*_members(solutions)))


def phase_margin(system):
    _w,w_expr = _frequency_expression(system)
    mag = 20*S.log(S.Abs(w_expr),10)
    mag_sol = _positive_roots(mag,_w,lambda:S.numer(S.cancel(S.expand(w_expr*S.conjugate(w_expr))-1)))
    if len(mag_sol) == 0:pm = S.Integer(-180)
    else:pm = ((S.arg(w_expr)*S.Integer(180)/S.pi).subs({_w:mag_sol[0]})+S.Integer(180)) % 360
    if _at_least_180(pm):pm = pm-360
    return pm


def _at_least_180(pm):
    # Reference: ``pm >= 180`` decided by SymPy's numeric sign evaluation of
    # an algebraic constant.  The native relational does not evaluate it, so
    # decide it from a high-precision value and refuse when it is ambiguous.
    relation = pm >= 180
    if relation in (S.true,S.false) or isinstance(relation,bool):
        return bool(relation)
    difference = (pm-180).evalf(50)
    if difference.is_Number and difference.is_real and abs(difference) > S.Float('1e-40'):
        return bool(difference > 0)
    raise TypeError('cannot determine truth value of Relational')


def gain_margin(system):
    _w,w_expr = _frequency_expression(system)
    mag = 20*S.log(S.Abs(w_expr),10)
    phase_sol = _positive_roots(S.numer(w_expr.as_real_imag()[1].cancel()),_w,None)
    if len(phase_sol) == 0:return S.oo
    return -mag.subs({_w:phase_sol[0]})


__all__=['phase_margin','gain_margin']
