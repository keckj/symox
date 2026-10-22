"""Hill musculotendon force balance with rigid or compliant tendon."""
from abc import ABC, abstractmethod
from enum import IntEnum
import symox as S
from ..mechanics import ForceActuator
from ..vector import dynamicsymbols
from .activation import ActivationBase, _column
from .curve import (CharacteristicCurveCollection, TendonForceLengthDeGroote2016,
    TendonForceLengthInverseDeGroote2016, FiberForceLengthPassiveDeGroote2016,
    FiberForceLengthPassiveInverseDeGroote2016, FiberForceLengthActiveDeGroote2016,
    FiberForceVelocityDeGroote2016, FiberForceVelocityInverseDeGroote2016)


class MusculotendonFormulation(IntEnum):
    RIGID_TENDON = 0
    FIBER_LENGTH_EXPLICIT = 1
    TENDON_FORCE_EXPLICIT = 2
    FIBER_LENGTH_IMPLICIT = 3
    TENDON_FORCE_IMPLICIT = 4
class MusculotendonBase(ForceActuator, ABC):
    @property
    @abstractmethod
    def curves(self):raise NotImplementedError
    def __init__(self, name, pathway, activation_dynamics, *, musculotendon_dynamics=MusculotendonFormulation.RIGID_TENDON,
                 tendon_slack_length=None, peak_isometric_force=None, optimal_fiber_length=None,
                 maximal_fiber_velocity=None, optimal_pennation_angle=None, fiber_damping_coefficient=None, with_defaults=False):
        if not isinstance(activation_dynamics, ActivationBase):raise TypeError('activation_dynamics must be an activation model')
        self.name, self.activation_dynamics = name, activation_dynamics
        self.musculotendon_dynamics = MusculotendonFormulation(musculotendon_dynamics)
        supplied = (tendon_slack_length, peak_isometric_force, optimal_fiber_length, maximal_fiber_velocity, optimal_pennation_angle, fiber_damping_coefficient)
        labels = ('l_T_slack', 'F_M_max', 'l_M_opt', 'v_M_max', 'alpha_opt', 'beta')
        parameters = []
        for label, value in zip(labels, supplied):
            # SymPy names the default parameters without assumptions.
            if value is None:value = S.Symbol(label+'_'+name)
            value = S.sympify(value);setattr(self, label, value);parameters.append(value)
        self._parameters = parameters;self._with_defaults = with_defaults
        # SymPy lists the constants of the four curves in this fixed order for
        # every formulation, including the unused tendon curve of a rigid tendon.
        self._curve_constants = [] if with_defaults else [S.Symbol(f'c_{i}_{group}_{name}')
            for group, count in (('fl_T', 4), ('fl_M_pas', 2), ('fl_M_act', 12), ('fv_M', 4)) for i in range(count)]
        self._states, self._inputs, self._equations = [], [], []
        self.l_MT, self.v_MT = pathway.length, pathway.extension_velocity
        height = self.l_M_opt*S.sin(self.alpha_opt)
        form = self.musculotendon_dynamics
        if form == MusculotendonFormulation.RIGID_TENDON:
            self.l_T = self.l_T_slack
            self.l_M = S.sqrt((self.l_MT-self.l_T)**2+height**2)
            self.v_M = self.v_MT*(self.l_MT-self.l_T)/self.l_M
            # Match the rigid-tendon model's fixed pennation convention.
            self.cos_alpha = S.cos(self.alpha_opt)
        elif form in (MusculotendonFormulation.FIBER_LENGTH_EXPLICIT, MusculotendonFormulation.FIBER_LENGTH_IMPLICIT):
            lm = dynamicsymbols('l_M_tilde_'+name);self._states.append(lm)
            self.l_M = lm*self.l_M_opt
            self.cos_alpha = S.sqrt(1-(height/self.l_M)**2)
            self.l_T = self.l_MT-self.l_M*self.cos_alpha
        else:
            ft = dynamicsymbols('F_T_tilde_'+name);self._states.append(ft)
            self.l_T = self.l_T_slack*self._curve('tendon_force_length_inverse', ft)
            self.l_M = S.sqrt((self.l_MT-self.l_T)**2+height**2)
            self.cos_alpha = (self.l_MT-self.l_T)/self.l_M
        self.l_M_tilde = self.l_M/self.l_M_opt
        self.fl_M_act = self._curve('fiber_force_length_active', self.l_M_tilde)
        self.fl_M_pas = self._curve('fiber_force_length_passive', self.l_M_tilde)
        if form != MusculotendonFormulation.RIGID_TENDON:
            self.F_T_tilde = self._curve('tendon_force_length', self.l_T/self.l_T_slack) if form in (1, 3) else self._states[0]
            if form in (1, 2):
                # Like SymPy, the explicit formulations invert the force
                # balance without the fiber damping term.
                fv = (self.F_T_tilde/self.cos_alpha-self.fl_M_pas)/(self.a*self.fl_M_act)
                self.v_M = self.v_M_max*self._curve('fiber_force_velocity_inverse', fv)
            else:
                vm = dynamicsymbols('v_M_tilde_'+name);self._inputs.append(vm)
                self.v_M = vm*self.v_M_max
            if form in (1, 3):self._equations.append(self.v_M/self.l_M_opt)
            else:
                lt = S.Dummy('l_T_tilde')
                derivative = S.diff(self._curve('tendon_force_length', lt), lt).subs(lt, self.l_T/self.l_T_slack)
                self._equations.append(derivative*(self.v_MT-self.v_M/self.cos_alpha)/self.l_T_slack)
        self.v_M_tilde = self.v_M/self.v_M_max
        self.fv_M = self._curve('fiber_force_velocity', self.v_M_tilde)
        if form in (1, 2):self.F_M_tilde = self.F_T_tilde/self.cos_alpha
        else:self.F_M_tilde = self.a*self.fl_M_act*self.fv_M+self.fl_M_pas+self.beta*self.v_M_tilde
        if form == 0:self.F_T_tilde = self.F_M_tilde*self.cos_alpha
        self.force_balance = self.F_M_tilde*self.cos_alpha-self.F_T_tilde
        self.F_M, self.F_T = self.F_M_max*self.F_M_tilde, self.F_M_max*self.F_T_tilde
        super().__init__(-self.F_T, pathway)
    def _curve(self, kind, value):
        curve = getattr(self.curves, kind)
        if self._with_defaults:return curve.with_defaults(value).doit()
        group = {'tendon_force_length':'fl_T', 'tendon_force_length_inverse':'fl_T',
            'fiber_force_length_passive':'fl_M_pas', 'fiber_force_length_passive_inverse':'fl_M_pas',
            'fiber_force_length_active':'fl_M_act', 'fiber_force_velocity':'fv_M',
            'fiber_force_velocity_inverse':'fv_M'}[kind]
        constants = [S.Symbol(f'c_{i}_{group}_{self.name}') for i in range(len(curve.defaults))]
        return curve(value, *constants).doit()
    @classmethod
    def with_defaults(cls, name, pathway, activation_dynamics, *, musculotendon_dynamics=MusculotendonFormulation.RIGID_TENDON,
                      tendon_slack_length=None, peak_isometric_force=None, optimal_fiber_length=None,
                      maximal_fiber_velocity=S.Float('10.0'), optimal_pennation_angle=S.Float('0.0'), fiber_damping_coefficient=S.Float('0.1')):
        return cls(name, pathway, activation_dynamics, musculotendon_dynamics=musculotendon_dynamics,
            tendon_slack_length=tendon_slack_length, peak_isometric_force=peak_isometric_force,
            optimal_fiber_length=optimal_fiber_length, maximal_fiber_velocity=maximal_fiber_velocity,
            optimal_pennation_angle=optimal_pennation_angle, fiber_damping_coefficient=fiber_damping_coefficient, with_defaults=True)
    @property
    def a(self):return self.activation_dynamics.a
    @property
    def e(self):return self.activation_dynamics.e
    @property
    def x(self):return _column(self._states+list(self.activation_dynamics.x))
    @property
    def r(self):return _column(self._inputs+list(self.activation_dynamics.r))
    @property
    def p(self):return _column([c for c in self._parameters if not c.is_number]+list(self.activation_dynamics.p)+self._curve_constants)
    state_vars = property(lambda self:self.x)
    input_vars = property(lambda self:self.r)
    constants = property(lambda self:self.p)
    M = property(lambda self:S.eye(self.x.rows))
    F = property(lambda self:self.rhs())
    def rhs(self):return _column(self._equations+list(self.activation_dynamics.rhs()))
    @property
    def algebraic_equations(self):return _column([self.force_balance]) if self.musculotendon_dynamics in (3, 4) else S.zeros(0, 1)


class MusculotendonDeGroote2016(MusculotendonBase):
    curves = CharacteristicCurveCollection(TendonForceLengthDeGroote2016, TendonForceLengthInverseDeGroote2016,
        FiberForceLengthPassiveDeGroote2016, FiberForceLengthPassiveInverseDeGroote2016,
        FiberForceLengthActiveDeGroote2016, FiberForceVelocityDeGroote2016, FiberForceVelocityInverseDeGroote2016)
