"""Musculotendon models, ported from SymPy 1.14."""
from abc import abstractmethod
from enum import IntEnum, unique
from symox.core.numbers import Float, Integer
from symox.core.symbol import Symbol, symbols
from symox.functions.elementary.miscellaneous import sqrt
from symox.functions.elementary.trigonometric import cos, sin
from symox.matrices.dense import MutableDenseMatrix as Matrix, diag, eye, zeros
from symox.physics.biomechanics.activation import ActivationBase
from symox.physics.biomechanics.curve import CharacteristicCurveCollection, FiberForceLengthActiveDeGroote2016, FiberForceLengthPassiveDeGroote2016, FiberForceLengthPassiveInverseDeGroote2016, FiberForceVelocityDeGroote2016, FiberForceVelocityInverseDeGroote2016, TendonForceLengthDeGroote2016, TendonForceLengthInverseDeGroote2016
from symox.physics.biomechanics._mixin import _NamedMixin
from symox.physics.mechanics.actuator import ForceActuator
from symox.physics.vector.functions import dynamicsymbols
__all__ = ['MusculotendonBase', 'MusculotendonDeGroote2016', 'MusculotendonFormulation']

@unique
class MusculotendonFormulation(IntEnum):
    RIGID_TENDON = 0
    FIBER_LENGTH_EXPLICIT = 1
    TENDON_FORCE_EXPLICIT = 2
    FIBER_LENGTH_IMPLICIT = 3
    TENDON_FORCE_IMPLICIT = 4

    def __str__(self):
        return str(self.value)
_DEFAULT_MUSCULOTENDON_FORMULATION = MusculotendonFormulation.RIGID_TENDON

class MusculotendonBase(ForceActuator, _NamedMixin):

    def __init__(self, name, pathway, activation_dynamics, *, musculotendon_dynamics=_DEFAULT_MUSCULOTENDON_FORMULATION, tendon_slack_length=None, peak_isometric_force=None, optimal_fiber_length=None, maximal_fiber_velocity=None, optimal_pennation_angle=None, fiber_damping_coefficient=None, with_defaults=False):
        self.name = name
        super().__init__(Symbol('F'), pathway)
        if not isinstance(activation_dynamics, ActivationBase):
            msg = f"Can't set attribute `activation_dynamics` to {activation_dynamics} as it must be of type `ActivationBase`, not {type(activation_dynamics)}."
            raise TypeError(msg)
        self._activation_dynamics = activation_dynamics
        self._child_objects = (self._activation_dynamics,)
        if tendon_slack_length is not None:
            self._l_T_slack = tendon_slack_length
        else:
            self._l_T_slack = Symbol(f'l_T_slack_{self.name}')
        if peak_isometric_force is not None:
            self._F_M_max = peak_isometric_force
        else:
            self._F_M_max = Symbol(f'F_M_max_{self.name}')
        if optimal_fiber_length is not None:
            self._l_M_opt = optimal_fiber_length
        else:
            self._l_M_opt = Symbol(f'l_M_opt_{self.name}')
        if maximal_fiber_velocity is not None:
            self._v_M_max = maximal_fiber_velocity
        else:
            self._v_M_max = Symbol(f'v_M_max_{self.name}')
        if optimal_pennation_angle is not None:
            self._alpha_opt = optimal_pennation_angle
        else:
            self._alpha_opt = Symbol(f'alpha_opt_{self.name}')
        if fiber_damping_coefficient is not None:
            self._beta = fiber_damping_coefficient
        else:
            self._beta = Symbol(f'beta_{self.name}')
        self._with_defaults = with_defaults
        if musculotendon_dynamics == MusculotendonFormulation.RIGID_TENDON:
            self._rigid_tendon_musculotendon_dynamics()
        elif musculotendon_dynamics == MusculotendonFormulation.FIBER_LENGTH_EXPLICIT:
            self._fiber_length_explicit_musculotendon_dynamics()
        elif musculotendon_dynamics == MusculotendonFormulation.TENDON_FORCE_EXPLICIT:
            self._tendon_force_explicit_musculotendon_dynamics()
        elif musculotendon_dynamics == MusculotendonFormulation.FIBER_LENGTH_IMPLICIT:
            self._fiber_length_implicit_musculotendon_dynamics()
        elif musculotendon_dynamics == MusculotendonFormulation.TENDON_FORCE_IMPLICIT:
            self._tendon_force_implicit_musculotendon_dynamics()
        else:
            msg = f'Musculotendon dynamics {repr(musculotendon_dynamics)} passed to `musculotendon_dynamics` was of type {type(musculotendon_dynamics)}, must be {MusculotendonFormulation}.'
            raise TypeError(msg)
        self._musculotendon_dynamics = musculotendon_dynamics
        self._force = -self._F_T

    @classmethod
    def with_defaults(cls, name, pathway, activation_dynamics, *, musculotendon_dynamics=_DEFAULT_MUSCULOTENDON_FORMULATION, tendon_slack_length=None, peak_isometric_force=None, optimal_fiber_length=None, maximal_fiber_velocity=Float('10.0'), optimal_pennation_angle=Float('0.0'), fiber_damping_coefficient=Float('0.1')):
        return cls(name, pathway, activation_dynamics=activation_dynamics, musculotendon_dynamics=musculotendon_dynamics, tendon_slack_length=tendon_slack_length, peak_isometric_force=peak_isometric_force, optimal_fiber_length=optimal_fiber_length, maximal_fiber_velocity=maximal_fiber_velocity, optimal_pennation_angle=optimal_pennation_angle, fiber_damping_coefficient=fiber_damping_coefficient, with_defaults=True)

    @abstractmethod
    def curves(cls):
        pass

    @property
    def tendon_slack_length(self):
        return self._l_T_slack

    @property
    def l_T_slack(self):
        return self._l_T_slack

    @property
    def peak_isometric_force(self):
        return self._F_M_max

    @property
    def F_M_max(self):
        return self._F_M_max

    @property
    def optimal_fiber_length(self):
        return self._l_M_opt

    @property
    def l_M_opt(self):
        return self._l_M_opt

    @property
    def maximal_fiber_velocity(self):
        return self._v_M_max

    @property
    def v_M_max(self):
        return self._v_M_max

    @property
    def optimal_pennation_angle(self):
        return self._alpha_opt

    @property
    def alpha_opt(self):
        return self._alpha_opt

    @property
    def fiber_damping_coefficient(self):
        return self._beta

    @property
    def beta(self):
        return self._beta

    @property
    def activation_dynamics(self):
        return self._activation_dynamics

    @property
    def excitation(self):
        return self._activation_dynamics._e

    @property
    def e(self):
        return self._activation_dynamics._e

    @property
    def activation(self):
        return self._activation_dynamics._a

    @property
    def a(self):
        return self._activation_dynamics._a

    @property
    def musculotendon_dynamics(self):
        return self._musculotendon_dynamics

    def _rigid_tendon_musculotendon_dynamics(self):
        self._l_MT = self.pathway.length
        self._v_MT = self.pathway.extension_velocity
        self._l_T = self._l_T_slack
        self._l_T_tilde = Integer(1)
        self._l_M = sqrt((self._l_MT - self._l_T) ** 2 + (self._l_M_opt * sin(self._alpha_opt)) ** 2)
        self._l_M_tilde = self._l_M / self._l_M_opt
        self._v_M = self._v_MT * (self._l_MT - self._l_T_slack) / self._l_M
        self._v_M_tilde = self._v_M / self._v_M_max
        if self._with_defaults:
            self._fl_T = self.curves.tendon_force_length.with_defaults(self._l_T_tilde)
            self._fl_M_pas = self.curves.fiber_force_length_passive.with_defaults(self._l_M_tilde)
            self._fl_M_act = self.curves.fiber_force_length_active.with_defaults(self._l_M_tilde)
            self._fv_M = self.curves.fiber_force_velocity.with_defaults(self._v_M_tilde)
        else:
            fl_T_constants = symbols(f'c_0:4_fl_T_{self.name}')
            self._fl_T = self.curves.tendon_force_length(self._l_T_tilde, *fl_T_constants)
            fl_M_pas_constants = symbols(f'c_0:2_fl_M_pas_{self.name}')
            self._fl_M_pas = self.curves.fiber_force_length_passive(self._l_M_tilde, *fl_M_pas_constants)
            fl_M_act_constants = symbols(f'c_0:12_fl_M_act_{self.name}')
            self._fl_M_act = self.curves.fiber_force_length_active(self._l_M_tilde, *fl_M_act_constants)
            fv_M_constants = symbols(f'c_0:4_fv_M_{self.name}')
            self._fv_M = self.curves.fiber_force_velocity(self._v_M_tilde, *fv_M_constants)
        self._F_M_tilde = self.a * self._fl_M_act * self._fv_M + self._fl_M_pas + self._beta * self._v_M_tilde
        self._F_T_tilde = self._F_M_tilde
        self._F_M = self._F_M_tilde * self._F_M_max
        self._cos_alpha = cos(self._alpha_opt)
        self._F_T = self._F_M * self._cos_alpha
        self._state_vars = zeros(0, 1)
        self._input_vars = zeros(0, 1)
        self._state_eqns = zeros(0, 1)
        self._curve_constants = Matrix(fl_T_constants + fl_M_pas_constants + fl_M_act_constants + fv_M_constants) if not self._with_defaults else zeros(0, 1)

    def _fiber_length_explicit_musculotendon_dynamics(self):
        self._l_M_tilde = dynamicsymbols(f'l_M_tilde_{self.name}')
        self._l_MT = self.pathway.length
        self._v_MT = self.pathway.extension_velocity
        self._l_M = self._l_M_tilde * self._l_M_opt
        self._l_T = self._l_MT - sqrt(self._l_M ** 2 - (self._l_M_opt * sin(self._alpha_opt)) ** 2)
        self._l_T_tilde = self._l_T / self._l_T_slack
        self._cos_alpha = (self._l_MT - self._l_T) / self._l_M
        if self._with_defaults:
            self._fl_T = self.curves.tendon_force_length.with_defaults(self._l_T_tilde)
            self._fl_M_pas = self.curves.fiber_force_length_passive.with_defaults(self._l_M_tilde)
            self._fl_M_act = self.curves.fiber_force_length_active.with_defaults(self._l_M_tilde)
        else:
            fl_T_constants = symbols(f'c_0:4_fl_T_{self.name}')
            self._fl_T = self.curves.tendon_force_length(self._l_T_tilde, *fl_T_constants)
            fl_M_pas_constants = symbols(f'c_0:2_fl_M_pas_{self.name}')
            self._fl_M_pas = self.curves.fiber_force_length_passive(self._l_M_tilde, *fl_M_pas_constants)
            fl_M_act_constants = symbols(f'c_0:12_fl_M_act_{self.name}')
            self._fl_M_act = self.curves.fiber_force_length_active(self._l_M_tilde, *fl_M_act_constants)
        self._F_T_tilde = self._fl_T
        self._F_T = self._F_T_tilde * self._F_M_max
        self._F_M = self._F_T / self._cos_alpha
        self._F_M_tilde = self._F_M / self._F_M_max
        self._fv_M = (self._F_M_tilde - self._fl_M_pas) / (self.a * self._fl_M_act)
        if self._with_defaults:
            self._v_M_tilde = self.curves.fiber_force_velocity_inverse.with_defaults(self._fv_M)
        else:
            fv_M_constants = symbols(f'c_0:4_fv_M_{self.name}')
            self._v_M_tilde = self.curves.fiber_force_velocity_inverse(self._fv_M, *fv_M_constants)
        self._dl_M_tilde_dt = self._v_M_max / self._l_M_opt * self._v_M_tilde
        self._state_vars = Matrix([self._l_M_tilde])
        self._input_vars = zeros(0, 1)
        self._state_eqns = Matrix([self._dl_M_tilde_dt])
        self._curve_constants = Matrix(fl_T_constants + fl_M_pas_constants + fl_M_act_constants + fv_M_constants) if not self._with_defaults else zeros(0, 1)

    def _tendon_force_explicit_musculotendon_dynamics(self):
        self._F_T_tilde = dynamicsymbols(f'F_T_tilde_{self.name}')
        self._l_MT = self.pathway.length
        self._v_MT = self.pathway.extension_velocity
        self._fl_T = self._F_T_tilde
        if self._with_defaults:
            self._fl_T_inv = self.curves.tendon_force_length_inverse.with_defaults(self._fl_T)
        else:
            fl_T_constants = symbols(f'c_0:4_fl_T_{self.name}')
            self._fl_T_inv = self.curves.tendon_force_length_inverse(self._fl_T, *fl_T_constants)
        self._l_T_tilde = self._fl_T_inv
        self._l_T = self._l_T_tilde * self._l_T_slack
        self._l_M = sqrt((self._l_MT - self._l_T) ** 2 + (self._l_M_opt * sin(self._alpha_opt)) ** 2)
        self._l_M_tilde = self._l_M / self._l_M_opt
        if self._with_defaults:
            self._fl_M_pas = self.curves.fiber_force_length_passive.with_defaults(self._l_M_tilde)
            self._fl_M_act = self.curves.fiber_force_length_active.with_defaults(self._l_M_tilde)
        else:
            fl_M_pas_constants = symbols(f'c_0:2_fl_M_pas_{self.name}')
            self._fl_M_pas = self.curves.fiber_force_length_passive(self._l_M_tilde, *fl_M_pas_constants)
            fl_M_act_constants = symbols(f'c_0:12_fl_M_act_{self.name}')
            self._fl_M_act = self.curves.fiber_force_length_active(self._l_M_tilde, *fl_M_act_constants)
        self._cos_alpha = (self._l_MT - self._l_T) / self._l_M
        self._F_T = self._F_T_tilde * self._F_M_max
        self._F_M = self._F_T / self._cos_alpha
        self._F_M_tilde = self._F_M / self._F_M_max
        self._fv_M = (self._F_M_tilde - self._fl_M_pas) / (self.a * self._fl_M_act)
        if self._with_defaults:
            self._fv_M_inv = self.curves.fiber_force_velocity_inverse.with_defaults(self._fv_M)
        else:
            fv_M_constants = symbols(f'c_0:4_fv_M_{self.name}')
            self._fv_M_inv = self.curves.fiber_force_velocity_inverse(self._fv_M, *fv_M_constants)
        self._v_M_tilde = self._fv_M_inv
        self._v_M = self._v_M_tilde * self._v_M_max
        self._v_T = self._v_MT - self._v_M / self._cos_alpha
        self._v_T_tilde = self._v_T / self._l_T_slack
        if self._with_defaults:
            self._fl_T = self.curves.tendon_force_length.with_defaults(self._l_T_tilde)
        else:
            self._fl_T = self.curves.tendon_force_length(self._l_T_tilde, *fl_T_constants)
        self._dF_T_tilde_dt = self._fl_T.diff(dynamicsymbols._t).subs({self._l_T_tilde.diff(dynamicsymbols._t): self._v_T_tilde})
        self._state_vars = Matrix([self._F_T_tilde])
        self._input_vars = zeros(0, 1)
        self._state_eqns = Matrix([self._dF_T_tilde_dt])
        self._curve_constants = Matrix(fl_T_constants + fl_M_pas_constants + fl_M_act_constants + fv_M_constants) if not self._with_defaults else zeros(0, 1)

    # SymPy 1.14 raises NotImplementedError for the two implicit
    # formulations. Symox keeps them (an extension): the normalized fiber
    # velocity becomes an input and the force balance an algebraic equation.
    def _implicit_fiber_velocity(self):
        self._v_M_tilde = dynamicsymbols(f'v_M_tilde_{self.name}')
        self._v_M = self._v_M_tilde * self._v_M_max
        if self._with_defaults:
            self._fv_M = self.curves.fiber_force_velocity.with_defaults(self._v_M_tilde)
            fv_M_constants = ()
        else:
            fv_M_constants = symbols(f'c_0:4_fv_M_{self.name}')
            self._fv_M = self.curves.fiber_force_velocity(self._v_M_tilde, *fv_M_constants)
        self._F_M_tilde = self.a * self._fl_M_act * self._fv_M + self._fl_M_pas + self._beta * self._v_M_tilde
        self._F_M = self._F_M_tilde * self._F_M_max
        self._force_balance = self._F_M_tilde * self._cos_alpha - self._F_T_tilde
        self._input_vars = Matrix([self._v_M_tilde])
        return fv_M_constants

    def _fiber_length_implicit_musculotendon_dynamics(self):
        self._l_M_tilde = dynamicsymbols(f'l_M_tilde_{self.name}')
        self._l_MT = self.pathway.length
        self._v_MT = self.pathway.extension_velocity
        self._l_M = self._l_M_tilde * self._l_M_opt
        self._l_T = self._l_MT - sqrt(self._l_M ** 2 - (self._l_M_opt * sin(self._alpha_opt)) ** 2)
        self._l_T_tilde = self._l_T / self._l_T_slack
        self._cos_alpha = (self._l_MT - self._l_T) / self._l_M
        fl_T_constants = fl_M_pas_constants = fl_M_act_constants = ()
        if self._with_defaults:
            self._fl_T = self.curves.tendon_force_length.with_defaults(self._l_T_tilde)
            self._fl_M_pas = self.curves.fiber_force_length_passive.with_defaults(self._l_M_tilde)
            self._fl_M_act = self.curves.fiber_force_length_active.with_defaults(self._l_M_tilde)
        else:
            fl_T_constants = symbols(f'c_0:4_fl_T_{self.name}')
            self._fl_T = self.curves.tendon_force_length(self._l_T_tilde, *fl_T_constants)
            fl_M_pas_constants = symbols(f'c_0:2_fl_M_pas_{self.name}')
            self._fl_M_pas = self.curves.fiber_force_length_passive(self._l_M_tilde, *fl_M_pas_constants)
            fl_M_act_constants = symbols(f'c_0:12_fl_M_act_{self.name}')
            self._fl_M_act = self.curves.fiber_force_length_active(self._l_M_tilde, *fl_M_act_constants)
        self._F_T_tilde = self._fl_T
        self._F_T = self._F_T_tilde * self._F_M_max
        fv_M_constants = self._implicit_fiber_velocity()
        self._dl_M_tilde_dt = self._v_M_max / self._l_M_opt * self._v_M_tilde
        self._state_vars = Matrix([self._l_M_tilde])
        self._state_eqns = Matrix([self._dl_M_tilde_dt])
        self._curve_constants = Matrix(fl_T_constants + fl_M_pas_constants + fl_M_act_constants + fv_M_constants) if not self._with_defaults else zeros(0, 1)

    def _tendon_force_implicit_musculotendon_dynamics(self):
        self._F_T_tilde = dynamicsymbols(f'F_T_tilde_{self.name}')
        self._l_MT = self.pathway.length
        self._v_MT = self.pathway.extension_velocity
        fl_T_constants = fl_M_pas_constants = fl_M_act_constants = ()
        if self._with_defaults:
            self._fl_T_inv = self.curves.tendon_force_length_inverse.with_defaults(self._F_T_tilde)
        else:
            fl_T_constants = symbols(f'c_0:4_fl_T_{self.name}')
            self._fl_T_inv = self.curves.tendon_force_length_inverse(self._F_T_tilde, *fl_T_constants)
        self._l_T_tilde = self._fl_T_inv
        self._l_T = self._l_T_tilde * self._l_T_slack
        self._l_M = sqrt((self._l_MT - self._l_T) ** 2 + (self._l_M_opt * sin(self._alpha_opt)) ** 2)
        self._l_M_tilde = self._l_M / self._l_M_opt
        if self._with_defaults:
            self._fl_M_pas = self.curves.fiber_force_length_passive.with_defaults(self._l_M_tilde)
            self._fl_M_act = self.curves.fiber_force_length_active.with_defaults(self._l_M_tilde)
        else:
            fl_M_pas_constants = symbols(f'c_0:2_fl_M_pas_{self.name}')
            self._fl_M_pas = self.curves.fiber_force_length_passive(self._l_M_tilde, *fl_M_pas_constants)
            fl_M_act_constants = symbols(f'c_0:12_fl_M_act_{self.name}')
            self._fl_M_act = self.curves.fiber_force_length_active(self._l_M_tilde, *fl_M_act_constants)
        self._cos_alpha = (self._l_MT - self._l_T) / self._l_M
        self._F_T = self._F_T_tilde * self._F_M_max
        fv_M_constants = self._implicit_fiber_velocity()
        self._v_T = self._v_MT - self._v_M / self._cos_alpha
        self._v_T_tilde = self._v_T / self._l_T_slack
        if self._with_defaults:
            self._fl_T = self.curves.tendon_force_length.with_defaults(self._l_T_tilde)
        else:
            self._fl_T = self.curves.tendon_force_length(self._l_T_tilde, *fl_T_constants)
        self._dF_T_tilde_dt = self._fl_T.diff(dynamicsymbols._t).subs({self._l_T_tilde.diff(dynamicsymbols._t): self._v_T_tilde})
        self._state_vars = Matrix([self._F_T_tilde])
        self._state_eqns = Matrix([self._dF_T_tilde_dt])
        self._curve_constants = Matrix(fl_T_constants + fl_M_pas_constants + fl_M_act_constants + fv_M_constants) if not self._with_defaults else zeros(0, 1)

    @property
    def algebraic_equations(self):
        """Symox extension: the force balance of the implicit formulations."""
        if self._musculotendon_dynamics in (MusculotendonFormulation.FIBER_LENGTH_IMPLICIT, MusculotendonFormulation.TENDON_FORCE_IMPLICIT):
            return Matrix([self._force_balance])
        return zeros(0, 1)

    @property
    def state_vars(self):
        state_vars = [self._state_vars]
        for child in self._child_objects:
            state_vars.append(child.state_vars)
        return Matrix.vstack(*state_vars)

    @property
    def x(self):
        state_vars = [self._state_vars]
        for child in self._child_objects:
            state_vars.append(child.state_vars)
        return Matrix.vstack(*state_vars)

    @property
    def input_vars(self):
        input_vars = [self._input_vars]
        for child in self._child_objects:
            input_vars.append(child.input_vars)
        return Matrix.vstack(*input_vars)

    @property
    def r(self):
        input_vars = [self._input_vars]
        for child in self._child_objects:
            input_vars.append(child.input_vars)
        return Matrix.vstack(*input_vars)

    @property
    def constants(self):
        musculotendon_constants = [self._l_T_slack, self._F_M_max, self._l_M_opt, self._v_M_max, self._alpha_opt, self._beta]
        musculotendon_constants = [c for c in musculotendon_constants if not c.is_number]
        constants = [Matrix(musculotendon_constants) if musculotendon_constants else zeros(0, 1)]
        for child in self._child_objects:
            constants.append(child.constants)
        constants.append(self._curve_constants)
        return Matrix.vstack(*constants)

    @property
    def p(self):
        musculotendon_constants = [self._l_T_slack, self._F_M_max, self._l_M_opt, self._v_M_max, self._alpha_opt, self._beta]
        musculotendon_constants = [c for c in musculotendon_constants if not c.is_number]
        constants = [Matrix(musculotendon_constants) if musculotendon_constants else zeros(0, 1)]
        for child in self._child_objects:
            constants.append(child.constants)
        constants.append(self._curve_constants)
        return Matrix.vstack(*constants)

    @property
    def M(self):
        M = [eye(len(self._state_vars))]
        for child in self._child_objects:
            M.append(child.M)
        return diag(*M)

    @property
    def F(self):
        F = [self._state_eqns]
        for child in self._child_objects:
            F.append(child.F)
        return Matrix.vstack(*F)

    def rhs(self):
        is_explicit = (MusculotendonFormulation.FIBER_LENGTH_EXPLICIT, MusculotendonFormulation.TENDON_FORCE_EXPLICIT)
        if self.musculotendon_dynamics is MusculotendonFormulation.RIGID_TENDON:
            child_rhs = [child.rhs() for child in self._child_objects]
            return Matrix.vstack(*child_rhs)
        elif self.musculotendon_dynamics in is_explicit:
            rhs = self._state_eqns
            child_rhs = [child.rhs() for child in self._child_objects]
            return Matrix.vstack(rhs, *child_rhs)
        return self.M.solve(self.F)

    def __repr__(self):
        return f'{self.__class__.__name__}({self.name!r}, pathway={self.pathway!r}, activation_dynamics={self.activation_dynamics!r}, musculotendon_dynamics={self.musculotendon_dynamics}, tendon_slack_length={self._l_T_slack!r}, peak_isometric_force={self._F_M_max!r}, optimal_fiber_length={self._l_M_opt!r}, maximal_fiber_velocity={self._v_M_max!r}, optimal_pennation_angle={self._alpha_opt!r}, fiber_damping_coefficient={self._beta!r})'

    def __str__(self):
        return str(self.force)

class MusculotendonDeGroote2016(MusculotendonBase):
    curves = CharacteristicCurveCollection(tendon_force_length=TendonForceLengthDeGroote2016, tendon_force_length_inverse=TendonForceLengthInverseDeGroote2016, fiber_force_length_passive=FiberForceLengthPassiveDeGroote2016, fiber_force_length_passive_inverse=FiberForceLengthPassiveInverseDeGroote2016, fiber_force_length_active=FiberForceLengthActiveDeGroote2016, fiber_force_velocity=FiberForceVelocityDeGroote2016, fiber_force_velocity_inverse=FiberForceVelocityInverseDeGroote2016)
