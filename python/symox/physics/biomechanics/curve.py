"""Characteristic musculotendon curves (De Groote et al., 2016), ported from SymPy 1.14."""
from dataclasses import dataclass
from symox.core.expr import UnevaluatedExpr
from symox.core.function import ArgumentIndexError, Function
from symox.core.numbers import Float, Integer
from symox.functions.elementary.exponential import exp, log
from symox.functions.elementary.hyperbolic import cosh, sinh
from symox.functions.elementary.miscellaneous import sqrt
from symox.printing.precedence import PRECEDENCE
__all__ = ['CharacteristicCurveCollection', 'CharacteristicCurveFunction', 'FiberForceLengthActiveDeGroote2016', 'FiberForceLengthPassiveDeGroote2016', 'FiberForceLengthPassiveInverseDeGroote2016', 'FiberForceVelocityDeGroote2016', 'FiberForceVelocityInverseDeGroote2016', 'TendonForceLengthDeGroote2016', 'TendonForceLengthInverseDeGroote2016']

class CharacteristicCurveFunction(Function):
    _symox_user_hooks = True

    @classmethod
    def eval(cls):
        msg = f'Cannot directly instantiate {cls.__name__!r}, instances of characteristic curves must be of a concrete subclass.'
        raise TypeError(msg)

    def _print_code(self, printer):
        return printer._print(printer.parenthesize(self.doit(deep=False, evaluate=False), PRECEDENCE['Atom']))
    _ccode = _print_code
    _cupycode = _print_code
    _cxxcode = _print_code
    _fcode = _print_code
    _jaxcode = _print_code
    _lambdacode = _print_code
    _mpmathcode = _print_code
    _octave = _print_code
    _pythoncode = _print_code
    _numpycode = _print_code
    _scipycode = _print_code



def _evalf_at(expr, prec):
    """SymPy's ``expr._eval_evalf(prec)`` on the evaluated curve formula."""
    from mpmath.libmp.libmpf import prec_to_dps
    return expr.evalf(prec_to_dps(prec))

def _doit_argument(expr, deep, hints):
    """``expr.doit(deep=deep, **hints)`` as in SymPy, where ``Basic.doit``
    forwards hints it does not use. The native ``doit`` only takes its own
    hints, so the curves are evaluated with the hints and the rest without."""
    try:
        return expr.doit(deep=deep, **hints)
    except NotImplementedError:
        pass
    from symox.core.traversal import preorder_traversal
    curves = {e: e.doit(deep=deep, **hints) for e in preorder_traversal(expr)
              if isinstance(e, CharacteristicCurveFunction)}
    return expr.xreplace(curves).doit(deep=deep) if curves else expr.doit(deep=deep)

class TendonForceLengthDeGroote2016(CharacteristicCurveFunction):
    _symox_user_hooks = True

    @classmethod
    def with_defaults(cls, l_T_tilde):
        c0 = Float('0.2')
        c1 = Float('0.995')
        c2 = Float('0.25')
        c3 = Float('33.93669377311689')
        return cls(l_T_tilde, c0, c1, c2, c3)

    @classmethod
    def eval(cls, l_T_tilde, c0, c1, c2, c3):
        pass

    def _eval_evalf(self, prec):
        return _evalf_at(self.doit(deep=False, evaluate=False), prec)

    def doit(self, deep=True, evaluate=True, **hints):
        (l_T_tilde, *constants) = self.args
        if deep:
            hints['evaluate'] = evaluate
            l_T_tilde = _doit_argument(l_T_tilde, deep, hints)
            (c0, c1, c2, c3) = [_doit_argument(c, deep, hints) for c in constants]
        else:
            (c0, c1, c2, c3) = constants
        if evaluate:
            return c0 * exp(c3 * (l_T_tilde - c1)) - c2
        return c0 * exp(c3 * UnevaluatedExpr(l_T_tilde - c1)) - c2

    def fdiff(self, argindex=1):
        (l_T_tilde, c0, c1, c2, c3) = self.args
        if argindex == 1:
            return c0 * c3 * exp(c3 * UnevaluatedExpr(l_T_tilde - c1))
        elif argindex == 2:
            return exp(c3 * UnevaluatedExpr(l_T_tilde - c1))
        elif argindex == 3:
            return -c0 * c3 * exp(c3 * UnevaluatedExpr(l_T_tilde - c1))
        elif argindex == 4:
            return Integer(-1)
        elif argindex == 5:
            return c0 * (l_T_tilde - c1) * exp(c3 * UnevaluatedExpr(l_T_tilde - c1))
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return TendonForceLengthInverseDeGroote2016

    def _latex(self, printer):
        l_T_tilde = self.args[0]
        _l_T_tilde = printer._print(l_T_tilde)
        return '\\operatorname{fl}^T \\left( %s \\right)' % _l_T_tilde

class TendonForceLengthInverseDeGroote2016(CharacteristicCurveFunction):
    _symox_user_hooks = True

    @classmethod
    def with_defaults(cls, fl_T):
        c0 = Float('0.2')
        c1 = Float('0.995')
        c2 = Float('0.25')
        c3 = Float('33.93669377311689')
        return cls(fl_T, c0, c1, c2, c3)

    @classmethod
    def eval(cls, fl_T, c0, c1, c2, c3):
        pass

    def _eval_evalf(self, prec):
        return _evalf_at(self.doit(deep=False, evaluate=False), prec)

    def doit(self, deep=True, evaluate=True, **hints):
        (fl_T, *constants) = self.args
        if deep:
            hints['evaluate'] = evaluate
            fl_T = _doit_argument(fl_T, deep, hints)
            (c0, c1, c2, c3) = [_doit_argument(c, deep, hints) for c in constants]
        else:
            (c0, c1, c2, c3) = constants
        if evaluate:
            return log((fl_T + c2) / c0) / c3 + c1
        return log(UnevaluatedExpr((fl_T + c2) / c0)) / c3 + c1

    def fdiff(self, argindex=1):
        (fl_T, c0, c1, c2, c3) = self.args
        if argindex == 1:
            return 1 / (c3 * (fl_T + c2))
        elif argindex == 2:
            return -1 / (c0 * c3)
        elif argindex == 3:
            return Integer(1)
        elif argindex == 4:
            return 1 / (c3 * (fl_T + c2))
        elif argindex == 5:
            return -log(UnevaluatedExpr((fl_T + c2) / c0)) / c3 ** 2
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return TendonForceLengthDeGroote2016

    def _latex(self, printer):
        fl_T = self.args[0]
        _fl_T = printer._print(fl_T)
        return '\\left( \\operatorname{fl}^T \\right)^{-1} \\left( %s \\right)' % _fl_T

class FiberForceLengthPassiveDeGroote2016(CharacteristicCurveFunction):
    _symox_user_hooks = True

    @classmethod
    def with_defaults(cls, l_M_tilde):
        c0 = Float('0.6')
        c1 = Float('4.0')
        return cls(l_M_tilde, c0, c1)

    @classmethod
    def eval(cls, l_M_tilde, c0, c1):
        pass

    def _eval_evalf(self, prec):
        return _evalf_at(self.doit(deep=False, evaluate=False), prec)

    def doit(self, deep=True, evaluate=True, **hints):
        (l_M_tilde, *constants) = self.args
        if deep:
            hints['evaluate'] = evaluate
            l_M_tilde = _doit_argument(l_M_tilde, deep, hints)
            (c0, c1) = [_doit_argument(c, deep, hints) for c in constants]
        else:
            (c0, c1) = constants
        if evaluate:
            return (exp(c1 * (l_M_tilde - 1) / c0) - 1) / (exp(c1) - 1)
        return (exp(c1 * UnevaluatedExpr(l_M_tilde - 1) / c0) - 1) / (exp(c1) - 1)

    def fdiff(self, argindex=1):
        (l_M_tilde, c0, c1) = self.args
        if argindex == 1:
            return c1 * exp(c1 * UnevaluatedExpr(l_M_tilde - 1) / c0) / (c0 * (exp(c1) - 1))
        elif argindex == 2:
            return -c1 * exp(c1 * UnevaluatedExpr(l_M_tilde - 1) / c0) * UnevaluatedExpr(l_M_tilde - 1) / (c0 ** 2 * (exp(c1) - 1))
        elif argindex == 3:
            return -exp(c1) * (-1 + exp(c1 * UnevaluatedExpr(l_M_tilde - 1) / c0)) / (exp(c1) - 1) ** 2 + exp(c1 * UnevaluatedExpr(l_M_tilde - 1) / c0) * (l_M_tilde - 1) / (c0 * (exp(c1) - 1))
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return FiberForceLengthPassiveInverseDeGroote2016

    def _latex(self, printer):
        l_M_tilde = self.args[0]
        _l_M_tilde = printer._print(l_M_tilde)
        return '\\operatorname{fl}^M_{pas} \\left( %s \\right)' % _l_M_tilde

class FiberForceLengthPassiveInverseDeGroote2016(CharacteristicCurveFunction):
    _symox_user_hooks = True

    @classmethod
    def with_defaults(cls, fl_M_pas):
        c0 = Float('0.6')
        c1 = Float('4.0')
        return cls(fl_M_pas, c0, c1)

    @classmethod
    def eval(cls, fl_M_pas, c0, c1):
        pass

    def _eval_evalf(self, prec):
        return _evalf_at(self.doit(deep=False, evaluate=False), prec)

    def doit(self, deep=True, evaluate=True, **hints):
        (fl_M_pas, *constants) = self.args
        if deep:
            hints['evaluate'] = evaluate
            fl_M_pas = _doit_argument(fl_M_pas, deep, hints)
            (c0, c1) = [_doit_argument(c, deep, hints) for c in constants]
        else:
            (c0, c1) = constants
        if evaluate:
            return c0 * log(fl_M_pas * (exp(c1) - 1) + 1) / c1 + 1
        return c0 * log(UnevaluatedExpr(fl_M_pas * (exp(c1) - 1)) + 1) / c1 + 1

    def fdiff(self, argindex=1):
        (fl_M_pas, c0, c1) = self.args
        if argindex == 1:
            return c0 * (exp(c1) - 1) / (c1 * (fl_M_pas * (exp(c1) - 1) + 1))
        elif argindex == 2:
            return log(fl_M_pas * (exp(c1) - 1) + 1) / c1
        elif argindex == 3:
            return c0 * fl_M_pas * exp(c1) / (c1 * (fl_M_pas * (exp(c1) - 1) + 1)) - c0 * log(fl_M_pas * (exp(c1) - 1) + 1) / c1 ** 2
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return FiberForceLengthPassiveDeGroote2016

    def _latex(self, printer):
        fl_M_pas = self.args[0]
        _fl_M_pas = printer._print(fl_M_pas)
        return '\\left( \\operatorname{fl}^M_{pas} \\right)^{-1} \\left( %s \\right)' % _fl_M_pas

class FiberForceLengthActiveDeGroote2016(CharacteristicCurveFunction):
    _symox_user_hooks = True

    @classmethod
    def with_defaults(cls, l_M_tilde):
        c0 = Float('0.814')
        c1 = Float('1.06')
        c2 = Float('0.162')
        c3 = Float('0.0633')
        c4 = Float('0.433')
        c5 = Float('0.717')
        c6 = Float('-0.0299')
        c7 = Float('0.2')
        c8 = Float('0.1')
        c9 = Float('1.0')
        c10 = Float('0.354')
        c11 = Float('0.0')
        return cls(l_M_tilde, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11)

    @classmethod
    def eval(cls, l_M_tilde, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11):
        pass

    def _eval_evalf(self, prec):
        return _evalf_at(self.doit(deep=False, evaluate=False), prec)

    def doit(self, deep=True, evaluate=True, **hints):
        (l_M_tilde, *constants) = self.args
        if deep:
            hints['evaluate'] = evaluate
            l_M_tilde = _doit_argument(l_M_tilde, deep, hints)
            constants = [_doit_argument(c, deep, hints) for c in constants]
        (c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) = constants
        if evaluate:
            return c0 * exp(-((l_M_tilde - c1) / (c2 + c3 * l_M_tilde)) ** 2 / 2) + c4 * exp(-((l_M_tilde - c5) / (c6 + c7 * l_M_tilde)) ** 2 / 2) + c8 * exp(-((l_M_tilde - c9) / (c10 + c11 * l_M_tilde)) ** 2 / 2)
        return c0 * exp(-(UnevaluatedExpr(l_M_tilde - c1) / (c2 + c3 * l_M_tilde)) ** 2 / 2) + c4 * exp(-(UnevaluatedExpr(l_M_tilde - c5) / (c6 + c7 * l_M_tilde)) ** 2 / 2) + c8 * exp(-(UnevaluatedExpr(l_M_tilde - c9) / (c10 + c11 * l_M_tilde)) ** 2 / 2)

    def fdiff(self, argindex=1):
        (l_M_tilde, c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) = self.args
        if argindex == 1:
            return c0 * (c3 * (l_M_tilde - c1) ** 2 / (c2 + c3 * l_M_tilde) ** 3 + (c1 - l_M_tilde) / (c2 + c3 * l_M_tilde) ** 2) * exp(-(l_M_tilde - c1) ** 2 / (2 * (c2 + c3 * l_M_tilde) ** 2)) + c4 * (c7 * (l_M_tilde - c5) ** 2 / (c6 + c7 * l_M_tilde) ** 3 + (c5 - l_M_tilde) / (c6 + c7 * l_M_tilde) ** 2) * exp(-(l_M_tilde - c5) ** 2 / (2 * (c6 + c7 * l_M_tilde) ** 2)) + c8 * (c11 * (l_M_tilde - c9) ** 2 / (c10 + c11 * l_M_tilde) ** 3 + (c9 - l_M_tilde) / (c10 + c11 * l_M_tilde) ** 2) * exp(-(l_M_tilde - c9) ** 2 / (2 * (c10 + c11 * l_M_tilde) ** 2))
        elif argindex == 2:
            return exp(-(l_M_tilde - c1) ** 2 / (2 * (c2 + c3 * l_M_tilde) ** 2))
        elif argindex == 3:
            return c0 * (l_M_tilde - c1) / (c2 + c3 * l_M_tilde) ** 2 * exp(-(l_M_tilde - c1) ** 2 / (2 * (c2 + c3 * l_M_tilde) ** 2))
        elif argindex == 4:
            return c0 * (l_M_tilde - c1) ** 2 / (c2 + c3 * l_M_tilde) ** 3 * exp(-(l_M_tilde - c1) ** 2 / (2 * (c2 + c3 * l_M_tilde) ** 2))
        elif argindex == 5:
            return c0 * l_M_tilde * (l_M_tilde - c1) ** 2 / (c2 + c3 * l_M_tilde) ** 3 * exp(-(l_M_tilde - c1) ** 2 / (2 * (c2 + c3 * l_M_tilde) ** 2))
        elif argindex == 6:
            return exp(-(l_M_tilde - c5) ** 2 / (2 * (c6 + c7 * l_M_tilde) ** 2))
        elif argindex == 7:
            return c4 * (l_M_tilde - c5) / (c6 + c7 * l_M_tilde) ** 2 * exp(-(l_M_tilde - c5) ** 2 / (2 * (c6 + c7 * l_M_tilde) ** 2))
        elif argindex == 8:
            return c4 * (l_M_tilde - c5) ** 2 / (c6 + c7 * l_M_tilde) ** 3 * exp(-(l_M_tilde - c5) ** 2 / (2 * (c6 + c7 * l_M_tilde) ** 2))
        elif argindex == 9:
            return c4 * l_M_tilde * (l_M_tilde - c5) ** 2 / (c6 + c7 * l_M_tilde) ** 3 * exp(-(l_M_tilde - c5) ** 2 / (2 * (c6 + c7 * l_M_tilde) ** 2))
        elif argindex == 10:
            return exp(-(l_M_tilde - c9) ** 2 / (2 * (c10 + c11 * l_M_tilde) ** 2))
        elif argindex == 11:
            return c8 * (l_M_tilde - c9) / (c10 + c11 * l_M_tilde) ** 2 * exp(-(l_M_tilde - c9) ** 2 / (2 * (c10 + c11 * l_M_tilde) ** 2))
        elif argindex == 12:
            return c8 * (l_M_tilde - c9) ** 2 / (c10 + c11 * l_M_tilde) ** 3 * exp(-(l_M_tilde - c9) ** 2 / (2 * (c10 + c11 * l_M_tilde) ** 2))
        elif argindex == 13:
            return c8 * l_M_tilde * (l_M_tilde - c9) ** 2 / (c10 + c11 * l_M_tilde) ** 3 * exp(-(l_M_tilde - c9) ** 2 / (2 * (c10 + c11 * l_M_tilde) ** 2))
        raise ArgumentIndexError(self, argindex)

    def _latex(self, printer):
        l_M_tilde = self.args[0]
        _l_M_tilde = printer._print(l_M_tilde)
        return '\\operatorname{fl}^M_{act} \\left( %s \\right)' % _l_M_tilde

class FiberForceVelocityDeGroote2016(CharacteristicCurveFunction):
    _symox_user_hooks = True

    @classmethod
    def with_defaults(cls, v_M_tilde):
        c0 = Float('-0.318')
        c1 = Float('-8.149')
        c2 = Float('-0.374')
        c3 = Float('0.886')
        return cls(v_M_tilde, c0, c1, c2, c3)

    @classmethod
    def eval(cls, v_M_tilde, c0, c1, c2, c3):
        pass

    def _eval_evalf(self, prec):
        return _evalf_at(self.doit(deep=False, evaluate=False), prec)

    def doit(self, deep=True, evaluate=True, **hints):
        (v_M_tilde, *constants) = self.args
        if deep:
            hints['evaluate'] = evaluate
            v_M_tilde = _doit_argument(v_M_tilde, deep, hints)
            (c0, c1, c2, c3) = [_doit_argument(c, deep, hints) for c in constants]
        else:
            (c0, c1, c2, c3) = constants
        if evaluate:
            return c0 * log(c1 * v_M_tilde + c2 + sqrt((c1 * v_M_tilde + c2) ** 2 + 1)) + c3
        return c0 * log(c1 * v_M_tilde + c2 + sqrt(UnevaluatedExpr(c1 * v_M_tilde + c2) ** 2 + 1)) + c3

    def fdiff(self, argindex=1):
        (v_M_tilde, c0, c1, c2, c3) = self.args
        if argindex == 1:
            return c0 * c1 / sqrt(UnevaluatedExpr(c1 * v_M_tilde + c2) ** 2 + 1)
        elif argindex == 2:
            return log(c1 * v_M_tilde + c2 + sqrt(UnevaluatedExpr(c1 * v_M_tilde + c2) ** 2 + 1))
        elif argindex == 3:
            return c0 * v_M_tilde / sqrt(UnevaluatedExpr(c1 * v_M_tilde + c2) ** 2 + 1)
        elif argindex == 4:
            return c0 / sqrt(UnevaluatedExpr(c1 * v_M_tilde + c2) ** 2 + 1)
        elif argindex == 5:
            return Integer(1)
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return FiberForceVelocityInverseDeGroote2016

    def _latex(self, printer):
        v_M_tilde = self.args[0]
        _v_M_tilde = printer._print(v_M_tilde)
        return '\\operatorname{fv}^M \\left( %s \\right)' % _v_M_tilde

class FiberForceVelocityInverseDeGroote2016(CharacteristicCurveFunction):
    _symox_user_hooks = True

    @classmethod
    def with_defaults(cls, fv_M):
        c0 = Float('-0.318')
        c1 = Float('-8.149')
        c2 = Float('-0.374')
        c3 = Float('0.886')
        return cls(fv_M, c0, c1, c2, c3)

    @classmethod
    def eval(cls, fv_M, c0, c1, c2, c3):
        pass

    def _eval_evalf(self, prec):
        return _evalf_at(self.doit(deep=False, evaluate=False), prec)

    def doit(self, deep=True, evaluate=True, **hints):
        (fv_M, *constants) = self.args
        if deep:
            hints['evaluate'] = evaluate
            fv_M = _doit_argument(fv_M, deep, hints)
            (c0, c1, c2, c3) = [_doit_argument(c, deep, hints) for c in constants]
        else:
            (c0, c1, c2, c3) = constants
        if evaluate:
            return (sinh((fv_M - c3) / c0) - c2) / c1
        return (sinh(UnevaluatedExpr(fv_M - c3) / c0) - c2) / c1

    def fdiff(self, argindex=1):
        (fv_M, c0, c1, c2, c3) = self.args
        if argindex == 1:
            return cosh((fv_M - c3) / c0) / (c0 * c1)
        elif argindex == 2:
            return (c3 - fv_M) * cosh((fv_M - c3) / c0) / (c0 ** 2 * c1)
        elif argindex == 3:
            return (c2 - sinh((fv_M - c3) / c0)) / c1 ** 2
        elif argindex == 4:
            return -1 / c1
        elif argindex == 5:
            return -cosh((fv_M - c3) / c0) / (c0 * c1)
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        return FiberForceVelocityDeGroote2016

    def _latex(self, printer):
        fv_M = self.args[0]
        _fv_M = printer._print(fv_M)
        return '\\left( \\operatorname{fv}^M \\right)^{-1} \\left( %s \\right)' % _fv_M

@dataclass(frozen=True)
class CharacteristicCurveCollection:
    tendon_force_length: CharacteristicCurveFunction
    tendon_force_length_inverse: CharacteristicCurveFunction
    fiber_force_length_passive: CharacteristicCurveFunction
    fiber_force_length_passive_inverse: CharacteristicCurveFunction
    fiber_force_length_active: CharacteristicCurveFunction
    fiber_force_velocity: CharacteristicCurveFunction
    fiber_force_velocity_inverse: CharacteristicCurveFunction

    def __iter__(self):
        yield self.tendon_force_length
        yield self.tendon_force_length_inverse
        yield self.fiber_force_length_passive
        yield self.fiber_force_length_passive_inverse
        yield self.fiber_force_length_active
        yield self.fiber_force_velocity
        yield self.fiber_force_velocity_inverse
