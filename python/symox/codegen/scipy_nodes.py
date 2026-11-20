"""Symbolic nodes for SciPy special functions (port of ``sympy.codegen.scipy_nodes``)."""
import symox as S
from .numpy_nodes import _function_new


def _cosm1(x, *, evaluate=True):
    return S.Add(S.cos(x, evaluate=evaluate), -S.S.One, evaluate=evaluate)


def _argindex_error(self, argindex):
    from ..core.function import ArgumentIndexError
    raise ArgumentIndexError(self, argindex)


class cosm1(S.Function):
    """ Minus one plus cosine of x, i.e. cos(x) - 1. For use when x is close to zero.

    Helper class for use with e.g. scipy.special.cosm1
    See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.cosm1.html
    """
    nargs = 1
    _symox_user_hooks = True

    def __new__(cls, *args, **options):
        return _function_new(super(), cls, args, options)

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return -S.sin(*self.args)
        else:
            _argindex_error(self, argindex)

    def _eval_rewrite_as_cos(self, x, **kwargs):
        return _cosm1(x)

    def _eval_evalf(self, *args, **kwargs):
        return self.rewrite(S.cos).evalf(*args, **kwargs)

    def _eval_simplify(self, **kwargs):
        x, = self.args
        candidate = _cosm1(x.simplify(**kwargs))
        if candidate != _cosm1(x, evaluate=False):
            return candidate
        else:
            return cosm1(x)


def _powm1(x, y, *, evaluate=True):
    return S.Add(S.Pow(x, y, evaluate=evaluate), -S.S.One, evaluate=evaluate)


class powm1(S.Function):
    """ Minus one plus x to the power of y, i.e. x**y - 1. For use when x is close to one or y is close to zero.

    Helper class for use with e.g. scipy.special.powm1
    See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.powm1.html
    """
    nargs = 2
    _symox_user_hooks = True

    def __new__(cls, *args, **options):
        return _function_new(super(), cls, args, options)

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            return S.Pow(self.args[0], self.args[1])*self.args[1]/self.args[0]
        elif argindex == 2:
            return S.log(self.args[0])*S.Pow(*self.args)
        else:
            _argindex_error(self, argindex)

    def _eval_rewrite_as_Pow(self, x, y, **kwargs):
        return _powm1(x, y)

    def _eval_evalf(self, *args, **kwargs):
        return self.rewrite(S.Pow).evalf(*args, **kwargs)

    def _eval_simplify(self, **kwargs):
        x, y = self.args
        candidate = _powm1(x.simplify(**kwargs), y.simplify(**kwargs))
        if candidate != _powm1(x, y, evaluate=False):
            return candidate
        else:
            return powm1(x, y)


__all__ = ['cosm1', 'powm1']
