"""NumPy and SciPy expression printers (``sympy.printing.numpy``)."""
from .pycode import ArrayPrinter


class NumPyPrinter(ArrayPrinter):
    """Python with NumPy: vectorized piecewise functions, logical operators,
    and the ``sympy.codegen.numpy_nodes`` functions."""


# SymPy's _known_functions_scipy_special and _known_constants_scipy_constants.
_known_functions_scipy_special = {
    'Ei': 'expi',
    'erf': 'erf',
    'erfc': 'erfc',
    'besselj': 'jv',
    'bessely': 'yv',
    'besseli': 'iv',
    'besselk': 'kv',
    'cosm1': 'cosm1',
    'powm1': 'powm1',
    'factorial': 'factorial',
    'gamma': 'gamma',
    'loggamma': 'gammaln',
    'digamma': 'psi',
    'polygamma': 'polygamma',
    'RisingFactorial': 'poch',
    'jacobi': 'eval_jacobi',
    'gegenbauer': 'eval_gegenbauer',
    'chebyshevt': 'eval_chebyt',
    'chebyshevu': 'eval_chebyu',
    'legendre': 'eval_legendre',
    'hermite': 'eval_hermite',
    'laguerre': 'eval_laguerre',
    'assoc_laguerre': 'eval_genlaguerre',
    'beta': 'beta',
    'LambertW': 'lambertw',
}
_known_constants_scipy_constants = {
    'GoldenRatio': 'golden_ratio',
    'Pi': 'pi',
}


def _unpack_integral_limits(integral_expr):
    """ helper function for _print_Integral that
        - accepts an Integral expression
        - returns a tuple of
           - a list variables of integration
           - a list of tuples of the upper and lower limits of integration
    """
    integration_vars = []
    limits = []
    for integration_range in integral_expr.limits:
        if len(integration_range) == 3:
            integration_var, lower_limit, upper_limit = integration_range
        else:
            raise NotImplementedError("Only definite integrals are supported")
        integration_vars.append(integration_var)
        limits.append((lower_limit, upper_limit))
    return integration_vars, limits


class SciPyPrinter(NumPyPrinter):
    """Python with SciPy and NumPy (port of SymPy's ``SciPyPrinter``)."""
    language = 'Python with SciPy and NumPy'

    def _function(self, name):
        if name in _known_functions_scipy_special:
            return 'scipy.special.' + _known_functions_scipy_special[name]
        return super()._function(name)

    # The scipy_nodes functions are printed by name (not as unsupported).
    def _print_cosm1(self, expr):
        return 'scipy.special.cosm1(%s)' % ', '.join(self._print(a) for a in expr.args)

    def _print_powm1(self, expr):
        return 'scipy.special.powm1(%s)' % ', '.join(self._print(a) for a in expr.args)

    def _print_Pi(self, expr):
        return 'scipy.constants.pi'

    def _print_GoldenRatio(self, expr):
        return 'scipy.constants.golden_ratio'

    # SciPy's lpmv has a different order of arguments from assoc_legendre
    def _print_assoc_legendre(self, expr):
        return "{0}({2}, {1}, {3})".format(
            'scipy.special.lpmv',
            self._print(expr.args[0]),
            self._print(expr.args[1]),
            self._print(expr.args[2]))

    def _print_lowergamma(self, expr):
        return "{0}({2})*{1}({2}, {3})".format(
            'scipy.special.gamma',
            'scipy.special.gammainc',
            self._print(expr.args[0]),
            self._print(expr.args[1]))

    def _print_uppergamma(self, expr):
        return "{0}({2})*{1}({2}, {3})".format(
            'scipy.special.gamma',
            'scipy.special.gammaincc',
            self._print(expr.args[0]),
            self._print(expr.args[1]))

    def _print_betainc(self, expr):
        betainc = 'scipy.special.betainc'
        beta = 'scipy.special.beta'
        args = [self._print(arg) for arg in expr.args]
        return f"({betainc}({args[0]}, {args[1]}, {args[3]}) - {betainc}({args[0]}, {args[1]}, {args[2]})) \
            * {beta}({args[0]}, {args[1]})"

    def _print_betainc_regularized(self, expr):
        return "{0}({1}, {2}, {4}) - {0}({1}, {2}, {3})".format(
            'scipy.special.betainc',
            self._print(expr.args[0]),
            self._print(expr.args[1]),
            self._print(expr.args[2]),
            self._print(expr.args[3]))

    def _print_fresnels(self, expr):
        return "{}({})[0]".format('scipy.special.fresnel', self._print(expr.args[0]))

    def _print_fresnelc(self, expr):
        return "{}({})[1]".format('scipy.special.fresnel', self._print(expr.args[0]))

    def _print_airyai(self, expr):
        return "{}({})[0]".format('scipy.special.airy', self._print(expr.args[0]))

    def _print_airyaiprime(self, expr):
        return "{}({})[1]".format('scipy.special.airy', self._print(expr.args[0]))

    def _print_airybi(self, expr):
        return "{}({})[2]".format('scipy.special.airy', self._print(expr.args[0]))

    def _print_airybiprime(self, expr):
        return "{}({})[3]".format('scipy.special.airy', self._print(expr.args[0]))

    def _print_Integral(self, e):
        integration_vars, limits = _unpack_integral_limits(e)

        if len(limits) == 1:
            # nicer (but not necessary) to prefer quad over nquad for 1D case
            module_str = 'scipy.integrate.quad'
            limit_str = "%s, %s" % tuple(map(self._print, limits[0]))
        else:
            module_str = 'scipy.integrate.nquad'
            limit_str = "({})".format(", ".join(
                "(%s, %s)" % tuple(map(self._print, l)) for l in limits))

        return "{}(lambda {}: {}, {})[0]".format(
            module_str,
            ", ".join(map(self._print, integration_vars)),
            self._print(e.args[0]),
            limit_str)

    def _print_Si(self, expr):
        return "{}({})[0]".format('scipy.special.sici', self._print(expr.args[0]))

    def _print_Ci(self, expr):
        return "{}({})[1]".format('scipy.special.sici', self._print(expr.args[0]))


__all__ = ['NumPyPrinter', 'SciPyPrinter']
