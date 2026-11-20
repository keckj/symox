"""Symbolic nodes for NumPy functions (port of ``sympy.codegen.numpy_nodes``)."""
import symox as S
from .ast import Token, none


def _logaddexp(x1, x2, *, evaluate=True):
    return S.log(S.Add(S.exp(x1, evaluate=evaluate), S.exp(x2, evaluate=evaluate), evaluate=evaluate))


_two = S.Integer(2)
_ln2 = S.log(_two)


def _lb(x, *, evaluate=True):
    return S.log(x, evaluate=evaluate)/_ln2


def _exp2(x, *, evaluate=True):
    return S.Pow(_two, x, evaluate=evaluate)


def _logaddexp2(x1, x2, *, evaluate=True):
    return _lb(S.Add(_exp2(x1, evaluate=evaluate),
                     _exp2(x2, evaluate=evaluate), evaluate=evaluate))


def _argindex_error(self, argindex):
    from ..core.function import ArgumentIndexError
    raise ArgumentIndexError(self, argindex)


def _sorted_args(args):
    from ..core.sorting import default_sort_key
    return sorted(args, key=default_sort_key)


def _should_evalf(arg):
    # SymPy's Function._should_evalf.
    if arg.is_Float:
        return arg._prec
    if not arg.is_Add:
        return -1
    from ..core.evalf import pure_complex
    m = pure_complex(arg)
    if m is None:
        return -1
    precs = [i._prec for i in m if i.is_Float]
    precs.append(-1)
    return max(precs)


def _function_new(base, cls, args, options):
    """SymPy's ``Function.__new__``: evaluated applications on Float
    arguments (or pure complex numbers with Float parts) are evaluated
    numerically at the largest argument precision."""
    from ..core.parameters import global_parameters
    evaluate = options.get('evaluate', global_parameters.evaluate)
    result = base.__new__(cls, *args, **options)
    if evaluate and isinstance(result, cls) and result.args:
        precs = [_should_evalf(a) for a in result.args]
        if min(precs) > 0:
            from mpmath.libmp import prec_to_dps
            result = result.evalf(prec_to_dps(max(precs)))
    return result


class logaddexp(S.Function):
    """ Logarithm of the sum of exponentiations of the inputs.

    Helper class for use with e.g. numpy.logaddexp

    See Also
    ========

    https://numpy.org/doc/stable/reference/generated/numpy.logaddexp.html
    """
    nargs = 2
    _symox_user_hooks = True

    def __new__(cls, *args):
        # SymPy's signature takes no options (evaluate=False is a TypeError).
        return _function_new(super(), cls, _sorted_args(args), {})

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            wrt, other = self.args
        elif argindex == 2:
            other, wrt = self.args
        else:
            _argindex_error(self, argindex)
        return S.S.One/(S.S.One + S.exp(other-wrt))

    def _eval_rewrite_as_log(self, x1, x2, **kwargs):
        return _logaddexp(x1, x2)

    def _eval_evalf(self, *args, **kwargs):
        return self.rewrite(S.log).evalf(*args, **kwargs)

    def _eval_simplify(self, *args, **kwargs):
        a, b = (x.simplify(**kwargs) for x in self.args)
        candidate = _logaddexp(a, b)
        if candidate != _logaddexp(a, b, evaluate=False):
            return candidate
        else:
            return logaddexp(a, b)


class logaddexp2(S.Function):
    """ Logarithm of the sum of exponentiations of the inputs in base-2.

    Helper class for use with e.g. numpy.logaddexp2

    See Also
    ========

    https://numpy.org/doc/stable/reference/generated/numpy.logaddexp2.html
    """
    nargs = 2
    _symox_user_hooks = True

    def __new__(cls, *args):
        # SymPy's signature takes no options (evaluate=False is a TypeError).
        return _function_new(super(), cls, _sorted_args(args), {})

    def fdiff(self, argindex=1):
        """
        Returns the first derivative of this function.
        """
        if argindex == 1:
            wrt, other = self.args
        elif argindex == 2:
            other, wrt = self.args
        else:
            _argindex_error(self, argindex)
        return S.S.One/(S.S.One + _exp2(other-wrt))

    def _eval_rewrite_as_log(self, x1, x2, **kwargs):
        return _logaddexp2(x1, x2)

    def _eval_evalf(self, *args, **kwargs):
        return self.rewrite(S.log).evalf(*args, **kwargs)

    def _eval_simplify(self, *args, **kwargs):
        a, b = (x.simplify(**kwargs).factor() for x in self.args)
        candidate = _logaddexp2(a, b)
        if candidate != _logaddexp2(a, b, evaluate=False):
            return candidate
        else:
            return logaddexp2(a, b)


class _AxisReduction(Token):
    _fields = ('array', 'axis')
    defaults = {'axis': none}

    def __init__(self, *args, **kwargs):
        # SymPy's Token.__new__ checks.
        if len(args) > len(self._fields):
            raise ValueError("Too many arguments (%d), expected at most %d" % (len(args), len(self._fields)))
        unknown = set(kwargs) - set(self._fields)
        if unknown:
            raise ValueError("Unknown keyword arguments: %s" % ' '.join(unknown))
        for field in self._fields[:len(args)]:
            if field in kwargs:
                raise TypeError("Got multiple values for attribute %r" % field)
        for field in self._fields[len(args):]:
            if field not in kwargs and field not in self.defaults:
                raise TypeError("No value for %r given and attribute has no default" % field)
        super().__init__(*args, **kwargs)
        # SymPy's Token keeps a field without a constructor as given (a list
        # stays a list).
        array = args[0] if args else kwargs['array']
        if array is not None:
            object.__setattr__(self, 'array', array)

    def _coerce(self, field, value):
        # SymPy's Token._construct: None gives the default, the axis is
        # sympified (``_construct_axis = staticmethod(sympify)``).
        if value is None or value is none:
            return self.defaults.get(field, none)
        if field == 'axis':
            return S.sympify(value)
        return value

    def __hash__(self):
        # Like SymPy's Basic hash: an unhashable array (a list) raises TypeError.
        return hash((type(self), self.array, self.axis))

    def __str__(self):
        from ..printing.str import sstr
        return '%s(%s)' % (type(self).__name__, ', '.join(
            sstr(getattr(self, k)) if k == 'array' else '%s=%s' % (k, sstr(getattr(self, k)))
            for k in self._fields if not (k in self.defaults and getattr(self, k) == self.defaults[k])))

    def __repr__(self):
        from ..printing.repr import srepr
        return '%s(%s)' % (type(self).__name__, ', '.join(
            srepr(getattr(self, k)) if k == 'array' else '%s=%s' % (k, srepr(getattr(self, k)))
            for k in self._fields if not (k in self.defaults and getattr(self, k) == self.defaults[k])))


class amin(_AxisReduction):
    """ Minimum value along an axis.

    Helper class for use with e.g. numpy.amin


    See Also
    ========

    https://numpy.org/doc/stable/reference/generated/numpy.amin.html
    """


class amax(_AxisReduction):
    """ Maximum value along an axis.

    Helper class for use with e.g. numpy.amax


    See Also
    ========

    https://numpy.org/doc/stable/reference/generated/numpy.amax.html
    """


class maximum(S.Function):
    """ Element-wise maximum of array elements.

    Helper class for use with e.g. numpy.maximum


    See Also
    ========

    https://numpy.org/doc/stable/reference/generated/numpy.maximum.html
    """
    _symox_user_hooks = True

    def __new__(cls, *args, **options):
        return _function_new(super(), cls, args, options)

    def _eval_rewrite_as_Max(self, *args, **kwargs):
        return S.Max(*self.args)


class minimum(S.Function):
    """ Element-wise minimum of array elements.

    Helper class for use with e.g. numpy.minimum


    See Also
    ========

    https://numpy.org/doc/stable/reference/generated/numpy.minimum.html
    """
    _symox_user_hooks = True

    def __new__(cls, *args, **options):
        return _function_new(super(), cls, args, options)

    def _eval_rewrite_as_Min(self, *args, **kwargs):
        return S.Min(*self.args)


__all__ = ['logaddexp', 'logaddexp2', 'amin', 'amax', 'maximum', 'minimum']
