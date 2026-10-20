"""NumPy ufunc evaluation with the native real-array result conventions."""
import numpy as np
from functools import lru_cache


@lru_cache(maxsize=32)
def _evaluator(arguments,expression):
    from .utilities.lambdify import _compile_object
    return _compile_object(arguments,expression)


def compile(arguments,expression):
    try:
        evaluator = _evaluator(tuple(arguments),expression)
    except (NotImplementedError,TypeError,ValueError,SyntaxError,AttributeError,RecursionError):
        return lambda *values:NotImplemented

    def evaluate(*values):
        arrays = []
        for value in values:
            if type(value) is np.ndarray:
                # The VM promotes other real dtypes to float64. Never gain
                # speed by evaluating float32 inputs at lower precision.
                if value.dtype != np.dtype('float64'):
                    return NotImplemented
                arrays.append(value)
            elif type(value) not in (int,float):
                return NotImplemented
        if not arrays or not any(value.ndim for value in arrays):
            return NotImplemented
        try:
            np.broadcast_shapes(*(value.shape for value in arrays))
        except ValueError:
            # Keep the VM's established equal-flat-size behavior as well.
            return NotImplemented
        non_scalar = [value for value in arrays if value.ndim]
        legacy_shape = non_scalar[0].shape
        if any(value.shape != legacy_shape for value in non_scalar[1:]):
            legacy_shape = None
        try:
            # Scalar operands must use the same IEEE float64 rules as the VM,
            # rather than Python integer powers/division at singular points.
            operands = tuple(np.float64(v) if type(v) in (int,float) else v for v in values)
            with np.errstate(all='ignore'):
                result = np.asarray(evaluator(*operands),dtype=np.float64)
            # For newly accepted broadcasting shapes, an unused argument must
            # not enlarge the reference result. Equal-shape VM calls retain
            # their established constant/scalar-output broadcasting extension.
            shape = result.shape if legacy_shape is None else legacy_shape
            if (result.shape == shape and result.flags.c_contiguous
                    and result.flags.writeable and result.flags.owndata
                    and all(result is not value for value in arrays)):
                return result
            # Constants and identity functions retain writable, independent
            # float64 arrays instead of exposing input aliases or readonly views.
            return np.array(np.broadcast_to(result,shape),dtype=np.float64,copy=True,order='C')
        except (NotImplementedError,TypeError,ValueError,ArithmeticError,AttributeError):
            return NotImplemented

    return evaluate
