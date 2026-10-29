"""Dispatch public scalar operations and finite matrix derivatives to arrays."""
from functools import wraps

from ... import _core
from ..._dispatch_arguments import FirstArgument
from ..._node_types import is_native_node
from ._array import NDimArray
from .arrayop import derive_by_array


def _array_variable(value):
    return isinstance(value, (NDimArray, _core.Matrix, _core._MatrixExpr)) or bool(
        getattr(value, '_symox_sparse', False))


def _check_derivative_options(options):
    if any(key != 'evaluate' or value is not True for key, value in options.items()):
        raise NotImplementedError('Array-variable derivatives require explicit evaluation')


def _diff_dispatch(original):
    first=FirstArgument(original)
    @wraps(original)
    def diff(*args,**options):
        expression=first.get(args,options)
        if isinstance(expression, NDimArray):
            expression,variables,options=first.split(args,options)
            return expression.diff(*variables,**options)
        variables=args[1:]
        if variables and options.get('evaluate', True) is True and (
                isinstance(expression, (_core._MatrixExpr, _core.Matrix))
                or getattr(expression, '_symox_sparse', False)):
            # Symbolic orders by scalar Symbols: SymPy's held ArrayDerivative.
            from .array_derivatives import symbolic_order_derivative
            held = symbolic_order_derivative(expression, variables)
            if held is not None:
                return held
        if variables and isinstance(expression, _core._MatrixExpr) and any(
                isinstance(v, _core._MatrixExpr) or (isinstance(v, tuple) and v and isinstance(v[0], _core._MatrixExpr))
                for v in variables):
            # ArrayDerivative._call_derive_default: a matrix expression
            # without the matrix variable keeps Derivative(expr, A).
            from .array_derivatives import held_matrix_derivative
            _check_derivative_options(options)
            held = held_matrix_derivative(expression, variables)
            if held is not None:
                return held
        if len(variables) == 1 and _array_variable(variables[0]):
            _check_derivative_options(options)
            variable = variables[0]
            if (isinstance(variable, _core._MatrixExpr) and variable._kind_tag == 'MatrixSymbol'
                    and not isinstance(expression, (NDimArray, _core.Matrix))
                    and not getattr(expression, '_symox_sparse', False)):
                # SymPy's _matrix_derivative (matrix-line algorithm).
                from ..._matrix_derivatives import _Unsupported, matrix_derivative
                try:
                    return matrix_derivative(expression, variable)
                except _Unsupported:
                    pass
            if (isinstance(variable, _core.Matrix) or getattr(variable, '_symox_sparse', False)) and not (
                    isinstance(expression, (NDimArray, _core.Matrix, _core._MatrixExpr))
                    or getattr(expression, '_symox_sparse', False) or getattr(expression, 'is_Matrix', False)):
                # SymPy's ArrayDerivative._call_derive_scalar_by_matrix: a
                # scalar derived by a matrix stays an (immutable) matrix.
                if getattr(variable, '_symox_sparse', False):
                    # Sparse applyfunc also maps the implicit zero, which
                    # SymPy's dok-based applyfunc never differentiates by.
                    return variable.as_immutable().applyfunc(
                        lambda entry: entry if entry == 0 else original(expression, entry))
                return variable.as_immutable().applyfunc(lambda entry: original(expression, entry))
            result = derive_by_array(expression, variable)
            # SymPy's MatrixCalculus.diff: a mutable matrix gives a mutable array.
            if (isinstance(expression, _core.Matrix) and not getattr(expression, '_immutable', False)
                    and hasattr(result, 'as_mutable')):
                result = result.as_mutable()
            return result
        return original(*args, **options)
    return diff


def _entrywise(original):
    first=FirstArgument(original)
    @wraps(original)
    def operation(*args, **options):
        if args and is_native_node(args[0]):
            return original(*args, **options)
        expression=first.get(args,options)
        if (isinstance(expression, (NDimArray, _core.Matrix))
                or getattr(expression, '_symox_sparse', False)):
            expression,rest,options=first.split(args,options)
            if isinstance(expression, NDimArray) or hasattr(expression, 'as_immutable'):
                # SymPy sympifies the argument first: mutable arrays and
                # matrices (dense or sparse) become their immutable
                # counterpart.
                expression=expression.as_immutable()
            return expression.applyfunc(lambda entry: original(entry, *rest, **options))
        return original(*args, **options)
    from ... import _BuiltinFunctionHead
    from ..._function_heads import is_builtin_head
    if is_builtin_head(original):
        return _BuiltinFunctionHead(original.__name__, operation)
    return operation


def install(namespace):
    if namespace.get('_array_adapters_installed', False):
        return
    namespace['_array_adapters_installed'] = True
    for name in ('simplify', 'expand', 'factor', 'trigsimp', 'powsimp', 'ratsimp',
                 'radsimp', 'cancel', 'together', 'conjugate', 'N'):
        if name in namespace:
            namespace[name] = _entrywise(namespace[name])
    namespace['diff'] = _diff_dispatch(namespace['diff'])
    for cls in (_core.Expr, _core.Matrix, _core._MatrixExpr):
        cls.diff = _diff_dispatch(cls.diff)
    from ..._sparse import _SparseBase
    _SparseBase.diff = _diff_dispatch(_SparseBase.diff)
