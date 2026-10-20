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
        if len(variables) == 1 and _array_variable(variables[0]):
            _check_derivative_options(options)
            return derive_by_array(expression, variables[0])
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
            if isinstance(expression, NDimArray):
                # SymPy sympifies the argument first: mutable arrays become
                # their immutable counterpart.
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
