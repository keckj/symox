"""Standalone JAX source generation with broadcasting, autodiff, and JIT."""

from __future__ import annotations

import builtins
import itertools
import keyword
import linecache

from .pytorch import _arguments, _contains_conditionals, _cse_preserving_assumptions, _expressions, _functions


_RUNTIME = '''import jax
import jax.numpy as jnp
import jax.scipy.special as jsp


def _sx_tensor(value, reference):
    if isinstance(value, (jax.Array, jax.core.Tracer)):
        return value
    dtype = reference.dtype if reference is not None else jnp.asarray(0.0).dtype
    if isinstance(value, bool):
        dtype = jnp.bool_
    elif isinstance(value, complex):
        dtype = jnp.promote_types(dtype, jnp.asarray(1j).dtype)
    elif dtype == jnp.bool_:
        dtype = jnp.asarray(0.0).dtype
    elif isinstance(value, float) and not jnp.issubdtype(dtype, jnp.inexact):
        dtype = jnp.asarray(0.0).dtype
    value = jnp.asarray(value, dtype=dtype)
    if isinstance(reference, jax.Array) and not isinstance(reference, jax.core.Tracer):
        value = jax.device_put(value, reference.sharding)
    return value


def _sx_prepare(values, dtype, device):
    if isinstance(device, str):
        devices = jax.devices(device)
        if not devices:
            raise ValueError("No JAX device found for backend " + device)
        device = devices[0]
    converted = tuple(jnp.asarray(value, dtype=dtype) for value in values)
    target_dtype = jnp.dtype(dtype) if dtype is not None else (jnp.result_type(*converted) if converted else jnp.asarray(0.0).dtype)
    converted = tuple(value if dtype is None and value.dtype == jnp.bool_ else value.astype(target_dtype) for value in converted)
    if device is not None:
        converted = tuple(jax.device_put(value, device) for value in converted)
    if converted:
        converted = tuple(jnp.broadcast_arrays(*converted))
        shape = converted[0].shape
    else:
        shape = ()
    reference = converted[0] if converted and converted[0].dtype == target_dtype else jnp.zeros((), dtype=target_dtype)
    if device is not None:
        reference = jax.device_put(reference, device)
    return converted, reference, shape


def _sx_output(value, reference, shape):
    return jnp.broadcast_to(_sx_tensor(value, reference), shape)


def _sx_dense(values, reference, batch_shape, expression_shape):
    if not values:
        result = jnp.empty(tuple(batch_shape) + expression_shape, dtype=reference.dtype)
        if isinstance(reference, jax.Array) and not isinstance(reference, jax.core.Tracer):
            result = jax.device_put(result, reference.sharding)
        return result
    return jnp.stack(values, axis=-1).reshape(tuple(batch_shape) + expression_shape)


def _sx_group(value, length):
    if not isinstance(value, (tuple, list)) or len(value) != length:
        raise ValueError("Nested JAX function arguments must match the declared list/tuple structure")
    return value
'''

_GENERATED = itertools.count(1)


def _generate(name, args, exprs, cse, user_functions, allow_callables=False):
    from symox import Symbol, _core

    if not isinstance(name, str) or not name.isidentifier() or keyword.iskeyword(name) or name.startswith("_sx_") or name in {"jax", "jnp", "jsp"} or name in vars(builtins):
        raise ValueError("Function name must be a Python identifier outside the reserved '_sx_' namespace")
    if not isinstance(cse, bool):
        raise TypeError("cse must be a bool")
    top, symbols, unpacking, groups = _arguments(args, backend="JAX")
    expressions, structure = _expressions(exprs)
    free = set().union(*(expr.free_symbols for expr in expressions)) if expressions else set()
    unknown = free.difference(symbols)
    if unknown:
        raise ValueError("Expression contains unbound symbols: " + ", ".join(sorted(map(str, unknown))))
    dummy_symbols = [Symbol(f"_sx_input_{index}", **symbol.assumptions0) for index, symbol in enumerate(symbols)]
    dummies = dict(zip(symbols, dummy_symbols))
    expressions = [expr.xreplace(dummies) for expr in expressions]
    functions, imports, namespace = _functions(user_functions, allow_callables)
    if name in {module.split(".")[0] for module in imports}:
        raise ValueError("Function name conflicts with a custom function module")
    replacements = []
    if cse and expressions and not _contains_conditionals(expressions):
        replacements, expressions = _cse_preserving_assumptions(expressions, "_sx_cse")
    names = {symbol: f"_sx_v_{index}" for index, symbol in enumerate(dummy_symbols)}
    names.update({symbol: f"_sx_t_{index}" for index, (symbol, _) in enumerate(replacements)})

    def emit(expr):
        return _core._jaxcode(expr, symbol_names=names, user_functions=functions, reference="_sx_ref")

    rendered = [f"_sx_output({emit(expr)}, _sx_ref, _sx_shape)" for expr in expressions]

    def output(node):
        kind = node[0]
        if kind == "leaf":
            return rendered[node[1]]
        if kind == "dense":
            items = ", ".join(output(child) for child in node[2])
            return f"_sx_dense([{items}], _sx_ref, _sx_shape, {node[1]!r})"
        items = ", ".join(output(child) for child in node[1])
        if kind == "list":
            return "[" + items + "]"
        return "(" + items + ("," if len(node[1]) == 1 else "") + ")"

    parameters = ", ".join(f"_sx_arg_{index}" for index in range(len(top)))
    parameters += (", *, " if parameters else "*, ") + "dtype=None, device=None"
    lines = [f"def {name}({parameters}):"]
    for path, length in groups:
        lines.append(f"    _sx_group({path}, {length})")
    for safe, path in unpacking:
        lines.append(f"    {safe} = {path}")
    variables = ", ".join(safe for safe, _ in unpacking)
    input_tuple = "(" + variables + ("," if len(unpacking) == 1 else "") + ")"
    lines.append(f"    _sx_values, _sx_ref, _sx_shape = _sx_prepare({input_tuple}, dtype, device)")
    for index, (safe, _) in enumerate(unpacking):
        lines.append(f"    {safe} = _sx_values[{index}]")
    for symbol, expression in replacements:
        lines.append(f"    {names[symbol]} = {emit(expression)}")
    lines.append("    return " + output(structure))
    source = _RUNTIME + "\n" + "".join(f"import {module}\n" for module in sorted(imports)) + "\n" + "\n".join(lines) + "\n"
    return source, namespace


def jaxcode_function(name, args, exprs, *, cse=True, user_functions=None):
    """Return standalone JAX Python source, with batch axes before matrix axes.

    The generated function accepts keyword-only ``dtype`` and ``device``.
    Devices can be JAX Device objects or backend names such as ``'cpu'``.
    Dtypes obey JAX's configuration, including ``jax_enable_x64``. Lists and
    tuples retain their structure; matrices and NumPy arrays become dense arrays.
    """
    return _generate(name, args, exprs, cse, user_functions)[0]


def jax_lambdify(args, exprs, *, cse=True, jit=False, jit_options=None, user_functions=None):
    """Build an inspectable JAX function, optionally transformed with jax.jit.

    ``dtype`` and ``device`` are static keyword arguments under JIT. Other JIT
    settings pass through ``jit_options``. Errors propagate without a fallback.
    ``user_functions`` accepts dotted Python names or callable implementations.
    CSE is disabled for expressions containing Piecewise to preserve its masks.
    """
    if not isinstance(jit, bool):
        raise TypeError("jit must be a bool")
    if jit_options is not None and not isinstance(jit_options, dict):
        raise TypeError("jit_options must be a dict")
    if jit_options and not jit:
        raise ValueError("jit_options require jit=True")
    name = "_symox_jaxgenerated"
    source, namespace = _generate(name, args, exprs, cse, user_functions, allow_callables=True)
    filename = f"<symox-jax-{next(_GENERATED)}>"
    exec(builtins.compile(source, filename, "exec"), namespace)
    function = namespace[name]
    linecache.cache[filename] = (len(source), None, source.splitlines(keepends=True), filename)
    function.__source__ = source
    function.__doc__ = "Generated JAX function. Expression dimensions follow broadcast batch dimensions."
    if jit:
        import jax
        options = dict(jit_options or {})
        static = options.get("static_argnames") or ()
        static = (static,) if isinstance(static, str) else tuple(static)
        options["static_argnames"] = tuple(dict.fromkeys((*static, "dtype", "device")))
        original = function
        function = jax.jit(original, **options)
        function.__source__ = source
        function.__wrapped__ = original
    return function


__all__ = ["jaxcode_function", "jax_lambdify"]
