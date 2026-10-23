"""Standalone, differentiable PyTorch functions generated from native expressions.

PyTorch is imported only when a callable is built or exported source is executed.
Lists and tuples of expressions retain their structure. Matrices and NumPy arrays become dense
tensors with the expression dimensions following the broadcast batch dimensions.
"""

from __future__ import annotations

import builtins
import itertools
import keyword
import linecache
import re


# These helpers are included in exported source: evaluating that source requires
# torch, but never symox, SymPy, or NumPy. Keep tensor operations traceable by
# Dynamo and avoid extracting tensor data into Python.
_RUNTIME = '''import torch


def _sx_tensor(value, reference):
    if isinstance(value, torch.Tensor):
        return value
    dtype = reference.dtype if reference is not None else torch.get_default_dtype()
    device = reference.device if reference is not None else None
    if isinstance(value, bool):
        dtype = torch.bool
    elif isinstance(value, complex):
        dtype = torch.promote_types(dtype, torch.complex128 if dtype in (torch.float64, torch.complex128) else torch.complex64)
    elif dtype == torch.bool:
        dtype = torch.get_default_dtype()
    elif isinstance(value, float) and not (dtype.is_floating_point or dtype.is_complex):
        dtype = torch.get_default_dtype()
    if dtype.is_complex:
        # This also covers real constants promoted by complex inputs. Real
        # components avoid a constant buffer that CUDA Inductor cannot load
        # when materializing a broadcast complex output.
        real_dtype = torch.float64 if dtype == torch.complex128 else (torch.float32 if dtype == torch.complex64 else torch.float16)
        return torch.complex(torch.as_tensor(value.real, dtype=real_dtype, device=device),
                             torch.as_tensor(value.imag, dtype=real_dtype, device=device))
    return torch.as_tensor(value, dtype=dtype, device=device)


def _sx_prepare(values, dtype, device):
    tensor_values = [value for value in values if isinstance(value, torch.Tensor)]
    target_device = device if device is not None else (tensor_values[0].device if tensor_values else None)
    target_dtype = dtype
    if target_dtype is None:
        for value in tensor_values:
            target_dtype = value.dtype if target_dtype is None else torch.promote_types(target_dtype, value.dtype)
        if target_dtype is None:
            target_dtype = torch.get_default_dtype()
        for value in values:
            if isinstance(value, complex):
                target_dtype = torch.promote_types(target_dtype, torch.complex128 if target_dtype in (torch.float64, torch.complex128) else torch.complex64)
    converted = []
    for value in values:
        if isinstance(value, torch.Tensor):
            if device is None and value.device != target_device:
                raise ValueError("PyTorch function inputs must be on the same device")
            value_dtype = torch.bool if value.dtype == torch.bool and dtype is None else target_dtype
            converted.append(value.to(dtype=value_dtype, device=target_device))
        else:
            value_dtype = torch.bool if isinstance(value, bool) and dtype is None else target_dtype
            converted.append(torch.as_tensor(value, dtype=value_dtype, device=target_device))
    if converted:
        converted = tuple(torch.broadcast_tensors(*converted))
        shape = converted[0].shape
    else:
        converted = ()
        shape = ()
    reference = converted[0] if converted and converted[0].dtype == target_dtype else torch.empty((), dtype=target_dtype, device=target_device)
    return converted, reference, shape


def _sx_output(value, reference, shape):
    return _sx_tensor(value, reference).expand(shape)


def _sx_dense(values, reference, batch_shape, expression_shape):
    if not values:
        return torch.empty(tuple(batch_shape) + expression_shape, dtype=reference.dtype, device=reference.device)
    return torch.stack(values, dim=-1).reshape(tuple(batch_shape) + expression_shape)


def _sx_group(value, length):
    if not isinstance(value, (tuple, list)) or len(value) != length:
        raise ValueError("Nested PyTorch function arguments must match the declared list/tuple structure")
    return value


def _sx_matrix(value, shape):
    if not isinstance(value, torch.Tensor):
        value = torch.as_tensor(value)
    if value.ndim < 2 or tuple(value.shape[-2:]) != shape:
        raise ValueError("Matrix input trailing dimensions must match the declared MatrixSymbol shape")
    return value
'''

_GENERATED = itertools.count(1)
_DOTTED_NAME = re.compile(r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*\Z", re.ASCII)


def _arguments(args, backend="PyTorch"):
    from symox import MatrixSymbol, Symbol

    top = list(args) if isinstance(args, (list, tuple)) else [args]
    leaves = []
    unpacking = []
    groups = []

    def visit(value, path):
        if isinstance(value, MatrixSymbol) and backend in ("PyTorch", "JAX"):
            if any(not dim.is_Integer for dim in value.shape):
                raise NotImplementedError("MatrixSymbol code generation requires concrete integer dimensions")
            shape = tuple(int(dim) for dim in value.shape)
            if 0 in shape:
                raise NotImplementedError("MatrixSymbol inputs require nonempty dimensions")
            alias = f"_sx_matrix_{len(groups)}"
            groups.append((alias, (path, shape)))
            for row in range(shape[0]):
                for column in range(shape[1]):
                    entry = value[row, column]
                    if entry in leaves:
                        raise ValueError(f"Duplicate input symbol: {value}")
                    unpacking.append((f"_sx_v_{len(leaves)}", f"{alias}[..., {row}, {column}]"))
                    leaves.append(entry)
            return
        if isinstance(value, Symbol):
            if value in leaves:
                raise ValueError(f"Duplicate input symbol: {value}")
            unpacking.append((f"_sx_v_{len(leaves)}", path))
            leaves.append(value)
        elif isinstance(value, (tuple, list)):
            groups.append((path, len(value)))
            for index, item in enumerate(value):
                visit(item, f"{path}[{index}]")
        else:
            raise TypeError(f"{backend} function arguments must be Symbols or nested lists/tuples of Symbols")

    for index, value in enumerate(top):
        visit(value, f"_sx_arg_{index}")
    return top, leaves, unpacking, groups


def _expressions(exprs):
    import numpy as np
    from symox import DenseMatrix as Matrix, sympify, _core

    leaves = []

    def visit(value):
        if getattr(value, "_symox_matrixexpr", False):
            value = value.as_explicit()
        if isinstance(value, Matrix) or getattr(value, "_symox_sparse", False):
            return "dense", tuple(value.shape), tuple(visit(item) for item in value)
        if isinstance(value, np.ndarray):
            return "dense", tuple(value.shape), tuple(visit(item) for item in value.flat)
        if isinstance(value, (tuple, list)):
            return "tuple" if isinstance(value, tuple) else "list", tuple(visit(item) for item in value)
        expression = _core._materialize_matrix_scalars(sympify(value))
        index = len(leaves)
        leaves.append(expression)
        return "leaf", index

    return leaves, visit(exprs)


def _contains_conditionals(expressions):
    stack = list(expressions)
    seen = set()
    while stack:
        value = stack.pop()
        if value in seen:
            continue
        seen.add(value)
        if value.class_name == "Piecewise":
            return True
        stack.extend(value.args)
    return False


def _check_bound(expressions, symbols):
    free = set().union(*(expr.free_symbols for expr in expressions)) if expressions else set()
    unknown = free.difference(symbols)
    if unknown:
        raise ValueError("Expression contains unbound symbols: " + ", ".join(sorted(map(str, unknown))))


def _cse_preserving_assumptions(expressions, prefix):
    """Keep proven domains when CSE introduces a fresh temporary symbol."""
    from symox import Symbol, _core

    replacements, expressions = _core._cse(expressions, True, None, prefix)
    mapping, annotated = {}, []
    for symbol, value in replacements:
        value = value.xreplace(mapping)
        assumptions = {field: True for field in (
            "positive", "negative", "real", "integer", "rational", "nonzero",
            "nonnegative", "nonpositive", "finite",
        ) if getattr(value, "is_" + field, None) is True}
        replacement = Symbol(symbol.name, **assumptions)
        mapping[symbol] = replacement
        annotated.append((replacement, value))
    return annotated, [value.xreplace(mapping) for value in expressions]


def _functions(user_functions, allow_callables):
    result, imports, namespace = {}, set(), {}
    for name, value in (user_functions or {}).items():
        if not isinstance(name, str):
            name = getattr(name, "__name__", None)
        if not isinstance(name, str):
            raise TypeError("Custom function names must be strings or named callables")
        if callable(value) and allow_callables:
            alias = f"_sx_user_{len(namespace)}"
            result[name] = alias
            namespace[alias] = value
        elif isinstance(value, str) and _DOTTED_NAME.fullmatch(value) and not any(keyword.iskeyword(p) for p in value.split(".")):
            if "." not in value and value not in vars(builtins):
                raise ValueError("Custom source functions must use a qualified import path, e.g. 'torch.sin'")
            result[name] = value
            if "." in value:
                # Import the entire containing module so torch.special and
                # external packages both work in a fresh interpreter.
                imports.add(value.rsplit(".", 1)[0])
        else:
            raise TypeError("Custom functions must be qualified Python names" + (" or callables" if allow_callables else ""))
    return result, imports, namespace


def _generate(name, args, exprs, cse, user_functions, allow_callables=False):
    from symox import Symbol, _core

    if not isinstance(name, str) or not name.isidentifier() or keyword.iskeyword(name) or name.startswith("_sx_") or name == "torch" or name in vars(builtins):
        raise ValueError("Function name must be a Python identifier outside the reserved '_sx_' namespace")
    if not isinstance(cse, bool):
        raise TypeError("cse must be a bool")
    top, symbols, unpacking, groups = _arguments(args)
    expressions, structure = _expressions(exprs)
    # Dummify before CSE so two distinct symbols named x (with different
    # assumptions) stay independent and user names cannot collide with helpers.
    dummy_symbols = [Symbol(f"_sx_input_{index}", **symbol.assumptions0) for index, symbol in enumerate(symbols)]
    dummies = dict(zip(symbols, dummy_symbols))
    expressions = [expr.xreplace(dummies) for expr in expressions]
    _check_bound(expressions, dummy_symbols)
    functions, imports, namespace = _functions(user_functions, allow_callables)
    if name in {module.split(".")[0] for module in imports}:
        raise ValueError("Function name conflicts with a custom function module")
    replacements = []
    # A shared CSE temporary would otherwise evaluate an inactive Piecewise
    # branch before the printer can apply its domain mask.
    if cse and expressions and not _contains_conditionals(expressions):
        prefix = "_sx_cse"
        names = {symbol.name for symbol in dummy_symbols}
        while any(symbol.startswith(prefix) for symbol in names):
            prefix += "_"
        replacements, expressions = _cse_preserving_assumptions(expressions, prefix)
    names = {symbol: f"_sx_v_{index}" for index, symbol in enumerate(dummy_symbols)}
    names.update({symbol: f"_sx_t_{index}" for index, (symbol, _) in enumerate(replacements)})

    def emit(expr):
        # Tensorize the imaginary unit before arithmetic. Inductor's complex
        # fallback kernels cannot carry a Python complex scalar through backward.
        from symox import I
        return _core._torchcode(expr, symbol_names=names, user_functions=functions,
                                reference="_sx_ref",
                                expression_overrides={I: "_sx_tensor(1j, _sx_ref)"})

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
        if isinstance(length, tuple):
            source_path, shape = length
            lines.append(f"    {path} = _sx_matrix({source_path}, {shape!r})")
        else:
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


def torchcode_function(name, args, exprs, *, cse=True, user_functions=None):
    """Return a complete Python module implementing a PyTorch function.

    The generated function accepts ``dtype=`` and ``device=`` per call. Tensor
    inputs infer both; Python scalar inputs use torch's default floating dtype.
    Input symbols are safely renamed, including Python keywords and Unicode.
    CSE is shared across all outputs, and disabled when Piecewise is present.
    ``user_functions`` maps symbolic names to qualified Python function names.
    """
    return _generate(name, args, exprs, cse, user_functions)[0]


def torch_lambdify(args, exprs, *, cse=True, compile=False, compile_options=None, user_functions=None):
    """Build an inspectable, broadcasting, autograd-compatible torch callable.

    ``compile=True`` invokes ``torch.compile`` with ``compile_options`` and
    propagates compilation errors. Custom functions may be callables or dotted
    Python names. The source is available as ``function.__source__`` and through
    ``inspect.getsource(function)``.
    """
    if not isinstance(compile, bool):
        raise TypeError("compile must be a bool")
    if compile_options is not None and not isinstance(compile_options, dict):
        raise TypeError("compile_options must be a dict")
    if compile_options and not compile:
        raise ValueError("compile_options require compile=True")
    name = "_symox_torchgenerated"
    source, namespace = _generate(name, args, exprs, cse, user_functions, allow_callables=True)
    filename = f"<symox-torch-{next(_GENERATED)}>"
    code = builtins.compile(source, filename, "exec")
    exec(code, namespace)
    function = namespace[name]
    linecache.cache[filename] = (len(source), None, source.splitlines(keepends=True), filename)
    function.__source__ = source
    function.__doc__ = "Generated PyTorch function. Expression dimensions follow broadcast batch dimensions."
    if compile:
        import torch
        original = function
        function = torch.compile(original, **(compile_options or {}))
        function.__source__ = source
        function.__wrapped__ = original
    return function


__all__ = ["torchcode_function", "torch_lambdify"]
