"""Standalone C++17 / SYCL 2020 functions and asynchronous USM kernels.

The generated source uses standard SYCL 2020 APIs and can be compiled by
AdaptiveCpp or another conforming SYCL compiler. No Python execution runtime is
required by source export.
"""

from __future__ import annotations

import re

from .pytorch import _arguments, _contains_conditionals, _cse_preserving_assumptions, _expressions


_HEADER = """#include <sycl/sycl.hpp>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
"""
_DTYPES = {"float32": "float", "float": "float", "float64": "double", "double": "double"}


def _generate(name, args, exprs, dtype, cse, user_functions):
    from symox import Symbol, _core
    from symox.extra.printing.sycl import _is_sycl_identifier

    if not _is_sycl_identifier(name) or name.startswith("_sx_"):
        raise ValueError("Function name must be a C++ identifier outside the reserved '_sx_' namespace")
    try:
        scalar_type = _DTYPES[dtype]
    except (KeyError, TypeError):
        raise ValueError("SYCL dtype must be 'float32' or 'float64'") from None
    if not isinstance(cse, bool):
        raise TypeError("cse must be a bool")
    _, symbols, _, _ = _arguments(args, backend="SYCL")
    expressions, structure = _expressions(exprs)
    free = set().union(*(expression.free_symbols for expression in expressions)) if expressions else set()
    unknown = free.difference(symbols)
    if unknown:
        raise ValueError("Expression contains unbound symbols: " + ", ".join(sorted(map(str, unknown))))
    if user_functions is not None and not isinstance(user_functions, dict):
        raise TypeError("user_functions must be a dict mapping symbolic names to C++ function names")
    functions = {}
    for key, value in (user_functions or {}).items():
        key = key if isinstance(key, str) else getattr(key, "__name__", None)
        if not isinstance(key, str) or not isinstance(value, str) or not re.fullmatch(
            r"(?:::)?[A-Za-z_][A-Za-z_0-9]*(?:::[A-Za-z_][A-Za-z_0-9]*)*", value
        ):
            raise ValueError("Custom functions must map names to C++ identifiers or namespace-qualified names")
        functions[key] = value
    dummy_symbols = [Symbol(f"_sx_input_{index}", **symbol.assumptions0) for index, symbol in enumerate(symbols)]
    expressions = [expression.xreplace(dict(zip(symbols, dummy_symbols))) for expression in expressions]
    # Store result types before CSE: a fresh Symbol does not remember that the
    # expression it replaced was boolean.
    output_types = ["bool" if _core._sycl_type(expression) == "bool" else scalar_type for expression in expressions]
    replacements = []
    if cse and expressions and not _contains_conditionals(expressions):
        replacements, expressions = _cse_preserving_assumptions(expressions, "_sx_cse")
    names = {symbol: f"_sx_v_{index}" for index, symbol in enumerate(dummy_symbols)}
    names.update({symbol: f"_sx_t_{index}" for index, (symbol, _) in enumerate(replacements)})
    temporary_types = {}

    def type_of(expression):
        return temporary_types.get(expression, "bool" if _core._sycl_type(expression) == "bool" else scalar_type)

    def emit(expression):
        return _core._syclcode(expression, dtype=dtype, symbol_names=names, user_functions=functions)

    scalar_output = structure[0] == "leaf"
    result_type = output_types[0] if scalar_output else f"_sx_result_{name}"
    lines = [_HEADER.rstrip(), "", "// Result fields follow depth-first row-major expression order."]
    if not scalar_output:
        lines.append(f"struct {result_type}")
        lines.append("{")
        for index, output_type in enumerate(output_types):
            lines.append(f"    {output_type} _sx_o_{index};")
        lines.extend(["};", ""])
    parameters = ", ".join(f"{scalar_type} {names[symbol]}" for symbol in dummy_symbols)
    lines.extend([f"inline {result_type} {name}({parameters})", "{"])
    for symbol, expression in replacements:
        temporary_types[symbol] = type_of(expression)
        lines.append(f"    {temporary_types[symbol]} {names[symbol]} = {emit(expression)};")
    rendered = [emit(expression) for expression in expressions]
    if scalar_output:
        lines.append(f"    return {rendered[0]};")
    else:
        lines.append("    return {" + ", ".join(rendered) + "};")
    lines.append("}")
    return "\n".join(lines) + "\n", symbols, structure, output_types, scalar_type, result_type


def syclcode_function(name, args, exprs, *, dtype="float64", cse=True, user_functions=None):
    """Return a standalone inline C++17 / SYCL 2020 scalar function.

    Arguments are safely renamed ``_sx_v_0``, ... in declared order. Scalar
    expressions return a scalar; containers, matrices and NumPy arrays return a
    struct with ``_sx_o_0``, ... fields in depth-first row-major order. Custom
    function mappings name existing C++ functions declared in the surrounding
    source. A SYCL device must support ``aspect::fp64`` to execute double code.
    """
    return _generate(name, args, exprs, dtype, cse, user_functions)[0]


def sycl_kernel(name, args, exprs, *, dtype="float64", cse=True, user_functions=None):
    """Export scalar ``name_eval`` and asynchronous USM kernel launcher ``name``.

    The launcher accepts a ``sycl::queue&``, a scalar element count, one constant
    USM input pointer per symbolic argument, one writable USM output pointer per
    expression leaf, and optional event dependencies. Buffers contain contiguous
    equal-length batches. The caller owns allocations and their lifetime, must
    use USM pointers valid in the queue's context, and must avoid overlapping
    input/output buffers. No allocation or wait is performed by the launcher.

    The returned event covers dependency execution even when the count is zero.
    Float64 execution checks the queue device's ``aspect::fp64`` before launch.
    """
    from symox.extra.printing.sycl import _is_sycl_identifier

    if not _is_sycl_identifier(name) or name.startswith("_sx_"):
        raise ValueError("Kernel name must be a C++ identifier outside the reserved '_sx_' namespace")
    source, symbols, structure, output_types, scalar_type, result_type = _generate(
        name + "_eval", args, exprs, dtype, cse, user_functions
    )
    parameters = ["    sycl::queue& _sx_queue", "    std::size_t _sx_count"]
    parameters += [f"    const {scalar_type}* _sx_in_{index}" for index in range(len(symbols))]
    parameters += [f"    {output_type}* _sx_out_{index}" for index, output_type in enumerate(output_types)]
    parameters.append("    const std::vector<sycl::event>& _sx_dependencies = {}")
    lines = [source.rstrip(), "", f"inline sycl::event {name}(", ",\n".join(parameters), ")", "{"]
    if scalar_type == "double":
        lines.extend(["    if (!_sx_queue.get_device().has(sycl::aspect::fp64))", '        throw std::runtime_error("Generated SYCL float64 kernel requires aspect::fp64");'])
    lines.extend([
        "    return _sx_queue.submit([&](sycl::handler& _sx_handler)",
        "    {",
        "        _sx_handler.depends_on(_sx_dependencies);",
        "        if (_sx_count == 0)",
        "        {",
        "            _sx_handler.single_task([]() {});",
        "            return;",
        "        }",
        "        _sx_handler.parallel_for(sycl::range<1>(_sx_count), [=](sycl::id<1> _sx_id)",
        "        {",
        "            const std::size_t _sx_index = _sx_id[0];",
    ])
    inputs = ", ".join(f"_sx_in_{index}[_sx_index]" for index in range(len(symbols)))
    lines.append(f"            const {result_type} _sx_value = {name}_eval({inputs});")
    for index in range(len(output_types)):
        value = "_sx_value" if structure[0] == "leaf" else f"_sx_value._sx_o_{index}"
        lines.append(f"            _sx_out_{index}[_sx_index] = {value};")
    lines.extend(["        });", "    });", "}"])
    return "\n".join(lines) + "\n"


adaptivecppcode_function = syclcode_function
adaptivecpp_kernel = sycl_kernel

__all__ = ["syclcode_function", "sycl_kernel", "adaptivecppcode_function", "adaptivecpp_kernel"]
