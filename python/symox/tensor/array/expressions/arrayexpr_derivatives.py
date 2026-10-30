"""``sympy.tensor.array.expressions.arrayexpr_derivatives`` import path.

``array_derive`` runs in Rust (``_native``); ``matrix_derive`` composes the
native ``convert_matrix_to_array`` and ``array_derive`` with
``convert_array_to_matrix`` and ``doit`` as SymPy does.
"""
from ._native import array_derive


def matrix_derive(expr, x):
    from . import convert_array_to_matrix
    from ._canon import convert_matrix_to_array
    ce = convert_matrix_to_array(expr)
    dce = array_derive(ce, x)
    return convert_array_to_matrix(dce).doit()


__all__ = ['array_derive', 'matrix_derive']
