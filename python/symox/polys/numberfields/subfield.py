"""Native primitive elements with exact constituent embeddings."""
from ..._core import primitive_element, to_number_field

__all__ = ['primitive_element', 'to_number_field']


def field_isomorphism(a,b,*,fast=True):
    """Return an exact embedding in falling powers of the target generator."""
    from ... import AlgebraicNumber, IsomorphismFailed, sympify
    a,b=sympify(a),sympify(b)
    a=a.to_primitive_element() if getattr(a,"is_AlgebraicNumber",False) else AlgebraicNumber(a)
    b=b.to_primitive_element() if getattr(b,"is_AlgebraicNumber",False) else AlgebraicNumber(b)
    if b.minpoly.degree()%a.minpoly.degree():return None
    try:
        return to_number_field(a.root,b.root).coeffs()
    except IsomorphismFailed:
        return None


__all__.append('field_isomorphism')
