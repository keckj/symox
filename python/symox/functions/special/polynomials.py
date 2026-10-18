"""Compatible imports backed by existing native implementations."""
from symox._api_exports import export as _export
_export(globals(), 'symox', 'assoc_laguerre assoc_legendre chebyshevt chebyshevt_root chebyshevu chebyshevu_root gegenbauer hermite hermite_prob jacobi jacobi_normalized laguerre legendre')

from ...core.function import DefinedFunction

class OrthogonalPolynomial(DefinedFunction):
    _matches=staticmethod(lambda expr:expr.class_name in ('assoc_laguerre','assoc_legendre','chebyshevt','chebyshevu',
        'gegenbauer','hermite','hermite_prob','jacobi','laguerre','legendre'))
    @classmethod
    def _eval_at_order(cls,n,x):
        from symox import sympify,Dummy
        n=sympify(n)
        if n.is_Integer and n>=0:
            symbol=Dummy('orthogonal_x');value=cls._ortho_poly(int(n),symbol)
            return (value.as_expr() if hasattr(value,'as_expr') else value).subs(symbol,x)
__all__.append('OrthogonalPolynomial')
