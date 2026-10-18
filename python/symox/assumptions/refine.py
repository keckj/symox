"""Assumption-dependent scalar rewrites."""
from .. import _core
from .assume import global_assumptions


def refine(expr, assumptions=True):
    """Rewrite expressions using supplied and global assumptions."""
    return _core._refine(expr, assumptions, global_assumptions)


def _handler(name):
    def operation(expr,assumptions):
        result=refine(expr,assumptions)
        return result if result!=expr else None
    operation.__name__=name;return operation
for _name in 'Pow abs arg atan2 im matrixelement re sign'.split():globals()['refine_'+_name]=_handler('refine_'+_name)
