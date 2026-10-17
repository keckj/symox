"""Assumption-dependent scalar rewrites."""
from .. import _core
from .assume import global_assumptions


def refine(expr, assumptions=True):
    """Rewrite expressions using supplied and global assumptions."""
    return _core._refine(expr, assumptions, global_assumptions)
