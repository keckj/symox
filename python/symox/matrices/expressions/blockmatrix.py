"""Native block diagonal matrices and their typed block grids."""
from ..._matrixexpr import BlockDiagMatrix, BlockMatrix, block_collapse

__all__ = ["BlockDiagMatrix", "BlockMatrix", "block_collapse"]

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'blockcut')
