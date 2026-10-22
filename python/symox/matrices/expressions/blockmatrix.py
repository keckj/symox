"""Native block diagonal matrices and their typed block grids."""
from ..._matrixexpr import BlockDiagMatrix, BlockMatrix, block_collapse
from ..._kind import matrix_tag

__all__ = ["BlockDiagMatrix", "BlockMatrix", "block_collapse"]

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox', 'blockcut')


def bounds(sizes):
    lower=0;result=[]
    for size in sizes:result.append((lower,lower+size));lower+=size
    return result


def bc_unpack(expr):
    return expr.blocks[0,0] if isinstance(expr,BlockMatrix) and expr.blockshape==(1,1) else expr


def _collapse_kind(expr,kinds):
    return block_collapse(expr) if matrix_tag(expr) in kinds else expr

def bc_matadd(expr):return _collapse_kind(expr,('MatAdd',))
def bc_block_plus_ident(expr):return _collapse_kind(expr,('MatAdd',))
def bc_dist(expr):return _collapse_kind(expr,('MatMul',))
def bc_matmul(expr):return _collapse_kind(expr,('MatMul','MatPow'))
def bc_transpose(expr):return _collapse_kind(expr,('Transpose',))


def blockinverse_1x1(expr):
    if getattr(expr,'_kind_tag',None)=='Inverse' and isinstance(expr.args[0],BlockMatrix) and expr.args[0].blockshape==(1,1):
        return BlockMatrix([[expr.args[0].blocks[0,0].inv()]])
    return expr


def blockinverse_2x2(expr):
    if getattr(expr,'_kind_tag',None)!='Inverse' or not isinstance(expr.args[0],BlockMatrix) or expr.args[0].blockshape!=(2,2):return expr
    from symox import ShapeError,NonInvertibleMatrixError
    grid=expr.args[0].blocks;A,B,C,D=[grid[i,j] for i in range(2) for j in range(2)]
    try:
        ai=A.inv();schur=(D-C*ai*B).inv()
        return BlockMatrix([[ai+ai*B*schur*C*ai,-ai*B*schur],[-schur*C*ai,schur]])
    except (ValueError,ShapeError,NonInvertibleMatrixError):
        try:
            di=D.inv();schur=(A-B*di*C).inv()
            return BlockMatrix([[schur,-schur*B*di],[-di*C*schur,di+di*C*schur*B*di]])
        except (ValueError,ShapeError,NonInvertibleMatrixError):return expr


def bc_inverse(expr):
    result=blockinverse_1x1(expr)
    if result!=expr:return result
    return blockinverse_2x2(expr)


def deblock(B):
    if not isinstance(B,BlockMatrix):return B
    grid=B.blocks;rows=[]
    for i in range(grid.rows):
        items=[grid[i,j] for j in range(grid.cols)]
        nested=[item if isinstance(item,BlockMatrix) else BlockMatrix([[item]]) for item in items]
        if len({item.blockshape[0] for item in nested})!=1:return B
        for k in range(nested[0].blockshape[0]):
            rows.append([item.blocks[k,j] for item in nested for j in range(item.blockshape[1])])
    try:return BlockMatrix(rows)
    except (ValueError,TypeError):return B


def reblock_2x2(expr):
    if not isinstance(expr,BlockMatrix) or not all(n>2 for n in expr.blockshape):return expr
    grid=expr.blocks;r,c=grid.shape
    def region(i0,i1,j0,j1):return bc_unpack(BlockMatrix([[grid[i,j] for j in range(j0,j1)] for i in range(i0,i1)]))
    return BlockMatrix([[region(0,1,0,1),region(0,1,1,c)],[region(1,r,0,1),region(1,r,1,c)]])

__all__ += ['bounds','bc_unpack','bc_matadd','bc_block_plus_ident','bc_dist','bc_matmul','bc_transpose',
 'bc_inverse','blockinverse_1x1','blockinverse_2x2','deblock','reblock_2x2']
