"""Graph decompositions composed from native graph and matrix operations."""


def _permutation(matrix, components):
    from . import PermutationMatrix
    from .combinatorics import Permutation
    indices = [i for component in components for i in component]
    return PermutationMatrix(Permutation(indices))


def connected_components_decomposition(self):
    """Return P, B with self = P.T * B * P and B block diagonal.

    P is a native PermutationMatrix; B is a BlockDiagMatrix.
    Both factors support as_explicit(), transpose and native matrix algebra.
    """
    from . import BlockDiagMatrix
    components = self.connected_components()
    return (_permutation(self, components),
            BlockDiagMatrix(*(self.extract(c, c) for c in components)))


def strongly_connected_components_decomposition(self, lower=True):
    """Return P, B with self = P.T * B * P and B block triangular."""
    from . import BlockMatrix
    components = self.strongly_connected_components()
    if not lower:
        components = list(reversed(components))
    blocks = [[self.extract(a, b) for b in components] for a in components]
    return _permutation(self, components), BlockMatrix(blocks)


def to_DM(self, domain=None, **kwargs):
    """Convert to a native DomainMatrix, inferring its domain if omitted."""
    from .polys.matrices import DomainMatrix
    if domain is not None:
        if kwargs:
            raise TypeError('Options cannot be used with domain parameter')
        return DomainMatrix(self.todod(), self.shape, domain)
    return DomainMatrix.from_Matrix(self, **kwargs)


def install(matrix, sparse_base):
    for cls in (matrix, sparse_base):
        cls.is_Matrix = True
        cls.connected_components_decomposition = connected_components_decomposition
        cls.strongly_connected_components_decomposition = strongly_connected_components_decomposition
        cls.to_DM = to_DM
