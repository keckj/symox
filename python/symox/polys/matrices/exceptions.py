"""Domain-matrix error names, independent of an optional SymPy installation."""


class DMError(Exception):
    """Base class for errors in exact matrix-domain operations."""


class DMBadInputError(DMError):
    pass


class DMDomainError(DMError):
    pass


class DMShapeError(DMError):
    pass


class DMValueError(DMError):
    pass


class DMRankError(DMError):
    pass


class DMNotAField(DMDomainError):
    pass


class DMNonSquareMatrixError(DMShapeError):
    pass


class DMNonInvertibleMatrixError(DMError):
    pass
