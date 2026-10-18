"""Matrices over native coefficient domains."""
from .exceptions import (DMError, DMBadInputError, DMDomainError, DMShapeError,
                         DMValueError, DMRankError, DMNotAField,
                         DMNonSquareMatrixError, DMNonInvertibleMatrixError)
from .domainmatrix import DomainMatrix
from .domainscalar import DomainScalar


def DM(rows, domain):
    return DomainMatrix.from_list(rows, domain)


__all__ = ['DomainMatrix', 'DomainScalar', 'DM', 'DMError', 'DMBadInputError',
           'DMDomainError', 'DMShapeError', 'DMValueError', 'DMRankError',
           'DMNotAField', 'DMNonSquareMatrixError', 'DMNonInvertibleMatrixError']
