"""Exceptions for the native polynomial and number-field APIs.

The class names and the hierarchy of the Python-defined classes follow
``sympy.polys.polyerrors``. ``PolynomialError`` and ``GeneratorsError`` also
derive from ``ValueError`` because native code historically raised plain
``ValueError`` for these conditions. The native exceptions imported from
``_core`` (``DomainError``, ``CoercionFailed``, ...) are ``ValueError``
subclasses that do not derive from ``BasePolynomialError``.
"""
from .._core import IsomorphismFailed, NotInvertible, DomainError, CoercionFailed, ExactQuotientFailed, OptionError, FlagError, UnificationFailed, RefinementFailed


class BasePolynomialError(Exception):
    """Base class for polynomial related exceptions."""

    def new(self, *args):
        raise NotImplementedError("abstract base class")


class PolynomialError(BasePolynomialError, ValueError):
    """Raised for invalid polynomial operations."""


class PolynomialDivisionFailed(BasePolynomialError):

    def __init__(self, f, g, domain):
        self.f = f
        self.g = g
        self.domain = domain
        super().__init__(f, g, domain)

    def __str__(self):
        if getattr(self.domain, "is_EX", False):
            msg = "You may want to use a different simplification algorithm. Note " \
                  "that in general it's not possible to guarantee to detect zero " \
                  "in this domain."
        elif not getattr(self.domain, "is_Exact", True):
            msg = "Your working precision or tolerance of computations may be set " \
                  "improperly. Adjust those parameters of the coefficient domain " \
                  "and try again."
        else:
            msg = "Zero detection is guaranteed in this coefficient domain. This " \
                  "may indicate a bug in SymPy or the domain is user defined and " \
                  "doesn't implement zero detection properly."
        return "couldn't reduce degree in a polynomial division algorithm when " \
               "dividing %s by %s. This can happen when it's not possible to " \
               "detect zero in the coefficient domain. The domain of computation " \
               "is %s. %s" % (self.f, self.g, self.domain, msg)


class OperationNotSupported(BasePolynomialError):

    def __init__(self, poly, func):
        self.poly = poly
        self.func = func
        super().__init__(poly, func)

    def __str__(self):
        return "`%s` operation not supported by %s representation" % (
            self.func, type(self.poly).__name__)


class HeuristicGCDFailed(BasePolynomialError):
    pass


class ModularGCDFailed(BasePolynomialError):
    pass


class HomomorphismFailed(BasePolynomialError):
    pass


class ExtraneousFactors(BasePolynomialError):
    pass


class EvaluationFailed(BasePolynomialError):
    pass


class NotReversible(BasePolynomialError):
    pass


class NotAlgebraic(BasePolynomialError):
    pass


class UnsolvableFactorError(BasePolynomialError):
    """Raised if a ``roots`` call with ``strict=True`` cannot find all roots."""


class GeneratorsError(BasePolynomialError, ValueError):
    """Raised when an operation refers to an unknown polynomial generator."""


class GeneratorsNeeded(GeneratorsError):
    pass


class ComputationFailed(BasePolynomialError):

    def __init__(self, func, nargs, exc):
        self.func = func
        self.nargs = nargs
        self.exc = exc
        super().__init__(func, nargs, exc)

    def __str__(self):
        return "%s(%s) failed without generators" % (
            self.func, ', '.join(map(str, self.exc.exprs[:self.nargs])))


class UnivariatePolynomialError(PolynomialError):
    pass


class MultivariatePolynomialError(PolynomialError):
    pass


class PolificationFailed(PolynomialError):

    def __init__(self, opt, origs, exprs, seq=False):
        if not seq:
            self.orig = origs
            self.expr = exprs
            self.origs = [origs]
            self.exprs = [exprs]
        else:
            self.origs = origs
            self.exprs = exprs
        self.opt = opt
        self.seq = seq
        super().__init__(opt, origs, exprs, seq)

    def __str__(self):
        if not self.seq:
            return "Cannot construct a polynomial from %s" % str(self.orig)
        return "Cannot construct polynomials from %s" % ', '.join(map(str, self.origs))


__all__ = [
    "BasePolynomialError", "ExactQuotientFailed", "PolynomialDivisionFailed",
    "OperationNotSupported", "HeuristicGCDFailed", "ModularGCDFailed",
    "HomomorphismFailed", "IsomorphismFailed", "ExtraneousFactors",
    "EvaluationFailed", "RefinementFailed", "CoercionFailed", "NotInvertible",
    "NotReversible", "NotAlgebraic", "DomainError", "PolynomialError",
    "UnificationFailed", "UnsolvableFactorError", "GeneratorsError",
    "GeneratorsNeeded", "ComputationFailed", "UnivariatePolynomialError",
    "MultivariatePolynomialError", "PolificationFailed", "OptionError",
    "FlagError",
]
