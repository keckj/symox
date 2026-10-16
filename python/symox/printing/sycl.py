"""SYCL 2020 expression printing for AdaptiveCpp and other SYCL compilers.

Generated snippets require ``<sycl/sycl.hpp>`` and ``<limits>``. Their inputs
are real scalars; ``float64`` is the default, with ``float32`` available.
Device execution of doubles requires the SYCL ``fp64`` aspect.
"""

from .slang import SlangCodePrinter, _normalize_slang_dtype


def _is_sycl_identifier(name):
    from symox import _core
    return isinstance(name, str) and _core._is_sycl_identifier(name)


def _normalize_sycl_dtype(dtype):
    try:
        return _normalize_slang_dtype(dtype)
    except ValueError:
        raise ValueError("SYCL dtype must be float32 or float64") from None


class SYCLCodePrinter(SlangCodePrinter):
    """Extensible native SYCL expression printer, suitable for AdaptiveCpp.

    Accepts ``dtype``, ``symbol_names`` and ``user_functions`` settings.
    ``loggamma`` requires a provably positive real argument: ``sycl::lgamma``
    computes log(abs(gamma)) on negative inputs. Custom function targets can
    be namespace-qualified C++ names, such as ``optics::distort``.
    """

    language = "SYCL"
    printmethod = "_syclcode"
    _native_method = "_syclcode"
    _default_settings = {"order": None, "dtype": "float64", "symbol_names": {}, "user_functions": {}}

    def __init__(self, settings=None, **kwargs):
        super().__init__(settings, **kwargs)
        self._custom_methods = {
            name: getattr(type(self), name)
            for name in dir(type(self))
            if name.startswith("_print_")
            and getattr(type(self), name) is not getattr(SYCLCodePrinter, name, None)
        }

    def _is_identifier(self, name):
        return _is_sycl_identifier(name)


SyclCodePrinter = SYCLCodePrinter
AdaptiveCppCodePrinter = SYCLCodePrinter


def syclcode(expr, assign_to=None, **settings):
    """Return a scalar SYCL expression using real floating point arithmetic."""
    return SYCLCodePrinter(settings).doprint(expr, assign_to)


adaptivecppcode = syclcode

__all__ = ["SYCLCodePrinter", "SyclCodePrinter", "AdaptiveCppCodePrinter", "syclcode", "adaptivecppcode"]
