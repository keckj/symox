"""Admission of finite Markov chain computations to the native kernels
(symox::stats::markov, binding crates/symox-py/src/markov.rs).

Only stock dense matrices of native expressions are admitted. Every function
returns None when the inputs are outside the admitted domain or the kernel
declines (undecidable step, non-unique solution, inexact linear algebra);
the caller then keeps its original path."""
import symox as S

_STOCK = (S.ImmutableDenseMatrix, S.MutableDenseMatrix)


def _stock(T):
    return type(T) in _STOCK


def communication_classes(T):
    return S._core._markov_communication_classes(T) if _stock(T) else None


def decompose(T):
    return S._core._markov_decompose(T) if _stock(T) else None


def canonical_form(T):
    return S._core._markov_canonical_form(T) if _stock(T) else None


def stationary_distribution(T):
    return S._core._markov_stationary(T) if _stock(T) else None


def ctmc_limiting_distribution(G):
    return S._core._markov_ctmc_limiting(G) if _stock(G) else None


def fundamental_matrix(T):
    """Raises the reference ValueError for a singular ``I - Q``."""
    return S._core._markov_fundamental_matrix(T) if _stock(T) else None


def absorbing_probabilities(T):
    return S._core._markov_absorbing_probabilities(T) if _stock(T) else None


def path_sums(kernels, initial, coded, labels, atoms, gate, event, expectation):
    """``(numerator, denominator)`` of the path enumeration of a query."""
    if not all(_stock(k) for k in kernels):
        return None
    # The random atoms are host Symbol subclasses: the kernel sees plain
    # dummies in their place (every atom is substituted on every path).
    plain = {atom: S.Dummy() for atom, _ in atoms}
    gate = None if gate is None else gate.xreplace(plain)
    return S._core._markov_path_sums(list(kernels), initial, list(coded), list(labels),
                                     [(plain[a], i) for a, i in atoms], gate,
                                     event.xreplace(plain), expectation)
