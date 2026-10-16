"""Immutable code generation nodes understood by the C++ printer."""

from dataclasses import dataclass


def _expression(value):
    from symox import sympify
    return value if isinstance(value, (Indexed, IndexedBase)) else sympify(value)


@dataclass(frozen=True, init=False)
class Assignment:
    lhs: object
    rhs: object

    def __init__(self, lhs, rhs):
        object.__setattr__(self, "lhs", _expression(lhs))
        object.__setattr__(self, "rhs", _expression(rhs))

    @property
    def args(self):
        return self.lhs, self.rhs

    @property
    def free_symbols(self):
        return self.lhs.free_symbols | self.rhs.free_symbols


@dataclass(frozen=True, init=False)
class CodeBlock:
    args: tuple

    def __init__(self, *args):
        object.__setattr__(self, "args", tuple(args))

    @property
    def free_symbols(self):
        return set().union(*(getattr(arg, "free_symbols", set()) for arg in self.args))


@dataclass(frozen=True)
class Comment:
    text: str

    def __post_init__(self):
        object.__setattr__(self, "text", str(self.text))

    @property
    def args(self):
        return (self.text,)


class _IndexedArithmetic:
    """Allow indexed references in scalar symbolic expressions."""

    is_Symbol = False
    is_Number = False
    is_negative = None
    is_Pow = False
    is_Atom = True
    is_commutative = True

    def _symox_(self):
        from symox import Symbol
        return Symbol(str(self))

    def __add__(self, other): return self._symox_() + other
    def __radd__(self, other): return other + self._symox_()
    def __sub__(self, other): return self._symox_() - other
    def __rsub__(self, other): return other - self._symox_()
    def __mul__(self, other): return self._symox_() * other
    def __rmul__(self, other): return other * self._symox_()
    def __truediv__(self, other): return self._symox_() / other
    def __rtruediv__(self, other): return other / self._symox_()
    def __pow__(self, other): return self._symox_() ** other
    def __neg__(self): return -self._symox_()


@dataclass(frozen=True, init=False)
class IndexedBase(_IndexedArithmetic):
    label: object
    shape: object

    def __init__(self, label, shape=None, **assumptions):
        from symox import Symbol
        label = Symbol(label, **assumptions) if isinstance(label, str) else _expression(label)
        object.__setattr__(self, "label", label)
        object.__setattr__(self, "shape", None if shape is None else tuple(_expression(d) for d in shape))

    @property
    def args(self):
        return (self.label,) if self.shape is None else (self.label, self.shape)

    @property
    def free_symbols(self):
        return self.label.free_symbols

    def __getitem__(self, indices):
        return Indexed(self, *(indices if isinstance(indices, tuple) else (indices,)))

    def __str__(self):
        return str(self.label)


@dataclass(frozen=True, init=False)
class Indexed(_IndexedArithmetic):
    base: IndexedBase
    indices: tuple

    def __init__(self, base, *indices):
        if not indices:
            raise ValueError("Indexed needs at least one index")
        if not isinstance(base, IndexedBase):
            base = IndexedBase(base)
        if base.shape is not None and len(indices) != len(base.shape):
            raise ValueError("Rank mismatch")
        object.__setattr__(self, "base", base)
        object.__setattr__(self, "indices", tuple(_expression(i) for i in indices))

    @property
    def args(self):
        return (self.base,) + self.indices

    @property
    def free_symbols(self):
        return self.base.free_symbols | set().union(*(i.free_symbols for i in self.indices)) | {self}

    def __str__(self):
        return f"{self.base}[{', '.join(map(str, self.indices))}]"

    def _symox_(self):
        # Flatten multidimensional references consistently when they enter
        # native scalar arithmetic (whose atoms are scalar symbols).
        from symox import Symbol
        from symox.printing.cxx import CXX17CodePrinter
        return Symbol(CXX17CodePrinter().doprint(self))
