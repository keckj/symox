"""Seed-independent serialization of stock code for reload authentication."""
from types import CodeType

_TYPE=type
_LEN=len
_SORTED=sorted
_REPLACE=CodeType.replace
_FROZENSET=frozenset
_TUPLE=tuple
_LIMIT=64


def code_key(value,dumps,depth=0):
    """Return bytes identifying a code object, independent of the hash seed.

    ``dumps`` is the caller's authenticated ``marshal.dumps``; callers keep
    checking that alias before each use.

    Stock CPython 3.10 marshal writes a frozenset constant (``x in {'a',
    'b'}`` or a set display of constants) in its iteration order, which
    depends on the hash seed and on the insertion history. A module compiled
    from source and the same module reloaded from a .pyc written by another
    process can therefore marshal differently although their code is
    identical. Here frozenset items are serialized in sorted order of their
    own keys.

    Only exact code, tuple and frozenset containers are walked. Every other
    constant goes through marshal, which rejects foreign objects with
    ValueError and never calls their equality, hashing or reduction.
    """
    if depth>_LIMIT:raise ValueError('code constants are nested too deeply')
    kind=_TYPE(value)
    if kind is CodeType:
        tag=b'c';parts=[dumps(_REPLACE(value,co_consts=()),0)]
        parts+=[code_key(item,dumps,depth+1) for item in value.co_consts]
    elif kind is _FROZENSET:
        tag=b'f';parts=_SORTED(code_key(item,dumps,depth+1) for item in value)
    elif kind is _TUPLE:
        tag=b't';parts=[code_key(item,dumps,depth+1) for item in value]
    else:
        return b'm'+dumps(value,0)
    return tag+b''.join(_LEN(part).to_bytes(8,'little')+part for part in parts)
