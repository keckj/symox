"""Admission for the exact Smith decomposition in Rust; no arithmetic here.

The admitted computation is the one of ``_smith_normal_decomp`` for an exact
``symox._core.Domain`` that the kernel admits (ZZ with stock ZZ elements;
QQ, GF(p) with p prime, QQ<a>, QQ_I, rational function fields and the
univariate rings QQ[x] with native elements of exactly that domain, checked in
Rust).
At registration the transitive protocol actually read by the reference body
is recorded: every global/builtin name and every class attribute (of Domain,
DomainElement, DomainMatrix, DDM, SDM) named by the reference code and by the
Python functions it reaches, plus the implicit operator slots. Each call
proves by identity that none of those bindings, nor any reached function's
code, globals, defaults, keyword defaults or closure, has changed (function
fields are read by ``_core._ddm_function_metadata``, without Python attribute
access), with plain dict lookups done in Rust. Any change, unverified platform or other domain runs the original
Python algorithm unchanged; a native refusal does the same.
"""
from types import ModuleType as _MODULE
from symox import _core

_MISSING = object()
_TYPE = type
_FUNCTION = (lambda: None).__class__
_CODE = (lambda: None).__code__.__class__
_TUPLE, _DICT, _STR = tuple, dict, str
_STATIC, _CLASS, _PROPERTY = staticmethod, classmethod, property
_READ = _core.__dict__.get('_ddm_function_metadata')
_CHECK = _core.__dict__.get('_snapshot_unchanged')


def _code_names(code, names):
    names.update(code.co_names)
    for constant in code.co_consts:
        if _TYPE(constant) is _CODE:
            _code_names(constant, names)


class Snapshot:
    """Identity proof of a transitive Python protocol captured at registration."""

    __slots__ = ('entries', 'functions', 'classes', 'ok')

    def __init__(self, roots, classes=(), implicit=()):
        self.entries = []
        self.functions = []
        self.classes = []
        self.ok = False
        if _TYPE(_READ) is not [].append.__class__:
            return
        names = set(implicit)
        visited = {}
        pending = list(roots)
        recorded = set()
        scanned = set()
        proxies = []
        if _TYPE(_CHECK) is not [].append.__class__:
            return
        import gc
        for cls in classes:
            if not isinstance(cls, _TYPE):
                return
            proxy = cls.__dict__
            # The class namespace itself (the mapping behind the read-only
            # proxy), so the per-call proof uses plain dict lookups.
            referents = gc.get_referents(proxy)
            if len(referents) != 1 or _TYPE(referents[0]) is not _DICT:
                return
            namespace = referents[0]
            if len(namespace) != len(proxy) or any(namespace[key] is not value for key, value in proxy.items()):
                return
            proxies.append(namespace)
            self.classes.append((cls, cls.__mro__, _TYPE(cls), namespace, len(namespace)))
        while True:
            while pending:
                function = pending.pop()
                if id(function) in visited:
                    continue
                metadata = _READ(function)
                if metadata is None or _TYPE(metadata) is not _TUPLE or len(metadata) != 5:
                    return
                visited[id(function)] = function
                self.functions.append((function, metadata))
                code, namespace = metadata[0], metadata[1]
                own = set()
                _code_names(code, own)
                names.update(own)
                if _TYPE(namespace) is not _DICT:
                    return
                builtins = namespace.get('__builtins__', _MISSING)
                self._record(namespace, '__builtins__', builtins, recorded, pending)
                if _TYPE(builtins) is _MODULE:
                    builtins = builtins.__dict__
                if _TYPE(builtins) is not _DICT:
                    return
                for name in own:
                    value = namespace.get(name, _MISSING)
                    self._record(namespace, name, value, recorded, pending)
                    if value is _MISSING:
                        self._record(builtins, name, builtins.get(name, _MISSING), recorded, pending)
                for cell in metadata[4] or ():
                    try:
                        value = cell.cell_contents
                    except ValueError:
                        continue
                    if _TYPE(value) is _FUNCTION:
                        pending.append(value)
            # Each name is looked up once in every class namespace: the
            # namespaces do not change during the capture, so a name already
            # scanned in an earlier pass yields nothing new.
            fresh = names - scanned
            scanned.update(fresh)
            for proxy in proxies:
                key = id(proxy)
                for name in fresh:
                    # Absent class names are covered by the namespace size
                    # check below: a binding added to a class changes it.
                    if name in proxy and (key, name) not in recorded:
                        self._record(proxy, name, proxy[name], recorded, pending)
            if not pending:
                break
        self.entries = _TUPLE(self.entries)
        self.functions = _TUPLE(self.functions)
        for container, name, value in self.entries:
            if _TYPE(container) is not _DICT or _TYPE(name) is not _STR:
                return
        self.ok = True

    def _record(self, container, name, value, recorded, pending):
        key = (id(container), name)
        if key in recorded:
            return
        recorded.add(key)
        self.entries.append((container, name, value))
        kind = _TYPE(value)
        if kind is _FUNCTION:
            pending.append(value)
        elif kind is _STATIC or kind is _CLASS:
            if _TYPE(value.__func__) is _FUNCTION:
                pending.append(value.__func__)
        elif kind is _PROPERTY:
            for accessor in (value.fget, value.fset, value.fdel):
                if _TYPE(accessor) is _FUNCTION:
                    pending.append(accessor)

    def unchanged(self):
        if not self.ok or _core.__dict__.get('_ddm_function_metadata') is not _READ:
            return False
        for cls, mro, meta, namespace, size in self.classes:
            if cls.__mro__ is not mro or _TYPE(cls) is not meta or len(namespace) != size:
                return False
        return _core.__dict__.get('_snapshot_unchanged') is _CHECK and _CHECK(self.entries, self.functions, _MISSING)


_ELEMENT_SLOTS = ('__eq__', '__ne__', '__mul__', '__rmul__', '__imul__', '__add__', '__radd__',
                  '__iadd__', '__sub__', '__rsub__', '__isub__', '__neg__', '__bool__', '__hash__',
                  '__index__', '__int__', '__getattribute__', '__getattr__', '__class__')
_MATRIX_SLOTS = ('__new__', '__init__', '__mul__', '__rmul__', '__matmul__', '__getattribute__',
                 '__getattr__', '__class__', 'matmul', 'to_list', 'unify', 'rep')

_STATE = None
_PENDING = None


def register(namespace):
    """Record the reference Smith module; the protocol is captured at first use.

    Package initialisation still installs element/domain methods after this
    module is imported, so the snapshot is taken lazily, once, from the
    completed runtime (the reference functions' metadata is fixed here).
    """
    global _PENDING
    if _PENDING is not None:
        return
    roots = [namespace.get('_smith_normal_decomp'), namespace.get('add_columns')]
    if _TYPE(_READ) is not [].append.__class__ or any(_TYPE(root) is not _FUNCTION for root in roots):
        _PENDING = False
        return
    _PENDING = (namespace, roots, [_READ(root) for root in roots])


def _capture():
    global _STATE
    _STATE = False
    if not _PENDING:
        return
    namespace, roots, metadata = _PENDING
    for root, fields in zip(roots, metadata):
        current = _READ(root)
        if fields is None or current is None or any(a is not b for a, b in zip(current, fields)):
            return
    from .domainmatrix import DomainMatrix
    from .ddm import DDM
    from .sdm import SDM
    import sys
    snapshot = Snapshot(roots, (_core.Domain, _core.DomainElement, DomainMatrix, DDM, SDM),
                        _ELEMENT_SLOTS + _MATRIX_SLOTS)
    if snapshot.ok:
        _STATE = (snapshot, namespace, roots[0], _core.ZZ, _core.__dict__.get('_ddm_smith'),
                  sys.getrecursionlimit)


def try_smith(rows, domain, shape, full):
    """Native ``_smith_normal_decomp(rows, domain, shape, full)``, or None."""
    if _STATE is None:
        _capture()
    state = _STATE
    if not state:
        return None
    snapshot, namespace, reference, zz, kernel, limit = state
    if _TYPE(domain) is not _core.Domain:
        return None
    if (namespace.get('_smith_normal_decomp') is not reference
            or _core.__dict__.get('_ddm_smith') is not kernel or not snapshot.unchanged()):
        return None
    if _TYPE(shape) is not _TUPLE or len(shape) != 2:
        return None
    # The reference recursion depth is min(shape); leave its RecursionError
    # (or success) to the reference when it could approach the limit.
    if min(shape) * 4 > limit():
        return None
    if domain is zz:
        return kernel(rows, shape, full)
    return kernel(rows, shape, full, domain)
