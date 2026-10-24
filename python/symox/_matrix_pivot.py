"""Host side of the native pivot search.

SymPy's ``_find_reasonable_pivot`` runs the full ``simplify`` on the
candidates its zero test leaves undecided, then ``equals(0)``.  The native
pivot search (``pivot_simp``) only expands and cancels, so it asks this hook
about entries that are not rational functions of their symbols.  The result
is returned only when it decides the zero test; errors count as undecided.
"""

_active = False


def pivot_simplify(value):
    global _active
    if _active:
        return None
    _active = True
    try:
        from . import S, simplify
        try:
            simplified = simplify(value)
        except Exception:
            simplified = value
        zero = getattr(simplified, 'is_zero', None)
        if zero is not None:
            return simplified
        try:
            if value.equals(0):
                return S.Zero
        except Exception:
            pass
        return None
    finally:
        _active = False
