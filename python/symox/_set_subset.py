"""Emptiness, finiteness and subset decisions for sets.

This is a port of the decision path of the reference ``Set.is_subset``
(trivial cases, ``_eval_is_subset``/``_eval_is_superset`` hooks, the
pairwise ``is_subset_sets`` rules and the intersection fallback) together
with the ``is_empty`` and ``is_finite_set`` properties it relies on.
"""
from itertools import product as _product

from .core.logic import fuzzy_and, fuzzy_bool, fuzzy_not, fuzzy_or

_INF_NAMES = frozenset(("Naturals", "Naturals0", "Integers", "Rationals", "Complexes", "UniversalSet"))


def _mod():
    import symox
    return symox


def _name(s):
    return getattr(s, "class_name", None)


def _is_interval(s):
    """Interval in the dispatch sense: the reference ``Reals`` is an Interval."""
    m = _mod()
    return isinstance(s, m.Interval) or _name(s) == "Reals"


def _interval_parts(s):
    m = _mod()
    if _name(s) == "Reals":
        return -m.oo, m.oo, True, True
    return s.start, s.end, bool(s.left_open), bool(s.right_open)


# emptiness / finiteness -----------------------------------------------------

def is_empty(s):
    m = _mod()
    name = _name(s)
    if s.is_EmptySet:
        return True
    if name == "Reals" or name in _INF_NAMES:
        return False
    if isinstance(s, m.Interval):
        start, end, lo, ro = _interval_parts(s)
        cond = (start >= end) if (lo or ro) else (start > end)
        return fuzzy_bool(cond)
    if isinstance(s, m.FiniteSet):
        return len(s.args) == 0
    if isinstance(s, m.Union):
        return fuzzy_and(a.is_empty for a in s.args)
    if isinstance(s, m.ProductSet):
        return fuzzy_or(a.is_empty for a in s.args)
    if isinstance(s, m.Range):
        size = _range_size(s)
        if size is None:
            return None
        return fuzzy_bool(m.sympify(size).is_zero)
    return None


def is_finite_set(s):
    m = _mod()
    name = _name(s)
    if s.is_EmptySet:
        return True
    if name == "Reals" or name in _INF_NAMES:
        return False
    if isinstance(s, m.FiniteSet):
        return True
    if isinstance(s, m.Interval):
        return s.measure.is_zero
    if isinstance(s, m.Union):
        return fuzzy_and(a.is_finite_set for a in s.args)
    if isinstance(s, m.Intersection):
        if fuzzy_or(a.is_finite_set for a in s.args):
            return True
        return None
    if isinstance(s, m.Complement):
        a, b = s.args
        a_finite = a.is_finite_set
        if a_finite is True:
            return True
        elif a_finite is False and b.is_finite_set:
            return False
        return None
    if isinstance(s, m.ProductSet):
        all_finite = fuzzy_and(a.is_finite_set for a in s.args)
        return fuzzy_or([s.is_empty, all_finite])
    if isinstance(s, m.Range):
        if s.start.is_integer and s.stop.is_integer:
            return True
        size = _range_size(s)
        return None if size is None else m.sympify(size).is_finite
    return None


def _range_size(r):
    """Size of a Range, or None when it cannot be computed."""
    m = _mod()
    start, stop, step = r.start, r.stop, r.step
    if all(v.is_Integer for v in (start, stop, step)):
        return len(range(int(start), int(stop), int(step)))
    if start.is_infinite or stop.is_infinite:
        if start == stop:
            return 0
        return m.oo
    return None


# pairwise rules (is_subset_sets) ------------------------------------------

def _measure_flag(s, attr):
    try:
        return getattr(s.measure, attr)
    except Exception:  # noqa: BLE001
        return None


def _interval_interval(a, b):
    a_start, a_end, a_lo, a_ro = _interval_parts(a)
    b_start, b_end, b_lo, b_ro = _interval_parts(b)
    m = _mod()
    if fuzzy_bool(a_start < b_start):
        return False
    if fuzzy_bool(a_end > b_end):
        return False
    if b_lo and not a_lo and fuzzy_bool(m.Eq(a_start, b_start)):
        return False
    if b_ro and not a_ro and fuzzy_bool(m.Eq(a_end, b_end)):
        return False
    return None


def _interval_finiteset(a, b):
    if fuzzy_not(_measure_flag(a, "is_zero")):
        return False
    return None


def _interval_union(a, b):
    m = _mod()
    if all(isinstance(s, (m.Interval, m.FiniteSet)) or _name(s) == "Reals" for s in b.args):
        intervals = [_interval_parts(s) for s in b.args if _is_interval(s)]
        a_start, a_end, _, _ = _interval_parts(a)
        if all(fuzzy_bool(a_start < s[0]) for s in intervals):
            return False
        if all(fuzzy_bool(a_end > s[1]) for s in intervals):
            return False
        if _measure_flag(a, "is_nonzero"):
            def no_overlap(s1, s2):
                return fuzzy_or([fuzzy_bool(s1[1] <= s2[0]), fuzzy_bool(s1[0] >= s2[1])])
            ap = (a_start, a_end)
            if all(no_overlap(s, ap) for s in intervals):
                return False
    return None


def _range_range(a, b):
    if a.step == b.step == 1:
        return fuzzy_and([fuzzy_bool(a.start >= b.start), fuzzy_bool(a.stop <= b.stop)])
    return None


def _range_interval(a, b):
    if a.step.is_positive:
        b_start, b_end, b_lo, b_ro = _interval_parts(b)
        if b_lo and a.inf.is_finite:
            cond_left = a.inf > b_start
        else:
            cond_left = a.inf >= b_start
        if b_ro and a.sup.is_finite:
            cond_right = a.sup < b_end
        else:
            cond_right = a.sup <= b_end
        return fuzzy_and([cond_left, cond_right])
    return None


def _range_finiteset(a, b):
    m = _mod()
    a_size = _range_size(a)
    if a_size is None:
        # symbolic Range of unknown size
        return None
    if fuzzy_bool(m.sympify(a_size) > len(b.args)):
        return False
    a_set = set(a)
    b_remaining = len(b.args)
    cnt_candidate = 0
    for x in b.args:
        if x.is_Integer:
            a_set.discard(x)
        elif fuzzy_not(x.is_integer):
            pass
        else:
            cnt_candidate += 1
        b_remaining -= 1
        if len(a_set) > b_remaining + cnt_candidate:
            return False
        if len(a_set) == 0:
            return True
    return None


def _interval_nonzero_false(a, b):
    if _measure_flag(a, "is_extended_nonzero"):
        return False
    return None


def _productset_finiteset(a, b):
    m = _mod()
    if not all(isinstance(s, m.FiniteSet) or s.is_EmptySet for s in a.args):
        return None
    return fuzzy_and(b.contains(m.Tuple(*t)) for t in _product(*(s.args for s in a.args)))


def is_subset_sets(a, b):
    m = _mod()
    an, bn = _name(a), _name(b)
    a_int, b_int = _is_interval(a), _is_interval(b)
    a_rng, b_rng = isinstance(a, m.Range), isinstance(b, m.Range)
    # most specific signatures first (Reals is an Interval subclass)
    if an == "Rationals" and bn == "Reals":
        return True
    if an == "Rationals" and b_rng:
        return False
    if an == "Complexes" and (b_int or b_rng or bn == "Rationals"):
        return False
    if a_rng and bn == "Complexes":
        return True
    if a_int:
        if b_int:
            return _interval_interval(a, b)
        if isinstance(b, m.FiniteSet):
            return _interval_finiteset(a, b)
        if isinstance(b, m.Union):
            return _interval_union(a, b)
        if b_rng or bn == "Rationals":
            return _interval_nonzero_false(a, b)
        return None
    if a_rng:
        if b_rng:
            return _range_range(a, b)
        if b_int:
            return _range_interval(a, b)
        if isinstance(b, m.FiniteSet):
            return _range_finiteset(a, b)
        return None
    if isinstance(a, m.ProductSet) and isinstance(b, m.FiniteSet):
        return _productset_finiteset(a, b)
    return None


# hooks -----------------------------------------------------------------------

def eval_is_subset(a, b):
    m = _mod()
    an = _name(a)
    if isinstance(a, m.FiniteSet):
        return fuzzy_and(b.contains(e) for e in a.args)
    if an == "Naturals":
        return m.Range(1, m.oo).is_subset(b)
    if an == "Naturals0":
        return m.Range(0, m.oo).is_subset(b)
    if an == "Integers":
        return m.Range(-m.oo, m.oo).is_subset(b)
    return None


def eval_is_superset(a, b):
    m = _mod()
    an = _name(a)
    if an == "Naturals":
        return m.Range(1, m.oo).is_superset(b)
    if an == "Naturals0":
        return m.Range(0, m.oo).is_superset(b)
    if an == "Integers":
        return m.Range(-m.oo, m.oo).is_superset(b)
    return None


def is_subset(a, b):
    m = _mod()
    if not isinstance(b, m.Set):
        raise ValueError("Unknown argument '%s'" % b)
    if isinstance(a, m.Union):
        return fuzzy_and(s.is_subset(b) for s in a.args)
    # Handle the trivial cases
    if a == b:
        return True
    empty = a.is_empty
    if empty is True:
        return True
    elif fuzzy_not(empty) and b.is_empty:
        return False
    if a.is_finite_set is False and b.is_finite_set:
        return False
    # Dispatch on subclass rules
    ret = a._eval_is_subset(b)
    if ret is not None:
        return ret
    ret = b._eval_is_superset(a)
    if ret is not None:
        return ret
    # Use pairwise rules
    ret = is_subset_sets(a, b)
    if ret is not None:
        return ret
    # Fall back on computing the intersection
    if a.intersect(b) == a:
        return True
    return None
