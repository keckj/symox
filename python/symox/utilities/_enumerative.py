"""Multiset partition algorithms adapted from SymPy; see LICENSE-SYMPY."""
class PartComponent:
    """Internal class used in support of the multiset partitions
    enumerators and the associated visitor functions."""
    __slots__ = ('c', 'u', 'v')

    def __init__(self):
        self.c = 0
        self.u = 0
        self.v = 0

    def __repr__(self):
        """for debug/algorithm animation purposes"""
        return 'c:%d u:%d v:%d' % (self.c, self.u, self.v)

    def __eq__(self, other):
        """Define  value oriented equality, which is useful for testers"""
        return isinstance(other, self.__class__) and self.c == other.c and (self.u == other.u) and (self.v == other.v)

    def __ne__(self, other):
        """Defined for consistency with __eq__"""
        return not self == other

def multiset_partitions_taocp(multiplicities):
    """Enumerates partitions of a multiset."""
    m = len(multiplicities)
    n = sum(multiplicities)
    pstack = [PartComponent() for i in range(n * m + 1)]
    f = [0] * (n + 1)
    for j in range(m):
        ps = pstack[j]
        ps.c = j
        ps.u = multiplicities[j]
        ps.v = multiplicities[j]
    f[0] = 0
    a = 0
    lpart = 0
    f[1] = m
    b = m
    while True:
        while True:
            k = b
            x = False
            for j in range(a, b):
                pstack[k].u = pstack[j].u - pstack[j].v
                if pstack[k].u == 0:
                    x = True
                elif not x:
                    pstack[k].c = pstack[j].c
                    pstack[k].v = min(pstack[j].v, pstack[k].u)
                    x = pstack[k].u < pstack[j].v
                    k = k + 1
                else:
                    pstack[k].c = pstack[j].c
                    pstack[k].v = pstack[k].u
                    k = k + 1
            if k > b:
                a = b
                b = k
                lpart = lpart + 1
                f[lpart + 1] = b
            else:
                break
        state = [f, lpart, pstack]
        yield state
        while True:
            j = b - 1
            while pstack[j].v == 0:
                j = j - 1
            if j == a and pstack[j].v == 1:
                if lpart == 0:
                    return
                lpart = lpart - 1
                b = a
                a = f[lpart]
            else:
                pstack[j].v = pstack[j].v - 1
                for k in range(j + 1, b):
                    pstack[k].v = pstack[k].u
                break

def factoring_visitor(state, primes):
    """Use with multiset_partitions_taocp to enumerate the ways a
    number can be expressed as a product of factors.  For this usage,
    the exponents of the prime factors of a number are arguments to
    the partition enumerator, while the corresponding prime factors
    are input here."""
    (f, lpart, pstack) = state
    factoring = []
    for i in range(lpart + 1):
        factor = 1
        for ps in pstack[f[i]:f[i + 1]]:
            if ps.v > 0:
                factor *= primes[ps.c] ** ps.v
        factoring.append(factor)
    return factoring

def list_visitor(state, components):
    """Return a list of lists to represent the partition."""
    (f, lpart, pstack) = state
    partition = []
    for i in range(lpart + 1):
        part = []
        for ps in pstack[f[i]:f[i + 1]]:
            if ps.v > 0:
                part.extend([components[ps.c]] * ps.v)
        partition.append(part)
    return partition

class MultisetPartitionTraverser:
    """Has methods to ``enumerate`` and ``count`` the partitions of a multiset."""

    def __init__(self):
        self.debug = False
        self.k1 = 0
        self.k2 = 0
        self.p1 = 0
        self.pstack = None
        self.f = None
        self.lpart = 0
        self.discarded = 0
        self.dp_stack = []
        if not hasattr(self, 'dp_map'):
            self.dp_map = {}

    def db_trace(self, msg):
        """Useful for understanding/debugging the algorithms.  Not
        generally activated in end-user code."""
        if self.debug:
            raise RuntimeError

    def _initialize_enumeration(self, multiplicities):
        """Allocates and initializes the partition stack."""
        num_components = len(multiplicities)
        cardinality = sum(multiplicities)
        self.pstack = [PartComponent() for i in range(num_components * cardinality + 1)]
        self.f = [0] * (cardinality + 1)
        for j in range(num_components):
            ps = self.pstack[j]
            ps.c = j
            ps.u = multiplicities[j]
            ps.v = multiplicities[j]
        self.f[0] = 0
        self.f[1] = num_components
        self.lpart = 0

    def decrement_part(self, part):
        """Decrements part (a subrange of pstack), if possible, returning
        True iff the part was successfully decremented."""
        plen = len(part)
        for j in range(plen - 1, -1, -1):
            if j == 0 and part[j].v > 1 or (j > 0 and part[j].v > 0):
                part[j].v -= 1
                for k in range(j + 1, plen):
                    part[k].v = part[k].u
                return True
        return False

    def decrement_part_small(self, part, ub):
        """Decrements part (a subrange of pstack), if possible, returning
        True iff the part was successfully decremented."""
        if self.lpart >= ub - 1:
            self.p1 += 1
            return False
        plen = len(part)
        for j in range(plen - 1, -1, -1):
            if j == 0 and (part[0].v - 1) * (ub - self.lpart) < part[0].u:
                self.k1 += 1
                return False
            if j == 0 and part[j].v > 1 or (j > 0 and part[j].v > 0):
                part[j].v -= 1
                for k in range(j + 1, plen):
                    part[k].v = part[k].u
                if plen > 1 and part[1].v == 0 and (part[0].u - part[0].v == (ub - self.lpart - 1) * part[0].v):
                    self.k2 += 1
                    self.db_trace('Decrement fails test 3')
                    return False
                return True
        return False

    def decrement_part_large(self, part, amt, lb):
        """Decrements part, while respecting size constraint."""
        if amt == 1:
            if not self.decrement_part(part):
                return False
        min_unalloc = lb - self.lpart
        if min_unalloc <= 0:
            return True
        total_mult = sum((pc.u for pc in part))
        total_alloc = sum((pc.v for pc in part))
        if total_mult <= min_unalloc:
            return False
        deficit = min_unalloc - (total_mult - total_alloc)
        if deficit <= 0:
            return True
        for i in range(len(part) - 1, -1, -1):
            if i == 0:
                if part[0].v > deficit:
                    part[0].v -= deficit
                    return True
                else:
                    return False
            elif part[i].v >= deficit:
                part[i].v -= deficit
                return True
            else:
                deficit -= part[i].v
                part[i].v = 0

    def decrement_part_range(self, part, lb, ub):
        """Decrements part (a subrange of pstack), if possible, returning
        True iff the part was successfully decremented."""
        return self.decrement_part_small(part, ub) and self.decrement_part_large(part, 0, lb)

    def spread_part_multiplicity(self):
        """Returns True if a new part has been created, and
        adjusts pstack, f and lpart as needed."""
        j = self.f[self.lpart]
        k = self.f[self.lpart + 1]
        base = k
        changed = False
        for j in range(self.f[self.lpart], self.f[self.lpart + 1]):
            self.pstack[k].u = self.pstack[j].u - self.pstack[j].v
            if self.pstack[k].u == 0:
                changed = True
            else:
                self.pstack[k].c = self.pstack[j].c
                if changed:
                    self.pstack[k].v = self.pstack[k].u
                elif self.pstack[k].u < self.pstack[j].v:
                    self.pstack[k].v = self.pstack[k].u
                    changed = True
                else:
                    self.pstack[k].v = self.pstack[j].v
                k = k + 1
        if k > base:
            self.lpart = self.lpart + 1
            self.f[self.lpart + 1] = k
            return True
        return False

    def top_part(self):
        """Return current top part on the stack, as a slice of pstack."""
        return self.pstack[self.f[self.lpart]:self.f[self.lpart + 1]]

    def enum_all(self, multiplicities):
        """Enumerate the partitions of a multiset."""
        self._initialize_enumeration(multiplicities)
        while True:
            while self.spread_part_multiplicity():
                pass
            state = [self.f, self.lpart, self.pstack]
            yield state
            while not self.decrement_part(self.top_part()):
                if self.lpart == 0:
                    return
                self.lpart -= 1

    def enum_small(self, multiplicities, ub):
        """Enumerate multiset partitions with no more than ``ub`` parts."""
        self.discarded = 0
        if ub <= 0:
            return
        self._initialize_enumeration(multiplicities)
        while True:
            while self.spread_part_multiplicity():
                self.db_trace('spread 1')
                if self.lpart >= ub:
                    self.discarded += 1
                    self.db_trace('  Discarding')
                    self.lpart = ub - 2
                    break
            else:
                state = [self.f, self.lpart, self.pstack]
                yield state
            while not self.decrement_part_small(self.top_part(), ub):
                self.db_trace('Failed decrement, going to backtrack')
                if self.lpart == 0:
                    return
                self.lpart -= 1
                self.db_trace('Backtracked to')
            self.db_trace('decrement ok, about to expand')

    def enum_large(self, multiplicities, lb):
        """Enumerate the partitions of a multiset with lb < num(parts)"""
        self.discarded = 0
        if lb >= sum(multiplicities):
            return
        self._initialize_enumeration(multiplicities)
        self.decrement_part_large(self.top_part(), 0, lb)
        while True:
            good_partition = True
            while self.spread_part_multiplicity():
                if not self.decrement_part_large(self.top_part(), 0, lb):
                    self.discarded += 1
                    good_partition = False
                    break
            if good_partition:
                state = [self.f, self.lpart, self.pstack]
                yield state
            while not self.decrement_part_large(self.top_part(), 1, lb):
                if self.lpart == 0:
                    return
                self.lpart -= 1

    def enum_range(self, multiplicities, lb, ub):
        """Enumerate the partitions of a multiset with
        ``lb < num(parts) <= ub``."""
        self.discarded = 0
        if ub <= 0 or lb >= sum(multiplicities):
            return
        self._initialize_enumeration(multiplicities)
        self.decrement_part_large(self.top_part(), 0, lb)
        while True:
            good_partition = True
            while self.spread_part_multiplicity():
                self.db_trace('spread 1')
                if not self.decrement_part_large(self.top_part(), 0, lb):
                    self.db_trace('  Discarding (large cons)')
                    self.discarded += 1
                    good_partition = False
                    break
                elif self.lpart >= ub:
                    self.discarded += 1
                    good_partition = False
                    self.db_trace('  Discarding small cons')
                    self.lpart = ub - 2
                    break
            if good_partition:
                state = [self.f, self.lpart, self.pstack]
                yield state
            while not self.decrement_part_range(self.top_part(), lb, ub):
                self.db_trace('Failed decrement, going to backtrack')
                if self.lpart == 0:
                    return
                self.lpart -= 1
                self.db_trace('Backtracked to')
            self.db_trace('decrement ok, about to expand')

    def count_partitions_slow(self, multiplicities):
        """Returns the number of partitions of a multiset whose elements
        have the multiplicities given in ``multiplicities``."""
        self.pcount = 0
        self._initialize_enumeration(multiplicities)
        while True:
            while self.spread_part_multiplicity():
                pass
            self.pcount += 1
            while not self.decrement_part(self.top_part()):
                if self.lpart == 0:
                    return self.pcount
                self.lpart -= 1

    def count_partitions(self, multiplicities):
        """Returns the number of partitions of a multiset whose components
        have the multiplicities given in ``multiplicities``."""
        self.pcount = 0
        self.dp_stack = []
        self._initialize_enumeration(multiplicities)
        pkey = part_key(self.top_part())
        self.dp_stack.append([(pkey, 0)])
        while True:
            while self.spread_part_multiplicity():
                pkey = part_key(self.top_part())
                if pkey in self.dp_map:
                    self.pcount += self.dp_map[pkey] - 1
                    self.lpart -= 1
                    break
                else:
                    self.dp_stack.append([(pkey, self.pcount)])
            self.pcount += 1
            while not self.decrement_part(self.top_part()):
                for (key, oldcount) in self.dp_stack.pop():
                    self.dp_map[key] = self.pcount - oldcount
                if self.lpart == 0:
                    return self.pcount
                self.lpart -= 1
            pkey = part_key(self.top_part())
            self.dp_stack[-1].append((pkey, self.pcount))

def part_key(part):
    """Helper for MultisetPartitionTraverser.count_partitions that
    creates a key for ``part``, that only includes information which can
    affect the count for that part.  (Any irrelevant information just
    reduces the effectiveness of dynamic programming.)"""
    rval = []
    for ps in part:
        rval.append(ps.u)
        rval.append(ps.v)
    return tuple(rval)
