"""Combinatorial iterators adapted from SymPy 1.14 (see LICENSE-SYMPY).

The algorithms run on Python containers and native Symox expressions.
"""
from collections import Counter, defaultdict, OrderedDict
from itertools import chain, combinations, combinations_with_replacement, cycle, islice, permutations, product, groupby
from operator import gt
from .iterables import group, variations, has_dups, has_variety, sift
from ._enumerative import multiset_partitions_taocp, list_visitor, MultisetPartitionTraverser

def as_int(value):
    from operator import index
    try:
        return index(value)
    except TypeError:
        raise ValueError('Expected an integer') from None

def is_palindromic(s, i=0, j=None):
    """Return True if the sequence is the same from left to right as it
    is from right to left in the whole sequence (default) or in the
    Python slice ``s[i: j]``; else False."""
    (i, j, _) = slice(i, j).indices(len(s))
    m = (j - i) // 2
    return all((s[i + k] == s[j - 1 - k] for k in range(m)))

def _iproduct2(iterable1, iterable2):
    """Cartesian product of two possibly infinite iterables"""
    it1 = iter(iterable1)
    it2 = iter(iterable2)
    elems1 = []
    elems2 = []
    sentinel = object()

    def append(it, elems):
        e = next(it, sentinel)
        if e is not sentinel:
            elems.append(e)
    n = 0
    append(it1, elems1)
    append(it2, elems2)
    if not elems1 or not elems2:
        return
    while n <= len(elems1) + len(elems2):
        for m in range(n - len(elems1) + 1, len(elems2)):
            yield (elems1[n - m], elems2[m])
        n += 1
        append(it1, elems1)
        append(it2, elems2)

def iproduct(*iterables):
    """Cartesian product of iterables."""
    if len(iterables) == 0:
        yield ()
        return
    elif len(iterables) == 1:
        for e in iterables[0]:
            yield (e,)
    elif len(iterables) == 2:
        yield from _iproduct2(*iterables)
    else:
        (first, others) = (iterables[0], iterables[1:])
        for (ef, eo) in _iproduct2(first, iproduct(*others)):
            yield ((ef,) + eo)

def multiset(seq):
    """Return the hashable sequence in multiset form with values being the
    multiplicity of the item in the sequence."""
    return dict(Counter(seq).items())

def ibin(n, bits=None, str=False):
    """Return a list of length ``bits`` corresponding to the binary value
    of ``n`` with small bits to the right (last). If bits is omitted, the
    length will be the number required to represent ``n``. If the bits are
    desired in reversed order, use the ``[::-1]`` slice of the returned list."""
    if n < 0:
        raise ValueError('negative numbers are not allowed')
    n = as_int(n)
    if bits is None:
        bits = 0
    else:
        try:
            bits = as_int(bits)
        except ValueError:
            bits = -1
        else:
            if n.bit_length() > bits:
                raise ValueError('`bits` must be >= {}'.format(n.bit_length()))
    if not str:
        if bits >= 0:
            return [1 if i == '1' else 0 for i in bin(n)[2:].rjust(bits, '0')]
        else:
            return variations(range(2), n, repetition=True)
    elif bits >= 0:
        return bin(n)[2:].rjust(bits, '0')
    else:
        return (bin(i)[2:].rjust(n, '0') for i in range(2 ** n))

def filter_symbols(iterator, exclude):
    """Only yield elements from `iterator` that do not occur in `exclude`."""
    exclude = set(exclude)
    for s in iterator:
        if s not in exclude:
            yield s

def numbered_symbols(prefix='x', cls=None, start=0, exclude=(), *args, **assumptions):
    """Generate an infinite stream of Symbols consisting of a prefix and
    increasing subscripts provided that they do not occur in ``exclude``."""
    exclude = set(exclude or [])
    if cls is None:
        from symox import Symbol
        cls = Symbol
    while True:
        name = '%s%s' % (prefix, start)
        s = cls(name, *args, **assumptions)
        if s not in exclude:
            yield s
        start += 1

def capture(func):
    """Return the printed output of func()."""
    from io import StringIO
    import sys
    stdout = sys.stdout
    sys.stdout = file = StringIO()
    try:
        func()
    finally:
        sys.stdout = stdout
    return file.getvalue()

def common_prefix(*seqs):
    """Return the subsequence that is a common start of sequences in ``seqs``."""
    if not all(seqs):
        return []
    elif len(seqs) == 1:
        return seqs[0]
    i = 0
    for i in range(min((len(s) for s in seqs))):
        if not all((seqs[j][i] == seqs[0][i] for j in range(len(seqs)))):
            break
    else:
        i += 1
    return seqs[0][:i]

def common_suffix(*seqs):
    """Return the subsequence that is a common ending of sequences in ``seqs``."""
    if not all(seqs):
        return []
    elif len(seqs) == 1:
        return seqs[0]
    i = 0
    for i in range(-1, -min((len(s) for s in seqs)) - 1, -1):
        if not all((seqs[j][i] == seqs[0][i] for j in range(len(seqs)))):
            break
    else:
        i -= 1
    if i == -1:
        return []
    else:
        return seqs[0][i + 1:]

def strongly_connected_components(G):
    """Strongly connected components of a directed graph in reverse topological
    order."""
    (V, E) = G
    Gmap = {vi: [] for vi in V}
    for (v1, v2) in E:
        Gmap[v1].append(v2)
    return _strongly_connected_components(V, Gmap)

def _strongly_connected_components(V, Gmap):
    """More efficient internal routine for strongly_connected_components"""
    lowlink = {}
    indices = {}
    stack = OrderedDict()
    callstack = []
    components = []
    nomore = object()

    def start(v):
        index = len(stack)
        indices[v] = lowlink[v] = index
        stack[v] = None
        callstack.append((v, iter(Gmap[v])))

    def finish(v1):
        if lowlink[v1] == indices[v1]:
            component = [stack.popitem()[0]]
            while component[-1] is not v1:
                component.append(stack.popitem()[0])
            components.append(component[::-1])
        (v2, _) = callstack.pop()
        if callstack:
            (v1, _) = callstack[-1]
            lowlink[v1] = min(lowlink[v1], lowlink[v2])
    for v in V:
        if v in indices:
            continue
        start(v)
        while callstack:
            (v1, it1) = callstack[-1]
            v2 = next(it1, nomore)
            if v2 is nomore:
                finish(v1)
            elif v2 not in indices:
                start(v2)
            elif v2 in stack:
                lowlink[v1] = min(lowlink[v1], indices[v2])
    return components

def connected_components(G):
    """Connected components of an undirected graph or weakly connected components
    of a directed graph."""
    (V, E) = G
    E_undirected = []
    for (v1, v2) in E:
        E_undirected.extend([(v1, v2), (v2, v1)])
    return strongly_connected_components((V, E_undirected))

def rotate_left(x, y):
    """Left rotates a list x by the number of steps specified
    in y."""
    if len(x) == 0:
        return []
    y = y % len(x)
    return x[y:] + x[:y]

def rotate_right(x, y):
    """Right rotates a list x by the number of steps specified
    in y."""
    if len(x) == 0:
        return []
    y = len(x) - y % len(x)
    return x[y:] + x[:y]

def least_rotation(x, key=None):
    """Returns the number of steps of left rotation required to
    obtain lexicographically minimal string/list/tuple, etc."""
    Id = lambda x: x
    if key is None:
        key = Id
    S = x + x
    f = [-1] * len(S)
    k = 0
    for j in range(1, len(S)):
        sj = S[j]
        i = f[j - k - 1]
        while i != -1 and sj != S[k + i + 1]:
            if key(sj) < key(S[k + i + 1]):
                k = j - i - 1
            i = f[i]
        if sj != S[k + i + 1]:
            if key(sj) < key(S[k]):
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return k

def multiset_combinations(m, n, g=None):
    """Return the unique combinations of size ``n`` from multiset ``m``."""
    from symox.utilities.sorting import ordered
    if g is None:
        if isinstance(m, dict):
            if any((as_int(v) < 0 for v in m.values())):
                raise ValueError('counts cannot be negative')
            N = sum(m.values())
            if n > N:
                return
            g = [[k, m[k]] for k in ordered(m)]
        else:
            m = list(m)
            N = len(m)
            if n > N:
                return
            try:
                m = multiset(m)
                g = [(k, m[k]) for k in ordered(m)]
            except TypeError:
                m = list(ordered(m))
                g = [list(i) for i in group(m, multiple=False)]
        del m
    else:
        N = sum((v for (k, v) in g))
    if n > N or not n:
        yield []
    else:
        for (i, (k, v)) in enumerate(g):
            if v >= n:
                yield ([k] * n)
                v = n - 1
            for v in range(min(n, v), 0, -1):
                for j in multiset_combinations(None, n - v, g[i + 1:]):
                    rv = [k] * v + j
                    if len(rv) == n:
                        yield rv

def multiset_permutations(m, size=None, g=None):
    """Return the unique permutations of multiset ``m``."""
    from symox.utilities.sorting import ordered
    if g is None:
        if isinstance(m, dict):
            if any((as_int(v) < 0 for v in m.values())):
                raise ValueError('counts cannot be negative')
            g = [[k, m[k]] for k in ordered(m)]
        else:
            m = list(ordered(m))
            g = [list(i) for i in group(m, multiple=False)]
        del m
    do = [gi for gi in g if gi[1] > 0]
    SUM = sum((gi[1] for gi in do))
    if not do or (size is not None and (size > SUM or size < 1)):
        if not do and size is None or size == 0:
            yield []
        return
    elif size == 1:
        for (k, v) in do:
            yield [k]
    elif len(do) == 1:
        (k, v) = do[0]
        v = v if size is None else size if size <= v else 0
        yield [k for i in range(v)]
    elif all((v == 1 for (k, v) in do)):
        for p in permutations([k for (k, v) in do], size):
            yield list(p)
    else:
        size = size if size is not None else SUM
        for (i, (k, v)) in enumerate(do):
            do[i][1] -= 1
            for j in multiset_permutations(None, size - 1, do):
                if j:
                    yield ([k] + j)
            do[i][1] += 1

def _partition(seq, vector, m=None):
    """Return the partition of seq as specified by the partition vector."""
    if m is None:
        m = max(vector) + 1
    elif isinstance(vector, int):
        (vector, m) = (m, vector)
    p = [[] for i in range(m)]
    for (i, v) in enumerate(vector):
        p[v].append(seq[i])
    return p

def _set_partitions(n):
    """Cycle through all partitions of n elements, yielding the
    current number of partitions, ``m``, and a mutable list, ``q``
    such that ``element[i]`` is in part ``q[i]`` of the partition."""
    n = as_int(n)
    if n < 0:
        raise ValueError("Set size must be nonnegative")
    if n == 0:
        yield (0, [])
        return
    p = [0] * n
    q = [0] * n
    nc = 1
    yield (nc, q)
    while nc != n:
        m = n
        while 1:
            m -= 1
            i = q[m]
            if p[i] != 1:
                break
            q[m] = 0
        i += 1
        q[m] = i
        m += 1
        nc += m - n
        p[0] += n - m
        if i == nc:
            p[nc] = 0
            nc += 1
        p[i - 1] -= 1
        p[i] += 1
        yield (nc, q)

def multiset_partitions(multiset, m=None):
    """Return unique partitions of the given multiset (in list form).
    If ``m`` is None, all multisets will be returned, otherwise only
    partitions with ``m`` parts will be returned."""
    if isinstance(multiset, int):
        n = multiset
        if m and m > n:
            return
        multiset = list(range(n))
        if m == 1:
            yield [multiset[:]]
            return
        for (nc, q) in _set_partitions(n):
            if m is None or nc == m:
                rv = [[] for i in range(nc)]
                for i in range(n):
                    rv[q[i]].append(multiset[i])
                yield rv
        return
    if len(multiset) == 1 and isinstance(multiset, str):
        multiset = [multiset]
    if not has_variety(multiset):
        n = len(multiset)
        if m and m > n:
            return
        if m == 1:
            yield [multiset[:]]
            return
        x = multiset[:1]
        for (size, p) in partitions(n, m, size=True):
            if m is None or size == m:
                rv = []
                for k in sorted(p):
                    rv.extend([x * k] * p[k])
                yield rv
    else:
        from symox.utilities.sorting import ordered
        multiset = list(ordered(multiset))
        n = len(multiset)
        if m and m > n:
            return
        if m == 1:
            yield [multiset[:]]
            return
        (elements, multiplicities) = zip(*group(multiset, False))
        if len(elements) < len(multiset):
            if m:
                mpt = MultisetPartitionTraverser()
                for state in mpt.enum_range(multiplicities, m - 1, m):
                    yield list_visitor(state, elements)
            else:
                for state in multiset_partitions_taocp(multiplicities):
                    yield list_visitor(state, elements)
        else:
            for (nc, q) in _set_partitions(n):
                if m is None or nc == m:
                    rv = [[] for i in range(nc)]
                    for i in range(n):
                        rv[q[i]].append(i)
                    yield [[multiset[j] for j in i] for i in rv]

def partitions(n, m=None, k=None, size=False):
    """Generate all partitions of positive integer, n."""
    if n <= 0 or (m is not None and m < 1) or (k is not None and k < 1) or (m and k and (m * k < n)):
        if size:
            yield (0, {})
        else:
            yield {}
        return
    if m is None:
        m = n
    else:
        m = min(m, n)
    k = min(k or n, n)
    (n, m, k) = (as_int(n), as_int(m), as_int(k))
    (q, r) = divmod(n, k)
    ms = {k: q}
    keys = [k]
    if r:
        ms[r] = 1
        keys.append(r)
    room = m - q - bool(r)
    if size:
        yield (sum(ms.values()), ms.copy())
    else:
        yield ms.copy()
    while keys != [1]:
        if keys[-1] == 1:
            del keys[-1]
            reuse = ms.pop(1)
            room += reuse
        else:
            reuse = 0
        while 1:
            i = keys[-1]
            newcount = ms[i] = ms[i] - 1
            reuse += i
            if newcount == 0:
                del keys[-1], ms[i]
            room += 1
            i -= 1
            (q, r) = divmod(reuse, i)
            need = q + bool(r)
            if need > room:
                if not keys:
                    return
                continue
            ms[i] = q
            keys.append(i)
            if r:
                ms[r] = 1
                keys.append(r)
            break
        room -= need
        if size:
            yield (sum(ms.values()), ms.copy())
        else:
            yield ms.copy()

def ordered_partitions(n, m=None, sort=True):
    """Generates ordered partitions of integer *n*."""
    if n < 1 or (m is not None and m < 1):
        yield []
        return
    if m is None:
        a = [1] * n
        y = -1
        v = n
        while v > 0:
            v -= 1
            x = a[v] + 1
            while y >= 2 * x:
                a[v] = x
                y -= x
                v += 1
            w = v + 1
            while x <= y:
                a[v] = x
                a[w] = y
                yield a[:w + 1]
                x += 1
                y -= 1
            a[v] = x + y
            y = a[v] - 1
            yield a[:w]
    elif m == 1:
        yield [n]
    elif n == m:
        yield ([1] * n)
    else:
        for b in range(1, n // m + 1):
            a = [b] * m
            x = n - b * m
            if not x:
                if sort:
                    yield a
            elif not sort and x <= m:
                for ax in ordered_partitions(x, sort=False):
                    mi = len(ax)
                    a[-mi:] = [i + b for i in ax]
                    yield a
                    a[-mi:] = [b] * mi
            else:
                for mi in range(1, m):
                    for ax in ordered_partitions(x, mi, sort=True):
                        a[-mi:] = [i + b for i in ax]
                        yield a
                        a[-mi:] = [b] * mi

def binary_partitions(n):
    """Generates the binary partition of *n*."""
    from math import ceil, log2
    power = int(2 ** ceil(log2(n)))
    acc = 0
    partition = []
    while power:
        if acc + power <= n:
            partition.append(power)
            acc += power
        power >>= 1
    last_num = len(partition) - 1 - (n & 1)
    while last_num >= 0:
        yield partition
        if partition[last_num] == 2:
            partition[last_num] = 1
            partition.append(1)
            last_num -= 1
            continue
        partition.append(1)
        partition[last_num] >>= 1
        x = partition[last_num + 1] = partition[last_num]
        last_num += 1
        while x > 1:
            if x <= len(partition) - last_num - 1:
                del partition[-x + 1:]
                last_num += 1
                partition[last_num] = x
            else:
                x >>= 1
    yield ([1] * n)

def uniq(seq, result=None):
    """Yield unique elements from ``seq`` as an iterator. The second
    parameter ``result``  is used internally; it is not necessary
    to pass anything for this."""
    try:
        n = len(seq)
    except TypeError:
        n = None

    def check():
        if n is not None and len(seq) != n:
            raise RuntimeError('sequence changed size during iteration')
    try:
        seen = set()
        result = result or []
        for (i, s) in enumerate(seq):
            if not (s in seen or seen.add(s)):
                yield s
                check()
    except TypeError:
        if s not in result:
            yield s
            check()
            result.append(s)
        if hasattr(seq, '__getitem__'):
            yield from uniq(seq[i + 1:], result)
        else:
            yield from uniq(seq, result)

def generate_bell(n):
    """Return permutations of [0, 1, ..., n - 1] such that each permutation
    differs from the last by the exchange of a single pair of neighbors.
    The ``n!`` permutations are returned as an iterator. In order to obtain
    the next permutation from a random starting permutation, use the
    ``next_trotterjohnson`` method of the Permutation class (which generates
    the same sequence in a different manner)."""
    n = as_int(n)
    if n < 1:
        raise ValueError('n must be a positive integer')
    if n == 1:
        yield (0,)
    elif n == 2:
        yield (0, 1)
        yield (1, 0)
    elif n == 3:
        yield from [(0, 1, 2), (0, 2, 1), (2, 0, 1), (2, 1, 0), (1, 2, 0), (1, 0, 2)]
    else:
        m = n - 1
        op = [0] + [-1] * m
        l = list(range(n))
        while True:
            yield tuple(l)
            big = (None, -1)
            for i in range(n):
                if op[i] and l[i] > big[1]:
                    big = (i, l[i])
            (i, _) = big
            if i is None:
                break
            j = i + op[i]
            (l[i], l[j]) = (l[j], l[i])
            (op[i], op[j]) = (op[j], op[i])
            if j == 0 or j == m or l[j + op[j]] > l[j]:
                op[j] = 0
            for i in range(j):
                if l[i] > l[j]:
                    op[i] = 1
            for i in range(j + 1, n):
                if l[i] > l[j]:
                    op[i] = -1

def generate_involutions(n):
    """Generate self-inverse permutations in lexicographic order.

    Construct fixed points or pairs instead of filtering all n! permutations.
    """
    n = as_int(n)
    if n < 0:
        raise ValueError('Permutation size must be nonnegative')
    permutation = [-1] * n
    def visit():
        i = next((k for k, value in enumerate(permutation) if value < 0), n)
        if i == n:
            yield tuple(permutation)
            return
        for j in range(i, n):
            if permutation[j] >= 0:
                continue
            permutation[i] = j
            permutation[j] = i
            yield from visit()
            permutation[i] = permutation[j] = -1
    yield from visit()


def multiset_derangements(s):
    """Generate derangements of the elements of s *in place*."""
    from symox.utilities.sorting import ordered
    try:
        ms = multiset(s)
    except TypeError:
        key = dict(enumerate(ordered(uniq(s))))
        h = []
        for si in s:
            for k in key:
                if key[k] == si:
                    h.append(k)
                    break
        for i in multiset_derangements(h):
            yield [key[j] for j in i]
        return
    mx = max(ms.values())
    n = len(s)
    if mx * 2 > n:
        return
    if len(ms) == n:
        yield from _set_derangements(s)
        return
    for M in ms:
        if ms[M] == mx:
            break
    inonM = [i for i in range(n) if s[i] != M]
    iM = [i for i in range(n) if s[i] == M]
    rv = [None] * n
    if 2 * mx == n:
        for i in inonM:
            rv[i] = M
        for p in multiset_permutations([s[i] for i in inonM]):
            for (i, pi) in zip(iM, p):
                rv[i] = pi
            yield rv
        rv[:] = [None] * n
        return
    if n - 2 * mx == 1 and len(ms.values()) == n - mx + 1:
        for (i, i1) in enumerate(inonM):
            ifill = inonM[:i] + inonM[i + 1:]
            for j in ifill:
                rv[j] = M
            for p in permutations([s[j] for j in ifill]):
                rv[i1] = s[i1]
                for (j, pi) in zip(iM, p):
                    rv[j] = pi
                k = i1
                for j in iM:
                    (rv[j], rv[k]) = (rv[k], rv[j])
                    yield rv
                    k = j
        rv[:] = [None] * n
        return

    def finish_derangements():
        """Place the last two elements into the partially completed
        derangement, and yield the results."""
        a = take[1][0]
        a_ct = take[1][1]
        b = take[0][0]
        b_ct = take[0][1]
        forced_a = []
        forced_b = []
        open_free = []
        for i in range(len(s)):
            if rv[i] is None:
                if s[i] == a:
                    forced_b.append(i)
                elif s[i] == b:
                    forced_a.append(i)
                else:
                    open_free.append(i)
        if len(forced_a) > a_ct or len(forced_b) > b_ct:
            return
        for i in forced_a:
            rv[i] = a
        for i in forced_b:
            rv[i] = b
        for a_place in combinations(open_free, a_ct - len(forced_a)):
            for a_pos in a_place:
                rv[a_pos] = a
            for i in open_free:
                if rv[i] is None:
                    rv[i] = b
            yield rv
            for i in open_free:
                rv[i] = None
        for i in forced_a:
            rv[i] = None
        for i in forced_b:
            rv[i] = None

    def iopen(v):
        return [i for i in range(n) if rv[i] is None and s[i] != v]

    def do(j):
        if j == 1:
            yield from finish_derangements()
        else:
            (M, mx) = take[j]
            for i in combinations(iopen(M), mx):
                for ii in i:
                    rv[ii] = M
                yield from do(j - 1)
                for ii in i:
                    rv[ii] = None
    take = sorted(ms.items(), key=lambda x: (x[1], x[0]))
    yield from do(len(take) - 1)
    rv[:] = [None] * n

def _set_derangements(s):
    """yield derangements of items in ``s`` which are assumed to contain
    no repeated elements"""
    if len(s) < 2:
        return
    if len(s) == 2:
        yield [s[1], s[0]]
        return
    if len(s) == 3:
        yield [s[1], s[2], s[0]]
        yield [s[2], s[0], s[1]]
        return
    for p in permutations(s):
        if not any((i == j for (i, j) in zip(p, s))):
            yield list(p)

def generate_derangements(s):
    """Return unique derangements of the elements of iterable ``s``."""
    if not has_dups(s):
        yield from _set_derangements(s)
    else:
        for p in multiset_derangements(s):
            yield list(p)

def necklaces(n, k, free=False):
    """A routine to generate necklaces that may (free=True) or may not
    (free=False) be turned over to be viewed. The "necklaces" returned
    are comprised of ``n`` integers (beads) with ``k`` different
    values (colors). Only unique necklaces are returned."""
    if k == 0 and n > 0:
        return
    a = [0] * n
    yield tuple(a)
    if n == 0:
        return
    while True:
        i = n - 1
        while a[i] == k - 1:
            i -= 1
            if i == -1:
                return
        a[i] += 1
        for j in range(n - i - 1):
            a[j + i + 1] = a[j]
        if n % (i + 1) == 0 and (not free or all((a <= a[j::-1] + a[-1:j:-1] for j in range(n - 1)))):
            yield tuple(a)

def bracelets(n, k):
    """Wrapper to necklaces to return a free (unrestricted) necklace."""
    return necklaces(n, k, free=True)

def generate_oriented_forest(n):
    """This algorithm generates oriented forests."""
    n = as_int(n)
    if n < 0:
        raise ValueError("Forest size must be nonnegative")
    P = list(range(-1, n))
    while True:
        yield P[1:]
        if P[n] > 0:
            P[n] = P[P[n]]
        else:
            for p in range(n - 1, 0, -1):
                if P[p] != 0:
                    target = P[p] - 1
                    for q in range(p - 1, 0, -1):
                        if P[q] == target:
                            break
                    offset = p - q
                    for i in range(p, n + 1):
                        P[i] = P[i - offset]
                    break
            else:
                break

def minlex(seq, directed=True, key=None):
    """Return the rotation of the sequence in which the lexically smallest
    elements appear first, e.g. `cba \\rightarrow acb`."""
    Id = lambda x: x
    if key is None:
        key = Id
    best = rotate_left(seq, least_rotation(seq, key=key))
    if not directed:
        rseq = seq[::-1]
        rbest = rotate_left(rseq, least_rotation(rseq, key=key))
        best = min(best, rbest, key=key)
    return tuple(best) if not isinstance(seq, str) else best

def runs(seq, op=gt):
    """Group the sequence into lists in which successive elements
    all compare the same with the comparison operator, ``op``:
    op(seq[i + 1], seq[i]) is True from all elements in a run."""
    cycles = []
    seq = iter(seq)
    try:
        run = [next(seq)]
    except StopIteration:
        return []
    while True:
        try:
            ei = next(seq)
        except StopIteration:
            break
        if op(ei, run[-1]):
            run.append(ei)
            continue
        else:
            cycles.append(run)
            run = [ei]
    if run:
        cycles.append(run)
    return cycles

def sequence_partitions(l, n, /):
    """Returns the partition of sequence $l$ into $n$ bins"""
    if n == 1 and l:
        yield [l]
        return
    for i in range(1, len(l)):
        for part in sequence_partitions(l[i:], n - 1):
            yield ([l[:i]] + part)

def sequence_partitions_empty(l, n, /):
    """Returns the partition of sequence $l$ into $n$ bins with
    empty sequence"""
    if n < 1:
        return
    if n == 1:
        yield [l]
        return
    for i in range(0, len(l) + 1):
        for part in sequence_partitions_empty(l[i:], n - 1):
            yield ([l[:i]] + part)

def kbins(l, k, ordered=None):
    """Return sequence ``l`` partitioned into ``k`` bins."""
    if ordered is None:
        yield from sequence_partitions(l, k)
    elif ordered == 11:
        for pl in multiset_permutations(l):
            pl = list(pl)
            yield from sequence_partitions(pl, k)
    elif ordered == 0:
        yield from multiset_partitions(l, k)
    elif ordered == 10:
        for p in multiset_partitions(l, k):
            for perm in permutations(p):
                yield list(perm)
    elif ordered == 1:
        for (kgot, p) in partitions(len(l), k, size=True):
            if kgot != k:
                continue
            for li in multiset_permutations(l):
                rv = []
                i = j = 0
                li = list(li)
                for (size, multiplicity) in sorted(p.items()):
                    for m in range(multiplicity):
                        j = i + size
                        rv.append(li[i:j])
                        i = j
                yield rv
    else:
        raise ValueError('ordered must be one of 00, 01, 10 or 11, not %s' % ordered)

def permute_signs(t):
    """Return iterator in which the signs of non-zero elements
    of t are permuted."""
    for signs in product(*[(1, -1)] * (len(t) - t.count(0))):
        signs = list(signs)
        yield type(t)([i * signs.pop() if i else i for i in t])

def signed_permutations(t):
    """Return iterator in which the signs of non-zero elements
    of t and the order of the elements are permuted and all
    returned values are unique."""
    return (type(t)(i) for j in multiset_permutations(t) for i in permute_signs(j))

def roundrobin(*iterables):
    """roundrobin recipe taken from itertools documentation:
    https://docs.python.org/3/library/itertools.html#itertools-recipes"""
    nexts = cycle((iter(it).__next__ for it in iterables))
    pending = len(iterables)
    while pending:
        try:
            for nxt in nexts:
                yield nxt()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))

class NotIterable:
    """Use this as mixin when creating a class which is not supposed to
    return true when iterable() is called on its instances because
    calling list() on the instance, for example, would result in
    an infinite loop."""
    pass

def iterable(i, exclude=(str, dict, NotIterable)):
    """Return a boolean indicating whether ``i`` is SymPy iterable.
    True also indicates that the iterator is finite, e.g. you can
    call list(...) on the instance."""
    if hasattr(i, '_iterable'):
        flag = i._iterable
        # On a class such as ``sin`` the instance property itself is found;
        # like the reference, a class is then iterable only if iter() works.
        if not (isinstance(i, type) and isinstance(flag, property)):
            return flag
    try:
        iter(i)
    except TypeError:
        return False
    if exclude:
        return not isinstance(i, exclude)
    return True

def is_sequence(i, include=None):
    """Return a boolean indicating whether ``i`` is a sequence in the SymPy
    sense. If anything that fails the test below should be included as
    being a sequence for your application, set 'include' to that object's
    type; multiple types should be passed as a tuple of types."""
    return hasattr(i, '__getitem__') and iterable(i) or (bool(include) and isinstance(i, include))

def postorder_traversal(node, keys=None):
    from symox.utilities.traversal import postorder_traversal as _postorder_traversal
    return _postorder_traversal(node, keys=keys)

def default_sort_key(*args, **kwargs):
    from symox import default_sort_key as _default_sort_key
    return _default_sort_key(*args, **kwargs)


def random_derangement(t, choice=None, strict=True):
    """Randomize a derangement, or minimize fixed points when strict=False.

    Uses a bounded rotation construction; sampling is not uniform over all
    possible derangements. ``choice`` selects from a nonempty sequence.
    """
    if choice is None:
        from secrets import choice
    values = list(t)
    if not values:
        return '' if isinstance(t, str) else []
    groups = {}
    for i, value in enumerate(values):
        groups.setdefault(value, []).append(i)
    largest = max(map(len, groups.values()))
    if strict and 2 * largest > len(values):
        raise ValueError('no derangement possible')
    pending = list(groups.values())
    indices = []
    while pending:
        selected = choice(pending)
        pending.remove(selected)
        indices.extend(selected)
    offset = choice(range(largest, len(values)-largest+1)) if 2*largest <= len(values) else largest
    result = [None] * len(values)
    for i, target in enumerate(indices):
        result[target] = values[indices[(i+offset) % len(values)]]
    return ''.join(result) if isinstance(t,str) else result


def interactive_traversal(expr):
    """Select an expression subtree using indices, f/l/r, or d to finish."""
    from .traversal import _children
    from random import choice
    while True:
        print('Current expression:', expr)
        children = tuple(_children(expr, None))
        if not children:
            return expr
        for i, child in enumerate(children):
            print('[%d] %s' % (i, child))
        try:
            selected = input('Your choice [index,f,l,r,d,?]: ')
        except EOFError:
            return expr
        if selected in ('', 'd'):
            return expr
        if selected == '?':
            print('Index: select a child; f: first; l: last; r: random; d: done')
            continue
        if selected in ('f','l','r'):
            expr = children[0] if selected=='f' else children[-1] if selected=='l' else choice(children)
            continue
        try:
            index = int(selected)
            if not 0 <= index < len(children):
                raise ValueError
            expr = children[index]
        except ValueError:
            print('Choice must be a valid child index')
