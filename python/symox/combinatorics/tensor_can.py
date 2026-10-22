"""Tensor canonicalization with the Butler-Portugal algorithm (ported from the reference tensor_can module)."""
from symox.combinatorics.permutations import Permutation, _af_rmul, _af_invert, _af_new
from symox.combinatorics.perm_groups import PermutationGroup, _orbit, _orbit_transversal
from symox.combinatorics.util import _distribute_gens_by_base, _orbits_transversals_from_bsgs

def dummy_sgs(dummies, sym, n):
    if len(dummies) > n:
        raise ValueError('List too large')
    res = []
    if sym is not None:
        for j in dummies[::2]:
            a = list(range(n + 2))
            if sym == 1:
                a[n] = n + 1
                a[n + 1] = n
            (a[j], a[j + 1]) = (a[j + 1], a[j])
            res.append(a)
    for j in dummies[:-3:2]:
        a = list(range(n + 2))
        a[j:j + 4] = (a[j + 2], a[j + 3], a[j], a[j + 1])
        res.append(a)
    return res

def _min_dummies(dummies, sym, indices):
    num_types = len(sym)
    m = [min(dx) if dx else None for dx in dummies]
    res = indices[:]
    for i in range(num_types):
        for (c, i) in enumerate(indices):
            for j in range(num_types):
                if i in dummies[j]:
                    res[c] = m[j]
                    break
    return res

def _trace_S(s, j, b, S_cosets):
    for h in S_cosets[b]:
        if s[h[b]] == j:
            return h
    return None

def _trace_D(gj, p_i, Dxtrav):
    for h in Dxtrav:
        if h[gj] == p_i:
            return h
    return None

def _dumx_remove(dumx, dumx_flat, p0):
    res = []
    for dx in dumx:
        if p0 not in dx:
            res.append(dx)
            continue
        k = dx.index(p0)
        if k % 2 == 0:
            p0_paired = dx[k + 1]
        else:
            p0_paired = dx[k - 1]
        dx.remove(p0)
        dx.remove(p0_paired)
        dumx_flat.remove(p0)
        dumx_flat.remove(p0_paired)
        res.append(dx)

def transversal2coset(size, base, transversal):
    a = []
    j = 0
    for i in range(size):
        if i in base:
            a.append(sorted(transversal[j].values()))
            j += 1
        else:
            a.append([list(range(size))])
    j = len(a) - 1
    while a[j] == [list(range(size))]:
        j -= 1
    return a[:j + 1]

def double_coset_can_rep(dummies, sym, b_S, sgens, S_transversals, g):
    size = g.size
    g = g.array_form
    num_dummies = size - 2
    indices = list(range(num_dummies))
    all_metrics_with_sym = not any((_ is None for _ in sym))
    num_types = len(sym)
    dumx = dummies[:]
    dumx_flat = []
    for dx in dumx:
        dumx_flat.extend(dx)
    b_S = b_S[:]
    sgensx = [h._array_form for h in sgens]
    if b_S:
        S_transversals = transversal2coset(size, b_S, S_transversals)
    dsgsx = []
    for i in range(num_types):
        dsgsx.extend(dummy_sgs(dumx[i], sym[i], num_dummies))
    idn = list(range(size))
    TAB = [(idn, idn, g)]
    for i in range(size - 2):
        b = i
        testb = b in b_S and sgensx
        if testb:
            sgensx1 = [_af_new(_) for _ in sgensx]
            deltab = _orbit(size, sgensx1, b)
        else:
            deltab = {b}
        if all_metrics_with_sym:
            md = _min_dummies(dumx, sym, indices)
        else:
            md = [min(_orbit(size, [_af_new(ddx) for ddx in dsgsx], ii)) for ii in range(size - 2)]
        p_i = min((min((md[h[x]] for x in deltab)) for (s, d, h) in TAB))
        dsgsx1 = [_af_new(_) for _ in dsgsx]
        Dxtrav = _orbit_transversal(size, dsgsx1, p_i, False, af=True) if dsgsx else None
        if Dxtrav:
            Dxtrav = [_af_invert(x) for x in Dxtrav]
        for ii in range(num_types):
            if p_i in dumx[ii]:
                if sym[ii] is not None:
                    deltap = dumx[ii]
                else:
                    p_i_index = dumx[ii].index(p_i) % 2
                    deltap = dumx[ii][p_i_index::2]
                break
        else:
            deltap = [p_i]
        TAB1 = []
        while TAB:
            (s, d, h) = TAB.pop()
            if min((md[h[x]] for x in deltab)) != p_i:
                continue
            deltab1 = [x for x in deltab if md[h[x]] == p_i]
            dg = _af_rmul(d, g)
            dginv = _af_invert(dg)
            sdeltab = [s[x] for x in deltab1]
            gdeltap = [dginv[x] for x in deltap]
            NEXT = [x for x in sdeltab if x in gdeltap]
            for j in NEXT:
                if testb:
                    s1 = _trace_S(s, j, b, S_transversals)
                    if not s1:
                        continue
                    else:
                        s1 = [s[ix] for ix in s1]
                else:
                    s1 = s
                if Dxtrav:
                    d1 = _trace_D(dg[j], p_i, Dxtrav)
                    if not d1:
                        continue
                else:
                    if p_i != dg[j]:
                        continue
                    d1 = idn
                assert d1[dg[j]] == p_i
                d1 = [d1[ix] for ix in d]
                h1 = [d1[g[ix]] for ix in s1]
                TAB1.append((s1, d1, h1))
        TAB1.sort(key=lambda x: x[-1])
        prev = [0] * size
        while TAB1:
            (s, d, h) = TAB1.pop()
            if h[:-2] == prev[:-2]:
                if h[-1] != prev[-1]:
                    return 0
            else:
                TAB.append((s, d, h))
            prev = h
        sgensx = [h for h in sgensx if h[b] == b]
        if b in b_S:
            b_S.remove(b)
        _dumx_remove(dumx, dumx_flat, p_i)
        dsgsx = []
        for i in range(num_types):
            dsgsx.extend(dummy_sgs(dumx[i], sym[i], num_dummies))
    return TAB[0][-1]

def canonical_free(base, gens, g, num_free):
    g = g.array_form
    size = len(g)
    if not base:
        return g[:]
    transversals = get_transversals(base, gens)
    for x in sorted(g[:-2]):
        if x not in base:
            base.append(x)
    h = g
    for transv in transversals:
        h_i = [size] * num_free
        s = None
        for sk in transv.values():
            h1 = _af_rmul(h, sk)
            hi = [h1.index(ix) for ix in range(num_free)]
            if hi < h_i:
                h_i = hi
                s = sk
        if s:
            h = _af_rmul(h, s)
    return h

def _get_map_slots(size, fixed_slots):
    res = list(range(size))
    pos = 0
    for i in range(size):
        if i in fixed_slots:
            continue
        res[i] = pos
        pos += 1
    return res

def _lift_sgens(size, fixed_slots, free, s):
    a = []
    j = k = 0
    fd = [y for (_, y) in sorted(zip(fixed_slots, free))]
    num_free = len(free)
    for i in range(size):
        if i in fixed_slots:
            a.append(fd[k])
            k += 1
        else:
            a.append(s[j] + num_free)
            j += 1
    return a

def canonicalize(g, dummies, msym, *v):
    from symox.combinatorics.testutil import canonicalize_naive
    if not isinstance(msym, list):
        if msym not in (0, 1, None):
            raise ValueError('msym must be 0, 1 or None')
        num_types = 1
    else:
        num_types = len(msym)
        if not all((msymx in (0, 1, None) for msymx in msym)):
            raise ValueError('msym entries must be 0, 1 or None')
        if len(dummies) != num_types:
            raise ValueError('dummies and msym must have the same number of elements')
    size = g.size
    num_tensors = 0
    v1 = []
    for (base_i, gens_i, n_i, sym_i) in v:
        if not _is_minimal_bsgs(base_i, gens_i):
            mbsgs = get_minimal_bsgs(base_i, gens_i)
            if not mbsgs:
                can = canonicalize_naive(g, dummies, msym, *v)
                return can
            (base_i, gens_i) = mbsgs
        v1.append((base_i, gens_i, [[]] * n_i, sym_i))
        num_tensors += n_i
    if num_types == 1 and (not isinstance(msym, list)):
        dummies = [dummies]
        msym = [msym]
    flat_dummies = []
    for dumx in dummies:
        flat_dummies.extend(dumx)
    if flat_dummies and flat_dummies != list(range(flat_dummies[0], flat_dummies[-1] + 1)):
        raise ValueError('dummies is not valid')
    (size1, sbase, sgens) = gens_products(*v1)
    if size != size1:
        raise ValueError('g has size %d, generators have size %d' % (size, size1))
    free = [i for i in range(size - 2) if i not in flat_dummies]
    num_free = len(free)
    g1 = canonical_free(sbase, sgens, g, num_free)
    if not flat_dummies:
        return g1
    sign = 0 if g1[-1] == size - 1 else 1
    start = 0
    for (i, (base_i, gens_i, n_i, sym_i)) in enumerate(v):
        free_i = []
        len_tens = gens_i[0].size - 2
        for j in range(n_i):
            h = g1[start:start + len_tens]
            fr = []
            for k in free:
                if k in h:
                    fr.append(h.index(k))
            free_i.append(fr)
            start += len_tens
        v1[i] = (base_i, gens_i, free_i, sym_i)
    (size, sbase, sgens) = gens_products(*v1)
    pos_free = [g1.index(x) for x in range(num_free)]
    size_red = size - num_free
    g1_red = [x - num_free for x in g1 if x in flat_dummies]
    if sign:
        g1_red.extend([size_red - 1, size_red - 2])
    else:
        g1_red.extend([size_red - 2, size_red - 1])
    map_slots = _get_map_slots(size, pos_free)
    sbase_red = [map_slots[i] for i in sbase if i not in pos_free]
    sgens_red = [_af_new([map_slots[i] for i in y._array_form if i not in pos_free]) for y in sgens]
    dummies_red = [[x - num_free for x in y] for y in dummies]
    transv_red = get_transversals(sbase_red, sgens_red)
    g1_red = _af_new(g1_red)
    g2 = double_coset_can_rep(dummies_red, msym, sbase_red, sgens_red, transv_red, g1_red)
    if g2 == 0:
        return 0
    g3 = _lift_sgens(size, pos_free, free, g2)
    return g3

def perm_af_direct_product(gens1, gens2, signed=True):
    gens1 = [list(x) for x in gens1]
    gens2 = [list(x) for x in gens2]
    s = 2 if signed else 0
    n1 = len(gens1[0]) - s
    n2 = len(gens2[0]) - s
    start = list(range(n1))
    end = list(range(n1, n1 + n2))
    if signed:
        gens1 = [gen[:-2] + end + [gen[-2] + n2, gen[-1] + n2] for gen in gens1]
        gens2 = [start + [x + n1 for x in gen] for gen in gens2]
    else:
        gens1 = [gen + end for gen in gens1]
        gens2 = [start + [x + n1 for x in gen] for gen in gens2]
    res = gens1 + gens2
    return res

def bsgs_direct_product(base1, gens1, base2, gens2, signed=True):
    s = 2 if signed else 0
    n1 = gens1[0].size - s
    base = list(base1)
    base += [x + n1 for x in base2]
    gens1 = [h._array_form for h in gens1]
    gens2 = [h._array_form for h in gens2]
    gens = perm_af_direct_product(gens1, gens2, signed)
    size = len(gens[0])
    id_af = list(range(size))
    gens = [h for h in gens if h != id_af]
    if not gens:
        gens = [id_af]
    return (base, [_af_new(h) for h in gens])

def get_symmetric_group_sgs(n, antisym=False):
    if n == 1:
        return ([], [_af_new(list(range(3)))])
    gens = [Permutation(n - 1)(i, i + 1)._array_form for i in range(n - 1)]
    if antisym == 0:
        gens = [x + [n, n + 1] for x in gens]
    else:
        gens = [x + [n + 1, n] for x in gens]
    base = list(range(n - 1))
    return (base, [_af_new(h) for h in gens])
riemann_bsgs = ([0, 2], [Permutation(0, 1)(4, 5), Permutation(2, 3)(4, 5), Permutation(5)(0, 2)(1, 3)])

def get_transversals(base, gens):
    if not base:
        return []
    stabs = _distribute_gens_by_base(base, gens)
    (orbits, transversals) = _orbits_transversals_from_bsgs(base, stabs)
    transversals = [{x: h._array_form for (x, h) in y.items()} for y in transversals]
    return transversals

def _is_minimal_bsgs(base, gens):
    base1 = []
    sgs1 = gens[:]
    size = gens[0].size
    for i in range(size):
        if not all((h._array_form[i] == i for h in sgs1)):
            base1.append(i)
            sgs1 = [h for h in sgs1 if h._array_form[i] == i]
    return base1 == base

def get_minimal_bsgs(base, gens):
    G = PermutationGroup(gens)
    (base, gens) = G.schreier_sims_incremental()
    if not _is_minimal_bsgs(base, gens):
        return None
    return (base, gens)

def tensor_gens(base, gens, list_free_indices, sym=0):

    def _get_bsgs(G, base, gens, free_indices):
        if not free_indices:
            return (base[:], gens[:])
        else:
            H = G.pointwise_stabilizer(free_indices)
            (base, sgs) = H.schreier_sims_incremental()
            return (base, sgs)
    if not base and list_free_indices.count([]) < 2:
        n = len(list_free_indices)
        size = gens[0].size
        size = n * (size - 2) + 2
        return (size, [], [_af_new(list(range(size)))])
    if any(list_free_indices):
        G = PermutationGroup(gens)
    else:
        G = None
    no_free = []
    size = gens[0].size
    id_af = list(range(size))
    num_indices = size - 2
    if not list_free_indices[0]:
        no_free.append(list(range(num_indices)))
    (res_base, res_gens) = _get_bsgs(G, base, gens, list_free_indices[0])
    for i in range(1, len(list_free_indices)):
        (base1, gens1) = _get_bsgs(G, base, gens, list_free_indices[i])
        (res_base, res_gens) = bsgs_direct_product(res_base, res_gens, base1, gens1, 1)
        if not list_free_indices[i]:
            no_free.append(list(range(size - 2, size - 2 + num_indices)))
        size += num_indices
    nr = size - 2
    res_gens = [h for h in res_gens if h._array_form != id_af]
    if sym is None or not no_free:
        if not res_gens:
            res_gens = [_af_new(id_af)]
        return (size, res_base, res_gens)
    base_comm = []
    for i in range(len(no_free) - 1):
        ind1 = no_free[i]
        ind2 = no_free[i + 1]
        a = list(range(ind1[0]))
        a.extend(ind2)
        a.extend(ind1)
        base_comm.append(ind1[0])
        a.extend(list(range(ind2[-1] + 1, nr)))
        if sym == 0:
            a.extend([nr, nr + 1])
        else:
            a.extend([nr + 1, nr])
        res_gens.append(_af_new(a))
    res_base = list(res_base)
    for i in base_comm:
        if i not in res_base:
            res_base.append(i)
    res_base.sort()
    if not res_gens:
        res_gens = [_af_new(id_af)]
    return (size, res_base, res_gens)

def gens_products(*v):
    (res_size, res_base, res_gens) = tensor_gens(*v[0])
    for i in range(1, len(v)):
        (size, base, gens) = tensor_gens(*v[i])
        (res_base, res_gens) = bsgs_direct_product(res_base, res_gens, base, gens, 1)
    res_size = res_gens[0].size
    id_af = list(range(res_size))
    res_gens = [h for h in res_gens if h != id_af]
    if not res_gens:
        res_gens = [id_af]
    return (res_size, res_base, res_gens)
