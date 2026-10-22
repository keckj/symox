"""Base and strong generating set utilities (ported from the reference util module)."""
from symox.combinatorics.permutations import Permutation, _af_invert, _af_rmul
def isprime(n):
    from symox.ntheory import isprime as f;return f(n)
rmul = Permutation.rmul
_af_new = Permutation._af_new

def _base_ordering(base, degree):
    base_len = len(base)
    ordering = [0] * degree
    for i in range(base_len):
        ordering[base[i]] = i
    current = base_len
    for i in range(degree):
        if i not in base:
            ordering[i] = current
            current += 1
    return ordering

def _check_cycles_alt_sym(perm):
    n = perm.size
    af = perm.array_form
    current_len = 0
    total_len = 0
    used = set()
    for i in range(n // 2):
        if i not in used and i < n // 2 - total_len:
            current_len = 1
            used.add(i)
            j = i
            while af[j] != i:
                current_len += 1
                j = af[j]
                used.add(j)
            total_len += current_len
            if current_len > n // 2 and current_len < n - 2 and isprime(current_len):
                return True
    return False

def _distribute_gens_by_base(base, gens):
    base_len = len(base)
    degree = gens[0].size
    stabs = [[] for _ in range(base_len)]
    max_stab_index = 0
    for gen in gens:
        j = 0
        while j < base_len - 1 and gen._array_form[base[j]] == base[j]:
            j += 1
        if j > max_stab_index:
            max_stab_index = j
        for k in range(j + 1):
            stabs[k].append(gen)
    for i in range(max_stab_index + 1, base_len):
        stabs[i].append(_af_new(list(range(degree))))
    return stabs

def _handle_precomputed_bsgs(base, strong_gens, transversals=None, basic_orbits=None, strong_gens_distr=None):
    if strong_gens_distr is None:
        strong_gens_distr = _distribute_gens_by_base(base, strong_gens)
    if transversals is None:
        if basic_orbits is None:
            (basic_orbits, transversals) = _orbits_transversals_from_bsgs(base, strong_gens_distr)
        else:
            transversals = _orbits_transversals_from_bsgs(base, strong_gens_distr, transversals_only=True)
    elif basic_orbits is None:
        base_len = len(base)
        basic_orbits = [None] * base_len
        for i in range(base_len):
            basic_orbits[i] = list(transversals[i].keys())
    return (transversals, basic_orbits, strong_gens_distr)

def _orbits_transversals_from_bsgs(base, strong_gens_distr, transversals_only=False, slp=False):
    from symox.combinatorics.perm_groups import _orbit_transversal
    base_len = len(base)
    degree = strong_gens_distr[0][0].size
    transversals = [None] * base_len
    slps = [None] * base_len
    if transversals_only is False:
        basic_orbits = [None] * base_len
    for i in range(base_len):
        (transversals[i], slps[i]) = _orbit_transversal(degree, strong_gens_distr[i], base[i], pairs=True, slp=True)
        transversals[i] = dict(transversals[i])
        if transversals_only is False:
            basic_orbits[i] = list(transversals[i].keys())
    if transversals_only:
        return transversals
    else:
        if not slp:
            return (basic_orbits, transversals)
        return (basic_orbits, transversals, slps)

def _remove_gens(base, strong_gens, basic_orbits=None, strong_gens_distr=None):
    from symox.combinatorics.perm_groups import _orbit
    base_len = len(base)
    degree = strong_gens[0].size
    if strong_gens_distr is None:
        strong_gens_distr = _distribute_gens_by_base(base, strong_gens)
    if basic_orbits is None:
        basic_orbits = []
        for i in range(base_len):
            basic_orbit = _orbit(degree, strong_gens_distr[i], base[i])
            basic_orbits.append(basic_orbit)
    strong_gens_distr.append([])
    res = strong_gens[:]
    for i in range(base_len - 1, -1, -1):
        gens_copy = strong_gens_distr[i][:]
        for gen in strong_gens_distr[i]:
            if gen not in strong_gens_distr[i + 1]:
                temp_gens = gens_copy[:]
                temp_gens.remove(gen)
                if temp_gens == []:
                    continue
                temp_orbit = _orbit(degree, temp_gens, base[i])
                if temp_orbit == basic_orbits[i]:
                    gens_copy.remove(gen)
                    res.remove(gen)
    return res

def _strip(g, base, orbits, transversals):
    h = g._array_form
    base_len = len(base)
    for i in range(base_len):
        beta = h[base[i]]
        if beta == base[i]:
            continue
        if beta not in orbits[i]:
            return (_af_new(h), i + 1)
        u = transversals[i][beta]._array_form
        h = _af_rmul(_af_invert(u), h)
    return (_af_new(h), base_len + 1)

def _strip_af(h, base, orbits, transversals, j, slp=[], slps={}):
    base_len = len(base)
    for i in range(j + 1, base_len):
        beta = h[base[i]]
        if beta == base[i]:
            continue
        if beta not in orbits[i]:
            if not slp:
                return (h, i + 1)
            return (h, i + 1, slp)
        u = transversals[i][beta]
        if h == u:
            if not slp:
                return (False, base_len + 1)
            return (False, base_len + 1, slp)
        h = _af_rmul(_af_invert(u), h)
        if slp:
            u_slp = slps[i][beta][:]
            u_slp.reverse()
            u_slp = [(i, (g,)) for g in u_slp]
            slp = u_slp + slp
    if not slp:
        return (h, base_len + 1)
    return (h, base_len + 1, slp)

def _strong_gens_from_distr(strong_gens_distr):
    if len(strong_gens_distr) == 1:
        return strong_gens_distr[0][:]
    else:
        result = strong_gens_distr[0]
        for gen in strong_gens_distr[1]:
            if gen not in result:
                result.append(gen)
        return result
