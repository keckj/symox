"""Testing helpers for groups (ported from the reference testutil module)."""
from symox.combinatorics import Permutation
from symox.combinatorics.util import _distribute_gens_by_base
rmul = Permutation.rmul

def _cmp_perm_lists(first, second):
    return {tuple(a) for a in first} == {tuple(a) for a in second}

def _naive_list_centralizer(self, other, af=False):
    from symox.combinatorics.perm_groups import PermutationGroup
    from symox.combinatorics.permutations import _af_commutes_with
    if hasattr(other, 'generators'):
        elements = list(self.generate_dimino(af=True))
        gens = [x._array_form for x in other.generators]
        commutes_with_gens = lambda x: all((_af_commutes_with(x, gen) for gen in gens))
        centralizer_list = []
        if not af:
            for element in elements:
                if commutes_with_gens(element):
                    centralizer_list.append(Permutation._af_new(element))
        else:
            for element in elements:
                if commutes_with_gens(element):
                    centralizer_list.append(element)
        return centralizer_list
    elif hasattr(other, 'getitem'):
        return _naive_list_centralizer(self, PermutationGroup(other), af)
    elif hasattr(other, 'array_form'):
        return _naive_list_centralizer(self, PermutationGroup([other]), af)

def _verify_bsgs(group, base, gens):
    from symox.combinatorics.perm_groups import PermutationGroup
    strong_gens_distr = _distribute_gens_by_base(base, gens)
    current_stabilizer = group
    for i in range(len(base)):
        candidate = PermutationGroup(strong_gens_distr[i])
        if current_stabilizer.order() != candidate.order():
            return False
        current_stabilizer = current_stabilizer.stabilizer(base[i])
    if current_stabilizer.order() != 1:
        return False
    return True

def _verify_centralizer(group, arg, centr=None):
    if centr is None:
        centr = group.centralizer(arg)
    centr_list = list(centr.generate_dimino(af=True))
    centr_list_naive = _naive_list_centralizer(group, arg, af=True)
    return _cmp_perm_lists(centr_list, centr_list_naive)

def _verify_normal_closure(group, arg, closure=None):
    from symox.combinatorics.perm_groups import PermutationGroup
    if closure is None:
        closure = group.normal_closure(arg)
    conjugates = set()
    if hasattr(arg, 'generators'):
        subgr_gens = arg.generators
    elif hasattr(arg, '__getitem__'):
        subgr_gens = arg
    elif hasattr(arg, 'array_form'):
        subgr_gens = [arg]
    for el in group.generate_dimino():
        conjugates.update((gen ^ el for gen in subgr_gens))
    naive_closure = PermutationGroup(list(conjugates))
    return closure.is_subgroup(naive_closure)

def canonicalize_naive(g, dummies, sym, *v):
    from symox.combinatorics.perm_groups import PermutationGroup
    from symox.combinatorics.tensor_can import gens_products, dummy_sgs
    from symox.combinatorics.permutations import _af_rmul
    v1 = []
    for i in range(len(v)):
        (base_i, gens_i, n_i, sym_i) = v[i]
        v1.append((base_i, gens_i, [[]] * n_i, sym_i))
    (size, sbase, sgens) = gens_products(*v1)
    dgens = dummy_sgs(dummies, sym, size - 2)
    if isinstance(sym, int):
        num_types = 1
        dummies = [dummies]
        sym = [sym]
    else:
        num_types = len(sym)
    dgens = []
    for i in range(num_types):
        dgens.extend(dummy_sgs(dummies[i], sym[i], size - 2))
    S = PermutationGroup(sgens)
    D = PermutationGroup([Permutation(x) for x in dgens])
    dlist = list(D.generate(af=True))
    g = g.array_form
    st = set()
    for s in S.generate(af=True):
        h = _af_rmul(g, s)
        for d in dlist:
            q = tuple(_af_rmul(d, h))
            st.add(q)
    a = list(st)
    a.sort()
    prev = (0,) * size
    for h in a:
        if h[:-2] == prev[:-2]:
            if h[-1] != prev[-1]:
                return 0
        prev = h
    return list(a[0])

def graph_certificate(gr):
    from symox.combinatorics.permutations import _af_invert
    from symox.combinatorics.tensor_can import get_symmetric_group_sgs, canonicalize
    if sorted(gr) != list(range(len(gr))):
        # Extension: the reference needs the vertices 0, ..., n - 1 and
        # fails otherwise; other hashable labels are renumbered in sorted
        # order, which does not change the certificate of the graph.
        label = {v: i for (i, v) in enumerate(sorted(gr))}
        if any(w not in label for neigh in gr.values() for w in neigh):
            raise ValueError('all adjacency vertices must be declared')
        gr = {label[v]: [label[w] for w in neigh] for (v, neigh) in gr.items()}
    items = list(gr.items())
    items.sort(key=lambda x: len(x[1]), reverse=True)
    pvert = [x[0] for x in items]
    pvert = _af_invert(pvert)
    num_indices = 0
    for (v, neigh) in items:
        num_indices += len(neigh)
    vertices = [[] for i in items]
    i = 0
    for (v, neigh) in items:
        for v2 in neigh:
            if pvert[v] < pvert[v2]:
                vertices[pvert[v]].append(i)
                vertices[pvert[v2]].append(i + 1)
                i += 2
    g = []
    for v in vertices:
        g.extend(v)
    assert len(g) == num_indices
    g += [num_indices, num_indices + 1]
    size = num_indices + 2
    assert sorted(g) == list(range(size))
    g = Permutation(g)
    vlen = [0] * (len(vertices[0]) + 1)
    for neigh in vertices:
        vlen[len(neigh)] += 1
    v = []
    for i in range(len(vlen)):
        n = vlen[i]
        if n:
            (base, gens) = get_symmetric_group_sgs(i)
            v.append((base, gens, n, 0))
    v.reverse()
    dummies = list(range(num_indices))
    can = canonicalize(g, dummies, 0, *v)
    return can

rmul = Permutation.rmul
