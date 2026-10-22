"""Exact finite permutation groups.

Order and membership use the native incremental stabilizer chain
(``_schreier``). The base, strong generating set, transversals and every
structural algorithm (centralizers, normal closures, series, blocks, Sylow
subgroups, subgroup search, randomized Schreier-Sims, presentations) are ported
from the reference ``perm_groups`` module so that deterministic outputs agree.
"""
from collections import deque
from itertools import chain, product
from math import factorial as _factorial, log, prod
import math
import operator
import random

from .permutations import (Permutation, Cycle, _af_commutes_with, _af_invert,
    _af_rmul, _af_rmuln, _af_pow, _af_new)
from ._schreier import StabilizerChain, multiply, inverse
from .util import (_check_cycles_alt_sym, _distribute_gens_by_base,
    _orbits_transversals_from_bsgs, _handle_precomputed_bsgs, _base_ordering,
    _strong_gens_from_distr, _strip, _strip_af)
# The package is imported while symox initializes, so numeric helpers resolve
# lazily on first use.
def _randrange(seed=None):
    from ..core.random import _randrange as f;return f(seed)
def randrange(*args):
    from ..core import random as r;return r.randrange(*args)
def choice(seq):
    from ..core import random as r;return r.choice(seq)
def factorial(n):
    from ..functions.combinatorial.factorials import factorial as f;return f(n)
def factorint(n,*args,**kwargs):
    from ..ntheory.factor_ import factorint as f;return f(n,*args,**kwargs)
def multiplicity(p,n):
    from ..ntheory.factor_ import multiplicity as f;return f(p,n)
def primefactors(n,*args,**kwargs):
    from ..ntheory import primefactors as f;return f(n,*args,**kwargs)
def isprime(n):
    from ..ntheory.primetest import isprime as f;return f(n)
def _sympify(a):
    from ..core.sympify import _sympify as f;return f(a)
def _primes():
    from ..ntheory import sieve;return iter(sieve)


rmul = Permutation.rmul_with_af


class _ReferencePermutationGroupMethods:

    def __contains__(self, i):
        if not isinstance(i, Permutation):
            raise TypeError('A PermutationGroup contains only Permutations as elements, not elements of type %s' % type(i))
        return self.contains(i)

    def equals(self, other):
        if not isinstance(other, PermutationGroup):
            return False
        set_self_gens = set(self.generators)
        set_other_gens = set(other.generators)
        if set_self_gens == set_other_gens:
            return True
        for gen1 in set_self_gens:
            if not other.contains(gen1):
                return False
        for gen2 in set_other_gens:
            if not self.contains(gen2):
                return False
        return True

    def __mul__(self, other):
        if isinstance(other, Permutation):
            from .coset import Coset
            return Coset(other, self, dir='+')
        gens1 = [perm._array_form for perm in self.generators]
        gens2 = [perm._array_form for perm in other.generators]
        n1 = self._degree
        n2 = other._degree
        start = list(range(n1))
        end = list(range(n1, n1 + n2))
        for i in range(len(gens2)):
            gens2[i] = [x + n1 for x in gens2[i]]
        gens2 = [start + gen for gen in gens2]
        gens1 = [gen + end for gen in gens1]
        together = gens1 + gens2
        gens = [_af_new(x) for x in together]
        return PermutationGroup(gens)

    def _random_pr_init(self, r, n, _random_prec_n=None):
        deg = self.degree
        random_gens = [x._array_form for x in self.generators]
        k = len(random_gens)
        if k < r:
            for i in range(k, r):
                random_gens.append(random_gens[i - k])
        acc = list(range(deg))
        random_gens.append(acc)
        self._random_gens = random_gens
        if _random_prec_n is None:
            for i in range(n):
                self.random_pr()
        else:
            for i in range(n):
                self.random_pr(_random_prec=_random_prec_n[i])

    def _union_find_merge(self, first, second, ranks, parents, not_rep):
        rep_first = self._union_find_rep(first, parents)
        rep_second = self._union_find_rep(second, parents)
        if rep_first != rep_second:
            if ranks[rep_first] >= ranks[rep_second]:
                (new_1, new_2) = (rep_first, rep_second)
            else:
                (new_1, new_2) = (rep_second, rep_first)
            total_rank = ranks[new_1] + ranks[new_2]
            if total_rank > self.max_div:
                return -1
            parents[new_2] = new_1
            ranks[new_1] = total_rank
            not_rep.append(new_2)
            return 1
        return 0

    def _union_find_rep(self, num, parents):
        (rep, parent) = (num, parents[num])
        while parent != rep:
            rep = parent
            parent = parents[rep]
        (temp, parent) = (num, parents[num])
        while parent != rep:
            parents[temp] = rep
            temp = parent
            parent = parents[temp]
        return rep

    @property
    def base(self):
        if self._base == []:
            self.schreier_sims()
        return self._base

    def baseswap(self, base, strong_gens, pos, randomized=False, transversals=None, basic_orbits=None, strong_gens_distr=None):
        (transversals, basic_orbits, strong_gens_distr) = _handle_precomputed_bsgs(base, strong_gens, transversals, basic_orbits, strong_gens_distr)
        base_len = len(base)
        degree = self.degree
        size = len(basic_orbits[pos]) * len(basic_orbits[pos + 1]) // len(_orbit(degree, strong_gens_distr[pos], base[pos + 1]))
        if pos + 2 > base_len - 1:
            T = []
        else:
            T = strong_gens_distr[pos + 2][:]
        if randomized is True:
            stab_pos = PermutationGroup(strong_gens_distr[pos])
            schreier_vector = stab_pos.schreier_vector(base[pos + 1])
            while len(_orbit(degree, T, base[pos])) != size:
                new = stab_pos.random_stab(base[pos + 1], schreier_vector=schreier_vector)
                T.append(new)
        else:
            Gamma = set(basic_orbits[pos])
            Gamma.remove(base[pos])
            if base[pos + 1] in Gamma:
                Gamma.remove(base[pos + 1])
            while len(_orbit(degree, T, base[pos])) != size:
                gamma = next(iter(Gamma))
                x = transversals[pos][gamma]
                temp = x._array_form.index(base[pos + 1])
                if temp not in basic_orbits[pos + 1]:
                    Gamma = Gamma - _orbit(degree, T, gamma)
                else:
                    y = transversals[pos + 1][temp]
                    el = rmul(x, y)
                    if el(base[pos]) not in _orbit(degree, T, base[pos]):
                        T.append(el)
                        Gamma = Gamma - _orbit(degree, T, base[pos])
        strong_gens_new_distr = strong_gens_distr[:]
        strong_gens_new_distr[pos + 1] = T
        base_new = base[:]
        (base_new[pos], base_new[pos + 1]) = (base_new[pos + 1], base_new[pos])
        strong_gens_new = _strong_gens_from_distr(strong_gens_new_distr)
        for gen in T:
            if gen not in strong_gens_new:
                strong_gens_new.append(gen)
        return (base_new, strong_gens_new)

    @property
    def basic_orbits(self):
        if self._basic_orbits == []:
            self.schreier_sims()
        return self._basic_orbits

    @property
    def basic_stabilizers(self):
        if self._transversals == []:
            self.schreier_sims()
        strong_gens = self._strong_gens
        base = self._base
        if not base:
            return []
        strong_gens_distr = _distribute_gens_by_base(base, strong_gens)
        basic_stabilizers = []
        for gens in strong_gens_distr:
            basic_stabilizers.append(PermutationGroup(gens))
        return basic_stabilizers

    @property
    def basic_transversals(self):
        if self._transversals == []:
            self.schreier_sims()
        return self._transversals

    def composition_series(self):
        der = self.derived_series()
        if not all((g.is_identity for g in der[-1].generators)):
            raise NotImplementedError('Group should be solvable')
        series = []
        for i in range(len(der) - 1):
            H = der[i + 1]
            up_seg = []
            for g in der[i].generators:
                K = PermutationGroup([g] + H.generators)
                order = K.order() // H.order()
                down_seg = []
                for (p, e) in factorint(order).items():
                    for _ in range(e):
                        down_seg.append(PermutationGroup([g] + H.generators))
                        g = g ** p
                up_seg = down_seg + up_seg
                H = K
            up_seg[0] = der[i]
            series.extend(up_seg)
        series.append(der[-1])
        return series

    def coset_transversal(self, H):
        if not H.is_subgroup(self):
            raise ValueError('The argument must be a subgroup')
        if H.order() == 1:
            return self.elements
        self._schreier_sims(base=H.base)
        base = self.base
        base_ordering = _base_ordering(base, self.degree)
        identity = Permutation(self.degree - 1)
        transversals = self.basic_transversals[:]
        for (l, t) in enumerate(transversals):
            transversals[l] = sorted(t.values(), key=lambda x: base_ordering[base[l] ^ x])
        orbits = H.basic_orbits
        h_stabs = H.basic_stabilizers
        g_stabs = self.basic_stabilizers
        indices = [x.order() // y.order() for (x, y) in zip(g_stabs, h_stabs)]
        if len(g_stabs) > len(h_stabs):
            T = g_stabs[len(h_stabs)].elements
        else:
            T = [identity]
        l = len(h_stabs) - 1
        t_len = len(T)
        while l > -1:
            T_next = []
            for u in transversals[l]:
                if u == identity:
                    continue
                b = base_ordering[base[l] ^ u]
                for t in T:
                    p = t * u
                    if all((base_ordering[h ^ p] >= b for h in orbits[l])):
                        T_next.append(p)
                    if t_len + len(T_next) == indices[l]:
                        break
                if t_len + len(T_next) == indices[l]:
                    break
            T += T_next
            t_len += len(T_next)
            l -= 1
        T.remove(identity)
        T = [identity] + T
        return T

    def _coset_representative(self, g, H):
        if H.order() == 1:
            return g
        if not self.base[:len(H.base)] == H.base:
            self._schreier_sims(base=H.base)
        orbits = H.basic_orbits[:]
        h_transversals = [list(_.values()) for _ in H.basic_transversals]
        transversals = [list(_.values()) for _ in self.basic_transversals]
        base = self.base
        base_ordering = _base_ordering(base, self.degree)

        def step(l, x):
            gamma = min(orbits[l], key=lambda y: base_ordering[y ^ x])
            i = [base[l] ^ h for h in h_transversals[l]].index(gamma)
            x = h_transversals[l][i] * x
            if l < len(orbits) - 1:
                for u in transversals[l]:
                    if base[l] ^ u == base[l] ^ x:
                        break
                x = step(l + 1, x * u ** (-1)) * u
            return x
        return step(0, g)

    def coset_table(self, H):
        if not H.is_subgroup(self):
            raise ValueError('The argument must be a subgroup')
        T = self.coset_transversal(H)
        n = len(T)
        A = list(chain.from_iterable(((gen, gen ** (-1)) for gen in self.generators)))
        table = []
        for i in range(n):
            row = [self._coset_representative(T[i] * x, H) for x in A]
            row = [T.index(r) for r in row]
            table.append(row)
        A = range(len(A))
        gamma = 1
        for (alpha, a) in product(range(n), A):
            beta = table[alpha][a]
            if beta >= gamma:
                if beta > gamma:
                    for x in A:
                        z = table[gamma][x]
                        table[gamma][x] = table[beta][x]
                        table[beta][x] = z
                        for i in range(n):
                            if table[i][x] == beta:
                                table[i][x] = gamma
                            elif table[i][x] == gamma:
                                table[i][x] = beta
                gamma += 1
            if gamma >= n - 1:
                return table

    def center(self):
        if not self._center:
            self._center = self.centralizer(self)
        return self._center

    def centralizer(self, other):
        if hasattr(other, 'generators'):
            if other.is_trivial or self.is_trivial:
                return self
            degree = self.degree
            identity = _af_new(list(range(degree)))
            orbits = other.orbits()
            num_orbits = len(orbits)
            orbits.sort(key=lambda x: -len(x))
            long_base = []
            orbit_reps = [None] * num_orbits
            orbit_reps_indices = [None] * num_orbits
            orbit_descr = [None] * degree
            for i in range(num_orbits):
                orbit = list(orbits[i])
                orbit_reps[i] = orbit[0]
                orbit_reps_indices[i] = len(long_base)
                for point in orbit:
                    orbit_descr[point] = i
                long_base = long_base + orbit
            (base, strong_gens) = self.schreier_sims_incremental(base=long_base)
            strong_gens_distr = _distribute_gens_by_base(base, strong_gens)
            i = 0
            for i in range(len(base)):
                if strong_gens_distr[i] == [identity]:
                    break
            base = base[:i]
            base_len = i
            for j in range(num_orbits):
                if base[base_len - 1] in orbits[j]:
                    break
            rel_orbits = orbits[:j + 1]
            num_rel_orbits = len(rel_orbits)
            transversals = [None] * num_rel_orbits
            for j in range(num_rel_orbits):
                rep = orbit_reps[j]
                transversals[j] = dict(other.orbit_transversal(rep, pairs=True))
            trivial_test = lambda x: True
            tests = [None] * base_len
            for l in range(base_len):
                if base[l] in orbit_reps:
                    tests[l] = trivial_test
                else:

                    def test(computed_words, l=l):
                        g = computed_words[l]
                        rep_orb_index = orbit_descr[base[l]]
                        rep = orbit_reps[rep_orb_index]
                        im = g._array_form[base[l]]
                        im_rep = g._array_form[rep]
                        tr_el = transversals[rep_orb_index][base[l]]
                        return im == tr_el._array_form[im_rep]
                    tests[l] = test

            def prop(g):
                return [rmul(g, gen) for gen in other.generators] == [rmul(gen, g) for gen in other.generators]
            return self.subgroup_search(prop, base=base, strong_gens=strong_gens, tests=tests)
        elif hasattr(other, '__getitem__'):
            gens = list(other)
            return self.centralizer(PermutationGroup(gens))
        elif hasattr(other, 'array_form'):
            return self.centralizer(PermutationGroup([other]))

    def commutator(self, G, H):
        ggens = G.generators
        hgens = H.generators
        commutators = []
        for ggen in ggens:
            for hgen in hgens:
                commutator = rmul(hgen, ggen, ~hgen, ~ggen)
                if commutator not in commutators:
                    commutators.append(commutator)
        res = self.normal_closure(commutators)
        return res

    def coset_factor(self, g, factor_index=False):
        if isinstance(g, (Cycle, Permutation)):
            g = g.list()
        if len(g) != self._degree:
            raise ValueError('g should be the same size as permutations of G')
        I = list(range(self._degree))
        basic_orbits = self.basic_orbits
        transversals = self._transversals
        factors = []
        base = self.base
        h = g
        for i in range(len(base)):
            beta = h[base[i]]
            if beta == base[i]:
                factors.append(beta)
                continue
            if beta not in basic_orbits[i]:
                return []
            u = transversals[i][beta]._array_form
            h = _af_rmul(_af_invert(u), h)
            factors.append(beta)
        if h != I:
            return []
        if factor_index:
            return factors
        tr = self.basic_transversals
        factors = [tr[i][factors[i]] for i in range(len(base))]
        return factors

    def generator_product(self, g, original=False):
        product = []
        if g.is_identity:
            return []
        if g in self.strong_gens:
            if not original or g in self.generators:
                return [g]
            else:
                slp = self._strong_gens_slp[g]
                for s in slp:
                    product.extend(self.generator_product(s, original=True))
                return product
        elif g ** (-1) in self.strong_gens:
            g = g ** (-1)
            if not original or g in self.generators:
                return [g ** (-1)]
            else:
                slp = self._strong_gens_slp[g]
                for s in slp:
                    product.extend(self.generator_product(s, original=True))
                l = len(product)
                product = [product[l - i - 1] ** (-1) for i in range(l)]
                return product
        f = self.coset_factor(g, True)
        for (i, j) in enumerate(f):
            slp = self._transversal_slp[i][j]
            for s in slp:
                if not original:
                    product.append(self.strong_gens[s])
                else:
                    s = self.strong_gens[s]
                    product.extend(self.generator_product(s, original=True))
        return product

    def coset_rank(self, g):
        factors = self.coset_factor(g, True)
        if not factors:
            return None
        rank = 0
        b = 1
        transversals = self._transversals
        base = self._base
        basic_orbits = self._basic_orbits
        for i in range(len(base)):
            k = factors[i]
            j = basic_orbits[i].index(k)
            rank += b * j
            b = b * len(transversals[i])
        return rank

    def coset_unrank(self, rank, af=False):
        if rank < 0 or rank >= self.order():
            return None
        base = self.base
        transversals = self.basic_transversals
        basic_orbits = self.basic_orbits
        m = len(base)
        v = [0] * m
        for i in range(m):
            (rank, c) = divmod(rank, len(transversals[i]))
            v[i] = basic_orbits[i][c]
        a = [transversals[i][v[i]]._array_form for i in range(m)]
        h = _af_rmuln(*a)
        if af:
            return h
        else:
            return _af_new(h)

    def derived_series(self):
        res = [self]
        current = self
        nxt = self.derived_subgroup()
        while not current.is_subgroup(nxt):
            res.append(nxt)
            current = nxt
            nxt = nxt.derived_subgroup()
        return res

    def derived_subgroup(self):
        r = self._r
        gens = [p._array_form for p in self.generators]
        set_commutators = set()
        degree = self._degree
        rng = list(range(degree))
        for i in range(r):
            for j in range(r):
                p1 = gens[i]
                p2 = gens[j]
                c = list(range(degree))
                for k in rng:
                    c[p2[p1[k]]] = p1[p2[k]]
                ct = tuple(c)
                if ct not in set_commutators:
                    set_commutators.add(ct)
        cms = [_af_new(p) for p in set_commutators]
        G2 = self.normal_closure(cms)
        return G2

    def generate(self, method='coset', af=False):
        if method == 'coset':
            return self.generate_schreier_sims(af)
        elif method == 'dimino':
            return self.generate_dimino(af)
        else:
            raise NotImplementedError('No generation defined for %s' % method)

    def generate_dimino(self, af=False):
        idn = list(range(self.degree))
        order = 0
        element_list = [idn]
        set_element_list = {tuple(idn)}
        if af:
            yield idn
        else:
            yield _af_new(idn)
        gens = [p._array_form for p in self.generators]
        for i in range(len(gens)):
            D = element_list.copy()
            N = [idn]
            while N:
                A = N
                N = []
                for a in A:
                    for g in gens[:i + 1]:
                        ag = _af_rmul(a, g)
                        if tuple(ag) not in set_element_list:
                            for d in D:
                                order += 1
                                ap = _af_rmul(d, ag)
                                if af:
                                    yield ap
                                else:
                                    p = _af_new(ap)
                                    yield p
                                element_list.append(ap)
                                set_element_list.add(tuple(ap))
                                N.append(ap)
        self._order = len(element_list)

    def generate_schreier_sims(self, af=False):
        n = self._degree
        u = self.basic_transversals
        basic_orbits = self._basic_orbits
        if len(u) == 0:
            for x in self.generators:
                if af:
                    yield x._array_form
                else:
                    yield x
            return
        if len(u) == 1:
            for i in basic_orbits[0]:
                if af:
                    yield u[0][i]._array_form
                else:
                    yield u[0][i]
            return
        u = list(reversed(u))
        basic_orbits = basic_orbits[::-1]
        stg = [list(range(n))]
        posmax = [len(x) for x in u]
        n1 = len(posmax) - 1
        pos = [0] * n1
        h = 0
        while 1:
            if pos[h] >= posmax[h]:
                if h == 0:
                    return
                pos[h] = 0
                h -= 1
                stg.pop()
                continue
            p = _af_rmul(u[h][basic_orbits[h][pos[h]]]._array_form, stg[-1])
            pos[h] += 1
            stg.append(p)
            h += 1
            if h == n1:
                if af:
                    for i in basic_orbits[-1]:
                        p = _af_rmul(u[-1][i]._array_form, stg[-1])
                        yield p
                else:
                    for i in basic_orbits[-1]:
                        p = _af_rmul(u[-1][i]._array_form, stg[-1])
                        p1 = _af_new(p)
                        yield p1
                stg.pop()
                h -= 1

    @property
    def is_perfect(self):
        if self._is_perfect is None:
            self._is_perfect = self.equals(self.derived_subgroup())
        return self._is_perfect

    @property
    def is_abelian(self):
        if self._is_abelian is not None:
            return self._is_abelian
        self._is_abelian = True
        gens = [p._array_form for p in self.generators]
        for x in gens:
            for y in gens:
                if y <= x:
                    continue
                if not _af_commutes_with(x, y):
                    self._is_abelian = False
                    return False
        return True

    def abelian_invariants(self):
        if self.is_trivial:
            return []
        gns = self.generators
        inv = []
        G = self
        H = G.derived_subgroup()
        Hgens = H.generators
        for p in primefactors(G.order()):
            ranks = []
            while True:
                pows = []
                for g in gns:
                    elm = g ** p
                    if not H.contains(elm):
                        pows.append(elm)
                K = PermutationGroup(Hgens + pows) if pows else H
                r = G.order() // K.order()
                G = K
                gns = pows
                if r == 1:
                    break
                ranks.append(multiplicity(p, r))
            if ranks:
                pows = [1] * ranks[0]
                for i in ranks:
                    for j in range(i):
                        pows[j] = pows[j] * p
                inv.extend(pows)
        inv.sort()
        return inv

    def is_elementary(self, p):
        return self.is_abelian and all((g.order() == p for g in self.generators))

    def _eval_is_alt_sym_naive(self, only_sym=False, only_alt=False):
        if only_sym and only_alt:
            raise ValueError('Both {} and {} cannot be set to True'.format(only_sym, only_alt))
        n = self.degree
        sym_order = _factorial(n)
        order = self.order()
        if order == sym_order:
            self._is_sym = True
            self._is_alt = False
            return not only_alt
        if 2 * order == sym_order:
            self._is_sym = False
            self._is_alt = True
            return not only_sym
        return False

    def _eval_is_alt_sym_monte_carlo(self, eps=0.05, perms=None):
        if perms is None:
            n = self.degree
            if n < 17:
                c_n = 0.34
            else:
                c_n = 0.57
            d_n = c_n * log(2) / log(n)
            N_eps = int(-log(eps) / d_n)
            perms = (self.random_pr() for i in range(N_eps))
            return self._eval_is_alt_sym_monte_carlo(perms=perms)
        for perm in perms:
            if _check_cycles_alt_sym(perm):
                return True
        return False

    def is_alt_sym(self, eps=0.05, _random_prec=None):
        if _random_prec is not None:
            N_eps = _random_prec['N_eps']
            perms = (_random_prec[i] for i in range(N_eps))
            return self._eval_is_alt_sym_monte_carlo(perms=perms)
        if self._is_sym or self._is_alt:
            return True
        if self._is_sym is False and self._is_alt is False:
            return False
        n = self.degree
        if n < 8:
            return self._eval_is_alt_sym_naive()
        elif self.is_transitive():
            return self._eval_is_alt_sym_monte_carlo(eps=eps)
        (self._is_sym, self._is_alt) = (False, False)
        return False

    @property
    def is_nilpotent(self):
        if self._is_nilpotent is None:
            lcs = self.lower_central_series()
            terminator = lcs[len(lcs) - 1]
            gens = terminator.generators
            degree = self.degree
            identity = _af_new(list(range(degree)))
            if all((g == identity for g in gens)):
                self._is_solvable = True
                self._is_nilpotent = True
                return True
            else:
                self._is_nilpotent = False
                return False
        else:
            return self._is_nilpotent

    def is_normal(self, gr, strict=True):
        if not self.is_subgroup(gr, strict=strict):
            return False
        d_self = self.degree
        d_gr = gr.degree
        if self.is_trivial and (d_self == d_gr or not strict):
            return True
        if self._is_abelian:
            return True
        new_self = self.copy()
        if not strict and d_self != d_gr:
            if d_self < d_gr:
                new_self = PermGroup(new_self.generators + [Permutation(d_gr - 1)])
            else:
                gr = PermGroup(gr.generators + [Permutation(d_self - 1)])
        gens2 = [p._array_form for p in new_self.generators]
        gens1 = [p._array_form for p in gr.generators]
        for g1 in gens1:
            for g2 in gens2:
                p = _af_rmuln(g1, g2, _af_invert(g1))
                if not new_self.coset_factor(p, True):
                    return False
        return True

    def is_primitive(self, randomized=True):
        if self._is_primitive is not None:
            return self._is_primitive
        if self.is_transitive() is False:
            return False
        if randomized:
            random_stab_gens = []
            v = self.schreier_vector(0)
            for _ in range(len(self)):
                random_stab_gens.append(self.random_stab(0, v))
            stab = PermutationGroup(random_stab_gens)
        else:
            stab = self.stabilizer(0)
        orbits = stab.orbits()
        for orb in orbits:
            x = orb.pop()
            if x != 0 and any((e != 0 for e in self.minimal_block([0, x]))):
                self._is_primitive = False
                return False
        self._is_primitive = True
        return True

    def minimal_blocks(self, randomized=True):

        def _number_blocks(blocks):
            n = len(blocks)
            appeared = {}
            m = 0
            b = [None] * n
            for i in range(n):
                if blocks[i] not in appeared:
                    appeared[blocks[i]] = m
                    b[i] = m
                    m += 1
                else:
                    b[i] = appeared[blocks[i]]
            return (tuple(b), m)
        if not self.is_transitive():
            return False
        blocks = []
        num_blocks = []
        rep_blocks = []
        if randomized:
            random_stab_gens = []
            v = self.schreier_vector(0)
            for i in range(len(self)):
                random_stab_gens.append(self.random_stab(0, v))
            stab = PermutationGroup(random_stab_gens)
        else:
            stab = self.stabilizer(0)
        orbits = stab.orbits()
        for orb in orbits:
            x = orb.pop()
            if x != 0:
                block = self.minimal_block([0, x])
                (num_block, _) = _number_blocks(block)
                rep = {j for j in range(self.degree) if num_block[j] == 0}
                minimal = True
                blocks_remove_mask = [False] * len(blocks)
                for (i, r) in enumerate(rep_blocks):
                    if len(r) > len(rep) and rep.issubset(r):
                        blocks_remove_mask[i] = True
                    elif len(r) < len(rep) and r.issubset(rep):
                        minimal = False
                        break
                blocks = [b for (i, b) in enumerate(blocks) if not blocks_remove_mask[i]]
                num_blocks = [n for (i, n) in enumerate(num_blocks) if not blocks_remove_mask[i]]
                rep_blocks = [r for (i, r) in enumerate(rep_blocks) if not blocks_remove_mask[i]]
                if minimal and num_block not in num_blocks:
                    blocks.append(block)
                    num_blocks.append(num_block)
                    rep_blocks.append(rep)
        return blocks

    @property
    def is_solvable(self):
        if self._is_solvable is None:
            if self.order() % 2 != 0:
                return True
            ds = self.derived_series()
            terminator = ds[len(ds) - 1]
            gens = terminator.generators
            degree = self.degree
            identity = _af_new(list(range(degree)))
            if all((g == identity for g in gens)):
                self._is_solvable = True
                return True
            else:
                self._is_solvable = False
                return False
        else:
            return self._is_solvable

    def is_subgroup(self, G, strict=True):
        from .named_groups import SymmetricPermutationGroup
        if isinstance(G, SymmetricPermutationGroup):
            if self.degree != G.degree:
                return False
            return True
        if not isinstance(G, PermutationGroup):
            return False
        if self == G or self.generators[0] == Permutation():
            return True
        if G.order() % self.order() != 0:
            return False
        if self.degree == G.degree or (self.degree < G.degree and (not strict)):
            gens = self.generators
        else:
            return False
        return all((G.contains(g, strict=strict) for g in gens))

    @property
    def is_polycyclic(self):
        return self.is_solvable

    def is_transitive(self, strict=True):
        if self._is_transitive:
            return self._is_transitive
        if strict:
            if self._is_transitive is not None:
                return self._is_transitive
            ans = len(self.orbit(0)) == self.degree
            self._is_transitive = ans
            return ans
        got_orb = False
        for x in self.orbits():
            if len(x) > 1:
                if got_orb:
                    return False
                got_orb = True
        return got_orb

    @property
    def is_trivial(self):
        if self._is_trivial is None:
            self._is_trivial = len(self) == 1 and self[0].is_Identity
        return self._is_trivial

    def lower_central_series(self):
        res = [self]
        current = self
        nxt = self.commutator(self, current)
        while not current.is_subgroup(nxt):
            res.append(nxt)
            current = nxt
            nxt = self.commutator(self, current)
        return res

    @property
    def max_div(self):
        if self._max_div is not None:
            return self._max_div
        n = self.degree
        if n == 1:
            return 1
        for x in _primes():
            if n % x == 0:
                d = n // x
                self._max_div = d
                return d

    def minimal_block(self, points):
        if not self.is_transitive():
            return False
        n = self.degree
        gens = self.generators
        parents = list(range(n))
        ranks = [1] * n
        not_rep = []
        k = len(points)
        if k > self.max_div:
            return [0] * n
        for i in range(k - 1):
            parents[points[i + 1]] = points[0]
            not_rep.append(points[i + 1])
        ranks[points[0]] = k
        i = 0
        len_not_rep = k - 1
        while i < len_not_rep:
            gamma = not_rep[i]
            i += 1
            for gen in gens:
                delta = self._union_find_rep(gamma, parents)
                temp = self._union_find_merge(gen(gamma), gen(delta), ranks, parents, not_rep)
                if temp == -1:
                    return [0] * n
                len_not_rep += temp
        for i in range(n):
            self._union_find_rep(i, parents)
        new_reps = {}
        return [new_reps.setdefault(r, i) for (i, r) in enumerate(parents)]

    def conjugacy_class(self, x):
        new_class = {x}
        last_iteration = new_class
        while len(last_iteration) > 0:
            this_iteration = set()
            for y in last_iteration:
                for s in self.generators:
                    conjugated = s * y * ~s
                    if conjugated not in new_class:
                        this_iteration.add(conjugated)
            new_class.update(last_iteration)
            last_iteration = this_iteration
        return new_class

    def conjugacy_classes(self):
        identity = _af_new(list(range(self.degree)))
        known_elements = {identity}
        classes = [known_elements.copy()]
        for x in self.generate():
            if x not in known_elements:
                new_class = self.conjugacy_class(x)
                classes.append(new_class)
                known_elements.update(new_class)
        return classes

    def normal_closure(self, other, k=10):
        if hasattr(other, 'generators'):
            degree = self.degree
            identity = _af_new(list(range(degree)))
            if all((g == identity for g in other.generators)):
                return other
            Z = PermutationGroup(other.generators[:])
            (base, strong_gens) = Z.schreier_sims_incremental()
            strong_gens_distr = _distribute_gens_by_base(base, strong_gens)
            (basic_orbits, basic_transversals) = _orbits_transversals_from_bsgs(base, strong_gens_distr)
            self._random_pr_init(r=10, n=20)
            _loop = True
            while _loop:
                Z._random_pr_init(r=10, n=10)
                for _ in range(k):
                    g = self.random_pr()
                    h = Z.random_pr()
                    conj = h ^ g
                    res = _strip(conj, base, basic_orbits, basic_transversals)
                    if res[0] != identity or res[1] != len(base) + 1:
                        gens = Z.generators
                        gens.append(conj)
                        Z = PermutationGroup(gens)
                        strong_gens.append(conj)
                        (temp_base, temp_strong_gens) = Z.schreier_sims_incremental(base, strong_gens)
                        (base, strong_gens) = (temp_base, temp_strong_gens)
                        strong_gens_distr = _distribute_gens_by_base(base, strong_gens)
                        (basic_orbits, basic_transversals) = _orbits_transversals_from_bsgs(base, strong_gens_distr)
                _loop = False
                for g in self.generators:
                    for h in Z.generators:
                        conj = h ^ g
                        res = _strip(conj, base, basic_orbits, basic_transversals)
                        if res[0] != identity or res[1] != len(base) + 1:
                            _loop = True
                            break
                    if _loop:
                        break
            return Z
        elif hasattr(other, '__getitem__'):
            return self.normal_closure(PermutationGroup(other))
        elif hasattr(other, 'array_form'):
            return self.normal_closure(PermutationGroup([other]))

    def orbit_rep(self, alpha, beta, schreier_vector=None):
        if schreier_vector is None:
            schreier_vector = self.schreier_vector(alpha)
        if schreier_vector[beta] is None:
            return False
        k = schreier_vector[beta]
        gens = [x._array_form for x in self.generators]
        a = []
        while k != -1:
            a.append(gens[k])
            beta = gens[k].index(beta)
            k = schreier_vector[beta]
        if a:
            return _af_new(_af_rmuln(*a))
        else:
            return _af_new(list(range(self._degree)))

    def orbit_transversal(self, alpha, pairs=False):
        return _orbit_transversal(self._degree, self.generators, alpha, pairs)

    def index(self, H):
        if H.is_subgroup(self):
            return self.order() // H.order()

    @property
    def is_symmetric(self):
        _is_sym = self._is_sym
        if _is_sym is not None:
            return _is_sym
        n = self.degree
        if n >= 8:
            if self.is_transitive():
                _is_alt_sym = self._eval_is_alt_sym_monte_carlo()
                if _is_alt_sym:
                    if any((g.is_odd for g in self.generators)):
                        (self._is_sym, self._is_alt) = (True, False)
                        return True
                    (self._is_sym, self._is_alt) = (False, True)
                    return False
                return self._eval_is_alt_sym_naive(only_sym=True)
            (self._is_sym, self._is_alt) = (False, False)
            return False
        return self._eval_is_alt_sym_naive(only_sym=True)

    @property
    def is_alternating(self):
        _is_alt = self._is_alt
        if _is_alt is not None:
            return _is_alt
        n = self.degree
        if n >= 8:
            if self.is_transitive():
                _is_alt_sym = self._eval_is_alt_sym_monte_carlo()
                if _is_alt_sym:
                    if all((g.is_even for g in self.generators)):
                        (self._is_sym, self._is_alt) = (False, True)
                        return True
                    (self._is_sym, self._is_alt) = (True, False)
                    return False
                return self._eval_is_alt_sym_naive(only_alt=True)
            (self._is_sym, self._is_alt) = (False, False)
            return False
        return self._eval_is_alt_sym_naive(only_alt=True)

    @classmethod
    def _distinct_primes_lemma(cls, primes):
        primes = sorted(primes)
        l = len(primes)
        for i in range(l):
            for j in range(i + 1, l):
                if primes[j] % primes[i] == 1:
                    return None
        return True

    @property
    def is_cyclic(self):
        if self._is_cyclic is not None:
            return self._is_cyclic
        if len(self.generators) == 1:
            self._is_cyclic = True
            self._is_abelian = True
            return True
        if self._is_abelian is False:
            self._is_cyclic = False
            return False
        order = self.order()
        if order < 6:
            self._is_abelian = True
            if order != 4:
                self._is_cyclic = True
                return True
        factors = factorint(order)
        if all((v == 1 for v in factors.values())):
            if self._is_abelian:
                self._is_cyclic = True
                return True
            primes = list(factors.keys())
            if PermutationGroup._distinct_primes_lemma(primes) is True:
                self._is_cyclic = True
                self._is_abelian = True
                return True
        if not self.is_abelian:
            self._is_cyclic = False
            return False
        self._is_cyclic = all((any((g ** (order // p) != self.identity for g in self.generators)) for (p, e) in factors.items() if e > 1))
        return self._is_cyclic

    @property
    def is_dihedral(self):
        if self._is_dihedral is not None:
            return self._is_dihedral
        order = self.order()
        if order % 2 == 1:
            self._is_dihedral = False
            return False
        if order == 2:
            self._is_dihedral = True
            return True
        if order == 4:
            self._is_dihedral = not self.is_cyclic
            return self._is_dihedral
        if self.is_abelian:
            self._is_dihedral = False
            return False
        n = order // 2
        gens = self.generators
        if len(gens) == 2:
            (x, y) = gens
            (a, b) = (x.order(), y.order())
            if a < b:
                (x, y, a, b) = (y, x, b, a)
            if a == 2 == b:
                self._is_dihedral = True
                return True
            if a == n and b == 2 and (y * x * y == ~x):
                self._is_dihedral = True
                return True
        (order_2, order_n) = ([], [])
        for p in self.elements:
            k = p.order()
            if k == 2:
                order_2.append(p)
            elif k == n:
                order_n.append(p)
        if len(order_2) != n + 1 - n % 2:
            self._is_dihedral = False
            return False
        if not order_n:
            self._is_dihedral = False
            return False
        x = order_n[0]
        y = order_2[0]
        if n % 2 == 0 and y == x ** (n // 2):
            y = order_2[1]
        self._is_dihedral = y * x * y == ~x
        return self._is_dihedral

    def pointwise_stabilizer(self, points, incremental=True):
        if incremental:
            (base, strong_gens) = self.schreier_sims_incremental(base=points)
            stab_gens = []
            degree = self.degree
            for gen in strong_gens:
                if [gen(point) for point in points] == points:
                    stab_gens.append(gen)
            if not stab_gens:
                stab_gens = _af_new(list(range(degree)))
            return PermutationGroup(stab_gens)
        else:
            gens = self._generators
            degree = self.degree
            for x in points:
                gens = _stabilizer(degree, gens, x)
        return PermutationGroup(gens)

    def random(self, af=False):
        rank = randrange(self.order())
        return self.coset_unrank(rank, af)

    def random_pr(self, gen_count=11, iterations=50, _random_prec=None):
        if self._random_gens == []:
            self._random_pr_init(gen_count, iterations)
        random_gens = self._random_gens
        r = len(random_gens) - 1
        if _random_prec is None:
            s = randrange(r)
            t = randrange(r - 1)
            if t == s:
                t = r - 1
            x = choice([1, 2])
            e = choice([-1, 1])
        else:
            s = _random_prec['s']
            t = _random_prec['t']
            if t == s:
                t = r - 1
            x = _random_prec['x']
            e = _random_prec['e']
        if x == 1:
            random_gens[s] = _af_rmul(random_gens[s], _af_pow(random_gens[t], e))
            random_gens[r] = _af_rmul(random_gens[r], random_gens[s])
        else:
            random_gens[s] = _af_rmul(_af_pow(random_gens[t], e), random_gens[s])
            random_gens[r] = _af_rmul(random_gens[s], random_gens[r])
        return _af_new(random_gens[r])

    def random_stab(self, alpha, schreier_vector=None, _random_prec=None):
        if schreier_vector is None:
            schreier_vector = self.schreier_vector(alpha)
        if _random_prec is None:
            rand = self.random_pr()
        else:
            rand = _random_prec['rand']
        beta = rand(alpha)
        h = self.orbit_rep(alpha, beta, schreier_vector)
        return rmul(~h, rand)

    def schreier_sims(self):
        if self._transversals:
            return
        self._schreier_sims()
        return

    def _schreier_sims(self, base=None):
        schreier = self.schreier_sims_incremental(base=base, slp_dict=True)
        (base, strong_gens) = schreier[:2]
        self._base = base
        self._strong_gens = strong_gens
        self._strong_gens_slp = schreier[2]
        if not base:
            self._transversals = []
            self._basic_orbits = []
            return
        strong_gens_distr = _distribute_gens_by_base(base, strong_gens)
        (basic_orbits, transversals, slps) = _orbits_transversals_from_bsgs(base, strong_gens_distr, slp=True)
        for (i, slp) in enumerate(slps):
            gens = strong_gens_distr[i]
            for k in slp:
                slp[k] = [strong_gens.index(gens[s]) for s in slp[k]]
        self._transversals = transversals
        self._basic_orbits = [sorted(x) for x in basic_orbits]
        self._transversal_slp = slps

    def schreier_sims_incremental(self, base=None, gens=None, slp_dict=False):
        if base is None:
            base = []
        if gens is None:
            gens = self.generators[:]
        degree = self.degree
        id_af = list(range(degree))
        if len(gens) == 1 and gens[0].is_Identity:
            if slp_dict:
                return (base, gens, {gens[0]: [gens[0]]})
            return (base, gens)
        (_base, _gens) = (base[:], gens[:])
        _gens = [x for x in _gens if not x.is_Identity]
        if not _gens:
            # The reference indexes the empty filtered list (IndexError);
            # the trivial group has the empty strong generating set.
            if slp_dict:
                return (_base, [], {})
            return (_base, [])
        for gen in _gens:
            if all((x == gen._array_form[x] for x in _base)):
                for new in id_af:
                    if gen._array_form[new] != new:
                        break
                else:
                    assert None
                _base.append(new)
        strong_gens_distr = _distribute_gens_by_base(_base, _gens)
        strong_gens_slp = []
        orbs = {}
        transversals = {}
        slps = {}
        base_len = len(_base)
        for i in range(base_len):
            (transversals[i], slps[i]) = _orbit_transversal(degree, strong_gens_distr[i], _base[i], pairs=True, af=True, slp=True)
            transversals[i] = dict(transversals[i])
            orbs[i] = list(transversals[i].keys())
        i = base_len - 1
        while i >= 0:
            continue_i = False
            db = {}
            for (beta, u_beta) in list(transversals[i].items()):
                for (j, gen) in enumerate(strong_gens_distr[i]):
                    gb = gen._array_form[beta]
                    u1 = transversals[i][gb]
                    g1 = _af_rmul(gen._array_form, u_beta)
                    slp = [(i, g) for g in slps[i][beta]]
                    slp = [(i, j)] + slp
                    if g1 != u1:
                        y = True
                        try:
                            u1_inv = db[gb]
                        except KeyError:
                            u1_inv = db[gb] = _af_invert(u1)
                        schreier_gen = _af_rmul(u1_inv, g1)
                        u1_inv_slp = slps[i][gb][:]
                        u1_inv_slp.reverse()
                        u1_inv_slp = [(i, (g,)) for g in u1_inv_slp]
                        slp = u1_inv_slp + slp
                        (h, j, slp) = _strip_af(schreier_gen, _base, orbs, transversals, i, slp=slp, slps=slps)
                        if j <= base_len:
                            y = False
                        elif h:
                            y = False
                            moved = 0
                            while h[moved] == moved:
                                moved += 1
                            _base.append(moved)
                            base_len += 1
                            strong_gens_distr.append([])
                        if y is False:
                            h = _af_new(h)
                            strong_gens_slp.append((h, slp))
                            for l in range(i + 1, j):
                                strong_gens_distr[l].append(h)
                                (transversals[l], slps[l]) = _orbit_transversal(degree, strong_gens_distr[l], _base[l], pairs=True, af=True, slp=True)
                                transversals[l] = dict(transversals[l])
                                orbs[l] = list(transversals[l].keys())
                            i = j - 1
                            continue_i = True
                    if continue_i is True:
                        break
                if continue_i is True:
                    break
            if continue_i is True:
                continue
            i -= 1
        strong_gens = _gens[:]
        if slp_dict:
            for (k, slp) in strong_gens_slp:
                strong_gens.append(k)
                for i in range(len(slp)):
                    s = slp[i]
                    if isinstance(s[1], tuple):
                        slp[i] = strong_gens_distr[s[0]][s[1][0]] ** (-1)
                    else:
                        slp[i] = strong_gens_distr[s[0]][s[1]]
            strong_gens_slp = dict(strong_gens_slp)
            for g in _gens:
                strong_gens_slp[g] = [g]
            return (_base, strong_gens, strong_gens_slp)
        strong_gens.extend([k for (k, _) in strong_gens_slp])
        return (_base, strong_gens)

    def schreier_sims_random(self, base=None, gens=None, consec_succ=10, _random_prec=None):
        if base is None:
            base = []
        if gens is None:
            gens = self.generators
        base_len = len(base)
        n = self.degree
        for gen in gens:
            if all((gen(x) == x for x in base)):
                new = 0
                while gen._array_form[new] == new:
                    new += 1
                base.append(new)
                base_len += 1
        strong_gens_distr = _distribute_gens_by_base(base, gens)
        transversals = {}
        orbs = {}
        for i in range(base_len):
            transversals[i] = dict(_orbit_transversal(n, strong_gens_distr[i], base[i], pairs=True))
            orbs[i] = list(transversals[i].keys())
        c = 0
        while c < consec_succ:
            if _random_prec is None:
                g = self.random_pr()
            else:
                g = _random_prec['g'].pop()
            (h, j) = _strip(g, base, orbs, transversals)
            y = True
            if j <= base_len:
                y = False
            elif not h.is_Identity:
                y = False
                moved = 0
                while h(moved) == moved:
                    moved += 1
                base.append(moved)
                base_len += 1
                strong_gens_distr.append([])
            if y is False:
                for l in range(1, j):
                    strong_gens_distr[l].append(h)
                    transversals[l] = dict(_orbit_transversal(n, strong_gens_distr[l], base[l], pairs=True))
                    orbs[l] = list(transversals[l].keys())
                c = 0
            else:
                c += 1
        strong_gens = strong_gens_distr[0][:]
        # The reference indexes level 1 unconditionally and raises an
        # IndexError when the base has a single point.
        for gen in (strong_gens_distr[1] if len(strong_gens_distr) > 1 else []):
            if gen not in strong_gens:
                strong_gens.append(gen)
        return (base, strong_gens)

    def schreier_vector(self, alpha):
        n = self.degree
        v = [None] * n
        v[alpha] = -1
        orb = [alpha]
        used = [False] * n
        used[alpha] = True
        gens = self.generators
        r = len(gens)
        for b in orb:
            for i in range(r):
                temp = gens[i]._array_form[b]
                if used[temp] is False:
                    orb.append(temp)
                    used[temp] = True
                    v[temp] = i
        return v

    def stabilizer(self, alpha):
        return PermGroup(_stabilizer(self._degree, self._generators, alpha))

    @property
    def strong_gens(self):
        if self._strong_gens == []:
            self.schreier_sims()
        return self._strong_gens

    def subgroup(self, gens):
        if not all((g in self for g in gens)):
            raise ValueError('The group does not contain the supplied generators')
        G = PermutationGroup(gens)
        return G

    def subgroup_search(self, prop, base=None, strong_gens=None, tests=None, init_subgroup=None):

        def get_reps(orbits):
            return [min(orbit, key=lambda x: base_ordering[x]) for orbit in orbits]

        def update_nu(l):
            temp_index = len(basic_orbits[l]) + 1 - len(res_basic_orbits_init_base[l])
            if temp_index >= len(sorted_orbits[l]):
                nu[l] = base_ordering[degree]
            else:
                nu[l] = sorted_orbits[l][temp_index]
        if base is None:
            (base, strong_gens) = self.schreier_sims_incremental()
        base_len = len(base)
        degree = self.degree
        identity = _af_new(list(range(degree)))
        base_ordering = _base_ordering(base, degree)
        base_ordering.append(degree)
        base_ordering.append(-1)
        strong_gens_distr = _distribute_gens_by_base(base, strong_gens)
        (basic_orbits, transversals) = _orbits_transversals_from_bsgs(base, strong_gens_distr)
        if init_subgroup is None:
            init_subgroup = PermutationGroup([identity])
        if tests is None:
            trivial_test = lambda x: True
            tests = []
            for i in range(base_len):
                tests.append(trivial_test)
        res = init_subgroup
        f = base_len - 1
        l = base_len - 1
        res_base = base[:]
        (res_base, res_strong_gens) = res.schreier_sims_incremental(base=res_base)
        res_strong_gens_distr = _distribute_gens_by_base(res_base, res_strong_gens)
        res_generators = res.generators
        res_basic_orbits_init_base = [_orbit(degree, res_strong_gens_distr[i], res_base[i]) for i in range(base_len)]
        orbit_reps = [None] * base_len
        orbits = _orbits(degree, res_strong_gens_distr[f])
        orbit_reps[f] = get_reps(orbits)
        orbit_reps[f].remove(base[f])
        c = [0] * base_len
        u = [identity] * base_len
        sorted_orbits = [None] * base_len
        for i in range(base_len):
            sorted_orbits[i] = basic_orbits[i][:]
            sorted_orbits[i].sort(key=lambda point: base_ordering[point])
        mu = [None] * base_len
        nu = [None] * base_len
        mu[l] = degree + 1
        update_nu(l)
        computed_words = [identity] * base_len
        while True:
            while l < base_len - 1 and computed_words[l](base[l]) in orbit_reps[l] and (base_ordering[mu[l]] < base_ordering[computed_words[l](base[l])] < base_ordering[nu[l]]) and tests[l](computed_words):
                new_point = computed_words[l](base[l])
                res_base[l] = new_point
                new_stab_gens = _stabilizer(degree, res_strong_gens_distr[l], new_point)
                res_strong_gens_distr[l + 1] = new_stab_gens
                orbits = _orbits(degree, new_stab_gens)
                orbit_reps[l + 1] = get_reps(orbits)
                l += 1
                temp_orbit = [computed_words[l - 1](point) for point in basic_orbits[l]]
                temp_orbit.sort(key=lambda point: base_ordering[point])
                sorted_orbits[l] = temp_orbit
                new_mu = degree + 1
                for i in range(l):
                    if base[l] in res_basic_orbits_init_base[i]:
                        candidate = computed_words[i](base[i])
                        if base_ordering[candidate] > base_ordering[new_mu]:
                            new_mu = candidate
                mu[l] = new_mu
                update_nu(l)
                c[l] = 0
                temp_point = sorted_orbits[l][c[l]]
                gamma = computed_words[l - 1]._array_form.index(temp_point)
                u[l] = transversals[l][gamma]
                computed_words[l] = rmul(computed_words[l - 1], u[l])
            g = computed_words[l]
            temp_point = g(base[l])
            if l == base_len - 1 and base_ordering[mu[l]] < base_ordering[temp_point] < base_ordering[nu[l]] and (temp_point in orbit_reps[l]) and tests[l](computed_words) and prop(g):
                res_generators.append(g)
                res_base = base[:]
                res_strong_gens.append(g)
                res_strong_gens_distr = _distribute_gens_by_base(res_base, res_strong_gens)
                res_basic_orbits_init_base = [_orbit(degree, res_strong_gens_distr[i], res_base[i]) for i in range(base_len)]
                orbit_reps[f] = get_reps(orbits)
                l = f
            while l >= 0 and c[l] == len(basic_orbits[l]) - 1:
                l = l - 1
            if l == -1:
                return PermutationGroup(res_generators)
            if l < f:
                f = l
                c[l] = 0
                temp_orbits = _orbits(degree, res_strong_gens_distr[f])
                orbit_reps[f] = get_reps(temp_orbits)
                mu[l] = degree + 1
                temp_index = len(basic_orbits[l]) + 1 - len(res_basic_orbits_init_base[l])
                if temp_index >= len(sorted_orbits[l]):
                    nu[l] = base_ordering[degree]
                else:
                    nu[l] = sorted_orbits[l][temp_index]
            c[l] += 1
            if l == 0:
                gamma = sorted_orbits[l][c[l]]
            else:
                gamma = computed_words[l - 1]._array_form.index(sorted_orbits[l][c[l]])
            u[l] = transversals[l][gamma]
            if l == 0:
                computed_words[l] = u[l]
            else:
                computed_words[l] = rmul(computed_words[l - 1], u[l])

    @property
    def transitivity_degree(self):
        if self._transitivity_degree is None:
            n = self.degree
            G = self
            for i in range(n):
                orb = G.orbit(i)
                if len(orb) != n - i:
                    self._transitivity_degree = i
                    return i
                G = G.stabilizer(i)
            self._transitivity_degree = n
            return n
        else:
            return self._transitivity_degree

    def _p_elements_group(self, p):
        gens = self.generators[:]
        gens = sorted(gens, key=lambda x: x.order(), reverse=True)
        gens_p = [g ** (g.order() / p) for g in gens]
        gens_r = []
        for i in range(len(gens)):
            x = gens[i]
            x_order = x.order()
            x_p = x ** (x_order / p)
            if i > 0:
                P = PermutationGroup(gens_p[:i])
            else:
                P = PermutationGroup(self.identity)
            if x ** (x_order / p) not in P:
                gens_r.append(x ** (x_order / p))
            else:
                g = P.generator_product(x_p, original=True)
                for s in g:
                    x = x * s ** (-1)
                x_order = x_order / p
                del gens[i]
                del gens_p[i]
                j = i - 1
                while j < len(gens) and gens[j].order() >= x_order:
                    j += 1
                gens = gens[:j] + [x] + gens[j:]
                gens_p = gens_p[:j] + [x] + gens_p[j:]
        return PermutationGroup(gens_r)

    def _sylow_alt_sym(self, p):
        n = self.degree
        gens = []
        identity = Permutation(n - 1)
        alt = p == 2 and all((g.is_even for g in self.generators))
        coeffs = []
        m = n
        while m > 0:
            coeffs.append(m % p)
            m = m // p
        power = len(coeffs) - 1
        for i in range(1, power + 1):
            if i == 1 and alt:
                continue
            gen = Permutation([(j + p ** (i - 1)) % p ** i for j in range(p ** i)])
            gens.append(identity * gen)
            if alt:
                gen = Permutation(0, 1) * gen * Permutation(0, 1) * gen
                gens.append(gen)
        start = 0
        while power > 0:
            a = coeffs[power]
            for _ in range(a):
                shift = Permutation()
                if start > 0:
                    for i in range(p ** power):
                        shift = shift(i, start + i)
                    if alt:
                        gen = Permutation(0, 1) * shift * Permutation(0, 1) * shift
                        gens.append(gen)
                        j = 2 * (power - 1)
                    else:
                        j = power
                    for (i, gen) in enumerate(gens[:j]):
                        if alt and i % 2 == 1:
                            continue
                        gen = shift * gen * shift
                        gens.append(gen)
                start += p ** power
            power = power - 1
        return gens

    def sylow_subgroup(self, p):
        from symox.combinatorics.homomorphisms import orbit_homomorphism, block_homomorphism
        if not isprime(p):
            raise ValueError('p must be a prime')

        def is_p_group(G):
            m = G.order()
            n = 0
            while m % p == 0:
                m = m / p
                n += 1
                if m == 1:
                    return (True, n)
            return (False, n)

        def _sylow_reduce(mu, nu):
            Q = mu.image().sylow_subgroup(p)
            Q = mu.invert_subgroup(Q)
            nu = nu.restrict_to(Q)
            R = nu.image().sylow_subgroup(p)
            return nu.invert_subgroup(R)
        order = self.order()
        if order % p != 0:
            return PermutationGroup([self.identity])
        (p_group, n) = is_p_group(self)
        if p_group:
            return self
        if self.is_alt_sym():
            return PermutationGroup(self._sylow_alt_sym(p))
        orbits = self.orbits()
        non_p_orbits = [o for o in orbits if len(o) % p != 0 and len(o) != 1]
        if non_p_orbits:
            G = self.stabilizer(list(non_p_orbits[0]).pop())
            return G.sylow_subgroup(p)
        if not self.is_transitive():
            orbits = sorted(orbits, key=len)
            omega1 = orbits.pop()
            omega2 = orbits[0].union(*orbits)
            mu = orbit_homomorphism(self, omega1)
            nu = orbit_homomorphism(self, omega2)
            return _sylow_reduce(mu, nu)
        blocks = self.minimal_blocks()
        if len(blocks) > 1:
            mu = block_homomorphism(self, blocks[0])
            nu = block_homomorphism(self, blocks[1])
            return _sylow_reduce(mu, nu)
        elif len(blocks) == 1:
            block = list(blocks)[0]
            if any((e != 0 for e in block)):
                mu = block_homomorphism(self, block)
                if not is_p_group(mu.image())[0]:
                    S = mu.image().sylow_subgroup(p)
                    return mu.invert_subgroup(S).sylow_subgroup(p)
        g = self.random()
        g_order = g.order()
        while g_order % p != 0 or g_order == 0:
            g = self.random()
            g_order = g.order()
        g = g ** (g_order // p)
        if order % p ** 2 != 0:
            return PermutationGroup(g)
        C = self.centralizer(g)
        while C.order() % p ** n != 0:
            S = C.sylow_subgroup(p)
            s_order = S.order()
            Z = S.center()
            P = Z._p_elements_group(p)
            h = P.random()
            C_h = self.centralizer(h)
            while C_h.order() % p * s_order != 0:
                h = P.random()
                C_h = self.centralizer(h)
            C = C_h
        return C.sylow_subgroup(p)

    def _block_verify(self, L, alpha):
        delta = sorted(self.orbit(alpha))
        p = [-1] * len(delta)
        blocks = [-1] * len(delta)
        B = [[]]
        u = [0] * len(delta)
        t = L.orbit_transversal(alpha, pairs=True)
        for (a, beta) in t:
            B[0].append(a)
            i_a = delta.index(a)
            p[i_a] = 0
            blocks[i_a] = alpha
            u[i_a] = beta
        rho = 0
        m = 0
        while rho <= m:
            beta = B[rho][0]
            for g in self.generators:
                d = beta ^ g
                i_d = delta.index(d)
                sigma = p[i_d]
                if sigma < 0:
                    m += 1
                    sigma = m
                    u[i_d] = u[delta.index(beta)] * g
                    p[i_d] = sigma
                    rep = d
                    blocks[i_d] = rep
                    newb = [rep]
                    for gamma in B[rho][1:]:
                        i_gamma = delta.index(gamma)
                        d = gamma ^ g
                        i_d = delta.index(d)
                        if p[i_d] < 0:
                            u[i_d] = u[i_gamma] * g
                            p[i_d] = sigma
                            blocks[i_d] = rep
                            newb.append(d)
                        else:
                            s = u[i_gamma] * g * u[i_d] ** (-1)
                            return (False, s)
                    B.append(newb)
                else:
                    for h in B[rho][1:]:
                        if h ^ g not in B[sigma]:
                            s = u[delta.index(beta)] * g * u[i_d] ** (-1)
                            return (False, s)
            rho += 1
        return (True, blocks)

    def _verify(H, K, phi, z, alpha):
        orbit = H.orbit(alpha)
        beta = alpha ^ z ** (-1)
        K_beta = K.stabilizer(beta)
        gammas = [alpha, beta]
        orbits = list({tuple(K_beta.orbit(o)) for o in orbit})
        orbit_reps = [orb[0] for orb in orbits]
        for rep in orbit_reps:
            if rep not in gammas:
                gammas.append(rep)
        betas = [alpha, beta]
        transversal = {alpha: phi.invert(H.identity), beta: phi.invert(z ** (-1))}
        for (s, g) in K.orbit_transversal(beta, pairs=True):
            if s not in transversal:
                transversal[s] = transversal[beta] * phi.invert(g)
        union = K.orbit(alpha).union(K.orbit(beta))
        while len(union) < len(orbit):
            for gamma in gammas:
                if gamma in union:
                    r = gamma ^ z
                    if r not in union:
                        betas.append(r)
                        transversal[r] = transversal[gamma] * phi.invert(z)
                        for (s, g) in K.orbit_transversal(r, pairs=True):
                            if s not in transversal:
                                transversal[s] = transversal[r] * phi.invert(g)
                        union = union.union(K.orbit(r))
                        break
        rels = []
        for b in betas:
            k_gens = K.stabilizer(b).generators
            for y in k_gens:
                new_rel = transversal[b]
                gens = K.generator_product(y, original=True)
                for g in gens[::-1]:
                    new_rel = new_rel * phi.invert(g)
                new_rel = new_rel * transversal[b] ** (-1)
                perm = phi(new_rel)
                try:
                    gens = K.generator_product(perm, original=True)
                except ValueError:
                    return (False, perm)
                for g in gens:
                    new_rel = new_rel * phi.invert(g) ** (-1)
                if new_rel not in rels:
                    rels.append(new_rel)
        for gamma in gammas:
            new_rel = transversal[gamma] * phi.invert(z) * transversal[gamma ^ z] ** (-1)
            perm = phi(new_rel)
            try:
                gens = K.generator_product(perm, original=True)
            except ValueError:
                return (False, perm)
            for g in gens:
                new_rel = new_rel * phi.invert(g) ** (-1)
            if new_rel not in rels:
                rels.append(new_rel)
        return (True, rels)

    def strong_presentation(self):
        from symox.combinatorics.fp_groups import FpGroup, simplify_presentation
        from symox.combinatorics.free_groups import free_group
        from symox.combinatorics.homomorphisms import block_homomorphism, homomorphism, GroupHomomorphism
        strong_gens = self.strong_gens[:]
        stabs = self.basic_stabilizers[:]
        base = self.base[:]
        gen_syms = ['x_%d' % i for i in range(len(strong_gens))]
        F = free_group(', '.join(gen_syms))[0]
        phi = homomorphism(F, self, F.generators, strong_gens)
        H = PermutationGroup(self.identity)
        while stabs:
            alpha = base.pop()
            K = H
            H = stabs.pop()
            new_gens = [g for g in H.generators if g not in K]
            if K.order() == 1:
                z = new_gens.pop()
                rels = [F.generators[-1] ** z.order()]
                intermediate_gens = [z]
                K = PermutationGroup(intermediate_gens)
            while new_gens:
                z = new_gens.pop()
                intermediate_gens = [z] + intermediate_gens
                K_s = PermutationGroup(intermediate_gens)
                orbit = K_s.orbit(alpha)
                orbit_k = K.orbit(alpha)
                if orbit_k == orbit:
                    if z in K:
                        rel = phi.invert(z)
                        perm = z
                    else:
                        t = K.orbit_rep(alpha, alpha ^ z)
                        rel = phi.invert(z) * phi.invert(t) ** (-1)
                        perm = z * t ** (-1)
                    for g in K.generator_product(perm, original=True):
                        rel = rel * phi.invert(g) ** (-1)
                    new_rels = [rel]
                elif len(orbit_k) == 1:
                    (success, new_rels) = K_s._verify(K, phi, z, alpha)
                else:
                    (check, block) = K_s._block_verify(K, alpha)
                    if check:
                        t = block_homomorphism(K_s, block)
                        m = t.codomain.degree
                        d = K_s.degree
                        p = Permutation()
                        for i in range(m):
                            p *= Permutation(i, i + d)
                        t_img = t.images
                        images = {g: g * p * t_img[g] * p for g in t_img}
                        for g in self.strong_gens[:-len(K_s.generators)]:
                            images[g] = g
                        K_s_act = PermutationGroup(list(images.values()))
                        f = GroupHomomorphism(self, K_s_act, images)
                        K_act = PermutationGroup([f(g) for g in K.generators])
                        (success, new_rels) = K_s_act._verify(K_act, f.compose(phi), f(z), d)
                for n in new_rels:
                    if n not in rels:
                        rels.append(n)
                K = K_s
        group = FpGroup(F, rels)
        return simplify_presentation(group)

    def presentation(self, eliminate_gens=True):
        from symox.combinatorics.fp_groups import FpGroup, simplify_presentation
        from symox.combinatorics.coset_table import CosetTable
        from symox.combinatorics.free_groups import free_group
        from symox.combinatorics.homomorphisms import homomorphism
        if self._fp_presentation:
            return self._fp_presentation

        def _factor_group_by_rels(G, rels):
            if isinstance(G, FpGroup):
                rels.extend(G.relators)
                return FpGroup(G.free_group, list(set(rels)))
            return FpGroup(G, rels)
        gens = self.generators
        len_g = len(gens)
        if len_g == 1:
            order = gens[0].order()
            if order == 1:
                return free_group([])[0]
            (F, x) = free_group('x')
            return FpGroup(F, [x ** order])
        if self.order() > 20:
            half_gens = self.generators[0:(len_g + 1) // 2]
        else:
            half_gens = []
        H = PermutationGroup(half_gens)
        H_p = H.presentation()
        len_h = len(H_p.generators)
        C = self.coset_table(H)
        n = len(C)
        gen_syms = ['x_%d' % i for i in range(len(gens))]
        F = free_group(', '.join(gen_syms))[0]
        images = [F.generators[i] for i in range(len_h)]
        R = homomorphism(H_p, F, H_p.generators, images, check=False)
        rels = R(H_p.relators)
        G_p = FpGroup(F, rels)
        T = homomorphism(G_p, self, G_p.generators, gens)
        C_p = CosetTable(G_p, [])
        C_p.table = [[None] * (2 * len_g) for i in range(n)]
        transversal = [None] * n
        transversal[0] = G_p.identity
        for i in range(2 * len_h):
            C_p.table[0][i] = 0
        gamma = 1
        for (alpha, x) in product(range(n), range(2 * len_g)):
            beta = C[alpha][x]
            if beta == gamma:
                gen = G_p.generators[x // 2] ** (-1) ** (x % 2)
                transversal[beta] = transversal[alpha] * gen
                C_p.table[alpha][x] = beta
                C_p.table[beta][x + (-1) ** (x % 2)] = alpha
                gamma += 1
                if gamma == n:
                    break
        C_p.p = list(range(n))
        beta = x = 0
        while not C_p.is_complete():
            while C_p.table[beta][x] == C[beta][x]:
                x = (x + 1) % (2 * len_g)
                if x == 0:
                    beta = (beta + 1) % n
            gen = G_p.generators[x // 2] ** (-1) ** (x % 2)
            new_rel = transversal[beta] * gen * transversal[C[beta][x]] ** (-1)
            perm = T(new_rel)
            nxt = G_p.identity
            for s in H.generator_product(perm, original=True):
                nxt = nxt * T.invert(s) ** (-1)
            new_rel = new_rel * nxt
            G_p = _factor_group_by_rels(G_p, [new_rel])
            C_p.scan_and_fill(0, new_rel)
            C_p = G_p.coset_enumeration([], strategy='coset_table', draft=C_p, max_cosets=n, incomplete=True)
        self._fp_presentation = simplify_presentation(G_p)
        return self._fp_presentation

    def polycyclic_group(self):
        from symox.combinatorics.pc_groups import PolycyclicGroup
        if not self.is_polycyclic:
            raise ValueError('The group must be solvable')
        der = self.derived_series()
        pc_series = []
        pc_sequence = []
        relative_order = []
        pc_series.append(der[-1])
        der.reverse()
        for i in range(len(der) - 1):
            H = der[i]
            for g in der[i + 1].generators:
                if g not in H:
                    H = PermutationGroup([g] + H.generators)
                    pc_series.insert(0, H)
                    pc_sequence.insert(0, g)
                    G1 = pc_series[0].order()
                    G2 = pc_series[1].order()
                    relative_order.insert(0, G1 // G2)
        return PolycyclicGroup(pc_sequence, pc_series, relative_order, collector=None)

def _orbit(degree, generators, alpha, action='tuples'):
    if not hasattr(alpha, '__getitem__'):
        alpha = [alpha]
    gens = [x._array_form for x in generators]
    if len(alpha) == 1 or action == 'union':
        orb = alpha
        used = [False] * degree
        for el in alpha:
            used[el] = True
        for b in orb:
            for gen in gens:
                temp = gen[b]
                if used[temp] == False:
                    orb.append(temp)
                    used[temp] = True
        return set(orb)
    elif action == 'tuples':
        alpha = tuple(alpha)
        orb = [alpha]
        used = {alpha}
        for b in orb:
            for gen in gens:
                temp = tuple([gen[x] for x in b])
                if temp not in used:
                    orb.append(temp)
                    used.add(temp)
        return set(orb)
    elif action == 'sets':
        alpha = frozenset(alpha)
        orb = [alpha]
        used = {alpha}
        for b in orb:
            for gen in gens:
                temp = frozenset([gen[x] for x in b])
                if temp not in used:
                    orb.append(temp)
                    used.add(temp)
        return {tuple(x) for x in orb}

def _orbits(degree, generators):
    orbs = []
    sorted_I = list(range(degree))
    I = set(sorted_I)
    while I:
        i = sorted_I[0]
        orb = _orbit(degree, generators, i)
        orbs.append(orb)
        I -= orb
        sorted_I = [i for i in sorted_I if i not in orb]
    return orbs

def _orbit_transversal(degree, generators, alpha, pairs, af=False, slp=False):
    tr = [(alpha, list(range(degree)))]
    slp_dict = {alpha: []}
    used = [False] * degree
    used[alpha] = True
    gens = [x._array_form for x in generators]
    for (x, px) in tr:
        px_slp = slp_dict[x]
        for gen in gens:
            temp = gen[x]
            if used[temp] == False:
                slp_dict[temp] = [gens.index(gen)] + px_slp
                tr.append((temp, _af_rmul(gen, px)))
                used[temp] = True
    if pairs:
        if not af:
            tr = [(x, _af_new(y)) for (x, y) in tr]
        if not slp:
            return tr
        return (tr, slp_dict)
    if af:
        tr = [y for (_, y) in tr]
        if not slp:
            return tr
        return (tr, slp_dict)
    tr = [_af_new(y) for (_, y) in tr]
    if not slp:
        return tr
    return (tr, slp_dict)

def _stabilizer(degree, generators, alpha):
    orb = [alpha]
    table = {alpha: list(range(degree))}
    table_inv = {alpha: list(range(degree))}
    used = [False] * degree
    used[alpha] = True
    gens = [x._array_form for x in generators]
    stab_gens = []
    for b in orb:
        for gen in gens:
            temp = gen[b]
            if used[temp] is False:
                gen_temp = _af_rmul(gen, table[b])
                orb.append(temp)
                table[temp] = gen_temp
                table_inv[temp] = _af_invert(gen_temp)
                used[temp] = True
            else:
                schreier_gen = _af_rmuln(table_inv[temp], gen, table[b])
                if schreier_gen not in stab_gens:
                    stab_gens.append(schreier_gen)
    return [_af_new(x) for x in stab_gens]


class PermutationGroup(_ReferencePermutationGroupMethods):
    is_group = True

    def __init__(self, *generators,dups=True,**kwargs):
        if generators and isinstance(generators[0], (list, tuple)):
            generators = tuple(generators[0])
        if not generators:generators=(Permutation([]),)
        if any(isinstance(g,Cycle) for g in generators):generators=tuple(Permutation(g) for g in generators)
        if any(not isinstance(g, Permutation) for g in generators):
            raise ValueError('Expected one or more permutation generators')
        self._degree = max(g.size for g in generators)
        if len({g.size for g in generators})>1:
            degree=kwargs.pop('degree',None)
            if degree is not None:self._degree=operator.index(degree)
        if kwargs:raise TypeError('unexpected group constructor keyword: '+next(iter(kwargs)))
        generators=[Permutation(g,size=self._degree) for g in generators]
        if dups:generators=list(dict.fromkeys(generators))
        if len(generators)>1:generators=[g for g in generators if not g.is_Identity]
        self._generators=list(generators or [Permutation(range(self._degree))])
        self._elements = None
        self._chain_cache=None
        self._order=None
        self._center = None
        self._is_abelian = None
        self._is_transitive = None
        self._is_sym = None
        self._is_alt = None
        self._is_primitive = None
        self._is_nilpotent = None
        self._is_solvable = None
        self._is_trivial = None
        self._transitivity_degree = None
        self._max_div = None
        self._is_perfect = None
        self._is_cyclic = None
        self._is_dihedral = None
        self._r = len(self._generators)
        # Assigned by schreier_sims (reference base and strong generators).
        self._base = []
        self._strong_gens = []
        self._strong_gens_slp = []
        self._basic_orbits = []
        self._transversals = []
        self._transversal_slp = []
        # Assigned by _random_pr_init.
        self._random_gens = []
        self._fp_presentation = None

    @property
    def degree(self):return self._degree
    @degree.setter
    def degree(self,value):self._degree=value

    @property
    def args(self):return tuple(self._generators)

    @property
    def generators(self):
        return list(self._generators)

    @property
    def identity(self):
        return Permutation(range(self._degree))

    def _enumerate(self):
        if self._elements is None:
            self._elements=tuple(self.generate_schreier_sims())
        return self._elements

    def _chain(self):
        if self._chain_cache is None:
            self._chain_cache=StabilizerChain([g._images for g in self._generators],self._degree)
        return self._chain_cache

    @property
    def elements(self):
        return list(self._enumerate())

    def order(self):
        if self._order is None:
            if self._is_sym:self._order=factorial(self._degree)
            elif self._is_alt:self._order=factorial(self._degree)/2
            else:
                self._order=(self._generators[0].order() if len(self._generators)==1 else self._chain().order)
        return self._order

    def contains(self, g, strict=True):
        permutation=g
        if not isinstance(permutation, Permutation):
            return False
        if permutation.size != self._degree:
            if strict:
                return False
            permutation = Permutation(permutation, size=self._degree)
        if permutation.is_Identity or permutation in self._generators:return True
        if self._is_sym:return True
        if self._is_alt:return permutation.is_even
        if len(self._generators)==1:
            # Powers of one generator act by a uniform rotation on each
            # disjoint cycle. Generalized CRT checks that the shifts describe
            # one common exponent, including noncoprime cycle lengths.
            residue,modulus=0,1
            for cycle in self._generators[0].full_cyclic_form:
                positions={point:index for index,point in enumerate(cycle)}
                shift=positions.get(permutation(cycle[0]))
                if shift is None:return False
                size=len(cycle)
                if any(permutation(point)!=cycle[(index+shift)%size] for index,point in enumerate(cycle)):
                    return False
                divisor=math.gcd(modulus,size)
                if (shift-residue)%divisor:return False
                reduced=size//divisor
                multiplier=((shift-residue)//divisor*pow(modulus//divisor,-1,reduced))%reduced
                residue+=modulus*multiplier;modulus*=reduced
            return True
        return self._chain().contains(permutation._images)

    def has(self,*patterns):
        """Structural containment over the generator tree, as ``Basic.has``."""
        def walk(node):
            yield node
            if isinstance(node,PermutationGroup):
                for g in node._generators:yield from walk(g)
            elif isinstance(node,Permutation):
                yield from (_sympify(i) for i in node.array_form)
        nodes=list(walk(self))
        for pattern in patterns:
            if isinstance(pattern,type):
                if any(isinstance(node,pattern) for node in nodes):return True
            elif any(type(node)==type(pattern) and node==pattern for node in nodes):return True
        return False

    def __eq__(self, other):
        if not isinstance(other, PermutationGroup):
            return NotImplemented
        return self._generators == other._generators

    def __hash__(self):
        return hash((PermutationGroup,tuple(self._generators)))

    def __str__(self):
        from symox.printing.str import sstr
        return sstr(self)

    __repr__ = __str__

    def copy(self):
        return PermutationGroup(*self._generators)

    def __getitem__(self,index):
        return self._generators[index]

    def __len__(self):return len(self._generators)

    def _point_orbit(self,points):
        orbit={operator.index(point) for point in points}
        if any(not 0<=point<self._degree for point in orbit):raise IndexError('point outside the permutation domain')
        queue=deque(orbit)
        while queue:
            value=queue.popleft()
            for g in self._generators:
                image=g(value)
                if image not in orbit:orbit.add(image);queue.append(image)
        return orbit

    def orbit(self,alpha,action='tuples'):
        try:point=operator.index(alpha)
        except TypeError:pass
        else:return self._point_orbit((point,))
        if action=='union':return self._point_orbit(alpha)
        if action not in ('tuples','sets'):raise ValueError('unknown orbit action')
        normalize=frozenset if action=='sets' else tuple
        state=normalize(alpha);orbit={state};queue=deque([state])
        while queue:
            state=queue.popleft()
            for g in self._generators:
                image=normalize(g(point) for point in state)
                if image not in orbit:orbit.add(image);queue.append(image)
        return {tuple(state) for state in orbit} if action=='sets' else orbit

    def orbits(self,rep=False):
        unseen=set(range(self._degree));result=[]
        while unseen:
            orbit=self.orbit(min(unseen));unseen.difference_update(orbit)
            result.append(min(orbit) if rep else orbit)
        return result

    def make_perm(self,n,seed=None):
        if isinstance(n,(tuple,list,range)):
            if seed is not None:raise ValueError('a sequence of generator indices requires seed=None')
            n,seed=len(n),n
        else:
            try:n=int(n)
            except TypeError:raise ValueError('n must be an integer or sequence') from None
        if isinstance(seed,(list,tuple)):
            stream=iter(seed)
            def choose():
                for index in stream:
                    if 0<=index<len(self):return index
                raise ValueError('generator-index sequence is too short')
        elif seed is None:
            choose=lambda:randrange(len(self))
        else:
            rng=random.Random(seed);choose=lambda:rng.randrange(len(self))
        result=self.identity
        for _ in range(n):result=self._generators[choose()]*result
        return result


PermGroup=PermutationGroup

# Compatible imports reuse the existing native objects.
from symox._api_exports import export as _export
_export(globals(), 'symox.combinatorics.coset', 'Coset')
_export(globals(), 'symox.combinatorics.named_groups', 'SymmetricPermutationGroup')
