"""Polycyclic groups and collectors (ported from the reference pc_groups module)."""
from symox.ntheory.primetest import isprime
from symox.combinatorics.perm_groups import PermutationGroup
from symox.printing.defaults import DefaultPrinting
from symox.combinatorics.free_groups import free_group

class PolycyclicGroup(DefaultPrinting):
    is_group = True
    is_solvable = True

    def __init__(self, pc_sequence, pc_series, relative_order, collector=None):
        self.pcgs = pc_sequence
        self.pc_series = pc_series
        self.relative_order = relative_order
        self.collector = Collector(self.pcgs, pc_series, relative_order) if not collector else collector

    def is_prime_order(self):
        return all((isprime(order) for order in self.relative_order))

    def length(self):
        return len(self.pcgs)

class Collector(DefaultPrinting):

    def __init__(self, pcgs, pc_series, relative_order, free_group_=None, pc_presentation=None):
        self.pcgs = pcgs
        self.pc_series = pc_series
        self.relative_order = relative_order
        self.free_group = free_group('x:{}'.format(len(pcgs)))[0] if not free_group_ else free_group_
        self.index = {s: i for (i, s) in enumerate(self.free_group.symbols)}
        if not pc_series:
            # Extension: the reference indexes an empty series and fails;
            # the tail subgroups of the sequence form its polycyclic series.
            identity = [pcgs[0] ** 0] if pcgs else []
            self.pc_series = [PermutationGroup(list(pcgs[i:]) or identity)
                              for i in range(len(pcgs) + 1)] if pcgs else []
            tails = self.pc_series
            for (i, g) in enumerate(pcgs):
                if (tails[i].order() != relative_order[i] * tails[i + 1].order()
                        or not tails[i + 1].is_normal(tails[i])):
                    raise ValueError('the generators are not a polycyclic sequence')
        self.pc_presentation = self.pc_relators()

    def _permutation(self, word):
        """Evaluate a word (or the empty relator ``()``) on the pcgs."""
        value = self.pcgs[0] ** 0
        if isinstance(word, tuple) and not word:
            return value
        if word.group != self.free_group:
            raise ValueError('word belongs to another free group')
        for (s, e) in word.array_form:
            value = value * self.pcgs[self.index[s]] ** e
        return value

    def minimal_uncollected_subword(self, word):
        if not word:
            return None
        array = word.array_form
        re = self.relative_order
        index = self.index
        for i in range(len(array)):
            (s1, e1) = array[i]
            if re[index[s1]] and (e1 < 0 or e1 > re[index[s1]] - 1):
                return ((s1, e1),)
        for i in range(len(array) - 1):
            (s1, e1) = array[i]
            (s2, e2) = array[i + 1]
            if index[s1] > index[s2]:
                e = 1 if e2 > 0 else -1
                return ((s1, e1), (s2, e))
        return None

    def relations(self):
        power_relators = {}
        conjugate_relators = {}
        for (key, value) in self.pc_presentation.items():
            if len(key.array_form) == 1:
                power_relators[key] = value
            else:
                conjugate_relators[key] = value
        return (power_relators, conjugate_relators)

    def subword_index(self, word, w):
        low = -1
        high = -1
        for i in range(len(word) - len(w) + 1):
            if word.subword(i, i + len(w)) == w:
                low = i
                high = i + len(w)
                break
        return (low, high)

    def map_relation(self, w):
        array = w.array_form
        s1 = array[0][0]
        s2 = array[1][0]
        key = ((s2, -1), (s1, 1), (s2, 1))
        key = self.free_group.dtype(key)
        return self.pc_presentation[key]

    def collected_word(self, word):
        free_group = self.free_group
        while True:
            w = self.minimal_uncollected_subword(word)
            if not w:
                break
            (low, high) = self.subword_index(word, free_group.dtype(w))
            if low == -1:
                continue
            (s1, e1) = w[0]
            if len(w) == 1:
                re = self.relative_order[self.index[s1]]
                q = e1 // re
                r = e1 - q * re
                key = ((w[0][0], re),)
                key = free_group.dtype(key)
                if self.pc_presentation[key]:
                    # The reference keeps only the first syllable of the
                    # power relator; x**e = x**r * (x**re)**q needs all of it.
                    word_ = free_group.dtype(((w[0][0], 1),)) ** r * self.pc_presentation[key] ** q
                elif r != 0:
                    word_ = ((w[0][0], r),)
                    word_ = free_group.dtype(word_)
                else:
                    word_ = None
                word = word.eliminate_word(free_group.dtype(w), word_)
            if len(w) == 2 and w[1][1] > 0:
                (s2, e2) = w[1]
                s2 = ((s2, 1),)
                s2 = free_group.dtype(s2)
                word_ = self.map_relation(free_group.dtype(w))
                word_ = s2 * word_ ** e1
                word_ = free_group.dtype(word_)
                word = word.substituted_word(low, high, word_)
            elif len(w) == 2 and w[1][1] < 0:
                (s2, e2) = w[1]
                s2 = ((s2, 1),)
                s2 = free_group.dtype(s2)
                word_ = self.map_relation(free_group.dtype(w))
                word_ = s2 ** (-1) * word_ ** e1
                word_ = free_group.dtype(word_)
                word = word.substituted_word(low, high, word_)
        return word

    def pc_relators(self):
        free_group = self.free_group
        rel_order = self.relative_order
        pc_relators = {}
        perm_to_free = {}
        pcgs = self.pcgs
        for (gen, s) in zip(pcgs, free_group.generators):
            perm_to_free[gen ** (-1)] = s ** (-1)
            perm_to_free[gen] = s
        pcgs = pcgs[::-1]
        series = self.pc_series[::-1]
        rel_order = rel_order[::-1]
        collected_gens = []
        for (i, gen) in enumerate(pcgs):
            re = rel_order[i]
            relation = perm_to_free[gen] ** re
            G = series[i]
            l = G.generator_product(gen ** re, original=True)
            l.reverse()
            word = free_group.identity
            for g in l:
                word = word * perm_to_free[g]
            word = self.collected_word(word)
            pc_relators[relation] = word if word else ()
            self.pc_presentation = pc_relators
            collected_gens.append(gen)
            if len(collected_gens) > 1:
                conj = collected_gens[len(collected_gens) - 1]
                conjugator = perm_to_free[conj]
                for j in range(len(collected_gens) - 1):
                    conjugated = perm_to_free[collected_gens[j]]
                    relation = conjugator ** (-1) * conjugated * conjugator
                    gens = conj ** (-1) * collected_gens[j] * conj
                    l = G.generator_product(gens, original=True)
                    l.reverse()
                    word = free_group.identity
                    for g in l:
                        word = word * perm_to_free[g]
                    word = self.collected_word(word)
                    pc_relators[relation] = word if word else ()
                    self.pc_presentation = pc_relators
        return pc_relators

    def exponent_vector(self, element):
        free_group = self.free_group
        G = PermutationGroup()
        for g in self.pcgs:
            G = PermutationGroup([g] + G.generators)
        gens = G.generator_product(element, original=True)
        gens.reverse()
        perm_to_free = {}
        for (sym, g) in zip(free_group.generators, self.pcgs):
            perm_to_free[g ** (-1)] = sym ** (-1)
            perm_to_free[g] = sym
        w = free_group.identity
        for g in gens:
            w = w * perm_to_free[g]
        word = self.collected_word(w)
        index = self.index
        exp_vector = [0] * len(free_group)
        word = word.array_form
        for t in word:
            exp_vector[index[t[0]]] = t[1]
        return exp_vector

    def depth(self, element):
        exp_vector = self.exponent_vector(element)
        return next((i + 1 for (i, x) in enumerate(exp_vector) if x), len(self.pcgs) + 1)

    def leading_exponent(self, element):
        exp_vector = self.exponent_vector(element)
        depth = self.depth(element)
        if depth != len(self.pcgs) + 1:
            return exp_vector[depth - 1]
        return None

    def _inverse_exponent(self, value, depth):
        # The reference multiplies by the float ``value**-1``; the leading
        # exponent has to be inverted modulo the relative order instead.
        order = self.relative_order[depth - 1]
        try:
            return pow(int(value), -1, int(order))
        except ValueError:
            raise NotImplementedError('leading exponent %s is not invertible '
                'modulo the relative order %s' % (value, order)) from None

    def _sift(self, z, g):
        h = g
        d = self.depth(h)
        # The reference stops at depth len(pcgs) and never sifts elements of
        # the last depth.
        while d <= len(self.pcgs) and z[d - 1] != 1:
            k = z[d - 1]
            e = self.leading_exponent(h) * self._inverse_exponent(self.leading_exponent(k), d)
            e = e % self.relative_order[d - 1]
            h = k ** (-e) * h
            d = self.depth(h)
        return h

    def induced_pcgs(self, gens):
        z = [1] * len(self.pcgs)
        # The reference pops from (and so empties) the caller's list.
        G = list(gens)
        while G:
            g = G.pop(0)
            h = self._sift(z, g)
            d = self.depth(h)
            if d <= len(self.pcgs):
                for gen in z:
                    if gen != 1:
                        G.append(h ** (-1) * gen ** (-1) * h * gen)
                # The reference omits the power of the new element, so it
                # misses (0 2)(1 3) in the subgroup generated by (0 1 2 3).
                G.append(h ** self.relative_order[d - 1])
                z[d - 1] = h
        z = [gen for gen in z if gen != 1]
        return z

    def constructive_membership_test(self, ipcgs, g):
        e = [0] * len(ipcgs)
        h = g
        d = self.depth(h)
        for (i, gen) in enumerate(ipcgs):
            while self.depth(gen) == d:
                f = self.leading_exponent(h) * self._inverse_exponent(self.leading_exponent(gen), d)
                f = f % self.relative_order[d - 1]
                h = gen ** (-f) * h
                e[i] = f
                d = self.depth(h)
        # The reference compares a Permutation with 1, which is never equal,
        # so it reports every member as outside the subgroup.
        if h.is_Identity:
            return e
        return False
