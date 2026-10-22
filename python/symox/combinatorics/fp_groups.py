"""Finitely presented groups (ported from the reference fp_groups module)."""
from symox.core.singleton import S
from symox.core.symbol import symbols
from symox.combinatorics.free_groups import FreeGroup, FreeGroupElement, free_group
from symox.combinatorics.rewritingsystem import RewritingSystem
from symox.combinatorics.coset_table import CosetTable, coset_enumeration_r, coset_enumeration_c
from symox.combinatorics import PermutationGroup
from symox.matrices import Matrix
from symox.polys.polytools import gcd
from symox.printing.defaults import DefaultPrinting
from symox.utilities import public
from symox.utilities.magic import pollute
from itertools import product

@public
def fp_group(fr_grp, relators=()):
    _fp_group = FpGroup(fr_grp, relators)
    return (_fp_group,) + tuple(_fp_group._generators)

@public
def xfp_group(fr_grp, relators=()):
    _fp_group = FpGroup(fr_grp, relators)
    return (_fp_group, _fp_group._generators)

@public
def vfp_group(fr_grpm, relators):
    _fp_group = FpGroup(symbols, relators)
    pollute([sym.name for sym in _fp_group.symbols], _fp_group.generators)
    return _fp_group

def _parse_relators(rels):
    return rels

class FpGroup(DefaultPrinting):
    is_group = True
    is_FpGroup = True
    is_PermutationGroup = False

    def __init__(self, fr_grp, relators):
        relators = _parse_relators(relators)
        self.free_group = fr_grp
        self.relators = relators
        self.generators = self._generators()
        self.dtype = type('FpGroupElement', (FpGroupElement,), {'group': self})
        self._coset_table = None
        self._is_standardized = False
        self._order = None
        self._center = None
        self._rewriting_system = RewritingSystem(self)
        self._perm_isomorphism = None
        return

    def _generators(self):
        return self.free_group.generators

    def make_confluent(self):
        self._rewriting_system.make_confluent()
        return

    def reduce(self, word):
        return self._rewriting_system.reduce(word)

    def equals(self, word1, word2):
        if self.reduce(word1 * word2 ** (-1)) == self.identity:
            return True
        elif self._rewriting_system.is_confluent:
            return False
        return None

    @property
    def identity(self):
        return self.free_group.identity

    def __contains__(self, g):
        return g in self.free_group

    def subgroup(self, gens, C=None, homomorphism=False):
        if not all((isinstance(g, FreeGroupElement) for g in gens)):
            raise ValueError('Generators must be `FreeGroupElement`s')
        if not all((g.group == self.free_group for g in gens)):
            raise ValueError('Given generators are not members of the group')
        if homomorphism:
            (g, rels, _gens) = reidemeister_presentation(self, gens, C=C, homomorphism=True)
        else:
            (g, rels) = reidemeister_presentation(self, gens, C=C)
        if g:
            g = FpGroup(g[0].group, rels)
        else:
            g = FpGroup(free_group('')[0], [])
        if homomorphism:
            from symox.combinatorics.homomorphisms import homomorphism
            return (g, homomorphism(g, self, g.generators, _gens, check=False))
        return g

    def coset_enumeration(self, H, strategy='relator_based', max_cosets=None, draft=None, incomplete=False):
        if not max_cosets:
            max_cosets = CosetTable.coset_table_max_limit
        if strategy == 'relator_based':
            C = coset_enumeration_r(self, H, max_cosets=max_cosets, draft=draft, incomplete=incomplete)
        else:
            C = coset_enumeration_c(self, H, max_cosets=max_cosets, draft=draft, incomplete=incomplete)
        if C.is_complete():
            C.compress()
        return C

    def standardize_coset_table(self):
        self._coset_table.standardize()
        self._is_standardized = True

    def coset_table(self, H, strategy='relator_based', max_cosets=None, draft=None, incomplete=False):
        if not H:
            if self._coset_table is not None:
                if not self._is_standardized:
                    self.standardize_coset_table()
            else:
                C = self.coset_enumeration([], strategy, max_cosets=max_cosets, draft=draft, incomplete=incomplete)
                self._coset_table = C
                self.standardize_coset_table()
            return self._coset_table.table
        else:
            C = self.coset_enumeration(H, strategy, max_cosets=max_cosets, draft=draft, incomplete=incomplete)
            C.standardize()
            return C.table

    def order(self, strategy='relator_based'):
        if self._order is not None:
            return self._order
        if self._coset_table is not None:
            self._order = len(self._coset_table.table)
        elif len(self.relators) == 0:
            self._order = self.free_group.order()
        elif len(self.generators) == 1:
            self._order = abs(gcd([r.array_form[0][1] for r in self.relators]))
        elif self._is_infinite():
            self._order = S.Infinity
        else:
            (gens, C) = self._finite_index_subgroup()
            if C:
                ind = len(C.table)
                self._order = ind * self.subgroup(gens, C=C).order()
            else:
                self._order = self.index([])
        return self._order

    def _is_infinite(self):
        used_gens = set()
        for r in self.relators:
            used_gens.update(r.contains_generators())
        if not set(self.generators) <= used_gens:
            return True
        abelian_rels = []
        for rel in self.relators:
            abelian_rels.append([rel.exponent_sum(g) for g in self.generators])
        if 0 in _integer_invariant_factors(abelian_rels):
            return True
        else:
            return None

    def _finite_index_subgroup(self, s=None):
        gen = self.most_frequent_generator()
        rels = list(self.generators)
        rels.extend(self.relators)
        if not s:
            if len(self.generators) == 2:
                s = [gen] + [g for g in self.generators if g != gen]
            else:
                rand = self.free_group.identity
                i = 0
                while (rand in rels or rand ** (-1) in rels or rand.is_identity) and i < 10:
                    rand = self.random()
                    i += 1
                s = [gen, rand] + [g for g in self.generators if g != gen]
        mid = (len(s) + 1) // 2
        half1 = s[:mid]
        half2 = s[mid:]
        draft1 = None
        draft2 = None
        m = 200
        C = None
        while not C and m / 2 < CosetTable.coset_table_max_limit:
            m = min(m, CosetTable.coset_table_max_limit)
            draft1 = self.coset_enumeration(half1, max_cosets=m, draft=draft1, incomplete=True)
            if draft1.is_complete():
                C = draft1
                half = half1
            else:
                draft2 = self.coset_enumeration(half2, max_cosets=m, draft=draft2, incomplete=True)
                if draft2.is_complete():
                    C = draft2
                    half = half2
            if not C:
                m *= 2
        if not C:
            return (None, None)
        C.compress()
        return (half, C)

    def most_frequent_generator(self):
        gens = self.generators
        rels = self.relators
        freqs = [sum((r.generator_count(g) for r in rels)) for g in gens]
        return gens[freqs.index(max(freqs))]

    def random(self):
        import random
        r = self.free_group.identity
        for i in range(random.randint(2, 3)):
            r = r * random.choice(self.generators) ** random.choice([1, -1])
        return r

    def index(self, H, strategy='relator_based'):
        if H == []:
            return self.order()
        else:
            C = self.coset_enumeration(H, strategy)
            return len(C.table)

    def __str__(self):
        if self.free_group.rank > 30:
            str_form = '<fp group with %s generators>' % self.free_group.rank
        else:
            str_form = '<fp group on the generators %s>' % str(self.generators)
        return str_form
    __repr__ = __str__

    def _to_perm_group(self):
        from symox.combinatorics import Permutation
        from symox.combinatorics.homomorphisms import homomorphism
        if self.order() is S.Infinity:
            raise NotImplementedError('Permutation presentation of infinite groups is not implemented')
        if self._perm_isomorphism:
            T = self._perm_isomorphism
            P = T.image()
        else:
            C = self.coset_table([])
            gens = self.generators
            images = [[C[i][2 * gens.index(g)] for i in range(len(C))] for g in gens]
            images = [Permutation(i) for i in images]
            P = PermutationGroup(images)
            T = homomorphism(self, P, gens, images, check=False)
            self._perm_isomorphism = T
        return (P, T)

    def _perm_group_list(self, method_name, *args):
        (P, T) = self._to_perm_group()
        perm_result = getattr(P, method_name)(*args)
        single = False
        if isinstance(perm_result, PermutationGroup):
            (perm_result, single) = ([perm_result], True)
        result = []
        for group in perm_result:
            gens = group.generators
            result.append(T.invert(gens))
        return result[0] if single else result

    def derived_series(self):
        return self._perm_group_list('derived_series')

    def lower_central_series(self):
        return self._perm_group_list('lower_central_series')

    def center(self):
        return self._perm_group_list('center')

    def derived_subgroup(self):
        return self._perm_group_list('derived_subgroup')

    def centralizer(self, other):
        T = self._to_perm_group()[1]
        other = T(other)
        return self._perm_group_list('centralizer', other)

    def normal_closure(self, other):
        T = self._to_perm_group()[1]
        other = T(other)
        return self._perm_group_list('normal_closure', other)

    def _perm_property(self, attr):
        P = self._to_perm_group()[0]
        return getattr(P, attr)

    @property
    def is_abelian(self):
        return self._perm_property('is_abelian')

    @property
    def is_nilpotent(self):
        return self._perm_property('is_nilpotent')

    @property
    def is_solvable(self):
        return self._perm_property('is_solvable')

    @property
    def elements(self):
        (P, T) = self._to_perm_group()
        return T.invert(P.elements)

    @property
    def is_cyclic(self):
        if len(self.generators) <= 1:
            return True
        try:
            (P, T) = self._to_perm_group()
        except NotImplementedError:
            raise NotImplementedError('Check for infinite Cyclic group is not implemented')
        return P.is_cyclic

    def abelian_invariants(self):
        try:
            (P, T) = self._to_perm_group()
        except NotImplementedError:
            raise NotImplementedError('abelian invariants is not implementedfor infinite group')
        return P.abelian_invariants()

    def composition_series(self):
        try:
            (P, T) = self._to_perm_group()
        except NotImplementedError:
            raise NotImplementedError('composition series is not implementedfor infinite group')
        return P.composition_series()

class FpSubgroup(DefaultPrinting):

    def __init__(self, G, gens, normal=False):
        super().__init__()
        self.parent = G
        self.generators = list({g for g in gens if g != G.identity})
        self._min_words = None
        self.C = None
        self.normal = normal

    def __contains__(self, g):
        if isinstance(self.parent, FreeGroup):
            if self._min_words is None:

                def _process(w):
                    (p, r) = w.cyclic_reduction(removed=True)
                    if not r.is_identity:
                        return [(r, p)]
                    else:
                        return [w, w ** (-1)]
                gens = []
                for w in self.generators:
                    if self.normal:
                        w = w.cyclic_reduction()
                    gens.extend(_process(w))
                for w1 in gens:
                    for w2 in gens:
                        if w1 == w2 or (not isinstance(w1, tuple) and w1 ** (-1) == w2):
                            continue
                        if isinstance(w1, tuple):
                            (s1, s2) = (w1[0][0], w1[0][0] ** (-1))
                        else:
                            (s1, s2) = (w1[0], w1[len(w1) - 1])
                        if isinstance(w2, tuple):
                            (r1, r2) = (w2[0][0], w2[0][0] ** (-1))
                        else:
                            (r1, r2) = (w2[0], w2[len(w1) - 1])
                        (p1, p2) = (w1, w2)
                        if isinstance(w1, tuple):
                            p1 = w1[0] * w1[1] * w1[0] ** (-1)
                        if isinstance(w2, tuple):
                            p2 = w2[0] * w2[1] * w2[0] ** (-1)
                        if r1 ** (-1) == s2 and (not (p1 * p2).is_identity):
                            new = _process(p1 * p2)
                            if new not in gens:
                                gens.extend(new)
                        if r2 ** (-1) == s1 and (not (p2 * p1).is_identity):
                            new = _process(p2 * p1)
                            if new not in gens:
                                gens.extend(new)
                self._min_words = gens
            min_words = self._min_words

            def _is_subword(w):
                (w, r) = w.cyclic_reduction(removed=True)
                if r.is_identity or self.normal:
                    return w in min_words
                else:
                    t = [s[1] for s in min_words if isinstance(s, tuple) and s[0] == r]
                    return [s for s in t if w.power_of(s)] != []
            known = {}

            def _word_break(w):
                if len(w) == 0:
                    return True
                i = 0
                while i < len(w):
                    i += 1
                    prefix = w.subword(0, i)
                    if not _is_subword(prefix):
                        continue
                    rest = w.subword(i, len(w))
                    if rest not in known:
                        known[rest] = _word_break(rest)
                    if known[rest]:
                        return True
                return False
            if self.normal:
                g = g.cyclic_reduction()
            return _word_break(g)
        else:
            if self.C is None:
                C = self.parent.coset_enumeration(self.generators)
                self.C = C
            i = 0
            C = self.C
            for j in range(len(g)):
                i = C.table[i][C.A_dict[g[j]]]
            return i == 0

    def order(self):
        if not self.generators:
            return S.One
        if isinstance(self.parent, FreeGroup):
            return S.Infinity
        if self.C is None:
            C = self.parent.coset_enumeration(self.generators)
            self.C = C
        return self.parent.order() / len(self.C.table)

    def to_FpGroup(self):
        if isinstance(self.parent, FreeGroup):
            gen_syms = ['x_%d' % i for i in range(len(self.generators))]
            return free_group(', '.join(gen_syms))[0]
        return self.parent.subgroup(C=self.C)

    def __str__(self):
        if len(self.generators) > 30:
            str_form = '<fp subgroup with %s generators>' % len(self.generators)
        else:
            str_form = '<fp subgroup on the generators %s>' % str(self.generators)
        return str_form
    __repr__ = __str__

def low_index_subgroups(G, N, Y=()):
    C = CosetTable(G, [])
    R = G.relators
    len_short_rel = 5
    R2 = {rel for rel in R if len(rel) > len_short_rel}
    R1 = {rel.identity_cyclic_reduction() for rel in set(R) - R2}
    R1_c_list = C.conjugates(R1)
    S = []
    descendant_subgroups(S, C, R1_c_list, C.A[0], R2, N, Y)
    return S

def descendant_subgroups(S, C, R1_c_list, x, R2, N, Y):
    A_dict = C.A_dict
    A_dict_inv = C.A_dict_inv
    if C.is_complete():
        for (w, alpha) in product(R2, C.omega):
            if not C.scan_check(alpha, w):
                return
        S.append(C)
    else:
        for (alpha, x) in product(range(len(C.table)), C.A):
            if C.table[alpha][A_dict[x]] is None:
                (undefined_coset, undefined_gen) = (alpha, x)
                break
        reach = C.omega + [C.n]
        for beta in reach:
            if beta < N:
                if beta == C.n or C.table[beta][A_dict_inv[undefined_gen]] is None:
                    try_descendant(S, C, R1_c_list, R2, N, undefined_coset, undefined_gen, beta, Y)

def try_descendant(S, C, R1_c_list, R2, N, alpha, x, beta, Y):
    D = C.copy()
    if beta == D.n and beta < N:
        D.table.append([None] * len(D.A))
        D.p.append(beta)
    D.table[alpha][D.A_dict[x]] = beta
    D.table[beta][D.A_dict_inv[x]] = alpha
    D.deduction_stack.append((alpha, x))
    if not D.process_deductions_check(R1_c_list[D.A_dict[x]], R1_c_list[D.A_dict_inv[x]]):
        return
    for w in Y:
        if not D.scan_check(0, w):
            return
    if first_in_class(D, Y):
        descendant_subgroups(S, D, R1_c_list, x, R2, N, Y)

def first_in_class(C, Y=()):
    n = C.n
    lamda = -1
    nu = [None] * n
    mu = [None] * n
    next_alpha = False
    for alpha in range(1, n):
        for beta in range(lamda + 1):
            nu[mu[beta]] = None
        for w in Y:
            if C.table[alpha][C.A_dict[w]] != alpha:
                next_alpha = True
                break
        if next_alpha:
            next_alpha = False
            continue
        mu[0] = alpha
        nu[alpha] = 0
        lamda = 0
        for beta in range(n):
            for x in C.A:
                gamma = C.table[beta][C.A_dict[x]]
                delta = C.table[mu[beta]][C.A_dict[x]]
                if gamma is None or delta is None:
                    next_alpha = True
                    break
                if nu[delta] is None:
                    lamda += 1
                    nu[delta] = lamda
                    mu[lamda] = delta
                if nu[delta] < gamma:
                    return False
                if nu[delta] > gamma:
                    next_alpha = True
                    break
            if next_alpha:
                next_alpha = False
                break
    return True

def simplify_presentation(*args, change_gens=False):
    if len(args) == 1:
        if not isinstance(args[0], FpGroup):
            raise TypeError('The argument must be an instance of FpGroup')
        G = args[0]
        (gens, rels) = simplify_presentation(G.generators, G.relators, change_gens=change_gens)
        if gens:
            return FpGroup(gens[0].group, rels)
        return FpGroup(FreeGroup([]), [])
    elif len(args) == 2:
        (gens, rels) = (args[0][:], args[1][:])
        if not gens:
            return (gens, rels)
        identity = gens[0].group.identity
    else:
        if len(args) == 0:
            m = 'Not enough arguments'
        else:
            m = 'Too many arguments'
        raise RuntimeError(m)
    prev_gens = []
    prev_rels = []
    while not set(prev_rels) == set(rels):
        prev_rels = rels
        while change_gens and (not set(prev_gens) == set(gens)):
            prev_gens = gens
            (gens, rels) = elimination_technique_1(gens, rels, identity)
        rels = _simplify_relators(rels)
    if change_gens:
        syms = [g.array_form[0][0] for g in gens]
        F = free_group(syms)[0]
        identity = F.identity
        gens = F.generators
        subs = dict(zip(syms, gens))
        for (j, r) in enumerate(rels):
            a = r.array_form
            rel = identity
            for (sym, p) in a:
                rel = rel * subs[sym] ** p
            rels[j] = rel
    return (gens, rels)

def _simplify_relators(rels):
    rels = rels[:]
    if not rels:
        return []
    identity = rels[0].group.identity
    exps = {}
    for i in range(len(rels)):
        rel = rels[i]
        if rel.number_syllables() == 1:
            g = rel[0]
            exp = abs(rel.array_form[0][1])
            if rel.array_form[0][1] < 0:
                rels[i] = rels[i] ** (-1)
                g = g ** (-1)
            if g in exps:
                exp = gcd(exp, exps[g].array_form[0][1])
            exps[g] = g ** exp
    one_syllables_words = list(exps.values())
    for (i, rel) in enumerate(rels):
        if rel in one_syllables_words:
            continue
        rel = rel.eliminate_words(one_syllables_words, _all=True)
        for g in rel.contains_generators():
            if g in exps:
                exp = exps[g].array_form[0][1]
                max_exp = (exp + 1) // 2
                rel = rel.eliminate_word(g ** max_exp, g ** (max_exp - exp), _all=True)
                rel = rel.eliminate_word(g ** (-max_exp), g ** (-(max_exp - exp)), _all=True)
        rels[i] = rel
    rels = [r.identity_cyclic_reduction() for r in rels]
    rels += one_syllables_words
    rels = list(set(rels))
    rels.sort()
    try:
        rels.remove(identity)
    except ValueError:
        pass
    return rels

def elimination_technique_1(gens, rels, identity):
    rels = rels[:]
    rels.sort()
    gens = gens[:]
    redundant_gens = {}
    redundant_rels = []
    used_gens = set()
    for rel in rels:
        contained_gens = rel.contains_generators()
        if any((g in contained_gens for g in redundant_gens)):
            continue
        contained_gens = list(contained_gens)
        contained_gens.sort(reverse=True)
        for gen in contained_gens:
            if rel.generator_count(gen) == 1 and gen not in used_gens:
                k = rel.exponent_sum(gen)
                gen_index = rel.index(gen ** k)
                bk = rel.subword(gen_index + 1, len(rel))
                fw = rel.subword(0, gen_index)
                chi = bk * fw
                redundant_gens[gen] = chi ** (-1 * k)
                used_gens.update(chi.contains_generators())
                redundant_rels.append(rel)
                break
    rels = [r for r in rels if r not in redundant_rels]
    rels = [r.eliminate_words(redundant_gens, _all=True).identity_cyclic_reduction() for r in rels]
    rels = list(set(rels))
    try:
        rels.remove(identity)
    except ValueError:
        pass
    gens = [g for g in gens if g not in redundant_gens]
    return (gens, rels)

def define_schreier_generators(C, homomorphism=False):
    y = []
    gamma = 1
    f = C.fp_group
    X = f.generators
    if homomorphism:
        _gens = {}
        tau = {}
        tau[0] = f.identity
    C.P = [[None] * len(C.A) for i in range(C.n)]
    for (alpha, x) in product(C.omega, C.A):
        beta = C.table[alpha][C.A_dict[x]]
        if beta == gamma:
            C.P[alpha][C.A_dict[x]] = '<identity>'
            C.P[beta][C.A_dict_inv[x]] = '<identity>'
            gamma += 1
            if homomorphism:
                tau[beta] = tau[alpha] * x
        elif x in X and C.P[alpha][C.A_dict[x]] is None:
            y_alpha_x = '%s_%s' % (x, alpha)
            y.append(y_alpha_x)
            C.P[alpha][C.A_dict[x]] = y_alpha_x
            if homomorphism:
                _gens[y_alpha_x] = tau[alpha] * x * tau[beta] ** (-1)
    grp_gens = list(free_group(', '.join(y)))
    C._schreier_free_group = grp_gens.pop(0)
    C._schreier_generators = grp_gens
    if homomorphism:
        C._schreier_gen_elem = _gens
    for (i, j) in product(range(len(C.P)), range(len(C.A))):
        if C.P[i][j] == '<identity>':
            C.P[i][j] = C._schreier_free_group.identity
        elif isinstance(C.P[i][j], str):
            r = C._schreier_generators[y.index(C.P[i][j])]
            C.P[i][j] = r
            beta = C.table[i][j]
            C.P[beta][j + 1] = r ** (-1)

def reidemeister_relators(C):
    R = C.fp_group.relators
    rels = [rewrite(C, coset, word) for word in R for coset in range(C.n)]
    order_1_gens = {i for i in rels if len(i) == 1}
    rels = list(filter(lambda rel: rel not in order_1_gens, rels))
    for i in range(len(rels)):
        w = rels[i]
        w = w.eliminate_words(order_1_gens, _all=True)
        rels[i] = w
    C._schreier_generators = [i for i in C._schreier_generators if not (i in order_1_gens or i ** (-1) in order_1_gens)]
    i = 0
    while i < len(rels):
        w = rels[i]
        j = i + 1
        while j < len(rels):
            if w.is_cyclic_conjugate(rels[j]):
                del rels[j]
            else:
                j += 1
        i += 1
    C._reidemeister_relators = rels

def rewrite(C, alpha, w):
    v = C._schreier_free_group.identity
    for i in range(len(w)):
        x_i = w[i]
        v = v * C.P[alpha][C.A_dict[x_i]]
        alpha = C.table[alpha][C.A_dict[x_i]]
    return v

def elimination_technique_2(C):
    rels = C._reidemeister_relators
    rels.sort(reverse=True)
    gens = C._schreier_generators
    for i in range(len(gens) - 1, -1, -1):
        rel = rels[i]
        for j in range(len(gens) - 1, -1, -1):
            gen = gens[j]
            if rel.generator_count(gen) == 1:
                k = rel.exponent_sum(gen)
                gen_index = rel.index(gen ** k)
                bk = rel.subword(gen_index + 1, len(rel))
                fw = rel.subword(0, gen_index)
                rep_by = (bk * fw) ** (-1 * k)
                del rels[i]
                del gens[j]
                rels = [rel.eliminate_word(gen, rep_by) for rel in rels]
                break
    C._reidemeister_relators = rels
    C._schreier_generators = gens
    return (C._schreier_generators, C._reidemeister_relators)

def reidemeister_presentation(fp_grp, H, C=None, homomorphism=False):
    if not C:
        C = coset_enumeration_r(fp_grp, H)
    C.compress()
    C.standardize()
    define_schreier_generators(C, homomorphism=homomorphism)
    reidemeister_relators(C)
    (gens, rels) = (C._schreier_generators, C._reidemeister_relators)
    (gens, rels) = simplify_presentation(gens, rels, change_gens=True)
    C.schreier_generators = tuple(gens)
    C.reidemeister_relators = tuple(rels)
    if homomorphism:
        _gens = [C._schreier_gen_elem[str(gen)] for gen in gens]
        return (C.schreier_generators, C.reidemeister_relators, _gens)
    return (C.schreier_generators, C.reidemeister_relators)
FpGroupElement = FreeGroupElement


def _integer_invariant_factors(rows):
    """Diagonal of the integer Smith normal form, as ``invariant_factors``."""
    A = [[int(x) for x in row] for row in rows]
    m = len(A)
    n = len(A[0]) if m else 0
    k = min(m, n)
    diag = []
    t = 0
    while t < k:
        entries = [(abs(A[i][j]), i, j) for i in range(t, m) for j in range(t, n) if A[i][j]]
        if not entries:
            break
        _, i, j = min(entries)
        A[t], A[i] = A[i], A[t]
        for row in A:
            row[t], row[j] = row[j], row[t]
        while True:
            p = A[t][t]
            changed = False
            for i in range(t + 1, m):
                q = A[i][t] // p
                if q:
                    A[i] = [a - q*b for a, b in zip(A[i], A[t])]
                if A[i][t]:
                    changed = True
            for j in range(t + 1, n):
                q = A[t][j] // p
                if q:
                    for row in A:
                        row[j] -= q*row[t]
                if A[t][j]:
                    changed = True
            if not changed:
                bad = next(((i, j) for i in range(t + 1, m) for j in range(t + 1, n) if A[i][j] % p), None)
                if bad is None:
                    break
                A[t] = [a + b for a, b in zip(A[t], A[bad[0]])]
                continue
            entries = [(abs(A[i][t]), i, t) for i in range(t, m) if A[i][t]] + [(abs(A[t][j]), t, j) for j in range(t, n) if A[t][j]]
            _, i, j = min(entries)
            A[t], A[i] = A[i], A[t]
            for row in A:
                row[t], row[j] = row[j], row[t]
        diag.append(abs(A[t][t]))
        t += 1
    return tuple(diag + [0]*(k - len(diag)))
