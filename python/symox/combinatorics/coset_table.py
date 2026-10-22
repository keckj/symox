"""Coset tables and Todd-Coxeter coset enumeration (ported from the reference module)."""
from symox.combinatorics.free_groups import free_group
from symox.printing.defaults import DefaultPrinting
from itertools import chain, product
from bisect import bisect_left

class CosetTable(DefaultPrinting):
    coset_table_max_limit = 4096000
    coset_table_limit = None
    max_stack_size = 100

    def __init__(self, fp_grp, subgroup, max_cosets=None):
        if not max_cosets:
            max_cosets = CosetTable.coset_table_max_limit
        self.fp_group = fp_grp
        self.subgroup = subgroup
        self.coset_table_limit = max_cosets
        self.p = [0]
        self.A = list(chain.from_iterable(((gen, gen ** (-1)) for gen in self.fp_group.generators)))
        self.P = [[None] * len(self.A)]
        self.table = [[None] * len(self.A)]
        self.A_dict = {x: self.A.index(x) for x in self.A}
        self.A_dict_inv = {}
        for (x, index) in self.A_dict.items():
            if index % 2 == 0:
                self.A_dict_inv[x] = self.A_dict[x] + 1
            else:
                self.A_dict_inv[x] = self.A_dict[x] - 1
        self.deduction_stack = []
        H = self.subgroup
        self._grp = free_group(', '.join(['a_%d' % i for i in range(len(H))]))[0]
        self.P = [[None] * len(self.A)]
        self.p_p = {}

    @property
    def omega(self):
        return [coset for coset in range(len(self.p)) if self.p[coset] == coset]

    def copy(self):
        self_copy = self.__class__(self.fp_group, self.subgroup)
        self_copy.table = [list(perm_rep) for perm_rep in self.table]
        self_copy.p = list(self.p)
        self_copy.deduction_stack = list(self.deduction_stack)
        return self_copy

    def __str__(self):
        return 'Coset Table on %s with %s as subgroup generators' % (self.fp_group, self.subgroup)
    __repr__ = __str__

    @property
    def n(self):
        if not self.table:
            return 0
        return max(self.omega) + 1

    def is_complete(self):
        return not any((None in self.table[coset] for coset in self.omega))

    def define(self, alpha, x, modified=False):
        A = self.A
        table = self.table
        len_table = len(table)
        if len_table >= self.coset_table_limit:
            raise ValueError('the coset enumeration has defined more than %s cosets. Try with a greater value max number of cosets ' % self.coset_table_limit)
        table.append([None] * len(A))
        self.P.append([None] * len(self.A))
        beta = len_table
        self.p.append(beta)
        table[alpha][self.A_dict[x]] = beta
        table[beta][self.A_dict_inv[x]] = alpha
        if modified:
            self.P[alpha][self.A_dict[x]] = self._grp.identity
            self.P[beta][self.A_dict_inv[x]] = self._grp.identity
            self.p_p[beta] = self._grp.identity

    def define_c(self, alpha, x):
        A = self.A
        table = self.table
        len_table = len(table)
        if len_table >= self.coset_table_limit:
            raise ValueError('the coset enumeration has defined more than %s cosets. Try with a greater value max number of cosets ' % self.coset_table_limit)
        table.append([None] * len(A))
        beta = len_table
        self.p.append(beta)
        table[alpha][self.A_dict[x]] = beta
        table[beta][self.A_dict_inv[x]] = alpha
        self.deduction_stack.append((alpha, x))

    def scan_c(self, alpha, word):
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        while i <= j and table[f][A_dict[word[i]]] is not None:
            f = table[f][A_dict[word[i]]]
            i += 1
        if i > j:
            if f != b:
                self.coincidence_c(f, b)
            return
        while j >= i and table[b][A_dict_inv[word[j]]] is not None:
            b = table[b][A_dict_inv[word[j]]]
            j -= 1
        if j < i:
            self.coincidence_c(f, b)
        elif j == i:
            table[f][A_dict[word[i]]] = b
            table[b][A_dict_inv[word[i]]] = f
            self.deduction_stack.append((f, word[i]))

    def coincidence_c(self, alpha, beta):
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        q = []
        self.merge(alpha, beta, q)
        while len(q) > 0:
            gamma = q.pop(0)
            for x in A_dict:
                delta = table[gamma][A_dict[x]]
                if delta is not None:
                    table[delta][A_dict_inv[x]] = None
                    self.deduction_stack.append((delta, x ** (-1)))
                    mu = self.rep(gamma)
                    nu = self.rep(delta)
                    if table[mu][A_dict[x]] is not None:
                        self.merge(nu, table[mu][A_dict[x]], q)
                    elif table[nu][A_dict_inv[x]] is not None:
                        self.merge(mu, table[nu][A_dict_inv[x]], q)
                    else:
                        table[mu][A_dict[x]] = nu
                        table[nu][A_dict_inv[x]] = mu

    def scan(self, alpha, word, y=None, fill=False, modified=False):
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        b_p = y
        if modified:
            f_p = self._grp.identity
        flag = 0
        while fill or flag == 0:
            flag = 1
            while i <= j and table[f][A_dict[word[i]]] is not None:
                if modified:
                    f_p = f_p * self.P[f][A_dict[word[i]]]
                f = table[f][A_dict[word[i]]]
                i += 1
            if i > j:
                if f != b:
                    if modified:
                        self.modified_coincidence(f, b, f_p ** (-1) * y)
                    else:
                        self.coincidence(f, b)
                return
            while j >= i and table[b][A_dict_inv[word[j]]] is not None:
                if modified:
                    b_p = b_p * self.P[b][self.A_dict_inv[word[j]]]
                b = table[b][A_dict_inv[word[j]]]
                j -= 1
            if j < i:
                if modified:
                    self.modified_coincidence(f, b, f_p ** (-1) * b_p)
                else:
                    self.coincidence(f, b)
            elif j == i:
                table[f][A_dict[word[i]]] = b
                table[b][A_dict_inv[word[i]]] = f
                if modified:
                    self.P[f][self.A_dict[word[i]]] = f_p ** (-1) * b_p
                    self.P[b][self.A_dict_inv[word[i]]] = b_p ** (-1) * f_p
                return
            elif fill:
                self.define(f, word[i], modified=modified)

    def scan_check(self, alpha, word):
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        f = alpha
        i = 0
        r = len(word)
        b = alpha
        j = r - 1
        while i <= j and table[f][A_dict[word[i]]] is not None:
            f = table[f][A_dict[word[i]]]
            i += 1
        if i > j:
            return f == b
        while j >= i and table[b][A_dict_inv[word[j]]] is not None:
            b = table[b][A_dict_inv[word[j]]]
            j -= 1
        if j < i:
            return False
        elif j == i:
            table[f][A_dict[word[i]]] = b
            table[b][A_dict_inv[word[i]]] = f
        return True

    def merge(self, k, lamda, q, w=None, modified=False):
        p = self.p
        rep = self.rep
        phi = rep(k, modified=modified)
        psi = rep(lamda, modified=modified)
        if phi != psi:
            mu = min(phi, psi)
            v = max(phi, psi)
            p[v] = mu
            if modified:
                if v == phi:
                    self.p_p[phi] = self.p_p[k] ** (-1) * w * self.p_p[lamda]
                else:
                    self.p_p[psi] = self.p_p[lamda] ** (-1) * w ** (-1) * self.p_p[k]
            q.append(v)

    def rep(self, k, modified=False):
        p = self.p
        lamda = k
        rho = p[lamda]
        if modified:
            s = p[:]
        while rho != lamda:
            if modified:
                s[rho] = lamda
            lamda = rho
            rho = p[lamda]
        if modified:
            rho = s[lamda]
            while rho != k:
                mu = rho
                rho = s[mu]
                p[rho] = lamda
                self.p_p[rho] = self.p_p[rho] * self.p_p[mu]
        else:
            mu = k
            rho = p[mu]
            while rho != lamda:
                p[mu] = lamda
                mu = rho
                rho = p[mu]
        return lamda

    def coincidence(self, alpha, beta, w=None, modified=False):
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        q = []
        if modified:
            self.modified_merge(alpha, beta, w, q)
        else:
            self.merge(alpha, beta, q)
        while len(q) > 0:
            gamma = q.pop(0)
            for x in A_dict:
                delta = table[gamma][A_dict[x]]
                if delta is not None:
                    table[delta][A_dict_inv[x]] = None
                    mu = self.rep(gamma, modified=modified)
                    nu = self.rep(delta, modified=modified)
                    if table[mu][A_dict[x]] is not None:
                        if modified:
                            v = self.p_p[delta] ** (-1) * self.P[gamma][self.A_dict[x]] ** (-1)
                            v = v * self.p_p[gamma] * self.P[mu][self.A_dict[x]]
                            self.modified_merge(nu, table[mu][self.A_dict[x]], v, q)
                        else:
                            self.merge(nu, table[mu][A_dict[x]], q)
                    elif table[nu][A_dict_inv[x]] is not None:
                        if modified:
                            v = self.p_p[gamma] ** (-1) * self.P[gamma][self.A_dict[x]]
                            v = v * self.p_p[delta] * self.P[mu][self.A_dict_inv[x]]
                            self.modified_merge(mu, table[nu][self.A_dict_inv[x]], v, q)
                        else:
                            self.merge(mu, table[nu][A_dict_inv[x]], q)
                    else:
                        table[mu][A_dict[x]] = nu
                        table[nu][A_dict_inv[x]] = mu
                        if modified:
                            v = self.p_p[gamma] ** (-1) * self.P[gamma][self.A_dict[x]] * self.p_p[delta]
                            self.P[mu][self.A_dict[x]] = v
                            self.P[nu][self.A_dict_inv[x]] = v ** (-1)

    def scan_and_fill(self, alpha, word):
        self.scan(alpha, word, fill=True)

    def scan_and_fill_c(self, alpha, word):
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        r = len(word)
        f = alpha
        i = 0
        b = alpha
        j = r - 1
        while True:
            while i <= j and table[f][A_dict[word[i]]] is not None:
                f = table[f][A_dict[word[i]]]
                i += 1
            if i > j:
                if f != b:
                    self.coincidence_c(f, b)
                return
            while j >= i and table[b][A_dict_inv[word[j]]] is not None:
                b = table[b][A_dict_inv[word[j]]]
                j -= 1
            if j < i:
                self.coincidence_c(f, b)
            elif j == i:
                table[f][A_dict[word[i]]] = b
                table[b][A_dict_inv[word[i]]] = f
                self.deduction_stack.append((f, word[i]))
            else:
                self.define_c(f, word[i])

    def look_ahead(self):
        R = self.fp_group.relators
        p = self.p
        for beta in self.omega:
            for w in R:
                self.scan(beta, w)
                if p[beta] < beta:
                    break

    def process_deductions(self, R_c_x, R_c_x_inv):
        p = self.p
        table = self.table
        while len(self.deduction_stack) > 0:
            if len(self.deduction_stack) >= CosetTable.max_stack_size:
                self.look_ahead()
                del self.deduction_stack[:]
                continue
            else:
                (alpha, x) = self.deduction_stack.pop()
                if p[alpha] == alpha:
                    for w in R_c_x:
                        self.scan_c(alpha, w)
                        if p[alpha] < alpha:
                            break
            beta = table[alpha][self.A_dict[x]]
            if beta is not None and p[beta] == beta:
                for w in R_c_x_inv:
                    self.scan_c(beta, w)
                    if p[beta] < beta:
                        break

    def process_deductions_check(self, R_c_x, R_c_x_inv):
        table = self.table
        while len(self.deduction_stack) > 0:
            (alpha, x) = self.deduction_stack.pop()
            if not all((self.scan_check(alpha, w) for w in R_c_x)):
                return False
            beta = table[alpha][self.A_dict[x]]
            if beta is not None:
                if not all((self.scan_check(beta, w) for w in R_c_x_inv)):
                    return False
        return True

    def switch(self, beta, gamma):
        A = self.A
        A_dict = self.A_dict
        table = self.table
        for x in A:
            z = table[gamma][A_dict[x]]
            table[gamma][A_dict[x]] = table[beta][A_dict[x]]
            table[beta][A_dict[x]] = z
            for alpha in range(len(self.p)):
                if self.p[alpha] == alpha:
                    if table[alpha][A_dict[x]] == beta:
                        table[alpha][A_dict[x]] = gamma
                    elif table[alpha][A_dict[x]] == gamma:
                        table[alpha][A_dict[x]] = beta

    def standardize(self):
        A = self.A
        A_dict = self.A_dict
        gamma = 1
        for (alpha, x) in product(range(self.n), A):
            beta = self.table[alpha][A_dict[x]]
            if beta >= gamma:
                if beta > gamma:
                    self.switch(gamma, beta)
                gamma += 1
                if gamma == self.n:
                    return

    def compress(self):
        gamma = -1
        A = self.A
        A_dict = self.A_dict
        A_dict_inv = self.A_dict_inv
        table = self.table
        chi = tuple([i for i in range(len(self.p)) if self.p[i] != i])
        for alpha in self.omega:
            gamma += 1
            if gamma != alpha:
                for x in A:
                    beta = table[alpha][A_dict[x]]
                    table[gamma][A_dict[x]] = beta
                    table[beta][A_dict_inv[x]] == gamma
        self.p = list(range(gamma + 1))
        del table[len(self.p):]
        for row in table:
            for j in range(len(self.A)):
                row[j] -= bisect_left(chi, row[j])

    def conjugates(self, R):
        R_c = list(chain.from_iterable(((rel.cyclic_conjugates(), (rel ** (-1)).cyclic_conjugates()) for rel in R)))
        R_set = set()
        for conjugate in R_c:
            R_set = R_set.union(conjugate)
        R_c_list = []
        for x in self.A:
            r = {word for word in R_set if word[0] == x}
            R_c_list.append(r)
            R_set.difference_update(r)
        return R_c_list

    def coset_representative(self, coset):
        for x in self.A:
            gamma = self.table[coset][self.A_dict[x]]
            if coset == 0:
                return self.fp_group.identity
            if gamma < coset:
                return self.coset_representative(gamma) * x ** (-1)

    def modified_define(self, alpha, x):
        self.define(alpha, x, modified=True)

    def modified_scan(self, alpha, w, y, fill=False):
        self.scan(alpha, w, y=y, fill=fill, modified=True)

    def modified_scan_and_fill(self, alpha, w, y):
        self.modified_scan(alpha, w, y, fill=True)

    def modified_merge(self, k, lamda, w, q):
        self.merge(k, lamda, q, w=w, modified=True)

    def modified_rep(self, k):
        self.rep(k, modified=True)

    def modified_coincidence(self, alpha, beta, w):
        self.coincidence(alpha, beta, w=w, modified=True)

def coset_enumeration_r(fp_grp, Y, max_cosets=None, draft=None, incomplete=False, modified=False):
    C = CosetTable(fp_grp, Y, max_cosets=max_cosets)
    if modified:
        _scan_and_fill = C.modified_scan_and_fill
        _define = C.modified_define
    else:
        _scan_and_fill = C.scan_and_fill
        _define = C.define
    if draft:
        C.table = draft.table[:]
        C.p = draft.p[:]
    R = fp_grp.relators
    A_dict = C.A_dict
    p = C.p
    for i in range(len(Y)):
        if modified:
            _scan_and_fill(0, Y[i], C._grp.generators[i])
        else:
            _scan_and_fill(0, Y[i])
    alpha = 0
    while alpha < C.n:
        if p[alpha] == alpha:
            try:
                for w in R:
                    if modified:
                        _scan_and_fill(alpha, w, C._grp.identity)
                    else:
                        _scan_and_fill(alpha, w)
                    if p[alpha] < alpha:
                        break
                if p[alpha] == alpha:
                    for x in A_dict:
                        if C.table[alpha][A_dict[x]] is None:
                            _define(alpha, x)
            except ValueError as e:
                if incomplete:
                    return C
                raise e
        alpha += 1
    return C

def modified_coset_enumeration_r(fp_grp, Y, max_cosets=None, draft=None, incomplete=False):
    return coset_enumeration_r(fp_grp, Y, max_cosets=max_cosets, draft=draft, incomplete=incomplete, modified=True)

def coset_enumeration_c(fp_grp, Y, max_cosets=None, draft=None, incomplete=False):
    X = fp_grp.generators
    R = fp_grp.relators
    C = CosetTable(fp_grp, Y, max_cosets=max_cosets)
    if draft:
        C.table = draft.table[:]
        C.p = draft.p[:]
        C.deduction_stack = draft.deduction_stack
        for (alpha, x) in product(range(len(C.table)), X):
            if C.table[alpha][C.A_dict[x]] is not None:
                C.deduction_stack.append((alpha, x))
    A = C.A
    R_cyc_red = [rel.identity_cyclic_reduction() for rel in R]
    R_c = list(chain.from_iterable(((rel.cyclic_conjugates(), (rel ** (-1)).cyclic_conjugates()) for rel in R_cyc_red)))
    R_set = set()
    for conjugate in R_c:
        R_set = R_set.union(conjugate)
    R_c_list = []
    for x in C.A:
        r = {word for word in R_set if word[0] == x}
        R_c_list.append(r)
        R_set.difference_update(r)
    for w in Y:
        C.scan_and_fill_c(0, w)
    for x in A:
        C.process_deductions(R_c_list[C.A_dict[x]], R_c_list[C.A_dict_inv[x]])
    alpha = 0
    while alpha < len(C.table):
        if C.p[alpha] == alpha:
            try:
                for x in C.A:
                    if C.p[alpha] != alpha:
                        break
                    if C.table[alpha][C.A_dict[x]] is None:
                        C.define_c(alpha, x)
                        C.process_deductions(R_c_list[C.A_dict[x]], R_c_list[C.A_dict_inv[x]])
            except ValueError as e:
                if incomplete:
                    return C
                raise e
        alpha += 1
    return C
