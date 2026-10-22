"""Free groups and reduced words (ported from the reference free_groups module)."""
from __future__ import annotations
from symox.core import S
from symox.core.expr import Expr
from symox.core.symbol import Symbol, symbols as _symbols
from symox.core.sympify import CantSympify
from symox.printing.defaults import DefaultPrinting
from symox.utilities import public
from symox.utilities.iterables import flatten, is_sequence
from symox.utilities.magic import pollute
from symox.utilities.misc import as_int

@public
def free_group(symbols):
    _free_group = FreeGroup(symbols)
    return (_free_group,) + tuple(_free_group.generators)

@public
def xfree_group(symbols):
    _free_group = FreeGroup(symbols)
    return (_free_group, _free_group.generators)

@public
def vfree_group(symbols):
    _free_group = FreeGroup(symbols)
    pollute([sym.name for sym in _free_group.symbols], _free_group.generators)
    return _free_group

def _parse_symbols(symbols):
    if not symbols:
        return ()
    if isinstance(symbols, str):
        return _symbols(symbols, seq=True)
    elif isinstance(symbols, (Expr, FreeGroupElement)):
        return (symbols,)
    elif is_sequence(symbols):
        if all((isinstance(s, str) for s in symbols)):
            return _symbols(symbols)
        elif all((isinstance(s, Expr) for s in symbols)):
            return symbols
    raise ValueError('The type of `symbols` must be one of the following: a str, Symbol/Expr or a sequence of one of these types')
_free_group_cache: dict[int, FreeGroup] = {}

class FreeGroup(DefaultPrinting):
    is_associative = True
    is_group = True
    is_FreeGroup = True
    is_PermutationGroup = False
    relators: list[Expr] = []

    def __new__(cls, symbols):
        symbols = tuple(_parse_symbols(symbols))
        rank = len(symbols)
        _hash = hash((cls.__name__, symbols, rank))
        obj = _free_group_cache.get(_hash)
        if obj is None:
            obj = object.__new__(cls)
            obj._hash = _hash
            obj._rank = rank
            obj.dtype = type('FreeGroupElement', (FreeGroupElement,), {'group': obj})
            obj.symbols = symbols
            obj.generators = obj._generators()
            obj._gens_set = set(obj.generators)
            for (symbol, generator) in zip(obj.symbols, obj.generators):
                if isinstance(symbol, Symbol):
                    name = symbol.name
                    if hasattr(obj, name):
                        setattr(obj, name, generator)
            _free_group_cache[_hash] = obj
        return obj

    def __getnewargs__(self):
        return (self.symbols,)

    def __getstate__(self):
        return None

    def _generators(group):
        gens = []
        for sym in group.symbols:
            elm = ((sym, 1),)
            gens.append(group.dtype(elm))
        return tuple(gens)

    def clone(self, symbols=None):
        return self.__class__(symbols or self.symbols)

    def __contains__(self, i):
        if not isinstance(i, FreeGroupElement):
            return False
        group = i.group
        return self == group

    def __hash__(self):
        return self._hash

    def __len__(self):
        return self.rank

    def __str__(self):
        if self.rank > 30:
            str_form = '<free group with %s generators>' % self.rank
        else:
            str_form = '<free group on the generators '
            gens = self.generators
            str_form += str(gens) + '>'
        return str_form
    __repr__ = __str__

    def __getitem__(self, index):
        symbols = self.symbols[index]
        return self.clone(symbols=symbols)

    def __eq__(self, other):
        return self is other

    def index(self, gen):
        if isinstance(gen, self.dtype):
            return self.generators.index(gen)
        else:
            raise ValueError('expected a generator of Free Group %s, got %s' % (self, gen))

    def order(self):
        if self.rank == 0:
            return S.One
        else:
            return S.Infinity

    @property
    def elements(self):
        if self.rank == 0:
            return {self.identity}
        else:
            raise ValueError('Group contains infinitely many elements, hence cannot be represented')

    @property
    def rank(self):
        return self._rank

    @property
    def is_abelian(self):
        return self.rank in (0, 1)

    @property
    def identity(self):
        return self.dtype()

    def contains(self, g):
        if not isinstance(g, FreeGroupElement):
            return False
        elif self != g.group:
            return False
        else:
            return True

    def center(self):
        return {self.identity}

class FreeGroupElement(CantSympify, DefaultPrinting, tuple):
    __slots__ = ()
    is_assoc_word = True

    def new(self, init):
        return self.__class__(init)
    _hash = None

    def __hash__(self):
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.group, frozenset(tuple(self))))
        return _hash

    def copy(self):
        return self.new(self)

    @property
    def is_identity(self):
        return not self.array_form

    @property
    def array_form(self):
        return tuple(self)

    @property
    def letter_form(self):
        return tuple(s for (i, j) in self.array_form for s in ((i,) * j if j > 0 else (-i,) * -j))

    def __getitem__(self, i):
        group = self.group
        r = self.letter_form[i]
        if r.is_Symbol:
            return group.dtype(((r, 1),))
        else:
            return group.dtype(((-r, -1),))

    def index(self, gen):
        if len(gen) != 1:
            raise ValueError()
        return self.letter_form.index(gen.letter_form[0])

    @property
    def letter_form_elm(self):
        group = self.group
        r = self.letter_form
        return [group.dtype(((elm, 1),)) if elm.is_Symbol else group.dtype(((-elm, -1),)) for elm in r]

    @property
    def ext_rep(self):
        return tuple(x for pair in self.array_form for x in pair)

    def __contains__(self, gen):
        return gen.array_form[0][0] in tuple([r[0] for r in self.array_form])

    def __str__(self):
        if self.is_identity:
            return '<identity>'
        str_form = ''
        array_form = self.array_form
        for i in range(len(array_form)):
            if i == len(array_form) - 1:
                if array_form[i][1] == 1:
                    str_form += str(array_form[i][0])
                else:
                    str_form += str(array_form[i][0]) + '**' + str(array_form[i][1])
            elif array_form[i][1] == 1:
                str_form += str(array_form[i][0]) + '*'
            else:
                str_form += str(array_form[i][0]) + '**' + str(array_form[i][1]) + '*'
        return str_form
    __repr__ = __str__

    def __pow__(self, n):
        n = as_int(n)
        result = self.group.identity
        if n == 0:
            return result
        if n < 0:
            n = -n
            x = self.inverse()
        else:
            x = self
        while True:
            if n % 2:
                result *= x
            n >>= 1
            if not n:
                break
            x *= x
        return result

    def __mul__(self, other):
        group = self.group
        if not isinstance(other, group.dtype):
            raise TypeError('only FreeGroup elements of same FreeGroup can be multiplied')
        if self.is_identity:
            return other
        if other.is_identity:
            return self
        r = list(self.array_form + other.array_form)
        zero_mul_simp(r, len(self.array_form) - 1)
        return group.dtype(tuple(r))

    def __truediv__(self, other):
        group = self.group
        if not isinstance(other, group.dtype):
            raise TypeError('only FreeGroup elements of same FreeGroup can be multiplied')
        return self * other.inverse()

    def __rtruediv__(self, other):
        group = self.group
        if not isinstance(other, group.dtype):
            raise TypeError('only FreeGroup elements of same FreeGroup can be multiplied')
        return other * self.inverse()

    def __add__(self, other):
        return NotImplemented

    def inverse(self):
        group = self.group
        r = tuple([(i, -j) for (i, j) in self.array_form[::-1]])
        return group.dtype(r)

    def order(self):
        if self.is_identity:
            return S.One
        else:
            return S.Infinity

    def commutator(self, other):
        group = self.group
        if not isinstance(other, group.dtype):
            raise ValueError('commutator of only FreeGroupElement of the same FreeGroup exists')
        else:
            return self.inverse() * other.inverse() * self * other

    def eliminate_words(self, words, _all=False, inverse=True):
        again = True
        new = self
        if isinstance(words, dict):
            while again:
                again = False
                for sub in words:
                    prev = new
                    new = new.eliminate_word(sub, words[sub], _all=_all, inverse=inverse)
                    if new != prev:
                        again = True
        else:
            while again:
                again = False
                for sub in words:
                    prev = new
                    new = new.eliminate_word(sub, _all=_all, inverse=inverse)
                    if new != prev:
                        again = True
        return new

    def eliminate_word(self, gen, by=None, _all=False, inverse=True):
        if by is None:
            by = self.group.identity
        if self.is_independent(gen) or gen == by:
            return self
        if gen == self:
            return by
        if gen ** (-1) == by:
            _all = False
        word = self
        l = len(gen)
        try:
            i = word.subword_index(gen)
            k = 1
        except ValueError:
            if not inverse:
                return word
            try:
                i = word.subword_index(gen ** (-1))
                k = -1
            except ValueError:
                return word
        word = word.subword(0, i) * by ** k * word.subword(i + l, len(word)).eliminate_word(gen, by)
        if _all:
            return word.eliminate_word(gen, by, _all=True, inverse=inverse)
        else:
            return word

    def __len__(self):
        return sum((abs(j) for (i, j) in self))

    def __eq__(self, other):
        group = self.group
        if not isinstance(other, group.dtype):
            return False
        return tuple.__eq__(self, other)

    def __lt__(self, other):
        group = self.group
        if not isinstance(other, group.dtype):
            raise TypeError('only FreeGroup elements of same FreeGroup can be compared')
        l = len(self)
        m = len(other)
        if l < m:
            return True
        elif l > m:
            return False
        for i in range(l):
            a = self[i].array_form[0]
            b = other[i].array_form[0]
            p = group.symbols.index(a[0])
            q = group.symbols.index(b[0])
            if p < q:
                return True
            elif p > q:
                return False
            elif a[1] < b[1]:
                return True
            elif a[1] > b[1]:
                return False
        return False

    def __le__(self, other):
        return self == other or self < other

    def __gt__(self, other):
        group = self.group
        if not isinstance(other, group.dtype):
            raise TypeError('only FreeGroup elements of same FreeGroup can be compared')
        return not self <= other

    def __ge__(self, other):
        return not self < other

    def exponent_sum(self, gen):
        if len(gen) != 1:
            raise ValueError('gen must be a generator or inverse of a generator')
        s = gen.array_form[0]
        return s[1] * sum((i[1] for i in self.array_form if i[0] == s[0]))

    def generator_count(self, gen):
        if len(gen) != 1 or gen.array_form[0][1] < 0:
            raise ValueError('gen must be a generator')
        s = gen.array_form[0]
        return s[1] * sum((abs(i[1]) for i in self.array_form if i[0] == s[0]))

    def subword(self, from_i, to_j, strict=True):
        group = self.group
        if not strict:
            from_i = max(from_i, 0)
            to_j = min(len(self), to_j)
        if from_i < 0 or to_j > len(self):
            raise ValueError('`from_i`, `to_j` must be positive and no greater than the length of associative word')
        if to_j <= from_i:
            return group.identity
        else:
            letter_form = self.letter_form[from_i:to_j]
            array_form = letter_form_to_array_form(letter_form, group)
            return group.dtype(array_form)

    def subword_index(self, word, start=0):
        l = len(word)
        self_lf = self.letter_form
        word_lf = word.letter_form
        index = None
        for i in range(start, len(self_lf) - l + 1):
            if self_lf[i:i + l] == word_lf:
                index = i
                break
        if index is not None:
            return index
        else:
            raise ValueError('The given word is not a subword of self')

    def is_dependent(self, word):
        try:
            return self.subword_index(word) is not None
        except ValueError:
            pass
        try:
            return self.subword_index(word ** (-1)) is not None
        except ValueError:
            return False

    def is_independent(self, word):
        return not self.is_dependent(word)

    def contains_generators(self):
        group = self.group
        gens = {group.dtype(((syllable[0], 1),)) for syllable in self.array_form}
        return gens

    def cyclic_subword(self, from_i, to_j):
        group = self.group
        l = len(self)
        letter_form = self.letter_form
        period1 = int(from_i / l)
        if from_i >= l:
            from_i -= l * period1
            to_j -= l * period1
        diff = to_j - from_i
        word = letter_form[from_i:to_j]
        period2 = int(to_j / l) - 1
        word += letter_form * period2 + letter_form[:diff - l + from_i - l * period2]
        word = letter_form_to_array_form(word, group)
        return group.dtype(word)

    def cyclic_conjugates(self):
        return {self.cyclic_subword(i, i + len(self)) for i in range(len(self))}

    def is_cyclic_conjugate(self, w):
        l1 = len(self)
        l2 = len(w)
        if l1 != l2:
            return False
        w1 = self.identity_cyclic_reduction()
        w2 = w.identity_cyclic_reduction()
        letter1 = w1.letter_form
        letter2 = w2.letter_form
        str1 = ' '.join(map(str, letter1))
        str2 = ' '.join(map(str, letter2))
        if len(str1) != len(str2):
            return False
        return str1 in str2 + ' ' + str2

    def number_syllables(self):
        return len(self.array_form)

    def exponent_syllable(self, i):
        return self.array_form[i][1]

    def generator_syllable(self, i):
        return self.array_form[i][0]

    def sub_syllables(self, from_i, to_j):
        if not isinstance(from_i, int) or not isinstance(to_j, int):
            raise ValueError('both arguments should be integers')
        group = self.group
        if to_j <= from_i:
            return group.identity
        else:
            r = tuple(self.array_form[from_i:to_j])
            return group.dtype(r)

    def substituted_word(self, from_i, to_j, by):
        lw = len(self)
        if from_i >= to_j or from_i > lw or to_j > lw:
            raise ValueError('values should be within bounds')
        if from_i == 0 and to_j == lw:
            return by
        elif from_i == 0:
            return by * self.subword(to_j, lw)
        elif to_j == lw:
            return self.subword(0, from_i) * by
        else:
            return self.subword(0, from_i) * by * self.subword(to_j, lw)

    def is_cyclically_reduced(self):
        if not self:
            return True
        return self[0] != self[-1] ** (-1)

    def identity_cyclic_reduction(self):
        word = self.copy()
        group = self.group
        while not word.is_cyclically_reduced():
            exp1 = word.exponent_syllable(0)
            exp2 = word.exponent_syllable(-1)
            r = exp1 + exp2
            if r == 0:
                rep = word.array_form[1:word.number_syllables() - 1]
            else:
                rep = ((word.generator_syllable(0), exp1 + exp2),) + word.array_form[1:word.number_syllables() - 1]
            word = group.dtype(rep)
        return word

    def cyclic_reduction(self, removed=False):
        word = self.copy()
        g = self.group.identity
        while not word.is_cyclically_reduced():
            exp1 = abs(word.exponent_syllable(0))
            exp2 = abs(word.exponent_syllable(-1))
            exp = min(exp1, exp2)
            start = word[0] ** abs(exp)
            end = word[-1] ** abs(exp)
            word = start ** (-1) * word * end ** (-1)
            g = g * start
        if removed:
            return (word, g)
        return word

    def power_of(self, other):
        if self.is_identity:
            return True
        l = len(other)
        if l == 1:
            gens = self.contains_generators()
            s = other in gens or other ** (-1) in gens
            return len(gens) == 1 and s
        (reduced, r1) = self.cyclic_reduction(removed=True)
        if not r1.is_identity:
            (other, r2) = other.cyclic_reduction(removed=True)
            if r1 == r2:
                return reduced.power_of(other)
            return False
        if len(self) < l or len(self) % l:
            return False
        prefix = self.subword(0, l)
        if prefix == other or prefix ** (-1) == other:
            rest = self.subword(l, len(self))
            return rest.power_of(other)
        return False

def letter_form_to_array_form(array_form, group):
    a = list(array_form[:])
    new_array = []
    n = 1
    symbols = group.symbols
    for i in range(len(a)):
        if i == len(a) - 1:
            if a[i] == a[i - 1]:
                if -a[i] in symbols:
                    new_array.append((-a[i], -n))
                else:
                    new_array.append((a[i], n))
            elif -a[i] in symbols:
                new_array.append((-a[i], -1))
            else:
                new_array.append((a[i], 1))
            return new_array
        elif a[i] == a[i + 1]:
            n += 1
        else:
            if -a[i] in symbols:
                new_array.append((-a[i], -n))
            else:
                new_array.append((a[i], n))
            n = 1

def zero_mul_simp(l, index):
    while index >= 0 and index < len(l) - 1 and (l[index][0] == l[index + 1][0]):
        exp = l[index][1] + l[index + 1][1]
        base = l[index][0]
        l[index] = (base, exp)
        del l[index + 1]
        if l[index][1] == 0:
            del l[index]
            index -= 1


def _free_element_to_sympy(self):
    """Convert a word to the corresponding word of a reference free group."""
    from sympy.combinatorics import free_group
    group, *gens = free_group(','.join(map(str, self.group.symbols)))
    value = group.identity
    for s, e in self.array_form:
        value *= gens[self.group.symbols.index(s)]**e
    return value


FreeGroupElement.to_sympy = _free_element_to_sympy
