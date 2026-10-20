"""Lazy three-valued Boolean operations used by symbolic assumptions."""


def fuzzy_bool(x):
    if x is not None and x in (True, False):
        return bool(x)
    return None


def fuzzy_not(v):
    return None if v is None else not v


def fuzzy_and(args):
    unknown = False
    for value in args:
        flag = fuzzy_bool(value)
        if flag is False:
            return False
        unknown |= flag is None
    return None if unknown else True


def fuzzy_or(args):
    unknown = False
    for value in args:
        flag = fuzzy_bool(value)
        if flag is True:
            return True
        unknown |= flag is None
    return None if unknown else False


def fuzzy_xor(args):
    parity = False
    for value in args:
        flag = fuzzy_bool(value)
        if flag is None:
            return None
        parity ^= flag
    return parity


def fuzzy_nand(args):
    return fuzzy_not(fuzzy_and(args))


def _fuzzy_group(args, quick_exit=False):
    false_seen = False
    for value in args:
        if value is True:
            continue
        if value is None or quick_exit and false_seen:
            return None
        false_seen = True
    return not false_seen


def _torf(args):
    common = None
    for value in args:
        if value is not True and value is not False:
            return None
        if common is not None and value is not common:
            return None
        common = value
    return common is True



# Propositional rule expressions used by the assumption fact compiler
# (sympy/core/logic.py, SymPy 1.14).
class Logic:
    """Logical expression"""
    # {} 'op' -> LogicClass
    op_2class = {}

    def __new__(cls, *args):
        obj = object.__new__(cls)
        obj.args = args
        return obj

    def __getnewargs__(self):
        return self.args

    def __hash__(self):
        return hash((type(self).__name__,) + tuple(self.args))

    def __eq__(a, b):
        if not isinstance(b, type(a)):
            return False
        else:
            return a.args == b.args

    def __ne__(a, b):
        if not isinstance(b, type(a)):
            return True
        else:
            return a.args != b.args

    def __lt__(self, other):
        if self.__cmp__(other) == -1:
            return True
        return False

    def __cmp__(self, other):
        if type(self) is not type(other):
            a = str(type(self))
            b = str(type(other))
        else:
            a = self.args
            b = other.args
        return (a > b) - (a < b)

    def __str__(self):
        return '%s(%s)' % (self.__class__.__name__,
                           ', '.join(str(a) for a in self.args))

    __repr__ = __str__

    @staticmethod
    def fromstring(text):
        """Logic from string with space around & and | but none after !.

           e.g.

           !a & b | c
        """
        lexpr = None  # current logical expression
        schedop = None  # scheduled operation
        for term in text.split():
            # operation symbol
            if term in '&|':
                if schedop is not None:
                    raise ValueError(
                        'double op forbidden: "%s %s"' % (term, schedop))
                if lexpr is None:
                    raise ValueError(
                        '%s cannot be in the beginning of expression' % term)
                schedop = term
                continue
            if '&' in term or '|' in term:
                raise ValueError('& and | must have space around them')
            if term[0] == '!':
                if len(term) == 1:
                    raise ValueError('do not include space after "!"')
                term = Not(term[1:])

            # already scheduled operation, e.g. '&'
            if schedop:
                lexpr = Logic.op_2class[schedop](lexpr, term)
                schedop = None
                continue

            # this should be atom
            if lexpr is not None:
                raise ValueError(
                    'missing op between "%s" and "%s"' % (lexpr, term))

            lexpr = term

        # let's check that we ended up in correct state
        if schedop is not None:
            raise ValueError('premature end-of-expression in "%s"' % text)
        if lexpr is None:
            raise ValueError('"%s" is empty' % text)

        # everything looks good now
        return lexpr


class AndOr_Base(Logic):

    def __new__(cls, *args):
        bargs = []
        for a in args:
            if a == cls.op_x_notx:
                return a
            elif a == (not cls.op_x_notx):
                continue    # skip this argument
            bargs.append(a)

        args = sorted(set(cls.flatten(bargs)), key=hash)

        for a in args:
            if Not(a) in args:
                return cls.op_x_notx

        if len(args) == 1:
            return args.pop()
        elif len(args) == 0:
            return not cls.op_x_notx

        return Logic.__new__(cls, *args)

    @classmethod
    def flatten(cls, args):
        # quick-n-dirty flattening for And and Or
        args_queue = list(args)
        res = []

        while True:
            try:
                arg = args_queue.pop(0)
            except IndexError:
                break
            if isinstance(arg, Logic):
                if isinstance(arg, cls):
                    args_queue.extend(arg.args)
                    continue
            res.append(arg)

        args = tuple(res)
        return args


class And(AndOr_Base):
    op_x_notx = False

    def _eval_propagate_not(self):
        # !(a&b&c ...) == !a | !b | !c ...
        return Or(*[Not(a) for a in self.args])

    # (a|b|...) & c == (a&c) | (b&c) | ...
    def expand(self):

        # first locate Or
        for i, arg in enumerate(self.args):
            if isinstance(arg, Or):
                arest = self.args[:i] + self.args[i + 1:]

                orterms = [And(*(arest + (a,))) for a in arg.args]
                for j in range(len(orterms)):
                    if isinstance(orterms[j], Logic):
                        orterms[j] = orterms[j].expand()

                res = Or(*orterms)
                return res

        return self


class Or(AndOr_Base):
    op_x_notx = True

    def _eval_propagate_not(self):
        # !(a|b|c ...) == !a & !b & !c ...
        return And(*[Not(a) for a in self.args])


class Not(Logic):

    def __new__(cls, arg):
        if isinstance(arg, str):
            return Logic.__new__(cls, arg)

        elif isinstance(arg, bool):
            return not arg
        elif isinstance(arg, Not):
            return arg.args[0]

        elif isinstance(arg, Logic):
            # XXX this is a hack to expand right from the beginning
            arg = arg._eval_propagate_not()
            return arg

        else:
            raise ValueError('Not: unknown argument %r' % (arg,))

    @property
    def arg(self):
        return self.args[0]


Logic.op_2class['&'] = And
Logic.op_2class['|'] = Or
Logic.op_2class['!'] = Not


__all__ = ['fuzzy_bool', 'fuzzy_not', 'fuzzy_and', 'fuzzy_or', 'fuzzy_xor', 'fuzzy_nand']
