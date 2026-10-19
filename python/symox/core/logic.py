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


__all__ = ['fuzzy_bool', 'fuzzy_not', 'fuzzy_and', 'fuzzy_or', 'fuzzy_xor', 'fuzzy_nand']
