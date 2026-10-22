"""Static helpers of SymPy 1.14's ``motzkin`` class.

Ported from ``sympy/functions/combinatorial/numbers.py`` (BSD-3-Clause, see
LICENSE-SYMPY). The arithmetic is kept identical, including the true
divisions of the reference recurrences.
"""


def _as_int(n):
    from .utilities.misc import as_int
    return as_int(n)


def is_motzkin(n):
    try:
        n = _as_int(n)
    except ValueError:
        return False
    if n > 0:
        if n in (1, 2):
            return True
        tn1 = 1
        tn = 2
        i = 3
        while tn < n:
            a = ((2*i + 1)*tn + (3*i - 3)*tn1)/(i + 2)
            i += 1
            tn1 = tn
            tn = a
        return tn == n
    return False


def find_motzkin_numbers_in_range(x, y):
    if 0 <= x <= y:
        motzkins = []
        if x <= 1 <= y:
            motzkins.append(1)
        tn1 = 1
        tn = 2
        i = 3
        while tn <= y:
            if tn >= x:
                motzkins.append(tn)
            a = ((2*i + 1)*tn + (3*i - 3)*tn1)/(i + 2)
            i += 1
            tn1 = tn
            tn = int(a)
        return motzkins
    raise ValueError('The provided range is not valid. This condition should satisfy x <= y')


def find_first_n_motzkins(n):
    try:
        n = _as_int(n)
    except ValueError:
        raise ValueError('The provided number must be a positive integer')
    if n < 0:
        raise ValueError('The provided number must be a positive integer')
    motzkins = [1]
    if n >= 1:
        motzkins.append(1)
    tn1 = 1
    tn = 2
    i = 3
    while i <= n:
        motzkins.append(tn)
        a = ((2*i + 1)*tn + (3*i - 3)*tn1)/(i + 2)
        i += 1
        tn1 = tn
        tn = int(a)
    return motzkins


def install(head):
    for function in (is_motzkin, find_motzkin_numbers_in_range, find_first_n_motzkins):
        setattr(head, function.__name__, staticmethod(function))
