from .. import sympify,SympifyError
from .._sympify import converter
__all__ = ['sympify','SympifyError']

class CantSympify:
    """Marker for objects which must not be implicitly converted to scalars."""

def kernS(s):
    """Port of the reference kernS: a kern symbol keeps the 2-arg Muls of
    negated parenthesised expressions from distributing."""
    from random import choice
    import string
    from symox import Symbol, SympifyError
    hit = False
    quoted = '"' in s or "'" in s
    if '(' in s and not quoted:
        if s.count('(') != s.count(")"):
            raise SympifyError('unmatched left parenthesis')
        s = ''.join(s.split())
        olds = s
        s = s.replace('*(', '* *(')
        s = s.replace('** *', '**')
        target = '-( *('
        s = s.replace('-(', target)
        i = nest = 0
        while True:
            j = s.find(target, i)
            if j == -1:
                break
            j += len(target) - 1
            for j in range(j, len(s)):
                if s[j] == "(":
                    nest += 1
                elif s[j] == ")":
                    nest -= 1
                if nest == 0:
                    break
            s = s[:j] + ")" + s[j:]
            i = j + 2
        if ' ' in s:
            kern = '_'
            while kern in s:
                kern += choice(string.ascii_letters + string.digits)
            s = s.replace(' ', kern)
            hit = kern in s
        else:
            hit = False
    for i in range(2):
        try:
            expr = sympify(s)
            break
        except TypeError:
            if hit:
                s = olds
                hit = False
                continue
            expr = sympify(s)
    if not hit:
        return expr
    k = Symbol(kern)

    def _subs(e):
        # Basic._subs fallback with the reference hack2 2-arg Mul rule.
        if e == k:
            return sympify(1)
        if not getattr(e, 'args', ()) or not e.has(k):
            return e
        args = [_subs(a) for a in e.args]
        rv = e.func(*args)
        if e.is_Mul and not rv.is_Mul:
            coeff = sympify(1)
            nonnumber = []
            for a in args:
                if a.is_Number:
                    coeff *= a
                else:
                    nonnumber.append(a)
            nonnumber = e.func(*nonnumber)
            if coeff == 1:
                return nonnumber
            return e.func(coeff, nonnumber, evaluate=False)
        return rv

    def _clear(expr):
        if isinstance(expr, (list, tuple, set)):
            return type(expr)([_clear(e) for e in expr])
        if hasattr(expr, 'subs'):
            return _subs(expr)
        return expr
    return _clear(expr)
__all__ += ['CantSympify','kernS']
__all__ += ['converter']


def _sympify(a):
    """Short version of sympify for internal usage: ``sympify(a, strict=True)``.

    Only objects with a registered or implicit symbolic conversion are
    accepted; strings and other non-Basic objects raise SympifyError.
    """
    return sympify(a, strict=True)
__all__ += ['_sympify']
