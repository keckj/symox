"""Decimal strings of mpmath ``_mpf_`` tuples without mpmath.

``prec_to_dps`` and ``to_str`` follow mpmath/libmp/libmpf.py, which the
reference ``StrPrinter._print_Float`` calls: the binary value is truncated
to a decimal fixed-point number with ``dps + 3`` digits, then rounded half
up to ``dps`` digits. Values whose binary exponent exceeds 3500 in
magnitude are first scaled by a power of ten computed with truncated
binary exponentiation, as mpmath does.
"""
import math

_LOG2_10 = math.log(10, 2)


def prec_to_dps(n):
    """Number of accurate decimals that can be represented with a
    precision of ``n`` bits."""
    return max(1, int(round(int(n) / 3.3219280948873626) - 1))


def _atanh_inverse(k, bits):
    """``atanh(1/k)`` as a fixed-point integer with ``bits`` fraction bits."""
    one = 1 << bits
    total, power, n, k2 = 0, one // k, 1, k * k
    while power:
        total += power // n
        power //= k2
        n += 2
    return total


def _log10_2(bits):
    """``log(2)/log(10)`` as a fixed-point integer with ``bits`` fraction
    bits, accurate to a few units."""
    work = bits + 20
    ln2 = 2 * _atanh_inverse(3, work)
    ln10 = 3 * ln2 + 2 * _atanh_inverse(9, work)
    return (ln2 << bits) // ln10


def _truncated_pow10(n, prec):
    """``10**n`` (``n > 0``) as ``(man, exp)``, truncated to ``prec`` bits at
    every step."""
    pm, pe, man, exp = 1, 0, 10, 0
    while True:
        if n & 1:
            pm, pe = pm * man, pe + exp
            extra = pm.bit_length() - prec
            if extra > 0:
                pm, pe = pm >> extra, pe + extra
            n -= 1
            if not n:
                return pm, pe
        man, exp = man * man, exp + exp
        extra = man.bit_length() - prec
        if extra > 0:
            man, exp = man >> extra, exp + extra
        n //= 2


def _scaled(man, exp, bitprec):
    """``(man, exp, exponent)`` with ``man*2**exp*10**exponent`` close to the
    input and a moderate binary exponent."""
    bits = abs(exp).bit_length() + 5
    b = abs(exp) * _log10_2(bits) >> bits
    if exp < 0:
        b = -b
    if b == 0:
        return man, exp, 0
    work = bitprec + 4 * abs(b).bit_length() + 10
    tm, te = _truncated_pow10(abs(b), work)
    if b > 0:
        shift = max(0, work + tm.bit_length() - man.bit_length())
        return (man << shift) // tm, exp - shift - te, b
    return man * tm, exp + te, b


def to_digits_exp(s, dps):
    """``(sign, digits, exponent)`` of ``s`` with about ``dps`` digits,
    rounded toward zero."""
    sign, man, exp, bc = s
    sign = '-' if sign else ''
    if not man:
        return '', '0', 0
    bitprec = int(dps * _LOG2_10) + 10
    exponent = 0
    if abs(exp + bc) > 3500:
        man, exp, exponent = _scaled(man, exp, bitprec)
        bc = man.bit_length()
    fixprec = max(bitprec - exp - bc, 0)
    fixdps = int(fixprec / _LOG2_10 + 0.5)
    offset = exp + fixprec
    fixed = man << offset if offset >= 0 else man >> -offset
    digits = str(fixed * 10**fixdps >> fixprec)
    exponent += len(digits) - fixdps - 1
    return sign, digits, exponent


def to_str(s, dps, strip_zeros=True, min_fixed=None, max_fixed=None, show_zero_exponent=False):
    """mpmath's decimal string of ``s`` with ``dps`` significant digits."""
    if not s[1]:
        if s == (0, 0, 0, 0):
            t = '0.0' if dps else '.0'
            return t + 'e+0' if show_zero_exponent else t
        if s == (0, 0, -456, -2):
            return '+inf'
        if s == (1, 0, -789, -3):
            return '-inf'
        if s == (0, 0, -123, -1):
            return 'nan'
        raise ValueError
    if min_fixed is None:
        min_fixed = min(-(dps // 3), -5)
    if max_fixed is None:
        max_fixed = dps
    sign, digits, exponent = to_digits_exp(s, dps + 3)
    if not dps:
        if digits[0] in '56789':
            exponent += 1
        digits = '.0'
    else:
        if len(digits) > dps and digits[dps] in '56789':
            digits = digits[:dps]
            i = dps - 1
            while i >= 0 and digits[i] == '9':
                i -= 1
            if i >= 0:
                digits = digits[:i] + str(int(digits[i]) + 1) + '0' * (dps - i - 1)
            else:
                digits = '1' + '0' * (dps - 1)
                exponent += 1
        else:
            digits = digits[:dps]
        if min_fixed < exponent < max_fixed:
            if exponent < 0:
                digits = '0' * int(-exponent) + digits
                split = 1
            else:
                split = exponent + 1
                if split > dps:
                    digits += '0' * (split - dps)
            exponent = 0
        else:
            split = 1
        digits = digits[:split] + '.' + digits[split:]
        if strip_zeros:
            digits = digits.rstrip('0')
            if digits[-1] == '.':
                digits += '0'
    if exponent == 0 and dps and not show_zero_exponent:
        return sign + digits
    if exponent >= 0:
        return sign + digits + 'e+' + str(exponent)
    return sign + digits + 'e' + str(exponent)
