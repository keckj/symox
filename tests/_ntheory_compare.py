"""Shared helpers: evaluate one expression string with SymPy and symox."""
import itertools
from collections import defaultdict

import sympy as sp
import symox as sx


def _namespace(mod):
    ns = {'itertools': itertools}
    exec('from %s import *' % mod.__name__, ns)
    exec('from %s.ntheory import *' % mod.__name__, ns)
    exec('from %s.ntheory.modular import crt, crt1, crt2, solve_congruence'
         % mod.__name__, ns)
    exec('from %s.ntheory.factor_ import (udivisors, udivisor_count, '
         'antidivisors, antidivisor_count, core, smoothness, '
         'find_carmichael_numbers_in_range, find_first_n_carmichaels)'
         % mod.__name__, ns)
    exec('from %s.discrete import *' % mod.__name__, ns)
    exec('from %s.functions.combinatorial.numbers import udivisor_sigma'
         % mod.__name__ if mod is sp else
         'from %s import udivisor_sigma' % mod.__name__, ns)
    ns['M'] = mod
    for name in 'xyzn':
        ns[name] = mod.Symbol(name)
    ns['k'] = mod.Symbol('k', integer=True)
    ns['kp'] = mod.Symbol('kp', integer=True, positive=True)
    ns['kn'] = mod.Symbol('kn', integer=True, negative=True)
    ns['po'] = mod.Symbol('po', positive=True)
    ns['pr'] = mod.Symbol('pr', prime=True)
    ns['ev'] = mod.Symbol('ev', even=True)
    ns['od'] = mod.Symbol('od', odd=True)
    return ns


NS = {'sympy': _namespace(sp), 'symox': _namespace(sx)}


def norm(v):
    """A library-independent description of a result."""
    if isinstance(v, bool) or v is None:
        return ('py', type(v).__name__, v)
    if isinstance(v, int):
        return ('py', 'int', v)
    if isinstance(v, float):
        return ('py', 'float', repr(v))
    if isinstance(v, (list, tuple)):
        return (type(v).__name__, [norm(i) for i in v])
    if isinstance(v, defaultdict):
        return ('defaultdict', [(norm(a), norm(b)) for a, b in v.items()])
    if isinstance(v, dict):
        return ('dict', [(norm(a), norm(b)) for a, b in v.items()])
    if isinstance(v, (set, frozenset)):
        return ('set', sorted(str(norm(i)) for i in v))
    if hasattr(v, '__next__'):
        return ('iter', [norm(i) for i in itertools.islice(v, 50)])
    kind = 'Expr'
    for k in ('Integer', 'Rational', 'Float', 'Symbol', 'Add', 'Mul', 'Pow',
              'Function'):
        if getattr(v, 'is_' + k, False) is True:
            kind = k
            break
    return ('expr', kind, str(v))


def run(lib, src):
    try:
        return norm(eval(src, NS[lib]))
    except Exception as e:  # noqa: BLE001
        msg = str(e).splitlines()[0] if str(e) else ''
        return ('exc', type(e).__name__, msg)


def check(src):
    expected = run('sympy', src)
    got = run('symox', src)
    assert got == expected, src
