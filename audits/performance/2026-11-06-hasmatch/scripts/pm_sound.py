"""Soundness of the old ``pattern_match`` bindings logged for one caller.

For each logged call (see pm_audit.py) where the old matcher returned a
binding, check with SymPy 1.14 that the pattern with that binding
substituted equals the expression (``expand(p.xreplace(m) - e) == 0``,
then ``simplify``).  A sound binding makes the caller's rewrite an
identity, whatever binding the reference would have chosen.

Usage: pm_sound.py LOG CALLER [LIMIT]
"""
import collections
import sys

import sympy


def parse(text):
    return sympy.parse_expr(text, evaluate=True)


def main(log, caller, limit=400):
    limit = int(limit)
    rows, seen = [], set()
    for line in open(log):
        parts = line.rstrip('\n').split('\t')
        if len(parts) != 5 or parts[0] != caller or parts[3] in ('None', 'ERROR'):
            continue
        if parts[3] == parts[4]:
            continue
        key = tuple(parts[1:4])
        if key not in seen:
            seen.add(key)
            rows.append(parts)
    counts = collections.Counter()
    bad = []
    for _, expr, pattern, old, _ in rows[:limit]:
        try:
            e, p = parse(expr), parse(pattern)
            binding = {}
            for item in old.split('; '):
                name, value = item.split('=', 1)
                binding[sympy.Symbol(name)] = parse(value)
            d = p.xreplace(binding) - e
            ok = sympy.expand(d) == 0 or sympy.simplify(d) == 0
        except Exception as error:
            counts['unparsed'] += 1
            continue
        counts['sound' if ok else 'UNSOUND'] += 1
        if not ok:
            bad.append((expr, pattern, old))
    print(caller, 'distinct differing bindings', len(rows), dict(counts))
    for row in bad[:20]:
        print('UNSOUND', *row, sep='\n  ')


if __name__ == '__main__':
    main(*sys.argv[1:])
