"""Print summary.md (medians over runs) from raw_sympy/raw_before/raw_after.tsv."""
import os
import statistics
import sys

here = os.path.dirname(os.path.abspath(__file__))


def load(name):
    rows = {}
    order = []
    with open(os.path.join(here, 'raw_%s.tsv' % name)) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            case = parts[1]
            if case not in rows:
                order.append(case)
            rows.setdefault(case, []).append(parts[2:])
    return rows, order


def med(rows, i):
    vals = [float(r[i]) for r in rows if r[0] != 'TIMEOUT']
    return statistics.median(vals) if vals else float('nan')


def digest(rows):
    tags = {r[2] for r in rows if len(r) > 2}
    return tags.pop() if len(tags) == 1 else 'MIXED %s' % sorted(tags)


sympy, order = load('sympy')
before, _ = load('before')
after, _ = load('after')
out = ['| call | SymPy 1.14 first | before first | after first | after repeat | after/SymPy | after/before | digest = SymPy |',
       '|---|---|---|---|---|---|---|---|']
for case in order:
    s, b, a = med(sympy[case], 0), med(before[case], 0), med(after[case], 0)
    same = digest(after[case]) == digest(sympy[case]) == digest(before[case])
    out.append('| `%s` | %.4f s | %.4f s | %.4f s | %.4f s | %.3g | %.3g | %s |' % (
        case, s, b, a, med(after[case], 1), a / s if s else float('nan'),
        a / b if b else float('nan'), 'yes' if same else 'NO'))
text = '\n'.join(out) + '\n'
open(os.path.join(here, 'summary.md'), 'w').write(text)
sys.stdout.write(text)
