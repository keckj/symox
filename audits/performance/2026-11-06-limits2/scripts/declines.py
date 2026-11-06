"""Group the Gruntz declines of a trace run: declines.py TAG [RAW_DIR].

Reads raw/tl/TAG_*.trace.json (statement harness) and raw/pt/TAG_*.json
(pytest run); prints the count per reason (the reason's prefix up to ':'),
and the total of distinct (reason, limit) pairs per source."""
import collections, glob, json, os, sys

tag = sys.argv[1]
raw = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(__file__), '..', 'raw')
for src in ('tl', 'pt'):
    pat = os.path.join(raw, src, '%s_*.json' % tag)
    files = [f for f in glob.glob(pat) if (src == 'pt') != f.endswith('.trace.json') or src == 'tl' and f.endswith('.trace.json')]
    files = [f for f in files if (f.endswith('.trace.json') if src == 'tl' else not f.endswith('.trace.json'))]
    total = 0
    groups = collections.Counter()
    per_file = {}
    for f in sorted(files):
        d = json.load(open(f))
        per_file[os.path.basename(f)] = len(d)
        total += len(d)
        for reason, _ in d:
            groups[reason.split(':')[0][:70]] += 1
    print('== %s %s: %d declines' % (src, tag, total), per_file)
    for k, v in groups.most_common():
        print('  %4d  %s' % (v, k))
