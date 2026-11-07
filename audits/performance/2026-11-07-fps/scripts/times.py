"""times.py TAG...: total and median recorded time (s) of each raw replay.

Only the keys present in every tag are summed, so the totals compare the
same statements.
"""
import json, os, statistics, sys

raw = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'raw')
data = {t: json.load(open(os.path.join(raw, t + '.json'))) for t in sys.argv[1:]}
keys = set.intersection(*(set(d) for d in data.values()))
for t, d in data.items():
    times = [d[k][2] for k in keys]
    print('%-14s keys=%d total=%.3f median=%.4f max=%.3f' % (
        t, len(keys), sum(times), statistics.median(times), max(times)))
