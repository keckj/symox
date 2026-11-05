"""Sum over integrands of the minimum time over runs (prof_targets.py
outputs), excluding integrands that timed out in any run.
usage: sum_min.py RUN1.txt RUN2.txt ..."""
import sys
best, timeouts = {}, set()
for path in sys.argv[1:]:
    for line in open(path):
        if line.startswith('TIME '):
            _, t, e = line.rstrip('\n').split(' ', 2)
            best[e] = min(best.get(e, float('inf')), float(t))
        elif line.startswith('EXC ') and line.rstrip().endswith('TimeoutError'):
            timeouts.add(line.split(' ', 1)[1].rsplit(' ', 1)[0])
print('%d integrands, %d timed out (excluded): %.3f s' % (
    len(best), len(timeouts), sum(t for e, t in best.items() if e not in timeouts)))
