"""Run the full tests/ suite in N disjoint file shards (default 24).

usage: run_shards.py WORKTREE OUTDIR [N]

Collects tests/ once, assigns each test file to the least loaded shard by
collected test count (largest files first), runs one pytest process per
shard with PYTHONPATH=WORKTREE/python (capped by prlimit --as=$SHARD_AS when
set) and writes OUTDIR/shardNN.txt plus
OUTDIR/summary.json (per-shard summary line and the failing node ids).
"""
import collections
import concurrent.futures
import json
import os
import re
import subprocess
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
n = int(sys.argv[3]) if len(sys.argv) > 3 else 24
out.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, PYTHONPATH=str(root / 'python'), PYTEST_DISABLE_PLUGIN_AUTOLOAD='1',
           OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', MPLBACKEND='Agg')
col = subprocess.run([sys.executable, '-m', 'pytest', '--collect-only', '-q', '-p', 'no:cacheprovider', 'tests'],
                     cwd=root, env=env, capture_output=True, text=True)
nodes = [line for line in col.stdout.splitlines() if '::' in line]
weights = collections.Counter(node.split('::', 1)[0] for node in nodes)
files = sorted(weights, key=lambda f: (-weights[f], f))
shards = [[] for _ in range(n)]
loads = [0] * n
for f in files:
    i = min(range(n), key=loads.__getitem__)
    shards[i].append(f)
    loads[i] += weights[f]
(out / 'shards.json').write_text(json.dumps(shards, indent=1))
print('collected', len(nodes), 'tests in', len(files), 'files', flush=True)


# SHARD_AS: address-space cap of each shard process (prlimit --as); the
# collection of the whole suite in one process needs more than 6 GB and runs
# uncapped.
cap = ['prlimit', '--as=' + os.environ['SHARD_AS']] if os.environ.get('SHARD_AS') else []


def run(i):
    log = out / ('shard%02d.txt' % i)
    with log.open('w') as fh:
        p = subprocess.run(cap + [sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider', '--tb=short',
                            '-rfE'] + shards[i], cwd=root, env=env, stdout=fh, stderr=subprocess.STDOUT)
    text = log.read_text()
    tail = [line for line in text.splitlines() if re.search(r'\d+ (passed|failed)', line)]
    failed = re.findall(r'^(?:FAILED|ERROR) (\S+)', text, re.M)
    return {'shard': i, 'returncode': p.returncode, 'summary': tail[-1] if tail else '', 'failed': failed}


with concurrent.futures.ThreadPoolExecutor(n) as ex:
    results = list(ex.map(run, range(n)))
tot = collections.Counter()
for r in results:
    for num, kind in re.findall(r'(\d+) (passed|failed|skipped|xfailed|xpassed|errors?)', r['summary']):
        tot[kind] += int(num)
(out / 'summary.json').write_text(json.dumps({'collected': len(nodes), 'totals': tot, 'shards': results}, indent=1))
print(json.dumps(tot))
for r in results:
    for f in r['failed']:
        print('FAILED', f)
