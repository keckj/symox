"""`run_suite.py <repo> <outdir> [shards]`: runs the whole tests/ suite in
disjoint file shards (default 24) and sums the junit results.

Every file matched by pytest.ini (test_*.py, compare_printing.py) under
tests/ goes to exactly one shard; tests/extra stays in its own shard.
Memory is limited to 6 GB per shard (prlimit --as), except for tests/extra:
jax and the CUDA runtimes reserve more address space than that and abort
under the cap."""
import collections
import concurrent.futures
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

root, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
nshards = int(sys.argv[3]) if len(sys.argv) > 3 else 24
out.mkdir(parents=True, exist_ok=True)
files = sorted(str(p.relative_to(root)) for p in (root / 'tests').rglob('*.py')
               if (p.name.startswith('test_') or p.name == 'compare_printing.py') and '__pycache__' not in p.parts)
extra = [f for f in files if f.startswith('tests/extra/')]
rest = [f for f in files if f not in extra]
shards = [[] for _ in range(nshards - 1)]
loads = [0] * (nshards - 1)
for f in sorted(rest, key=lambda f: -(root / f).stat().st_size):
    i = min(range(nshards - 1), key=loads.__getitem__)
    shards[i].append(f)
    loads[i] += (root / f).stat().st_size
shards.append(extra)
assert sorted(f for s in shards for f in s) == files
(out / 'shards.json').write_text(json.dumps(shards, indent=1))
env = dict(os.environ, PYTHONPATH=str(root / 'python'), OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           MKL_NUM_THREADS='1', MPLBACKEND='Agg')


def run(i, shard):
    junit = out / ('shard-%02d.xml' % i)
    cap = [] if shard is extra else ['prlimit', '--as=6000000000']
    cmd = cap + [sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider',
           '--tb=short', '--junitxml=' + str(junit)] + shard
    t0 = time.time()
    with open(out / ('shard-%02d.log' % i), 'w') as fh:
        rc = subprocess.run(cmd, cwd=root, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
    print('shard', i, 'rc', rc, round(time.time() - t0, 1), 's', flush=True)
    return i, rc, junit


counts = collections.Counter()
failures = []
t0 = time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=nshards) as pool:
    results = list(pool.map(lambda a: run(*a), [(i, s) for i, s in enumerate(shards) if s]))
for i, rc, junit in results:
    if not junit.exists():
        failures.append('shard %d: no junit (rc %d)' % (i, rc))
        continue
    for case in ET.parse(junit).getroot().iter('testcase'):
        name = case.get('classname', '') + '::' + case.get('name', '')
        kinds = [c.tag for c in case]
        if 'failure' in kinds or 'error' in kinds:
            counts['failed'] += 1
            failures.append(name)
        elif 'skipped' in kinds:
            counts['skipped'] += 1
        else:
            counts['passed'] += 1
summary = {'files': len(files), 'shards': sum(1 for s in shards if s), 'counts': dict(counts),
           'failures': failures, 'elapsed_seconds': round(time.time() - t0, 1)}
(out / 'summary.json').write_text(json.dumps(summary, indent=1))
print(json.dumps({k: v for k, v in summary.items() if k != 'failures'}))
print('\n'.join(failures[:200]))
