"""Full tests/ suite in 24 disjoint file shards with per-file coverage
proof (adapted from audits/sympy-1.14.0/native-sdm-irref-2026-10-28).

usage: run_suite.py <worktree> <evidence dir>
Prints "SUITE CLEAN" only when every collected test of every file ran and
none failed or errored.
"""
import concurrent.futures
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

root = Path(sys.argv[1]).resolve()
evidence = Path(sys.argv[2]).resolve()
evidence.mkdir(parents=True, exist_ok=True)
here = Path(__file__).resolve().parent
SHARDS = 24
started = time.time()
env = dict(os.environ, PYTHONPATH=str(root/'python') + os.pathsep + str(here),
           SYMOX_CI_PY=str(root/'python'),
           PYTEST_DISABLE_PLUGIN_AUTOLOAD='1', OPENBLAS_NUM_THREADS='1',
           OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', MPLBACKEND='Agg')
so = root/'python/symox/_core.abi3.so'
identity = {'so_sha256': hashlib.sha256(so.read_bytes()).hexdigest(),
            'head': subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=root, capture_output=True, text=True).stdout.strip(),
            'diff_sha256': hashlib.sha256(subprocess.run(['git', 'diff', 'HEAD'], cwd=root, capture_output=True).stdout).hexdigest()}
(evidence/'identity.json').write_text(json.dumps(identity, indent=2))


def run(name, files, collect=False):
    folder = evidence/name
    folder.mkdir(exist_ok=True)
    local = dict(env, SYMOX_VALIDATION_RECORD=str(folder/'execution.json'),
                 MPLCONFIGDIR=str(folder/'matplotlib'))
    # The collection and the tests/extra shard import torch, jax and cuda,
    # which reserve more address space than --as allows (MemoryError inside
    # pytest), so they are capped by a cgroup on real memory instead; the
    # other shards keep the address-space limit.
    cgroup = collect or any(f.startswith('tests/extra/') for f in files)
    cap = (['systemd-run', '--user', '--scope', '-q', '-p', 'MemoryMax=6G'] if cgroup
           else ['prlimit', '--as=6000000000'])
    command = cap + [sys.executable, '-m', 'pytest', '-p', 'validation_plugin', '-q',
                     '--tb=short', '-o', 'cache_dir='+str(folder/'pytest-cache')]
    if collect:
        command.append('--collect-only')
    else:
        command.extend(['--junitxml='+str(folder/'junit.xml'), '--durations=8'])
    command.extend(files)
    (folder/'command.json').write_text(json.dumps(command, indent=2))
    start = time.time()
    with (folder/'pytest.log').open('w') as output:
        process = subprocess.run(command, cwd=root, env=local, stdout=output, stderr=subprocess.STDOUT)
    result = {'name': name, 'returncode': process.returncode,
              'elapsed_seconds': time.time()-start, 'files': files}
    (folder/'result.json').write_text(json.dumps(result, indent=2))
    print(name, process.returncode, round(result['elapsed_seconds'], 2), flush=True)
    return result


collection_result = run('collection', ['tests'], True)
if collection_result['returncode'] != 0:
    raise SystemExit('Collection failed; see collection/pytest.log')
collection = json.loads((evidence/'collection/execution.json').read_text())
files = collection['files']
weights = Counter(node.split('::', 1)[0] for node in collection['collected'])
shards = [[] for _ in range(SHARDS)]
loads = [0]*(SHARDS-1)
for path in sorted(files, key=lambda f: (-weights[f], f)):
    if path.startswith('tests/extra/'):
        shards[SHARDS-1].append(path)
    else:
        index = min(range(SHARDS-1), key=loads.__getitem__)
        shards[index].append(path)
        loads[index] += max(weights[path], 1)
assert Counter(p for shard in shards for p in shard) == Counter(files)
(evidence/'shards.json').write_text(json.dumps(shards, indent=2))
print('Collected', len(collection['collected']), 'tests in', len(files), 'files', flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=SHARDS) as pool:
    futures = [pool.submit(run, 'shard-%02d' % i, shard) for i, shard in enumerate(shards) if shard]
    results = [future.result() for future in concurrent.futures.as_completed(futures)]
collected, finished = [], []
counts = Counter()
failures = []
skips = Counter()
per_file = Counter()
for result in results:
    folder = evidence/result['name']
    report = json.loads((folder/'execution.json').read_text())
    collected.extend(report['collected'])
    finished.extend(report['finished'])
    xml = ET.parse(folder/'junit.xml').getroot()
    for case in xml.iter('testcase'):
        if case.find('failure') is not None or case.find('error') is not None:
            counts['failed'] += 1
            failures.append({'shard': result['name'], 'class': case.get('classname'), 'name': case.get('name')})
        elif case.find('skipped') is not None:
            child = case.find('skipped')
            key = 'xfailed' if child.get('type') == 'pytest.xfail' else 'skipped'
            counts[key] += 1
            if key == 'skipped':
                skips[child.get('message', '')] += 1
        else:
            counts['passed'] += 1
expected = Counter(collection['collected'])
fin = Counter(finished)
files_incomplete = sorted({n.split('::', 1)[0] for n in (expected - fin).elements()})
coverage = {'collection_matches': Counter(collected) == expected,
            'execution_matches': fin == expected,
            'expected': len(collection['collected']), 'collected_in_shards': len(collected),
            'finished': len(finished), 'files': len(files), 'files_incomplete': files_incomplete,
            'missing': list((expected-fin).elements())[:50],
            'unexpected': list((fin-expected).elements())[:50]}
summary = {'identity': identity, 'counts': dict(counts), 'coverage': coverage,
           'shards': sorted(results, key=lambda r: r['name']),
           'failures': failures, 'skip_reasons': dict(skips),
           'elapsed_seconds': time.time()-started}
(evidence/'summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps({key: summary[key] for key in ('counts', 'coverage', 'failures', 'elapsed_seconds')}, indent=2), flush=True)
success = (not failures and all(r['returncode'] == 0 for r in results)
           and coverage['collection_matches'] and coverage['execution_matches'])
print('SUITE CLEAN' if success else 'SUITE NOT CLEAN', flush=True)
raise SystemExit(0 if success else 1)
