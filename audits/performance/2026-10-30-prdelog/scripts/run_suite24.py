"""Full pytest collection and disjoint file shards, with execution coverage proof.

Copy of audits/sympy-1.14.0/native-sdm-irref-2026-10-28/run_suite.py with 24
shards, evidence in ~/symox-ci/validation/prdelog-20261030k (identity.json
instead of a source manifest), shards under prlimit --as=6000000000 and the
collection under systemd-run MemoryMax=6G.
usage: python run_suite24.py WORKTREE
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
evidence = Path('/home/keckj/symox-ci/validation/prdelog-20261030k')
evidence.mkdir(exist_ok=True)
started = time.time()
env = dict(os.environ, PYTHONPATH=str(root/'python') + os.pathsep + str(evidence),
           SYMOX_CI_PY=str(root/'python'),
           PYTEST_DISABLE_PLUGIN_AUTOLOAD='1', OPENBLAS_NUM_THREADS='1',
           OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', MPLBACKEND='Agg')
manifest={'identity':json.loads((evidence/'identity.json').read_text())}


def run(name, files, collect=False):
    folder = evidence/name
    folder.mkdir(exist_ok=True)
    local = dict(env, SYMOX_VALIDATION_RECORD=str(folder/'execution.json'),
                 MPLCONFIGDIR=str(folder/'matplotlib'))
    command = ['prlimit','--as=6000000000',sys.executable, '-m', 'pytest', '-p', 'validation_plugin', '-q',
               '--tb=short', '-o', 'cache_dir='+str(folder/'pytest-cache')]
    if collect:
        # The collection JSON exceeds a 6 GB address-space cap; cap RSS instead.
        command[:2] = ['systemd-run', '--user', '--scope', '-q', '-p', 'MemoryMax=6G']
        command.append('--collect-only')
    else:
        command.extend(['--junitxml='+str(folder/'junit.xml'), '--durations=8'])
    command.extend(files)
    (folder/'command.json').write_text(json.dumps(command, indent=2))
    start = time.time()
    with (folder/'pytest.log').open('w') as output:
        process = subprocess.run(command, cwd=root, env=local,
                                 stdout=output, stderr=subprocess.STDOUT)
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
shards = [[] for _ in range(24)]
loads = [0]*23
# Keep optional runtime tests together; every collected file remains included.
for path in sorted(files, key=lambda f: (-weights[f], f)):
    if path.startswith('tests/extra/'):
        shards[23].append(path)
    else:
        index = min(range(23), key=loads.__getitem__)
        shards[index].append(path)
        loads[index] += max(weights[path], 1)
assert Counter(p for shard in shards for p in shard) == Counter(files)
(evidence/'shards.json').write_text(json.dumps(shards, indent=2))
print('Collected', len(collection['collected']), 'tests in', len(files), 'files', flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=24) as pool:
    futures = [pool.submit(run, 'shard-%02d'%i, shard) for i, shard in enumerate(shards) if shard]
    results = [future.result() for future in concurrent.futures.as_completed(futures)]
collected, finished = [], []
counts = Counter()
failures = []
skips = Counter()
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
coverage = {'collection_matches': Counter(collected) == expected,
            'execution_matches': Counter(finished) == expected,
            'expected': len(collection['collected']), 'collected_in_shards': len(collected),
            'finished': len(finished), 'missing': list((expected-Counter(finished)).elements()),
            'unexpected': list((Counter(finished)-expected).elements())}
summary = {'source_identity': manifest['identity'], 'counts': dict(counts),
           'coverage': coverage, 'shards': sorted(results, key=lambda r:r['name']),
           'failures': failures, 'skip_reasons': dict(skips),
           'elapsed_seconds': time.time()-started}
(evidence/'summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps({key: summary[key] for key in ('counts', 'coverage', 'failures', 'elapsed_seconds')}, indent=2), flush=True)
success = (not failures and all(r['returncode'] == 0 for r in results)
           and coverage['collection_matches'] and coverage['execution_matches'])
raise SystemExit(0 if success else 1)
