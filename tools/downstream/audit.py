#!/usr/bin/env python3
"""Run pinned upstream tests first with SymPy 1.14, then strict native imports.

All downloads, environments and upstream working trees stay in --cache.
No upstream assertions are patched, and a dirty reference suite is never
reported as a Symox regression or a qualified integration.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix=path.name+'.', delete=False) as stream:
        stream.write(json.dumps(value, indent=2, sort_keys=True)+'\n')
        pending = Path(stream.name)
    pending.replace(path)


def prepare_sources(projects, cache):
    """Fetch immutable published archives or pinned Git revisions, not HEAD."""
    inventory_path = cache/'source-inventory.json'
    inventory = read(inventory_path) if inventory_path.exists() else {}
    for project in projects:
        name = project.get('distribution')
        if not name or project.get('reference_incompatible'):
            continue
        current = inventory.get(name)
        if isinstance(current, dict) and Path(current.get('root', '/nonexistent')).is_dir():
            continue
        git = project.get('source_git')
        if git:
            target = cache/'git'/project['id']
            target.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(['git', 'init', str(target)], check=True, capture_output=True)
            subprocess.run(['git', '-C', str(target), 'fetch', '--depth', '1', git['url'], git['commit']], check=True)
            subprocess.run(['git', '-C', str(target), 'checkout', '--detach', 'FETCH_HEAD'], check=True)
            inventory[name] = dict(root=str(target), commit=git['commit'])
        else:
            artifact = project.get('source_artifact')
            if not artifact:
                continue
            archive = cache/'archives'/artifact['filename']
            archive.parent.mkdir(parents=True, exist_ok=True)
            if not archive.exists():
                archive.write_bytes(urllib.request.urlopen(artifact['url'], timeout=120).read())
            if hashlib.sha256(archive.read_bytes()).hexdigest() != artifact['digests']['sha256']:
                raise ValueError('Published source digest mismatch: '+name)
            destination = cache/'sources'/name
            destination.mkdir(parents=True, exist_ok=True)
            with tarfile.open(archive) as package:
                for member in package.getmembers():
                    target = (destination/member.name).resolve()
                    if not target.is_relative_to(destination.resolve()):
                        raise ValueError('Unsafe archive path: '+member.name)
                    if member.issym() or member.islnk():
                        link = (target.parent/member.linkname if member.issym() else destination/member.linkname).resolve()
                        if not link.is_relative_to(destination.resolve()):
                            raise ValueError('Unsafe archive link: '+member.name)
                    if not (member.isfile() or member.isdir() or member.issym() or member.islnk()):
                        raise ValueError('Unsupported archive member: '+member.name)
                    package.extract(member, destination)
            inventory[name] = dict(root=str(next(destination.iterdir())), sha256=artifact['digests']['sha256'])
    write(inventory_path, inventory)


def command(argv, cwd, env, log, timeout):
    start = time.monotonic()
    log.parent.mkdir(parents=True, exist_ok=True)
    timed_out = False
    with log.open('w') as stream:
        process = subprocess.Popen(list(map(str, argv)), cwd=cwd, env=env,
            stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            code = process.returncode
    return dict(command=list(map(str, argv)), exit_code=code, timed_out=timed_out,
                seconds=round(time.monotonic()-start, 3), log=str(log),
                tail=log.read_text(errors='replace')[-16000:])


def phase_summary(report):
    counts = dict(passed=0, failed=0, errors=0, skipped=0, xfailed=0, xpassed=0)
    for item in report.get('reports', []):
        when, outcome = item['phase'], item['outcome']
        if item.get('wasxfail') is not None:
            if outcome == 'skipped': counts['xfailed'] += 1
            elif when == 'call' and outcome == 'passed': counts['xpassed'] += 1
        elif when == 'call':
            counts[outcome] += 1
        elif outcome == 'failed':
            counts['errors'] += 1
        elif outcome == 'skipped':
            counts['skipped'] += 1
    counts['skipped'] += len(report.get('collection_skips', []))
    counts['collection_errors'] = len(report.get('collection_errors', []))
    return counts


def native_status(execution):
    proof = execution.get('provenance', {})
    if proof.get('blocked_reference_execution'):
        return 'reference_execution_blocked'
    if execution.get('timed_out'):
        return 'native_timeout'
    if execution['exit_code']:
        return 'native_failure'
    if not proof.get('routed_modules'):
        return 'no_sympy_code_exercised'
    if execution.get('reference_passes_not_reproduced', {}).get('count'):
        return 'native_does_not_reproduce_reference_passes'
    if execution.get('counts', {}).get('passed', 0) < 1:
        return 'native_no_tests_passed'
    return 'native_suite_passes'


def archive_evidence(args, name, paths):
    directory = args.out/'evidence'/name
    directory.mkdir(parents=True, exist_ok=True)
    evidence = {}
    for path in paths:
        if not path.is_file(): continue
        target = directory/(path.name+'.gz')
        target.write_bytes(gzip.compress(path.read_bytes(), mtime=0))
        evidence[str(target.relative_to(args.out))] = hashlib.sha256(target.read_bytes()).hexdigest()
    return evidence


def run_project(project, args):
    name = project['id']
    destination = args.out/'results'/f'{name}.json'
    result = dict(project=project, started_at=datetime.now(timezone.utc).isoformat(),
                  reference_version='1.14.0', status='preparing',
                  candidate=dict(wheel=str(args.wheel), sha256=hashlib.sha256(args.wheel.read_bytes()).hexdigest()),
                  harness_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in SUPPORT.glob('*.py')})
    def save(status):
        result['status'] = status
        result['evidence'] = archive_evidence(args, name, (args.cache/'logs'/name).glob('*'))
        write(destination, result)
        print(name+': '+status, flush=True)
    if project.get('kind') in ('native', 'torch', 'publication'):
        if destination.exists() and project.get('kind') != 'publication':
            previous = read(destination)
            if (previous.get('candidate', {}).get('sha256') == result['candidate']['sha256']
                    and previous.get('project') == project):
                return previous
            requested = result['candidate']['sha256']
            result = previous
            result['requested_candidate_sha256'] = requested
            save('external_adapter_needs_candidate_refresh')
            return result
        save('special_adapter_required' if project['kind'] != 'publication' else 'publication_not_a_testable_distribution')
        return result
    if project.get('reference_incompatible'):
        save('excluded_reference_requirement')
        return result
    inventory = read(args.cache/'source-inventory.json')
    source_info = inventory.get(project['distribution'], {})
    if not isinstance(source_info, dict) or not source_info.get('root'):
        save('source_unavailable')
        return result
    source = Path(source_info['root'])
    result['source'] = str(source)
    result['source_identity'] = source_info
    tests = [path for path in project['tests'] if (source/path).exists()]
    if not tests:
        result['test_candidates'] = source_info
        save('upstream_tests_missing_from_source_archive')
        return result
    environment = args.cache/'envs'/name
    python = environment/'bin/python'
    env = os.environ.copy()
    env.update(MPLBACKEND='Agg', QT_QPA_PLATFORM='offscreen', OMP_NUM_THREADS='2',
               OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2', CMAKE_BUILD_PARALLEL_LEVEL='2',
               PYTHONHASHSEED='0', DEVITO_LANGUAGE='openmp', DEVITO_ARCH='gcc')
    env.update(project.get('environment', {}))
    env['SPYDER_CONFDIR'] = str(args.cache/'configs'/name)
    env['PATH'] = str(environment/'bin')+os.pathsep+env.get('PATH', '')
    # Keep other projects and the development environment out of resolution.
    env.pop('PYTHONPATH', None)
    constraints = args.cache/'constraints'/f'{name}.txt'
    constraints.parent.mkdir(parents=True, exist_ok=True)
    constraints.write_text('sympy==1.14.0\n'+(project.get('numpy_constraint', 'numpy<2') or '')+'\n')
    if not python.exists():
        result['venv'] = command(['uv', 'venv', '--python', project.get('python', '3.12'), environment],
                                ROOT, env, args.cache/'logs'/name/'venv.log', 120)
        if result['venv']['exit_code']:
            save('environment_creation_failed')
            return result
    dependencies = ['sympy==1.14.0', 'pytest>=8,<9', 'pytest-cov', 'pytest-timeout', 'pytest-xdist',
                    'pytest-mock', 'hypothesis', 'packaging', 'setuptools', 'wheel', 'mpmath', 'numpy',
                    *project.get('extra_dependencies', [])]
    package = [project['distribution']+'=='+project['version']] if project.get('install_mode') == 'wheel' else ['-e', source]
    result['install'] = command(['uv', 'pip', 'install', '--python', python, '--constraint', constraints,
                                  '--reinstall-package', 'symox', *dependencies, *package, args.wheel],
                                 ROOT, env, args.cache/'logs'/name/'install.log', args.install_timeout)
    if result['install']['exit_code']:
        save('installation_failed')
        return result
    site = next((environment/'lib').glob('python*/site-packages'))
    (site/'symox-downstream-audit.pth').write_text(str(SUPPORT)+'\n')
    preflight = 'import json,sympy,symox; from importlib.metadata import version; assert sympy.__version__=="1.14.0"; print(json.dumps({"sympy":sympy.__version__,"symox":symox.__version__,"numpy":version("numpy")}))'
    result['preflight'] = command([python, '-c', preflight], source, env,
                                  args.cache/'logs'/name/'preflight.log', 30)
    if result['preflight']['exit_code']:
        save('reference_preflight_failed')
        return result
    result['freeze'] = command(['uv', 'pip', 'freeze', '--python', python], ROOT, env,
                               args.cache/'logs'/name/'freeze.txt', 30)
    save('running_reference')
    test_cwd = args.cache if project.get('install_mode') == 'wheel' else source
    if project.get('install_mode') == 'wheel':
        tests = [str((site if project.get('installed_test_package') else source)/path) for path in tests]
    selections = [Path(t) if Path(t).is_absolute() else source/t for t in tests]
    selected_files = [p for selected in selections for p in
                      ([selected] if selected.is_file() else sorted(selected.rglob('*.py')))]
    result['selected_test_sha256'] = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in selected_files}
    if project.get('installed_test_package'):
        differences = [str(p) for p in selected_files if not (source/p.relative_to(site)).is_file()
                       or p.read_bytes() != (source/p.relative_to(site)).read_bytes()]
        result['installed_tests_match_pinned_source'] = not differences
        if differences:
            result['modified_test_files'] = differences
            save('installed_tests_differ_from_pinned_source')
            return result
    argv = [str(python), '-m', 'pytest', '-q', '-ra', '-p', 'audit_plugin', '--timeout='+str(project.get('per_test_timeout', 60)),
            *(['--import-mode=importlib'] if project.get('install_mode') == 'wheel' else []),
            '--junitxml='+str(args.cache/'logs'/name/'reference-junit.xml'),
            *tests, *project.get('pytest_args', [])]
    reference_passed = set()
    for backend in ('sympy', 'symox'):
        stage = 'reference' if backend == 'sympy' else 'native'
        report_path = args.cache/'logs'/name/f'{stage}-tests.json'
        proof_path = args.cache/'logs'/name/f'{stage}-provenance.json'
        report_path.unlink(missing_ok=True)
        proof_path.unlink(missing_ok=True)
        for old in proof_path.parent.glob(proof_path.stem+'.*.json'):
            old.unlink()
        run_env = dict(env, SYMOX_DOWNSTREAM_BACKEND=backend, PYTHONPATH=str(SUPPORT),
                       SYMOX_DOWNSTREAM_TEST_REPORT=str(report_path),
                       SYMOX_DOWNSTREAM_PROVENANCE=str(proof_path))
        stage_argv = [word.replace('reference-junit.xml',stage+'-junit.xml') for word in argv]
        execution = command(stage_argv, test_cwd, run_env, args.cache/'logs'/name/f'{stage}.log', args.timeout)
        if report_path.exists():
            report = read(report_path)
            execution['counts'] = phase_summary(report)
            evidence = args.out/'evidence'/name/f'{stage}-tests.json.gz'
            evidence.parent.mkdir(parents=True, exist_ok=True)
            evidence.write_bytes(gzip.compress(report_path.read_bytes(), mtime=0))
            execution['test_evidence'] = str(evidence.relative_to(args.out))
            execution['test_evidence_sha256'] = hashlib.sha256(evidence.read_bytes()).hexdigest()
            execution['first_failures'] = [item for item in report.get('reports', []) if item['outcome']=='failed'][:3]
            execution['collection_errors'] = report.get('collection_errors', [])[:4]
            execution['collection_skips'] = report.get('collection_skips', [])
            passed = {row['nodeid'] for row in report.get('reports', [])
                      if row['phase'] == 'call' and row['outcome'] == 'passed'}
            if backend == 'sympy':
                reference_passed = passed
            else:
                missing = sorted(reference_passed-passed)
                execution['reference_passes_not_reproduced'] = dict(count=len(missing), examples=missing[:30])
        proofs = [read(p) for p in proof_path.parent.glob(proof_path.stem+'.*.json')]
        if proof_path.exists(): proofs.append(read(proof_path))
        if proofs:
            execution['process_provenance'] = proofs
            execution['provenance'] = dict(
                routed_modules={k:v for proof in proofs for k,v in proof['routed_modules'].items()},
                blocked_reference_execution=[value for proof in proofs for value in proof['blocked_reference_execution']])
        result[stage] = execution
        counts = execution.get('counts', {})
        clean = execution['exit_code'] == 0 and counts.get('passed', 0) > 0
        if backend == 'sympy':
            result['reference_qualified'] = clean
            if not clean:
                save('reference_not_clean' if not execution['timed_out'] else 'reference_timeout')
                return result
            save('running_native')
        else:
            save(native_status(execution))
    return result


def render(args):
    rows = [read(path) for path in sorted((args.out/'results').glob('*.json'))]
    manifest = read(SUPPORT/'projects.json')
    projects = {p['id']:p for p in manifest['projects']}
    acceptable = {'native_suite_passes', 'publication_not_a_testable_distribution', 'excluded_reference_requirement'}
    blockers = [r['project']['id'] for r in rows if r['status'] not in acceptable]
    statuses = {}
    for row in rows: statuses[row['status']] = statuses.get(row['status'], 0)+1
    readiness = dict(reference_version='1.14.0', universal_parity_certified=False,
        featured_drop_in_ready=not blockers and len(rows) == len(projects),
        statuses=statuses, blockers=blockers)
    lines = ['# Featured-project integration audit', '',
        'Reference: **SymPy 1.14.0**. **'+('The candidate passes the selected integration gate.' if readiness['featured_drop_in_ready'] else 'The candidate does not pass the featured-project replacement gate.')+'**', '',
        'These results measure the stated upstream test selections, not universal compatibility or mathematical parity.',
        'Each native run requires a clean reference run on the same unmodified upstream selection and environment.',
        'A passing reference case that becomes skipped, expected to fail, or absent prevents native success.', '',
        'The experimental process-local importer routes `sympy.*` to native `symox.*` and blocks execution of installed reference SymPy Python files.',
        'Its facade advertises the **target namespace version** 1.14.0 for dependency version gates; the actual Symox version and wheel digest are recorded separately.',
        'This audit shim is not a supported production replacement installer.', '',
        'Sources come from the [official featured list](https://www.sympy.org/en/index.html). The yt unit-system component, yt unit tests and Devito symbolic tests have separate rows.',
        'Commands, source revisions, wheel digests, dependency freezes, test reports and process provenance are in [the audit artifacts](../audits/downstream/summary.json).',
        'Recorded assertions are checked against [the pinned upstream sources](../audits/downstream/source-integrity.json).',
        'Local regression checks, their exact scope and the production fix commit are recorded in [the validation evidence](../audits/downstream/validation.json).', '',
        '| Project | Version | SymPy 1.14 reference | Native Symox | Outcome |', '|---|---|---|---|---|']
    labels = dict(native_suite_passes='Selected tests pass', native_failure='Native incompatibility',
        reference_not_clean='Reference fails; native not run', reference_timeout='Reference timed out; native not run',
        no_sympy_code_exercised='No symbolic integration exercised', excluded_reference_requirement='Requires an older SymPy',
        publication_not_a_testable_distribution='Publication, not a package',
        special_adapter_required='Adapter not yet run', running_reference='Reference running', running_native='Native running')
    def count(row, stage):
        data = row.get(stage, {}).get('counts')
        if not data:
            execution = row.get(stage)
            if execution and execution.get('exit_code'):
                return 'startup failure (exit '+str(execution['exit_code'])+')'
            return 'not run/completed'
        parts = [str(data['passed'])+' pass']
        for key, label in [('failed','fail'),('errors','error'),('collection_errors','collection error'),
                           ('skipped','skip'),('xfailed','xfail'),('xpassed','xpass')]:
            if data.get(key): parts.append(str(data[key])+' '+label)
        return ', '.join(parts)
    for row in rows:
        project = projects.get(row['project']['id'], row['project'])
        link = '../audits/downstream/results/'+project['id']+'.json'
        lines.append('| ['+project['name']+']('+link+') | '+project.get('version','—')+' | '+count(row,'reference')+' | '+count(row,'native')+' | '+labels.get(row['status'], row['status'])+' |')
    lines += ['', '## Scope and first blockers', '']
    for row in rows:
        project = projects.get(row['project']['id'], row['project'])
        lines += ['### '+project['name'], '', project.get('scope', 'Upstream test tree selected in the pinned manifest.'), '']
        if project.get('reference_incompatible'):
            lines += ['Excluded before native testing: '+project['reference_incompatible']+'.', '']
        execution = row.get('native',row.get('reference',row.get('install',{})))
        failures = execution.get('collection_errors') or execution.get('first_failures')
        if failures:
            lines += ['```text', str(failures[0].get('detail',''))[-1700:], '```', '']
        elif execution.get('exit_code'):
            lines += ['```text',execution.get('tail','')[-1700:], '```', '']
        if row.get('reference_confirmation'):
            confirm = row['reference_confirmation'].get('counts', {})
            lines += ['The four MPI cases that failed in the full parallel reference run were rerun without pytest workers: '+str(confirm.get('passed', 0))+' passed. The full-suite reference remains unqualified; this is not evidence of a Symox regression.', '']
        if row['status'] == 'native_suite_passes':
            lines += ['All passing reference node IDs were reproduced. The skip/xfail counts above remain outside the proven scope.', '']
    lines += ['## Reproduction', '',
        'See [the harness instructions](../tools/downstream/README.md) for runtime setup, dependency pins and native adapters.', '',
        '```bash', 'python tools/downstream/audit.py --wheel /path/to/symox.whl --cache /tmp/symox-featured --out audits/downstream',
        '# A CI gate for a specific integration:',
        'python tools/downstream/audit.py --wheel /path/to/symox.whl --projects pyneqsys --require-native-pass',
        '```', '',
        'The ordinary audit exits successfully when evidence collection finishes, including negative results. `--require-native-pass` makes an unsuccessful integration fail the command.',
        'Native runtimes use `tools/downstream/external.py`. Publications and incompatible dependency requirements are never counted as tested native packages.', '']
    body = re.sub(r'\x1b\[[0-?]*[ -/]*[@-~]', '', '\n'.join(lines))
    (ROOT/'docs/downstream-integration.md').write_text('\n'.join(line.rstrip() for line in body.splitlines())+'\n')
    index = [dict(project=r['project']['id'], status=r['status'],
                  reference=r.get('reference', {}).get('counts'), native=r.get('native', {}).get('counts'),
                  result='results/'+r['project']['id']+'.json') for r in rows]
    write(args.out/'summary.json', dict(**readiness, results=index))
    write(args.out/'readiness.json', readiness)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=SUPPORT/'projects.json')
    parser.add_argument('--cache', type=Path, default=Path('/tmp/symox-featured'))
    parser.add_argument('--out', type=Path, default=ROOT/'audits/downstream')
    parser.add_argument('--wheel', type=Path)
    parser.add_argument('--projects', nargs='+')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--timeout', type=int, default=900)
    parser.add_argument('--install-timeout', type=int, default=900)
    parser.add_argument('--report-only', action='store_true')
    parser.add_argument('--require-native-pass', action='store_true',
                        help='Exit unsuccessfully unless every selected project reproduces its reference passes')
    args = parser.parse_args()
    args.cache=args.cache.resolve();args.out=args.out.resolve()
    if args.report_only:
        render(args);return
    if args.wheel is None:parser.error('--wheel is required')
    if not args.wheel.is_file():parser.error('--wheel must name an existing file')
    if args.jobs < 1:parser.error('--jobs must be positive')
    args.wheel=args.wheel.resolve()
    projects = read(args.manifest)['projects']
    if args.projects:
        unknown = set(args.projects)-{p['id'] for p in projects}
        if unknown: parser.error('Unknown projects: '+', '.join(sorted(unknown)))
        projects=[p for p in projects if p['id'] in args.projects]
    prepare_sources(projects, args.cache)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_project,p,args):p for p in projects}
        for future in as_completed(futures):
            try:future.result()
            except Exception as error:
                project=futures[future]
                write(args.out/'results'/(project['id']+'.json'),dict(project=project,status='harness_error',error=repr(error)))
                print(project['id']+': '+repr(error),flush=True)
    render(args)
    if args.require_native_pass:
        if any(read(args.out/'results'/(p['id']+'.json'))['status'] != 'native_suite_passes' for p in projects):
            raise SystemExit(1)


if __name__ == '__main__':
    main()
