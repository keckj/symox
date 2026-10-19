#!/usr/bin/env python3
"""Audit real upstream tests in an existing native/Conda runtime.

Use --prefix for the installed Sage/Octave/Cadabra or PyTorch environment.
The Python wheel is force installed into that isolated prefix. No reference
failure permits a native run. Commands, raw logs and import provenance are
retained just as for the ordinary Python projects.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re

from audit import ROOT, SUPPORT, archive_evidence, command, phase_summary, read, render, write


def parse_counts(kind, execution, report_path):
    counts = dict(passed=0, failed=0, errors=0, collection_errors=0,
                  skipped=0, xfailed=0, xpassed=0)
    if kind == 'pytest' and report_path.exists():
        report = read(report_path)
        execution['first_failures'] = [r for r in report['reports'] if r['outcome'] == 'failed'][:3]
        execution['collection_errors'] = report.get('collection_errors', [])[:4]
        execution['passed_nodeids'] = sorted(r['nodeid'] for r in report['reports']
            if r['phase'] == 'call' and r['outcome'] == 'passed')
        return phase_summary(report)
    output = Path(execution['log']).read_text(errors='replace')
    if kind == 'sage':
        match = re.search(r'\[(\d+) tests(?:, (\d+) failures)?[,\]]', output)
        if match:
            counts.update(passed=int(match[1])-int(match[2] or 0), failed=int(match[2] or 0))
    elif kind == 'octave':
        match = re.search(r'^SYMOX_BIST_RESULTS=(.+)$', output, re.M)
        if match:
            execution['bist_files'] = json.loads(match[1])
            for row in execution['bist_files']:
                for key in ('passed', 'failed', 'skipped', 'xfailed'):
                    counts[key] += row[key]
    elif kind == 'generator' and execution['exit_code'] == 0:
        # One real script execution, not a claim about the C++ test suite.
        counts['passed'] = 1
    elif kind == 'cadabra':
        execution['upstream_pass_markers'] = len(re.findall(r'Test\s*\d+[a-z]?\s+passed', output, re.I))
        if execution['exit_code'] == 0: counts['passed'] = 1
    if execution['exit_code'] and not counts['failed']:
        counts['errors'] = 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project', choices=['sagemath', 'octsympy', 'cadabra', 'qmcpack', 'torchinductor'])
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--wheel', type=Path, required=True)
    parser.add_argument('--cache', type=Path, default=Path('/tmp/symox-featured'))
    parser.add_argument('--out', type=Path, default=ROOT/'audits/downstream')
    parser.add_argument('--micromamba', type=Path)
    parser.add_argument('--timeout', type=int, default=900)
    args = parser.parse_args()
    args.prefix = args.prefix.resolve(); args.wheel = args.wheel.resolve()
    args.cache = args.cache.resolve(); args.out = args.out.resolve()
    project = next(p for p in read(SUPPORT/'projects.json')['projects'] if p['id'] == args.project)
    name = args.project
    source = args.cache/'git'/('cadabra' if name == 'cadabra' else name)
    python = args.prefix/'bin/python'
    site = next((args.prefix/'lib').glob('python*/site-packages'))
    env = dict(os.environ, MPLBACKEND='Agg', QT_QPA_PLATFORM='offscreen',
               OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', PYTHONHASHSEED='0')
    env.pop('PYTHONPATH', None); env.pop('SYMOX_DOWNSTREAM_BACKEND', None)
    prefix = [str(args.micromamba), 'run', '-p', str(args.prefix)] if args.micromamba else []
    logs = args.cache/'logs'/name
    result = dict(project=project, started_at=datetime.now(timezone.utc).isoformat(),
        reference_version='1.14.0', candidate=dict(wheel=str(args.wheel),
        sha256=hashlib.sha256(args.wheel.read_bytes()).hexdigest()),
        harness_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in SUPPORT.glob('*.py')})
    destination = args.out/'results'/(name+'.json')
    def save(status):
        result['status'] = status
        result['evidence'] = archive_evidence(args, name, logs.glob('*'))
        write(destination, result)
        print(name+': '+status, flush=True)
    result['install'] = command(['uv', 'pip', 'install', '--python', python,
        '--reinstall-package', 'symox', args.wheel], ROOT, env, logs/'install.log', 120)
    if result['install']['exit_code']:
        save('installation_failed'); return
    (site/'symox-downstream-audit.pth').write_text(str(SUPPORT)+'\n')
    result['preflight'] = command(prefix+[str(python), '-c',
        'import sympy,symox; assert sympy.__version__=="1.14.0"; print(sympy.__version__,symox.__version__)'],
        args.cache, env, logs/'preflight.log', 60)
    if result['preflight']['exit_code']:
        save('reference_preflight_failed'); return
    result['freeze'] = command(['uv', 'pip', 'freeze', '--python', python], ROOT, env, logs/'freeze.txt', 30)
    if args.micromamba:
        result['conda_explicit'] = command([args.micromamba, 'list', '-p', args.prefix, '--explicit'],
                                          ROOT, env, logs/'conda-explicit.txt', 60)
    if name == 'sagemath':
        kind = 'sage'; selected = [site/'sage/interfaces/sympy.py']
        argv = [python, '-m', 'sage.doctest', '--random-seed=0', '--optional=sage,sympy', '--timeout=600', *selected]
    elif name == 'octsympy':
        kind = 'octave'; selected = [source/'inst/@sym'/(n+'.m') for n in ('sym','diff','int','solve')]
        env.update(SYMOX_OCTSYMPY_SOURCE=str(source), SYMOX_AUDIT_PYTHON=str(python))
        argv = [args.prefix/'bin/octave', '--no-gui', '--no-window-system', '--quiet', SUPPORT/'octave_bist.m']
    elif name == 'qmcpack':
        kind = 'generator'; selected = [source/'src/QMCWaveFunctions/tests/gen_matrix_ops.py']
        argv = [python, *selected]
    elif name == 'torchinductor':
        kind = 'pytest'; selected = [args.cache/'git/torch/test/test_sympy_utils.py']
        argv = [python, '-m', 'pytest', '-q', '-ra', '-p', 'audit_plugin', '--timeout=60', *selected]
    else:
        kind = 'cadabra'; selected = [source/'tests/sympy_cdb.cdb']
        argv = [args.prefix/'bin/cadabra2', *selected]
    result['test_kind'] = kind
    result['selected_tests'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in selected}
    reference_nodes = set()
    reference_output = None
    for stage, backend in [('reference', 'sympy'), ('native', 'symox')]:
        save('running_'+stage)
        report_path = logs/(stage+'-tests.json')
        proof_path = logs/(stage+'-provenance.json')
        report_path.unlink(missing_ok=True)
        for path in logs.glob(proof_path.stem+'*.json'): path.unlink()
        run_env = dict(env, SYMOX_DOWNSTREAM_BACKEND=backend, PYTHONPATH=str(SUPPORT),
            SYMOX_DOWNSTREAM_TEST_REPORT=str(report_path), SYMOX_DOWNSTREAM_PROVENANCE=str(proof_path))
        execution = command(prefix+list(map(str,argv)), args.cache, run_env, logs/(stage+'.log'), args.timeout)
        execution['counts'] = parse_counts(kind, execution, report_path)
        if kind == 'pytest':
            nodes = set(execution.pop('passed_nodeids', []))
            if stage == 'reference': reference_nodes = nodes
            else: execution['reference_passes_not_reproduced'] = sorted(reference_nodes-nodes)
        elif kind == 'generator':
            output = Path(execution['log']).read_bytes()
            execution['output_sha256'] = hashlib.sha256(output).hexdigest()
            if stage == 'reference': reference_output = output
            else: execution['reference_output_reproduced'] = output == reference_output
        proofs = [read(p) for p in logs.glob(proof_path.stem+'*.json')]
        execution['process_provenance'] = proofs
        result[stage] = execution
        artifacts = args.out/'evidence'/name; artifacts.mkdir(parents=True, exist_ok=True)
        for path in [Path(execution['log']), report_path]:
            if path.exists(): (artifacts/(path.name+'.gz')).write_bytes(gzip.compress(path.read_bytes(), mtime=0))
        clean = execution['exit_code'] == 0 and execution['counts']['passed'] > 0
        if stage == 'reference':
            result['reference_qualified'] = clean
            if not clean:
                save('reference_timeout' if execution['timed_out'] else 'reference_not_clean'); break
        else:
            active = any(p['routed_modules'] for p in proofs)
            blocked = any(p['blocked_reference_execution'] for p in proofs)
            omissions = execution.get('reference_passes_not_reproduced')
            same = execution.get('reference_output_reproduced', True)
            if blocked: save('reference_execution_blocked')
            elif not active: save('no_sympy_code_exercised')
            elif clean and same and not omissions: save('native_suite_passes')
            else: save('native_timeout' if execution['timed_out'] else 'native_failure')
    render(args)


if __name__ == '__main__':
    main()
