"""Verify that downstream evidence cannot silently use the reference backend."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import hashlib
from types import SimpleNamespace

import pytest


SUPPORT = Path(__file__).resolve().parents[1]/'tools/downstream'
spec = importlib.util.spec_from_file_location('downstream_audit', SUPPORT/'audit.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def run_backend(tmp_path, backend, source):
    env = dict(os.environ, PYTHONPATH=str(SUPPORT), SYMOX_DOWNSTREAM_BACKEND=backend,
               SYMOX_DOWNSTREAM_PROVENANCE=str(tmp_path/'proof.json'))
    completed = subprocess.run([sys.executable, '-c', source], env=env, cwd=tmp_path,
                               capture_output=True, text=True, timeout=30)
    proofs = [json.loads(p.read_text()) for p in tmp_path.glob('proof.*.json')]
    return completed, proofs


def test_native_imports_keep_actual_symbol_identity_and_no_reference_files(tmp_path):
    process, proofs = run_backend(tmp_path, 'symox', '''
import sympy, symox
from sympy.core.symbol import Symbol
assert sympy.Symbol is symox.Symbol is Symbol
assert sympy.__version__ == '1.14.0'
assert symox.__version__ == sympy.__symox_version__
assert (Symbol('x')+1).subs(Symbol('x'), 2) == 3
''')
    assert process.returncode == 0, process.stderr
    assert proofs[0]['routed_modules']['sympy.core.symbol'] == 'symox.core.symbol'
    assert not proofs[0]['blocked_reference_execution']
    assert all('/site-packages/sympy/' not in (path or '') for path in proofs[0]['loaded_symbolic_modules'].values())


def test_missing_native_module_does_not_fall_back_to_an_installed_reference(tmp_path):
    process, proofs = run_backend(tmp_path, 'symox', '''
import importlib.abc, sys
class MissingNative(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'symox.diffgeom':
            raise ModuleNotFoundError('deliberately unavailable native module')
sys.meta_path.insert(0, MissingNative())
import sympy.diffgeom
''')
    assert process.returncode != 0
    assert 'requires missing native module symox.diffgeom' in process.stderr
    assert 'sympy.diffgeom' not in proofs[0]['loaded_symbolic_modules']


def test_preloaded_package_aliases_preserve_public_simplify_callable(tmp_path):
    process, proofs = run_backend(tmp_path, 'symox', '''
import sympy, symox
import importlib
for name in ['sympy.simplify', 'sympy.simplify.simplify', 'sympy.core.expr']:
    importlib.import_module(name)
assert sympy.simplify is symox.simplify
assert callable(sympy.simplify)
x = sympy.Symbol('x')
assert sympy.simplify((x+1)**2-x*x-2*x) == 1
''')
    assert process.returncode == 0, process.stderr
    assert proofs[0]['routed_modules']['sympy.simplify'] == 'symox.simplify'
    assert not proofs[0]['blocked_reference_execution']


def test_reference_code_execution_is_blocked_even_without_an_import(tmp_path):
    process, proofs = run_backend(tmp_path, 'symox', "exec(compile('x=1', '/env/site-packages/sympy/core/basic.py', 'exec'))")
    assert process.returncode != 0
    assert 'Reference SymPy code is forbidden' in process.stderr
    assert proofs[0]['blocked_reference_execution'] == ['/env/site-packages/sympy/core/basic.py']


def test_reference_control_really_imports_sympy_114(tmp_path):
    process, proofs = run_backend(tmp_path, 'sympy', "import sympy; assert sympy.__version__=='1.14.0'; assert sympy.Symbol.__module__.startswith('sympy.')")
    assert process.returncode == 0, process.stderr
    assert not proofs[0]['routed_modules']


def test_counts_keep_collection_skips_and_empty_reason_xfails():
    report = dict(reports=[
        dict(phase='call', outcome='passed', wasxfail=None),
        dict(phase='setup', outcome='skipped', wasxfail=None),
        dict(phase='call', outcome='skipped', wasxfail=''),
        dict(phase='setup', outcome='failed', wasxfail=None)],
        collection_skips=[dict(nodeid='optional.py')], collection_errors=[])
    assert audit.phase_summary(report) == dict(passed=1, failed=0, errors=1,
        skipped=2, xfailed=1, xpassed=0, collection_errors=0)


@pytest.mark.parametrize('change,status', [
    ({}, 'native_suite_passes'),
    ({'reference_passes_not_reproduced': {'count': 1}}, 'native_does_not_reproduce_reference_passes'),
    ({'exit_code': 1}, 'native_failure'),
    ({'provenance': {'routed_modules': {}}}, 'no_sympy_code_exercised'),
    ({'provenance': {'blocked_reference_execution': ['reference.py']}}, 'reference_execution_blocked'),
])
def test_release_gate_rejects_additional_skips_or_backend_contamination(change, status):
    execution = dict(exit_code=0, counts={'passed': 1},
        provenance={'routed_modules': {'sympy': 'symox'}, 'blocked_reference_execution': []})
    execution.update(change)
    assert audit.native_status(execution) == status


def test_plugin_records_skipped_modules_without_changing_test_outcomes(tmp_path):
    (tmp_path/'test_skip.py').write_text("import pytest\npytest.skip('optional backend absent', allow_module_level=True)\n")
    (tmp_path/'test_pass.py').write_text('def test_pass(): assert 2+2==4\n')
    report = tmp_path/'tests.json'
    env = dict(os.environ, PYTHONPATH=str(SUPPORT), SYMOX_DOWNSTREAM_TEST_REPORT=str(report))
    process = subprocess.run([sys.executable, '-m', 'pytest', '-q', '-p', 'audit_plugin', str(tmp_path)],
                             env=env, cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert process.returncode == 0, process.stdout+process.stderr
    counts = audit.phase_summary(json.loads(report.read_text()))
    assert counts['passed'] == counts['skipped'] == 1


def test_source_integrity_detects_changed_assertions_even_with_updated_run_hashes(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(SUPPORT))
    import verify_sources
    source = tmp_path/'sources/pkg'
    source.mkdir(parents=True)
    test = source/'test_math.py'
    test.write_text('def test_math(): assert 2+2 == 4\n')
    archive = tmp_path/'archives/pkg.tar.gz'
    archive.parent.mkdir()
    with tarfile.open(archive, 'w:gz') as package:
        package.add(test, arcname='pkg/test_math.py')
    project = dict(source_artifact=dict(filename=archive.name,
        digests={'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()}))
    result = dict(source=str(source), selected_test_sha256={str(test):verify_sources.sha(test.read_bytes())})
    assert verify_sources.verify(result, project, tmp_path)['status'] == 'pinned_upstream_tests_match'
    test.write_text('def test_math(): pass\n')
    result['selected_test_sha256'][str(test)] = verify_sources.sha(test.read_bytes())
    assert verify_sources.verify(result, project, tmp_path)['status'] == 'mismatch'


def test_external_success_cannot_be_reused_for_a_different_wheel(tmp_path):
    wheel = tmp_path/'candidate.whl'
    wheel.write_bytes(b'different candidate')
    project = dict(id='native_example', kind='native')
    args = SimpleNamespace(out=tmp_path/'out', cache=tmp_path/'cache', wheel=wheel)
    audit.write(args.out/'results/native_example.json', dict(project=project,
        candidate={'sha256':'old-wheel-digest'}, status='native_suite_passes'))
    result = audit.run_project(project, args)
    assert result['status'] == 'external_adapter_needs_candidate_refresh'
    assert result['candidate']['sha256'] == 'old-wheel-digest'
    assert result['requested_candidate_sha256'] == hashlib.sha256(wheel.read_bytes()).hexdigest()
