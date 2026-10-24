"""Run the unmodified Sage SymPy bridge doctests on both symbolic backends.

    python examples/showcase_v4/verify_sage.py --prefix /path/to/sage-prefix

The prefix must already contain Sage and SymPy 1.14.0. The repository package
and experimental routing hook are selected for this process only. No packages
or upstream tests in the prefix are changed.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate(prefix, destination):
    root = Path(__file__).resolve().parents[2]
    prefix = Path(prefix).resolve()
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    python = prefix/'bin/python'
    site = next((prefix/'lib').glob('python*/site-packages'))
    test = site/'sage/interfaces/sympy.py'
    reference = json.loads((root/'audits/downstream/results/sagemath.json').read_text())
    assert sha(test) in reference['selected_tests'].values(), 'Upstream test differs from the audited selection'
    extension = next((root/'python/symox').glob('_core*.so'))
    metadata = dict(started_at=datetime.now(timezone.utc).isoformat(),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        binary_sha256=sha(extension), test=str(test), test_sha256=sha(test),
        scope='Sage sage/interfaces/sympy.py doctests only; experimental backend routing', records={})
    mamba = Path.home()/'.cache/symox-validation/bin/micromamba'
    runtime = [str(mamba),'run','-p',str(prefix)] if mamba.exists() else []
    command = runtime+[str(python),'-m','sage.doctest','--random-seed=0',
        '--optional=sage,sympy','--timeout=600',str(test)]
    for lib in ('sympy','symox'):
        for old in destination.glob(f'{lib}-provenance.*.json'):
            old.unlink()
        env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(root/'tools/downstream'),str(root/'python')]),
            SYMOX_DOWNSTREAM_BACKEND=lib,
            SYMOX_DOWNSTREAM_PROVENANCE=str(destination/f'{lib}-provenance.json'),
            MPLBACKEND='Agg', QT_QPA_PLATFORM='offscreen', PYTHONHASHSEED='0',
            OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2')
        preflight = subprocess.run(runtime+[str(python),'-c',
            'import sympy,symox,sage.version; assert sympy.__version__=="1.14.0"; print(sage.version.version); print(sympy.__version__,symox.__version__); print(symox.__file__)'],
            env=env,cwd=root,capture_output=True,text=True,timeout=60)
        assert preflight.returncode == 0, preflight.stderr
        (destination/f'{lib}-preflight.log').write_text(preflight.stdout+preflight.stderr)
        print('running Sage bridge tests:',lib,flush=True)
        started = time.perf_counter()
        with (destination/f'{lib}.log').open('w') as log:
            run = subprocess.run(command,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=900)
        elapsed = time.perf_counter()-started
        output = (destination/f'{lib}.log').read_text()
        match = re.search(r'\[(\d+) tests(?:, (\d+) failures)?[,\]]',output)
        assert match, output[-3000:]
        total,failed = int(match[1]),int(match[2] or 0)
        proofs = [json.loads(p.read_text()) for p in destination.glob(f'{lib}-provenance.*.json')]
        assert proofs and not any(p['blocked_reference_execution'] for p in proofs)
        if lib == 'symox':
            assert any(p['routed_modules'].get('sympy') == 'symox' for p in proofs)
            assert all('/site-packages/sympy/' not in str(path)
                       for p in proofs for path in p['loaded_symbolic_modules'].values())
        metadata['records'][lib] = dict(command=command,exit_code=run.returncode,
            total=total,passed=total-failed,failed=failed,wall_seconds=elapsed,
            log=f'{lib}.log',preflight=preflight.stdout,process_provenance=proofs)
        (destination/'result.json').write_text(json.dumps(metadata,indent=2)+'\n')
        print(lib, total-failed,'/',total,'passed;',elapsed,'seconds',flush=True)
        if lib == 'sympy':
            assert not failed and run.returncode == 0, 'Reference did not qualify'
    assert sha(extension) == metadata['binary_sha256'], 'Native binary changed during the run'
    assert sha(test) == metadata['test_sha256'], 'Upstream tests changed during the run'
    metadata['verified'] = all(r['failed']==0 and r['exit_code']==0 for r in metadata['records'].values())
    (destination/'result.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return metadata


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix',required=True,type=Path)
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('data-refresh')/'sage')
    args=parser.parse_args()
    generate(args.prefix,args.out)
