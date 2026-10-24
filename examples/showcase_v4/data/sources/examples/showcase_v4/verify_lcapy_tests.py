"""Run Lcapy's unmodified upstream tests with the two symbolic backends."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import xml.etree.ElementTree as ET


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate(prefix,source,destination):
    root=Path(__file__).resolve().parents[2]
    prefix,source,destination=map(lambda p:Path(p).resolve(),(prefix,source,destination))
    destination.mkdir(parents=True,exist_ok=True)
    tests={str(p.relative_to(source)):sha(p) for p in sorted((source/'lcapy/tests').rglob('*.py'))}
    extension=next((root/'python/symox').glob('_core*.so'))
    result=dict(started_at=datetime.now(timezone.utc).isoformat(),binary_sha256=sha(extension),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
        upstream_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),
        test_sha256=tests,scope='Lcapy upstream lcapy/tests selection; experimental backend routing',records={})
    project=next(p for p in json.loads((root/'tools/downstream/projects.json').read_text())['projects']
                 if p['id']=='lcapy')
    assert result['upstream_commit']==project['source_git']['commit'],'Different upstream revision'
    assert not subprocess.check_output(['git','diff','HEAD','--','lcapy'],cwd=source),'Modified upstream package or tests'
    for lib in ('sympy','symox'):
        for old in destination.glob(f'{lib}-provenance.*.json'):
            old.unlink()
        env=dict(os.environ,PYTHONPATH=os.pathsep.join([str(root/'tools/downstream'),str(root/'python'),str(source)]),
            SYMOX_DOWNSTREAM_BACKEND=lib,SYMOX_DOWNSTREAM_PROVENANCE=str(destination/f'{lib}-provenance.json'),
            MPLBACKEND='Agg',PYTHONHASHSEED='0',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
        preflight=subprocess.run([str(prefix/'bin/python'),'-c',
            'import lcapy,sympy,symox; assert lcapy.__version__=="1.26"; '
            'assert sympy.__version__=="1.14.0"; '
            'print(lcapy.__version__,sympy.__version__,symox.__version__); print(symox.__file__)'],
            cwd=source,env=env,capture_output=True,text=True,timeout=60)
        assert preflight.returncode==0,preflight.stderr
        (destination/f'{lib}-preflight.log').write_text(preflight.stdout+preflight.stderr)
        command=[str(prefix/'bin/python'),'-m','pytest','-q','-ra','-p','audit_plugin','--timeout=60',
                 '--junitxml='+str(destination/f'{lib}-junit.xml'),'lcapy/tests','--tb=native']
        print('running Lcapy upstream tests:',lib,flush=True)
        started=time.perf_counter()
        with (destination/f'{lib}.log').open('w') as log:
            run=subprocess.run(command,cwd=source,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=900)
        elapsed=time.perf_counter()-started
        xml=ET.parse(destination/f'{lib}-junit.xml').getroot()
        cases=list(xml.iter('testcase'))
        counts={key:sum(node.find(tag) is not None for node in cases)
                for key,tag in [('failed','failure'),('errors','error'),('skipped','skipped')]}
        counts['total']=len(cases)
        counts['passed']=len(cases)-sum(counts[k] for k in ('failed','errors','skipped'))
        proofs=[json.loads(p.read_text()) for p in destination.glob(f'{lib}-provenance.*.json')]
        assert proofs and not any(p['blocked_reference_execution'] for p in proofs)
        if lib=='symox':
            assert any(p['routed_modules'].get('sympy')=='symox' for p in proofs)
            assert all('/site-packages/sympy/' not in str(path)
                       for p in proofs for path in p['loaded_symbolic_modules'].values())
        result['records'][lib]=dict(command=command,exit_code=run.returncode,wall_seconds=elapsed,
                                   **counts,preflight=preflight.stdout,process_provenance=proofs,
                                   test_ids=sorted((node.get('classname'),node.get('name')) for node in cases))
        (destination/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print(lib,counts,elapsed,'seconds',flush=True)
        if lib=='sympy':
            assert run.returncode==0 and counts['passed']==317 and counts['skipped']==0,'Reference did not qualify'
    assert tests=={str(p.relative_to(source)):sha(p) for p in sorted((source/'lcapy/tests').rglob('*.py'))}
    assert sha(extension)==result['binary_sha256'],'Native binary changed during run'
    assert result['records']['sympy']['test_ids']==result['records']['symox']['test_ids'],'Different test selections'
    result['verified']=all(row['exit_code']==0 and row['passed']==317 and row['skipped']==0 for row in result['records'].values())
    (destination/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prefix',required=True,type=Path)
    p.add_argument('--source',required=True,type=Path)
    p.add_argument('--out',type=Path,default=Path(__file__).with_name('data-refresh')/'lcapy-tests')
    args=p.parse_args()
    generate(args.prefix,args.source,args.out)
