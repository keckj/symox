"""Prepared replay. Run only after full R6 gates and a root-authorized quiet window."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,sys

p=argparse.ArgumentParser()
p.add_argument('--baseline-python',type=Path,required=True)
p.add_argument('--r6-python',type=Path,required=True)
p.add_argument('--r6-manifest',type=Path,required=True)
p.add_argument('--r6-native',required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
assert not a.out.exists(),'Keep earlier raw results immutable'
a.out.mkdir(parents=True)
manifest=json.loads(a.r6_manifest.read_text())
assert manifest['identity']['native_sha256']==a.r6_native
r6root=a.r6_python.parent
for name,h in manifest['files'].items():assert hashlib.sha256((r6root/name).read_bytes()).hexdigest()==h,name
worker=Path(__file__).with_name('worker.py')
(a.out/'identity.json').write_text(json.dumps({'r6_identity':manifest['identity'],'baseline_python':str(a.baseline_python),
    'r6_python':str(a.r6_python),'worker_sha256':hashlib.sha256(worker.read_bytes()).hexdigest()},indent=2)+'\n')
commands=[];reference=None
for repeat in range(3):
    arms=['R5','private-disabled','private-enabled'] if repeat%2==0 else ['private-enabled','private-disabled','R5']
    for arm in arms:
        package=a.baseline_python if arm=='R5' else a.r6_python
        expected_native='532386bda46046d3d21c99f49c1b863fac62bdee5c357d03cf5123bde6913d4b' if arm=='R5' else a.r6_native
        cmd=['taskset','-c','24',sys.executable,str(worker),arm]
        commands.append({'round':repeat,'arm':arm,'meaning':'R5' if arm=='R5' else 'R6 '+arm.removeprefix('private-'),
                         'command':cmd,'PYTHONPATH':str(package),'expected_native':expected_native})
        result=subprocess.run(cmd,env=dict(os.environ,PYTHONPATH=str(package)),capture_output=True,text=True,timeout=30)
        prefix=a.out/f'{repeat}-{arm}'
        prefix.with_suffix('.json').write_text(result.stdout);prefix.with_suffix('.stderr').write_text(result.stderr)
        assert result.returncode==0,result.stderr
        data=json.loads(result.stdout)
        assert Path(data['package']).resolve()==(package/'symox/__init__.py').resolve()
        assert data['native_sha256']==expected_native and data['affinity']==[24]
        signatures=[(c['case'],c['choices'],c['digests']) for c in data['cases']]
        if reference is None:reference=signatures
        assert signatures==reference
        print(repeat,arm,data['certificates'],'certificates',flush=True)
(a.out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
for name,h in manifest['files'].items():assert hashlib.sha256((r6root/name).read_bytes()).hexdigest()==h,name
print('ALL CHOICES AND EXACT CERTIFICATES IDENTICAL; SOURCES UNCHANGED')
