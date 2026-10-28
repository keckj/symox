from pathlib import Path
import os,json,subprocess,sys,time
root=Path('/home/keckj/symox-ci/wt/domainmatrix-auto-benchmark-20261027');root.mkdir(exist_ok=True)
paths={'R5':Path('/home/keckj/symox-ci/wt/native-polymatrix-real-log-r1-20261027/python'),
       'private-disabled':root/'after-release/python','private-enabled':root/'after-release/python'}
commands=[]
for repeat in range(3):
    arms=['R5','private-disabled','private-enabled'] if repeat%2==0 else ['private-enabled','private-disabled','R5']
    for arm in arms:
        cmd=['taskset','-c','24',sys.executable,str(root/'worker.py'),arm]
        env=dict(os.environ,PYTHONPATH=str(paths[arm]))
        commands.append({'round':repeat,'arm':arm,'command':cmd,'PYTHONPATH':str(paths[arm])})
        p=subprocess.run(cmd,env=env,cwd=root,capture_output=True,text=True,timeout=30)
        (root/f'{repeat}-{arm}.json').write_text(p.stdout)
        (root/f'{repeat}-{arm}.stderr').write_text(p.stderr)
        assert p.returncode==0,(arm,p.stderr)
        d=json.loads(p.stdout);print(repeat,arm,d['certificates'],'certificates',flush=True)
(root/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
reference=json.loads((root/'0-R5.json').read_text())
for entry in commands:
 d=json.loads((root/f"{entry['round']}-{entry['arm']}.json").read_text())
 assert [c['case'] for c in d['cases']]==[c['case'] for c in reference['cases']]
 assert [(c['choices'],c['digests']) for c in d['cases']]==[(c['choices'],c['digests']) for c in reference['cases']]
print('ALL CHOICES AND EXACT CERTIFICATES IDENTICAL',flush=True)
