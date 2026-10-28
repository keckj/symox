from pathlib import Path
import hashlib,json,subprocess,time
root=Path('/home/keckj/symox-ci/wt/native-lrt-limits-20261027')
output=Path('/home/keckj/symox-ci/wt/native-lrt-limits-20261027/.validation')
expected=json.loads((output/'rust-build-input.json').read_text())
result={}
for name,command in [('tests',['cargo','test','-p','symox']),('clippy',['cargo','clippy','-p','symox','--','-D','clippy::correctness'])]:
    start=time.time()
    with (output/(name+'.log')).open('w') as log:
        process=subprocess.run(command,cwd=root,stdout=log,stderr=subprocess.STDOUT)
    result[name]={'command':command,'returncode':process.returncode,'elapsed_seconds':time.time()-start}
result['rust_inputs_unchanged']=all(hashlib.sha256((root/p).read_bytes()).hexdigest()==h for p,h in expected.items())
(output/'rustchecks.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
raise SystemExit(0 if result['rust_inputs_unchanged'] and all(result[n]['returncode']==0 for n in ['tests','clippy']) else 1)
