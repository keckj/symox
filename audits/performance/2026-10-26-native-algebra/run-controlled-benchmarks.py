from pathlib import Path
import subprocess
cache=Path('/home/keckj/.cache/symox-validation/native-algebra-20261026')
packages={'before':'/home/keckj/symox-ci/wt/native-kernels-20261026/python','after':'/home/keckj/symox-ci/wt/native-algebra-20261026/python'}
for index in range(3):
 for backend in (['after'] if index==0 else ['before','after']):
  for kind,script in [('count','operation-count-algebra-benchmark.py'),('simplex','simplex-algebra-benchmark.py')]:
   env=packages[backend]
   command=f'SYMOX_CI_PY={env} PYTHONPATH={env} taskset -c 24 /home/keckj/symox-ci/venv/bin/python /home/keckj/symox-ci/{script}'
   target=cache/f'{kind}-{backend}-{index}.json'
   assert not target.exists(),target
   with target.open('w') as out,(cache/f'{kind}-{backend}-{index}.stderr.log').open('w') as err:
    subprocess.run(['ssh','-o','BatchMode=yes','demeter',command],stdout=out,stderr=err,check=True)
   print(kind,backend,index,flush=True)
