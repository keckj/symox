import subprocess
from pathlib import Path
output=Path(__file__).parent/'performance';output.mkdir(exist_ok=True)
for runner in ['log-to-atan','sparse-calculus','fraction-free']:
 for trial in range(3):
  pairs=[('before','native-calculus-20261026'),('after','native-derivative-r1b-20261027')]
  if trial%2:pairs.reverse()
  for label,snapshot in pairs:
   remote='/home/keckj/symox-ci/bench/native-derivative-r1b-20261027';name=f'{runner}-{label}-{trial}'
   cmd=f'export SYMOX_CI_PY=/home/keckj/symox-ci/wt/{snapshot}/python; export PYTHONPATH="$SYMOX_CI_PY"; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1; taskset -c 24 /home/keckj/symox-ci/venv/bin/python {remote}/benchmark-{runner}.py {remote}/{name}.json'
   with (output/(name+'.log')).open('w') as log:subprocess.run(['ssh','demeter',cmd],stdout=log,stderr=subprocess.STDOUT,check=True)
   subprocess.run(['scp',f'demeter:{remote}/{name}.json',str(output/(name+'.json'))],check=True)
   print(name,flush=True)
