import subprocess
from pathlib import Path
output=Path(__file__).parent/'performance';output.mkdir(exist_ok=True)
for trial in range(3):
 for label,snapshot in ([('before','native-products-20261026'),('after','native-rational-20261026')] if trial%2==0 else [('after','native-rational-20261026'),('before','native-products-20261026')]):
  remote='/home/keckj/symox-ci/bench/rational-reduction-20261026';name=f'{label}-{trial}'
  cmd=f'export SYMOX_CI_PY=/home/keckj/symox-ci/wt/{snapshot}/python; export PYTHONPATH="$SYMOX_CI_PY"; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1; taskset -c 24 /home/keckj/symox-ci/venv/bin/python {remote}/benchmark-rational-reduction.py {remote}/{name}.json'
  with (output/(name+'.log')).open('w') as log:subprocess.run(['ssh','demeter',cmd],stdout=log,stderr=subprocess.STDOUT,check=True)
  subprocess.run(['scp',f'demeter:{remote}/{name}.json',str(output/(name+'.json'))],check=True)
  print(name,flush=True)
