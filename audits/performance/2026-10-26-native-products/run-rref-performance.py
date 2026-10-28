import subprocess,json
from pathlib import Path
output=Path(__file__).parent/'rref-performance'
output.mkdir(exist_ok=True)
for trial in range(3):
    for label,snapshot in ([('before','native-corrections-20261026-r1'),('after','native-products-20261026')] if trial%2==0 else [('after','native-products-20261026'),('before','native-corrections-20261026-r1')]):
        remote='/home/keckj/symox-ci/bench/ddm-matmul-20261026'
        name=f'rref-{label}-{trial}'
        command=f'export SYMOX_CI_PY=/home/keckj/symox-ci/wt/{snapshot}/python; export PYTHONPATH="$SYMOX_CI_PY"; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1; taskset -c 24 /home/keckj/symox-ci/venv/bin/python {remote}/benchmark-ddm-rref.py {remote}/{name}.json'
        with (output/(name+'.log')).open('w') as log:
            subprocess.run(['ssh','demeter',command],stdout=log,stderr=subprocess.STDOUT,check=True)
        subprocess.run(['scp',f'demeter:{remote}/{name}.json',str(output/(name+'.json'))],check=True)
        print(name,flush=True)
