"""Run bench.py for every case: 7 fresh processes (median) and repeated calls.

python drive.py <label> <PYTHONPATH> <out.json>
"""
import json
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
label, path, out = sys.argv[1:4]
env = dict(os.environ, PYTHONPATH=path)
import importlib.util
spec = importlib.util.spec_from_file_location('cases', os.path.join(HERE, 'bench.py'))
names = ['sin(A)', 'cos(-A)', 'exp(log(A))', 'sqrt(A)', 'sin(x).subs(x,A)', 'Function(f)(A)',
         '2*x*A', 'A**2*B', 'ZeroMatrix**(1/2)', 'sin(x)*cos(x)+x**2', 'srepr(Tuple)',
         'Markov Str P', 'Markov int P']
result = {}
for name in names:
    def call(mode):
        text = subprocess.run([sys.executable, os.path.join(HERE, 'bench.py'), label, mode, name],
                              env=env, capture_output=True, text=True, check=True).stdout
        return json.loads(text)
    fresh = [call('fresh') for _ in range(7)]
    repeat = call('repeat')
    fresh = fresh[0] if isinstance(fresh[0], str) else statistics.median(fresh)
    result[name] = {'fresh': fresh, 'repeat': repeat}
    print(label, name, result[name], flush=True)
json.dump(result, open(out, 'w'), indent=1)
