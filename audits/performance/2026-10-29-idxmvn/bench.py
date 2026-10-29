"""Before/after timings of the matrix-head calls of the idxmvn wave.

Usage: PYTHONPATH=<tree>/python python bench.py <label>
Prints one JSON line per call: first call in a fresh process (median of
5 processes) and repeated call (median of 7 runs of 200 calls).
"""
import json
import statistics
import subprocess
import sys
import time

SETUP = ("from symox import *\n"
         "A=MatrixSymbol('A',3,3);B=MatrixSymbol('B',3,3);C=MatrixSymbol('C',3,3);x=Symbol('x')\n")
CALLS = {
    'conjugate(A)': 'conjugate(A)',
    'conjugate(x*A)': 'conjugate(x*A)',
    'conjugate(A**2)': 'conjugate(A**2)',
    'str(sin(A - A.T))': 'str(sin(A - A.T))',
    'srepr(sin(2*A - B))': 'srepr(sin(2*A - B))',
    're(A)': 're(A)',
    'im(x*A)': 'im(x*A)',
    'Abs(A)': 'Abs(A)',
    'Abs(x*A*B)': 'Abs(x*A*B)',
    're(log(A))': 're(log(A))',
    'Heaviside(A) error': 'Heaviside(A)',
    'x + A error': 'x + A',
}


def guarded(code):
    return "try:\n    " + code + "\nexcept Exception:\n    pass\n"


def first_call(code):
    script = (SETUP + "import time\nt=time.perf_counter()\n" + guarded(code)
              + "print(time.perf_counter()-t)\n")
    samples = []
    for _ in range(5):
        out = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True,
                             check=True).stdout
        samples.append(float(out.strip().splitlines()[-1]))
    return statistics.median(samples)


def repeated(code):
    namespace = {}
    exec(SETUP, namespace)
    body = compile(guarded(code), '<bench>', 'exec')
    runs = []
    for _ in range(7):
        t = time.perf_counter()
        for _ in range(200):
            exec(body, namespace)
        runs.append((time.perf_counter() - t)/200)
    return statistics.median(runs)


def outcome(code):
    namespace = {}
    exec(SETUP, namespace)
    try:
        value = eval(code, namespace)
        return str(value)[:80]
    except Exception as error:
        return '!! %s: %s' % (type(error).__name__, str(error)[:60])


if __name__ == '__main__':
    label = sys.argv[1]
    for name, code in CALLS.items():
        print(json.dumps({'label': label, 'call': name, 'first_s': first_call(code),
                          'repeated_s': repeated(code), 'result': outcome(code)}), flush=True)
