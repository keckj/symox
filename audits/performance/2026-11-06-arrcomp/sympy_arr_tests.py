"""sympy_arr_tests.py <worktree> <outdir>: runs the 21 SymPy 1.14 array /
array-expression / matrix-expression test files with `sympy` imports
rewritten to `symox`; prints one summary line per file and the total."""
import glob, os, re, subprocess, sys
wt, out = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
import sympy
base = os.path.dirname(sympy.__file__)
groups = [('arr', 'tensor/array/tests', None), ('expr', 'tensor/array/expressions/tests', None),
          ('mat', 'matrices/expressions/tests', ['test_applyfunc.py', 'test_derivatives.py', 'test_hadamard.py',
                                                  'test_kronecker.py', 'test_matexpr.py', 'test_trace.py'])]
os.makedirs(out, exist_ok=True)
tot_p = tot_n = 0
for pre, d, only in groups:
    for f in sorted(glob.glob(os.path.join(base, d, 'test_*.py'))):
        name = os.path.basename(f)
        if only and name not in only:
            continue
        src = open(f).read()
        src = re.sub(r'\bsympy\b', 'symox', src)
        dst = os.path.join(out, '%s_%s' % (pre, name))
        open(dst, 'w').write(src)
        env = dict(os.environ, PYTHONPATH=os.path.join(wt, 'python'))
        r = subprocess.run([sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider', '--rootdir', out,
                            '-o', 'addopts=', dst, '-rfE', '--tb=line'], cwd=out, env=env, capture_output=True, text=True)
        open(dst + '.log', 'w').write(r.stdout + r.stderr)
        last = [l for l in r.stdout.splitlines() if l.strip()][-1] if r.stdout.strip() else 'no output'
        p = int(m.group(1)) if (m := re.search(r'(\d+) passed', last)) else 0
        n = sum(int(x) for x in re.findall(r'(\d+) (?:passed|failed|xfailed|error)', last))
        print('%s_%s | %s' % (pre, name, last.strip('= ')))
        tot_p += p
print('TOTAL passed', tot_p)
