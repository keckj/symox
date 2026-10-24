"""Refresh every displayed measurement, verify results, then publish the data.

    .venv/bin/python examples/showcase_v4/refresh.py [--resume]

Build Symox in release mode first. Workers are strictly sequential so benchmark
processes do not contend with each other. Old measurements/video are archived.
Partial runs stay in data-refresh/ and cannot silently feed the renderer.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib
from itertools import permutations
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys

import bench_modules
from bench_current import INTEGRANDS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LIBS = ('symox','sympy')
CORE = ['integrals']+[f'det{n}' for n in range(3,8)]+['taylor_exp','solve','throughput','pendulum6']
MODULES = ['normal_moment','normal_var','binom_var','binom_skew','dice','exp_P','gamma_mgf',
           'sum_normals','beta_E','normal_entropy','beam','hydrogen','tf','pauli','lagrange2',
           'wigner9j','cg','berlekamp','hill','lfsr','vigenere','rsa']


def write(path, data):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(data,indent=1)+'\n')
    temp.replace(path)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def measure(lib, case):
    samples, best = [], None
    for _ in range(1 if case in ('integrals','throughput') else 3):
        try:
            r = subprocess.run([sys.executable,str(HERE/'bench_current.py'),lib,case],
                               capture_output=True,text=True,timeout=150)
        except subprocess.TimeoutExpired:
            return dict(status='process_timeout', budget_seconds=150,
                        detail='Worker exceeded setup + operation + serialization budget; operation timing unavailable')
        if r.returncode:
            raise RuntimeError(f'{lib} {case}: {r.stderr[-2000:]}')
        row = json.loads(r.stdout)
        if 'timeout' in row:
            return row
        elapsed = row.get('time',row.get('elapsed',sum(x['time'] for x in row.get('items',[]))))
        samples.append(elapsed)
        if best is None or elapsed < best[0]:
            best = elapsed,row
        if elapsed > 5:
            break
    best[1]['samples_seconds'] = samples
    return best[1]


def verify(core, modules):
    import sympy as sp
    import mpmath as mp
    x = sp.Symbol('x', real=True)
    checks = {}
    for lib in LIBS:
        items = core[lib]['integrals']['items']
        assert len(items) == len(INTEGRANDS)
        for item in items:
            f = sp.sympify(item['input'], locals={'x':x})
            F = sp.sympify(item['result'], locals={'x':x})
            # Classical antiderivatives are checked on smooth real intervals.
            # A floor of a real argument is locally constant away from its
            # jumps; its distributional derivative is not being asserted zero.
            floors = F.atoms(sp.floor)
            assert all(v.args[0].is_real is True for v in floors)
            smooth = F.xreplace({v:sp.Dummy('locally_constant',real=True) for v in floors})
            residual = sp.diff(smooth,x)-f
            zero = sp.simplify(sp.trigsimp(residual)) == 0
            if not zero:
                zero = sp.cancel(sp.expand(residual.rewrite(sp.exp))) == 0
            assert zero, (lib,item['input'])
            item['verified'] = 'exact symbolic differentiation on smooth real intervals, away from poles and branch boundaries'
        print(lib, '63 antiderivatives verified',flush=True)
        for n in range(3,8):
            row = core[lib][f'det{n}']
            if 'time' not in row:
                assert 'timeout' in row
                continue
            symbols = [[sp.Symbol(f'a{i}{j}') for j in range(n)] for i in range(n)]
            expected = {}
            for p in permutations(range(n)):
                monomial = math.prod(symbols[i][p[i]] for i in range(n))
                sign = (-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
                expected[monomial] = sign
            # Expanded determinants are flat sums; parse their monomials one
            # at a time to avoid Python's compiler recursion limit at 7! terms.
            actual = {}
            for word in re.findall(r'[+-]?[^+-]+',row['result'].replace(' ','')):
                coefficient,monomial = sp.sympify(word).as_coeff_Mul()
                actual[monomial] = actual.get(monomial,0)+coefficient
            assert actual == expected, (lib,n)
            row['verified'] = 'all monomials and signs against the Leibniz formula'
            row.pop('result')
        series = core[lib]['taylor_exp']
        if 'time' in series:
            actual = sp.sympify(series['result'],locals={'x':x}).removeO()
            # Independent coefficient recurrence: f'=cos(x)*f, f(0)=1.
            coeffs = [sp.S.One]
            for n in range(24):
                coeffs.append(sum((-1)**(j//2)*coeffs[n-j]/sp.factorial(j)
                                  for j in range(0,n+1,2))/(n+1))
            assert sp.expand(actual-sum(c*x**i for i,c in enumerate(coeffs))) == 0
            series['verified'] = 'exact recurrence f prime = cos(x) f through order 24'
        roots = core[lib]['solve']
        if roots.get('status') == 'process_timeout':
            print(lib, 'solve incomplete: process budget, no operation ratio',flush=True)
            continue
        assert roots['count'] == 8 and len(roots['real']) == 2
        with mp.workdps(70):
            nums = [[mp.mpc(str(sp.re(sp.sympify(a)).evalf(65)),str(sp.im(sp.sympify(a)).evalf(65)))
                     for a in pair] for pair in roots['numeric']]
            assert min(abs(a[0]-b[0]) for i,a in enumerate(nums) for b in nums[i+1:]) > mp.mpf('1e-20')
            assert all(abs(a*a+b*b-4)<mp.mpf('1e-48') and abs(b-a**4+1)<mp.mpf('1e-48') for a,b in nums)
        roots['verified'] = '8 distinct roots; both residuals < 1e-48 at 60 digits; elimination degree 8'
    # Compare entire returned sets, rather than checking only that some roots work.
    if all('numeric' in core[lib]['solve'] for lib in LIBS):
        a,b = [core[lib]['solve']['numeric'] for lib in LIBS]
        ac = [[complex(sp.N(sp.sympify(v),16)) for v in pair] for pair in a]
        bc = [[complex(sp.N(sp.sympify(v),16)) for v in pair] for pair in b]
        assert all(min(max(abs(u-v) for u,v in zip(p,q)) for q in bc)<1e-12 for p in ac)
        checks['solve'] = 'whole returned sets matched; numerical residual validation, not a general exact solver certificate'
    else:
        checks['solve'] = 'Native worker exceeded 150 s process budget; no native solution or operation timing claim. Completed reference result checked separately.'
    for a,b in zip(core['symox']['pendulum6']['probes'],core['sympy']['pendulum6']['probes']):
        assert max(abs(u-v) for u,v in zip(a,b)) < 1e-8
    checks['pendulum'] = 'accelerations agree at two states, tolerance 1e-8; no trajectory-equivalence claim'
    a,b = [core[lib]['throughput'] for lib in LIBS]
    assert len(a['results']) == len(b['results']) == 40
    for u,v,srca,srcb in zip(a['results'],b['results'],a['srcs'],b['srcs']):
        assert sp.expand(sp.sympify(srca)-sp.sympify(srcb)) == 0
        assert sp.expand(sp.sympify(u)-sp.sympify(v)) == 0
    checks['throughput'] = 'identical deterministic input stream; first 40 outputs exactly matched'
    for case,row in modules.items():
        if case == 'normal_entropy' and 'timeout' in row['sympy']:
            sigma = sp.Symbol('sigma',positive=True)
            h = sp.sympify(row['symox']['result'],locals={'sigma':sigma})
            assert sp.simplify(h-(1+sp.log(2*sp.pi*sigma**2))/2) == 0
            row['verified'] = 'independent Gaussian entropy formula in exactness.json'
        else:
            assert row['same'], (case,'different results')
            row['verified'] = 'matching strings or symbolic difference simplified to zero'
    # The timed LFSR returns a sum and the timed RSA case returns its last
    # round trip. Check the full streams/round trips separately as well.
    streams = {}
    for lib in LIBS:
        sp = importlib.import_module(lib)
        crypto = importlib.import_module(lib+'.crypto.crypto')
        F = sp.FF(2)
        key = [F(v) for v in (1,0,0,1,0,1,1)]
        streams[lib] = [int(v) for v in crypto.lfsr_sequence(key,[F(1)]*7,3000)]
        P,Q,E = sp.nextprime(10**60),sp.nextprime(10**61),65537
        public = crypto.rsa_public_key(P,Q,E)
        private = crypto.rsa_private_key(P,Q,E)
        assert len(str(int(P*Q))) == 122
        assert all(int(crypto.decipher_rsa(crypto.encipher_rsa(m,public),private)) == m for m in range(1,200))
    assert streams['symox'] == streams['sympy']
    checks['crypto'] = dict(lfsr='all 3000 bits agree exactly',
                            stream_sha256=hashlib.sha256(bytes(streams['symox'])).hexdigest(),
                            rsa='all 199 round trips recover the plaintext; modulus has 122 digits')
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--lcapy-source',type=Path)
    parser.add_argument('--sage-prefix',type=Path,required=True)
    parser.add_argument('--lcapy-test-prefix',type=Path,required=True)
    args = parser.parse_args()
    import symox, sympy
    from symox import _core
    stage = HERE/'data-refresh'
    stage.mkdir(exist_ok=True)
    identity = dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                    binary_sha256=digest(Path(_core.__file__)),
                    working_code_diff_sha256=hashlib.sha256(subprocess.check_output(['git','diff','--binary','HEAD','--','crates','python'],cwd=ROOT)).hexdigest(),
                    inputs_sha256={name:digest(HERE/name) for name in ('bench_current.py','bench_modules.py')})
    if args.resume:
        metadata = json.loads((stage/'metadata.json').read_text())
        assert metadata['identity'] == identity, 'Code/build changed; start a fresh measurement run'
    else:
        metadata = dict(identity=identity, measured_at=datetime.now(timezone.utc).isoformat(),
            versions=dict(symox=symox.__version__,sympy=sympy.__version__),python=sys.version,
            cpu=next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                      if line.startswith('model name')),platform.machine()),
            protocol='Sequential fresh processes; smaller warm-up; best of up to 3; >5 s once; 120 s timed-operation alarm; separate 150 s process budget. Integration wall and 6 s throughput once.',
            build='maturin develop --release', selection='Selected examples; not an aggregate performance or global parity claim')
        write(stage/'metadata.json',metadata)
    core = json.loads((stage/'core-records.json').read_text()) if args.resume and (stage/'core-records.json').exists() else {lib:{} for lib in LIBS}
    modules = json.loads((stage/'cases4.json').read_text()) if args.resume and (stage/'cases4.json').exists() else {}
    for case in CORE:
        for lib in LIBS:
            if case not in core[lib]:
                print('measuring',lib,case,flush=True)
                core[lib][case] = measure(lib,case)
                write(stage/'core-records.json',core)
                print('done',lib,case,core[lib][case].get('time',core[lib][case].get('timeout','wall/stream')),flush=True)
    for case in MODULES:
        if case in modules:
            continue
        records = {}
        for lib in LIBS:
            samples=[]
            for _ in range(3):
                print('measuring',lib,case,flush=True)
                r=bench_modules.run(lib,case)
                if 'time' not in r:
                    records[lib]=r
                    break
                samples.append(r['time'])
                if lib not in records or r['time']<records[lib]['time']:
                    records[lib]=r
                if r['time']>5:
                    break
            records[lib]['samples_seconds']=samples
        a,b=records['symox'],records['sympy']
        records['same']='result' in a and 'result' in b and bench_modules.same(a['result'],b['result'])
        modules[case]=records
        write(stage/'cases4.json',modules)
        print('done',case,'same',records['same'],flush=True)
    print('verifying results independently',flush=True)
    checks=verify(core,modules)
    assert all('time' in core[lib]['solve'] for lib in LIBS), 'Required checked solver regression did not complete'
    from verify_exactness import generate as exactness
    from verify_downstream import generate as downstream
    write(stage/'exactness.json',exactness())
    write(stage/'downstream.json',downstream(args.lcapy_source))
    from verify_sage import generate as sage_tests
    from verify_lcapy_tests import generate as lcapy_tests
    sage=sage_tests(args.sage_prefix,stage/'sage')
    assert sage['verified'], 'Sage bridge tests did not reproduce the reference'
    assert args.lcapy_source, 'Provide the Lcapy source checkout for upstream tests'
    lcapy=lcapy_tests(args.lcapy_test_prefix,args.lcapy_source,stage/'lcapy-tests')
    assert lcapy['verified'], 'Lcapy tests did not reproduce the reference'
    write(stage/'sage.json',sage)
    write(stage/'lcapy-tests.json',lcapy)
    write(stage/'verification.json',checks)
    write(stage/'cases4.json',modules)
    cases={lib:{k:v for k,v in rows.items() if k not in ('throughput','pendulum6')} for lib,rows in core.items()}
    cases['versions']=metadata['versions']
    write(stage/'cases3.json',cases)
    for lib in LIBS:
        write(stage/f'tp_{lib}.json',core[lib]['throughput'])
    data=json.loads((HERE/'data/data.json').read_text())
    data.update(versions=metadata['versions'],cpu=metadata['cpu'],measured_at=metadata['measured_at'],
                race={lib:{'6':core[lib]['pendulum6']['time']} for lib in LIBS})
    # Frames are an illustrative simulation, not timing data. Keep the existing
    # unit-pendulum trajectory, explicitly separate from this refreshed race.
    data.pop('taylor',None)
    data['frames_note']='Illustrative unit-pendulum trajectory retained from original render; readiness times freshly measured.'
    write(stage/'data.json',data)
    archive=HERE/'archive/2026-10-18'
    if not archive.exists():
        archive.mkdir(parents=True)
        shutil.copytree(HERE/'data',archive/'data')
        shutil.copy2(HERE/'symox_showcase_v4.mp4',archive/'symox_showcase_v4.mp4')
    # All validation has passed before any renderer inputs are replaced.
    assert digest(Path(_core.__file__))==identity['binary_sha256'],'Native binary changed during measurements'
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==identity['commit'],'Source commit changed during measurements'
    patch=subprocess.check_output(['git','diff','--binary','HEAD','--','crates','python'],cwd=ROOT)
    assert hashlib.sha256(patch).hexdigest()==identity['working_code_diff_sha256'],'Native/Python source changed during measurements'
    for name in ['cases3.json','cases4.json','tp_symox.json','tp_sympy.json','data.json','metadata.json',
                 'exactness.json','downstream.json','verification.json','sage.json','lcapy-tests.json']:
        shutil.copy2(stage/name,HERE/'data'/name)
    evidence=HERE/'data/evidence'
    evidence.mkdir(exist_ok=True)
    (evidence/'working-code.patch').write_bytes(patch)
    for subdir in ['sage','lcapy-tests']:
        shutil.copytree(stage/subdir,evidence/subdir,dirs_exist_ok=True)
    subprocess.run([sys.executable,str(HERE/'prepare_modules.py')],check=True)
    freeze_sources()
    print('Verified data published. Original measurements and video archived.',flush=True)


def freeze_sources():
    data = HERE/'data'
    sources = data/'sources'
    sources.mkdir(exist_ok=True)
    paths = list(HERE.glob('*.py')) + [HERE.parent/'pendulum_race.py', ROOT/'tools/downstream/sitecustomize.py',
                                    ROOT/'crates/symox-core/src/poly/roots.rs',ROOT/'crates/symox-core/src/solve.rs',
                                    ROOT/'tests/test_solve_decomposed_radicals.py',
                                    ROOT/'tests/test_solve_checks.py',ROOT/'tests/test_solve_sets_parity.py',
                                    ROOT/'python/symox/_general_simplify.py',ROOT/'tests/test_general_simplify_controls.py']
    metadata = json.loads((data/'metadata.json').read_text())
    metadata['source_sha256'] = {str(path.relative_to(ROOT)):digest(path) for path in sorted(paths)}
    for path in paths:
        destination = sources/path.relative_to(ROOT)
        destination.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,destination)
    write(data/'metadata.json',metadata)
    write(data/'manifest.json',{str(p.relative_to(data)):digest(p) for p in sorted(data.rglob('*'))
                               if p.is_file() and p.name!='manifest.json' and '__pycache__' not in p.parts})


if __name__=='__main__':
    main()
