"""Check an actual Lcapy RC circuit with both backends in isolated processes.

The native worker uses the repository's experimental import-routing hook. It
is a demonstration of this circuit, not a supported production Lcapy installer.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

NETLIST = 'V1 1 0 step\nR1 1 2 R\nC1 2 0 C'


def worker():
    import lcapy
    import sympy as sp
    warm = lcapy.Circuit('V1 1 0 step\nR1 1 2 1\nR2 2 0 1')
    warm[2].V.laplace()
    warm[2].V.time()
    started = time.perf_counter()
    c = lcapy.Circuit(NETLIST)
    laplace = c[2].V.laplace().sympy
    response = c[2].V.time().sympy
    elapsed = time.perf_counter()-started
    names = {str(v): v for v in laplace.free_symbols | response.free_symbols}
    s, t, R, C, V = [names[k] for k in ('s','t','R','C','V1')]
    assert sp.simplify(laplace - V/(s*(1+R*C*s))) == 0
    expected = V*(1-sp.exp(-t/(R*C)))*sp.Heaviside(t)
    assert sp.simplify(response-expected) == 0
    return dict(time=elapsed,lcapy=lcapy.__version__, laplace=str(laplace), response=str(response),
                verified=True, values=[float(response.subs({R:1000,C:sp.Rational(1,1000000),V:1,t:sp.Rational(k,10000)}))
                                       for k in range(61)])


def generate(source=None, repeat=3):
    root = Path(__file__).resolve().parents[2]
    records = {}
    with tempfile.TemporaryDirectory(prefix='symox-lcapy-') as tmp:
        for lib in ('sympy','symox'):
            env = os.environ.copy()
            paths = [str(root/'tools/downstream'), str(root/'python')]
            if source:
                paths.append(str(Path(source).resolve()))
            paths.append(env.get('PYTHONPATH',''))
            env['PYTHONPATH'] = os.pathsep.join(paths)
            env['SYMOX_DOWNSTREAM_BACKEND'] = lib
            runs=[]
            for index in range(repeat):
                env['SYMOX_DOWNSTREAM_PROVENANCE'] = str(Path(tmp)/f'{lib}-{index}.json')
                r = subprocess.run([sys.executable,__file__,'--worker'],env=env,capture_output=True,text=True,timeout=120)
                if r.returncode:
                    raise RuntimeError(r.stderr)
                row = json.loads(r.stdout)
                evidence = list(Path(tmp).glob(f'{lib}-{index}.*.json'))
                assert len(evidence) == 1
                proof = json.loads(evidence[0].read_text())
                assert not proof['blocked_reference_execution']
                if lib == 'symox':
                    assert proof['native_version']
                    assert proof['routed_modules']['sympy'] == 'symox'
                    assert all('/site-packages/sympy/' not in str(path) for path in proof['loaded_symbolic_modules'].values())
                row['provenance']=proof
                runs.append(row)
            records[lib]=min(runs,key=lambda row:row['time'])
            records[lib]['samples_seconds']=[row['time'] for row in runs]
            records[lib]['all_process_provenance']=[row['provenance'] for row in runs]
    assert max(abs(a-b) for a,b in zip(records['symox']['values'],records['sympy']['values'])) < 1e-12
    extension=next((root/'python/symox').glob('_core*.so'))
    return dict(netlist=NETLIST, records=records, verified=True,
                measured_at=datetime.now(timezone.utc).isoformat(),
                commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),
                binary_sha256=hashlib.sha256(extension.read_bytes()).hexdigest(),
                protocol='Best of 3 fresh processes after a smaller resistive-circuit warmup',
                timing_scope='Circuit construction, Laplace voltage and time response; imports, warmup, verification and samples excluded',
                scope='One RC circuit; experimental process-local backend routing',
                milliseconds=[k/10 for k in range(61)])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--lcapy-source', type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker()))
    else:
        out = generate(args.lcapy_source)
        Path(__file__).with_name('data').joinpath('downstream.json').write_text(json.dumps(out,indent=2)+'\n')
        print(out['scope'], 'verified with both backends')
