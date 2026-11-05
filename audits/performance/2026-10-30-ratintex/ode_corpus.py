"""Replay the 315-case ODE corpus on symox only and score each method against
the stored SymPy 1.14 observations of audits/sympy-1.14.0/corpus-2026-10-26/ode-570f.
usage: (from worktree) PYTHONPATH=python python ode_corpus.py OUT.json [--jobs 16] [--timeout 60]
"""
import argparse, concurrent.futures, glob, gzip, json, os, re, select, subprocess, sys
from pathlib import Path

ROOT = Path.cwd()
CORPUS = ROOT / 'audits/sympy-1.14.0/corpus-2026-10-26'
WORKER = re.search(r"WORKER=r'''(.*?)'''", (CORPUS / 'replay-ode.py').read_text(), re.S).group(1)


def stored():
    ref = {}
    for p in sorted(glob.glob(str(CORPUS / 'ode-570f/[0-9][0-9].json.gz'))):
        for row in json.loads(gzip.open(p).read()):
            ref[(row['id'], row['method'])] = row['sympy']
    return ref


def same(a, b):
    if a.get('error') in ('DEADLINE', 'MemoryError', 'WORKER_EXIT', 'PROTOCOL'):
        return False
    if a['ok'] and b['ok']:
        return a['text'] == b['text']
    return not a['ok'] and not b['ok'] and (a['error'], a['text']) == (b['error'], b['text'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--jobs', type=int, default=16)
    ap.add_argument('--timeout', type=float, default=60)
    args = ap.parse_args()
    cases = json.loads(gzip.open(CORPUS / 'ode-cases.json.gz').read())['cases']
    ref = stored()

    def shard(index):
        rows = []
        proc = None

        def start():
            p = subprocess.Popen([sys.executable, '-u', '-c', WORKER, 'symox'], stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, env=os.environ.copy())
            json.loads(p.stdout.readline())
            return p
        for case in cases[index::args.jobs]:
            for method in ('classify', 'hint', 'default'):
                if proc is None or proc.poll() is not None:
                    proc = start()
                proc.stdin.write(json.dumps(dict(case, method=method)) + '\n')
                proc.stdin.flush()
                if not select.select([proc.stdout], [], [], args.timeout)[0]:
                    proc.kill(); proc.wait(); proc = None
                    res = {'ok': False, 'error': 'DEADLINE', 'text': ''}
                else:
                    line = proc.stdout.readline()
                    res = json.loads(line) if line else {'ok': False, 'error': 'WORKER_EXIT', 'text': ''}
                    if not line:
                        proc = None
                rows.append({'id': case['id'], 'method': method, 'same': same(res, ref[(case['id'], method)]), **res})
        if proc is not None:
            proc.kill(); proc.wait()
        return rows
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as ex:
        rows = sum(ex.map(shard, range(args.jobs)), [])
    key = {'classify': 'cls', 'hint': 'hint', 'default': 'def'}
    score = {key[m]: sum(r['same'] for r in rows if r['method'] == m) for m in key}
    score['n'] = len(cases)
    Path(args.out).write_text(json.dumps({'score': score, 'rows': rows}, indent=1))
    print(json.dumps(score))


if __name__ == '__main__':
    main()
