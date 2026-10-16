"""Measure API presence and a curated differential corpus, not global parity.

Run with a development installation and SymPy installed::

    python tools/sympy_coverage.py --corpus --json /tmp/symox-coverage.json

API names say nothing about signatures, domains or behavior. The expression
corpus comes from symox's tests and is intentionally biased toward implemented
features. Neither percentage estimates the fraction of SymPy implemented.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import inspect
import json
from pathlib import Path
import runpy
import subprocess


ROOT = Path(__file__).resolve().parents[1]


def surface(reference, actual, names):
    names = sorted(set(names))
    missing = [name for name in names if not hasattr(actual, name)]
    return {
        "present": len(names) - len(missing),
        "total": len(names),
        "presence_percent": round(100 * (len(names) - len(missing)) / len(names), 2) if names else None,
        "missing": missing,
    }


def methods(value):
    for name in dir(type(value)):
        if not name.startswith('_'):
            member = inspect.getattr_static(value, name)
            if callable(member) or isinstance(member, (classmethod, staticmethod)):
                yield name


def converted(value):
    if hasattr(value, 'to_sympy'):
        return value.to_sympy()
    if isinstance(value, (list, tuple)):
        return type(value)(converted(item) for item in value)
    if isinstance(value, dict):
        return {converted(key): converted(item) for key, item in value.items()}
    return value


def corpus_report(sp, sx):
    suite = runpy.run_path(str(ROOT / 'tests' / 'compare_sympy.py'))
    reference, native = suite['namespace'](sp), suite['namespace'](sx)
    results, details = Counter(), []
    for case in suite['CASES']:
        try:
            expected = eval(case, reference)
        except Exception as error:
            status, detail = 'reference_error', f'{type(error).__name__}: {error}'
        else:
            try:
                actual = eval(case, native)
                back = converted(actual)
                if expected == back:
                    status, detail = 'same_tree', None
                elif isinstance(expected, sp.Expr) and isinstance(back, sp.Expr) and sp.simplify(expected-back) == 0:
                    status, detail = 'proved_equivalent', None
                else:
                    # Do not count numerical agreement or matching text as proof.
                    status, detail = 'unverified_difference', {'symox': str(actual), 'sympy': str(expected)}
            except Exception as error:
                status, detail = 'native_error', f'{type(error).__name__}: {error}'
        results[status] += 1
        if detail is not None:
            details.append({'case': case, 'status': status, 'detail': detail,
                            'known_deviation': case in suite['KNOWN_DEVIATIONS']})
    total = sum(results.values()) - results['reference_error']
    confirmed = results['same_tree'] + results['proved_equivalent']
    return {'total': total, 'confirmed': confirmed,
            'confirmed_percent': round(100*confirmed/total, 2) if total else None,
            'counts': dict(results), 'details': details}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--corpus', action='store_true', help='also run exact differential expression checks')
    parser.add_argument('--json', type=Path, help='write the complete report, including missing API names')
    options = parser.parse_args()
    import sympy as sp
    import symox as sx

    x, X = sx.Symbol('x'), sp.Symbol('x')
    samples = {
        'Expr methods': (X, x),
        'Poly methods': (sp.Poly(X*X+1, X), sx.Poly(x*x+1, x)),
        'Matrix methods': (sp.Matrix([[1, 2], [3, 4]]), sx.Matrix([[1, 2], [3, 4]])),
        'MatrixExpr methods': (sp.MatrixSymbol('A', 2, 2), sx.MatrixSymbol('A', 2, 2)),
        'AlgebraicNumber methods': (sp.AlgebraicNumber(sp.sqrt(2)), sx.AlgebraicNumber(sx.sqrt(2))),
        'Rational domain methods': (sp.QQ, sx.QQ),
    }
    report = {
        'measured_at': datetime.now(timezone.utc).isoformat(),
        'sympy_version': sp.__version__, 'symox_version': sx.__version__,
        'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'source_dirty': bool(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()),
        'limitations': [
            'API presence does not verify signatures, domains or behavior.',
            'The expression corpus is curated for symox; its percentage is not global SymPy coverage.',
            'There is no combined global parity percentage.',
        ],
        'surfaces': {'SymPy top-level exports': surface(sp, sx, sp.__all__)},
    }
    for name, (reference, actual) in samples.items():
        report['surfaces'][name] = surface(reference, actual, methods(reference))
    if options.corpus:
        report['corpus'] = corpus_report(sp, sx)
    for name, values in report['surfaces'].items():
        print(f"{name}: {values['present']}/{values['total']} names present ({values['presence_percent']}%)")
    if options.corpus:
        values = report['corpus']
        print(f"Curated corpus: {values['confirmed']}/{values['total']} confirmed ({values['confirmed_percent']}%)")
        print('Corpus counts:', json.dumps(values['counts'], sort_keys=True))
    print('These are separate indicators, not a percentage of global SymPy parity.')
    if options.json:
        options.json.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
