"""Freeze the small public contracts; use --baseline only on wave45."""
import argparse
import importlib.util
import json
from pathlib import Path

root = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('polytools_wave46', root / 'tests/test_polytools_wave46.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
parser = argparse.ArgumentParser()
parser.add_argument('output', type=Path)
parser.add_argument('--baseline', action='store_true')
parser.add_argument('--initial-only', action='store_true', help='Only the original 91 contracts')
parser.add_argument('--generators-only', action='store_true', help='Only the additional generator/domain regressions')
args = parser.parse_args()
if args.baseline:
    # Admission equivalent to the pre-wave46 public wrappers. The baseline
    # extension has no gcd_list/lcm_list exports; its gcd/lcm do the work.
    for name, op in [('gcd_list', 'gcd'), ('lcm_list', 'lcm')]:
        setattr(module.symox._core, name, lambda seq, *g, _op=op, **k: getattr(module.symox, _op)(list(seq), *g, **k))
rows = [dict(name=n, reference=module.run(module.sympy, n), actual=module.run(module.symox, n))
        for n in sorted(module.cases(module.sympy)) if (not args.generators_only or 'Poly generators ' in n)
        and (not args.initial_only or 'Poly generators ' not in n)]
args.output.write_text(json.dumps(rows, indent=2) + '\n')
print(json.dumps(dict(total=len(rows), same=sum(r['reference'] == r['actual'] for r in rows))))
for row in rows:
    if row['reference'] != row['actual']:
        print(row)
