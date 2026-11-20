"""Small differential contracts; invoke once on each frozen extension."""
import argparse
import importlib.util
import json
from pathlib import Path

root = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location('polytools_wave47', root / 'tests/test_polytools_wave47.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
parser = argparse.ArgumentParser()
parser.add_argument('output', type=Path)
args = parser.parse_args()
rows = [dict(name=n, reference=module.run(module.sympy, n), actual=module.run(module.symox, n))
        for n in sorted(module.cases(module.sympy))]
args.output.write_text(json.dumps(rows, indent=2) + '\n')
print(json.dumps(dict(total=len(rows), same=sum(r['reference'] == r['actual'] for r in rows))))
for row in rows:
    if row['reference'] != row['actual']:
        print(row)
