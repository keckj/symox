from pathlib import Path
import json,random,sympy as S
root=Path(__file__).resolve().parent
cases={}
for n,kind in [(4,'ZZ'),(8,'ZZ'),(12,'ZZ'),(4,'QQ'),(8,'QQ')]:
 pairs=[]
 for variant in range(3):
  rng=random.Random(673210+n*13+variant*31+(1000 if kind=='QQ' else 0))
  rows=[[S.Rational(rng.randrange(-6,7),rng.randrange(1,8) if kind=='QQ' else 1) for _ in range(n)] for _ in range(n)]
  coefficients=S.Matrix(rows).charpoly().all_coeffs()
  # Independent exact Cayley-Hamilton certificate, beyond output matching.
  value=S.zeros(n)
  for c in coefficients:value=value*S.Matrix(rows)+c*S.eye(n)
  assert value==S.zeros(n)
  pairs.append({'matrix':[[[str(c.p),str(c.q)] for c in row] for row in rows],'coefficients':[[str(c.p),str(c.q)] for c in coefficients],'cayley_hamilton_residual':'zero matrix'})
 cases[f'{kind}_{n}']={'domain':kind,'n':n,'inputs':pairs}
(root/'fixtures.json').write_text(json.dumps({'reference_version':S.__version__,'cases':cases},indent=2)+'\n')
print(len(cases),'families,',sum(len(v['inputs']) for v in cases.values()),'independent characteristic-polynomial and Cayley-Hamilton certificates')
