"""data/cases4.json -> data/scene4.json: TeX of the symox results and plot samples.

Run with the venv python (needs SymPy); the ManimGX environment has none.
"""
import json
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from sympy.physics.hydrogen import R_nl

HERE = Path(__file__).parent / "data"
C = json.loads((HERE / "cases4.json").read_text())
x, r = sp.symbols("x r")

STATS = ["normal_moment", "normal_var", "binom_var", "binom_skew", "dice", "exp_P",
         "gamma_mgf", "sum_normals", "beta_E", "normal_entropy"]
tex = {}
for case in STATS:
    e = sp.sympify(C[case]["symox"]["result"], locals={"lam": sp.Symbol("lambda", positive=True)})
    tex[case] = sp.latex(e)

defl = sp.sympify(C["beam"]["symox"]["result"]).rewrite(sp.Piecewise)
f = sp.lambdify(x, defl, "numpy")
xs = np.linspace(0, 10, 121)
beam = [float(f(v)) for v in xs]

rho = sp.lambdify(r, (r*R_nl(5, 2, r, 1))**2, "numpy")
rs = np.linspace(0, 70, 281)
hydrogen = [float(rho(v)) for v in rs]

dice = {}
for roll in product(range(1, 7), repeat=3):
    dice[sum(roll)] = dice.get(sum(roll), 0) + 1

out = {"stats_tex": tex, "beam": {"x": xs.tolist(), "y": beam}, "hydrogen": {"r": rs.tolist(), "rho": hydrogen},
       "dice": dice}
(HERE / "scene4.json").write_text(json.dumps(out, indent=1))
for k, v in tex.items():
    print(k, v)
print("beam min", min(beam), "at", xs[int(np.argmin(beam))])

# Typeset the verified current Taylor result rather than retaining an old
# typesetting artifact when the input measurements are refreshed.
cases = json.loads((HERE / 'cases3.json').read_text())
polynomial = sp.sympify(cases['symox']['taylor_exp']['result']).removeO()
terms = [polynomial.coeff(x, k)*x**k for k in range(25) if polynomial.coeff(x, k) != 0]
series_lines = []
for offset in range(0, len(terms), 4):
    line = ' + '.join(sp.latex(term) for term in terms[offset:offset+4]).replace('+ -', '- ')
    if offset and not line.startswith('-'):
        line = '+ ' + line
    series_lines.append(line)
series_lines[-1] += r' + O(x^{25})'
(HERE / 'series_tex.json').write_text(json.dumps(series_lines, indent=1)+'\n')
