"""Independent PLDU identity checks beyond the frozen small-shape corpus."""
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path
from random import Random
import sys

stage = Path('/home/keckj/.cache/symox-validation/native-ddm-fflu-r11-20261028/final-work')
sys.path.insert(0, str(stage / 'python'))
import symox as sx
from symox.polys.matrices.ddm import DDM
import symox.polys.matrices._fflu_native as admission

rng = Random(20261028)
names = ('_ddm_fflu', '_ddm_fflu_p', '_ddm_fflu_l', '_ddm_fflu_d')
calls = {name: 0 for name in names}
for name in names:
    original = getattr(sx._core, name)
    def counted(*args, _name=name, _original=original):
        result = _original(*args)
        assert result is True, _name
        calls[_name] += 1
        return result
    setattr(sx._core, name, counted)

def rank(a, n):
    a = [row[:] for row in a]
    r = 0
    for j in range(n):
        pivot = next((i for i in range(r, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        divisor = a[r][j]
        a[r] = [v / divisor for v in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][j]
            a[i] = [v - scale*w for v, w in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r

cases = []
for index in range(320):
    m, n = rng.randrange(11), rng.randrange(13)
    domain = sx.QQ if index % 2 else sx.ZZ
    a = [[F(rng.randrange(-9, 10), rng.randrange(1, 12) if index % 2 else 1)
          if rng.randrange(4) else F(0) for _ in range(n)] for _ in range(m)]
    if index % 3 == 0:
        leading = rng.randrange(n + 1)
        for row in a:
            row[:leading] = [F(0)] * leading
    if index % 5 == 0 and m > 1:
        a[-1] = [F(3)*v for v in a[0]]
    expected_rank = rank(a, n)
    outputs = []
    for fallback in (False, True):
        saved = admission._READY
        if fallback:
            admission._READY = False
        try:
            source = DDM([[domain(v.numerator, v.denominator) if index % 2 else domain(v.numerator)
                           for v in row] for row in a], (m, n), domain)
            factors = source.fflu()
        finally:
            admission._READY = saved
        p, l, d, u = [[[F(str(v)) for v in row] for row in factor] for factor in factors]
        assert all(d[i][i] for i in range(m))
        assert all(sum(p[i]) == 1 for i in range(m))
        assert all(sum(p[i][j] for i in range(m)) == 1 for j in range(m))
        assert all(v in (0, 1) for row in p for v in row)
        assert all(l[i][j] == 0 for i in range(m) for j in range(i+1, m))
        assert all(d[i][j] == 0 for i in range(m) for j in range(m) if i != j)
        pivots = [next((j for j, v in enumerate(row) if v), None) for row in u]
        nonzero = [j for j in pivots if j is not None]
        assert pivots == nonzero + [None] * (m-len(nonzero))
        assert nonzero == sorted(set(nonzero))
        assert len(nonzero) == expected_rank == rank(u, n)
        for i in range(m):
            for j in range(n):
                assert sum(p[i][k]*a[k][j] for k in range(m)) == sum(l[i][k]*u[k][j]/d[k][k] for k in range(m))
        outputs.append((p,l,d,u))
    assert outputs[0] == outputs[1]
    cases.append({'index': index, 'shape': [m,n], 'domain': 'QQ' if index % 2 else 'ZZ', 'rank': expected_rank})

assert not any(name == 'sympy' or name.startswith('sympy.') for name in sys.modules)
binary = Path(sx._core.__file__)
result = {'cases': cases, 'certified_factorizations': len(cases)*2, 'native_admissions': calls,
          'native_sha256': hashlib.sha256(binary.read_bytes()).hexdigest(),
          'source': str(stage), 'reference_imported': False,
          'certificate': 'independent Fraction rank, echelon, invertible diagonal, permutation, and P*A=L*D^-1*U; exact native/fallback factors'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({key: value for key, value in result.items() if key != 'cases'}))
