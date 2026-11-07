"""The algebraic composite cases of SymPy 1.14's
`polys/domains/tests/test_domains.py` (`test_Domain_unify`,
`test_Domain_convert`), run live on one library.

Usage: python upstream.py sympy|symox OUT.jsonl

The upstream file cannot be imported against symox as a whole (it needs
`sympy.external.gmpy`, `polys.rings`, `polys.fields`...), so the loops are
transcribed: for every pair of domains of `test_Domain_convert` that
involves an algebraic ground, the record keeps `str(K3.unify(K4))` and
whether the four conversions of one and zero of the upstream
`check_domains` hold (or the exception class). The `test_Domain_unify`
assertions that involve `ALG` are recorded the same way."""
import json, sys

mode, out = sys.argv[1], sys.argv[2]
M = __import__(mode)
sqrt, I, symbols = M.sqrt, M.I, M.symbols
ZZ, QQ, RR, CC, EX = M.ZZ, M.QQ, M.RR, M.CC, M.EX
ZZ_I, QQ_I = M.ZZ_I, M.QQ_I
x, y, z = symbols('x y z')
ALG = QQ.algebraic_field(sqrt(2), sqrt(3))
QQ2 = QQ.algebraic_field(sqrt(2))
QQ3 = QQ.algebraic_field(sqrt(3))
records = []


def record(case, fn):
    try:
        r = fn()
        res = str(r)
    except Exception as e:
        res = 'error ' + type(e).__name__
    records.append({'case': case, 'result': res})


for case, fn in [
    ('unify(ZZ[x], ALG)', lambda: ZZ[x].unify(ALG)),
    ('unify(ZZ.frac_field(x), ALG)', lambda: ZZ.frac_field(x).unify(ALG)),
    ('unify(ALG, ZZ[x])', lambda: ALG.unify(ZZ[x])),
    ('unify(ALG[x], ZZ_I)', lambda: ALG[x].unify(ZZ_I)),
    ('unify(EX, ALG[x])', lambda: EX.unify(ALG[x])),
    ('unify(ALG[x], RR)', lambda: ALG[x].unify(RR)),
    ('ALG[x] == unify', lambda: ZZ[x].unify(ALG) == ALG[x]),
    ('ALG.frac_field(x) == unify', lambda: ZZ.frac_field(x).unify(ALG) == ALG.frac_field(x)),
]:
    record(case, fn)


def composite_domains(K):
    return [K, K[y], K[z], K[y, z], K.frac_field(y), K.frac_field(z), K.frac_field(y, z)]


def check(K3, K4):
    K5 = K3.unify(K4)
    ok = (K5.convert_from(K3.one, K3) == K5.one and K5.convert_from(K4.one, K4) == K5.one
          and K5.convert_from(K3.zero, K3) == K5.zero and K5.convert_from(K4.zero, K4) == K5.zero)
    return '%s %s' % (K5, ok)


doms = [('ZZ', ZZ), ('QQ', QQ), ('QQ2', QQ2), ('QQ3', QQ3), ('QQ_I', QQ_I), ('ZZ_I', ZZ_I), ('RR', RR), ('CC', CC)]
for i, (n1, K1) in enumerate(doms):
    for n2, K2 in doms[i:]:
        if 'QQ2' not in (n1, n2) and 'QQ3' not in (n1, n2):
            continue
        for a, K3 in enumerate(composite_domains(K1)):
            for b, K4 in enumerate(composite_domains(K2)):
                if a == 0 and b == 0:
                    continue  # no composite domain involved
                record('convert %s/%d %s/%d' % (n1, a, n2, b), lambda: check(K3, K4))

with open(out, 'w') as fh:
    for r in records:
        fh.write(json.dumps(r) + '\n')
print(len(records), 'records')
