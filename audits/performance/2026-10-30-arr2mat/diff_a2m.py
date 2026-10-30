"""Differential corpus for convert_array_to_matrix, its passes and the
_EditArrayContraction editor: `diff_a2m.py <sympy|symox> <out.jsonl> [n]`.

The cases are generated from a fixed seed with the same constructor calls in
both modes; each line records str(result) or the exception type."""
import importlib
import json
import random
import sys

mode, out = sys.argv[1], sys.argv[2]
NRAND = int(sys.argv[3]) if len(sys.argv) > 3 else 400
lib = importlib.import_module(mode)
ae = importlib.import_module(mode + '.tensor.array.expressions.array_expressions')
a2m = importlib.import_module(mode + '.tensor.array.expressions.from_array_to_matrix')
m2a = importlib.import_module(mode + '.tensor.array.expressions.from_matrix_to_array')
cm, ca = m2a.convert_matrix_to_array, a2m.convert_array_to_matrix
tp, con, dia, pd, add = ae._array_tensor_product, ae._array_contraction, ae._array_diagonal, ae._permute_dims, ae._array_add
S = lib
FALLBACKS = []
if mode == 'symox':
    # Count the inputs on which the native port declines and the stock
    # Python algorithm (symox._arrayexpr) runs instead.
    _st = importlib.import_module('symox.tensor.array.expressions._stock')
    _orig_to_stock = _st.to_stock

    def _counting_to_stock(e, _top=[True]):
        if _top[0]:
            FALLBACKS.append(1)
        _top[0] = False
        try:
            return _orig_to_stock(e)
        finally:
            _top[0] = True
    _st.to_stock = _counting_to_stock
k = S.Symbol('k')
M, N, P = [S.MatrixSymbol(s, k, k) for s in 'MNP']
a, b = S.MatrixSymbol('a', k, 1), S.MatrixSymbol('b', k, 1)
c = S.MatrixSymbol('c', 1, k)
X, Y = S.MatrixSymbol('X', 3, 3), S.MatrixSymbol('Y', 3, 3)
u = S.MatrixSymbol('u', 3, 1)
I = S.Identity(k)
I3 = S.Identity(3)
x, y = S.symbols('x y')
A = ae.ArraySymbol('A', (k, k))
OK = ae.OneArray(k)
FACTORS = [M, N, P, a, b, c, I, X, Y, u, I3, A, OK, ae.OneArray(1, k)]


def dims(e):
    return list(e.shape)


def rand_array(rng):
    nf = rng.randint(1, 4)
    fs = [rng.choice(FACTORS) for _ in range(nf)]
    if rng.random() < 0.2:
        fs.insert(0, rng.choice([x, 2, y]))
    expr = tp(*fs) if len(fs) > 1 else fs[0]
    shape = dims(expr)
    axes = list(range(len(shape)))
    rng.shuffle(axes)
    groups, used = [], set()
    for i in axes:
        if i in used or rng.random() < 0.35:
            continue
        partners = [j for j in axes if j != i and j not in used and shape[j] == shape[i]]
        if partners:
            j = rng.choice(partners)
            g = [i, j]
            if rng.random() < 0.1:
                more = [q for q in partners if q != j]
                if more:
                    g.append(more[0])
            used.update(g)
            groups.append(tuple(sorted(g)))
    kind = rng.random()
    if groups:
        if kind < 0.6:
            expr = con(expr, *groups)
        elif kind < 0.85:
            expr = dia(expr, *groups)
        else:
            expr = con(dia(expr, groups[0]), *[]) if len(groups) == 1 else con(dia(expr, groups[0]), *[])
    r = len(expr.shape)
    if r > 1 and rng.random() < 0.3:
        p = list(range(r))
        rng.shuffle(p)
        expr = pd(expr, p)
    if rng.random() < 0.15:
        try:
            expr = add(expr, expr)
        except Exception:
            pass
    return expr


MATS = [M, N, P, I]


def rand_mat(rng, depth=0):
    r = rng.random()
    if depth > 2 or r < 0.3:
        return rng.choice(MATS)
    if r < 0.5:
        return S.MatMul(*[rand_mat(rng, depth + 1) for _ in range(rng.randint(2, 3))])
    if r < 0.6:
        return rand_mat(rng, depth + 1).T
    if r < 0.7:
        return S.MatAdd(rand_mat(rng, depth + 1), rand_mat(rng, depth + 1))
    if r < 0.8:
        return S.HadamardProduct(rand_mat(rng, depth + 1), rand_mat(rng, depth + 1))
    if r < 0.9:
        return S.Trace(rand_mat(rng, depth + 1)) * rand_mat(rng, depth + 1)
    return a * b.T * rand_mat(rng, depth + 1) * (c * rand_mat(rng, depth + 1) * a)


def fmt(v):
    if isinstance(v, (list, tuple)):
        return '[' + ', '.join(fmt(i) for i in v) + ']'
    return str(v)


def run(fn):
    # ' [stock]' marks a result (value or exception) computed by the stock
    # fallback.
    n0 = len(FALLBACKS)
    try:
        r = fmt(fn())
    except RecursionError:
        r = 'EXC:RecursionError'
    except Exception as e:
        r = 'EXC:' + type(e).__name__
    return r + (' [stock]' if len(FALLBACKS) > n0 else '')


PASSES = [
    ('convert', lambda e: ca(e)),
    ('a2m', lambda e: a2m._array2matrix(e)),
    ('rtd', lambda e: a2m._remove_trivial_dims(e)),
    ('hadamard', lambda e: a2m.identify_hadamard_products(e)),
    ('removable', lambda e: a2m.identify_removable_identity_matrices(e)),
    ('remove_id', lambda e: a2m.remove_identity_matrices(e)),
    ('multi_id', lambda e: a2m._array_contraction_to_diagonal_multiple_identity(e)),
    ('diag2contr', lambda e: a2m._array_diag2contr_diagmatrix(e)),
]


def editor(e):
    ed = ae._EditArrayContraction(e)
    out = [ed.number_of_contraction_indices, [str(x.element) + str(x.indices) for x in ed.args_with_ind],
           ed.get_contraction_indices(), str(ed.get_contraction_indices_to_ind_rel_pos()),
           ed.number_of_diagonal_indices]
    ed.merge_scalars()
    out.append(str(ed.to_array_contraction()))
    return out


rng = random.Random(20261030)
cases = []
for i in range(NRAND):
    try:
        cases.append(('arr%d' % i, rand_array(rng)))
    except Exception as e:
        cases.append(('arr%d' % i, 'EXC:' + type(e).__name__))
for i in range(NRAND // 2):
    try:
        cases.append(('mat%d' % i, cm(rand_mat(rng))))
    except Exception as e:
        cases.append(('mat%d' % i, 'EXC:' + type(e).__name__))

with open(out, 'w') as fh:
    for name, e in cases:
        rec = {'case': name, 'input': str(e)}
        if isinstance(e, str):
            fh.write(json.dumps(rec) + '\n')
            continue
        for pname, fn in PASSES:
            rec[pname] = run(lambda: fn(e))
        if isinstance(e, (ae.ArrayContraction, ae.ArrayDiagonal, ae.ArrayTensorProduct)):
            rec['editor'] = run(lambda: editor(e))
        fh.write(json.dumps(rec) + '\n')
