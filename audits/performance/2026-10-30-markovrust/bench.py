"""Markov chain public-call benchmark.

usage: python bench.py sympy|symox CALL [repeats]
Prints one JSON line: first-call time and the median of the repeated calls
(seconds), plus the srepr of the result (for the cross-check)."""
import importlib, json, random, statistics, sys, time, warnings

lib = importlib.import_module(sys.argv[1])
st = importlib.import_module(sys.argv[1] + '.stats')
R, Eq = lib.Rational, lib.Eq


def rational_chain(n, seed, density=0.6, absorbing=()):
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        if i in absorbing:
            rows.append([1 if j == i else 0 for j in range(n)])
            continue
        w = [rng.randint(1, 9) if rng.random() < density else 0 for _ in range(n)]
        if sum(w) == 0:
            w[(i + 1) % n] = 1
        rows.append([R(x, sum(w)) for x in w])
    return lib.Matrix(rows)


def generator(n, seed):
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        r = [R(rng.randint(0, 5), rng.randint(1, 3)) if j != i else 0 for j in range(n)]
        r[i] = -sum(r)
        rows.append(r)
    return lib.Matrix(rows)


D8 = st.DiscreteMarkovChain('X', list(range(8)), rational_chain(8, 1, 1.0))
A8 = st.DiscreteMarkovChain('Y', list(range(8)), rational_chain(8, 2, 0.5, absorbing=(0, 5)))
C5 = st.ContinuousMarkovChain('C', list(range(5)), generator(5, 3))
D6 = st.DiscreteMarkovChain('Z', list(range(6)), rational_chain(6, 4, 0.7))
X, Y, Z = D8, A8, D6
CALLS = {
    'classes8': lambda: Y.communication_classes(),
    'decompose8': lambda: Y.decompose(),
    'canonical8': lambda: Y.canonical_form(),
    'stationary8': lambda: X.stationary_distribution(),
    'fundamental8': lambda: Y.fundamental_matrix(),
    'absorbing8': lambda: Y.absorbing_probabilities(),
    'regular8': lambda: X.is_regular(),
    'ctmc_limiting5': lambda: C5.limiting_distribution(),
    'P_path6x4': lambda: st.P(Eq(Z[1], 0) & Eq(Z[2], 1) & Eq(Z[3], 2) & Eq(Z[4], 3), Eq(Z[0], 0)),
    'E_path8x3': lambda: st.E(X[1] + X[2] + X[3], Eq(X[0], 0)),
    'P_path8x4': lambda: st.P(lib.Ne(X[2], 0) & lib.Ne(X[3], 1) & lib.Ne(X[4], 2) & lib.Ne(X[5], 3), Eq(X[0], 0)),
}
name = sys.argv[2]
repeats = int(sys.argv[3]) if len(sys.argv) > 3 else 5
f = CALLS[name]
warnings.simplefilter('ignore')
t0 = time.perf_counter()
try:
    r = f()
    out = lib.srepr(r) if isinstance(r, lib.Basic) else repr(r)
except Exception as e:  # noqa
    out = 'ERR:%s:%s' % (type(e).__name__, e)
first = time.perf_counter() - t0
times = []
for _ in range(repeats):
    t0 = time.perf_counter()
    try:
        f()
    except Exception:  # noqa
        pass
    times.append(time.perf_counter() - t0)
print(json.dumps({'call': name, 'lib': sys.argv[1], 'first': first, 'median': statistics.median(times),
                  'result': out[:20000]}))
