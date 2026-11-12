"""Constructors on deep shared expression graphs stay polynomial in depth.

``e_{k+1} = exp((h(e_k) - e_k)/y) + sin(e_k) + x`` reuses ``e_k`` several
times per level: the expanded tree grows exponentially while the graph grows
linearly. ``Min``, ``Abs`` and relationals on such graphs used to walk the
expanded tree in a sign query (``reference_numeric_sign``); a nested ``Min``
(SymPy's Python path, ``_is_connected``) in ``factor_terms``.
"""
import subprocess
import sys
import textwrap

import pytest

SCRIPT = textwrap.dedent('''
    import sys, time
    import symox as sx
    x, y = sx.symbols("x y")
    kind, depth = sys.argv[1], int(sys.argv[2])
    make = {
        "Min": lambda e: sx.Min(e, y),
        "Max": lambda e: sx.Max(e, y),
        "Abs": lambda e: (e + y - sx.Abs(e - y)) / 2,
        "Lt": lambda e: (sx.StrictLessThan(e, y), e)[1],
        "MinMin": lambda e: sx.Min(sx.Min(e, y), x + y),
    }[kind]
    e = x
    start = time.perf_counter()
    for _ in range(depth):
        m = make(e)
        e = sx.exp((m - e) / y) + sx.sin(e) + x
    print(time.perf_counter() - start)
''')


@pytest.mark.parametrize('kind', ['Min', 'Max', 'Abs', 'Lt', 'MinMin'])
def test_constructors_on_deep_shared_graphs_are_not_exponential(kind):
    # Before the fix, depth 14 already took minutes (x5 to x25 per level);
    # depth 40 would never finish. The subprocess bounds a regression.
    result = subprocess.run([sys.executable, '-c', SCRIPT, kind, '40'],
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert float(result.stdout) < 60, result.stdout


SOFTDTW = textwrap.dedent('''
    import sys, time
    import symox as sx
    n = int(sys.argv[1])
    X, Y, g = sx.symbols(f"x0:{n}"), sx.symbols(f"y0:{n}"), sx.Symbol("g")
    d2 = lambda i, j: (X[i] - Y[j]) ** 2
    start = time.perf_counter()
    prev, acc = [], None
    for i in range(n):
        acc = d2(i, 0) if i == 0 else d2(i, 0) + acc
        prev.append(acc)
    for j in range(1, n):
        cur = [d2(0, j) + prev[0]]
        for i in range(1, n):
            a, b, c = prev[i - 1], cur[i - 1], prev[i]
            mu = sx.Min(sx.Min(a, c), b) if sys.argv[2] == "nested" else sx.Min(a, c, b)
            soft = sx.log(sx.exp((mu - a) / g) + sx.exp((mu - b) / g) + sx.exp((mu - c) / g))
            cur.append(d2(i, j) + mu - g * soft)
        prev = cur
    print(time.perf_counter() - start)
''')


@pytest.mark.parametrize('form', ['nested', 'flat'])
def test_soft_dtw_grid_is_not_exponential(form):
    # Soft-DTW recursion unrolled on an n x n grid (KeOps SoftDTW_SqDist):
    # each cell reuses three neighbours. Before the fix, 4 x 5 took 48 s and
    # every extra cell multiplied the time. ``nested`` takes SymPy's Python
    # path (``_is_connected``, ``factor_terms``, ``_nodes``), ``flat`` the
    # native kernel (``ordered`` sort keys).
    result = subprocess.run([sys.executable, '-c', SOFTDTW, '8', form],
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert float(result.stdout) < 60, result.stdout
