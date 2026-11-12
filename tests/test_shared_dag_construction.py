"""Constructors on deep shared expression graphs stay polynomial in depth.

``e_{k+1} = exp((h(e_k) - e_k)/y) + sin(e_k) + x`` reuses ``e_k`` several
times per level: the expanded tree grows exponentially while the graph grows
linearly. ``Min``, ``Abs`` and relationals on such graphs used to walk the
expanded tree in a sign query (``reference_numeric_sign``).
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
    }[kind]
    e = x
    start = time.perf_counter()
    for _ in range(depth):
        m = make(e)
        e = sx.exp((m - e) / y) + sx.sin(e) + x
    print(time.perf_counter() - start)
''')


@pytest.mark.parametrize('kind', ['Min', 'Max', 'Abs', 'Lt'])
def test_constructors_on_deep_shared_graphs_are_not_exponential(kind):
    # Before the fix, depth 14 already took minutes (x5 to x25 per level);
    # depth 40 would never finish. The subprocess bounds a regression.
    result = subprocess.run([sys.executable, '-c', SCRIPT, kind, '40'],
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr
    assert float(result.stdout) < 60, result.stdout
