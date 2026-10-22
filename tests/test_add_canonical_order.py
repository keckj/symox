"""Canonical argument order and equality of Add and Mul over host classes.

Expected values are those of SymPy 1.14. A held (``evaluate=False``) sum or
product equals a canonical one exactly when their public arguments agree in
order; the internal hash order of a canonical node never decides it, so the
result does not depend on the import or class registration order.
"""
import os
import subprocess
import sys
import textwrap

from symox import Add, Expr, Function, Mul, S, Symbol, symbols
from symox.calculus.accumulationbounds import AccumBounds as B

x, y = symbols('x y')


class _E(Expr):
    pass


class _F(Function):
    pass


def _same(*values):
    first = values[0]
    assert all(v == first for v in values)
    assert len({hash(v) for v in values}) == 1
    assert len({v.args for v in values}) == 1


def test_accumbounds_add_both_construction_orders():
    b = B(1, 2)
    _same(x + b, b + x, Add(b, x), Add(x, b))
    assert (x + b).args == (B(1, 2), x)
    assert Add(b, x).args == (B(1, 2), x)
    assert x + B(1, 2) == Add(B(1, 2), x)
    assert B(1, 2) - x == Add(B(1, 2), -x)


def test_accumbounds_mul_both_construction_orders():
    b = B(1, 2)
    _same(x*b, b*x, Mul(b, x), Mul(x, b))
    assert Mul(b, x).args == (B(1, 2), x)


def test_accumbounds_coefficient_folding():
    # SymPy Add.flatten / Mul.flatten fold Numbers and bounds into the
    # coefficient, which leads the arguments.
    b = B(1, 2)
    assert Add(b, 1, x).args == (B(2, 3), x)
    assert Add(x, b, B(3, 4)).args == (B(4, 6), x)
    assert Mul(b, 2, x).args == (B(2, 4), x)
    assert Add(x + b, 1).args == (B(2, 3), x)
    assert Mul(x*b, 3).args == (B(3, 6), x)
    assert Add(b, x, y).args == (B(1, 2), x, y)
    assert (x - b).args == (B(-2, -1), x)
    assert Add(B(1, S.Infinity), S.Infinity, x) == S.Infinity + x


def test_accumbounds_operator_chains_fold():
    b = B(1, 2)
    assert (x + b + 1).args == (B(2, 3), x)
    assert (x + b - 1).args == (B(0, 1), x)
    assert (2*(x*b)).args == (B(2, 4), x)
    assert (x*b*x).args == (B(1, 2), x**2)
    assert ((b + x)*2).args == (B(2, 4), 2*x)
    assert (x + b) + S.Infinity == x + S.Infinity
    # Held as in the reference: x + bounds is not real, so no folding.
    assert (x + b + B(3, 4)).args == (B(3, 4), x + B(1, 2))
    assert ((x*b)*B(1, 3)).args == (B(1, 3), x*B(1, 2))


def test_host_expr_and_function_subclasses():
    for o in (_E(x), _F(x)):
        _same(x + o, o + x, Add(o, x), Add(x, o))
        _same(x*o, o*x, Mul(o, x), Mul(x, o))
        assert (x + o).args == (x, o)
        assert (x*o).args == (x, o)


def test_held_equality_follows_public_arguments():
    assert (Add(y, x, evaluate=False) == x + y) is False
    assert (Add(x, y, evaluate=False) == x + y) is True
    assert (Mul(y, x, evaluate=False) == x*y) is False
    assert (Mul(x, 2, evaluate=False) == 2*x) is False
    assert hash(Add(x, y, evaluate=False)) == hash(x + y)
    a = symbols('a0:30')
    for i in range(30):
        for j in range(30):
            if i != j:
                held = Add(a[j], a[i], evaluate=False)
                canonical = a[i] + a[j]
                assert (held == canonical) == (held.args == canonical.args)
                held = Mul(a[j], a[i], evaluate=False)
                canonical = a[i]*a[j]
                assert (held == canonical) == (held.args == canonical.args)


_IMPORT_ORDER_SCRIPT = textwrap.dedent('''
    import importlib, sys
    for name in sys.argv[1:]:
        importlib.import_module(name)
    from symox import Add, Mul, Symbol, Expr, Function
    from symox.calculus.accumulationbounds import AccumBounds as B
    class E(Expr):
        pass
    class F(Function):
        pass
    x = Symbol('x')
    out = []
    for o in (B(1, 2), E(x), F(x)):
        for O, op in ((Add, lambda a, b: a + b), (Mul, lambda a, b: a*b)):
            values = [op(x, o), op(o, x), O(o, x), O(x, o)]
            out.append((str(values[0].args),
                        all(v == values[0] for v in values),
                        len({hash(v) for v in values}) == 1,
                        len({v.args for v in values}) == 1))
    print(out)
''')


def test_import_order_independent():
    root = os.path.join(os.path.dirname(__file__), os.pardir, 'python')
    env = dict(os.environ, PYTHONPATH=os.path.abspath(root))
    orders = [
        [],
        ['symox.calculus.accumulationbounds'],
        ['symox.stats', 'symox.physics.quantum', 'symox.calculus.accumulationbounds'],
        ['symox.printing.codeprinter', 'symox.tensor', 'symox.calculus.accumulationbounds',
         'symox.stats'],
    ]
    results = set()
    for order in orders:
        run = subprocess.run([sys.executable, '-c', _IMPORT_ORDER_SCRIPT, *order],
                             env=env, capture_output=True, text=True, timeout=300)
        assert run.returncode == 0, run.stderr
        results.add(run.stdout.strip())
    assert len(results) == 1
    (result,) = results
    assert 'False' not in result
    assert "'(AccumBounds(1, 2), x)'" in result
