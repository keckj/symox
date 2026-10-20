"""Host instances of Symbol and Function subclasses live exactly as long as needed.

A Python instance bound to a native node keeps its identity and attributes
for as long as an expression references it, whatever the native algorithms
rebuild, and is collected (cycles included) once nothing references it.
"""
import gc
import pickle
import weakref

import symox as sp


x, y = sp.symbols('x y')


class Tagged(sp.Function):
    def __new__(cls, *args, tag=None):
        obj = super().__new__(cls, *args)
        obj.tag = tag
        return obj


class Stamped(sp.Function):
    @classmethod
    def eval(cls, arg):
        return None

    def __new__(cls, *args):
        obj = super().__new__(cls, *args)
        obj.stamp = str(args)
        return obj


class Colored(sp.Symbol):
    pass


class Picklable(sp.Symbol):
    def __new__(cls, name, color=None):
        obj = super().__new__(cls, name)
        obj.color = color
        return obj

    def __getnewargs_ex__(self):
        return (self.name, self.color), {}


def _churn(count=3000):
    for i in range(count):
        junk = [Tagged(sp.Integer(i), tag=i) + x, Colored('c%d' % (i % 7))]
        del junk
    gc.collect()


def _find(expr, cls):
    return [node for node in sp.preorder_traversal(expr) if isinstance(node, cls)]


def test_attributes_survive_while_referenced():
    a = Colored('a')
    a.color = 'red'
    expr = (Tagged(x, tag='t1') + a) * 2
    del a
    _churn()
    assert _find(expr, Tagged)[0].tag == 't1'
    assert [s for s in expr.free_symbols if isinstance(s, Colored)][0].color == 'red'


def test_identity_is_kept_through_traversal():
    g = Tagged(x, tag='same')
    expr = sp.sin(g) + g ** 2
    _churn(500)
    found = _find(expr, Tagged)
    assert found and all(node is g for node in found)


def test_instances_are_freed_when_unreferenced():
    g = Tagged(x, tag='gone')
    expr = g * y + 1
    ref = weakref.ref(g)
    del g, expr
    gc.collect()
    assert ref() is None


def test_self_referencing_cycles_are_collected():
    h = Tagged(y, tag=None)
    h.tag = h + 1
    ref = weakref.ref(h)
    del h
    gc.collect()
    assert ref() is None


def test_rebuilt_applications_keep_their_attributes():
    result = (Stamped(x) + y).subs(x, 5)
    _churn()
    assert _find(result, Stamped)[0].stamp == '(5,)'
    derived = sp.diff(Stamped(x ** 2) * x, x)
    _churn(500)
    assert all(hasattr(node, 'stamp') for node in _find(derived, Stamped))
    replaced = (Stamped(x) * Stamped(y)).xreplace({y: x + 1})
    assert sorted(node.stamp for node in _find(replaced, Stamped)) == ['(x + 1,)', '(x,)']


def test_simplify_and_expand_keep_attributes():
    g = Tagged(x, tag='kept')
    expr = sp.simplify(g * (sp.sin(x) ** 2 + sp.cos(x) ** 2))
    expanded = sp.expand((g + 1) ** 3)
    _churn(1000)
    assert _find(expr, Tagged)[0].tag == 'kept'
    assert {node.tag for node in _find(expanded, Tagged)} == {'kept'}


def test_matrix_entries_keep_attributes():
    M = sp.Matrix([[Tagged(x, tag='m'), 1], [0, Tagged(y, tag='n')]])
    product = M * M
    _churn()
    assert M[0, 0].tag == 'm'
    assert {node.tag for node in _find(product[0, 0], Tagged)} == {'m'}


def test_deep_expressions_are_collected():
    expr = Tagged(x, tag=0)
    refs = []
    for i in range(400):
        node = Tagged(expr, tag=i)
        refs.append(weakref.ref(node))
        expr = node + i
    _churn(100)
    assert _find(expr, Tagged)[0].tag == 399
    del expr, node
    gc.collect()
    assert all(ref() is None for ref in refs)


def test_pickle_round_trip_keeps_subclass_state():
    p = Picklable('p', color='blue')
    restored = pickle.loads(pickle.dumps(p * x + 1))
    _churn(500)
    symbol = [s for s in restored.free_symbols if isinstance(s, Picklable)][0]
    assert symbol.color == 'blue'
    assert restored == p * x + 1


def test_retained_instances_stay_bounded():
    # Native results may sit in bounded caches (as SymPy's cacheit does),
    # but the number of instances kept alive does not grow with the work.
    refs = []

    def run(start, count):
        for i in range(start, start + count):
            g = Tagged(sp.Integer(i), tag=[i] * 10)
            refs.append(weakref.ref(g))
            sp.diff((g * x + y).subs(x, y), y)
            del g
        gc.collect()
        return sum(ref() is not None for ref in refs)

    first = run(0, 10000)
    second = run(10000, 10000)
    assert second <= first + 64
    assert second < 20000 // 2
