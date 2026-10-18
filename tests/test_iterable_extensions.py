"""Container algorithms: exact small oracle comparisons and structural checks."""
import copy
import itertools
from collections import Counter
import pytest
from symox.utilities import iterables as sx
from sympy.utilities import iterables as sp


def snapshot(iterator):
    return [copy.deepcopy(v) for v in iterator]


CASES=[
    ('ibin',(10,)),('ibin',(10,8)),('common_prefix',([1,2,3],[1,2,4])),
    ('common_suffix',('abc','zbc')),('rotate_left',([1,2,3],7)),
    ('rotate_right',('abcd',-2)),('least_rotation',('baca',)),
    ('minlex',('cbab',)),('runs',([1,2,2,1,4],)),
    ('roundrobin',('abc','de','f')),('uniq',([1,2,1,3,2],)),
    ('sequence_partitions',('abcd',2)),('sequence_partitions_empty',('abc',3)),
    ('generate_derangements',('aabc',)),('multiset_derangements',('aabc',)),
    ('multiset_permutations',('aabc',)),('multiset_combinations',('aabc',2)),
    ('permute_signs',((0,1,2),)),('signed_permutations',((1,2),)),
    ('filter_symbols',([1,2,3,1],{2})),
]
CASES += [(name,(n,)) for name in ('partitions','ordered_partitions','binary_partitions','generate_involutions','generate_oriented_forest') for n in range(1,7)]
CASES += [('generate_bell',(n,)) for n in range(1,6)]
CASES += [(name,(n,2)) for name in ('necklaces','bracelets') for n in range(6)]
CASES += [('multiset_partitions',(v,k)) for v in [1,4,[1,1,2,2],[1,1,1],['a','b','a']] for k in [None,1,2,4]]
CASES += [('kbins',([1,2,3,4],2,flag)) for flag in [None,0,1,10,11]]
CASES += [('ordered_partitions',(7,k)) for k in range(1,5)]
CASES += [('partitions',(7,m,k,True)) for m in [None,2,4] for k in [None,2,4]]


@pytest.mark.parametrize('name,args',CASES)
def test_reference_order(name,args):
    actual=getattr(sx,name)(*copy.deepcopy(args))
    expected=getattr(sp,name)(*copy.deepcopy(args))
    if isinstance(actual,(list,tuple,str,int,dict)):
        assert actual==expected
    else:
        assert snapshot(actual)==snapshot(expected)


def test_fair_infinite_product():
    pairs=list(itertools.islice(sx.iproduct(itertools.count(),itertools.count()),100))
    assert (3,4) in pairs and (4,3) in pairs
    assert len(pairs)==len(set(pairs))
    assert list(sx.iproduct([],itertools.count()))==[]
    assert list(sx.iproduct())==[()]
    assert sorted(sx.iproduct([1,2],[3],[4,5]))==list(itertools.product([1,2],[3],[4,5]))


def test_graph_components():
    graph=([0,1,2,3,4],[(0,1),(1,0),(1,2),(2,3),(3,2)])
    assert sx.strongly_connected_components(graph)==sp.strongly_connected_components(graph)
    assert sx.connected_components(graph)==sp.connected_components(graph)
    graph=(list(range(3000)),[(i,i+1) for i in range(2999)])
    assert len(sx.strongly_connected_components(graph))==3000


@pytest.mark.parametrize('values',[[],[1,2],list(range(20)),list('aabbcc'),list('aaaab')])
def test_bounded_random_derangement(values):
    for pick in [lambda seq:seq[0],lambda seq:seq[-1]]:
        possible=not values or 2*max(Counter(values).values())<=len(values)
        if not possible:
            with pytest.raises(ValueError):sx.random_derangement(values,choice=pick)
        result=sx.random_derangement(values,choice=pick,strict=False)
        assert Counter(result)==Counter(values)
        fixed=sum(a==b for a,b in zip(values,result))
        assert fixed==max(0,2*max(Counter(values).values(),default=0)-len(values))


def test_symbol_generators_and_wrappers(monkeypatch):
    import symox as s
    x=s.Symbol('x0')
    symbols=sx.numbered_symbols(exclude=[x],real=True)
    values=list(itertools.islice(symbols,3))
    assert [v.name for v in values]==['x0','x1','x2'] # assumptions distinguish x0
    assert all(v.is_real is True for v in values)
    assert sx.capture(lambda:print('native'))=='native\n'
    assert list(sx.postorder_traversal(s.sin(x)))==[x,s.sin(x)]
    assert sx.default_sort_key(x)==s.default_sort_key(x)
    monkeypatch.setattr('builtins.input',lambda prompt:'f')
    assert sx.interactive_traversal(s.sin(x))==x
    assert sx.iterable(iter([1])) and not sx.is_sequence(iter([1]))
    assert sx.is_sequence([1]) and not sx.iterable('abc')
    assert not sx.iterable(sx.NotIterable())
    assert sx.is_palindromic('abcdcba')


def test_empty_multiset_partition():
    assert list(sx.multiset_partitions(0))==[[]]


def test_involutions_without_factorial_filter():
    assert sum(1 for _ in sx.generate_involutions(10))==9496
    assert list(sx.generate_involutions(0))==[()]
