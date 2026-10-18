"""Exact finite root systems: Cartan matrices, counts and reflection closure."""
import pytest
import symox as s
from symox.liealgebras import CartanType
from symox.liealgebras.root_system import RootSystem
from symox.liealgebras.cartan_matrix import CartanMatrix
from symox.liealgebras.dynkin_diagram import DynkinDiagram


@pytest.mark.parametrize('name,count,dimension',[
    ('A1',2,3),('A2',6,8),('A3',12,15),('A5',30,35),
    ('B2',8,10),('B3',18,21),('B5',50,55),('C2',8,10),('C3',18,21),
    ('D3',12,15),('D4',24,28),('D5',40,45),('E6',72,78),('E7',126,133),('E8',240,248),
    ('F4',48,52),('G2',12,14),
])
def test_root_counts_lengths_and_reflections(name,count,dimension):
    cartan = CartanType(name)
    assert cartan.roots() == count
    assert cartan.basis() == dimension
    roots = RootSystem(name).all_roots()
    values = {tuple(s.sympify(v) for v in root) for root in roots.values()}
    assert len(roots) == len(values) == count
    assert cartan.cartan_matrix().det() > 0
    for i in range(1,cartan.rank()+1):
        alpha = tuple(cartan.simple_root(i))
        norm = sum(v*v for v in alpha)
        assert alpha in values
        for root in values:
            factor = 2*sum(a*b for a,b in zip(root,alpha))/norm
            assert factor.is_Integer
            assert tuple(r-factor*a for r,a in zip(root,alpha)) in values
    assert sum(cartan.highest_root()[i]**2 for i in range(cartan.dimension())) != 0


@pytest.mark.parametrize('name',['A2','A4','B2','B4','C3','C4','D4','E6','E7','E8','F4','G2'])
def test_cartan_matrix_matches_reference(name):
    from sympy.liealgebras import CartanType as reference
    actual = CartanType(name)
    expected = reference(name)
    assert actual.cartan_matrix().to_sympy() == expected.cartan_matrix()
    assert actual.rank() == expected.rank()
    assert actual.dimension() == expected.dimension()
    assert CartanMatrix(name) == actual.cartan_matrix()
    assert DynkinDiagram(name) == actual.dynkin_diagram()
    if name[0] in 'ABCD':
        assert actual.positive_roots() == expected.positive_roots()


@pytest.mark.parametrize('value',['A0','B1','D2','E5','E9','F3','G3','H3','A2x',''])
def test_invalid_finite_types(value):
    with pytest.raises(ValueError):
        CartanType(value)


def test_constructor_and_root_system_operations():
    c = CartanType(['A',3])
    assert c == CartanType('A3') and CartanType(c) is c
    assert c.lie_algebra() == 'su(4)'
    root = RootSystem('A3')
    assert root.add_simple_roots(1,2) == [1,0,-1,0]
    assert root.add_as_roots([1,0,-1,0],[0,0,1,-1]) == [1,0,0,-1]
    assert root.root_space() == 'alpha[1] + alpha[2] + alpha[3]'
    with pytest.raises(ValueError):
        c.simple_root(0)
    with pytest.raises(ValueError):
        c.simple_root(4)
    with pytest.raises(TypeError):
        CartanType(3)
