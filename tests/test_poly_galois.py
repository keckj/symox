"""Exact Galois classification against an independent implementation."""
import pickle
import random
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx

from symox.combinatorics import Permutation, PermutationGroup
from symox.polys.numberfields import galois_group


@pytest.mark.parametrize('source,group', [
    ('2*x+3', 'S1'), ('x**2-2', 'S2'),
    ('x**3-2', 'S3'), ('x**3-3*x+1', 'A3'),
    ('x**4-2', 'D4'), ('x**4+x+1', 'S4'),
    ('x**4+8*x+12', 'A4'), ('x**4+x**3+x**2+x+1', 'C4'),
    ('x**4+1', 'V'), ('x**4-10*x**2+1', 'V'),
    ('x**4-5*x**2+5', 'C4'), ('x**4-2*x**2+2', 'D4'),
])
@pytest.mark.parametrize('rational', [False, True])
def test_galois_groups(source, group, rational):
    x = sx.Symbol('x')
    expr = eval(source, {'x': x})
    # Nonmonic/rational inputs must preserve the splitting field.
    expr = expr/3 if rational else expr
    p = sx.Poly(expr, x)
    reference = sp.Poly(p.as_expr().to_sympy(), x.to_sympy())
    expected_name, expected_alt = reference.galois_group(by_name=True)
    name, alternating = p.galois_group(by_name=True)
    assert name.name == expected_name.name == group
    assert alternating is expected_alt
    assert pickle.loads(pickle.dumps(name)) == name
    actual, alternating = p.galois_group()
    expected, expected_alt = reference.galois_group()
    assert isinstance(actual, PermutationGroup)
    assert [g.array_form for g in actual.generators] == [g.array_form for g in expected.generators]
    assert actual.order() == expected.order()
    assert actual.degree == expected.degree
    assert actual.is_transitive() == expected.is_transitive()
    assert actual.is_abelian == expected.is_abelian
    assert actual.is_solvable == expected.is_solvable
    assert alternating is expected_alt
    assert {tuple(g.array_form) for g in actual.elements} == {
        tuple(g.array_form) for g in expected.generate_schreier_sims()}
    assert pickle.loads(pickle.dumps(actual)) == actual
    assert galois_group(expr, x, by_name=True) == (name, alternating)
    assert sx.galois_group(p, by_name=True, randomize=True) == (name, alternating)


def test_seeded_irreducible_cubics_and_quartics():
    rng = random.Random(20261017)
    x = sx.Symbol('x')
    checked = 0
    for _ in range(100):
        degree = rng.choice([3, 4])
        expr = rng.choice([1, 2, 3])*x**degree + sum(rng.randint(-5, 5)*x**i for i in range(degree))
        p = sx.Poly(expr, x)
        expected = sp.Poly(expr.to_sympy(), x.to_sympy())
        assert p.is_irreducible == expected.is_irreducible
        if not p.is_irreducible:
            continue
        name, alternating = p.galois_group(by_name=True)
        expected_name, expected_alt = expected.galois_group(by_name=True)
        assert name.name == expected_name.name
        assert alternating is expected_alt
        checked += 1
    assert checked >= 40


@pytest.mark.parametrize('source', ['x**2-1', 'x+y', '0', '1', 'x**7-2'])
def test_invalid_polynomials(source):
    x, y = sx.symbols('x y')
    p = sx.Poly(eval(source, {'x': x, 'y': y}), x, y) if source == 'x+y' else sx.Poly(eval(source, {'x': x}), x)
    with pytest.raises(ValueError):
        p.galois_group()


def test_unsupported_degrees_and_domains():
    x, a = sx.symbols('x a')
    # Degrees 5 and 6 are supported (SymPy 1.14: M20 and D6).
    assert sx.Poly(x**5-2, x).galois_group(by_name=True)[0].name == 'M20'
    assert sx.Poly(x**6-2, x).galois_group(by_name=True)[0].name == 'D6'
    with pytest.raises(ValueError):
        sx.Poly(x**7-2, x).galois_group()
    for p in (sx.Poly(x**2+x+1, x, modulus=2), sx.Poly(x**2-a, x),
              sx.Poly(x**2-sx.sqrt(2), x)):
        with pytest.raises(ValueError):
            p.galois_group()
    with pytest.raises(ValueError, match='max_tries'):
        sx.Poly(x**2-2, x).galois_group(max_tries=0)


def test_permutation_arithmetic_and_membership():
    p = Permutation([1, 2, 3, 0])
    q = Permutation([1, 0, 2, 3])
    from sympy.combinatorics import Permutation as SP
    assert (p*q).array_form == (SP(p.array_form)*SP(q.array_form)).array_form
    assert (~p).array_form == (~SP(p.array_form)).array_form
    assert (p**-13).array_form == (SP(p.array_form)**-13).array_form
    assert p.order() == 4 and p.signature() == -1
    assert p.cyclic_form == [[0, 1, 2, 3]]
    group = PermutationGroup(p)
    assert p in group and q not in group
    assert p*~p == group.identity
    images = p.array_form
    images[0] = 99
    assert p(0) == 1
    for bad in ([0, 0], [-1], [1, 2]):
        with pytest.raises(ValueError):
            Permutation(bad)


def test_galois_without_sympy():
    subprocess.run([sys.executable, '-c', '''
import sys
class BlockSymPy:
    def find_spec(self, name, *args):
        if name == 'sympy' or name.startswith('sympy.'):
            raise RuntimeError('SymPy imported at runtime')
sys.meta_path.insert(0, BlockSymPy())
import symox as s
x = s.Symbol('x')
for f, name, order in [(x**3-2, 'S3', 6), (x**4-2, 'D4', 8),
                        (x**4+x**3+x*x+x+1, 'C4', 4)]:
    p = s.Poly(f, x)
    assert p.galois_group(by_name=True)[0].name == name
    assert p.galois_group()[0].order() == order
'''], check=True)
