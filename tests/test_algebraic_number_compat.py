"""Algebraic values retain exact power-basis metadata without SymPy at runtime."""
import pickle
import subprocess
import sys
from fractions import Fraction

import numpy as np
import pytest
import sympy as sp
import symox as sm


def native(text):
    return sm.sympify(text)


def same(a, b):
    delta = sp.expand(sp.sympify(str(a)) - sp.sympify(str(b)))
    assert delta == 0 or sp.simplify(delta) == 0, (a, b, delta)


ROOTS = ['0', '2/3', 'sqrt(2)', '2**(1/3)', 'sqrt(2)+I',
         '2**(1/3)+sqrt(2)', '2**(1/3)+I*sqrt(3)',
         'sqrt(2)+sqrt(3)+sqrt(5)']


@pytest.mark.parametrize('root', ROOTS)
@pytest.mark.parametrize('coefficients', [None, [2, 3], [1, 0, 1], []])
def test_constructor_coordinates_metadata_and_pickle(root, coefficients):
    x = sm.Symbol('z')
    a = sm.AlgebraicNumber(native(root), coefficients, alias='a', gen=x)
    expected = sp.AlgebraicNumber(sp.sympify(root), coefficients, alias='a', gen=sp.Symbol('z'))
    assert a.is_AlgebraicNumber and a.is_number and a.is_algebraic and a.is_finite
    assert a.is_commutative and not a.is_Number and not a.is_Atom
    assert a.free_symbols == set()
    assert a.alias == sm.Symbol('a') and a.root == native(root)
    assert a.coeffs() == [sm.sympify(str(c)) for c in expected.coeffs()]
    assert a.native_coeffs() == [Fraction(int(c.p), int(c.q)) for c in expected.coeffs()]
    same(a.as_expr(), expected.as_expr())
    same(a.minpoly.as_expr(), expected.minpoly.as_expr())
    same(a.as_poly(x).as_expr(), expected.as_poly(sp.Symbol('z')).as_expr())
    same(a.as_expr(x), expected.as_expr(sp.Symbol('z')))
    assert a.is_primitive_element == expected.is_primitive_element
    assert a.func(*a.args) == a
    restored = pickle.loads(pickle.dumps(a))
    assert restored == a and hash(restored) == hash(a)
    assert restored.alias == a.alias and restored.minpoly.gens == a.minpoly.gens
    same(restored.as_expr(), a.as_expr())
    assert complex(a) == pytest.approx(complex(expected.evalf()), rel=1e-12, abs=1e-12)


@pytest.mark.parametrize('root', ROOTS)
def test_composition_field_elements_and_element_minimal_polynomial(root):
    z = sm.Symbol('z')
    a = sm.AlgebraicNumber(native(root), [2, 3], alias='a', gen=z)
    b = sm.AlgebraicNumber(a, [2, -1])
    assert b.alias == a.alias and b.minpoly.gens == a.minpoly.gens
    assert b.root == a.root
    same(b.as_expr(), 2 * sp.sympify(str(a.as_expr())) - 1)
    c = a.field_element([2, -1])
    same(c.as_expr(), 2 * sp.sympify(root) - 1)
    same(a.primitive_element().as_expr(), root)
    d = a.to_primitive_element()
    same(d.root, a.as_expr())
    assert d.is_primitive_element
    if not a.is_primitive_element:
        assert d.alias is None
    same(a.minpoly_of_element().as_expr(), sp.minimal_polynomial(sp.sympify(str(a.as_expr())), sp.Symbol('z')))
    same(a.to_root(), a.as_expr())
    same(a.to_root(minpoly=a.minpoly_of_element()), a.as_expr())


@pytest.mark.parametrize('eta,theta', [
    ('sqrt(2)', 'sqrt(2)+sqrt(3)'),
    ('2**(1/3)', '2**(1/3)+sqrt(2)'),
    ('I*sqrt(3)', '2**(1/3)+I*sqrt(3)'),
    ('sqrt(30)', 'sqrt(2)+sqrt(3)+sqrt(5)'),
    ('sqrt(6)', '(sqrt(2)+sqrt(3)+sqrt(5))**2'),
    ('2/3', 'sqrt(2)'),
])
def test_to_number_field_exact_selected_embedding(eta, theta):
    a = sm.to_number_field(native(eta), native(theta), alias='theta')
    expected = sp.to_number_field(sp.sympify(eta), sp.sympify(theta))
    assert a.root == native(theta)
    assert a.coeffs() == [sm.sympify(str(c)) for c in expected.coeffs()]
    same(a.as_expr(), eta)
    same(a.minpoly_of_element().as_expr(), sp.minimal_polynomial(sp.sympify(eta), sp.Symbol('_x')))
    assert not a.is_primitive_element


@pytest.mark.parametrize('items', [
    ['sqrt(2)', 'sqrt(3)', 'sqrt(5)'],
    ['2**(1/3)', 'sqrt(2)'],
    ['2**(1/3)', 'I*sqrt(3)'],
    ['sqrt(6)', 'sqrt(15)', 'sqrt(7)', 'sqrt(70)'],
    ['sqrt(2)', 'sqrt(8)'], ['1/2', 'sqrt(2)'],
])
def test_to_number_field_joint_generators_and_membership(items):
    a = sm.to_number_field([native(item) for item in items])
    assert a.is_primitive_element
    for item in items:
        b = sm.to_number_field(native(item), a)
        assert b.root == a.root
        same(b.as_expr(), item)
    expected = sp.to_number_field([sp.sympify(item) for item in items])
    # A primitive element need not use the same coefficients as SymPy.
    assert sp.degree(sp.sympify(str(a.minpoly.as_expr()))) == expected.minpoly.degree()


@pytest.mark.parametrize('eta,theta', [
    ('sqrt(3)', 'sqrt(2)'),
    ('I', '2**(1/3)+I*sqrt(3)'),
    ('sqrt(3)', '2**(1/3)+I*sqrt(3)'),
    ('sqrt(2)', '(sqrt(2)+sqrt(3)+sqrt(5))**2'),
])
def test_outside_selected_subfields_fail_without_guessing(eta, theta):
    from symox.polys.polyerrors import IsomorphismFailed
    with pytest.raises(IsomorphismFailed):
        sm.to_number_field(native(eta), native(theta))


def test_constructor_pairs_are_certified_and_invalid_options_are_explicit():
    x = sm.Symbol('x')
    a = sm.AlgebraicNumber((x*x-2, sm.sqrt(2)), [2, 3])
    same(a.as_expr(), 2*sp.sqrt(2)+3)
    b = sm.AlgebraicNumber((3*x-2, sm.Rational(2, 3)))
    same(b.as_expr(), sp.Rational(2, 3))
    for polynomial, root in [(x*x-3, sm.sqrt(2)), (x**4-4, sm.sqrt(2)), (0, 0),
                             (x*x-4, 2), (x+sm.Symbol('y'), 1),
                             (x*x-sm.Float(2), sm.sqrt(2))]:
        with pytest.raises((TypeError, ValueError, NotImplementedError)):
            sm.AlgebraicNumber((polynomial, root))
    for coefficients in [[sm.Float(2), 1], [sm.sqrt(2), 1], [x, 1]]:
        with pytest.raises((TypeError, ValueError)):
            sm.AlgebraicNumber(sm.sqrt(2), coefficients)
    with pytest.raises(ValueError):
        sm.AlgebraicNumber(sm.sqrt(2), alias=2)
    with pytest.raises(ValueError):
        sm.to_number_field([])
    with pytest.raises(NotImplementedError):
        a.to_root(radicals=False)
    with pytest.raises(NotImplementedError):
        a.to_primitive_element(radicals=False)
    with pytest.raises(ValueError):
        a.to_root(minpoly=x*x-2)


def test_structural_identity_scalar_and_numpy_interop():
    a = sm.AlgebraicNumber(sm.sqrt(2))
    root = sm.sqrt(2)
    assert a != root and root != a
    assert not (a == root) and not (root == a)
    assert len({a, root}) == 2
    assert sm.AlgebraicNumber(root) == a
    x = sm.Symbol('x')
    same(a+1, sp.sqrt(2)+1)
    same(1+a, sp.sqrt(2)+1)
    same(a-1, sp.sqrt(2)-1)
    same(1-a, 1-sp.sqrt(2))
    same(a*x, sp.sqrt(2)*sp.Symbol('x'))
    same(x*a, sp.sqrt(2)*sp.Symbol('x'))
    same(a/2, sp.sqrt(2)/2)
    same(2/a, sp.sqrt(2))
    same(a**2, 2)
    reference = sp.AlgebraicNumber(sp.sqrt(2))
    assert sp.sympify(reference) is reference
    assert sm.sympify(a) is a
    assert sm.sympify(a).as_expr() == root
    assert sm.Poly(x+a,x,extension=a).degree() == 1
    same(np.dot(np.array([a,1],dtype=object),np.array([a,2],dtype=object)),4)
    zero=a.field_element([0])
    assert zero.is_zero and zero.is_rational and zero.is_integer
    with pytest.raises(ZeroDivisionError):
        1/zero
    assert float(a) == pytest.approx(2**0.5)
    assert complex(sm.AlgebraicNumber(sm.I)) == 1j
    assert a.is_real is True and a.is_positive is True
    assert sm.AlgebraicNumber(sm.I).is_real is False


def test_import_paths_and_native_only_process():
    from symox.polys.numberfields import to_number_field
    from symox.polys.numberfields.subfield import to_number_field as from_subfield
    from symox.polys.polyerrors import IsomorphismFailed
    assert to_number_field is from_subfield is sm.to_number_field
    error = IsomorphismFailed('not in the selected field')
    restored = pickle.loads(pickle.dumps(error))
    assert type(restored) is IsomorphismFailed and restored.args == error.args
    script = r'''
import sys, importlib.abc
class NoSymPy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise AssertionError('runtime SymPy import forbidden')
sys.meta_path.insert(0, NoSymPy())
import symox as s, pickle
x=s.Symbol('x')
theta=2**s.Rational(1,3)+s.sqrt(2)
a=s.to_number_field(s.sqrt(2),theta)
assert a.is_AlgebraicNumber
assert pickle.loads(pickle.dumps(a)) == a
assert s.Poly(x+a,x,extension=theta).degree() == 1
b=s.AlgebraicNumber(s.sqrt(2)+s.sqrt(3)+s.sqrt(5))
assert b.minpoly.degree() == 8
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable, '-c', script], check=True)


@pytest.mark.parametrize("source", ["0", "1", "-2", "sqrt(2)", "-sqrt(2)",
                                    "I", "-I", "1+I", "sqrt(2)+I*sqrt(3)"])
def test_nonzero_means_real_nonzero_for_field_elements(source):
    value = sm.AlgebraicNumber(sm.sympify(source), alias="a")
    assert value.is_nonzero == value.as_expr().is_nonzero
    if not value.is_real:
        assert value.is_zero is False
        assert value.is_nonzero is False
    assert pickle.loads(pickle.dumps(value)).is_nonzero == value.is_nonzero
