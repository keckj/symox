"""Field linear algebra and exact ideal/quotient identities."""
import pickle
import random

import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('kind',range(10))
def test_native_scalars_on_the_left_of_module_vectors(kind):
    x=sx.Symbol('x')
    domains=[sx.ZZ,sx.QQ,sx.GF(7),sx.QQ.algebraic_field(sx.sqrt(2)),
             sx.ZZ_I,sx.QQ_I,sx.RR,sx.CC,sx.EX,sx.QQ.poly_ring(x)]
    domain=domains[kind]
    scalar=domain.convert(x+1 if domain.is_Composite else 2)
    module=domain.free_module(3)
    vector=module([1,2,3])
    assert scalar*vector==vector*scalar
    assert (scalar*vector).module==module
    assert all(value.parent()==domain for value in scalar*vector)


def test_native_scalars_delegate_to_quotient_and_extension_arithmetic():
    from symox.polys.agca.extensions import MonogenicFiniteExtension
    ring=sx.ZZ/[15]
    a=ring(4)
    scalar=sx.ZZ(2)
    assert scalar+a==ring(6)
    assert scalar-a==ring(13)
    assert scalar*a==ring(8)
    assert scalar/a==ring(8)
    assert a/scalar==ring(2)
    x=sx.Symbol('x')
    extension=MonogenicFiniteExtension(sx.Poly(x*x-2,x,domain=sx.QQ))
    generator=extension.generator
    rational=sx.QQ(2)
    assert rational*generator==generator+generator
    assert rational+generator==extension(2)+generator
    assert rational-generator==extension(2)-generator
    assert rational/generator==generator
    with pytest.raises(sx.CoercionFailed):
        sx.sympify(generator)


@pytest.mark.parametrize('index', range(3))
def test_field_submodules(index):
    k = [sx.QQ, sx.GF(7), sx.QQ.algebraic_field(sx.sqrt(2))][index]
    module = k.free_module(3)
    e1, e2, e3 = module.basis
    assert e1+e2 == module([1, 1, 0])
    assert 2*e1-e2 == module([2, -1, 0])
    sub = module.submodule([1, 2, 0], [2, 4, 0], [0, 1, 3])
    assert not sub.is_full_module() and not sub.is_zero()
    rng = random.Random(32)
    for _ in range(60):
        a, b = k(rng.randrange(-40, 41)), k(rng.randrange(-40, 41))
        v = module([a, 2*a+b, 3*b])
        assert v in sub
        assert v+e3 not in sub
    assert module.submodule(*module.basis).is_full_module()
    assert module.submodule().is_zero()
    assert module.submodule().contains([0, 0, 0])
    assert sub.submodule([1, 2, 0]).contains([2, 4, 0])
    with pytest.raises(ValueError):
        sub.submodule([0, 0, 1])
    assert pickle.loads(pickle.dumps(e1)) == e1
    assert pickle.loads(pickle.dumps(sub)).contains([2, 5, 3])


def test_free_module_coordinate_validation():
    assert sx.ZZ.free_module(2)([1, 2])*3 == [3, 6]
    assert sx.QQ.free_module(0).basis == ()
    assert sx.QQ.free_module(0).submodule().is_full_module()
    with pytest.raises(ValueError):
        sx.QQ.free_module(-1)
    with pytest.raises(TypeError):
        sx.QQ.free_module(1.2)
    with pytest.raises(ValueError):
        sx.QQ.free_module(2)([1])
    with pytest.raises(ValueError):
        sx.QQ.free_module(2)(sx.QQ.free_module(3)([1, 2, 3]))
    assert [1, sx.sqrt(2)] not in sx.QQ.free_module(2)


@pytest.mark.parametrize('a,b', [(0, 0), (6, 15), (-12, 18), (2, 3), (13, 0)])
def test_integer_ideals_and_quotients(a, b):
    import math
    ideal = sx.ZZ.ideal(a, b)
    g = math.gcd(a, b)
    assert ideal == sx.ZZ.ideal(g)
    assert ideal.is_zero() == (g == 0)
    assert ideal.is_whole_ring() == (g == 1)
    ring = sx.ZZ / ideal
    for n in range(-20, 21):
        assert ideal.contains(n) == (n % g == 0 if g else n == 0)
        assert int(ring(n).data) == (n % g if g else n)
        for m in [-5, 0, 7]:
            assert ring(n)+ring(m) == ring(n+m)
            assert ring(n)*ring(m) == ring(n*m)
    if g > 1:
        for n in range(1, g):
            if math.gcd(n, g) == 1:
                assert ring(n)**-1 * ring(n) == ring.one
    assert pickle.loads(pickle.dumps(ring(2))) == ring(2)


@pytest.mark.parametrize('ground', [sx.QQ, sx.GF(7)])
def test_field_ideals(ground):
    assert ground.ideal().is_zero()
    assert ground.ideal(0).contains(0)
    assert not ground.ideal(0).contains(1)
    assert ground.ideal(0, 2).is_whole_ring()
    ring = ground.quotient_ring([])
    a = ring(3)
    assert a**-1*a == ring.one


@pytest.mark.parametrize('index', range(5))
def test_polynomial_quotients_against_independent_remainder(index):
    x = sx.Symbol('x')
    X = sp.Symbol('x')
    modulus = [x*x+1, x**3-x+1, x**4, x*x-2, x**3-2][index]
    k = sx.QQ.old_poly_ring(x)
    ideal = k.ideal(modulus)
    ring = k.quotient_ring(ideal)
    rng = random.Random(index)
    for _ in range(20):
        a = sum(rng.randrange(-10, 11)*x**j for j in range(7))
        b = sum(rng.randrange(-10, 11)*x**j for j in range(4))
        for actual, expression in [(ring(a)*ring(b), a*b),
                                   (ring(a)+ring(b), a+b),
                                   (ring(a)**3, a**3)]:
            reference = sp.rem(sx.sympify(expression).to_sympy(), modulus.to_sympy(), X)
            assert ring.to_sympy(actual).to_sympy() == reference
    if index != 2:
        u = ring(x)
        assert u**-1*u == ring.one
    assert pickle.loads(pickle.dumps(ring(x))) == ring(x)


def test_multivariate_ideal_membership_and_arithmetic():
    x, y = sx.symbols('x y')
    k = sx.QQ.old_poly_ring(x, y)
    i = k.ideal(x*x-y, y*y-1)
    j = k.ideal(x-y)
    assert i.contains(x**4-1)
    assert not i.contains(x-y)
    assert (i+j).contains(x-y)
    assert (i*j).contains((x*x-y)*(x-y))
    assert not (i*j).contains(x-y)
    assert (i*j).subset(i) and i.subset(i+j)
    quotient = k / i
    assert quotient(x)**4 == quotient.one
    assert quotient(x)**2 == quotient(y)
    with pytest.raises(ValueError):
        k.quotient_ring(sx.ZZ.ideal(3))


@pytest.mark.parametrize('n', [1, 2, 3, 4, 6, 8])
def test_supported_cyclotomic_fields(n):
    k = sx.QQ.cyclotomic_field(n, ss=True)
    if n <= 2:
        assert k is sx.QQ
    else:
        assert k.ext.alias == sx.Symbol('zeta'+str(n))
        assert k.unit**n == k.one
        for d in sp.divisors(n)[:-1]:
            assert k.unit**int(d) != k.one


@pytest.mark.parametrize('root_index', [0, 1, -1, -2])
def test_algebraic_field_from_polynomial(root_index):
    x = sx.Symbol('x')
    k = sx.QQ.alg_field_from_poly(sx.Poly(x*x-2, x), alias='a', root_index=root_index)
    assert k.unit**2 == k(2)
    assert k.ext.alias == sx.Symbol('a')
    assert (k.ext.root.is_positive is True) == (root_index in [1, -1])
    assert k.ext.minpoly.as_expr().subs(k.ext.minpoly.gens[0], k.ext.root).expand() == 0
