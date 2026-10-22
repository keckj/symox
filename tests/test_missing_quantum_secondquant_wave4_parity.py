"""Second quantization parity with SymPy 1.14 (expected values recorded from the reference)."""
from itertools import permutations

import pytest

import symox as S
from symox.physics.secondquant import (AntiSymmetricTensor, B, BBra, BKet, Bd, Commutator,
    Dagger, F, FBra, FKet, Fd, FixedBosonicBasis, InnerProduct, KroneckerDelta, NO,
    VarBosonicBasis, apply_operators, evaluate_deltas, matrix_rep, substitute_dummies, wicks)


def _fermi_symbols(cls=S.Symbol):
    i, j, k, l = S.symbols('i j k l', below_fermi=True, cls=cls)
    a, b, c, d = S.symbols('a b c d', above_fermi=True, cls=cls)
    return i, j, k, l, a, b, c, d


def test_bosonic_apply_operator():
    x, y = S.symbols('x y')
    assert repr(B(x).apply_operator(y)) == 'y*AnnihilateBoson(x)'
    assert repr(Dagger(B(x)).apply_operator(y)) == 'y*CreateBoson(x)'
    assert repr(B(0).apply_operator(BKet((3,)))) == 'sqrt(3)*FockStateBosonKet((2,))'
    assert repr(Dagger(B(0)).apply_operator(BKet((3,)))) == '2*FockStateBosonKet((4,))'
    n = S.Symbol('n')
    assert B(0).apply_operator(BKet((n,))) == S.sqrt(n)*BKet((n - 1,))
    assert B(n).apply_operator(BKet([n])) == B(n)*BKet([n])


def test_fermionic_apply_operator():
    p, q, n = S.symbols('p q n')
    assert repr(F(p).apply_operator(q)) == 'q*AnnihilateFermion(p)'
    assert repr(F(1).apply_operator(FKet((1, 2)))) == 'FockStateKet((2,))'
    assert repr(Fd(3).apply_operator(FKet((1, 2)))) == 'FockStateKet((1, 2, 3))'
    assert Fd(1).apply_operator(FKet([n])) == FKet([1, n]) == -FKet([n, 1])
    assert Fd(n).apply_operator(FKet([])) == FKet([n])


def test_bosonic_up_down_and_printing():
    b = BKet([1, 2])
    assert repr((b.up(1), b.down(1), b.up(0))) == '(FockState((1, 3)), FockState((1, 1)), FockState((2, 2)))'
    assert repr((BBra([1, 2]).up(1), BBra([1, 2]).down(1))) == '(FockState((1, 3)), FockState((1, 1)))'
    assert BKet([0, 2]).down(0) == 0
    assert str(BKet((1, 2))) == '|(1, 2)>' and str(BBra((1, 2))) == '<(1, 2)|'
    assert repr(2*BKet((1,))) == '2*FockStateBosonKet((1,))'
    assert len(b) == 2 and b[1] == 2


def test_fermionic_states_follow_the_fermi_sea():
    i, j, k, l, a, b, c, d = _fermi_symbols()
    one, two, five = S.Integer(1), S.Integer(2), S.Integer(5)
    # SymPy needs sympified integer labels; symox also accepts Python ints.
    assert repr((FKet([]).up(one), FKet([2]).up(one), FKet([1]).up(one), FKet([1]).down(one), FKet([2]).down(one))) == \
        '(FockStateKet((1,)), FockStateKet((1, 2)), 0, FockStateKet(()), 0)'
    assert repr((FBra([]).up(one), FBra([1]).down(one))) == '(FockStateKet((1,)), FockStateKet(()))'
    assert repr(FKet([]).up(a)) == 'FockStateKet((a,))'
    assert FKet([], fermi_level=4).up(two) == 0
    assert repr(FKet([], fermi_level=4).up(five)) == 'FockStateKet((5,), fermi_level=4)'
    assert str(FKet((1, 2), fermi_level=4)) == '|(-1, -2)>'
    assert repr(FKet((1, 5), fermi_level=4).down(five)) == '-FockStateFermionKet((1,))'
    assert repr(-FKet((2, 1))) == 'FockStateKet((1, 2))'
    # Holes beyond the Fermi level and Pauli violations vanish.
    assert FKet([0, 1]) == 0 and FKet([i, j]) == 0 and FKet([1, 1]) == 0
    assert Fd(i).apply_operator(FKet([j])) == 0
    assert F(i).apply_operator(FKet([])) == 0
    assert repr(Fd(a).apply_operator(FKet([]))) == 'FockStateKet((a,))'
    assert FBra([1])*FKet([1]) == 1 and InnerProduct(FBra([1]), FKet([2])) == 0
    p = S.Symbol('p')
    restricted = FKet([1], fermi_level=1).up(p)
    assert str(restricted).startswith('-KroneckerDelta(_a, p)*FockStateFermionKet((1, p))')


def test_bosonic_bases():
    b = VarBosonicBasis(5)
    assert repr((b, b.state(3), b.index(b.state(3)), len(b))) == \
        '([FockState((0,)), FockState((1,)), FockState((2,)), FockState((3,)), FockState((4,))], FockState((3,)), 3, 5)'
    b = FixedBosonicBasis(2, 2)
    assert repr((b, b.state(1), b.index(b.state(1)), len(b))) == \
        '([FockState((2, 0)), FockState((1, 1)), FockState((0, 2))], FockState((1, 1)), 1, 3)'
    b = FixedBosonicBasis(2, 3)
    assert repr((b, b.state(3), b.index(b.state(3)))) == \
        '([FockState((2, 0, 0)), FockState((1, 1, 0)), FockState((0, 2, 0)), FockState((1, 0, 1)), FockState((0, 1, 1)), FockState((0, 0, 2))], FockState((1, 0, 1)), 3)'
    assert matrix_rep(Bd(0)*B(0), VarBosonicBasis(4)) == S.diag(0, 1, 2, 3)


def test_commutator():
    i, j, k, l, a, b, c, d = _fermi_symbols()
    p, q = S.symbols('p q')
    m = S.Symbol('m', above_fermi=True)
    assert Commutator(B(0), Bd(0)) == 1 and Commutator(Bd(0), B(0)) == -1
    assert Commutator(B(m), Bd(0)) == KroneckerDelta(m, 0)
    assert Commutator(F(m), Fd(m)) == 1 - 2*NO(Fd(m)*F(m))
    c1 = Commutator(F(p), Fd(q))
    assert Commutator.eval(c1, c1) == 0
    assert repr(Commutator(Fd(a), Fd(i))) == '2*NO(CreateFermion(a)*CreateFermion(i))'
    assert Commutator(Fd(a), F(i)) == -2*NO(F(i)*Fd(a))
    assert repr(Commutator(Fd(p)*Fd(q), F(i)).doit()) == \
        '-AnnihilateFermion(i)*CreateFermion(p)*CreateFermion(q) + CreateFermion(p)*CreateFermion(q)*AnnihilateFermion(i)'
    assert Commutator(Fd(a)*F(i), Fd(b)*F(j)).doit(wicks=True) == 0
    assert Commutator(Fd(j), NO(Fd(a)*F(i))).doit(wicks=True) == -KroneckerDelta(j, i)*Fd(a)
    c = Commutator(Fd(p)*Fd(q), F(i))
    assert repr(c) == 'Commutator(CreateFermion(p)*CreateFermion(q),f(i))'
    assert str(c) == '[CreateFermion(p)*CreateFermion(q),f(i)]'
    assert Commutator(F(i), Fd(p)*Fd(q)) == -c
    X, Y, Z = S.symbols('X Y Z', commutative=False)
    assert Commutator(Commutator(X, Y), Z) != 0


def test_antisymmetric_tensor():
    i, j, k, l, a, b, c, d = _fermi_symbols()
    v = AntiSymmetricTensor('v', (a, i), (b, j))
    assert repr(v) == 'AntiSymmetricTensor(v, (a, i), (b, j))'
    assert repr(AntiSymmetricTensor('v', (i, a), (b, j))) == '-AntiSymmetricTensor(v, (a, i), (b, j))'
    assert repr(AntiSymmetricTensor('t', (a, b), (i, j)).upper) == '(a, b)'
    assert v.symbol == S.Symbol('v') and v.lower == S.Tuple(b, j)
    p, q, r, s = S.symbols('p q r s')
    assert repr(AntiSymmetricTensor('v', (q, p), (r, s))) == '-AntiSymmetricTensor(v, (p, q), (r, s))'
    assert AntiSymmetricTensor('v', (p, p), (r, s)) == 0
    # Substitution canonicalizes again, as the reference rebuild does.
    x = S.Dummy('x')
    assert AntiSymmetricTensor('v', (x, a), ()).subs(x, b) == -AntiSymmetricTensor('v', (a, b), ())


def test_evaluate_deltas_respect_the_fermi_level():
    i, j, k, l, a, b, c, d = _fermi_symbols(S.Dummy)
    p, q = S.symbols('p q', cls=S.Dummy)
    f = S.Function('f')
    assert evaluate_deltas(KroneckerDelta(i, a)*f(a)) == 0
    assert evaluate_deltas(KroneckerDelta(i, p)*f(p)) == f(i)
    assert evaluate_deltas(KroneckerDelta(a, p)*f(p)) == f(a)
    assert evaluate_deltas(KroneckerDelta(i, j)*f(j)*f(i)) == f(i)**2


def test_substitute_dummies_merges_equivalent_terms():
    i, j, k, l, a, b, c, d = _fermi_symbols(S.Dummy)
    p, q, r, s = S.symbols('p q r s', cls=S.Dummy)
    f = S.Function('f')
    assert str(substitute_dummies(f(a, b) + f(c, d))) == '2*f(_a, _b)'
    assert substitute_dummies(f(a, b, c, d) - f(b, a, c, d)) == 0
    expr = AntiSymmetricTensor('v', (a, b), (i, j))*AntiSymmetricTensor('t', (i, j), (a, b)) - \
        AntiSymmetricTensor('v', (c, d), (k, l))*AntiSymmetricTensor('t', (k, l), (c, d))
    assert substitute_dummies(expr) == 0
    assert str(substitute_dummies(f(p, q) + f(r, s) + f(p, s))) == '3*f(_p, _q)'
    assert substitute_dummies(f(p, q)*Fd(p)*F(q) - f(r, s)*Fd(r)*F(s)) == 0
    ccsd = wicks(NO(Fd(a)*F(i))*NO(Fd(j)*F(b)), keep_only_fully_contracted=True)
    assert substitute_dummies(ccsd) == 0


def test_substitute_dummies_ambiguous_orders():
    # Reference test_equivalent_internal_lines_VT2conjT2_ambiguous_order_AT.
    i, j = S.symbols('i j', below_fermi=True, cls=S.Dummy)
    a, b, c, d = S.symbols('a b c d', above_fermi=True, cls=S.Dummy)
    p1, p2, p3, p4 = S.symbols('p1 p2 p3 p4', above_fermi=True, cls=S.Dummy)
    v = lambda *x: AntiSymmetricTensor('v', x[:2], x[2:])
    t = lambda *x: AntiSymmetricTensor('t', x[:2], x[2:])
    for template in (v(p1, p2, p3, p4)*t(p1, p2, i, j)*t(p3, p4, i, j),
                     v(p1, p2, p3, p4)*t(p1, p2, j, i)*t(p3, p4, i, j)):
        orders = permutations([a, b, c, d])
        base = substitute_dummies(template.subs(list(zip([p1, p2, p3, p4], next(orders)))))
        for order in orders:
            assert substitute_dummies(template.subs(list(zip([p1, p2, p3, p4], order)))) == base


def test_wicks_fully_contracted_fock_and_two_body():
    # Reference test_fully_contracted.
    i, j, k, l, a, b, c, d = _fermi_symbols()
    p, q, r, s = S.symbols('p q r s', cls=S.Dummy)
    fock = AntiSymmetricTensor('f', (p,), (q,))*NO(Fd(p)*F(q))
    V = (AntiSymmetricTensor('v', (p, q), (r, s))*NO(Fd(p)*Fd(q)*F(s)*F(r)))/4
    assert wicks(NO(Fd(i)*F(a))*fock, keep_only_fully_contracted=True,
                 simplify_kronecker_deltas=True) == AntiSymmetricTensor('f', (a,), (i,))
    assert wicks(NO(Fd(i)*Fd(j)*F(b)*F(a))*V, keep_only_fully_contracted=True,
                 simplify_kronecker_deltas=True) == AntiSymmetricTensor('v', (a, b), (i, j))
