"""Wave-4 wrong-result regressions in sympy.physics.quantum, wigner and secondquant.

Expected values are the reference SymPy 1.14 results.
"""
import pytest

import symox as S
from symox.physics.quantum import TensorProduct, Dagger, Ket, Bra
from symox.physics.quantum.qapply import qapply
from symox.physics.quantum.represent import represent, enumerate_states
from symox.physics.quantum.state import TimeDepKet, TimeDepBra, Wavefunction
from symox.physics.quantum.qubit import Qubit, IntQubit, IntQubitBra
from symox.physics.quantum.grover import (OracleGate, apply_grover,
    grover_iteration, superposition_basis)
from symox.physics.quantum.spin import (Jx, Jy, Jz, JxKet, JyKet, JzKet,
    JzBra, JxBra, JzKetCoupled, JxKetCoupled)
from symox.physics.quantum.cg import CG, Wigner3j, Wigner6j, Wigner9j
from symox.physics.quantum.cartesian import XKet, XBra, PxKet
from symox.physics.wigner import wigner_3j, wigner_6j, clebsch_gordan

half = S.Rational(1, 2)
r2 = S.sqrt(2)


def test_tensor_product_of_noncommutative_symbols():
    A, B, C, D = S.symbols('A B C D', commutative=False)
    x, y = S.symbols('x y')
    assert str(TensorProduct(A, B)) == 'AxB'
    assert str(TensorProduct(A, B, C)) == 'AxBxC'
    assert str(TensorProduct(2*A, B)) == '2*AxB'
    assert str(TensorProduct(x*A, y*B)) == 'x*y*AxB'
    assert TensorProduct(A, B).args == (A, B)
    assert TensorProduct(A) == A
    assert str(TensorProduct(A, B)*TensorProduct(C, D)) == '(A*C)x(B*D)'
    assert str(TensorProduct(Ket('a'), Ket('b'))) == '|a>x|b>'
    assert str(Dagger(TensorProduct(Ket('a'), Ket('b')))) == '<a|x<b|'
    assert TensorProduct(S.Matrix([1, 2]), S.Matrix([3, 4])) == S.Matrix([3, 4, 6, 8])


@pytest.mark.parametrize('args,kwargs,value,bits', [
    ((1, 1), {}, 3, (1, 1)),
    ((3,), {}, 3, (1, 1)),
    ((5,), {}, 5, (1, 0, 1)),
    ((0,), {}, 0, (0,)),
    ((1,), {}, 1, (1,)),
    ((0, 5), {}, 0, (0, 0, 0, 0, 0)),
    ((3, 2), {}, 3, (1, 1)),
    ((3,), {'nqubits': 4}, 3, (0, 0, 1, 1)),
    ((1, 0, 1), {}, 5, (1, 0, 1)),
])
def test_int_qubit_argument_conventions(args, kwargs, value, bits):
    q = IntQubit(*args, **kwargs)
    assert str(q) == '|%d>' % value
    assert q.qubit_values == bits
    assert q.as_int() == value


def test_int_qubit_bra_qubit_argument_and_errors():
    assert str(IntQubitBra(5)) == '<5|'
    assert str(IntQubitBra(1, 1)) == '<3|'
    assert str(Dagger(IntQubit(5))) == '<5|'
    assert str(IntQubit(Qubit('101'))) == '|5>'
    assert IntQubit(5).nqubits == 3
    assert represent(IntQubit(1, 1)) == S.Matrix([0, 0, 0, 1])
    for args in ((2, 1), (4, 2), (-1,)):
        with pytest.raises(ValueError):
            IntQubit(*args)


def test_grover_oracle_marks_integer_states():
    assert qapply(apply_grover(lambda q: q == IntQubit(2), 2)) == IntQubit(2, nqubits=2)
    assert qapply(apply_grover(lambda q: q == IntQubit(3), 2)) == IntQubit(3, nqubits=2)
    oracle = OracleGate(2, lambda q: q == IntQubit(2))
    assert qapply(grover_iteration(superposition_basis(2), oracle)) == IntQubit(2, nqubits=2)
    assert qapply(oracle*IntQubit(2)) == -IntQubit(2)
    assert qapply(oracle*IntQubit(1, nqubits=2)) == IntQubit(1, nqubits=2)
    basis = superposition_basis(2)
    assert basis == sum((IntQubit(i, nqubits=2) for i in range(4)), S.Integer(0))/2


@pytest.mark.parametrize('state,basis,column', [
    (JzKet(1, -1), Jx, [half, r2/2, half]),
    (JzKet(1, 1), Jx, [half, -r2/2, half]),
    (JzKet(1, 0), Jx, [r2/2, 0, -r2/2]),
    (JzKet(1, -1), Jy, [-half, -r2*S.I/2, half]),
    (JzKet(1, -1), Jz, [0, 0, 1]),
    (JzKet(1, -1), None, [0, 0, 1]),
    (JzKet(half, half), Jx, [-r2/2, r2/2]),
    (JzKet(half, -half), Jy, [r2*S.I/2, -r2/2]),
    (JxKet(1, -1), Jz, [half, -r2/2, half]),
    (JxKet(1, 1), Jz, [half, r2/2, half]),
    (JxKet(1, 0), None, [0, 1, 0]),
    (JxKet(1, 0), Jy, [0, 1, 0]),
    (JyKet(1, 1), Jz, [half, r2*S.I/2, -half]),
    (JyKet(1, 1), Jx, [-S.I, 0, 0]),
    (JzKetCoupled(1, -1, (half, half)), Jx, [0, half, r2/2, half]),
    (JzKetCoupled(1, -1, (half, half)), None, [0, 0, 0, 1]),
    (JzKetCoupled(0, 0, (half, half)), None, [1, 0, 0, 0]),
    (JxKetCoupled(1, 0, (half, half)), Jz, [0, -r2/2, 0, r2/2]),
])
def test_spin_state_representation_basis(state, basis, column):
    options = {} if basis is None else {'basis': basis}
    assert (represent(state, **options) - S.Matrix(column)).applyfunc(S.simplify) == S.zeros(len(column), 1)


def test_spin_bra_representation_basis():
    assert represent(JzBra(1, -1), basis=Jx) == S.Matrix([[half, r2/2, half]])
    assert represent(JxBra(1, 1), basis=Jz) == S.Matrix([[half, r2/2, half]])


def test_wigner_coupling_classes_use_reference_argument_order():
    assert Wigner3j(6, 0, 4, 0, 2, 0).doit() == S.sqrt(715)/143
    assert Wigner3j(1, 1, 1, -1, 0, 0).doit() == S.sqrt(3)/3
    assert Wigner3j(half, half, half, -half, 1, 0).doit() == S.sqrt(6)/6
    assert Wigner3j(2, 1, 2, -1, 0, 0).doit() == -S.sqrt(5)/5
    assert wigner_3j(2, 6, 4, 0, 0, 0) == S.sqrt(715)/143
    assert CG(3, 1, 2, -1, 5, 0).doit() == S.sqrt(105)/21
    assert CG(1, 0, 1, 0, 2, 0).doit() == S.sqrt(6)/3
    assert Wigner6j(1, 2, 3, 2, 1, 2).doit() == S.sqrt(21)/105
    assert Wigner9j(1, 2, 1, 2, 1, 1, 1, 1, 0).doit() == S.Rational(1, 90)
    assert clebsch_gordan(S.Rational(3, 2), half, 1, S.Rational(3, 2), -half, 1) == S.sqrt(3)/2
    assert wigner_6j(5, 5, 5, 5, 5, 5) == S.Rational(1, 52)


def test_enumerate_states_reference_labels():
    assert [str(v) for v in enumerate_states(XKet('foo'), 1, 3)] == ['|foo_1>', '|foo_2>', '|foo_3>']
    assert [str(v) for v in enumerate_states(XKet('foo'), [1, 3, 5])] == ['|foo_1>', '|foo_3>', '|foo_5>']
    assert [str(v) for v in enumerate_states(PxKet('p'), 2, 2)] == ['|p_2>', '|p_3>']
    assert [str(v) for v in enumerate_states(XBra('foo'), 1, 2)] == ['<foo_1|', '<foo_2|']
    assert enumerate_states(Ket('psi'), 1, 3) == []
    with pytest.raises(TypeError):
        enumerate_states(XKet('foo'), 1)


def test_wavefunction_derivative_and_limits():
    x, y = S.symbols('x y')
    assert S.diff(Wavefunction(x**2, x), x) == Wavefunction(2*x, x)
    assert Wavefunction(x**2, x).diff(x) == Wavefunction(2*x, x)
    assert S.diff(Wavefunction(x**2*y, x, y), y) == Wavefunction(x**2, x, y)
    assert str(S.diff(Wavefunction(x**2, x), x)) == 'Wavefunction(2*x, x)'
    assert Wavefunction(x**2, (x, 0, 1)).limits == {x: (0, 1)}
    assert Wavefunction(x**2, (x, 0, 1)).norm == S.sqrt(5)/5


def test_time_dependent_state_label_excludes_time():
    t = S.Symbol('t')
    ket, bra = TimeDepKet('psi', t), TimeDepBra('psi', t)
    assert ket.label == (S.Symbol('psi'),) and ket.time == t
    assert bra.label == (S.Symbol('psi'),) and bra.time == t
    assert str(ket) == '|psi;t>' and str(bra) == '<psi;t|'
    assert str(Dagger(ket)) == '<psi;t|'
    assert TimeDepKet('psi', 'a', t).label == (S.Symbol('psi'), S.Symbol('a'))


# Second quantization with the Fermi vacuum (reference SymPy 1.14).
# Sums are compared as multisets of terms with sorted factors: the order of
# the terms of a sum of operator products is not part of the contract.
from symox.physics.secondquant import (F, Fd, B, NO, wicks, contraction, evaluate_deltas,
    PermutationOperator, simplify_index_permutations, KroneckerDelta,
    SubstitutionOfAmbigousOperatorFailed, ContractionAppliesOnlyToFermions)


def _split(text, seps):
    parts, depth, cur, k = [], 0, '', 0
    while k < len(text):
        ch = text[k]
        if ch in '([':
            depth += 1
        elif ch in ')]':
            depth -= 1
        if depth == 0:
            hit = next((s for s in seps if text.startswith(s, k)), None)
            if hit:
                parts.append(cur); cur = '-' if hit.strip() == '-' else ''; k += len(hit); continue
        cur += ch; k += 1
    parts.append(cur)
    return [p for p in parts if p]


def _terms(text):
    text = text.replace('-1*', '-').replace('+ -', '- ')
    result = []
    for term in _split(text, [' + ', ' - ']):
        negative = term.startswith('-')
        body = term[1:] if negative else term
        noncommutative = [f for f in _split(body, ['*']) if 'Fermion(' in f or f.startswith('NO(')]
        commutative = sorted(f for f in _split(body, ['*']) if f not in noncommutative)
        result.append(('-' if negative else '') + '*'.join(commutative + noncommutative))
    return sorted(result)


def same_sum(value, expected):
    assert _terms(str(value)) == _terms(expected), (str(value), expected)


p, q, r, s = S.symbols('p q r s', cls=S.Dummy)
P, Q, R = S.symbols('p q r')
f = S.Function('f')


def test_no_doit_splits_general_dummies_at_the_fermi_level():
    same_sum(NO(Fd(p)*F(q)).doit(),
             'KroneckerDelta(_a, _p)*KroneckerDelta(_a, _q)*CreateFermion(_a)*AnnihilateFermion(_a)'
             ' + KroneckerDelta(_a, _p)*KroneckerDelta(_i, _q)*CreateFermion(_a)*AnnihilateFermion(_i)'
             ' - KroneckerDelta(_a, _q)*KroneckerDelta(_i, _p)*AnnihilateFermion(_a)*CreateFermion(_i)'
             ' - KroneckerDelta(_i, _p)*KroneckerDelta(_i, _q)*AnnihilateFermion(_i)*CreateFermion(_i)')
    # -(...) distributed: NO(F(p)*Fd(q)) = -NO(Fd(q)*F(p)).
    same_sum(NO(F(p)*Fd(q)).doit(),
             '-KroneckerDelta(_a, _q)*KroneckerDelta(_a, _p)*CreateFermion(_a)*AnnihilateFermion(_a)'
             ' - KroneckerDelta(_a, _q)*KroneckerDelta(_i, _p)*CreateFermion(_a)*AnnihilateFermion(_i)'
             ' + KroneckerDelta(_a, _p)*KroneckerDelta(_i, _q)*AnnihilateFermion(_a)*CreateFermion(_i)'
             ' + KroneckerDelta(_i, _q)*KroneckerDelta(_i, _p)*AnnihilateFermion(_i)*CreateFermion(_i)')
    assert str(NO(Fd(p)*F(q)).doit(remove_brackets=False)) == ':CreateFermion(_p)*AnnihilateFermion(_q):'


def test_no_ordering_and_ambiguous_bracket_removal():
    same_sum(NO(F(P)*Fd(Q)), '-NO(CreateFermion(q)*AnnihilateFermion(p))')
    same_sum(NO(2*Fd(P)*F(Q)+F(P)*Fd(Q)),
             '2*NO(CreateFermion(p)*AnnihilateFermion(q)) - NO(CreateFermion(q)*AnnihilateFermion(p))')
    assert str(NO(F(p)*F(q)*Fd(r)*Fd(s))) == ':CreateFermion(_r)*CreateFermion(_s)*AnnihilateFermion(_p)*AnnihilateFermion(_q):'
    assert str(NO(Fd(p)*F(q)*F(r)).get_subNO(1)) == ':CreateFermion(_p)*AnnihilateFermion(_r):'
    assert NO(Fd(P)) == Fd(P)
    for expr in (NO(Fd(P)*F(Q)), NO(F(0)*Fd(1))):
        with pytest.raises(SubstitutionOfAmbigousOperatorFailed):
            expr.doit()
    no = NO(Fd(p)*F(q))
    assert (no.has_q_creators, no.has_q_annihilators) == (1, 1)
    assert list(no.iter_q_creators()) == [0, 1] and list(no.iter_q_annihilators()) == [1, 0]


def test_contraction_uses_the_fermi_vacuum():
    assert str(contraction(Fd(p), F(q))) == 'KroneckerDelta(_i, _q)*KroneckerDelta(_p, _q)'
    assert str(contraction(F(p), Fd(q))) == 'KroneckerDelta(_a, _q)*KroneckerDelta(_p, _q)'
    assert contraction(Fd(p), Fd(q)) == 0 and contraction(F(p), F(q)) == 0
    assert str(contraction(F(0), Fd(0))) == 'KroneckerDelta(0, _a)'
    assert str(contraction(Fd(0), F(0))) == 'KroneckerDelta(0, _i)'
    with pytest.raises(ContractionAppliesOnlyToFermions):
        contraction(F(P), B(Q))


def test_wicks_fermi_vacuum():
    assert str(wicks(Fd(p)*F(q), keep_only_fully_contracted=True)) == 'KroneckerDelta(_i, _q)*KroneckerDelta(_p, _q)'
    same_sum(wicks(Fd(p)*F(q)), 'KroneckerDelta(_i, _q)*KroneckerDelta(_p, _q) + NO(CreateFermion(_p)*AnnihilateFermion(_q))')
    same_sum(wicks(F(P)*Fd(Q)*F(R)),
             'KroneckerDelta(_a, q)*KroneckerDelta(p, q)*AnnihilateFermion(r)'
             ' + KroneckerDelta(_i, r)*KroneckerDelta(q, r)*AnnihilateFermion(p)'
             ' - NO(CreateFermion(q)*AnnihilateFermion(p)*AnnihilateFermion(r))')
    same_sum(wicks(F(p)*F(q)*Fd(r)*Fd(s), keep_only_fully_contracted=True),
             '-KroneckerDelta(_a, _r)*KroneckerDelta(_a, _s)*KroneckerDelta(_p, _r)*KroneckerDelta(_q, _s)'
             ' + KroneckerDelta(_a, _s)*KroneckerDelta(_a, _r)*KroneckerDelta(_p, _s)*KroneckerDelta(_q, _r)')
    assert str(wicks(Fd(p)*F(q), keep_only_fully_contracted=True, simplify_kronecker_deltas=True)) == 'KroneckerDelta(_i, _p)'
    same_sum(wicks(F(p)*Fd(q)*F(r)*Fd(s), keep_only_fully_contracted=True, simplify_kronecker_deltas=True),
             'KroneckerDelta(_a, _p)*KroneckerDelta(_a, _r) + KroneckerDelta(_a, _p)*KroneckerDelta(_i, _q)')
    same_sum(wicks(F(0)*Fd(0)+Fd(0)*F(0)), 'KroneckerDelta(0, _a) + KroneckerDelta(0, _i)')
    assert str(wicks(F(0)*Fd(0), keep_only_fully_contracted=True)) == 'KroneckerDelta(0, _a)'
    same_sum(wicks(F(0)*F(1)*Fd(0)*Fd(1), keep_only_fully_contracted=True), '-KroneckerDelta(0, _a)*KroneckerDelta(1, _a)')
    assert wicks(Fd(P), keep_only_fully_contracted=True) == 0 and wicks(Fd(P)) == Fd(P) and wicks(0) == 0
    assert wicks(NO(Fd(P)*F(Q)), keep_only_fully_contracted=True) == 0


def test_permutation_operator_get_permuted():
    assert PermutationOperator(p, q).get_permuted(f(p, q)) == -f(q, p)
    assert PermutationOperator(P, Q).get_permuted(f(P)) == f(P)
    assert PermutationOperator(Q, P).args == (P, Q)
    assert PermutationOperator(P, Q).get_permuted(f(P, R)*f(Q, R)) == -f(P, R)*f(Q, R)
    # Reference: f(p, q)*PermutationOperator(p, q).
    assert simplify_index_permutations(f(P, Q)-f(Q, P), [PermutationOperator(P, Q)]) == PermutationOperator(P, Q)*f(P, Q)


def test_evaluate_deltas_reference_cases():
    i, j = S.symbols('i j')
    assert evaluate_deltas(KroneckerDelta(i, j)*f(j)) == f(i)
    assert evaluate_deltas(KroneckerDelta(i, j)*f(R)) == KroneckerDelta(i, j)*f(R)
    assert evaluate_deltas(KroneckerDelta(P, Q)*KroneckerDelta(Q, R)*f(P)) == f(R)
    assert evaluate_deltas(KroneckerDelta(i, j)*f(i)*f(j)) == f(i)**2


def test_operator_printing():
    assert str(Fd(P)) == 'f+(p)' and str(F(P)) == 'f(p)' and str(B(P)) == 'b(p)'
    assert str(NO(Fd(P)*F(Q))) == ':CreateFermion(p)*AnnihilateFermion(q):'
    assert repr(NO(Fd(P)*F(Q))) == 'NO(CreateFermion(p)*AnnihilateFermion(q))'
    assert str(Fd(P)*F(Q)) == 'CreateFermion(p)*AnnihilateFermion(q)'


def test_rotation_applies_with_the_state_projection_as_first_index():
    from symox.physics.quantum.spin import Rotation
    value = qapply(Rotation(0, S.pi/2, 0)*JzKet(half, half))
    assert value == -r2*JzKet(half, -half)/2 + r2*JzKet(half, half)/2
    assert str(value) == '-sqrt(2)*|1/2,-1/2>/2 + sqrt(2)*|1/2,1/2>/2'
