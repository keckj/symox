"""Control (sympy.physics.control.lti) parity with SymPy 1.14.

Expected strings were produced by the reference SymPy 1.14 on the same
inputs.
"""

import pytest

from symox import I, Matrix, Rational, sqrt, symbols
from symox.physics.control.lti import (
    Feedback, MIMOFeedback, MIMOParallel, MIMOSeries, Parallel, Series,
    StateSpace, TransferFunction, TransferFunctionMatrix,
)

s, a, b, p = symbols('s a b p')
A = Matrix([[-1, 0], [0, -2]]); B = Matrix([[1], [1]]); C = Matrix([[1, 1]]); D = Matrix([[0]])
ss1 = StateSpace(A, B, C, D)
ss2 = StateSpace(Matrix([[0, 1], [-2, -3]]), Matrix([[0], [1]]), Matrix([[1, 0]]), Matrix([[0]]))
tf1 = TransferFunction(s + 1, s**2 + 2*s + 3, s)
tf2 = TransferFunction(s - 2, s + 4, s)
tf3 = TransferFunction(1, s, s)
t1 = TransferFunction(1, s + 1, s); t2 = TransferFunction(2, s + 3, s)
t3 = TransferFunction(s, s**2 + 1, s); t4 = TransferFunction(s - 4, s - 5, s)
tfm1 = TransferFunctionMatrix([[t1, t2], [t3, t4]])
tfm2 = TransferFunctionMatrix([[t2, t1], [t4, t3]])
A3 = Matrix([[1, 0, 0], [0, 2, 0], [0, 0, 3]]); B3 = Matrix([[1], [1], [0]]); C3 = Matrix([[1, 0, 0]]); D3 = Matrix([[0]])
ss3 = StateSpace(A3, B3, C3, D3)

CASES = {
    'ss1': (lambda: ss1, 'StateSpace(Matrix([\n[-1,  0],\n[ 0, -2]]), Matrix([\n[1],\n[1]]), Matrix([[1, 1]]), Matrix([[0]]))'),
    'series_ss': (lambda: Series(ss1, ss2), 'Series(StateSpace(\nMatrix([\n[-1,  0],\n[ 0, -2]]),\n\nMatrix([\n[1],\n[1]]),\n\nMatrix([[1, 1]]),\n\nMatrix([[0]])), StateSpace(\nMatrix([\n[ 0,  1],\n[-2, -3]]),\n\nMatrix([\n[0],\n[1]]),\n\nMatrix([[1, 0]]),\n\nMatrix([[0]])))'),
    'series_ss_doit': (lambda: Series(ss1, ss2).doit(), 'StateSpace(Matrix([\n[-1,  0,  0,  0],\n[ 0, -2,  0,  0],\n[ 0,  0,  0,  1],\n[ 1,  1, -2, -3]]), Matrix([\n[1],\n[1],\n[0],\n[0]]), Matrix([[0, 0, 1, 0]]), Matrix([[0]]))'),
    'parallel_ss': (lambda: Parallel(ss1, ss2), 'Parallel(StateSpace(\nMatrix([\n[-1,  0],\n[ 0, -2]]),\n\nMatrix([\n[1],\n[1]]),\n\nMatrix([[1, 1]]),\n\nMatrix([[0]])), StateSpace(\nMatrix([\n[ 0,  1],\n[-2, -3]]),\n\nMatrix([\n[0],\n[1]]),\n\nMatrix([[1, 0]]),\n\nMatrix([[0]])))'),
    'parallel_ss_doit': (lambda: Parallel(ss1, ss2).doit(), 'StateSpace(Matrix([\n[-1,  0,  0,  0],\n[ 0, -2,  0,  0],\n[ 0,  0,  0,  1],\n[ 0,  0, -2, -3]]), Matrix([\n[1],\n[1],\n[0],\n[1]]), Matrix([[1, 1, 1, 0]]), Matrix([[0]]))'),
    'feedback_ss': (lambda: Feedback(ss1, ss2), 'Feedback(StateSpace(Matrix([\n[-1,  0],\n[ 0, -2]]), Matrix([\n[1],\n[1]]), Matrix([[1, 1]]), Matrix([[0]])), StateSpace(Matrix([\n[ 0,  1],\n[-2, -3]]), Matrix([\n[0],\n[1]]), Matrix([[1, 0]]), Matrix([[0]])), -1)'),
    'feedback_ss_doit': (lambda: Feedback(ss1, ss2).doit(), 'StateSpace(Matrix([\n[-1,  0, -1,  0],\n[ 0, -2, -1,  0],\n[ 0,  0,  0,  1],\n[ 1,  1, -2, -3]]), Matrix([\n[1],\n[1],\n[0],\n[0]]), Matrix([[1, 1, 0, 0]]), Matrix([[0]]))'),
    'series_tf_ss': (lambda: Series(tf1, ss1), 'Series(TransferFunction(s + 1, s**2 + 2*s + 3, s), StateSpace(\nMatrix([\n[-1,  0],\n[ 0, -2]]),\n\nMatrix([\n[1],\n[1]]),\n\nMatrix([[1, 1]]),\n\nMatrix([[0]])))'),
    'series_tf_ss_doit': (lambda: Series(tf1, ss1).doit(), 'StateSpace(Matrix([\n[ 0,  1,  0,  0],\n[-3, -2,  0,  0],\n[ 1,  1, -1,  0],\n[ 1,  1,  0, -2]]), Matrix([\n[0],\n[1],\n[0],\n[0]]), Matrix([[0, 0, 1, 1]]), Matrix([[0]]))'),
    'series_rewrite_tf': (lambda: Series(tf1, tf2).rewrite(TransferFunction), 'TransferFunction((s - 2)*(s + 1), (s + 4)*(s**2 + 2*s + 3), s)'),
    'parallel_rewrite_tf': (lambda: Parallel(tf1, tf2).rewrite(TransferFunction), 'TransferFunction((s - 2)*(s**2 + 2*s + 3) + (s + 1)*(s + 4), (s + 4)*(s**2 + 2*s + 3), s)'),
    'expr_rewrite': (lambda: (tf1*tf2 - tf3*tf1).rewrite(TransferFunction), 'TransferFunction(s*(s - 2)*(s + 1)*(s**2 + 2*s + 3) + (-s - 1)*(s + 4)*(s**2 + 2*s + 3), s*(s + 4)*(s**2 + 2*s + 3)**2, s)'),
    'expr': (lambda: tf1*tf2 - tf3*tf1, 'Parallel(Series(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s)), Series(TransferFunction(-1, 1, s), TransferFunction(1, s, s), TransferFunction(s + 1, s**2 + 2*s + 3, s)))'),
    'nested': (lambda: Series(tf1, Series(tf2, tf3)), 'Series(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s), TransferFunction(1, s, s))'),
    'nested_par': (lambda: Parallel(tf1, Parallel(tf2, tf3)), 'Parallel(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s), TransferFunction(1, s, s))'),
    'nested_doit': (lambda: Series(tf1, Series(tf2, tf3)).doit(), 'TransferFunction((s - 2)*(s + 1), s*(s + 4)*(s**2 + 2*s + 3), s)'),
    'feedback_tf': (lambda: Feedback(tf1, tf2), 'Feedback(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s), -1)'),
    'feedback_tf_doit': (lambda: Feedback(tf1, tf2).doit(), 'TransferFunction((s + 1)*(s + 4)*(s**2 + 2*s + 3), ((s - 2)*(s + 1) + (s + 4)*(s**2 + 2*s + 3))*(s**2 + 2*s + 3), s)'),
    'feedback_tf_neg': (lambda: Feedback(tf1, tf2, 1).doit(), 'TransferFunction((s + 1)*(s + 4)*(s**2 + 2*s + 3), (-(s - 2)*(s + 1) + (s + 4)*(s**2 + 2*s + 3))*(s**2 + 2*s + 3), s)'),
    'tfm1': (lambda: tfm1, 'TransferFunctionMatrix(((TransferFunction(1, s + 1, s), TransferFunction(2, s + 3, s)), (TransferFunction(s, s**2 + 1, s), TransferFunction(s - 4, s - 5, s))))'),
    'mimoseries': (lambda: MIMOSeries(tfm1, tfm2), 'MIMOSeries(TransferFunctionMatrix(((TransferFunction(1, s + 1, s), TransferFunction(2, s + 3, s)), (TransferFunction(s, s**2 + 1, s), TransferFunction(s - 4, s - 5, s)))), TransferFunctionMatrix(((TransferFunction(2, s + 3, s), TransferFunction(1, s + 1, s)), (TransferFunction(s - 4, s - 5, s), TransferFunction(s, s**2 + 1, s)))))'),
    'mimoseries_doit': (lambda: MIMOSeries(tfm1, tfm2).doit(), 'TransferFunctionMatrix(((TransferFunction(s*(s + 1)*(s + 3) + 2*(s + 1)*(s**2 + 1), (s + 1)**2*(s + 3)*(s**2 + 1), s), TransferFunction(4*(s - 5)*(s + 1) + (s - 4)*(s + 3)**2, (s - 5)*(s + 1)*(s + 3)**2, s)), (TransferFunction(s**2*(s - 5)*(s + 1) + (s - 4)*(s**2 + 1)**2, (s - 5)*(s + 1)*(s**2 + 1)**2, s), TransferFunction(s*(s - 5)*(s - 4)*(s + 3) + (s - 5)*(2*s - 8)*(s**2 + 1), (s - 5)**2*(s + 3)*(s**2 + 1), s))))'),
    'mimoparallel': (lambda: MIMOParallel(tfm1, tfm2), 'MIMOParallel(TransferFunctionMatrix(((TransferFunction(1, s + 1, s), TransferFunction(2, s + 3, s)), (TransferFunction(s, s**2 + 1, s), TransferFunction(s - 4, s - 5, s)))), TransferFunctionMatrix(((TransferFunction(2, s + 3, s), TransferFunction(1, s + 1, s)), (TransferFunction(s - 4, s - 5, s), TransferFunction(s, s**2 + 1, s)))))'),
    'mimoparallel_doit': (lambda: MIMOParallel(tfm1, tfm2).doit(), 'TransferFunctionMatrix(((TransferFunction(3*s + 5, (s + 1)*(s + 3), s), TransferFunction(3*s + 5, (s + 1)*(s + 3), s)), (TransferFunction(s*(s - 5) + (s - 4)*(s**2 + 1), (s - 5)*(s**2 + 1), s), TransferFunction(s*(s - 5) + (s - 4)*(s**2 + 1), (s - 5)*(s**2 + 1), s))))'),
    'mimofeedback': (lambda: MIMOFeedback(tfm1, tfm2), 'MIMOFeedback(TransferFunctionMatrix(((TransferFunction(1, s + 1, s), TransferFunction(2, s + 3, s)), (TransferFunction(s, s**2 + 1, s), TransferFunction(s - 4, s - 5, s)))), TransferFunctionMatrix(((TransferFunction(2, s + 3, s), TransferFunction(1, s + 1, s)), (TransferFunction(s - 4, s - 5, s), TransferFunction(s, s**2 + 1, s)))), -1)'),
    'elem_zeros': (lambda: tfm1.elem_zeros(), '[[[], []], [[0], [4]]]'),
    'eval_frequency': (lambda: tfm1.eval_frequency(2), 'Matrix([[1/3, 2/5], [2/5, 2/3]])'),
    'from_Matrix': (lambda: TransferFunctionMatrix.from_Matrix(Matrix([[1/(s+1), s/(s+2)], [0, (s-1)/(s**2+1)]]), s), 'TransferFunctionMatrix(((TransferFunction(1, s + 1, s), TransferFunction(s, s + 2, s)), (TransferFunction(0, 1, s), TransferFunction(s - 1, s**2 + 1, s))))'),
    'from_Matrix_args': (lambda: TransferFunctionMatrix.from_Matrix(Matrix([[1/(s+1), s/(s+2)]]), s).args, '(((TransferFunction(1, s + 1, s), TransferFunction(s, s + 2, s)),),)'),
    'tfm_shape': (lambda: tfm1.shape, '(2, 2)'),
    'ctrl_sub': (lambda: ss3.controllable_subspace(), '[Matrix([\n[1],\n[1],\n[0]]), Matrix([\n[1],\n[2],\n[0]])]'),
    'obs_sub': (lambda: ss3.observable_subspace(), '[Matrix([\n[1],\n[1],\n[1]])]'),
    'ctrl_mat': (lambda: ss3.controllability_matrix(), 'Matrix([[1, 1, 1], [1, 2, 4], [0, 0, 0]])'),
    'obs_mat': (lambda: ss3.observability_matrix(), 'Matrix([[1, 0, 0], [1, 0, 0], [1, 0, 0]])'),
    'is_ctrl': (lambda: ss3.is_controllable(), 'False'),
    'zpk': (lambda: TransferFunction.from_zpk([1, 2], [3, 4], 5, s), 'TransferFunction(5*(s - 2)*(s - 1), (s - 4)*(s - 3), s)'),
    'zpk2': (lambda: TransferFunction.from_zpk([a], [b, -b], 2*a, s), 'TransferFunction(2*a*(-a + s), (-b + s)*(b + s), s)'),
    'zpk3': (lambda: TransferFunction.from_zpk([], [-1, -1], 1, s), 'TransferFunction(1, (s + 1)**2, s)'),
    'ss_rewrite_tf': (lambda: ss1.rewrite(TransferFunction), '[[TransferFunction(2*s + 3, (s + 1)*(s + 2), s)]]'),
    'tf_rewrite_ss': (lambda: tf1.rewrite(StateSpace), 'StateSpace(Matrix([\n[ 0,  1],\n[-3, -2]]), Matrix([\n[0],\n[1]]), Matrix([[1, 1]]), Matrix([[0]]))'),
    'repr_tf': (lambda: repr(tf1), 'TransferFunction(s + 1, s**2 + 2*s + 3, s)'),
    'repr_ss': (lambda: repr(ss1), 'StateSpace(\nMatrix([\n[-1,  0],\n[ 0, -2]]),\n\nMatrix([\n[1],\n[1]]),\n\nMatrix([[1, 1]]),\n\nMatrix([[0]]))'),
    'repr_series': (lambda: repr(Series(tf1, tf2)), 'Series(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s))'),
    'repr_fb': (lambda: repr(Feedback(tf1, tf2)), 'Feedback(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s), -1)'),
    'zeros': (lambda: tf2.zeros(), '[2]'),
    'dc_gain': (lambda: tf2.dc_gain(), '-1/2'),
    'is_stable': (lambda: tf1.is_stable(), 'True'),
    'tf_neg': (lambda: -tf1, 'TransferFunction(-s - 1, s**2 + 2*s + 3, s)'),
    'tf_pow': (lambda: tf2**2, 'TransferFunction((s - 2)**2, (s + 4)**2, s)'),
    'tf_add': (lambda: tf1 + tf2, 'Parallel(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s))'),
    'tf_mul': (lambda: tf1 * tf2, 'Series(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s))'),
    'tf_add_doit': (lambda: (tf1 + tf2).doit(), 'TransferFunction((s - 2)*(s**2 + 2*s + 3) + (s + 1)*(s + 4), (s + 4)*(s**2 + 2*s + 3), s)'),
    'tf_from_rat': (lambda: TransferFunction.from_rational_expression((s+1)/(s**2+2), s), 'TransferFunction(s + 1, s**2 + 2, s)'),
    'tf_to_expr': (lambda: tf2.to_expr(), '(s - 2)/(s + 4)'),
    'series_mul_parallel': (lambda: Series(tf1, tf2) * Parallel(tf2, tf3), 'Series(TransferFunction(s + 1, s**2 + 2*s + 3, s), TransferFunction(s - 2, s + 4, s), Parallel(TransferFunction(s - 2, s + 4, s), TransferFunction(1, s, s)))'),
    'ss_series_tf_rewrite': (lambda: Series(ss1, ss2).rewrite(TransferFunction), 'TransferFunction(2*s + 3, s**4 + 6*s**3 + 13*s**2 + 12*s + 4, s)'),
    'ss_add': (lambda: ss1 + ss2, 'StateSpace(Matrix([\n[-1,  0,  0,  0],\n[ 0, -2,  0,  0],\n[ 0,  0,  0,  1],\n[ 0,  0, -2, -3]]), Matrix([\n[1],\n[1],\n[0],\n[1]]), Matrix([[1, 1, 1, 0]]), Matrix([[0]]))'),
    'ss_mul': (lambda: ss1 * ss2, 'StateSpace(Matrix([\n[ 0,  1,  0,  0],\n[-2, -3,  0,  0],\n[ 1,  0, -1,  0],\n[ 1,  0,  0, -2]]), Matrix([\n[0],\n[1],\n[0],\n[0]]), Matrix([[0, 0, 1, 1]]), Matrix([[0]]))'),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_control_matches_reference(name):
    build, expected = CASES[name]
    assert str(build()) == expected


def test_mixed_nesting_is_kept():
    assert str(Series(tf1, Parallel(tf2, tf3))) == (
        'Series(TransferFunction(s + 1, s**2 + 2*s + 3, s), '
        'Parallel(TransferFunction(s - 2, s + 4, s), TransferFunction(1, s, s)))')
    assert str(Parallel(tf1, Series(tf2, tf3))) == (
        'Parallel(TransferFunction(s + 1, s**2 + 2*s + 3, s), '
        'Series(TransferFunction(s - 2, s + 4, s), TransferFunction(1, s, s)))')
    assert str(Series(tf1, Parallel(tf2, tf3)).doit()) == (
        'TransferFunction((s + 1)*(s*(s - 2) + s + 4), s*(s + 4)*(s**2 + 2*s + 3), s)')


def test_elem_poles_and_zeros_values():
    # Same sets as the reference; only the root order may differ.
    poles = tfm1.elem_poles()
    assert [[set(e) for e in row] for row in poles] == [[{-1}, {-3}], [{I, -I}, {5}]]
    assert [[set(e) for e in row] for row in tfm1.elem_zeros()] == [[set(), set()], [{0}, {4}]]


def test_mimofeedback_doit_equivalent_to_reference():
    expected = MIMOFeedback(tfm1, tfm2).sensitivity * tfm1._expr_mat
    got = MIMOFeedback(tfm1, tfm2).doit()._expr_mat
    assert (got - expected).applyfunc(lambda e: e.cancel()) == Matrix.zeros(2, 2)
