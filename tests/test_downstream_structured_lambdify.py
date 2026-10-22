"""Nested numerical signatures and matrix results preserve their structures."""
import numpy as np
import pytest
import sympy as sp
import symox as sx


def test_nested_signature_with_dynamic_and_derivative_parameters():
    t=sx.Symbol('t');q=sx.Function('q')(t);qd=q.diff(t);x=sx.Symbol('x')
    result=sx.lambdify([[q,qd],[x]],[sx.Matrix([[q,qd],[x,q+qd]]),sx.Matrix([qd*x])],modules='numpy',cse=True)
    actual=result([2,3],[4])
    reference=sp.lambdify([[q.to_sympy(),qd.to_sympy()],[x.to_sympy()]],
        [sp.Matrix([[q.to_sympy(),qd.to_sympy()],[x.to_sympy(),q.to_sympy()+qd.to_sympy()]]),sp.Matrix([qd.to_sympy()*x.to_sympy()])],modules='numpy',cse=True)([2,3],[4])
    assert isinstance(actual,list)
    for value,expected in zip(actual,reference):np.testing.assert_equal(value,expected)
    with pytest.raises(TypeError):result([2,3])
    with pytest.raises(ValueError):result([2],[4])
    with pytest.raises(ValueError):result([2,3,4],[4])
    complex_value=result([2+1j,3-1j],[4])
    np.testing.assert_equal(complex_value[0],[[2+1j,3-1j],[4,5]])
    np.testing.assert_equal(complex_value[1],[[12-4j]])


def test_nested_outputs_keep_tuple_list_and_array_shapes():
    x,y=sx.symbols('x y')
    expression=(x+y,[sx.Matrix([x,y]),(x-y,)])
    value=sx.lambdify((x,y),expression,'numpy')(2,3)
    assert isinstance(value,tuple) and isinstance(value[1],list) and isinstance(value[1][1],tuple)
    assert value[0]==5 and value[1][1]==(-1,)
    np.testing.assert_equal(value[1][0],[[2],[3]])
    vectorized=sx.lambdify([x],sx.Matrix([[x,x**2]]),'numpy')(np.array([2,3]))
    np.testing.assert_equal(vectorized,[[[2,3],[4,9]]])


@pytest.mark.parametrize('shape',[(0,0),(0,2),(2,0)])
def test_empty_matrix_output_matches_reference(shape):
    actual=sx.lambdify([],sx.zeros(*shape),'numpy')()
    expected=sp.lambdify([],sp.zeros(*shape),'numpy')()
    np.testing.assert_equal(actual,expected)


@pytest.mark.parametrize('modules', [['scipy', 'numpy'], ['numpy'], ['numpy', 'mpmath'], ['math'], ['sympy']])
def test_namespace_sinc_receives_the_printed_argument(modules):
    # SymPy's NumPy-family printers call a normalized sinc on x/pi, so a
    # user ``sinc`` entry (Lcapy maps it to numpy.sinc) sees that argument.
    def run(lib):
        n = lib.Symbol('n')
        seen = []
        function = lib.lambdify(n, lib.sinc(n) + 1, [{'sinc': lambda a: seen.append(float(a)) or 7.0}] + modules)
        return function(2), seen
    actual, expected = run(sx), run(sp)
    assert actual[0] == expected[0]
    assert actual[1] == pytest.approx(expected[1])


@pytest.mark.parametrize('source', ['sqrt(t)', '1/sqrt(t)', '2*sqrt(t + 1)', 'sqrt(t) + 1/sqrt(t + 1)', 'sqrt(2)*t', 't**(3/2)'])
@pytest.mark.parametrize('modules', [['numpy'], ['math'], ['mpmath'], ['scipy', 'numpy']])
def test_namespace_sqrt_receives_square_roots(source, modules):
    # SymPy prints b**(1/2) as sqrt(b) and b**(-1/2) as 1/sqrt(b); Lcapy
    # maps ``sqrt`` to a complex-aware square root.
    import cmath

    def run(lib):
        t = lib.Symbol('t')
        expr = lib.sympify(source, locals={'t': t})
        seen = []
        function = lib.lambdify(t, expr, [{'sqrt': lambda a: seen.append(a) or cmath.sqrt(a)}] + modules)
        return complex(function(2j)), sorted(seen, key=lambda v: (v.real, v.imag))
    actual, expected = run(sx), run(sp)
    assert actual == expected
