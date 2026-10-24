"""Typed matrix scalars lower to native entries for heterogeneous codegen."""
import pickle

import pytest
import sympy as sp
import symox as sm


def test_held_matrix_scalar_nodes_roundtrip_and_differentiate():
    a,b=sm.MatrixSymbol("A",2,2),sm.MatrixSymbol("B",2,2)
    x=sm.Symbol("x", real=True)
    for node in (sm.Trace(a+b),sm.Determinant(a*b),sm.Trace(x*a),sm.Determinant(x*a)):
        assert node.func(*node.args)==node
        assert pickle.loads(pickle.dumps(node))==node
        assert sm.sympify(node.to_sympy())==node
        assert sm.sympify(sm.srepr(node))==node
        assert node.xreplace({sm.Symbol("absent"):1})==node
    assert sm.Trace(x*a).diff(x)==sm.Trace(a)
    assert sm.Determinant(x*a).diff(x)==2*x*sm.Determinant(a)
    n=sm.Symbol("n", integer=True, positive=True)
    assert sm.Trace(sm.Identity(n)).diff(n)==1
    assert sm.Trace(a+b).doit()==sm.Trace(a)+sm.Trace(b)
    assert sm.Determinant(a*b).doit()==sm.Determinant(a)*sm.Determinant(b)


def test_matrix_scalar_substitutions_preserve_held_nodes():
    a,b=sm.MatrixSymbol("A",2,2),sm.MatrixSymbol("B",2,2)
    x=sm.Symbol("x")
    assert sm.Trace(x*a).xreplace({x:2})==sm.Trace(2*a)
    assert sm.Determinant(x*a).subs(x,2)==sm.Determinant(2*a)
    dense=sm.ImmutableMatrix([[1,2],[3,4]])
    trace=sm.Trace(a).subs(a,dense)
    assert trace==sm.Trace(dense)
    assert trace.doit()==5


def test_complex_scalar_derivatives_remain_held_when_unknown():
    z=sm.Symbol("z")
    x=sm.Symbol("x", real=True)
    for f,reference in [(sm.conjugate,sp.conjugate),(sm.re,sp.re),(sm.im,sp.im)]:
        assert sm.diff(f(z),z).to_sympy()==sp.diff(reference(sp.Symbol("z")),sp.Symbol("z"))
        assert sm.diff(f(x),x).to_sympy()==sp.diff(reference(sp.Symbol("x",real=True)),sp.Symbol("x",real=True))


def test_doit_releases_held_scalar_arithmetic():
    x=sm.Symbol("x")
    held=sm.Add(x,x,evaluate=False)
    assert held.doit()==2*x
    assert sm.UnevaluatedExpr(held).doit()==2*x
    assert sm.UnevaluatedExpr(held).doit(deep=False)==held
    assert sm.Derivative(x*x,x).doit()==2*x


def test_native_matrix_scalar_code_and_torch_jax_sources():
    a=sm.MatrixSymbol("A",2,2)
    expr=sm.Trace(a)+sm.Determinant(a)
    # As in SymPy, the C++ printer does not support matrix scalars.
    with pytest.raises(NotImplementedError, match="Unsupported by .*Determinant"):
        sm.cxxcode(expr)
    for exporter in (sm.extra.torchcode_function,sm.extra.jaxcode_function):
        source=exporter("matrix_scalar",a,expr)
        assert "Trace(" not in source and "Determinant(" not in source
        assert "symox" not in source and "sympy" not in source


def test_unevaluated_tensor_derivatives_remain_explicitly_unsupported():
    a=sm.MatrixSymbol("A",2,2)
    x=sm.Symbol("x")
    # SymPy 1.14: the matrix-line derivative gives Identity(2).
    assert sm.diff(sm.Trace(a),a)==sm.Identity(2)
    for operation in (lambda:sm.Derivative(a,x),lambda:sm.Derivative(sm.Trace(a),a),
                      lambda:sm.Trace(a).diff()):
        with pytest.raises(NotImplementedError,match="matrix|tensor"):
            operation()
