"""Matrix entries retain their typed parents through scalar operations."""
import pickle

import pytest
import sympy as sp
import symox as s


def test_matrix_entries_reconstruct_pickle_and_convert():
    n = s.Symbol("n", integer=True, positive=True)
    i, j = s.symbols("i j", integer=True)
    a = s.MatrixSymbol("A", n, n)
    entry = a[i, j]
    assert entry.args == (a, i, j)
    assert entry.func(*entry.args) == entry
    assert entry.free_symbols == {a, i, j}
    assert entry.has(a) and not entry.has(n)
    assert s.sympify(a) == a
    assert s.sympify(entry.to_sympy()) == entry
    expression = (entry + 1)**2
    assert pickle.loads(pickle.dumps(expression)) == expression
    assert expression.to_sympy() == (sp.MatrixSymbol("A", n.to_sympy(), n.to_sympy())[i.to_sympy(), j.to_sympy()] + 1)**2


def test_scalar_replacements_descend_matrix_metadata_once():
    n, m, p = s.symbols("n m p", integer=True, positive=True)
    a = s.MatrixSymbol("A", n, n)
    b = s.MatrixSymbol("B", n, n)
    entry = (a*b)[0, 0]
    replaced = entry.xreplace({n: m, m: p})
    assert replaced.args[0].shape == (m, m)
    finite = replaced.subs(m, 2)
    aa, bb = s.MatrixSymbol("A", 2, 2), s.MatrixSymbol("B", 2, 2)
    assert finite == aa[0, 0]*bb[0, 0] + aa[0, 1]*bb[1, 0]
    assert (a[0, 0]+a[1, 1]).xreplace({a: b}) == b[0, 0]+b[1, 1]
    concrete = s.Matrix([[2, 3], [5, 7]])
    assert aa[1, 0].subs(aa, concrete) == 5
    assert (aa[0, 1]+1).xreplace({aa: concrete}) == 4
    for method in ("subs", "xreplace"):
        with pytest.raises(TypeError, match="Matrix"):
            getattr(aa[0, 0]+1, method)({aa: 1})
        with pytest.raises(ValueError, match="dimension"):
            getattr(a[0, 0], method)({n: -1})


def test_matrix_entries_native_and_subclassable_code_printers():
    from symox.printing.cxx import CXX17CodePrinter
    a = s.MatrixSymbol("A", 3, 4)
    i, j = s.symbols("i j", integer=True)
    entry = a[i, j]
    assert s.cxxcode(entry) == "A[4*i + j]"
    assert CXX17CodePrinter().doprint(entry) == "A[4*i + j]"
    assert s.numpycode(entry) == "A[i, j]"
    assert s.torchcode(entry) == "A[i, j]"
    assert s.jaxcode(entry) == "A[i, j]"
    for printer in (s.slangcode, s.syclcode):
        with pytest.raises(NotImplementedError, match="matrix"):
            printer(entry)


def test_scalar_matrix_arithmetic_does_not_create_commutative_proxy():
    x = s.Symbol("x")
    a = s.MatrixSymbol("A", 2, 2)
    assert (x*a).shape == (2, 2)
    assert (x*a)[0, 1] == x*a[0, 1]
    for operation in (lambda: x+a, lambda: x-a, lambda: x/a, lambda: x**a):
        with pytest.raises((TypeError, ValueError)):
            operation()
    b = s.MatrixSymbol("B", 2, 2)
    assert s.expand(a*(a+b)) == a*a+a*b
    assert s.simplify(s.MatAdd(a, a)) == 2*a
    assert s.conjugate(a) == a.conjugate()
    assert s.transpose(a) == a.T
    assert s.diff(a, x) == s.ZeroMatrix(2, 2)
