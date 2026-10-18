"""Operation labels and weighted traversal of compact shared expression DAGs."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox import count_ops


@pytest.mark.parametrize('source',[
    'x','2','-2','Rational(-1,2)','Float(-2)','-oo','E','exp(-1)','exp(-x)',
    'x+y','x-y','-x-y','x+y-z','x*y','-x*y','2*x*y','x/y','-x/y',
    'x/2','2/x','-2/x','x/y/z','x**(-2)','x**(Rational(1,2))',
    'sin(x)+cos(y)','log(x,2)','Max(x,y,z)','Min(x,y)',
    'Eq(x,0)','x>0','And(x>0,x<1)','And(x>0,x<y,y<2)','Or(x>0,y>0)',
    'Piecewise((x*x,x>0),(0,True))','Order(x*x,x)',
    'Integral(x*x,x)','Sum(x,(x,0,5))','Lambda(x,x+1)',
    'Subs(Function("f")(x),x,1)','Derivative(Function("f")(x),x,2)',
])
def test_scalar_operation_labels_match_reference(source):
    x,y,z=sx.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    namespace=vars(sx)|dict(x=x,y=y,z=z);reference=vars(sp)|dict(x=X,y=Y,z=Z)
    actual=eval(source,namespace);expected=eval(source,reference)
    assert count_ops(actual,visual=True).to_sympy()==sp.count_ops(expected,visual=True)
    assert count_ops(actual)==sp.count_ops(expected)
    if isinstance(actual,sx.Expr):
        assert actual.count_ops(visual=True)==count_ops(actual,visual=True)
        assert actual.count_ops()==count_ops(actual)


def test_containers_matrices_arrays_and_sparse_large_shapes():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    pairs=[([x,sx.sin(x),None,True,x+2],[X,sp.sin(X),None,True,X+2]),
        ({x:sx.sin(x),x+2:y+1},{X:sp.sin(X),X+2:Y+1}),
        (sx.Matrix([[x*x,sx.sin(x)],[x+1,0]]),sp.Matrix([[X*X,sp.sin(X)],[X+1,0]])),
        (sx.Array([[x*x,sx.sin(x)],[x+1,0]]),sp.Array([[X*X,sp.sin(X)],[X+1,0]])),
        (sx.FiniteSet(x+1,sx.sin(x)),sp.FiniteSet(X+1,sp.sin(X)))]
    for native,reference in pairs:
        assert count_ops(native,True).to_sympy()==sp.count_ops(reference,True)
    sparse=sx.ImmutableSparseNDimArray({3:x*x},(10**12,))
    assert count_ops(sparse)==1
    matrix=sx.SparseMatrix(10**6,10**6,{(1,2):sx.sin(x)+x*x})
    assert count_ops(matrix,True)==sx.Symbol('ADD')+sx.Symbol('SIN')+sx.Symbol('POW')


def test_shared_graph_counts_tree_occurrences_without_unfolding():
    x=sx.Symbol('x');one=sx.Integer(1)
    value=sx._core._restore_arithmetic('Mul',sx.Integer(2),[(x,one)],False)
    for _ in range(80):
        value=sx._core._restore_arithmetic('Mul',one,[(value,one),(value,one)],False)
    expected=2**81-1
    assert count_ops(value)==expected
    assert count_ops(value,True)==expected*sx.Symbol('MUL')
    assert value.count_ops()==expected


def test_cyclic_container_inputs_are_rejected():
    values=[];values.append(values)
    with pytest.raises(ValueError,match='Cyclic'):count_ops(values)


def test_held_numeric_powers_are_not_evaluated_for_operation_counts():
    x=sx.Symbol('x')
    huge=sx.Pow(2,10**9,evaluate=False)
    value=sx.Mul(huge,x,evaluate=False)
    assert count_ops(value,True)==sx.Symbol('MUL')+sx.Symbol('POW')
    reciprocal=sx.Mul(sx.Pow(2,-10**9,evaluate=False),x,evaluate=False)
    with pytest.raises(NotImplementedError,match='power resource limit'):count_ops(reciprocal)


def test_nsimplify_method_accepts_native_constant_basis():
    value=sx.Float('3.141592653589793')
    assert value.nsimplify(constants=[sx.pi])==sx.pi
    assert sx.Float('0.3333333333333333').nsimplify()==sx.Rational(1,3)


def test_counting_without_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import Symbol,sin,Array
from symox import count_ops
x=Symbol('x')
assert count_ops(x*x+sin(x),True)==Symbol('ADD')+Symbol('POW')+Symbol('SIN')
assert (x*x).count_ops(visual=True)==Symbol('POW')
assert count_ops(Array([x*x,sin(x)]))==2
'''
    subprocess.run([sys.executable,'-c',source],check=True)
