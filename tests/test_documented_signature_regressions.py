"""Runtime regressions found by the versioned public API audit."""
from itertools import islice
import pytest
import sympy as sp
import symox as sx
from symox.utilities.iterables import iterable, is_sequence


@pytest.mark.parametrize('name', ['Integers','Naturals','Naturals0','Rationals'])
def test_countable_set_order_is_exact_and_has_no_duplicates(name):
    actual=list(islice(getattr(sx.S,name),200))
    expected=list(islice(getattr(sp.S,name),200))
    assert [x.to_sympy() for x in actual] == expected
    assert len(set(actual)) == 200
    assert getattr(sx.S,name).is_iterable is True


def test_tuple_protocol_and_scalar_noniterability():
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=sx.Tuple(1,x,[2,3]);expected=sp.Tuple(1,X,(2,3))
    assert len(actual)==len(expected)
    assert [v.to_sympy() for v in actual]==list(expected)
    for key in [0,-1,slice(None,None,-1),slice(1,3),slice(None,None,2)]:
        assert actual[key].to_sympy()==expected[key]
    assert iterable(actual) and is_sequence(actual)
    for scalar in [x,sx.Integer(2),sx.sin(x)]:
        assert not iterable(scalar) and not is_sequence(scalar)
        with pytest.raises(TypeError):iter(scalar)
        with pytest.raises(TypeError):len(scalar)
        with pytest.raises(TypeError):scalar[0]
    assert sx.Symbol.diff(x*x,x)==2*x


@pytest.mark.parametrize('num', ['1.23456789123456789','1e20','1.e20','12e3',
    '12.e3','1.230000000000000000000','0.00000001234567890123456789',
    '1_234.567_890_123_456_789','-1.234567890123456789e-900',10**80])
@pytest.mark.parametrize('options',[{}, {'dps':''}, {'precision':''}, {'dps':40}, {'precision':180}])
def test_float_infers_and_preserves_input_precision(num,options):
    actual=sx.Float(num=num,**options);expected=sp.Float(num=num,**options)
    assert actual._prec==expected._prec
    # SymPy 1.14's explicit-precision underscore parser misplaces the decimal
    # point. Compare the mathematical literal; do not reproduce that bug.
    if isinstance(num,str) and '_' in num:
        expected=sp.Float(num.replace('_',''),**options)
    assert actual.to_sympy()==expected


def test_named_math_arguments_and_root_branches():
    assert sx.Integer(i=12)==12
    assert sx.sqrt(arg=4)==2
    assert sx.cbrt(arg=8)==2
    assert sx.log(arg=8,base=2)==3
    assert sx.N(x=sx.Rational(1,3),n=30)._prec==sp.N(sp.Rational(1,3),30)._prec
    for arg in [-8,2,sx.Symbol('x')]:
        for n in [2,3,5,-3,sx.Rational(2,3)]:
            for k in range(4):
                actual=sx.root(arg=arg,n=n,k=k).to_sympy()
                expected=sp.root(sp.sympify(str(arg)),sp.sympify(str(n)),k)
                assert sp.simplify(actual-expected)==0


def test_documented_structural_helpers():
    x=sx.Symbol('x')
    assert sx.FiniteSet((1,2),frozenset({2,3})).to_sympy()==sp.FiniteSet((1,2),{2,3})
    assert sx.Dict({x:1}).has(sx.Symbol)
    assert sx.Dict({x:1}).has(x)
    assert not sx.Dict({x:1}).has(sx.sin)
    assert sx.DiracDelta(x).is_simple(x)
    assert not sx.DiracDelta(sx.cos(x)).is_simple(x)
    n=sx.Symbol('n',positive=True,integer=True)
    for shape in [(2,3),(3,2),(n+1,n),(n,n+1)]:
        for name in ['DiagonalMatrix','DiagonalOf']:
            actual=getattr(sx,name)(sx.MatrixSymbol('A',*shape)).diagonal_length
            reference_shape=[sp.sympify(str(v),locals={'n':sp.Symbol('n',positive=True,integer=True)}) for v in shape]
            expected=getattr(sp,name)(sp.MatrixSymbol('A',*reference_shape)).diagonal_length
            assert actual.to_sympy()==expected


@pytest.mark.parametrize('digits',[15,60,120])
@pytest.mark.parametrize('points',[1,(1,2)])
def test_nsolve_short_signature_and_secant_points(digits,points):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expected=sp.nsolve(X*X-2,points,prec=digits)
    for actual in [sx.nsolve(x*x-2,points,prec=digits),
                   sx.nsolve(x*x-2,x,points,prec=digits,solver='secant')]:
        assert actual._prec==expected._prec
        assert abs(actual.to_sympy()-expected)<sp.Rational(1,10**max(10,digits-3))
    root=sx.nsolve(x*x+1,1+1j,prec=digits)
    assert abs(root.to_sympy()-sp.I)<sp.Rational(1,10**max(10,digits-3))


def test_boolean_constants_after_substitution_and_number_truth():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    for name in ['And','Or','Xor','Equivalent','Implies']:
        actual=getattr(sx,name)(x,y)
        expected=getattr(sp,name)(X,Y)
        for value in [0,1]:
            assert actual.subs(x,value).to_sympy()==expected.subs(X,value)
    for name in ['Not','Xor','Equivalent','Implies']:
        for value in [0,1,2,sx.Rational(1,2),sx.Float(0),sx.oo,sx.pi,sx.I]:
            values=(value,) if name=='Not' else (value,x)
            reference=tuple(v.to_sympy() if hasattr(v,'to_sympy') else v for v in values)
            assert getattr(sx,name)(*values).to_sympy()==getattr(sp,name)(*reference)
    assert sx.DiracDelta(sx.pi)==0
    assert sx.Heaviside(sx.Symbol('positive',positive=True))==1
