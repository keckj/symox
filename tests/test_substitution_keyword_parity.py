"""Real named substitution calls must bind before native type dispatch."""
import pytest
import symox as sx
import sympy as sp


def equal(actual,expected):
    if isinstance(expected,sp.Poly):
        assert isinstance(actual,sx.Poly)
        assert tuple(g.to_sympy() for g in actual.gens)==expected.gens
        assert str(actual.domain)==str(expected.domain)
        assert actual.as_expr().to_sympy()==expected.as_expr()
    else:assert actual.to_sympy()==expected


def build(module,kind):
    x,y=module.symbols('x y')
    if kind=='scalar':return x+2*y
    if kind=='derivative':return module.Derivative(module.Function('f')(x),x,evaluate=False)
    if kind=='poly':return module.Poly(x**2+y,x,y)
    if kind=='abstract':return x*module.MatrixSymbol('A',2,2)
    return getattr(module,kind)([[x,2*y],[x+y,1]])


@pytest.mark.parametrize('kind',['scalar','derivative','poly','abstract','Matrix',
    'ImmutableMatrix','SparseMatrix','ImmutableSparseMatrix'])
@pytest.mark.parametrize('form',['both_named','second_named','mapping','pairs','none'])
def test_substitution_argument_names_are_runtime_contracts(kind,form):
    def result(module):
        x,y=module.symbols('x y');value=build(module,kind)
        if form=='both_named':return value.subs(arg1=x,arg2=y+1)
        if form=='second_named':
            if module is sp and getattr(value,'is_Matrix',False) and kind!='abstract':
                # MatrixBase's iterator precheck rejects this spelling even
                # though the scalar Basic contract accepts it. Native support
                # is compared with the equivalent accepted positional call.
                with pytest.raises(TypeError):value.subs(x,arg2=y+1)
                return value.subs(x,y+1)
            return value.subs(x,arg2=y+1)
        if form=='mapping':return value.subs(arg1={x:y+1})
        if form=='pairs':return value.subs(arg1=[(x,y+1),(y,2)])
        return value.subs(arg1={x:y+1},arg2=None)
    actual=result(sx);expected=result(sp)
    equal(actual,expected)


@pytest.mark.parametrize('kind',['scalar','poly','abstract','Matrix','SparseMatrix'])
def test_duplicate_substitution_arguments_raise_without_mutation(kind):
    x,y=sx.symbols('x y');value=build(sx,kind);before=str(value)
    for args,kwargs in [((x,),{'arg1':y}),((x,y),{'arg2':1}),
                        ((),{'arg2':1}),((x,y,1),{})]:
        with pytest.raises(TypeError):value.subs(*args,**kwargs)
        assert str(value)==before


@pytest.mark.parametrize('method,name,index',[('row_del','row',-1),('col_del','col',0)])
def test_dense_matrix_deletion_accepts_reference_keyword(method,name,index):
    actual=sx.Matrix([[1,2],[3,4]]);expected=sp.Matrix([[1,2],[3,4]])
    assert getattr(actual,method)(**{name:index}) is None
    assert getattr(expected,method)(**{name:index}) is None
    equal(actual,expected)


@pytest.mark.parametrize('kind',['scalar','poly','Matrix','ImmutableMatrix','SparseMatrix'])
def test_named_mapping_retains_simultaneous_substitution(kind):
    x,y=sx.symbols('x y');a,b=sp.symbols('x y')
    actual=build(sx,kind).subs(arg1={x:y,y:x},simultaneous=True)
    expected=build(sp,kind).subs(arg1={a:b,b:a},simultaneous=True)
    equal(actual,expected)


@pytest.mark.parametrize('kind',['scalar','derivative','poly','abstract','Matrix','SparseMatrix'])
def test_unused_substitution_hints_are_ignored(kind):
    # The reference forwards extra keywords to `_subs`, which ignores them
    # (Lcapy passes `safe=` and `causal=`).
    x,y=sx.symbols('x y');a,b=sp.symbols('x y')
    actual=build(sx,kind).subs(x,y+1,safe=True,causal=False)
    expected=build(sp,kind).subs(a,b+1,safe=True,causal=False)
    equal(actual,expected)
