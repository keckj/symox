"""Class selectors span machine and expression-count Derivative nodes."""
import pytest
import sympy as sp
import symox as sx


def tree(module,count):
    x,n=module.symbols('x n');f=module.Function('f')(x)
    order=eval(count,dict(n=n,Float=module.Float,Rational=module.Rational))
    machine=module.Derivative(f,x,evaluate=False)
    other=module.Derivative(f,(x,order),evaluate=False)
    expression=module.Add(machine,other,module.sin(x),evaluate=False)
    return x,n,f,machine,other,expression


def foreign(values):return {v.to_sympy() for v in values}


@pytest.mark.parametrize('count',['2','n','2**32','Float(1)','Rational(1,2)'])
def test_derivative_class_and_instance_atom_queries_select_both_heads(count):
    x,n,f,machine,other,expr=tree(sx,count)
    X,N,F,M,O,E=tree(sp,count)
    for actual,expected in [((sx.Derivative,),(sp.Derivative,)),((machine,),(M,)),((other,),(O,)),((machine,other),(M,O)),((sx.Derivative,sx.Symbol),(sp.Derivative,sp.Symbol))]:
        assert foreign(expr.atoms(*actual))==E.atoms(*expected)
    assert expr.atoms(sx.Derivative)=={machine,other}
    assert expr.atoms((sx.Derivative,sx.Symbol))==E.atoms((sp.Derivative,sp.Symbol))==set()
    assert not other.is_Function
    assert foreign(expr.atoms(sx.Function))==E.atoms(sp.Function)
    assert not (expr.atoms(sx.Function)&{machine,other})


@pytest.mark.parametrize('count',['2','n','2**32','Float(1)','Rational(1,2)'])
def test_has_find_and_count_keep_the_public_derivative_class_protocol(count):
    x,n,f,machine,other,expr=tree(sx,count)
    X,N,F,M,O,E=tree(sp,count)
    assert expr.has(sx.Derivative)==E.has(sp.Derivative)
    assert expr.has(machine,other)==E.has(M,O)
    assert expr.has(sx.Derivative,sx.Symbol)==E.has(sp.Derivative,sp.Symbol)
    # A tuple literal is skipped by the reference (strict ``_sympify``).
    assert E.has((sp.Derivative,sp.Symbol)) is False
    assert expr.has((sx.Derivative,sx.Symbol)) is False
    assert foreign(expr.find(sx.Derivative))==E.find(sp.Derivative)
    assert {k.to_sympy():v for k,v in expr.find(sx.Derivative,group=True).items()}==E.find(sp.Derivative,group=True)
    assert expr.count(sx.Derivative)==E.count(sp.Derivative)==2
    assert expr.find(lambda node:getattr(node,'is_Derivative',False))=={machine,other}
    assert sx.sin(x).atoms(sx.Derivative)==set()
    assert not sx.sin(x).has(sx.Derivative)


def test_custom_derivative_named_type_keeps_instancecheck_and_errors():
    _,_,_,machine,other,expr=tree(sx,'n');seen=[]
    class Meta(type):
        def __instancecheck__(cls,value):
            seen.append(value)
            return bool(getattr(value,'is_Derivative',False))
    class Derivative(metaclass=Meta):pass
    assert expr.atoms(Derivative)=={machine,other}
    assert machine in seen and other in seen
    assert expr.find(Derivative)=={machine,other}
    assert expr.has(Derivative)
    class MarkerError(RuntimeError):pass
    class BrokenMeta(type):
        def __instancecheck__(cls,value):raise MarkerError('selector callback')
    class Derivative(metaclass=BrokenMeta):pass
    with pytest.raises(MarkerError,match='selector callback'):expr.atoms(Derivative)
    with pytest.raises(MarkerError,match='selector callback'):expr.find(Derivative)
    assert expr.atoms(sx.Derivative)=={machine,other}


def test_callable_named_derivative_is_not_mistaken_for_a_native_function_head():
    _,_,_,machine,other,expr=tree(sx,'n')
    def Derivative(value):return bool(getattr(value,'is_Derivative',False))
    assert expr.atoms(Derivative)==set()
    assert expr.find(Derivative)=={machine,other}
    # A callable selector belongs to find; has already refused this literal
    # before the symbolic-order representation and retains that error here.
    with pytest.raises(TypeError,match='cannot convert object of type function'):
        expr.has(Derivative)


def test_ode_and_pde_auto_detection_uses_the_same_complete_derivative_selector():
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    f=sx.Function('f')(x);F=sp.Function('f')(X)
    equation=sx.diff(f,x)
    assert sx.classify_ode(equation)==sx.classify_ode(equation,f)==sp.classify_ode(sp.diff(F,X),F)
    u=sx.Function('u')(x,y);U=sp.Function('u')(X,Y)
    equation=sx.diff(u,x)+sx.diff(u,y)
    assert sx.classify_pde(equation)==sx.classify_pde(equation,u)==sp.classify_pde(sp.diff(U,X)+sp.diff(U,Y),U)
