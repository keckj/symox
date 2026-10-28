"""Native scalar operation histograms preserve reference and host protocols."""
import pytest
import sympy as sp
import symox as sx


def expressions(source):
    x,y,z = sx.symbols('x y z')
    X,Y,Z = sp.symbols('x y z')
    return (eval(source,vars(sx)|dict(x=x,y=y,z=z)),
            eval(source,vars(sp)|dict(x=X,y=Y,z=Z)))


@pytest.mark.parametrize('source',[
    '-Rational(7,3)', 'Float(-2)', 'Float(-2)*x', 'Float(-2)+x+y',
    '-x-y-z', 'x-y-z', '-x-y+z', '-x/(2*y)', '2*x/(3*y)',
    '(x+5)/(2*(8*x*x+14*x+3))', 'x*exp(-y)', 'x*exp(-2*y)',
    'x*y**(-z)', 'x*y**(-2*z)', 'x*y**(-sin(z))',
    'sin(x)+cos(y)+Function("f")(x,y)', 'Max(x,y,z)', 'Min(x,y,z)',
    'Heaviside(x)', 'Heaviside(x,Rational(1,3))', 'DiracDelta(x,2)',
    'hyper([Rational(1,3),Rational(2,3)],[2],x)',
    'Derivative(Function("f")(x),x,7)', 'Integral(x*x,(x,0,y+1))',
    'Sum(x*x,(x,1,5))', 'Product(x,(x,1,5))',
    'Subs(Function("f")(x),x,1+y)', 'Limit(sin(x)/x,x,0)',
    'UnevaluatedExpr(x+y+z)', 'Order(x*x,x)', 'Lambda(x,x+y+z)',
    'Tuple(x+y+z,Eq(x,y))', 'Eq(x+y+z,y)',
    'And(x>0,y>0,z>0)', 'Or(x>0,y>0)', 'Not(x>0,evaluate=False)',
    'Piecewise((x*x,And(x>0,y>0,z>0)),(0,True))',
    'Q.positive(x+y+z)', 'Q.eq(x+y,z)', 'Contains(x,Interval(-1,1))',
    'Mul(-oo,x,evaluate=False)', 'Mul(-2,x+y,evaluate=False)',
    'Add(-x,-y,-z,evaluate=False)', 'Pow(E,x,evaluate=False)',
])
def test_standard_scalar_histograms_match_reference(source):
    actual,expected=expressions(source)
    assert sx.count_ops(actual,True).to_sympy()==sp.count_ops(expected,True)
    assert sx.count_ops(actual)==sp.count_ops(expected)
    if isinstance(actual,sx.Expr):
        assert actual.count_ops(visual=True)==sx.count_ops(actual,True)


@pytest.mark.parametrize('name',['Order','Tuple','FiniteSet','Max','Min','And','Or','Not',
                                 'Integral','Derivative','sin','f','ß','\ua7d1'])
def test_undefined_names_keep_function_prefix_and_python_unicode_case(name):
    x=sx.Symbol('collision_x');X=sp.Symbol('collision_x')
    actual=sx.Function(name)(x)
    expected=sp.Function(name)(X)
    assert sx.count_ops(actual,True).to_sympy()==sp.count_ops(expected,True)
    assert sx.count_ops(actual)==1


def test_plain_native_scalars_do_not_enter_python_graph_walker(monkeypatch):
    import symox._operation_counts as implementation
    implementation._COUNTS.clear()
    def unavailable(*args,**kwargs):
        raise AssertionError('native scalar counting entered the Python graph walker')
    monkeypatch.setattr(implementation,'_count_ops_fallback',unavailable)
    x=sx.Symbol('native_only_count_x')
    for expression,expected in [(sx.sin(x)+x*x,3),(sx.Heaviside(x),2),
                                (sx.UnevaluatedExpr(x+1),1),
                                (sx.Function('Max')(x),1)]:
        assert sx.count_ops(expression)==expected
        assert sx.count_ops(expression,True).is_Add or expected==1


def test_host_function_protocol_retains_python_fallback(monkeypatch):
    import symox._operation_counts as implementation
    class Host(sx.Function):
        is_Function=False
    x=sx.Symbol('host_count_x')
    value=Host(x+1)
    assert sx._core._native_count_ops(value,False) is None
    original=implementation._count_ops_fallback
    seen=[]
    def recording(expression,visual=False):
        seen.append(expression)
        return original(expression,visual)
    monkeypatch.setattr(implementation,'_count_ops_fallback',recording)
    implementation._COUNTS.clear()
    assert sx.count_ops(value)==1
    assert sx.count_ops(value,True)==sx.Symbol('ADD')
    assert seen==[value,value]


def test_unhashable_host_keeps_counting_protocol():
    class Host(sx.Function):
        __hash__=None
    x=sx.Symbol('unhashable_count_x')
    value=Host(x+1)
    assert sx._core._native_count_ops(value,True) is None
    assert sx.count_ops(value,True)==sx.Symbol('ADD')+sx.Symbol('HOST')


def test_shared_held_counts_keep_big_integers_and_visual_multiplicity():
    x=sx.Symbol('weighted_count_x');one=sx.Integer(1)
    expression=sx._core._restore_arithmetic('Mul',sx.Integer(2),[(x,one)],False)
    for _ in range(100):
        expression=sx._core._restore_arithmetic('Mul',one,[(expression,one),(expression,one)],False)
    expected=2**101-1
    assert sx.count_ops(expression)==expected
    assert sx.count_ops(expression,True)==expected*sx.Symbol('MUL')


def test_huge_held_numeric_powers_are_inspected_without_materializing_them():
    x=sx.Symbol('held_count_x')
    positive=sx.Pow(2,10**9,evaluate=False)
    assert sx.count_ops(sx.Mul(positive,x,evaluate=False),True)==sx.Symbol('MUL')+sx.Symbol('POW')
    for power in [positive,sx.Pow(2,-10**9,evaluate=False)]:
        value=sx.Mul(power,1/x,evaluate=False)
        with pytest.raises(NotImplementedError,match='numeric power resource limit'):
            sx.count_ops(value)


@pytest.mark.parametrize('kind',['indexed','base','index','arithmetic'])
def test_indexed_frontend_metadata_keeps_shape_and_bound_operations(kind):
    def build(module):
        n=module.Symbol('metadata_count_n',integer=True)
        i=module.Idx('metadata_count_i',(1,n-1))
        base=module.IndexedBase('metadata_count_A',shape=(2*n+3,))
        return {'indexed':base[i], 'base':base, 'index':i,
                'arithmetic':base[i]+base[n+1]}[kind]
    actual,expected=build(sx),build(sp)
    assert sx.count_ops(actual,True).to_sympy()==sp.count_ops(expected,True)
    assert sx.count_ops(actual)==sp.count_ops(expected)


@pytest.mark.parametrize('source',[
    'Mul(2,2,evaluate=False)',
    'Mul(Rational(-2,3),1,evaluate=False)',
    'Mul(-1,Rational(2,3),2,evaluate=False)',
    'Mul(1,Mul(-1,1,evaluate=False),evaluate=False)',
    'Mul(0,Mul(-1,Float("2.5"),x,evaluate=False),evaluate=False)',
    'Add(Mul(-1,Float("1.25"),evaluate=False),Mul(1,Float("1.25"),evaluate=False),evaluate=False)',
    'Q.positive(Order(x*x,x))',
    'Q.positive(Piecewise((x*x,x>0),(0,True)))',
])
def test_held_normalization_and_basic_presentation_modes(source):
    actual,expected=expressions(source)
    assert sx.count_ops(actual,True).to_sympy()==sp.count_ops(expected,True)
    assert sx.count_ops(actual)==sp.count_ops(expected)


@pytest.mark.parametrize('error',[KeyboardInterrupt,SystemExit])
def test_host_protocol_interrupts_are_not_bypassed(error):
    class Host(sx.Function):
        @property
        def is_Function(self):
            raise error
    value=Host(sx.Symbol('interrupt_count_x'),evaluate=False)
    assert sx._core._native_count_ops(value,False) is None
    with pytest.raises(error):
        sx.count_ops(value,True)
