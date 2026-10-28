"""Exact native held orders; evaluation rules absent from the engine stay held."""
import pickle
import pytest
import sympy as sp
import symox as sx


def both(body='sin(x)', orders=('n',), assumptions=None, evaluate=False):
    out=[]
    for api in (sx,sp):
        x,y=api.symbols('x y')
        n=api.Symbol('n',**(assumptions or {}));m=api.Symbol('m')
        scope=dict(x=x,y=y,n=n,m=m,f=api.Function('f'),sin=api.sin,exp=api.exp,I=api.I,nan=api.nan,pi=api.pi,R=api.Rational,F=api.Float)
        expr=eval(body,scope)
        pairs=[(x,eval(order,scope)) for order in orders]
        out.append((api.Derivative(expr,*pairs,evaluate=evaluate),scope))
    return out


def equal(actual,expected):
    assert actual.to_sympy()==expected
    assert str(actual)==str(expected)


@pytest.mark.parametrize('order',['n','m+n','R(1,2)','pi','I','nan','F(0)','F(1)','F(2)','4294967296'])
@pytest.mark.parametrize('body',['sin(x)','exp(x)','f(x)','x**2'])
def test_held_orders_preserve_exact_args_and_head(body,order):
    (actual,_),(expected,_)=both(body,(order,))
    equal(actual,expected)
    assert actual.is_Derivative and not actual.is_Function
    assert actual.func is sx.Derivative
    assert actual.func(*actual.args)==actual
    assert actual.variable_count[0][1].to_sympy()==expected.variable_count[0][1]


@pytest.mark.parametrize('assumptions',[{},dict(integer=True),dict(integer=True,nonnegative=True),dict(real=True),dict(zero=True)])
def test_symbol_assumptions_do_not_coerce_count(assumptions):
    (actual,_),(expected,_)=both(assumptions=assumptions)
    equal(actual,expected)
    assert actual.derivative_count.to_sympy()==expected.derivative_count
    assert actual.free_symbols=={sx.Symbol('x'),sx.Symbol('n',**assumptions)}


@pytest.mark.parametrize('order',['n','R(1,2)','pi','F(0)','F(2)'])
@pytest.mark.parametrize('body',['sin(x)','exp(x)','f(x)','x**2'])
def test_evaluation_holds_when_no_symbolic_formula(body,order):
    (actual,_),(expected,_)=both(body,(order,),evaluate=True)
    equal(actual,expected)
    equal(actual.doit(),expected.doit())


@pytest.mark.parametrize('value',[0,1,2,3,sp.Rational(1,2),sp.Float(2),sp.I,sp.nan])
@pytest.mark.parametrize('simultaneous',[False,True])
def test_count_specialization_matches_held_reference(value,simultaneous):
    (actual,a),(expected,b)=both('f(x)')
    actual_value=sx.sympify(value)
    equal(actual.subs(a['n'],actual_value,simultaneous=simultaneous),expected.subs(b['n'],value,simultaneous=simultaneous))


@pytest.mark.parametrize('operation',['subs','xreplace'])
@pytest.mark.parametrize('value',[-1,-2])
def test_negative_specializations_error_and_next_call_is_healthy(operation,value):
    (actual,a),_=both()
    for _ in range(2):
        with pytest.raises(ValueError,match='order of differentiation must be nonnegative'):
            getattr(actual,operation)({a['n']:value})
    assert actual.subs(a['n'],0)==sx.sin(a['x'])
    assert sx.diff(sx.sin(a['x']),a['x'])==sx.cos(a['x'])


@pytest.mark.parametrize('bad',['I','nan'])
def test_nonreal_order_construction_and_evaluation_have_distinct_errors(bad):
    (actual,_),(expected,_)=both(orders=(bad,))
    equal(actual,expected)
    with pytest.raises(TypeError) as error:expected.doit()
    with pytest.raises(TypeError,match=str(error.value)):actual.doit()
    assert sx.Derivative(sx.sin(sx.Symbol('x')),(sx.Symbol('x'),sx.Symbol('n'))).subs(sx.Symbol('n'),0)==sx.sin(sx.Symbol('x'))


def test_invalid_counts_cannot_cancel_before_error_is_reported():
    x,y,n=sx.symbols('x y n')
    # Both held derivatives become the same invalid object; canonical Add must
    # not silently erase the invalid count and return zero.
    expr=sx.Derivative(sx.sin(x),(x,n))-sx.Derivative(sx.sin(x),(x,-n-2))
    with pytest.raises(ValueError):expr.subs(n,-1)
    assert sx.Derivative(sx.sin(y),(y,n)).subs(n,2).doit()==-sx.sin(y)


def test_grouping_commutation_and_canonical():
    results=[]
    for api in (sx,sp):
        x,y,n,m=api.symbols('x y n m');f=api.Function('f')
        held=api.Derivative(f(x,y),(x,n),(y,m),(x,2))
        results.append((held,held.canonical,api.Derivative(f(x,y),(x,n),(y,m),(x,2),evaluate=True),api.Derivative(api.Derivative(f(x),(x,n)),x)))
    for actual,expected in zip(*results):equal(actual,expected)


@pytest.mark.parametrize('operation',['diff','doit','subs_x_y','subs_x_one','subs_y_x','simultaneous','xreplace'])
def test_derivative_variables_substitution_and_formal_differentiation(operation):
    results=[]
    for api in (sx,sp):
        x,y,n=api.symbols('x y n');f=api.Function('f');e=api.Derivative(f(x),(x,n))
        result={'diff':lambda:e.diff(x),'doit':lambda:e.doit(),'subs_x_y':lambda:e.subs(x,y),'subs_x_one':lambda:e.subs(x,1),'subs_y_x':lambda:api.Derivative(f(y),(x,n)).subs(y,x),'simultaneous':lambda:e.subs({x:1,n:2},simultaneous=True),'xreplace':lambda:e.xreplace({x:y,n:2})}[operation]()
        results.append(result)
    equal(*results)


def test_order_metadata_is_free_but_not_differentiated():
    x,y,n=sx.symbols('x y n');f=sx.Function('f');e=sx.Derivative(f(y),(x,n))
    assert e.free_symbols=={y,n}
    assert e._wrt_variables==[x]
    assert e.diff(n)==0
    with pytest.raises(TypeError,match='Cannot give expansion for symbolic count'):e.variables
    assert e.subs(n,2).variables==(x,x)


@pytest.mark.parametrize('count',['n','n+1','R(1,2)'])
def test_printing_operation_counts_and_pickle(count):
    (actual,_),(expected,_)=both('f(x)',(count,))
    assert sx.srepr(actual)==sp.srepr(expected)
    assert sx.latex(actual)==sp.latex(expected)
    assert sx.pretty(actual)==sp.pretty(expected)
    assert sx.pretty(actual,use_unicode=False)==sp.pretty(expected,use_unicode=False)
    assert sx.count_ops(actual)==sp.count_ops(expected)
    assert str(sx.count_ops(actual,visual=True))==str(sp.count_ops(expected,visual=True))
    restored=pickle.loads(pickle.dumps(actual))
    assert restored==actual and restored.func is sx.Derivative
    equal(sx.sympify(str(actual)),expected)


def test_cse_never_replaces_a_derivative_body_with_a_temporary():
    x,n=sx.symbols('x n');f=sx.Function('f');e=sx.Derivative(f(x),(x,n))
    replacements,reduced=sx.cse([e+1,e+2])
    assert len(replacements)==1 and replacements[0][1]==e
    assert len(reduced)==2


@pytest.mark.parametrize('printer',['ccode','cxxcode','fcode','rust_code'])
def test_code_generation_explicitly_refuses_symbolic_order(printer):
    (actual,_),_=both()
    with pytest.raises((NotImplementedError,ValueError,TypeError)):
        getattr(sx,printer)(actual)


@pytest.mark.parametrize('source',['Derivative(sin(x),(x,n))','Derivative(sin(x),(x,1/2))','Derivative(sin(x),(x,4294967296))','Derivative(sin(x),x,4294967296)'])
def test_parser_retains_symbolic_and_wide_exact_counts(source):
    equal(sx.sympify(source),sp.sympify(source))


def test_nested_host_substitution_catches_its_own_error_without_poisoning_caller():
    x,n,t=sx.symbols('x n t')
    derivative=sx.Derivative(sx.sin(x),(x,n))
    class Caught(sx.Function):
        def _eval_subs(self,old,new):
            if old==t:
                with pytest.raises(ValueError):derivative.subs(n,-1)
                assert derivative.subs(n,0)==sx.sin(x)
                return sx.Integer(7)
    result=(Caught(t)+derivative).subs(t,3)
    assert result==7+derivative
    assert derivative.subs(n,1).doit()==sx.cos(x)


@pytest.mark.parametrize('exception',[RuntimeError('host sentinel'),KeyboardInterrupt('host interrupted')])
def test_host_exception_retains_priority_and_error_frame_is_always_cleared(exception):
    x,n,t=sx.symbols('x n t')
    derivative=sx.Derivative(sx.sin(x),(x,n))
    class Failing(sx.Function):
        def _eval_subs(self,old,new):
            if old==t:raise exception
    with pytest.raises(type(exception),match=str(exception)):
        (Failing(t)+derivative).subs(t,2)
    assert derivative.subs(n,0)==sx.sin(x)
    assert pickle.loads(pickle.dumps(derivative))==derivative


@pytest.mark.parametrize('body',['sin(x)','f(x)'])
def test_commutativity_and_diff_wrt_use_body_instead_of_count(body):
    x=sx.Symbol('x');n=sx.Symbol('n',commutative=False);f=sx.Function('f')
    expr=eval(body,dict(x=x,sin=sx.sin,f=f))
    derivative=sx.Derivative(expr,(x,n))
    assert derivative.is_commutative is True
    assert derivative.free_symbols=={x,n}
    if body=='f(x)':
        with pytest.raises(TypeError):derivative._diff_wrt
    else:assert derivative._diff_wrt is False


@pytest.mark.parametrize('order',['n','R(1,2)','F(0)','F(2)','I','nan'])
def test_differentiating_an_independent_body_returns_zero_before_merging(order):
    (actual,a),(expected,b)=both('f(y)',(order,))
    equal(actual.diff(a['x']),expected.diff(b['x']))


@pytest.mark.parametrize('order',['I','nan'])
def test_differentiating_a_nonreal_held_order_raises_at_evaluation(order):
    (actual,a),(expected,b)=both('sin(x)',(order,))
    with pytest.raises(TypeError):expected.diff(b['x'])
    with pytest.raises(TypeError):actual.diff(a['x'])
    assert actual.subs(a['x'],a['y']).is_Derivative


def test_float_zero_count_stays_inexact_when_merged_or_counted():
    (actual,a),(expected,b)=both('sin(x)',('F(0)',))
    assert actual.derivative_count.to_sympy()==expected.derivative_count
    equal(actual.diff(a['x']),expected.diff(b['x']))


@pytest.mark.parametrize('first,second',[('R(1,2)','n'),('R(1,2)','R(3,2)'),('17','n')])
def test_derivative_variables_use_exact_node_ranks(first,second):
    results=[]
    for api in (sx,sp):
        x,y,n,m=api.symbols('x y n m');f=api.Function('f');g=api.Function('g')
        env=dict(R=api.Rational,n=n)
        a=api.Derivative(f(x),(x,eval(first,env)))
        b=api.Derivative(g(y),(y,eval(second,env)))
        result=api.Derivative(a*b,(b,m),(a,n)).canonical
        results.append(result)
    equal(*results)


@pytest.mark.parametrize('raw_boolean',[False,True])
def test_paired_python_boolean_counts_follow_reference_preprocessing(raw_boolean):
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=sx.Derivative(sx.sin(x),(x,raw_boolean))
    expected=sp.Derivative(sp.sin(X),(X,raw_boolean))
    equal(actual,expected)
    if raw_boolean:
        with pytest.raises(TypeError):actual.doit()
        with pytest.raises(TypeError):actual.derivative_count
        assert pickle.loads(pickle.dumps(actual))==actual


@pytest.mark.parametrize('boolean_name',['true','false'])
def test_boolean_expressions_as_counts_and_substitutions_remain_held(boolean_name):
    results=[]
    for api in (sx,sp):
        x,n=api.symbols('x n');count=getattr(api,boolean_name)
        e=api.Derivative(api.sin(x),(x,count))
        result=api.Derivative(api.sin(x),(x,n)).subs(n,count)
        results.append((e,result))
    for actual,expected in zip(*results):equal(actual,expected)
    for actual in results[0]:
        with pytest.raises(TypeError):actual.doit()
        with pytest.raises(TypeError):sx.Derivative(actual,sx.Symbol('x'))


def test_python_integer_subclasses_keep_integer_protocol_and_exact_order():
    class Plain(int):pass
    class WithSympyHook(int):
        def _sympy_(self):return sp.Rational(1,2)
    class WithIntegerOverride(int):
        def __int__(self):return 7
    for count in [Plain(2),WithSympyHook(2),WithIntegerOverride(2),Plain(4294967296)]:
        x=sx.Symbol('x');X=sp.Symbol('x')
        equal(sx.Derivative(sx.sin(x),(x,count)),sp.Derivative(sp.sin(X),(X,count)))


def test_integer_count_hook_error_is_not_retried_or_replaced():
    x=sx.Symbol('x');n=sx.Symbol('n')
    class InvalidInteger(int):
        def __int__(self):raise RuntimeError('integer conversion sentinel')
    for args in [((x,InvalidInteger(2)),),(x,InvalidInteger(2))]:
        with pytest.raises(RuntimeError,match='integer conversion sentinel'):
            sx.Derivative(sx.sin(x),*args)
    assert sx.Derivative(sx.sin(x),(x,n)).subs(n,0)==sx.sin(x)


@pytest.mark.parametrize('simultaneous',[False,True])
def test_host_error_keeps_priority_over_an_invalid_count_in_the_same_call(simultaneous):
    x,n=sx.symbols('x n');derivative=sx.Derivative(sx.sin(x),(x,n))
    class Failing(sx.Function):
        def _eval_subs(self,old,new):
            if old==n:raise RuntimeError('same call host sentinel')
    for _ in range(2):
        with pytest.raises(RuntimeError,match='same call host sentinel'):
            (Failing(n)+derivative).subs(n,-1,simultaneous=simultaneous)
    assert derivative.subs(n,0)==sx.sin(x)


@pytest.mark.parametrize('body,assumptions',[('f(x)',{}),('f(y)',{}),('f(y)',dict(positive=True)),('f(f(x))',{}),('sin(x)',{})])
def test_symbolic_derivative_as_variable_uses_reference_doit_check(body,assumptions):
    (actual,_),(expected,_)=both(body,assumptions=assumptions)
    assert actual._diff_wrt is expected._diff_wrt


def test_a_host_returning_a_derivative_mid_call_still_reports_a_later_invalid_order():
    x,n,t=sx.symbols('x n t');derivative=sx.Derivative(sx.sin(x),(x,n))
    class Returned(sx.Function):
        def _eval_subs(self,old,new):
            if old==t:return derivative
    with pytest.raises(ValueError,match='order of differentiation must be nonnegative'):
        Returned(t).subs([(t,1),(n,-1)])
    assert derivative.subs(n,0)==sx.sin(x)


@pytest.mark.parametrize('order',['n','n+1','-n','n/2'])
def test_noncommutative_count_is_free_metadata_but_has_zero_formal_derivative(order):
    (actual,a),(expected,b)=both('sin(x)',(order,),assumptions=dict(commutative=False))
    assert a['n'] in actual.free_symbols
    equal(actual.diff(a['n']),expected.diff(b['n']))
    assert actual.is_commutative is True


@pytest.mark.parametrize('container',['product','sum','matrix','matrixexpr'])
def test_traversals_propagate_invalid_orders_instead_of_returning_a_candidate(container):
    for api in (sx,sp):
        x,n=api.symbols('x n');d=api.Derivative(api.sin(x),(x,api.I))
        e={'product':lambda:x*d,'sum':lambda:d+x,'matrix':lambda:api.Matrix([[d]]),'matrixexpr':lambda:api.MatMul(d,api.MatrixSymbol('A',2,2),evaluate=False)}[container]()
        with pytest.raises(TypeError):e.diff(x)
    x,n=sx.symbols('x n');d=sx.Derivative(sx.sin(x),(x,n));matrix=sx.Matrix([[d]])
    with pytest.raises(ValueError):matrix.xreplace({n:-1})
    assert matrix.subs(n,0)==sx.Matrix([[sx.sin(x)]])


@pytest.mark.parametrize('operation',['expand','simplify','factor_terms','evalf'])
def test_uncertified_counter_transformations_keep_the_original_operator_held(operation):
    a,x=sx.symbols('a x');count=(a+1)**2-a*a-2*a-2
    derivative=sx.Derivative(sx.sin(x),(x,count))
    if operation in ('expand','evalf'):result=getattr(derivative,operation)()
    else:result=getattr(sx,operation)(derivative)
    assert result==derivative
    assert result.is_Derivative and result.variable_count[0][1]==count
    # This intentionally held scope differs from reference expand, which
    # proves the order negative and raises. It must never become NaN or zero.
    A,X=sp.symbols('a x');reference=sp.Derivative(sp.sin(X),(X,(A+1)**2-A*A-2*A-2))
    if operation=='expand':
        with pytest.raises(ValueError):reference.expand()


@pytest.mark.parametrize('container',['matrix','matrixexpr'])
def test_matrix_substitution_and_pickle_preserve_an_unmodified_counter(container):
    x,n,t=sx.symbols('x n t');d=sx.Derivative(sx.sin(x),(x,n))
    matrix=sx.Matrix([[d,t]]) if container=='matrix' else d*sx.MatrixSymbol('A',2,2)
    unchanged=matrix.subs(t,3)
    assert unchanged.has(d)
    restored=pickle.loads(pickle.dumps(unchanged))
    assert restored==unchanged


@pytest.mark.parametrize('digits',[30,80])
def test_float_count_precision_survives_held_construction_merging_and_pickle(digits):
    results=[]
    for api in (sx,sp):
        x=api.Symbol('x');zero=api.Float('0.0',digits)
        e=api.Derivative(api.sin(x),(x,zero))
        results.append((e,e.diff(x),pickle.loads(pickle.dumps(e))))
    for actual,expected in zip(*results):equal(actual,expected)


@pytest.mark.parametrize('operation',['subs','xreplace'])
def test_a_nonfinite_context_keeps_typed_counter_substitutions_and_its_guard(operation):
    x,n,a=sx.symbols('x n a');derivative=sx.Derivative(sx.sin(x),(x,n))
    domain_guard=sx.Not(sx.Contains(a,sx.S.Complexes),evaluate=False)
    guarded=sx.And(domain_guard,sx.Eq(sx.exp(a*x)+derivative,0,evaluate=False),evaluate=False)
    with pytest.raises(ValueError,match='order of differentiation must be nonnegative'):
        getattr(guarded,operation)({n:-1})
    specialized=getattr(guarded,operation)({n:0})
    assert specialized.has(sx.sin(x))
    assert specialized.has(domain_guard)
    assert not specialized.has(derivative)


def test_nonreal_order_errors_also_cross_other_formal_variable_derivatives():
    for api in (sx,sp):
        x,y=api.symbols('x y');f=api.Function('f')
        derivative=api.Derivative(f(x,y),(x,api.I))
        for expr in [derivative,y*derivative,api.Matrix([[derivative]])]:
            with pytest.raises(TypeError):expr.diff(y)


def test_host_function_jet_protocol_is_explicitly_unsupported_not_assumed_true():
    x,n=sx.symbols('x n')
    class Host(sx.Function):pass
    jet=sx.Derivative(Host(x),(x,n))
    with pytest.raises(NotImplementedError,match='Host _diff_wrt'):
        jet._diff_wrt
    with pytest.raises(NotImplementedError,match='Host _diff_wrt'):
        sx.Derivative(jet,(jet,n))
    assert jet.subs(n,0)==Host(x)


def test_new_counts_reject_unknown_host_symbol_protocol_but_keep_machine_path():
    n=sx.Symbol('n')
    class HostSymbol(sx.Symbol):pass
    x=HostSymbol('x')
    # An unchanged literal u32 order stays on the established path.
    assert sx.Derivative(x*x,(x,1)).doit()==2*x
    with pytest.raises(NotImplementedError,match='host or specialized Symbol'):
        sx.Derivative(x*x,(x,n))
    ordinary=sx.Symbol('t');derivative=sx.Derivative(sx.sin(ordinary),(ordinary,n))
    with pytest.raises(NotImplementedError,match='host or specialized Symbol'):
        derivative.subs(ordinary,x)
    assert derivative.subs(n,0)==sx.sin(ordinary)


@pytest.mark.parametrize('value',[0.0,-0.0,1.0])
def test_python_float_pairs_are_preprocessed_before_becoming_float_expressions(value):
    x=sx.Symbol('x');X=sp.Symbol('x')
    equal(sx.Derivative(sx.sin(x),(x,value)),sp.Derivative(sp.sin(X),(X,value)))


def test_integer_subclass_zero_preprocessing_differs_between_pair_and_positional_order():
    class Converted(int):
        def __int__(self):return 7
    x=sx.Symbol('x');X=sp.Symbol('x');count=Converted(0)
    equal(sx.Derivative(sx.sin(x),(x,count)),sp.Derivative(sp.sin(X),(X,count)))
    equal(sx.Derivative(sx.sin(x),x,count),sp.Derivative(sp.sin(X),X,count))


def test_pair_count_equality_hook_error_is_preserved_and_frame_is_cleared():
    class EqualityError(int):
        def __eq__(self,other):raise RuntimeError('pair equality sentinel')
    x,n=sx.symbols('x n')
    with pytest.raises(RuntimeError,match='pair equality sentinel'):
        sx.Derivative(sx.sin(x),(x,EqualityError(2)))
    assert sx.Derivative(sx.sin(x),(x,n)).subs(n,0)==sx.sin(x)


def test_diff_wrt_property_preserves_host_count_errors_and_clears_the_frame():
    x=sx.Symbol('x');f=sx.Function('f');state={'raise':False}
    class Count(sx.Function):
        def _eval_is_negative(self):return None
        def _eval_is_extended_real(self):
            if state['raise']:raise RuntimeError('diff_wrt count sentinel')
            return None
    n=Count(x);derivative=sx.Derivative(f(x),(x,n))
    state['raise']=True
    with pytest.raises(RuntimeError,match='diff_wrt count sentinel'):
        derivative._diff_wrt
    state['raise']=False
    assert sx.Derivative(sx.sin(x),(x,sx.Symbol('m'))).subs(sx.Symbol('m'),0)==sx.sin(x)
