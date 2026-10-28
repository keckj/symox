"""Literal integer closure certificates; symbolic/floating counts stay held."""
import pickle
from collections import Counter
import pytest
import sympy as sp
import symox as sx

ORDERS=[0,1,2,3,4,63,64,2**31-1,2**31,2**32-1,2**32,2**32+1,2**64-1,2**100+2,2**128+3]


def run(body,x,order,mode):
    if mode=='diff':return sx.diff(body,x,order)
    if mode=='doit':return sx.Derivative(body,(x,order),evaluate=False).doit()
    if mode=='unsimplified':return sx.diff(body,x,order,simplify=False)
    return sx.Derivative(body,(x,order),evaluate=True)


@pytest.mark.parametrize('function',['sin','cos'])
@pytest.mark.parametrize('order',ORDERS)
@pytest.mark.parametrize('mode',['diff','doit','unsimplified','constructor'])
def test_complete_periodic_classes_across_integer_storage_boundaries(function,order,mode):
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=run(getattr(sx,function)(x),x,order,mode)
    # D^4 is the identity on both functions; the oracle derives at most 3 times.
    expected=sp.diff(getattr(sp,function)(X),X,order%4)
    assert actual.to_sympy()==expected


@pytest.mark.parametrize('order',ORDERS)
@pytest.mark.parametrize('mode',['diff','doit','unsimplified'])
def test_exponential_fixed_point_has_no_order_loop(order,mode):
    x=sx.Symbol('x')
    assert run(sx.exp(x),x,order,mode)==sx.exp(x)


@pytest.mark.parametrize('source',['x','1','7','I','pi','x**7-2*x**3+5','(x+1)**8','(x+1)*(x+2)*(x+3)'])
@pytest.mark.parametrize('order',[2**32-1,2**32,2**128+7])
@pytest.mark.parametrize('mode',['diff','doit','unsimplified'])
def test_finite_polynomials_vanish_only_beyond_certified_degree(source,order,mode):
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=run(eval(source,dict(x=x,I=sx.I,pi=sx.pi)),x,order,mode)
    # Exact degree is an independent certificate for all higher orders. The
    # reference Mul override otherwise enumerates N-way Leibniz partitions.
    polynomial=sp.Poly(eval(source,dict(x=X,I=sp.I,pi=sp.pi)),X)
    degree=int(polynomial.degree())
    assert order>degree
    expected=sp.diff(polynomial.as_expr(),X,degree+1)
    assert actual==0 and actual.to_sympy()==expected


@pytest.mark.parametrize('source',['sin(2*x+1)','cos(-2*x+1)','exp(3*x-2)','(2*x+3)**8','x**12','3*x**6+5*x**3-2','sin(2*x)+x**5','(1+I)*x**5+sqrt(2)*(x+1)**4'])
@pytest.mark.parametrize('order',range(10))
def test_closed_form_recurrences_match_reference_for_entire_small_families(source,order):
    x=sx.Symbol('x');X=sp.Symbol('x')
    actual=sx.diff(eval(source,dict(x=x,sin=sx.sin,cos=sx.cos,exp=sx.exp,sqrt=sx.sqrt,I=sx.I)),x,order)
    expected=sp.diff(eval(source,dict(x=X,sin=sp.sin,cos=sp.cos,exp=sp.exp,sqrt=sp.sqrt,I=sp.I)),X,order)
    assert sp.expand(actual.to_sympy()-expected)==0


@pytest.mark.parametrize('function',['sin','cos','exp'])
@pytest.mark.parametrize('slope',['0','1','-1','I'])
@pytest.mark.parametrize('order',[0,1,2**32-1,2**128+2])
def test_degenerate_and_unit_complex_affine_slopes(function,slope,order):
    x=sx.Symbol('x');X=sp.Symbol('x');a=eval(slope,dict(I=sx.I));A=eval(slope,dict(I=sp.I))
    body=getattr(sx,function)(a*x+sx.pi/3)
    actual=sx.diff(body,x,order)
    if order==0:expected=getattr(sp,function)(A*X+sp.pi/3)
    elif A==0:expected=sp.S.Zero
    else:
        power=A**(order%4)
        if function=='exp':expected=power*sp.exp(A*X+sp.pi/3)
        else:expected=power*sp.diff(getattr(sp,function)(A*sp.Symbol('t')+sp.pi/3),sp.Symbol('t'),order%4)/(A**(order%4))
        expected=expected.subs(sp.Symbol('t'),X)
    assert sp.expand(actual.to_sympy()-expected)==0


@pytest.mark.parametrize('slope',['2','2+I'])
@pytest.mark.parametrize('order',[2**32+1,2**128])
def test_huge_numeric_slope_powers_are_exact_and_export_without_materialization(slope,order):
    x=sx.Symbol('x');X=sp.Symbol('x');a=eval(slope,dict(I=sx.I));A=eval(slope,dict(I=sp.I))
    actual=sx.diff(sx.exp(a*x+1),x,order)
    coefficient,moving=actual.as_independent(x)
    assert coefficient==sx.Pow(a,order,evaluate=False)
    assert moving==sx.exp(a*x+1)
    expected=sp.Pow(A,order,evaluate=False)
    assert coefficient.to_sympy()==expected
    reference=actual.to_sympy()
    assert Counter(sp.Mul.make_args(reference))==Counter([expected,sp.exp(A*X+1)])
    assert pickle.loads(pickle.dumps(actual))==actual


@pytest.mark.parametrize('degree',[2**32-1,2**32,2**80])
@pytest.mark.parametrize('remaining',[0,1,3])
@pytest.mark.parametrize('base_source',['x','x+1','2*x+3'])
def test_large_polynomial_boundary_coefficients_are_held_exactly(degree,remaining,base_source):
    x=sx.Symbol('x');order=degree-remaining;base=eval(base_source,dict(x=x))
    actual=sx.diff(base**degree,x,order,simplify=False)
    assert not actual.is_Derivative
    assert not actual.has(sx.nan,sx.zoo,sx.oo)
    slope=2 if base_source=='2*x+3' else 1
    coefficient=sx.factorial(degree,evaluate=False) if remaining<=1 else sx.FallingFactorial(degree,order,evaluate=False)
    expected=[coefficient]
    if slope!=1:expected.append(sx.Pow(slope,order,evaluate=False))
    if remaining:expected.append(base**remaining)
    # Compare the COMPLETE factor multiset without activating the constructor
    # of a gigantic numeric power. Storage metadata is intentionally private.
    assert Counter(sx.Mul.make_args(actual))==Counter(expected)
    reference=actual.to_sympy()
    assert reference.has(sp.factorial(sp.Integer(degree),evaluate=False)) if remaining<=1 else reference.has(sp.FallingFactorial(degree,order,evaluate=False))
    assert pickle.loads(pickle.dumps(actual))==actual
    assert sx.diff(base**degree,x,degree+1)==0


def test_polynomial_roots_remain_regular_at_the_degree_boundary():
    x=sx.Symbol('x');m=2**80;base=x+1
    assert sx.diff(base**m,x,m-1,simplify=False).subs(x,-1)==0
    result=sx.diff(base**m,x,m,simplify=False)
    assert result==sx.factorial(m,evaluate=False)
    assert result.subs(x,-1)==result


@pytest.mark.parametrize('order',['Float(0)','Float(1)','Float(2**32)','Rational(1,2)','n'])
def test_nonliteral_counts_are_never_normalized_to_an_integer(order):
    x,n=sx.symbols('x n');count=eval(order,dict(Float=sx.Float,Rational=sx.Rational,n=n))
    held=sx.Derivative(sx.sin(x),(x,count))
    assert held._integer_order_value() is None
    assert held.doit()==held


@pytest.mark.parametrize('source',['a*x**2','f(1)*x','sin(x**2)','exp(x**2)','x**(-1)','oo*x**2','zoo*x'])
def test_unproved_domains_keep_the_original_big_order_fallback(source):
    x,a=sx.symbols('x a');f=sx.Function('f');body=eval(source,dict(x=x,a=a,f=f,sin=sx.sin,exp=sx.exp,oo=sx.oo,zoo=sx.zoo))
    held=sx.Derivative(body,(x,2**100))
    assert held._integer_order_value() is None
    assert held.doit()==held


def test_inexact_coefficients_preserve_the_existing_precision_fallback():
    x=sx.Symbol('x');body=sx.Float('1.234567890123456789',70)*x**5
    held=sx.Derivative(body,(x,2**100))
    assert held._integer_order_value() is None
    assert held.doit()==held
    assert held.variable_count[0][1]==2**100


def test_symbolic_formulas_stay_open_and_specializations_keep_the_count_protocol():
    x=sx.Symbol('x');n=sx.Symbol('n',integer=True,nonnegative=True)
    held=sx.Derivative(x,(x,n),evaluate=True)
    assert held.is_Derivative
    assert held.subs(n,sx.Rational(1,2)).is_Derivative
    assert held.subs(n,sx.Float(0)).is_Derivative
    with pytest.raises(ValueError):held.subs(n,-1)
    assert held.subs(n,0)==x
    assert held.subs(n,1).doit()==1


def test_machine_order_groups_that_exceed_u32_in_total_do_not_wrap():
    x=sx.Symbol('x');held=sx.Derivative(sx.Derivative(sx.sin(x),(x,2**32-1)),x)
    assert held.derivative_count==2**32
    assert held.doit()==sx.sin(x)


def test_known_nonfinite_variable_does_not_gain_a_new_literal_certificate():
    x=sx.Symbol('x',infinite=True);held=sx.Derivative(x**3,(x,2**100))
    assert held._integer_order_value() is None
    assert held.doit()==held


def test_host_derivative_hooks_are_left_on_the_established_route():
    x=sx.Symbol('x');calls=[]
    class F(sx.Function):
        def _eval_derivative(self,variable):calls.append(variable);return sx.Integer(7)
    assert sx.diff(F(x),x,3,simplify=False)==0
    assert calls
    calls.clear();held=sx.Derivative(F(x),(x,2**100))
    assert held._integer_order_value() is None
    assert held.doit()==held
    assert not calls


def test_export_retains_an_explicitly_held_numeric_power_in_a_product():
    x=sx.Symbol('x');order=2**32
    result=sx.diff(sx.exp(2*x),x,order).to_sympy()
    assert any(node==sp.Pow(2,order,evaluate=False) for node in sp.preorder_traversal(result))


@pytest.mark.parametrize('order',[-1,-2**100])
def test_invalid_orders_still_raise_before_a_closed_form_is_attempted(order):
    x=sx.Symbol('x')
    with pytest.raises(ValueError,match='order of differentiation must be nonnegative'):
        sx.Derivative(sx.sin(x),(x,order),evaluate=True)


def test_evaluate_false_retains_the_complete_literal_count_until_doit():
    x=sx.Symbol('x');n=2**128+1
    held=sx.Derivative(sx.sin(x),(x,n),evaluate=False)
    assert held.is_Derivative
    assert held.expr==sx.sin(x) and held.variable_count==((x,sx.Integer(n)),)
    assert held.doit()==sx.cos(x)
    assert pickle.loads(pickle.dumps(held))==held


def test_noncommutative_variable_keeps_the_existing_held_big_order_route():
    x=sx.Symbol('x',commutative=False)
    held=sx.Derivative(sx.sin(x),(x,2**100))
    assert held._integer_order_value() is None
    assert held.doit()==held


def test_entrywise_matrix_evaluation_and_multiple_exact_variable_groups():
    x,y=sx.symbols('x y');n=2**100+1
    matrix=sx.Matrix([[sx.sin(x),sx.exp(x)],[x**7,1]])
    expected=sx.Matrix([[sx.cos(x),sx.exp(x)],[0,0]])
    assert matrix.applyfunc(lambda entry:sx.diff(entry,x,n))==expected
    # The existing Matrix.diff count parser still accepts only machine orders;
    # this scalar tranche deliberately preserves its explicit refusal.
    with pytest.raises(ValueError,match='derivative order must be a nonnegative integer'):
        matrix.diff(x,n)
    held=sx.Derivative(sx.sin(x)*sx.exp(y),(x,n),(y,n),evaluate=False)
    # Parameters depending on another variable have no new coefficient
    # certificate: retain the complete original fallback, not a partial result.
    assert held._integer_order_value() is None
    assert held.doit()==held


@pytest.mark.parametrize('source',['sin(I*x)','cos(I*x)','sin(2*I*x)','cos(2*I*x)','sinh(2*x+1)','cosh(-2*x+1)'])
@pytest.mark.parametrize('order',[0,1,2,63,64,2**32-1,2**32,2**128+1])
@pytest.mark.parametrize('mode',['diff','doit'])
def test_canonical_hyperbolic_images_keep_the_complex_sine_cosine_domain(source,order,mode):
    x=sx.Symbol('x');X=sp.Symbol('x')
    env=lambda m,v:dict(x=v,I=m.I,sin=m.sin,cos=m.cos,sinh=m.sinh,cosh=m.cosh)
    body=eval(source,env(sx,x));reference=eval(source,env(sp,X))
    actual=run(body,x,order,mode)
    if order<=64:
        assert sp.expand(actual.to_sympy()-sp.diff(reference,X,order))==0
    else:
        # D^2 fixes the unit-slope hyperbolic images. For slope magnitude2,
        # replace the ENTIRE exact large power with a formal scalar before
        # doing any reference arithmetic, preventing integer materialization.
        remainder=order%2
        slope=2 if '2*' in source else 1
        oracle=sp.diff(reference,X,remainder)/(slope**remainder)
        converted=actual.to_sympy()
        if slope==2:
            A=sp.Symbol('A');power=sp.Pow(2,order,evaluate=False)
            assert converted.has(power)
            converted=converted.xreplace({power:A})
            oracle*=A
        assert sp.expand(converted-oracle)==0


@pytest.mark.parametrize('function',['sin','cos','exp'])
@pytest.mark.parametrize('order',[2**32-1,2**32,2**128+1])
@pytest.mark.parametrize('method',['subs','xreplace'])
def test_a_later_exact_count_specialization_evaluates_through_the_native_kernel(function,order,method):
    x,n=sx.symbols('x n');X=sp.Symbol('x')
    held=sx.Derivative(getattr(sx,function)(x),(x,n),evaluate=False)
    specialized=held.subs(n,order) if method=='subs' else held.xreplace({n:sx.Integer(order)})
    assert specialized.is_Derivative
    assert specialized.variable_count==((x,sx.Integer(order)),)
    actual=specialized.doit()
    expected=sp.exp(X) if function=='exp' else sp.diff(getattr(sp,function)(X),X,order%4)
    assert actual.to_sympy()==expected


def test_compact_machine_derivatives_remain_compact_in_validation_and_evalf():
    x,y=sx.symbols('x y');f=sx.Function('f');order=2**32-1
    held=sx.Derivative(f(x),(x,order),evaluate=False)
    assert sx.count_ops(held)==2
    assert held.evalf(subs={y:1})==held
    # An Order validator visits the derivative while replacing a symbol.
    assert (held+y).subs(y,x)==held+x
    assert sx.Matrix([[held]]).subs(y,x)==sx.Matrix([[held]])


@pytest.mark.parametrize('remaining',[1,3])
def test_held_numeric_slope_and_factorial_are_regular_at_the_polynomial_root(remaining):
    x=sx.Symbol('x');degree=2**80
    result=sx.diff((2*x+3)**degree,x,degree-remaining,simplify=False)
    assert result._has_deferred_numeric_power
    assert result.subs(x,-sx.Rational(3,2))==0


@pytest.mark.parametrize('order',[2**32-1,2**32,2**32+1,2**32+2,2**128+1])
@pytest.mark.parametrize('mode',['diff','doit','constructor','unsimplified'])
def test_negative_imaginary_unit_and_e_slopes_use_their_exact_canonical_images(order,mode):
    x=sx.Symbol('x')
    phase=[1,-sx.I,-1,sx.I][order%4]
    actual=run(sx.exp(-sx.I*x),x,order,mode)
    assert actual==phase*sx.exp(-sx.I*x)
    assert not actual._has_deferred_numeric_power
    actual=run(sx.exp(sx.E*x),x,order,mode)
    assert actual==sx.exp(order)*sx.exp(sx.E*x)
    assert actual.to_sympy()==sp.exp(order)*sp.exp(sp.E*sp.Symbol("x"))
    assert not actual._has_deferred_numeric_power


@pytest.mark.parametrize('degree',[4095,4096,4097,5000])
@pytest.mark.parametrize('remaining',[0,1,3])
def test_reasonable_factorial_and_falling_coefficients_still_materialize(degree,remaining):
    x=sx.Symbol('x');actual=sx.diff(x**degree,x,degree-remaining,simplify=False)
    expected=sx.Integer(int(sp.factorial(degree)/sp.factorial(remaining)))*x**remaining
    assert actual==expected
