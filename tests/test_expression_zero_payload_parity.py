"""EX's falsy Float payload survives only inside canonical nonzero DMP rows."""
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.polys.polyclasses import DMP
from sympy.polys.polyclasses import DMP as ReferenceDMP


def expression(value):
    if hasattr(value,'ex'):value=value.ex
    if not isinstance(value,(sx.Expr,sp.Basic)) and hasattr(value,'as_expr'):value=value.as_expr()
    return sx.to_sympy(value) if isinstance(value,sx.Expr) else sp.sympify(value)


def exact(value):
    e=expression(value)
    return e.xreplace({v:sp.Rational(v) for v in e.atoms(sp.Float)})


def scalar(actual,expected):
    a,b=expression(actual),expression(expected)
    assert exact(a)==exact(b)
    assert isinstance(a,sp.Float)==isinstance(b,sp.Float)
    if isinstance(b,sp.Float):assert a._prec==b._prec
    assert sorted(v._prec for v in a.atoms(sp.Float))==sorted(v._prec for v in b.atoms(sp.Float))


def nested(actual,expected):
    assert isinstance(actual,list)==isinstance(expected,list)
    if isinstance(expected,list):
        assert len(actual)==len(expected)
        for a,b in zip(actual,expected):nested(a,b)
    else:scalar(actual,expected)


def polynomial(actual,expected):
    assert str(actual.domain)==str(expected.domain)
    assert tuple(sx.to_sympy(g) for g in actual.gens)==expected.gens
    assert set(actual.as_dict())==set(expected.as_dict())
    for m,c in actual.as_dict().items():scalar(c,expected.as_dict()[m])
    scalar(actual.as_expr(),expected.as_expr())
    nested(actual.rep.to_sympy_list(),expected.rep.to_sympy_list())
    scalar(actual.LC(),expected.LC());scalar(actual.TC(),expected.TC())
    for m in [(0,)*len(actual.gens),(1,)*len(actual.gens)]:scalar(actual.nth(*m),expected.nth(*m))
    if len(actual.gens)==1:nested(actual.all_coeffs(),expected.all_coeffs())


def data(precision,kind,module):
    z=module.Float(0,precision=precision)
    choices=[{(5,):1,(0,):z},{(5,):1,(2,):z,(0,):3},
             {(7,):z,(0,):1},{(0,):z},
             {(2,2):1,(0,0):z},{(2,2):1,(2,0):z},
             {(2,2):1,(0,2):1,(0,0):z},
             {(2,2):1,(1,0):z,(0,0):1},
             {(2,2):1,(0,2):z,(0,0):1}]
    return choices[kind]


@pytest.mark.parametrize('precision',[1,53,139,140,211])
@pytest.mark.parametrize('kind',range(9))
def test_public_coefficients_raw_views_and_zero_evaluation(precision,kind):
    native=data(precision,kind,sx);reference=data(precision,kind,sp)
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    count=len(next(iter(native)));gens=(x,y)[:count];oracle_gens=(X,Y)[:count]
    p=sx.Poly.from_dict(native,*gens,domain=sx.EX)
    q=sp.Poly.from_dict(reference,*oracle_gens,domain=sp.EX)
    polynomial(p,q)
    for j in range(count):
        a=p.eval(gens[j],0);b=q.eval(oracle_gens[j],0)
        polynomial(a,b) if isinstance(b,sp.Poly) else scalar(a,b)
    for a,b in [(p.set_domain(sx.EX),q.set_domain(sp.EX)),(p.to_field(),q.to_field()),(p.to_ring(),q.to_ring()),(sx.Poly(p,*gens,domain=sx.EX),sp.Poly(q,*oracle_gens,domain=sp.EX))]:polynomial(a,b)


@pytest.mark.parametrize('precision',[1,53,139,140,211])
def test_lists_pickle_and_dmp_payload_protocols(precision):
    x=sx.Symbol('x');X=sp.Symbol('x');z=sx.Float(0,precision=precision);Z=sp.Float(0,precision=precision)
    p=sx.Poly.from_list([1,0,z],x,domain=sx.EX);q=sp.Poly.from_list([1,0,Z],X,domain=sp.EX)
    polynomial(p,q)
    for protocol in [0,pickle.HIGHEST_PROTOCOL]:
        # SymPy's Poly pickle also loses the EX domain. Preserve our established
        # stronger domain round trip and its coefficient payload deliberately.
        polynomial(pickle.loads(pickle.dumps(p,protocol=protocol)),q)
        nested(pickle.loads(pickle.dumps(p.rep,protocol=protocol)).to_sympy_list(),q.rep.to_sympy_list())
    for make in [lambda:DMP.from_sympy_list([1,0,z],0,sx.EX),lambda:DMP.from_dict({(2,):1,(0,):z},0,sx.EX),lambda:p.rep]:
        a=make();b=ReferenceDMP.from_sympy_list([1,0,Z],0,sp.EX)
        nested(a.to_sympy_list(),b.to_sympy_list());scalar(a.TC(),b.TC());scalar(a.nth(0),b.nth(0))
        assert set(a.to_dict())==set(b.to_dict())
        nested(a.neg().to_sympy_list(),b.neg().to_sympy_list())
        assert a.convert(sx.EX) is a


@pytest.mark.parametrize('precision',[1,53,139,140,211])
def test_leading_zero_does_not_force_huge_dense_allocation(precision):
    x=sx.Symbol('x');z=sx.Float(0,precision=precision)
    p=sx.Poly.from_dict({(1<<50,):z,(1<<28,):1},x,domain=sx.EX)
    assert p.degree()==1<<28 and p.as_dict()=={(1<<28,):1}
    assert p.eval(0)==0
    q=sx.Poly.from_dict({(1<<31,):z},x,domain=sx.EX)
    assert q.is_zero and q.LC()==q.TC()==q.eval(0)==0
    # Re-enter Dense with a retained Float zero, then remove it with the same
    # arithmetic that the reference uses to choose structural zeros.
    f=sx.Poly.from_dict({(5,):1,(0,):z},x,domain=sx.EX)
    X=sp.Symbol('x');g=sp.Poly.from_dict({(5,):1,(0,):sp.Float(0,precision=precision)},X,domain=sp.EX)
    for a,b in [(f+0,g+0),(f*1,g*1),(-f,-g),(f-f,g-g),(f.diff(),g.diff()),(f.integrate(),g.integrate())]:polynomial(a,b)


def test_input_callbacks_run_once_and_exceptions_remain_visible():
    x=sx.Symbol('x');z=sx.Float(0,precision=139);calls=[]
    class Input:
        def _sympy_(self):calls.append(1);return z
    p=sx.Poly.from_dict({(3,):1,(0,):Input()},x,domain=sx.EX)
    assert calls==[1];scalar(p.eval(0),sp.Float(0,precision=139))
    pickle.loads(pickle.dumps(p));p.rep.to_list();assert calls==[1]
    error=RuntimeError('payload input callback')
    class Bad:
        def _sympy_(self):raise error
    with pytest.raises(RuntimeError) as caught:sx.Poly.from_dict({(0,):Bad()},x,domain=sx.EX)
    assert caught.value is error


def test_payload_path_does_not_import_sympy_at_runtime():
    code='''import sys,importlib.abc,pickle
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,fullname,path=None,target=None):
  if fullname=='sympy' or fullname.startswith('sympy.'):raise AssertionError('runtime SymPy import')
sys.meta_path.insert(0,Block())
import symox as s
x=s.Symbol('x');z=s.Float(0,precision=139)
p=s.Poly.from_dict({(5,):1,(0,):z},x,domain=s.EX)
assert p.TC().is_Float and p.TC()._prec==139
assert p.eval(0).is_Float and p.eval(0)._prec==139
assert pickle.loads(pickle.dumps(p)).eval(0)._prec==139
assert p.rep.to_sympy_list()[-1]._prec==139
'''
    result=subprocess.run([sys.executable,'-B','-c',code],capture_output=True,text=True,timeout=20)
    assert result.returncode==0,result.stderr


def test_public_exact_zero_evaluation_uses_native_payload_without_expression_view(monkeypatch):
    from symox.polys import _evaluation
    x=sx.Symbol('x')
    p=sx.Poly.from_dict({(4000,):1,(0,):3},x,domain=sx.EX)
    def expression_view_was_called(*args,**kwargs):raise RuntimeError('expression view')
    monkeypatch.setattr(_evaluation,'_native_as_expr',expression_view_was_called)
    assert p.eval(0)==3 and p.eval([0])==3 and p.eval({x:0})==3
    with pytest.raises(RuntimeError,match='expression view'):p.eval(1)


@pytest.mark.parametrize('precision',[1,53,139,140,211])
def test_conversion_to_exact_domain_discards_only_the_converted_zero_payload(precision):
    x=sx.Symbol('x');X=sp.Symbol('x')
    p=sx.Poly.from_dict({(5,):1,(0,):sx.Float(0,precision=precision)},x,domain=sx.EX)
    q=sp.Poly.from_dict({(5,):1,(0,):sp.Float(0,precision=precision)},X,domain=sp.EX)
    polynomial(p.set_domain(sx.QQ),q.set_domain(sp.QQ))
    polynomial(p.set_domain(sx.QQ).set_domain(sx.EX),q.set_domain(sp.QQ).set_domain(sp.EX))
    plain=sx.Poly(x**5,x,domain=sx.EX);reference_plain=sp.Poly(X**5,X,domain=sp.EX)
    assert (p==plain)==(q==reference_plain)==False
    assert (p.rep==plain.rep)==(q.rep==reference_plain.rep)==False


def test_purepoly_pickle_keeps_class_and_coefficient_payload():
    x=sx.Symbol('x');X=sp.Symbol('x')
    p=sx.PurePoly.from_dict({(5,):1,(0,):sx.Float(0,precision=211)},x,domain=sx.EX)
    q=sp.PurePoly.from_dict({(5,):1,(0,):sp.Float(0,precision=211)},X,domain=sp.EX)
    assert type(p) is sx.PurePoly
    polynomial(p,q)
    restored=pickle.loads(pickle.dumps(p))
    assert type(restored) is sx.PurePoly
    polynomial(restored,q)


@pytest.mark.parametrize('precision',[1,53,139,140,211])
def test_three_variable_rows_and_partial_zero_evaluations(precision):
    x,y,z=sx.symbols('x y z');X,Y,Z=sp.symbols('x y z')
    native={(2,1,3):1,(2,1,0):sx.Float(0,precision=precision),(0,2,3):3,(0,2,0):sx.Float(0,precision=precision),(1,2,0):sx.Float(0,precision=precision)}
    oracle={m:expression(c) for m,c in native.items()}
    p=sx.Poly.from_dict(native,x,y,z,domain=sx.EX);q=sp.Poly.from_dict(oracle,X,Y,Z,domain=sp.EX)
    polynomial(p,q)
    for j in range(3):polynomial(p.eval(p.gens[j],0),q.eval(q.gens[j],0))


def test_public_as_expr_override_and_subclasses_keep_the_previous_fallback(monkeypatch):
    from symox.polys import _evaluation
    calls=[];original=_evaluation._native_as_expr
    def read_representation(poly):calls.append('internal');return original(poly)
    monkeypatch.setattr(_evaluation,'_native_as_expr',read_representation)
    x=sx.Symbol('x')
    class Custom(sx.Poly):
        def as_expr(self,*gens):calls.append('public');return sx.Integer(99)
    p=Custom.from_dict({(3,):1,(0,):sx.Float(0,precision=139)},x,domain=sx.EX)
    assert p.eval(0)==0
    assert calls==['internal']  # The old adapter used its captured native view.
    calls.clear()
    p=sx.Poly.from_dict({(3,):1,(0,):sx.Float(0,precision=139)},x,domain=sx.EX)
    def public(poly,*gens):calls.append('public');return sx.Integer(99)
    monkeypatch.setattr(sx.Poly,'as_expr',public)
    assert p.eval(0)==0 and calls==['internal']


def test_modified_domain_conversion_keeps_its_evaluation_callbacks(monkeypatch):
    x=sx.Symbol('x');p=sx.Poly.from_dict({(3,):1,(0,):sx.Float(0,precision=139)},x,domain=sx.EX)
    calls=[];domain_type=type(sx.EX);convert=domain_type.convert
    def observed(domain,value,*args,**kwargs):calls.append(value);return convert(domain,value,*args,**kwargs)
    monkeypatch.setattr(domain_type,'convert',observed)
    assert p.eval(0)==0
    assert len(calls)==2
    calls.clear()
    class Conversion:
        def __get__(self,instance,owner):
            calls.append('class' if instance is None else 'instance')
            if instance is None:return self
            return lambda value,*args,**kwargs:convert(instance,value,*args,**kwargs)
    monkeypatch.setattr(domain_type,'convert',Conversion())
    assert p.eval(0)==0
    assert calls==['instance','instance']


def test_custom_dictionary_views_keep_their_pickle_and_representation_protocols(monkeypatch):
    from symox.polys import _poly_methods
    x=sx.Symbol('x');calls=[]
    class Custom(sx.Poly):
        def as_dict(self,*args,**kwargs):calls.append('dict');return {(0,):7}
    p=Custom.from_dict({(3,):1,(0,):sx.Float(0,precision=139)},x,domain=sx.EX)
    restore,args=_poly_methods._reduce_polynomial(p)
    assert args[1]=={(0,):7} and calls==['dict']
    calls.clear()
    assert p.rep.to_sympy_list()==[7] and calls==['dict']
    calls.clear()
    p=sx.Poly.from_dict({(3,):1,(0,):sx.Float(0,precision=139)},x,domain=sx.EX)
    def dictionary(poly,*args,**kwargs):calls.append('dict');return {(0,):7}
    monkeypatch.setattr(sx.Poly,'as_dict',dictionary)
    restore,args=_poly_methods._reduce_polynomial(p)
    assert args[1]=={(0,):7} and calls==['dict']


def test_raw_dmp_constructor_keeps_unconverted_scalar_types():
    for value in [1,1.25]:
        p=DMP([value],sx.EX,0);q=ReferenceDMP([value],sp.EX,0)
        assert type(p.to_list()[0]) is type(q.to_list()[0]) is type(value)
        assert p.to_list()[0]==q.to_list()[0]==value


def test_observable_coefficients_and_list_iterators_keep_the_complete_fallback():
    from symox import _core
    x=sx.Symbol('x');c=sx.gamma(sx.Rational(1,3))
    assert _core.Poly._from_coefficient_terms({(1,):c},(x,),sx.EX) is None
    assert sx.Poly.from_dict({(1,):c},x,domain=sx.EX).as_expr()==c*x
    assert DMP.from_sympy_list([c,0],0,sx.EX).to_sympy_dict()=={(1,):c}
    calls=[]
    class Observed(list):
        def __iter__(self):calls.append('iter');return super().__iter__()
    DMP.from_sympy_list(Observed([1,0]),0,sx.EX)
    assert calls==['iter']
    calls.clear()
    class Input:
        def __init__(self,name,value):self.name,self.value=name,value
        def _sympy_(self):calls.append(self.name);return self.value
    p=DMP.from_list([Input('one',sx.Integer(1)),Input('zero',sx.Integer(0))],0,sx.EX)
    assert calls==['one','zero'] and set(p.to_dict())=={(1,)}
    # The old DMP holder accepts wider raw dictionaries independently of the
    # native Poly coordinate limit. Keep that established storage path.
    p=DMP.from_dict({(1<<32,):sx.EX.one},0,sx.EX)
    assert set(p.to_dict())=={(1<<32,)}
    for p in [DMP.from_dict({},1<<80,sx.EX),DMP.from_list([],1<<80,sx.EX)]:
        assert p.lev==1<<80 and p.is_zero


@pytest.mark.parametrize('precision',[1,53,139,140,211])
def test_nonzero_binary_coefficients_are_unchanged_next_to_zero_payloads(precision):
    x=sx.Symbol('x');X=sp.Symbol('x')
    c=sx.Float('1.234567890123456789012345678901234567890123456789',precision=precision)
    z=sx.Float(0,precision=precision)
    p=sx.Poly.from_dict({(3,):c,(1,):sx.Rational(2,3),(0,):z},x,domain=sx.EX)
    q=sp.Poly.from_dict({(3,):expression(c),(1,):sp.Rational(2,3),(0,):expression(z)},X,domain=sp.EX)
    polynomial(p,q)
    scalar(p.LC(),expression(c))
    a=DMP.from_sympy_list([c,sx.Rational(2,3),z],0,sx.EX)
    b=ReferenceDMP.from_sympy_list([expression(c),sp.Rational(2,3),expression(z)],0,sp.EX)
    nested(a.to_sympy_list(),b.to_sympy_list())


@pytest.mark.parametrize('precision',[1,53,139,140,211])
def test_numeric_zero_points_preserve_the_stored_coefficient_precision(precision):
    x=sx.Symbol('x');X=sp.Symbol('x')
    p=sx.Poly.from_dict({(5,):1,(0,):sx.Float(0,precision=precision)},x,domain=sx.EX)
    q=sp.Poly.from_dict({(5,):1,(0,):sp.Float(0,precision=precision)},X,domain=sp.EX)
    for point,reference_point in [(0,0),(0.0,0.0),(-0.0,-0.0),(sx.Float(0,precision=211),sp.Float(0,precision=211))]:
        scalar(p.eval(point),q.eval(reference_point))


@pytest.mark.parametrize('generator',[sx.Integer(0),sx.Integer(2),sx.Rational(1,2),sx.Float('1.25',precision=139)])
def test_zero_evaluation_uses_the_formal_generator_position(generator):
    p=sx.Poly.from_dict({(3,):1,(0,):sx.Float(0,precision=139)},generator,domain=sx.EX)
    q=sp.Poly.from_dict({(3,):1,(0,):sp.Float(0,precision=139)},expression(generator),domain=sp.EX)
    polynomial(p,q)
    scalar(p.eval(0),q.eval(0))
