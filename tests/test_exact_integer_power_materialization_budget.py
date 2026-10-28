"""Recipe-local exact powers; direct arithmetic keeps complete evaluation."""
import copy
import pickle
import pytest
import sympy as sp
import symox as sx


def recipe(base,order):
    x=sx.Symbol('x')
    return sx.diff(sx.exp(base*x),x,order).as_independent(x)[0]


@pytest.mark.parametrize('source',['2','-2','Rational(2,3)','Rational(-2,3)','Rational(2**200+1,2**200-1)'])
@pytest.mark.parametrize('order',[2**32,2**32+1,2**63-1,2**128+1])
def test_recipe_powers_preserve_exact_operands_and_implicit_reconstruction(source,order):
    x=sx.Symbol('x');base=eval(source,dict(Rational=sx.Rational))
    B=sp.sympify(eval(source,dict(Rational=sp.Rational)))
    power=recipe(base,order)
    expected=sx.Pow(abs(base),order,evaluate=False)
    if B<0 and order%2:expected=sx.Mul(-1,expected,evaluate=False)
    assert power==expected and power._has_deferred_numeric_power
    assert (power*x).subs(x,1)==power
    assert (power*x).subs(x,0)==0
    assert power/power==1
    expected_reference=sp.Pow(abs(B),order,evaluate=False)
    if B<0 and order%2:expected_reference=sp.Mul(-1,expected_reference,evaluate=False)
    assert power.to_sympy()==expected_reference
    for restored in [pickle.loads(pickle.dumps(power)),copy.copy(power),copy.deepcopy(power)]:
        assert restored==power and restored._has_deferred_numeric_power
        assert (restored*x).subs(x,0)==0


@pytest.mark.parametrize('source',['2','Rational(2,3)','2+I'])
@pytest.mark.parametrize('replacement',[0,1,-1])
def test_recipe_base_reconstruction_restores_exact_trivial_values(source,replacement):
    base=eval(source,dict(Rational=sx.Rational,I=sx.I));order=2**32+1
    value=recipe(base,order)
    for actual in [value.subs(base,replacement),value.xreplace({base:sx.Integer(replacement)})]:
        assert actual==replacement**order
        assert not actual._has_deferred_numeric_power


@pytest.mark.parametrize('replacement',[0,1,2])
def test_recipe_order_reconstruction_uses_the_exact_ordinary_small_power(replacement):
    order=2**32+1;value=recipe(sx.Integer(2),order)
    for actual in [value.subs(order,replacement),value.xreplace({sx.Integer(order):sx.Integer(replacement)})]:
        assert actual==2**replacement
        assert not actual._has_deferred_numeric_power


def test_recipe_base_becoming_float_keeps_the_precision_evaluator():
    # A unit Float makes this gigantic exponent safe to evaluate; it must not
    # be coerced back to a held exact-integer power.
    order=2**32;value=recipe(sx.Integer(2),order)
    assert value.subs(2,sx.Float(1,70))==sx.Float(1,70)
    assert not value.subs(2,sx.Float(1,70))._has_deferred_numeric_power


def test_full_materialization_is_still_available_by_direct_pow_int_and_doit():
    # Recipe storage is selected above its 1Mi-bit estimate; the actual
    # integer remains a reasonable 65KiB, so all complete routes are tested.
    order=524289;held=recipe(sx.Integer(2),order)
    assert held._has_deferred_numeric_power
    direct=sx.Pow(2,order)
    assert direct.is_Integer and not direct._has_deferred_numeric_power
    assert int(direct).bit_length()==order+1
    assert int(held)==int(direct)
    assert held.doit()==direct and held.doit(deep=False)==direct
    # Rebuilding the public constructor is an explicit evaluation request.
    assert held.func(*held.args)==direct


@pytest.mark.parametrize('base',[0,1,-1])
@pytest.mark.parametrize('order',[2**32,-2**32,2**63-1,-2**63])
def test_constant_size_integer_bases_keep_their_exact_evaluation(base,order):
    assert sx.Pow(base,order).to_sympy()==sp.Pow(base,order)


@pytest.mark.parametrize('source',['2','-2','Rational(2,3)','Rational(-2,3)','Rational(2**200+1,2**200-1)'])
@pytest.mark.parametrize('order',[-9,-2,-1,0,1,2,9])
def test_materializable_rational_powers_still_match_reference(source,order):
    base=eval(source,dict(Rational=sx.Rational));B=sp.sympify(eval(source,dict(Rational=sp.Rational)))
    assert sx.Pow(base,order).to_sympy()==B**order


def test_float_powers_keep_the_precision_evaluator():
    base=sx.Float('1.234567890123456789',70);B=sp.Float('1.234567890123456789',70)
    assert sx.Pow(base,-9).to_sympy()==B**(-9)
    assert sx.Pow(sx.Float(1),2**63-1)==sx.Float(1)


def test_expand_and_cse_preserve_complete_exact_recipe_coefficients():
    x,y=sx.symbols('x y');power=recipe(sx.Integer(2),2**32)
    original=power*(x+1)
    expanded=sx.expand(original)
    assert expanded==power*x+power
    assert expanded._has_deferred_numeric_power
    replacements,reduced=sx.cse([original,power*(y+1)])
    for index,value in enumerate(reduced):
        for symbol,body in reversed(replacements):value=value.xreplace({symbol:body})
        assert value==[original,power*(y+1)][index]
        assert value._has_deferred_numeric_power


@pytest.mark.parametrize('order',[2**32,2**32+1,2**128+1])
def test_exact_recipe_addition_products_and_reciprocals_do_not_change_values(order):
    power=recipe(sx.Integer(2),order);x=sx.Symbol('x')
    assert power-power==0
    assert power/power==1
    assert ((power+1)+(power+2))-2*power==3
    assert (power**2)/(power*power)==1
    inverse=1/power
    assert inverse*power==1
    assert (inverse*x).subs(x,0)==0
    assert pickle.loads(pickle.dumps(inverse))==inverse
    # Explicit ordinary arithmetic is complete even if it consumes a recipe.
    if order==2**32:
        small=recipe(sx.Integer(2),524289)
        assert (small*sx.Pow(2,524289,evaluate=False)).doit()==sx.Pow(2,1048578)


@pytest.mark.parametrize('customization',['constructor','func','matrix_dispatch','native_node','matrix_type','metaclass_call','constructor_code','power_hook'])
def test_cse_recipe_reconstruction_keeps_observable_custom_protocols(monkeypatch,customization):
    import importlib
    compat=importlib.import_module('symox._compat')
    module=importlib.import_module('symox.simplify.cse_main')
    core=importlib.import_module('symox._core')
    x,y=sx.symbols('x y');power=recipe(sx.Integer(2),524289)
    inputs=[power*(x+1),power*(y+1)];calls=[]
    if customization=='constructor':
        original=compat.Pow.__new__
        def constructor(cls,*args,**kwargs):
            calls.append('constructor');return original(cls,*args,**kwargs)
        monkeypatch.setattr(compat.Pow,'__new__',staticmethod(constructor))
    elif customization=='func':
        original=core.Expr.func
        def head(self):
            calls.append('func');return original.__get__(self)
        monkeypatch.setattr(type(power),'func',property(head))
    elif customization=='matrix_dispatch':
        original=compat._matrix_operand
        def matrix(value):
            calls.append('matrix');return original(value)
        monkeypatch.setattr(compat,'_matrix_operand',matrix)
    elif customization=='native_node':
        original=core._node
        def node(*args,**kwargs):
            calls.append('node');return original(*args,**kwargs)
        monkeypatch.setattr(core,'_node',node)
    elif customization=='matrix_type':
        class MatrixMeta(type):
            def __instancecheck__(cls,value):
                calls.append('matrix_type');return False
        monkeypatch.setattr(core,'Matrix',MatrixMeta('Matrix',(),{}))
    elif customization=='metaclass_call':
        meta=type(compat.Pow);original=meta.__call__
        def call(cls,*args,**kwargs):
            calls.append('call');return original(cls,*args,**kwargs)
        monkeypatch.setattr(meta,'__call__',call)
    elif customization=='constructor_code':
        original=compat.Pow.__new__
        # Mutate the actual captured function, preserving its object identity.
        monkeypatch.setitem(original.__globals__,'_recipe_protocol_calls',calls)
        def new_code(cls,b,e,evaluate=None):
            _recipe_protocol_calls.append('code')
            return _core._node('Pow',(b,e),True if evaluate is None else evaluate)
        monkeypatch.setattr(original,'__code__',new_code.__code__)
    else:
        kind=type(sx.Integer(2))
        monkeypatch.setattr(kind,'_symox_head',True,raising=False)
        def power_hook(self,e):
            calls.append('power_hook');return core.Pow(self,e)
        monkeypatch.setattr(kind,'_eval_power',power_hook,raising=False)
        host=importlib.import_module('symox._host');original=host.eval_power
        def dispatch(*args):
            calls.append('power_dispatch');return original(*args)
        monkeypatch.setattr(host,'eval_power',dispatch)
    assert not module._stock_recipe_reconstruction(power,list(power.args))
    actual=sx.cse(inputs);observed=list(calls);calls.clear()
    # The established CSE route is the authority for which callbacks occur.
    # Do not assume an operator is called merely because it was overridden.
    monkeypatch.setattr(module,'_stock_recipe_reconstruction',lambda *args:False)
    expected=sx.cse(inputs)
    assert actual==expected
    assert observed==calls and observed


@pytest.mark.parametrize('order',[524289,2**32+1,2**128+1])
def test_recipe_metadata_is_transparent_to_public_expression_identity(order):
    marked=recipe(sx.Integer(2),order)
    ordinary=sx.Pow(2,order,evaluate=False)
    assert marked==ordinary and hash(marked)==hash(ordinary)
    assert {ordinary:'ordinary'}[marked]=='ordinary'
    assert marked.func is ordinary.func and marked.args==ordinary.args
    assert marked.free_symbols==ordinary.free_symbols==set()
    assert sx.count_ops(marked)==sx.count_ops(ordinary)==1
    assert str(marked)==str(ordinary)
    assert sx.srepr(marked)==sx.srepr(ordinary)
    assert marked.to_sympy()==ordinary.to_sympy()


@pytest.mark.parametrize('replacement',[0,1,-1])
def test_reciprocal_recipe_reconstruction_keeps_the_ordinary_zero_pole(replacement):
    order=2**32+1;inverse=1/recipe(sx.Integer(2),order)
    actual=inverse.subs(2,replacement)
    expected=sx.zoo if replacement==0 else sx.Integer(replacement)
    assert actual==expected
    assert not actual._has_deferred_numeric_power


@pytest.mark.parametrize('reverse',[False,True])
def test_equal_raw_and_marked_sum_terms_keep_recipe_storage(reverse):
    marked=recipe(sx.Integer(2),2**128+1)
    raw=sx.Pow(2,2**128+1,evaluate=False)
    left,right=(raw+1,marked+2) if reverse else (marked+2,raw+1)
    combined=left+right
    assert combined==2*marked+3 and combined._has_deferred_numeric_power


@pytest.mark.parametrize('immutable',[False,True])
def test_matrix_pickle_keeps_each_recipe_and_raw_owner(immutable):
    x=sx.Symbol('x');power=recipe(sx.Integer(2),2**128+1)
    raw=sx.Pow(2,2**128+1,evaluate=False)
    cls=sx.ImmutableMatrix if immutable else sx.Matrix
    original=cls([[power,raw],[power*x,0]])
    restored=pickle.loads(pickle.dumps(original))
    assert type(restored) is type(original) and restored==original
    assert restored[0,0]._has_deferred_numeric_power
    assert not restored[0,1]._has_deferred_numeric_power
    assert restored.subs(x,0)[1,0]==0
