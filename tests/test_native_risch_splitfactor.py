"""Exact QQ splitfactor, independent certificates, and complete old protocols."""
from fractions import Fraction
import importlib
import sys

import pytest
import symox as sx
import sympy as sp

risch=importlib.import_module('symox.integrals.risch')
selector=importlib.import_module('symox.integrals._risch_splitfactor_native')
reference=importlib.import_module('sympy.integrals.risch')


def fixture(seed,package,module):
    variables=list(package.symbols('split_x split_t split_u'))[:1+seed%3]
    active=len(variables)-(1 if seed%4==0 and len(variables)>1 else 0)
    t=variables[active-1];mode=seed%5
    special=[package.S.One,t,t*t+1,t*(t-1),t*t+1][mode]
    dt=[package.S.One,t,t*t+1,2*t*(t-1),(t*t+1)/2][mode]
    scale=[package.S.One,-6,6,package.Rational(2,3),package.Rational(-5,7)][seed%5]
    expr=scale*(t+2)**((seed//5)%5)*special**((seed//7)%5)
    if seed%17==0:expr=package.S.Zero
    ds=[package.Poly(dt if i==active-1 else 1,v)for i,v in enumerate(variables)]
    de=module.DifferentialExtension(extension={'D':ds})
    if active<len(variables):de.decrement_level()
    return package.Poly(expr,t,domain='QQ'),de


def signature(pair):
    return tuple((tuple(t.name for t in p.gens),str(p.domain),tuple(sorted(
        (monomial,Fraction(int(v.p),int(v.q)))for monomial,v in p.as_dict().items())))for p in pair)


def trim(p):
    while p and not p[-1]:p.pop()
    return p


def dense(p):
    d={m[0]:Fraction(int(v.p),int(v.q))for m,v in p.as_dict().items()}
    return trim([d.get(i,Fraction())for i in range(max(d,default=-1)+1)])


def mul(a,b):
    if not a or not b:return []
    p=[Fraction()]*(len(a)+len(b)-1)
    for i,v in enumerate(a):
        for j,w in enumerate(b):p[i+j]+=v*w
    return trim(p)


def add(a,b):
    return trim([(a[i]if i<len(a)else Fraction())+(b[i]if i<len(b)else Fraction())
                 for i in range(max(len(a),len(b)))])


def sub(a,b):return add(a,[-v for v in b])


def diff(p):return trim([i*p[i]for i in range(1,len(p))])


def divrem(a,b):
    assert b
    a=list(a);q=[Fraction()]*max(0,len(a)-len(b)+1)
    while a and len(a)>=len(b):
        k=len(a)-len(b);v=a[-1]/b[-1];q[k]=v
        for i,w in enumerate(b):a[i+k]-=v*w
        trim(a)
    return trim(q),a


def gcd(a,b):
    while b:a,b=b,divrem(a,b)[1]
    return [v/a[-1]for v in a]if a else []


def certificate(p,de,pair):
    p,tau=dense(p),dense(de.d);normal,special=map(dense,pair)
    assert mul(normal,special)==p
    assert divrem(mul(diff(special),tau),special)[1]==[]
    if normal:
        sq,rem=divrem(normal,gcd(normal,diff(normal)))
        assert rem==[] and gcd(sq,mul(diff(sq),tau))==[Fraction(1)]


def trace(call):
    kernel=sx._core._risch_splitfactor;calls=[];previous=sys.getprofile()
    def profile(frame,event,arg):
        if event=='c_call' and arg is kernel:calls.append('native')
    sys.setprofile(profile)
    try:return call(),calls
    finally:sys.setprofile(previous)


def old(monkeypatch,p,de,**kwargs):
    with monkeypatch.context()as patch:
        patch.setattr(risch,'_native_splitfactor',lambda *args:None)
        return risch.splitfactor(p,de,**kwargs)


@pytest.mark.parametrize('seed',range(96))
def test_selected_splitfactor_matches_old_reference_and_independent_certificates(monkeypatch,seed):
    p,de=fixture(seed,sx,risch);rp,rde=fixture(seed,sp,reference)
    before=(de.level,de.t,de.d,tuple(de.T),tuple(de.D),tuple(p.as_dict().items()))
    actual,calls=trace(lambda:risch.splitfactor(p,de))
    assert calls==['native']
    selected=selector.reduce(p,de,sx.Poly(1,de.t,domain=p.get_domain()),False,None)
    assert selected is not None and signature(selected)==signature(actual)
    expected=old(monkeypatch,p,de)
    assert signature(actual)==signature(expected)==signature(reference.splitfactor(rp,rde))
    assert (actual[0]is p)==(expected[0]is p)
    certificate(p,de,actual)
    assert before==(de.level,de.t,de.d,tuple(de.T),tuple(de.D),tuple(p.as_dict().items()))
    assert all(a is b for a,b in zip(before[3],de.T))
    assert all(a is b for a,b in zip(before[4],de.D))


@pytest.mark.parametrize('levels',[1,2,3])
@pytest.mark.parametrize('value',[6,-6,sx.Rational(2,3),sx.Rational(-5,7)])
def test_constant_content_keeps_domains_units_representatives_and_product(monkeypatch,levels,value):
    variables=list(sx.symbols('split_x split_t split_u'))[:levels];t=variables[-1]
    de=risch.DifferentialExtension(extension={'D':[sx.Poly(1,v)for v in variables]})
    p=sx.Poly(value,t,domain='QQ');actual,calls=trace(lambda:risch.splitfactor(p,de))
    selected=selector.reduce(p,de,sx.Poly(1,t,domain='QQ'),False,None)
    assert selected is not None and signature(selected)==signature(actual)
    expected=old(monkeypatch,p,de)
    assert calls==['native'] and signature(actual)==signature(expected)
    assert actual[0]is not p and str(actual[0].domain)=='QQ'
    assert str(actual[1].domain)==('QQ'if levels==1 else 'ZZ')
    if levels>1 and (isinstance(value,int)or value.is_Integer):assert actual[1].as_expr()==abs(value)
    certificate(p,de,actual)


@pytest.mark.parametrize('coefficientD,z',[(True,None),(False,0),(False,False),(False,sx.Symbol('split_z'))])
def test_other_options_preserve_complete_original_computation(monkeypatch,coefficientD,z):
    p,de=fixture(22,sx,risch);kwargs={'coefficientD':coefficientD,'z':z}
    assert selector.reduce(p,de,sx.Poly(1,de.t,domain=p.get_domain()),coefficientD,z)is None
    actual,calls=trace(lambda:risch.splitfactor(p,de,**kwargs))
    assert signature(actual)==signature(old(monkeypatch,p,de,**kwargs))
    if coefficientD:assert calls==[]
    # The old recursive call intentionally drops z. Its canonical QQ
    # subproblem can independently select the exact native algorithm.


@pytest.mark.parametrize('domain',['ZZ','RR','EX','GF(7)','QQ[y]'])
def test_other_domains_keep_values_domains_precision_and_errors(monkeypatch,domain):
    t,y=sx.symbols('split_t y')
    expr=t+1+(y if domain=='QQ[y]'else 0)
    if domain=='RR':expr=t+sx.Float('1.234567890123456789',80)
    p=sx.Poly(expr,t,domain=domain)
    de=risch.DifferentialExtension(extension={'D':[sx.Poly(1,t)]})
    assert selector.reduce(p,de,sx.Poly(1,de.t,domain=p.get_domain()),False,None)is None
    try:actual,calls=trace(lambda:risch.splitfactor(p,de))
    except Exception as actual:
        with pytest.raises(type(actual))as expected:old(monkeypatch,p,de)
        assert str(actual)==str(expected.value)
    else:
        expected=old(monkeypatch,p,de)
        assert actual==expected
        assert all(a.gens==b.gens and str(a.domain)==str(b.domain)and a.as_dict()==b.as_dict()
                   for a,b in zip(actual,expected))


@pytest.mark.parametrize('name',['derivation','Poly','gcd','exquo','to_field','as_poly','diff','get_domain','has'])
def test_consumed_callbacks_match_complete_old_body(monkeypatch,name):
    p,de=fixture(22,sx,risch);events=[]
    owner=risch if name in ('Poly','derivation')else sx.Expr if name=='has'else sx.Poly
    original=getattr(owner,name)
    def callback(*args,**kwargs):events.append(name);return original(*args,**kwargs)
    monkeypatch.setattr(owner,name,callback)
    actual,calls=trace(lambda:risch.splitfactor(p,de));observed=events[:];events.clear()
    expected=old(monkeypatch,p,de)
    assert calls==[] and signature(actual)==signature(expected)and events==observed and events


@pytest.mark.parametrize('mode',['subclass','class_lookup','slot_property'])
def test_real_extension_reads_are_preserved(monkeypatch,mode):
    p,de=fixture(22,sx,risch);events=[]
    if mode=='subclass':
        class Custom(risch.DifferentialExtension):
            def __getattribute__(self,name):events.append(name);return super().__getattribute__(name)
        de=Custom(extension={'D':de.D})
    elif mode=='class_lookup':
        def callback(self,name):events.append(name);return object.__getattribute__(self,name)
        monkeypatch.setattr(risch.DifferentialExtension,'__getattribute__',callback)
    else:
        descriptor=risch.DifferentialExtension.t
        def callback(self):events.append('t');return descriptor.__get__(self,type(self))
        monkeypatch.setattr(risch.DifferentialExtension,'t',property(callback))
    events.clear();actual,calls=trace(lambda:risch.splitfactor(p,de));observed=events[:];events.clear()
    expected=old(monkeypatch,p,de)
    assert calls==[] and signature(actual)==signature(expected)and events==observed and events


@pytest.mark.parametrize('name',['_risch_splitfactor','_prove_fields'])
def test_foreign_native_or_shared_proof_delegate_is_never_called(monkeypatch,name):
    p,de=fixture(22,sx,risch);events=[]
    owner=sx._core if name=='_risch_splitfactor'else selector._shared
    def callback(*args):events.append(name);raise RuntimeError('foreign delegate')
    monkeypatch.setattr(owner,name,callback)
    actual=risch.splitfactor(p,de)
    assert signature(actual)==signature(old(monkeypatch,p,de))and events==[]
    if name=='_prove_fields':
        actual=risch.polynomial_reduce(p,de)
        with monkeypatch.context()as patch:
            patch.setattr(risch,'_native_polynomial_reduce',lambda *args:None)
            expected=risch.polynomial_reduce(p,de)
        assert signature(actual)==signature(expected)and events==[]


@pytest.mark.parametrize('mode',['code','defaults','foreign_globals','module_class'])
def test_shared_proof_origin_and_function_state_are_authenticated(monkeypatch,mode):
    from types import FunctionType,ModuleType
    p,de=fixture(22,sx,risch);events=[];shared=selector._shared
    if mode=='code':monkeypatch.setattr(shared._prove_fields,'__code__',(lambda p,de:None).__code__)
    elif mode=='defaults':monkeypatch.setattr(shared._prove_fields,'__defaults__',(None,))
    elif mode=='foreign_globals':
        function=shared._prove_fields
        monkeypatch.setattr(shared,'_prove_fields',FunctionType(function.__code__,dict(function.__globals__)))
    else:
        class Hooks(ModuleType):
            def __getattribute__(self,name):
                events.append(name);raise RuntimeError('foreign module lookup')
        monkeypatch.setattr(shared,'__class__',Hooks)
    actual,calls=trace(lambda:risch.splitfactor(p,de))
    assert calls==[] and signature(actual)==signature(old(monkeypatch,p,de))and events==[]


def test_missing_unused_d_member_still_runs_the_complete_old_algorithm(monkeypatch):
    t=sx.Symbol('split_t');p=sx.Poly(t*t+1,t,domain='QQ')
    de=risch.DifferentialExtension.__new__(risch.DifferentialExtension)
    de.T=[t];de.D=[sx.Poly(t*t+1,t)];de.t=t;de.level=-1
    actual,calls=trace(lambda:risch.splitfactor(p,de))
    assert calls==[] and signature(actual)==signature(old(monkeypatch,p,de))


@pytest.mark.parametrize('mode',['degree','bits','zero_derivation','tuple_variables','tuple_derivations','duplicate_variables','host_variable'])
def test_nonadmitted_payloads_retain_available_original_computations(monkeypatch,mode):
    t=sx.Symbol('split_t');p=sx.Poly(t*t+1,t,domain='QQ')
    de=risch.DifferentialExtension(extension={'D':[sx.Poly(t*t+1,t)]})
    if mode=='degree':p=sx.Poly((t+2)**65,t,domain='QQ')
    elif mode=='bits':p=sx.Poly((2**4100)*(t+2),t,domain='QQ')
    elif mode=='zero_derivation':de.D[0]=de.d=sx.Poly(0,t)
    elif mode=='tuple_variables':de.T=tuple(de.T)
    elif mode=='tuple_derivations':de.D=tuple(de.D)
    elif mode=='duplicate_variables':
        de.T=[t,t];de.D=[sx.Poly(1,t),de.d]
    else:
        class Host(sx.Symbol):pass
        t=Host('split_host_t');p=sx.Poly(t*t+1,t,domain='QQ')
        de=risch.DifferentialExtension(extension={'D':[sx.Poly(t*t+1,t)]})
    assert selector.reduce(p,de,sx.Poly(1,de.t,domain=p.get_domain()),False,None)is None
    actual=risch.splitfactor(p,de)
    assert signature(actual)==signature(old(monkeypatch,p,de))


def test_late_import_hook_observes_the_complete_original_import_log(monkeypatch):
    import builtins
    p,de=fixture(22,sx,risch);events=[];original=builtins.__import__
    def callback(name,globals=None,locals=None,fromlist=(),level=0):
        events.append((name,fromlist,level))
        return original(name,globals,locals,fromlist,level)
    monkeypatch.setattr(builtins,'__import__',callback)
    actual=risch.splitfactor(p,de);observed=events[:];events.clear()
    with monkeypatch.context()as patch:
        patch.setattr(risch,'_native_splitfactor',lambda *args:None)
        events.clear()
        expected=risch.splitfactor(p,de)
        legacy=events[:];events.clear()
    assert signature(actual)==signature(expected)and observed==legacy and legacy


def test_real_constructor_slot_is_proved_after_its_dictionary_is_restored():
    import os,subprocess
    code="""
import symox as s
from symox.integrals import risch
t=s.Symbol('split_t');p=s.Poly((t*t+1)**2,t,domain='QQ');one=s.Poly(1,t,domain='QQ')
de=risch.DifferentialExtension(extension={'D':[s.Poly(t*t+1,t)]})
original=s.Poly.__new__
def changed(cls,*args,**kwargs):return original(cls,*args,**kwargs)
s.Poly.__new__=staticmethod(changed)
s.Poly.__new__=original
assert s._core._risch_splitfactor(p,de.d,de.t,de.T,de.D,de.level,one)is None
def prefix(*args,**kwargs):
 risch.Poly=s.Poly
 return one
errors=[]
for mode in range(2):
 risch.Poly=prefix
 if mode:risch._native_splitfactor=lambda *args:None
 try:risch.splitfactor(p,de)
 except TypeError as error:errors.append(str(error))
 else:raise AssertionError('changed real constructor slot was omitted')
assert errors[0]==errors[1]
"""
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())


@pytest.mark.parametrize('levels',[1,2,3])
@pytest.mark.parametrize('caller',['canonical_representation','hermite_reduce'])
def test_actual_risch_callers_select_splitfactor_with_exact_decomposition(monkeypatch,levels,caller):
    def inputs(package,module):
        variables=list(package.symbols('split_x split_t split_u'))[:levels];t=variables[-1]
        de=module.DifferentialExtension(extension={'D':[package.Poly(1,v)for v in variables[:-1]]+[package.Poly(t*t+1,t)]})
        return package.Poly(t+1,t,domain='QQ'),package.Poly((t+2)**3*(t*t+1)**2,t,domain='QQ'),de
    a,d,de=inputs(sx,risch);ra,rd,rde=inputs(sp,reference)
    actual,calls=trace(lambda:getattr(risch,caller)(a,d,de))
    assert calls
    with monkeypatch.context()as patch:
        patch.setattr(risch,'_native_splitfactor',lambda *args:None)
        expected=getattr(risch,caller)(a,d,de)
    ref=getattr(reference,caller)(ra,rd,rde)
    flatten=lambda result:(result[0],*result[1],*result[2])if caller=='canonical_representation'else tuple(p for pair in result for p in pair)
    actual_signature=signature(flatten(actual));old_signature=signature(flatten(expected))
    reference_signature=signature(flatten(ref))
    assert actual_signature==old_signature
    if caller=='canonical_representation':assert actual_signature==reference_signature
    else:
        # Existing Poly.cancel domain inference in Symox returns this integer
        # g pair over ZZ, while SymPy retains QQ. Preserve that representative;
        # the baseline witness is archived separately, rather than weakening
        # the splitfactor domain comparison or claiming full Hermite parity.
        assert tuple(row[1]for row in actual_signature[:2])==('ZZ','ZZ')
        assert tuple(row[1]for row in reference_signature[:2])==('QQ','QQ')
        assert tuple((row[0],row[2])for row in actual_signature[:2])==tuple((row[0],row[2])for row in reference_signature[:2])
        assert actual_signature[2:]==reference_signature[2:]
    A,D,tau=dense(a),dense(d),dense(de.d)
    if caller=='canonical_representation':
        q,(b,ds),(c,dn)=actual
        Q,B,DS,C,DN=map(dense,(q,b,ds,c,dn))
        assert mul(A,mul(DS,DN))==mul(D,add(mul(Q,mul(DS,DN)),add(mul(B,DN),mul(C,DS))))
    else:
        (ga,gd),(ha,hd),(rr,dr)=actual
        GA,GD,HA,HD,R,DR=map(dense,(ga,gd,ha,hd,rr,dr))
        GD2=mul(GD,GD);DG=mul(sub(mul(diff(GA),GD),mul(GA,diff(GD))),tau)
        assert mul(A,mul(GD2,mul(HD,DR)))==add(mul(D,mul(DG,mul(HD,DR))),
            add(mul(D,mul(HA,mul(GD2,DR))),mul(D,mul(R,mul(GD2,HD)))))


@pytest.mark.parametrize('order',[(),('symox.diffgeom','symox.vector','symox.holonomic'),
    ('symox.diffgeom','symox.holonomic','symox.vector'),('symox.vector','symox.diffgeom','symox.holonomic'),
    ('symox.vector','symox.holonomic','symox.diffgeom'),('symox.holonomic','symox.diffgeom','symox.vector'),
    ('symox.holonomic','symox.vector','symox.diffgeom')])
def test_real_native_splitfactor_keeps_audit_log_empty_after_inputs(order):
    import os,subprocess
    code='''
import importlib,sys
import symox as s
from symox.integrals import risch
for name in ORDER:importlib.import_module(name)
selector=importlib.import_module('symox.integrals._risch_splitfactor_native')
x,t=s.symbols('audit_split_x audit_split_t')
de=risch.DifferentialExtension(extension={'D':[s.Poly(1,x),s.Poly(t*t+1,t)]})
p=s.Poly((t*t+1)**2*(t+2),t,domain='QQ');one=s.Poly(1,t,domain='QQ')
events=[]
def reject(event,args):
 if event=='object.__getattr__':
  events.append(args[1]);raise RuntimeError('unexpected metadata audit')
sys.addaudithook(reject)
selected=selector.reduce(p,de,one,False,None)
assert selected is not None and events==[]
actual=risch.splitfactor(p,de)
risch._native_splitfactor=lambda *args:None
expected=risch.splitfactor(p,de)
assert actual==expected and events==[]
assert [str(v.domain)for v in actual]==[str(v.domain)for v in expected]
'''.replace('ORDER',repr(order))
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())


@pytest.mark.parametrize('mode',['alias','module','collision'])
def test_c_metadata_reader_declines_foreign_origin_without_calls(monkeypatch,mode):
    from types import ModuleType
    p,de=fixture(22,sx,risch);one=sx.Poly(1,de.t,domain=p.domain);events=[]
    original=sx._core._ddm_function_metadata;key=None
    def foreign(*args):events.append('reader');raise RuntimeError('foreign reader')
    if mode=='alias':monkeypatch.setattr(sx._core,'_ddm_function_metadata',foreign)
    elif mode=='module':
        class Hooks(ModuleType):
            def __getattribute__(self,name):events.append(name);raise RuntimeError('module callback')
        monkeypatch.setattr(sx._core,'__class__',Hooks)
    else:
        class Key:
            def __hash__(self):return hash('_ddm_function_metadata')
            def __eq__(self,other):events.append('eq');raise RuntimeError('foreign key equality')
        monkeypatch.delattr(sx._core,'_ddm_function_metadata')
        key=Key();sx._core.__dict__[key]=None
    try:
        assert selector.reduce(p,de,one,False,None)is None and events==[]
        def outcome(call):
            try:return ('value',signature(call()))
            except Exception as error:return ('error',type(error),str(error))
        actual=outcome(lambda:risch.splitfactor(p,de));observed=events[:];events.clear()
        expected=outcome(lambda:old(monkeypatch,p,de))
        assert actual==expected and events==observed
        if mode=='module':assert observed
        else:assert observed==[]
    finally:
        if key is not None:del sx._core.__dict__[key]


@pytest.mark.parametrize('mode',['closure','defaults','kwdefaults','attributes','f_builtins'])
def test_metadata_payloads_are_proved_before_any_python_container_hook(mode):
    from types import FunctionType
    events=[];reader=sx._core._ddm_function_metadata;builtins=selector._ROOT['builtins']
    value=object()
    def factory():
        keep=value
        def function():return keep
        return function
    original=factory();code=original.__code__;scope=original.__globals__
    class Tuple(tuple):
        def __iter__(self):events.append('iter');raise RuntimeError('foreign tuple iteration')
        def __len__(self):events.append('len');raise RuntimeError('foreign tuple length')
        def __getitem__(self,index):events.append('getitem');raise RuntimeError('foreign tuple indexing')
    class Dict(dict):
        def __iter__(self):events.append('iter');raise RuntimeError('foreign dictionary iteration')
        def __len__(self):events.append('len');raise RuntimeError('foreign dictionary length')
    closure=Tuple(original.__closure__)if mode=='closure'else original.__closure__
    if mode=='f_builtins':
        saved=scope['__builtins__'];scope['__builtins__']=Dict(builtins)
        try:function=FunctionType(code,scope,closure=closure)
        finally:scope['__builtins__']=saved
    else:function=FunctionType(code,scope,closure=closure)
    defaults=Tuple((None,))if mode=='defaults'else None
    keywords=Dict(key=None)if mode=='kwdefaults'else None
    function.__defaults__=defaults;function.__kwdefaults__=keywords
    if mode=='attributes':function.__dict__=Dict()
    states=((function,code,scope,defaults,keywords,(value,),()),)
    assert not selector._states_match(states,reader,builtins)
    assert events==[]


@pytest.mark.parametrize('name',['derivation','Poly','gcd','exquo','to_field','as_poly','diff','get_domain','has'])
def test_raising_consumed_callbacks_keep_exact_old_error_and_journal(monkeypatch,name):
    p,de=fixture(22,sx,risch);events=[]
    owner=risch if name in ('Poly','derivation')else sx.Expr if name=='has'else sx.Poly
    def callback(*args,**kwargs):events.append(name);raise RuntimeError('splitfactor '+name)
    monkeypatch.setattr(owner,name,callback)
    with pytest.raises(RuntimeError)as actual:risch.splitfactor(p,de)
    observed=events[:];events.clear()
    with pytest.raises(RuntimeError)as expected:old(monkeypatch,p,de)
    assert str(actual.value)==str(expected.value)and events==observed and events


def test_preimport_reader_wrapper_is_not_blessed():
    import os,subprocess
    code='''
import builtins,sys
original=builtins.__import__;events=[]
def foreign(*args):events.append('reader');raise RuntimeError('foreign metadata')
def importing(name,globals=None,locals=None,fromlist=(),level=0):
 result=original(name,globals,locals,fromlist,level)
 if name=='_core' and globals and globals.get('__name__')=='symox.integrals._risch_splitfactor_native':
  sys.modules['symox._core']._ddm_function_metadata=foreign
 return result
builtins.__import__=importing
import symox as s
from symox.integrals import risch
selector=sys.modules['symox.integrals._risch_splitfactor_native']
t=s.Symbol('split_t');p=s.Poly((t*t+1)**2,t,domain='QQ')
de=risch.DifferentialExtension(extension={'D':[s.Poly(t*t+1,t)]})
events.clear()
assert selector._ROOT is None
actual=risch.splitfactor(p,de);risch._native_splitfactor=lambda *args:None
expected=risch.splitfactor(p,de)
assert actual==expected and events==[]
'''
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())


@pytest.mark.parametrize('owner_name,name',[('shared','_prove_fields'),('selector','_metadata'),
    ('rational','_stock_ratint_protocol_with_metadata')])
def test_consumed_namespace_keys_are_checked_before_collision_lookup(monkeypatch,owner_name,name):
    p,de=fixture(22,sx,risch);one=sx.Poly(1,de.t,domain=p.domain);events=[]
    owner={'shared':selector._shared,'selector':selector,
        'rational':importlib.import_module('symox.integrals.rationaltools')}[owner_name]
    class Key:
        def __hash__(self):return hash(name)
        def __eq__(self,other):events.append('eq');raise RuntimeError('foreign namespace equality')
    monkeypatch.delattr(owner,name);key=Key();owner.__dict__[key]=None
    try:
        assert selector.reduce(p,de,one,False,None)is None and events==[]
        actual=risch.splitfactor(p,de)
        assert signature(actual)==signature(old(monkeypatch,p,de))and events==[]
    finally:del owner.__dict__[key]


def test_stock_rationaltools_reload_keeps_real_native_admission_without_new_audits():
    import os,subprocess
    code='''
import importlib,sys
import symox as s
from symox.integrals import risch
selector=importlib.import_module('symox.integrals._risch_splitfactor_native')
rt=importlib.import_module('symox.integrals.rationaltools')
importlib.reload(rt)
t=s.Symbol('split_reload_t');p=s.Poly((t*t+1)**2*(t+2),t,domain='QQ');one=s.Poly(1,t,domain='QQ')
de=risch.DifferentialExtension(extension={'D':[s.Poly(t*t+1,t)]});events=[]
def reject(event,args):
 if event in ('object.__getattr__','marshal.dumps'):
  events.append(event);raise RuntimeError('unexpected metadata audit')
sys.addaudithook(reject)
assert selector.reduce(p,de,one,False,None)is not None
actual=risch.splitfactor(p,de);risch._native_splitfactor=lambda *args:None
expected=risch.splitfactor(p,de)
assert actual==expected and events==[]
'''
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())


@pytest.mark.parametrize('mode',['code','nested_constants','foreign_constant','globals','f_builtins',
    'defaults','kwdefaults','closure','attributes'])
def test_reloaded_helpers_are_authenticated_against_bootstrap_without_callbacks(mode):
    import os,subprocess
    code='''
import importlib,sys,types
import symox as s
from symox.integrals import risch
selector=importlib.import_module('symox.integrals._risch_splitfactor_native')
rt=importlib.import_module('symox.integrals.rationaltools');importlib.reload(rt)
t=s.Symbol('split_reload_t');p=s.Poly((t*t+1)**2*(t+2),t,domain='QQ');one=s.Poly(1,t,domain='QQ')
de=risch.DifferentialExtension(extension={'D':[s.Poly(t*t+1,t)]});events=[]
function=(rt._atan_stock_operators_with_metadata if MODE=='nested_constants'else rt._atan_vector_proof_with_metadata);namespace=function.__globals__
class Foreign:
 def __eq__(self,other):events.append('eq');raise RuntimeError('foreign constant equality')
 def __hash__(self):events.append('hash');raise RuntimeError('foreign constant hash')
class Tuple(tuple):
 def __iter__(self):events.append('iter');raise RuntimeError('closure iteration')
 def __len__(self):events.append('len');raise RuntimeError('closure length')
class Dict(dict):
 def __iter__(self):events.append('iter');raise RuntimeError('dictionary iteration')
 def __len__(self):events.append('len');raise RuntimeError('dictionary length')
mode=MODE
if mode=='code':function.__code__=(lambda *args:None).__code__
elif mode in ('nested_constants','foreign_constant'):
 constants=list(function.__code__.co_consts)
 if mode=='nested_constants':
  i=next(i for i,v in enumerate(constants)if type(v)is types.CodeType)
  child=constants[i];constants[i]=child.replace(co_consts=child.co_consts+('changed nested constant',))
 else:constants.append(Foreign())
 function.__code__=function.__code__.replace(co_consts=tuple(constants))
elif mode=='globals':
 rt._atan_vector_proof_with_metadata=types.FunctionType(function.__code__,dict(namespace))
elif mode=='f_builtins':
 saved=namespace['__builtins__'];namespace['__builtins__']=Dict(function.__builtins__)
 try:clone=types.FunctionType(function.__code__,namespace)
 finally:namespace['__builtins__']=saved
 rt._atan_vector_proof_with_metadata=clone
elif mode=='defaults':function.__defaults__=(None,)
elif mode=='kwdefaults':function.__kwdefaults__=Dict()
elif mode=='closure':
 rt._atan_vector_proof_with_metadata=types.FunctionType(function.__code__,namespace,closure=Tuple(()))
else:function.__dict__=Dict()
def reject(event,args):
 if event in ('object.__getattr__','marshal.dumps'):
  events.append(event);raise RuntimeError('unexpected metadata audit')
sys.addaudithook(reject)
assert selector.reduce(p,de,one,False,None)is None
actual=risch.splitfactor(p,de);risch._native_splitfactor=lambda *args:None
expected=risch.splitfactor(p,de)
assert actual==expected and events==[]
'''.replace('MODE',repr(mode))
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())


@pytest.mark.parametrize('raises',[False,True])
def test_fresh_bootstrap_proves_actual_math_function_builtins_after_global_restore(raises):
    import os,subprocess
    code='''
import builtins,sys,types,importlib.abc,importlib.machinery
original_len=builtins.len;events=[];installed=False
def length(value):
 events.append('len')
 if RAISES:raise RuntimeError('actual derivation f_builtins')
 return original_len(value)
class Finder(importlib.abc.MetaPathFinder):
 def find_spec(self,name,path=None,target=None):
  if name!='symox.integrals.risch':return None
  spec=importlib.machinery.PathFinder.find_spec(name,path);delegate=spec.loader
  class Loader:
   def create_module(self,spec):return delegate.create_module(spec)
   def exec_module(self,r):
    global installed
    delegate.exec_module(r);function=r.derivation;namespace=function.__globals__
    saved=namespace['__builtins__'];foreign=dict(function.__builtins__);foreign['len']=length
    namespace['__builtins__']=foreign
    try:clone=types.FunctionType(function.__code__,namespace,argdefs=function.__defaults__,closure=function.__closure__)
    finally:namespace['__builtins__']=saved
    assert clone.__globals__ is namespace and clone.__builtins__ is foreign
    r.derivation=clone;installed=True
  spec.loader=Loader();return spec
sys.meta_path.insert(0,Finder())
import symox as s
from symox.integrals import risch
selector=sys.modules['symox.integrals._risch_splitfactor_native']
t=s.Symbol('split_fb_t');p=s.Poly((t*t+1)**2*(t+2),t,domain='QQ')
de=risch.DifferentialExtension(extension={'D':[s.Poly(t*t+1,t)]})
assert installed and selector._ROOT is None and selector._shared._ROOT is not None
def outcome():
 try:return ('value',risch.splitfactor(p,de))
 except Exception as error:return ('error',type(error).__name__,str(error))
events.clear();actual=outcome();observed=events[:];events.clear()
risch._native_splitfactor=lambda *args:None
expected=outcome()
assert actual==expected and events==observed and events
if RAISES:assert actual==('error','RuntimeError','actual derivation f_builtins')
'''.replace('RAISES',repr(raises))
    subprocess.run([sys.executable,'-c',code],check=True,env=os.environ.copy())
