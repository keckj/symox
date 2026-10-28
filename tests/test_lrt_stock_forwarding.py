"""LRT uses authenticated scalar forwarding for ordinary Symbol/Dummy inputs."""
import importlib

import pytest
import sympy as sp
import symox as sx
import symox.holonomic.recurrence

native = importlib.import_module('symox.integrals.rationaltools')
reference = importlib.import_module('sympy.integrals.rationaltools')


@pytest.mark.parametrize('domain', ['ZZ', 'QQ'])
@pytest.mark.parametrize('dummy_x', [False, True])
@pytest.mark.parametrize('dummy_t', [False, True])
def test_official_forwarding_keeps_lrt_native(domain, dummy_x, dummy_t, monkeypatch):
    x = sx.Dummy('x') if dummy_x else sx.Symbol('x')
    t = sx.Dummy('t') if dummy_t else sx.Symbol('t')
    f = sx.Poly(1, x, domain=domain)
    g = sx.Poly(x*x+x+1, x, domain=domain)
    calls = []
    original = sx._core._ratint_logpart

    def selected(*args):
        result = original(*args)
        calls.append(result is not None)
        return result

    monkeypatch.setattr(sx._core, '_ratint_logpart', selected)
    result = native.ratint_logpart(f, g, x, t)
    X, T = x.to_sympy(), t.to_sympy()
    expected = reference.ratint_logpart(sp.Poly(1, X, domain=domain),
                                      sp.Poly(X*X+X+1, X, domain=domain), X, T)
    assert calls == [True]
    assert [(a.as_expr().to_sympy(), str(a.domain),
             b.as_expr().to_sympy(), str(b.domain)) for a, b in result] == [
                 (a.as_expr(), str(a.domain), b.as_expr(), str(b.domain))
                 for a, b in expected]


def test_wild_parameter_keeps_complete_fallback(monkeypatch):
    x = sx.Symbol('x')
    t = sx.Wild('t')
    f, g = sx.Poly(1, x), sx.Poly(x*x+x+1, x)

    def run():
        try:
            result = native.ratint_logpart(f, g, x, t)
        except Exception as error:
            return type(error), str(error)
        return [(a.as_expr(), b.as_expr(), a.domain, b.domain) for a, b in result]

    with monkeypatch.context() as patch:
        patch.setattr(sx._core, '_ratint_logpart', lambda *args: None)
        expected = run()

    def forbidden(*args):
        raise AssertionError('Wild parameter reached stock scalar forwarding')

    monkeypatch.setattr(sx._core, '_ratint_logpart', forbidden)
    assert run() == expected


def _lrt_pair():
    x,t=sx.symbols('x t')
    return x,t,sx.Poly(1,x),sx.Poly(x*x+x+1,x)


def test_reload_keeps_the_verified_official_resultant_wrapper(monkeypatch):
    import symox.polys as polys
    x,t,f,g=_lrt_pair()
    # The public generator-ordering wrapper is installed after rationaltools'
    # first bootstrap import. Reload must recognize that exact forwarder.
    for _ in range(3):
        importlib.reload(native)
        assert native.resultant is polys.resultant
        assert native._stock_lrt_resultant()
        calls=[];old=sx._core._ratint_logpart
        def selected(*args):
            value=old(*args);calls.append(value is not None);return value
        with monkeypatch.context() as patch:
            patch.setattr(sx._core,'_ratint_logpart',selected)
            native.ratint_logpart(f,g,x,t)
        assert calls==[True]


@pytest.mark.parametrize('owner',['base','public','symbol','one'])
def test_priority_descriptors_decline_without_calling_them_and_keep_old_errors(owner,monkeypatch):
    x,t,f,g=_lrt_pair();calls=[]
    cls={'base':sx._core.Expr,'public':sx.Expr,'symbol':type(x),'one':type(sx.S.One)}[owner]
    def observed(self):calls.append(type(self).__name__);raise KeyboardInterrupt('LRT priority callback')
    monkeypatch.setattr(cls,'_op_priority',property(observed))
    assert not native._atan_stock_operators(f,g,native._LRT_OPERATORS,native._LRT_BOOTSTRAP_STATES)
    assert calls==[]
    def forbidden(*args):raise AssertionError('custom priority reached native LRT')
    monkeypatch.setattr(sx._core,'_ratint_logpart',forbidden)
    with pytest.raises(KeyboardInterrupt,match='LRT priority callback'):native.ratint_logpart(f,g,x,t)
    trace=list(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    with pytest.raises(KeyboardInterrupt,match='LRT priority callback'):native.ratint_logpart(f,g,x,t)
    assert calls==trace and calls


def test_lrt_one_arithmetic_override_is_observed(monkeypatch):
    x,t,f,g=_lrt_pair();calls=[]
    def observed(self,other):calls.append(type(other).__name__);raise KeyboardInterrupt('LRT scalar multiply')
    monkeypatch.setattr(type(sx.S.One),'__mul__',observed)
    def forbidden(*args):raise AssertionError('custom scalar operator reached native LRT')
    monkeypatch.setattr(sx._core,'_ratint_logpart',forbidden)
    with pytest.raises(KeyboardInterrupt,match='LRT scalar multiply'):native.ratint_logpart(f,g,x,t)
    trace=list(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    with pytest.raises(KeyboardInterrupt,match='LRT scalar multiply'):native.ratint_logpart(f,g,x,t)
    assert calls==trace and calls


def test_user_wraps_resultant_does_not_acquire_official_forwarding_proof(monkeypatch):
    from functools import wraps
    x,t,f,g=_lrt_pair();calls=[];old=native.resultant
    @wraps(old)
    def observed(*args,**kwargs):calls.append(1);return old(*args,**kwargs)
    monkeypatch.setattr(native,'resultant',observed)
    assert not native._stock_lrt_resultant()
    native.ratint_logpart(f,g,x,t)
    assert calls==[1]


def test_changed_official_resultant_delegate_keeps_callback(monkeypatch):
    x,t,f,g=_lrt_pair()
    importlib.reload(native)
    wrapper=native.resultant
    cells=dict(zip(wrapper.__code__.co_freevars,wrapper.__closure__))
    calls=[]
    def observed(*args,**kwargs):calls.append(1);raise KeyboardInterrupt('resultant delegate')
    monkeypatch.setattr(cells['native'],'cell_contents',observed)
    assert not native._stock_lrt_resultant()
    with pytest.raises(KeyboardInterrupt,match='resultant delegate'):native.ratint_logpart(f,g,x,t)
    assert calls==[1]


def test_untrusted_fresh_module_capture_never_enables_native():
    import os,subprocess,sys
    script="""import importlib,sys
import symox as S
import symox.holonomic.recurrence
r=importlib.import_module('symox.integrals.rationaltools')
x,t=S.symbols('x t');f,g=S.Poly(1,x),S.Poly(x*x+x+1,x)
expected=[(h.as_expr(),q.as_expr()) for h,q in r.ratint_logpart(f,g,x,t)]
old=S._core.Expr.__mul__;calls=[]
def observed(self,other):calls.append(1);return old(self,other)
S._core.Expr.__mul__=observed
sys.modules.pop('symox.integrals.rationaltools')
r=importlib.import_module('symox.integrals.rationaltools')
assert not r._ATAN_BOOTSTRAP_CAPTURE
def forbidden(*args):raise AssertionError('untrusted capture reached native')
S._core._ratint_logpart=forbidden
assert [(h.as_expr(),q.as_expr()) for h,q in r.ratint_logpart(f,g,x,t)]==expected
assert calls
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


def test_restored_priority_and_resultant_contexts_recover_native_selection(monkeypatch):
    x,t,f,g=_lrt_pair();old=sx._core._ratint_logpart
    with monkeypatch.context() as patch:
        patch.setattr(type(x),'_op_priority',property(lambda self:12))
        assert not native._atan_stock_operators(f,g,native._LRT_OPERATORS,native._LRT_BOOTSTRAP_STATES)
    importlib.reload(native)
    calls=[]
    def selected(*args):
        value=old(*args);calls.append(value is not None);return value
    monkeypatch.setattr(sx._core,'_ratint_logpart',selected)
    native.ratint_logpart(f,g,x,t)
    assert calls==[True]


def test_atan_quartic_priority_callback_uses_the_same_common_guard(monkeypatch):
    x=sx.Symbol('x');f=sx.Poly(x**4+x+1,x,domain='QQ');g=sx.Poly(x**3+2,x,domain='QQ');calls=[]
    def observed(self):calls.append(type(self).__name__);raise KeyboardInterrupt('Atan priority callback')
    monkeypatch.setattr(sx._core.Expr,'_op_priority',property(observed))
    assert not native._atan_stock_operators(f,g) and calls==[]
    def forbidden(*args):raise AssertionError('priority callback reached native Atan')
    monkeypatch.setattr(sx._core,'_log_to_atan',forbidden)
    with pytest.raises(KeyboardInterrupt,match='Atan priority callback'):native.log_to_atan(f,g)
    trace=list(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_log_to_atan',lambda *args:None)
    with pytest.raises(KeyboardInterrupt,match='Atan priority callback'):native.log_to_atan(f,g)
    assert trace==calls and calls


def test_priority_guard_applies_without_importing_ore():
    import os,subprocess,sys
    script="""import importlib,sys
import symox as S
r=importlib.import_module('symox.integrals.rationaltools')
assert 'symox.holonomic.recurrence' not in sys.modules
x=S.Symbol('x');f=S.Poly(x**4+x+1,x,domain='QQ');g=S.Poly(x**3+2,x,domain='QQ');calls=[]
def observed(self):calls.append(type(self).__name__);raise KeyboardInterrupt('plain priority callback')
S._core.Expr._op_priority=property(observed)
assert not r._atan_stock_operators(f,g) and not calls
def forbidden(*args):raise AssertionError('priority reached native without Ore')
S._core._log_to_atan=forbidden
try:r.log_to_atan(f,g)
except KeyboardInterrupt as error:assert str(error)=='plain priority callback'
else:raise AssertionError('priority callback was hidden')
assert calls
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


def test_lrt_generator_power_override_is_observed(monkeypatch):
    x,t,f,g=_lrt_pair();calls=[]
    def observed(self,other):calls.append(other);raise KeyboardInterrupt('LRT scalar power')
    monkeypatch.setattr(type(x),'__pow__',observed)
    def forbidden(*args):raise AssertionError('custom power reached native LRT')
    monkeypatch.setattr(sx._core,'_ratint_logpart',forbidden)
    with pytest.raises(KeyboardInterrupt,match='LRT scalar power'):native.ratint_logpart(f,g,x,t)
    trace=list(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    with pytest.raises(KeyboardInterrupt,match='LRT scalar power'):native.ratint_logpart(f,g,x,t)
    assert trace==calls and calls


@pytest.mark.parametrize('imports', [
    ('symox.vector',),
    ('symox.vector', 'symox.holonomic.recurrence'),
    ('symox.holonomic.recurrence', 'symox.vector'),
    ('symox.tensor.tensor', 'symox.vector', 'symox.holonomic.recurrence'),
])
def test_vector_and_ore_import_orders_keep_both_native_routes(imports):
    import os, subprocess, sys
    script = """import importlib
import symox as S
r=importlib.import_module('symox.integrals.rationaltools')
for name in IMPORTS: importlib.import_module(name)
x,t=S.symbols('x t');f,g=S.Poly(1,x),S.Poly(x*x+x+1,x)
for reload in range(3):
 importlib.reload(r)
 assert r._atan_stock_operators(f,g,r._LRT_OPERATORS,r._LRT_BOOTSTRAP_STATES)
 calls=[];old=S._core._ratint_logpart
 def selected(*args):
  value=old(*args);calls.append(value is not None);return value
 S._core._ratint_logpart=selected
 try:result=r.ratint_logpart(f,g,x,t)
 finally:S._core._ratint_logpart=old
 assert calls==[True] and len(result)==1
 f,g=S.Poly(x**4+x+1,x,domain='QQ'),S.Poly(x**3+2,x,domain='QQ')
 assert r._stock_atan_protocol(f,g)
 calls=[];old=S._core._log_to_atan
 def selected(*args):
  value=old(*args);calls.append(value is not None);return value
 S._core._log_to_atan=selected
 try:result=r.log_to_atan(f,g)
 finally:S._core._log_to_atan=old
 assert calls==[True]
 f,g=S.Poly(1,x),S.Poly(x*x+x+1,x)
from symox.vector import CoordSys3D
N=CoordSys3D('N')
assert x*N.i==N.i*x
try:x/N.i
except TypeError as error:assert str(error)=='division by a vector or dyadic is undefined'
else:raise AssertionError('vector division lost its error')
"""
    script = 'IMPORTS=' + repr(imports) + '\n' + script
    subprocess.run([sys.executable, '-c', script], env=os.environ.copy(), check=True)


@pytest.mark.parametrize('mode', ['delegate', 'isinstance', 'wraps'])
def test_vector_forwarding_changes_keep_callbacks_and_complete_fallback(mode, monkeypatch):
    from functools import wraps
    import symox.vector._adapters as adapters
    x,t,f,g=_lrt_pair();calls=[]
    proof=next(p for p in adapters._SCALAR_DISPATCH_PROOFS if p[-1]=='__mul__')
    operation,delegate=proof[:2]
    def observed(*args,**kwargs):
        calls.append(1);raise KeyboardInterrupt('vector forwarding callback')
    if mode=='delegate':
        cells=dict(zip(operation.__code__.co_freevars,operation.__closure__))
        monkeypatch.setattr(cells['original'],'cell_contents',observed)
    elif mode=='isinstance':
        monkeypatch.setattr(adapters,'isinstance',observed,raising=False)
    else:
        monkeypatch.setattr(sx._core.Expr,'__mul__',wraps(sx._core.Expr.__mul__)(observed))
    assert not native._atan_stock_operators(f,g,native._LRT_OPERATORS,native._LRT_BOOTSTRAP_STATES)
    assert calls==[]
    def forbidden(*args):raise AssertionError('changed vector forwarding reached native LRT')
    monkeypatch.setattr(sx._core,'_ratint_logpart',forbidden)
    with pytest.raises(KeyboardInterrupt,match='vector forwarding callback'):
        native.ratint_logpart(f,g,x,t)
    trace=list(calls);calls.clear()
    monkeypatch.setattr(sx._core,'_ratint_logpart',lambda *args:None)
    with pytest.raises(KeyboardInterrupt,match='vector forwarding callback'):
        native.ratint_logpart(f,g,x,t)
    assert trace==calls and calls


def test_vector_forwarding_cycle_declines_without_recursing(monkeypatch):
    import symox.vector._adapters as adapters
    x,t,f,g=_lrt_pair()
    proof=adapters._SCALAR_DISPATCH_PROOFS[0]
    operation=proof[0]
    cells=dict(zip(operation.__code__.co_freevars,operation.__closure__))
    monkeypatch.setattr(cells['original'],'cell_contents',operation)
    assert not native._atan_stock_operators(f,g,native._LRT_OPERATORS,native._LRT_BOOTSTRAP_STATES)


def test_marker_checks_follow_the_actual_ore_layers_without_losing_native_vector_routes():
    import os,subprocess,sys
    script="""import sys,importlib
sys.modules['sympy']=None
import symox as S
import symox.vector
r=importlib.import_module('symox.integrals.rationaltools')
assert 'symox.holonomic.recurrence' not in sys.modules
x=S.Symbol('x');f,g=S.Poly(x*x+1,x,domain='QQ'),S.Poly(x-1,x,domain='QQ')
calls=[]
def marker(self):calls.append(1);raise KeyboardInterrupt('actual Ore marker')
S._core.Expr._symox_holonomic_function=property(marker)
assert r._stock_atan_protocol(f,g) and not calls
expected=r.log_to_atan(f,g)
old=S._core._log_to_atan
S._core._log_to_atan=lambda *args:None
try:assert r.log_to_atan(f,g)==expected and not calls
finally:S._core._log_to_atan=old
# The same property becomes observable once an actual Ore layer is installed.
import symox.holonomic.recurrence as ore
assert not r._stock_atan_protocol(f,g) and not calls
try:r.log_to_atan(f,g)
except KeyboardInterrupt as error:assert str(error)=='actual Ore marker'
else:raise AssertionError('real Ore marker callback was bypassed')
assert calls
calls.clear()
# A loaded module alone does not imply that its wrappers are still installed.
for proof in ore._SCALAR_DISPATCH_PROOFS:
 operation,delegate=proof[:2]
 name=dict(zip(operation.__code__.co_freevars,operation.__closure__))['name'].cell_contents
 setattr(S._core.Expr,name,delegate)
assert r._stock_atan_protocol(f,g) and not calls
selected=[]
def spy(*args):
 value=old(*args);selected.append(value is not None);return value
S._core._log_to_atan=spy
assert r.log_to_atan(f,g)==expected and selected==[True] and not calls
"""
    subprocess.run([sys.executable,'-c',script],env=os.environ.copy(),check=True)


def test_priority_sentinel_is_a_stored_override_not_a_missing_definition(monkeypatch):
    x,t,f,g=_lrt_pair()
    monkeypatch.setattr(type(x),'_op_priority',native._ATAN_MISSING)
    assert not native._atan_stock_operators(f,g,native._LRT_OPERATORS,native._LRT_BOOTSTRAP_STATES)
