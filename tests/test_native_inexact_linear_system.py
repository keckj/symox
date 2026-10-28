"""Bit-exact differential against the established inexact domain elimination."""
import importlib
import random
from types import SimpleNamespace

import pytest
import symox as sx
import sympy as sp

solvers = importlib.import_module('symox.polys.solvers')
P = (1, 53, 139, 140, 211)
x, y, t = sx.symbols('x y t')


def signature(expr):
    if hasattr(expr, 'as_expr'):
        expr = expr.as_expr()
    return (str(expr), sorted((f._prec, f._mpf_) for f in expr.atoms(sx.Float)))


def row_signature(rows):
    return [[(str(a.parent()), signature(a.parent().to_sympy(a))) for a in row] for row in rows]


def old_reduce(rows, k, width):
    """Unchanged reduction loop from the pre-port public implementation."""
    rows = [row[:] for row in rows]
    pivots = []; r = 0
    for col in range(width):
        piv = next((i for i in range(r, len(rows)) if rows[i][col]), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = k.one / rows[r][col]
        rows[r] = [v*inv for v in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col]:
                f = rows[i][col]
                rows[i] = [a-f*b for a,b in zip(rows[i], rows[r])]
        pivots.append(col); r += 1
        if r == len(rows):
            break
    return rows, pivots


def domain(p, complex_, composite=''):
    k = sx.ComplexField(p) if complex_ else sx.RealField(p)
    if composite == 'ring':
        return k.poly_ring(t)
    if composite == 'field':
        return k.frac_field(t)
    return k


@pytest.mark.parametrize('p', P)
@pytest.mark.parametrize('complex_', (False, True))
@pytest.mark.parametrize('seed', range(8))
def test_full_rref_stored_binary_values_and_pivot_order(p, complex_, seed):
    k = domain(p, complex_); rng = random.Random(10237+seed)
    raw = [[sx.Rational(rng.randrange(-23,24),rng.randrange(1,13)) +
            (sx.I*sx.Rational(rng.randrange(-7,8),11) if complex_ else 0)
            for j in range(4)] for i in range(3)]
    raw[0][0] = sx.Rational(1,2**71)  # first pivot, never the largest pivot
    if seed % 3 == 0:
        raw[-1] = raw[1][:]
    rows = [[k.convert(a) for a in row] for row in raw]
    before = row_signature(rows)
    expected, pivots = old_reduce(rows,k,4)
    actual = solvers._native_inexact_rref(rows,k,4)
    assert actual is not None
    assert actual[1] == pivots
    assert row_signature(actual[0]) == row_signature(expected)
    assert row_signature(rows) == before
    assert all(a.parent() == k for row in actual[0] for a in row)


@pytest.mark.parametrize('p', P)
@pytest.mark.parametrize('complex_', (False, True))
@pytest.mark.parametrize('composite', ('ring','field'))
def test_composite_domain_and_parameter_cancellation(p,complex_,composite):
    k = domain(p,complex_,composite)
    a = t + (sx.I if complex_ else 1)
    # Unit pivots over a polynomial ring; nonconstant pivots over its field.
    raw = [[1,a,2*a+1],[0,1,2]] if composite=='ring' else [[a,1,a+2],[1,0,1]]
    rows = [[k.convert(v) for v in row] for row in raw]
    expected,pivots = old_reduce(rows,k,3)
    actual = solvers._native_inexact_rref(rows,k,3)
    assert actual is not None and actual[1] == pivots
    assert row_signature(actual[0]) == row_signature(expected)
    assert all(a.parent() == k for row in actual[0] for a in row)


@pytest.mark.parametrize('complex_', (False,True))
@pytest.mark.parametrize('p', P)
@pytest.mark.parametrize('kind', ('unique','free','inconsistent','zero','tiny'))
def test_public_solver_is_bit_identical_to_previous_complete_path(monkeypatch,p,complex_,kind):
    k = domain(p,complex_)
    z = sx.I if complex_ else 0
    eqs = {'unique':[3*x+(7+z)*y-11,5*x+2*y-7],
           'free':[x+(3+z)*y-2], 'inconsistent':[x-1,x-2],
           'zero':[sx.S.Zero], 'tiny':[sx.Rational(1,10**60)*x-y, x-3]}[kind]
    ring = SimpleNamespace(symbols=(x,y),domain=k)
    actual = solvers.solve_lin_sys(eqs,ring,False)
    monkeypatch.setattr(solvers,'_native_inexact_rref',lambda *args:None)
    expected = solvers.solve_lin_sys(eqs,ring,False)
    assert (actual is None) == (expected is None)
    if actual is not None:
        assert actual.keys() == expected.keys()
        assert {str(a):signature(v) for a,v in actual.items()} == {str(a):signature(v) for a,v in expected.items()}


@pytest.mark.parametrize('k', (sx.RR,sx.CC,sx.RealField(140),sx.ComplexField(140)))
def test_raw_ring_result_and_exact_reconstruction(monkeypatch,k):
    from symox.polys.rings import ring
    r,a,b = ring('a,b',k)
    eqs=[a+2*b-3,3*a+7*b-11]
    actual=solvers.solve_lin_sys(eqs,r)
    monkeypatch.setattr(solvers,'_native_inexact_rref',lambda *args:None)
    expected=solvers.solve_lin_sys(eqs,r)
    assert actual.keys()==expected.keys()
    assert {str(a):signature(v) for a,v in actual.items()} == {str(a):signature(v) for a,v in expected.items()}
    assert actual[a]==k(-1) and actual[b]==k(2)
    assert actual[a]+k(2)*actual[b] == k(3)
    assert k(3)*actual[a]+k(7)*actual[b] == k(11)
    assert all(v.parent()==k for v in actual.values())


@pytest.mark.parametrize('attribute', ('__mul__','__sub__','__truediv__','__bool__'))
def test_operator_callbacks_keep_full_python_loop(monkeypatch,attribute):
    cls=sx._core.DomainElement; original=getattr(cls,attribute); calls=[]
    def custom(self,*args):
        calls.append(attribute)
        return original(self,*args)
    monkeypatch.setattr(cls,attribute,custom)
    rows=[[sx.RR(1),sx.RR(2)],[sx.RR(3),sx.RR(5)]]
    assert solvers._native_inexact_rref(rows,sx.RR,2) is None
    assert solvers.solve_lin_sys([x+y-2,3*x+2*y-5],SimpleNamespace(symbols=(x,y),domain=sx.RR),False)
    assert calls


@pytest.mark.parametrize('attribute', ('one','convert','to_sympy','is_QQ_I'))
def test_domain_callbacks_preserve_order_and_old_results(monkeypatch,attribute):
    cls=sx._core.Domain; original=cls.__dict__[attribute]; calls=[]
    if isinstance(original,property) or hasattr(original,'__get__') and attribute in ('one','is_QQ_I'):
        def getter(self):
            calls.append(attribute); return original.__get__(self,cls)
        replacement=property(getter)
    else:
        def replacement(self,*args,**kwargs):
            calls.append(attribute); return original(self,*args,**kwargs)
    monkeypatch.setattr(cls,attribute,replacement)
    native=solvers._native_inexact_rref
    r=SimpleNamespace(symbols=(x,y),domain=sx.RR)
    result=solvers.solve_lin_sys([x+y-2,3*x+2*y-5],r,False)
    actual_calls=calls[:]; calls.clear()
    monkeypatch.setattr(solvers,'_native_inexact_rref',lambda *args:None)
    expected=solvers.solve_lin_sys([x+y-2,3*x+2*y-5],r,False)
    assert result==expected and actual_calls==calls and calls


@pytest.mark.parametrize('name', ('next','range','zip'))
def test_module_builtin_callbacks_are_not_skipped(monkeypatch,name):
    import builtins
    calls=[];original=getattr(builtins,name)
    def custom(*args):
        calls.append(name); return original(*args)
    monkeypatch.setattr(solvers,name,custom,raising=False)
    r=SimpleNamespace(symbols=(x,y),domain=sx.RR)
    result=solvers.solve_lin_sys([x+y-2,3*x+2*y-5],r,False)
    actual_calls=calls[:];calls.clear()
    monkeypatch.setattr(solvers,'_native_inexact_rref',lambda *args:None)
    expected=solvers.solve_lin_sys([x+y-2,3*x+2*y-5],r,False)
    assert result==expected and actual_calls==calls and calls


def test_attribute_probe_callback_is_preserved(monkeypatch):
    calls=[]
    def getter(self):
        calls.append(1); return False
    monkeypatch.setattr(sx._core.DomainElement,'is_Matrix',property(getter),raising=False)
    rows=[[sx.RR(1),sx.RR(2)]]
    assert solvers._native_inexact_rref(rows,sx.RR,2) is None
    actual=solvers.solve_lin_sys([x-2],SimpleNamespace(symbols=(x,),domain=sx.RR),False)
    assert actual[x]==sx.Float(2) and calls


def test_direct_binding_declines_other_domains_and_invalid_storage():
    native=sx._core._inexact_linear_rref
    assert native([[sx.QQ(1)]],sx.QQ,1) is None
    assert native(((sx.RR(1),),),sx.RR,1) is None
    assert native([[sx.RR(1),sx.RR(2)]],sx.RR,1) is None
    assert native([[sx.CC(1)]],sx.RR,1) is None
    assert native([[1]],sx.RR,1) is None
    assert native([],sx.RR,2)==([],[])
    assert native([[sx.RR(1)]],sx.RR,1)[1]==[0]  # augmented pivot, not "unhandled"


@pytest.mark.parametrize('complex_', (False,True))
def test_polynomial_nonunit_error_matches_existing_path(complex_):
    k=domain(140,complex_,'ring');rows=[[k.convert(t),k.one]]
    def error(function):
        try:function()
        except Exception as exc:return type(exc),str(exc)
        pytest.fail('non-unit unexpectedly invertible in polynomial ring')
    assert error(lambda:solvers._native_inexact_rref(rows,k,2))==error(lambda:old_reduce(rows,k,2))


@pytest.mark.parametrize('p',P)
@pytest.mark.parametrize('complex_',(False,True))
def test_reference_exact_binary_coefficients_with_unique_dyadic_solution(p,complex_):
    from sympy.polys.rings import ring
    from sympy.polys.solvers import solve_lin_sys
    k=domain(p,complex_); rk=sp.ComplexField(p) if complex_ else sp.RealField(p)
    equations=[x+y-2,y-1]
    result=solvers.solve_lin_sys(equations,SimpleNamespace(symbols=(x,y),domain=k),False)
    # Both systems have an exact dyadic solution, independently certified.
    r,a,b=ring('a,b',rk)
    expected=solve_lin_sys([a+b-2,b-1],r)
    assert sp.Rational(str(sx.Rational(sx.re(result[x]))))==sp.Rational(sp.re(rk.to_sympy(expected[a])))==1
    assert sp.Rational(str(sx.Rational(sx.re(result[y]))))==sp.Rational(sp.re(rk.to_sympy(expected[b])))==1
    assert sx.im(result[x])==0 and sx.im(result[y])==0


@pytest.mark.parametrize('value',(sx.oo,-sx.oo,sx.nan,sx.zoo))
@pytest.mark.parametrize('k',(sx.RR,sx.CC))
def test_nonfinite_conversion_error_is_unchanged(monkeypatch,k,value):
    ring=SimpleNamespace(symbols=(x,),domain=k)
    matrix=sx.Matrix([[value]]);rhs=sx.Matrix([1])
    def outcome():
        try:
            return ('result',solvers._solve_inexact(matrix,rhs,[x],ring,False))
        except Exception as exc:
            return (type(exc),str(exc))
    actual=outcome()
    monkeypatch.setattr(solvers,'_native_inexact_rref',lambda *args:None)
    assert actual==outcome() and actual[0]!='result'


def test_custom_raw_ring_factories_and_metadata_keep_protocol(monkeypatch):
    calls=[]
    class Ring:
        domain=sx.RR
        symbols=(x,y)
        @property
        def gens(self):
            calls.append('gens');return ('X','Y')
        def __call__(self,value):
            calls.append(('ring',str(value)));return ('wrapped',str(value))
    ring=Ring()
    actual=solvers.solve_lin_sys([x+3*y-2],ring)
    before=calls[:];calls.clear()
    monkeypatch.setattr(solvers,'_native_inexact_rref',lambda *args:None)
    assert actual==solvers.solve_lin_sys([x+3*y-2],ring) and calls==before
    assert actual['X'][0]=='wrapped' and len(before)==2


def test_augmented_column_inconsistency_and_zeros():
    for rows,width in [([],3),([[sx.RR.zero,sx.RR.zero]],2),([[sx.RR.zero,sx.RR.one]],2)]:
        expected,pivots=old_reduce(rows,sx.RR,width)
        actual=solvers._native_inexact_rref(rows,sx.RR,width)
        assert actual is not None and actual[1]==pivots
        assert row_signature(actual[0])==row_signature(expected)


@pytest.mark.parametrize('p,complex_',[(p,False) for p in P]+[(p,True) for p in (1,53,139)])
def test_reference_full_rows_from_identical_stored_binary_inputs(p,complex_):
    k=domain(p,complex_)
    rk=sp.ComplexField(p) if complex_ else sp.RealField(p)
    raw=[[(2,1),(3,-2),(5,3)],[(7,-1),(11,2),(17,1)]]
    rows=[[k.convert(sx.Rational(a,7)+(sx.I*sx.Rational(b,11) if complex_ else 0))
           for a,b in row] for row in raw]
    def rational(value):
        return sp.Rational(str(sx.Rational(value)))
    rr=[[rk.convert(rational(sx.re(k.to_sympy(v)))+sp.I*rational(sx.im(k.to_sympy(v))))
         for v in row] for row in rows]
    actual,pivots=solvers._native_inexact_rref(rows,k,3)
    expected,ref_pivots=old_reduce(rr,rk,3)
    assert pivots==ref_pivots
    for ar,br in zip(actual,expected):
        for a,b in zip(ar,br):
            assert rational(sx.re(k.to_sympy(a)))==sp.Rational(sp.re(rk.to_sympy(b)))
            assert rational(sx.im(k.to_sympy(a)))==sp.Rational(sp.im(rk.to_sympy(b)))


@pytest.mark.parametrize('p',(53,140))
@pytest.mark.parametrize('complex_',(False,True))
def test_heurisch_calls_native_without_changing_primitive(monkeypatch,p,complex_):
    from symox.integrals.heurisch import heurisch
    a=sx.Float('0.5',precision=p)
    coefficient=a*sx.I if complex_ else a
    integrand=sx.exp(coefficient*x)
    native=sx._core._inexact_linear_rref;calls=[]
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_inexact_linear_rref',spy)
    actual=heurisch(integrand,x)
    assert calls and all(calls)
    monkeypatch.setattr(solvers,'_native_inexact_rref',lambda *args:None)
    expected=heurisch(integrand,x)
    assert signature(actual)==signature(expected)
    assert sp.simplify(sp.diff(sx.to_sympy(actual),sp.Symbol('x'))-sx.to_sympy(integrand))==0


@pytest.mark.parametrize('p',(53,140))
def test_heurisch_existing_unknown_is_preserved(monkeypatch,p):
    from symox.integrals.heurisch import heurisch
    integrand=sx.exp(sx.Float('0.5',precision=p)*(1+sx.I)*x)
    actual=heurisch(integrand,x)
    monkeypatch.setattr(solvers,'_native_inexact_rref',lambda *args:None)
    assert actual is None and heurisch(integrand,x) is None
