"""Exact ground PolyMatrix RREF/nullspace and complete protocol fallbacks."""
import importlib
import random

import pytest
import symox as sx
import sympy as sp
from sympy.polys.polymatrix import PolyMatrix as RefPolyMatrix

pm=importlib.import_module('symox.polys.polymatrix')
PolyMatrix=pm.PolyMatrix
x,y,z=sx.symbols('x y z')
rx,ry,rz=sp.symbols('x y z')


def fields(kind):
    if kind=='QQ':return sx.QQ,sp.QQ
    if kind=='ZRF':return sx.ZZ.frac_field(y),sp.ZZ.frac_field(ry)
    if kind=='RF':return sx.QQ.frac_field(y),sp.QQ.frac_field(ry)
    if kind=='multi':return sx.QQ.frac_field(y,z),sp.QQ.frac_field(ry,rz)
    if kind=='nested':return sx.QQ.frac_field(y).frac_field(z),sp.QQ.frac_field(ry).frac_field(rz)
    raise AssertionError(kind)


def exact_entries(matrix):
    return [sx.to_sympy(e) for e in matrix.flat()]


def same(actual,expected):
    assert actual.shape==expected.shape
    want=list(expected.to_Matrix())
    assert len(actual.flat())==len(want)
    assert all(sp.cancel(a-b)==0 for a,b in zip(exact_entries(actual),want))


@pytest.mark.parametrize('kind',('QQ','RF','multi','ZRF'))
@pytest.mark.parametrize('seed',range(12))
def test_reference_full_rref_canonical_basis_and_exact_certificates(monkeypatch,kind,seed):
    rng=random.Random(192581+seed);k,rk=fields(kind)
    rows=seed%4;cols=1+(3*seed)%4
    raw=[]
    for i in range(rows):
        for j in range(cols):
            a=sp.Rational(rng.randrange(-7,8),rng.randrange(1,8))
            if kind!='QQ':a=(a*ry+rng.randrange(-3,4))/(ry+2)
            if kind in ('multi','nested'):a+=(j+1)/(rz+3)
            raw.append(a)
    if raw:raw[0]=sp.S.Zero
    if rows>1 and seed%3==0:raw[-cols:]=raw[:cols]
    a=PolyMatrix(rows,cols,[sx.from_sympy(e) for e in raw],x,ring=k)
    reference=RefPolyMatrix(rows,cols,raw,ring=rk[rx])
    before=a.flat();ids=[id(e) for e in before]
    calls=[];native=sx._core._polymatrix_ground
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_polymatrix_ground',spy)
    reduced,pivots=a.rref();basis=a.nullspace()
    expected,ref_pivots=reference.rref();ref_basis=reference.nullspace()
    assert calls and all(calls)
    assert type(reduced) is PolyMatrix and type(pivots) is tuple
    assert reduced.domain==k and reduced.gens==(x,)
    assert pivots==ref_pivots
    same(reduced,expected)
    assert len(basis)==len(ref_basis)==cols-len(pivots)
    original=sp.Matrix(rows,cols,raw)
    for actual,want in zip(basis,ref_basis):
        assert type(actual) is PolyMatrix and actual.domain==k and actual.gens==(x,)
        same(actual,want)
        vector=sp.Matrix(exact_entries(actual))
        assert all(sp.cancel(v)==0 for v in original*vector)
    assert ids==[id(e) for e in a.flat()] and before==a.flat()
    # The original complete Python algorithm remains a second exact oracle.
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    old,old_pivots=a.rref();old_basis=a.nullspace()
    assert pivots==old_pivots and reduced==old
    assert len(basis)==len(old_basis) and all(b==c for b,c in zip(basis,old_basis))


@pytest.mark.parametrize('rows,cols',[(0,0),(0,4),(3,0),(1,1),(3,5)])
def test_empty_zero_and_full_rank_spaces_are_proven_not_unhandled(rows,cols):
    a=PolyMatrix(rows,cols,[sx.S.Zero]*(rows*cols),x,ring=sx.QQ)
    assert pm._native_ground(a) is not None and pm._native_ground(a,True) is not None
    result,pivots=a.rref()
    assert result.shape==(rows,cols) and pivots==()
    assert len(a.nullspace())==cols
    if rows==cols==1:
        full=PolyMatrix([[1]],x,ring=sx.QQ)
        assert pm._native_ground(full,True)==([], [0])
        assert full.nullspace()==[] and full.rank()==1


@pytest.mark.parametrize('cols',(2**63,2**80))
def test_zero_row_rref_preserves_huge_dimension_fallback(cols):
    a=PolyMatrix(0,cols,[],x,ring=sx.QQ)
    reduced,pivots=a.rref()
    assert reduced.shape==(0,cols) and reduced.domain==sx.QQ and pivots==()


@pytest.mark.parametrize('name',('_entry','_cancel','_like'))
def test_matrix_hooks_keep_the_old_number_and_order_of_calls(monkeypatch,name):
    a=PolyMatrix([[0,1,2],[1,2,3]],x,ring=sx.QQ);calls=[]
    target=pm if name=='_cancel' else PolyMatrix
    original=getattr(target,name)
    def hook(*args,**kwargs):
        calls.append((name,len(args)));return original(*args,**kwargs)
    monkeypatch.setattr(target,name,hook)
    result=a.nullspace();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.nullspace()
    assert actual==calls and calls
    assert result==expected


def test_public_cancel_callback_is_not_run_speculatively_or_skipped(monkeypatch):
    a=PolyMatrix([[1,2,3],[2,4,6]],x,ring=sx.QQ);calls=[];original=sx.cancel
    def hook(value,*args,**kwargs):
        calls.append(str(value));return original(value,*args,**kwargs)
    monkeypatch.setattr(sx,'cancel',hook)
    assert pm._native_ground(a) is None and calls==[]
    result=a.rref();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.rref()
    assert actual==calls and calls and result==expected


def test_wrapped_cancel_closure_callbacks_keep_the_full_route(monkeypatch):
    a=PolyMatrix([[1,2]],x,ring=sx.QQ);calls=[]
    function=sx.cancel
    index=function.__code__.co_freevars.index('native');cell=function.__closure__[index]
    original=cell.cell_contents
    def hook(*args,**kwargs):
        calls.append(1);return original(*args,**kwargs)
    monkeypatch.setattr(cell,'cell_contents',hook)
    assert pm._native_ground(a) is None and calls==[]
    result=a.rref();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.rref()
    assert actual==calls and calls
    assert result==expected


@pytest.mark.parametrize('name',('range','zip','sorted','set'))
def test_module_builtin_callbacks_keep_legacy_order(monkeypatch,name):
    import builtins
    a=PolyMatrix([[1,2,3],[2,4,6]],x,ring=sx.QQ);calls=[];original=getattr(builtins,name)
    def hook(*args):
        calls.append(name);return original(*args)
    monkeypatch.setattr(pm,name,hook,raising=False)
    assert pm._native_ground(a,True) is None and calls==[]
    result=a.nullspace();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.nullspace()
    assert actual==calls and calls
    assert result==expected


def test_expr_operator_callbacks_are_not_skipped(monkeypatch):
    a=PolyMatrix([[1,2],[3,7]],x,ring=sx.QQ);calls=[];original=sx.Expr.__truediv__
    def hook(self,other):
        calls.append(1);return original(self,other)
    monkeypatch.setattr(sx.Expr,'__truediv__',hook)
    assert pm._native_ground(a) is None and calls==[]
    result=a.rref();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.rref()
    assert actual==calls and calls
    assert result==expected


def test_subclasses_and_list_subclasses_keep_their_protocols(monkeypatch):
    calls=[]
    class Custom(PolyMatrix):
        def _rref_rows(self):
            calls.append('rref');return super()._rref_rows()
    a=Custom([[1,2,3]],x,ring=sx.QQ)
    assert pm._native_ground(a) is None and calls==[]
    result=a.nullspace();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.nullspace()
    assert actual==calls==['rref']
    assert result==expected


@pytest.mark.parametrize('k',(sx.RR,sx.CC,sx.RealField(139),sx.ComplexField(211),sx.QQ_I,sx.QQ.algebraic_field(sx.sqrt(2)),sx.EX))
def test_other_domains_retain_complete_existing_calculation(monkeypatch,k):
    a=PolyMatrix([[1,2,3],[2,4,6]],x,ring=k)
    assert pm._native_ground(a) is None and pm._native_ground(a,True) is None
    actual=a.rref();basis=a.nullspace()
    def signature(m):
        return (m.shape,str(m.domain),[(str(e),sorted((f._prec,f._mpf_) for f in e.atoms(sx.Float))) for e in m.flat()])
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.rref();expected_basis=a.nullspace()
    assert actual[1]==expected[1] and signature(actual[0])==signature(expected[0])
    assert [signature(m) for m in basis]==[signature(m) for m in expected_basis]


def test_floating_payload_inside_exact_domain_expression_is_not_rationalized(monkeypatch):
    k=sx.QQ.frac_field(y)
    a=PolyMatrix([[sx.Float('0.5',60)+y,1]],x,ring=k)
    assert pm._native_ground(a) is None
    actual=a.rref()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.rref()
    assert actual==expected
    assert [(f._prec,f._mpf_) for e in actual[0].flat() for f in e.atoms(sx.Float)]==[(f._prec,f._mpf_) for e in expected[0].flat() for f in e.atoms(sx.Float)]


@pytest.mark.parametrize('kind',('nonfield','nonground'))
def test_original_ground_check_errors_are_kept(kind):
    a=PolyMatrix([[1]],x,ring=sx.ZZ) if kind=='nonfield' else PolyMatrix([[x]],x,ring=sx.QQ)
    for name in ('rref','nullspace'):
        with pytest.raises(ValueError,match=f'PolyMatrix {name} is only for ground field elements'):
            getattr(a,name)()


def test_private_binding_refusal_is_not_a_false_empty_nullspace():
    from symox.polys import _polymatrix_native
    stock=_polymatrix_native._ROOT[4];fn=sx._core._polymatrix_ground
    assert fn([x],1,1,sx.QQ,True,stock) is None
    assert fn([sx.S.One],1,2,sx.QQ,True,stock) is None
    assert fn([sx.S.One],1,1,sx.RR,True,stock) is None
    assert fn([],0,2**80,sx.QQ,False,stock) is None
    assert fn([sx.S.One],1,1,sx.QQ,True,stock)==([], [0])



def test_prde_constant_system_uses_selected_integer_rational_function_field(monkeypatch):
    from symox.integrals.prde import constant_system
    from symox.integrals.risch import DifferentialExtension
    from sympy.integrals.prde import constant_system as reference
    from sympy.integrals.risch import DifferentialExtension as RefExtension
    a=PolyMatrix([[-(y+3)/(y-1),(y+1)/(y-1),1],[-y-3,y+1,y-1],[2*(y+3)/(y-1),0,0]],x)
    b=PolyMatrix([[(y+1)/(y-1)],[y+1],[0]],x)
    assert a.domain==sx.ZZ.frac_field(y)
    ar=RefPolyMatrix([[-(ry+3)/(ry-1),(ry+1)/(ry-1),1],[-ry-3,ry+1,ry-1],[2*(ry+3)/(ry-1),0,0]],rx)
    br=RefPolyMatrix([[(ry+1)/(ry-1)],[ry+1],[0]],rx)
    de=DifferentialExtension(extension={'D':[sx.Poly(1,y)]})
    dre=RefExtension(extension={'D':[sp.Poly(1,ry)]})
    calls=[];native=sx._core._polymatrix_ground
    def spy(*args):
        result=native(*args);calls.append((str(args[3]),result is not None));return result
    monkeypatch.setattr(sx._core,'_polymatrix_ground',spy)
    from symox.integrals import _prde_native
    before=_prde_native.STATS['constant_system_field']
    actual=constant_system(a,b,de)
    # Since prdeconst the whole system is computed by the native
    # fraction-field kernel, which does not go through PolyMatrix.rref.
    assert ('ZZ(y)',True) in calls or _prde_native.STATS['constant_system_field']==before+1
    expected=reference(ar,br,dre)
    for value,want in zip(actual,expected):same(value,want)
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    old=constant_system(a,b,de)
    for value,want in zip(actual,old):
        assert value==want and value.domain==want.domain


def test_prde_parametric_polynomial_solution_and_nullspace_certificate(monkeypatch):
    from symox.integrals.prde import param_poly_rischDE
    from symox.integrals.risch import DifferentialExtension
    from sympy.integrals.prde import param_poly_rischDE as reference
    from sympy.integrals.risch import DifferentialExtension as RefExtension
    de=DifferentialExtension(extension={'D':[sx.Poly(1,y)]})
    dre=RefExtension(extension={'D':[sp.Poly(1,ry)]})
    a=sx.Poly(y*y-y,y,field=True);b=sx.Poly(1,y,field=True)
    q=[sx.Poly(y,y,field=True),sx.Poly(y*y,y,field=True)]
    h,matrix=param_poly_rischDE(a,b,q,3,de)
    rh,reference_matrix=reference(sp.Poly(ry*ry-ry,ry,field=True),sp.Poly(1,ry,field=True),
                                   [sp.Poly(ry,ry,field=True),sp.Poly(ry*ry,ry,field=True)],3,dre)
    same(matrix,reference_matrix)
    assert len(h)==len(rh)
    assert all(sp.cancel(sx.to_sympy(f.as_expr())-r.as_expr())==0 for f,r in zip(h,rh))
    calls=[];native=sx._core._polymatrix_ground
    def spy(*args):
        result=native(*args);calls.append(result is not None);return result
    monkeypatch.setattr(sx._core,'_polymatrix_ground',spy)
    basis=matrix.nullspace()
    assert calls and all(calls) and len(basis)==1
    same(basis[0],reference_matrix.nullspace()[0])
    assert all(sp.cancel(v)==0 for v in sp.Matrix(matrix.rows,matrix.cols,exact_entries(matrix))*sp.Matrix(exact_entries(basis[0])))
    assert h[0]+h[1]==sx.Poly(y,y,domain=sx.QQ)


def test_list_subclass_storage_callbacks_are_not_speculated_or_skipped(monkeypatch):
    calls=[]
    class Items(list):
        def __getitem__(self,key):
            calls.append(str(key));return super().__getitem__(key)
    a=PolyMatrix([[1,2,3],[2,4,6]],x,ring=sx.QQ)
    a._items=Items(a._items)
    assert pm._native_ground(a) is None and calls==[]
    result=a.rref();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.rref()
    assert actual==calls and calls
    assert result==expected


def test_instance_reduction_override_is_respected(monkeypatch):
    a=PolyMatrix([[1,2,3]],x,ring=sx.QQ);calls=[]
    original=a._rref_rows
    def custom():
        calls.append(1);return original()
    a._rref_rows=custom
    assert pm._native_ground(a,True) is None and calls==[]
    result=a.nullspace()
    assert calls==[1]
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    calls.clear();expected=a.nullspace()
    assert calls==[1] and result==expected


def test_ignored_rref_options_do_not_call_user_functions():
    a=PolyMatrix([[1,2],[3,7]],x,ring=sx.QQ)
    def forbidden(*args):raise AssertionError('option is ignored by the established PolyMatrix API')
    reduced,pivots=a.rref(simplify=forbidden,normalize_last=forbidden)
    assert exact_entries(reduced)==[1,0,0,1] and pivots==(0,1)


def test_sympify_delegate_callbacks_used_for_basis_constants_are_preserved(monkeypatch):
    module=importlib.import_module('symox.core.sympify')
    a=PolyMatrix([[1,2,3]],x,ring=sx.QQ);calls=[];original=module.sympify
    def hook(*args,**kwargs):
        calls.append(1);return original(*args,**kwargs)
    monkeypatch.setattr(module,'sympify',hook)
    assert pm._native_ground(a,True) is None and calls==[]
    result=a.nullspace();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.nullspace()
    assert actual==calls and calls
    assert result==expected


def test_registered_integer_converter_callbacks_are_preserved(monkeypatch):
    from symox.core.sympify import converter
    a=PolyMatrix([[1,2,3]],x,ring=sx.QQ);calls=[]
    def convert(value):
        calls.append(value);return sx.Integer(value)
    monkeypatch.setitem(converter,int,convert)
    assert pm._native_ground(a,True) is None and calls==[]
    result=a.nullspace();actual=calls[:];calls.clear()
    monkeypatch.setattr(pm,'_native_ground',lambda *args:None)
    expected=a.nullspace()
    assert actual==calls and calls
    assert result==expected


@pytest.mark.parametrize('kind',('QQ','RF'))
def test_large_exact_coefficient_limbs_are_preserved_in_reduction_and_basis(kind):
    k,rk=fields(kind);large=(1<<521)+3**147
    coefficient=sp.Integer(large)+(ry if kind=='RF' else 0)
    raw=[coefficient,sp.S.One,sp.Integer(2),2*coefficient,sp.Integer(2),sp.Integer(4)]
    a=PolyMatrix(2,3,[sx.from_sympy(e) for e in raw],x,ring=k)
    reference=RefPolyMatrix(2,3,raw,ring=rk[rx])
    assert pm._native_ground(a) is not None
    reduced,pivots=a.rref();expected,want=reference.rref()
    assert pivots==want==(0,)
    same(reduced,expected)
    basis=a.nullspace()
    assert len(basis)==2
    for value,ref in zip(basis,reference.nullspace()):
        same(value,ref)
        assert all(sp.cancel(e)==0 for e in sp.Matrix(2,3,raw)*sp.Matrix(exact_entries(value)))
