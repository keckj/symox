"""Exact LU certificates, mutation conventions and complete protocol fallbacks."""
import ast
import builtins
import copy
import ctypes
import importlib
import random
import types
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp
import symox as sx
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix
from sympy.polys.matrices.dense import ddm_ilu as reference

dense=importlib.import_module('symox.polys.matrices.dense')
helper=importlib.import_module('symox.polys.matrices._lu_native')
algebraic=importlib.import_module('symox.polys.domains._algebraic_methods')


def qq(rows):
    return [[sx.QQ(Fraction(value).numerator,Fraction(value).denominator) for value in row] for row in rows]


def rational_rows(rows):
    return [[Fraction(str(value)) for value in row] for row in rows]


def certificate(source,combined,swaps):
    m=len(source);n=len(source[0]) if m else 0
    permuted=[list(map(Fraction,row)) for row in source]
    for i,j in swaps:permuted[i],permuted[j]=permuted[j],permuted[i]
    lower=[[Fraction(int(i==j)) if i<=j else combined[i][j] if j<n else Fraction(0) for j in range(m)] for i in range(m)]
    upper=[[combined[i][j] if i<=j else Fraction(0) for j in range(n)] for i in range(m)]
    assert [[sum(lower[i][k]*upper[k][j] for k in range(m)) for j in range(n)] for i in range(m)]==permuted


def old_function():
    # The original body is also available after native dispatch, so execute its
    # unchanged statements against the actual module globals for hook tests.
    module=ast.parse(Path(dense.__file__).read_text())
    function=next(node for node in module.body if isinstance(node,ast.FunctionDef) and node.name=='ddm_ilu')
    start=next(i for i,node in enumerate(function.body) if isinstance(node,ast.Assign)
               and any(isinstance(t,ast.Name) and t.id=='m' for t in node.targets))
    function.body=[node for node in function.body[start:]
                   if not (isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='native' for t in node.targets))
                   and not (isinstance(node,ast.If) and isinstance(node.test,ast.Compare)
                            and isinstance(node.test.left,ast.Name) and node.test.left.id=='native')]
    compiled=compile(ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[])),dense.__file__,'exec')
    code=next(c for c in compiled.co_consts if type(c) is types.CodeType)
    return types.FunctionType(code,dense.__dict__)


@pytest.mark.parametrize('shape',[(1,1),(1,5),(2,2),(3,2),(2,5),(4,4),(4,6),(6,4)])
@pytest.mark.parametrize('seed',range(8))
def test_native_exact_values_swaps_publication_and_lattice_certificate(shape,seed):
    m,n=shape;rng=random.Random(79139+seed*811+m*17+n)
    source=[[Fraction(rng.randrange(-8,9),rng.randrange(1,13)) for _ in range(n)] for _ in range(m)]
    if seed%3==0:source[0][0]=Fraction(0)
    if seed%4==0 and m>1:source[-1]=source[0][:]
    if seed%5==0:
        for row in source:row[0]=Fraction(0)
    a=qq(source);b=qq(source)
    a_rows={id(row):i for i,row in enumerate(a)};b_rows={id(row):i for i,row in enumerate(b)}
    a_cells=[[v for v in row] for row in a];b_cells=[[v for v in row] for row in b]
    expected=reference(b);actual=sx._core._ddm_ilu(a)
    assert actual==expected and a==b
    assert [a_rows[id(row)] for row in a]==[b_rows[id(row)] for row in b]
    for left,right in zip(a,b):
        origin=a_rows[id(left)]
        assert [(v is a_cells[origin][j]) for j,v in enumerate(left)]==[(v is b_cells[origin][j]) for j,v in enumerate(right)]
    assert all(type(value) is sx._core.DomainElement and value.parent()==sx.QQ for row in a for value in row)
    certificate(source,rational_rows(a),actual)


@pytest.mark.parametrize('shape',[(2,2),(3,3),(4,2),(2,5)])
@pytest.mark.parametrize('entry',['ddm','domainmatrix_dense','domainmatrix_sparse'])
def test_actual_public_lu_first_call_uses_native_and_preserves_input(shape,entry,monkeypatch):
    m,n=shape;raw=[[Fraction((i+1)*(j+2)+(i==j),i+j+1) for j in range(n)] for i in range(m)]
    calls=[];native=sx._core._ddm_ilu
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_ilu',record)
    a=DDM(qq(raw),shape,sx.QQ)
    if entry!='ddm':
        a=DomainMatrix.from_rep(a)
        if entry.endswith('sparse'):a=a.to_sparse()
    before=rational_rows(a if entry=='ddm' else a.to_list())
    lower,upper,swaps=a.lu()
    lower=sp.Matrix([[sp.Rational(str(v)) for v in row] for row in (lower if entry=='ddm' else lower.to_list())])
    upper=sp.Matrix([[sp.Rational(str(v)) for v in row] for row in (upper if entry=='ddm' else upper.to_list())])
    permuted=sp.Matrix(raw)
    for i,j in swaps:permuted.row_swap(i,j)
    assert lower*upper==permuted
    assert rational_rows(a if entry=='ddm' else a.to_list())==before
    assert calls and all(result is not None for result in calls)


def test_public_lu_solve_keeps_exact_two_sided_system_certificate(monkeypatch):
    calls=[];native=sx._core._ddm_ilu
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_ilu',record)
    a=DomainMatrix.from_list([[0,2,1],[1,3,4],[4,1,2]],sx.QQ)
    b=DomainMatrix.from_list([[2,3],[4,5],[6,7]],sx.QQ)
    x=a.lu_solve(b)
    assert a*x==b and calls and all(result is not None for result in calls)


@pytest.mark.parametrize('rows',[[],[[]],[[],[]]])
def test_empty_shapes_are_returned_without_row_replacement(rows):
    rows=copy.deepcopy(rows);old=list(rows)
    assert sx._core._ddm_ilu(rows)==[]
    assert all(left is right for left,right in zip(rows,old))


@pytest.mark.parametrize('source',[
    [[1,2],[3,4]],[[sx.ZZ(2),sx.ZZ(4)],[sx.ZZ(1),sx.ZZ(3)]],
    [[sx.QQ(2),0],[sx.QQ(3),sx.QQ(4)]],
    [[sx.QQ(2),sx.QQ(3)],[sx.QQ(4)]],
    [[],[sx.QQ(1)]],
])
def test_other_representations_and_ragged_shapes_keep_existing_behavior(source):
    a=copy.deepcopy(source);b=copy.deepcopy(source);before=copy.deepcopy(source)
    assert sx._core._ddm_ilu(a) is None and a==before
    try:expected=reference(b)
    except Exception as error:
        with pytest.raises(type(error)) as actual:dense.ddm_ilu(a)
        assert str(actual.value)==str(error)
    else:assert dense.ddm_ilu(a)==expected
    assert a==b
    assert [[type(v) for v in row] for row in a]==[[type(v) for v in row] for row in b]


def test_aliased_rows_are_declined_without_mutation_then_preserve_feedback():
    row=qq([[2,3]])[0];a=[row,row];oldrow=row[:];b=[oldrow,oldrow]
    assert sx._core._ddm_ilu(a) is None and a[0] is a[1] is row
    assert dense.ddm_ilu(a)==reference(b) and a==b and a[0] is a[1] is row


@pytest.mark.parametrize('name',['len','range','min'])
def test_actual_builtin_overrides_are_not_skipped(name,monkeypatch):
    calls=[];original=getattr(__import__('builtins'),name);old=old_function()
    def hook(*args):calls.append(args);return original(*args)
    monkeypatch.setattr(dense,name,hook,raising=False)
    a=qq([[0,2],[3,4]]);expected=old(a);events=calls[:];calls.clear()
    b=qq([[0,2],[3,4]]);actual=dense.ddm_ilu(b)
    assert actual==expected and a==b and calls==events


@pytest.mark.parametrize('name',['__bool__','__mul__','__sub__','parent','__truediv__','__isub__'])
def test_coefficient_callbacks_and_added_inplace_hook_are_preserved(name,monkeypatch):
    old=old_function();calls=[];original=getattr(sx._core.DomainElement,name,None)
    if original is None:original=sx._core.DomainElement.__sub__
    def hook(*args):calls.append(name);return original(*args)
    monkeypatch.setattr(sx._core.DomainElement,name,hook,raising=False)
    a=qq([[2,3],[4,7]]);expected=old(a);events=calls[:];calls.clear()
    b=qq([[2,3],[4,7]]);actual=dense.ddm_ilu(b)
    assert actual==expected and a==b and calls==events and calls


@pytest.mark.parametrize('name',['is_ZZ_I','is_QQ_I'])
def test_domain_flag_access_callbacks_are_preserved(name,monkeypatch):
    old=old_function();calls=[];original=sx._core.Domain.__dict__[name]
    def flag(value):calls.append(name);return original.__get__(value)
    monkeypatch.setattr(sx._core.Domain,name,property(flag))
    a=qq([[2,3],[4,7]]);expected=old(a);events=calls[:];calls.clear()
    b=qq([[2,3],[4,7]]);actual=dense.ddm_ilu(b)
    assert actual==expected and a==b and calls==events and calls


def test_native_division_delegate_cell_change_keeps_original_semantics(monkeypatch):
    divide=sx._core.DomainElement.__truediv__;cell=divide.__closure__[0];original=cell.cell_contents
    old=old_function();calls=[]
    def hook(*args):calls.append('divide');return original(*args)
    cell.cell_contents=hook
    try:
        a=qq([[2,3],[4,7]]);expected=old(a);events=calls[:];calls.clear()
        b=qq([[2,3],[4,7]]);actual=dense.ddm_ilu(b)
        assert actual==expected and a==b and calls==events and calls
    finally:cell.cell_contents=original


def test_changed_stock_ddm_list_assignment_and_custom_rows_are_not_bypassed(monkeypatch):
    a=DDM(qq([[0,2],[3,4]]),(2,2),sx.QQ)
    def hook(*args):raise KeyboardInterrupt('row publication')
    monkeypatch.setattr(DDM,'__setitem__',hook)
    with pytest.raises(KeyboardInterrupt,match='row publication'):dense.ddm_ilu(a)
    class Row(list):
        def __len__(self):raise KeyboardInterrupt('row length')
    row=Row(qq([[1]])[0])
    assert sx._core._ddm_ilu([row]) is None
    with pytest.raises(KeyboardInterrupt,match='row length'):dense.ddm_ilu([row])


def test_reload_does_not_authenticate_replaced_len(monkeypatch):
    def hook(*args):raise KeyboardInterrupt('len override')
    monkeypatch.setattr(dense,'len',hook,raising=False)
    importlib.reload(helper)
    with pytest.raises(KeyboardInterrupt,match='len override'):dense.ddm_ilu(qq([[1]]))


@pytest.mark.parametrize('digits',[100,500,1500])
def test_large_exact_rational_pivots_and_nearly_singular_inputs(digits):
    big=10**digits;tiny=Fraction(1,big)
    source=[[Fraction(big,3),Fraction(1,7),tiny],
            [Fraction(big,3),Fraction(1,7)+tiny,Fraction(5,9)],
            [Fraction(0),Fraction(4,11),Fraction(big,17)]]
    a=qq(source);b=qq(source)
    assert dense.ddm_ilu(a)==reference(b) and a==b
    certificate(source,rational_rows(a),reference(qq(source)))


def test_domain_and_coefficient_attribute_callbacks_keep_their_actual_order(monkeypatch):
    for cls in (sx._core.Domain,sx._core.DomainElement):
        old=old_function();calls=[];native=object.__getattribute__
        def observe(value,name):
            calls.append(name)
            return native(value,name)
        with monkeypatch.context() as context:
            context.setattr(cls,'__getattribute__',observe)
            a=qq([[2,3],[4,7]]);calls.clear();expected=old(a);events=calls[:]
            b=qq([[2,3],[4,7]]);calls.clear();actual=dense.ddm_ilu(b)
            assert actual==expected and a==b and calls==events and calls


def test_native_admission_does_not_add_module_globals_callback(monkeypatch):
    def forbidden():raise AssertionError('new globals callback')
    monkeypatch.setattr(dense,'globals',forbidden,raising=False)
    a=qq([[2,3],[4,7]])
    assert dense.ddm_ilu(a)==[]
    assert rational_rows(a)==[[Fraction(2),Fraction(3)],[Fraction(2),Fraction(1)]]


def test_presentation_hooks_and_real_units_import_do_not_disable_unrelated_lu(monkeypatch):
    import symox.physics.units
    calls=[];native=sx._core._ddm_ilu
    def record(*args):
        result=native(*args);calls.append(result);return result
    def forbidden(*args):raise AssertionError('not a consumed operation')
    monkeypatch.setattr(sx._core,'_ddm_ilu',record)
    monkeypatch.setattr(sx._core.DomainElement,'__repr__',forbidden)
    monkeypatch.setattr(sx._core.DomainElement,'__str__',forbidden)
    a=qq([[2,3],[4,7]])
    assert dense.ddm_ilu(a)==[] and calls==[[]]


class _MethodDef(ctypes.Structure):
    _fields_=[('name',ctypes.c_char_p),('function',ctypes.c_void_p),('flags',ctypes.c_int),('doc',ctypes.c_char_p)]


_method_defs=[]


def foreign_builtin(name,owner):
    original=getattr(builtins,name)
    get_function=ctypes.pythonapi.PyCFunction_GetFunction
    get_function.argtypes=[ctypes.py_object];get_function.restype=ctypes.c_void_p
    get_flags=ctypes.pythonapi.PyCFunction_GetFlags
    get_flags.argtypes=[ctypes.py_object];get_flags.restype=ctypes.c_int
    definition=_MethodDef(name.encode(),get_function(original),get_flags(original),None)
    _method_defs.append(definition)
    create=ctypes.pythonapi.PyCFunction_NewEx
    create.argtypes=[ctypes.POINTER(_MethodDef),ctypes.c_void_p,ctypes.c_void_p];create.restype=ctypes.py_object
    return create(ctypes.byref(definition),id(owner),id(None))


@pytest.mark.parametrize('name',['len','min'])
@pytest.mark.parametrize('owner_kind',['none','foreign','hostile'])
def test_foreign_builtin_at_helper_bootstrap_is_not_authenticated(name,owner_kind,monkeypatch):
    events=[]
    class Hostile(types.ModuleType):
        def __getattribute__(self,name):events.append(name);raise AssertionError('module attribute hook')
    owner=None if owner_kind=='none' else types.ModuleType('foreign') if owner_kind=='foreign' else Hostile('foreign')
    replaced=foreign_builtin(name,owner)
    fresh=types.ModuleType('symox.polys.matrices._lu_native_fresh');fresh.__package__='symox.polys.matrices'
    code=compile(Path(helper.__file__).read_text(),helper.__file__,'exec')
    with monkeypatch.context() as context:
        context.setattr(builtins,name,replaced)
        exec(code,fresh.__dict__)
    a=qq([[2,3],[4,7]]);before=copy.deepcopy(a)
    assert fresh.try_ilu(a,(range,min),None,(2,2)) is None and a==before and events==[]


@pytest.mark.parametrize('name',['len','min','range'])
def test_python_builtin_replacements_at_helper_bootstrap_are_not_authenticated(name,monkeypatch):
    original=getattr(builtins,name)
    if name=='range':
        def create(cls,*args):return original(*args)
        replaced=type('range',(),{'__module__':'builtins','__new__':create})
    else:
        def replaced(*args):return original(*args)
    fresh=types.ModuleType('symox.polys.matrices._lu_native_fresh');fresh.__package__='symox.polys.matrices'
    code=compile(Path(helper.__file__).read_text(),helper.__file__,'exec')
    with monkeypatch.context() as context:
        context.setattr(builtins,name,replaced)
        exec(code,fresh.__dict__)
    assert fresh.try_ilu(qq([[2,3],[4,7]]),(range,min),None,(2,2)) is None


@pytest.mark.parametrize('changed',['len','range','min'])
def test_cloned_function_uses_its_actual_globals_and_preserves_callbacks(changed):
    calls=[];real=getattr(builtins,changed)
    def hook(*args):calls.append(changed);return real(*args)
    namespace=dict(dense.__dict__);namespace[changed]=hook
    current=types.FunctionType(dense.ddm_ilu.__code__,namespace)
    old=types.FunctionType(old_function().__code__,namespace)
    a=qq([[2,3],[4,7]]);expected=old(a);events=calls[:];calls.clear()
    b=qq([[2,3],[4,7]]);actual=current(b)
    assert actual==expected and a==b and calls==events and calls


@pytest.mark.parametrize('wrong_dimension',['rows','columns','empty','negative'])
def test_real_length_prefix_is_not_replayed_and_false_dimensions_decline(wrong_dimension,monkeypatch):
    events=[]
    def hook(value):
        index=len(events);events.append(index)
        if wrong_dimension=='empty' and index==0:return 0
        if wrong_dimension=='negative' and index==0:return -1
        if wrong_dimension=='rows' and index==0:return len(value)-1
        if wrong_dimension=='columns' and index==1:return len(value)-1
        return len(value)
    monkeypatch.setattr(dense,'len',hook,raising=False)
    source=[[2,3,4],[5,6,7],[8,9,12]]
    a=qq(source);expected=old_function()(a);calls=events[:];events.clear()
    b=qq(source);actual=dense.ddm_ilu(b)
    assert actual==expected and a==b and events==calls
    assert len(events)==(1 if wrong_dimension=='empty' else 2)


def test_length_callback_can_remove_itself_without_an_extra_global_lookup():
    events=[];namespace=dict(dense.__dict__)
    def hook(value):
        events.append('len')
        if len(events)==2:namespace.pop('len')
        return len(value)
    namespace['len']=hook
    current=types.FunctionType(dense.ddm_ilu.__code__,namespace)
    a=qq([[2,3],[4,7]])
    assert current(a)==[] and events==['len','len']


@pytest.mark.parametrize('missing',['len','range','min'])
def test_missing_builtin_keeps_original_name_error_and_mutation_prefix(missing):
    namespace=dict(dense.__dict__);lookup=dict(vars(builtins));lookup.pop(missing)
    namespace['__builtins__']=lookup;namespace.pop(missing,None)
    current=types.FunctionType(dense.ddm_ilu.__code__,namespace)
    old=types.FunctionType(old_function().__code__,namespace)
    a=qq([[2,3],[4,7]]);b=qq([[2,3],[4,7]])
    with pytest.raises(NameError) as expected:old(a)
    with pytest.raises(NameError) as actual:current(b)
    assert str(actual.value)==str(expected.value) and a==b


@pytest.mark.parametrize('dimensions',[(1,2),(2,1),(-1,2),(2,-1),(True,2),(2,True)])
def test_native_refuses_false_prefix_dimensions_before_any_publication(dimensions):
    a=qq([[2,3],[4,7]]);before=copy.deepcopy(a)
    assert sx._core._ddm_ilu(a,None,dimensions) is None and a==before
