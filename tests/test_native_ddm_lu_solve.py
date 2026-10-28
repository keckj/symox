"""QQ triangular solve: exact systems, publication prefixes and real hooks."""
import ast
import builtins
import copy
import importlib
import inspect
import random
import types
from fractions import Fraction
from pathlib import Path

import pytest
import symox as sx
from symox.polys.matrices.ddm import DDM
from symox.polys.matrices.domainmatrix import DomainMatrix
from symox.polys.matrices.exceptions import DMNonInvertibleMatrixError
from sympy.polys.matrices.dense import ddm_ilu_solve as reference

dense=importlib.import_module('symox.polys.matrices.dense')


def qq(rows):
    return [[sx.QQ(Fraction(v).numerator,Fraction(v).denominator) for v in row] for row in rows]


def values(rows):
    return [[None if v is None else Fraction(str(v)) for v in row] for row in rows]


def old_function(namespace=None):
    tree=ast.parse(Path(dense.__file__).read_text())
    function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='ddm_ilu_solve')
    function.body=[n for n in function.body if not (isinstance(n,ast.If)
                   and isinstance(n.test,ast.Call) and isinstance(n.test.func,ast.Name)
                   and n.test.func.id=='_native_ilu_solve')]
    code=compile(ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[])),dense.__file__,'exec')
    return types.FunctionType(next(c for c in code.co_consts if type(c)is types.CodeType),dense.__dict__ if namespace is None else namespace)


def factors(source):
    m=len(source);n=len(source[0]) if m else 0
    return DDM(qq(source),(m,n),sx.QQ).lu()


@pytest.mark.parametrize('shape',[(1,1),(2,2),(3,3),(5,5),(4,2),(5,3)])
@pytest.mark.parametrize('rhs_columns',[0,1,3])
@pytest.mark.parametrize('seed',range(8))
def test_exact_solution_and_multiple_rhs_match_reference_with_independent_certificate(shape,rhs_columns,seed,monkeypatch):
    m,n=shape;rng=random.Random(90443+seed*29+m*17+n)
    a=[[Fraction((i==j)*20+rng.randrange(-5,6),rng.randrange(1,7)) for j in range(n)] for i in range(m)]
    known=[[Fraction(rng.randrange(-9,10),rng.randrange(1,9)) for _ in range(rhs_columns)] for _ in range(n)]
    b=[[sum(a[i][j]*known[j][k] for j in range(n)) for k in range(rhs_columns)] for i in range(m)]
    lower,upper,swaps=factors(a);rhs=qq(b)
    x=[[None]*rhs_columns for _ in range(n)];expected=[[None]*rhs_columns for _ in range(n)]
    row_ids=list(x);calls=[];native=sx._core._ddm_ilu_solve
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_ilu_solve',record)
    reference(expected,lower,upper,swaps,rhs)
    assert dense.ddm_ilu_solve(x,lower,upper,swaps,rhs) is None
    assert values(x)==values(expected)==known and calls==[True]
    assert all(a is b for a,b in zip(x,row_ids))
    assert [[sum(a[i][j]*values(x)[j][k] for j in range(n)) for k in range(rhs_columns)] for i in range(m)]==b


@pytest.mark.parametrize('entry',['ddm','dense','sparse'])
@pytest.mark.parametrize('shape',[(2,2),(3,3),(4,2)])
def test_real_public_solve_uses_native_substitution_and_preserves_operands(entry,shape,monkeypatch):
    m,n=shape;raw=[[int(i==j)*10+(i+1)*(j+2) for j in range(n)] for i in range(m)]
    rhs=[[sum(raw[i][j]*(j+k+1) for j in range(n)) for k in range(2)] for i in range(m)]
    a=DDM(qq(raw),shape,sx.QQ);b=DDM(qq(rhs),(m,2),sx.QQ)
    if entry!='ddm':
        a=DomainMatrix.from_rep(a);b=DomainMatrix.from_rep(b)
        if entry=='sparse':a=a.to_sparse();b=b.to_sparse()
    calls=[];native=sx._core._ddm_ilu_solve
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_ilu_solve',record)
    out=a.lu_solve(b)
    rows=out if entry=='ddm' else out.to_list()
    assert values(rows)==[[Fraction(j+k+1) for k in range(2)] for j in range(n)] and calls==[True]
    assert values(a if entry=='ddm' else a.to_list())==raw
    assert values(b if entry=='ddm' else b.to_list())==rhs


def test_combined_lu_storage_and_swapped_rhs_keep_the_original_conventions():
    combined=qq([[1,2],[3,-2]]);rhs=qq([[5],[6]]);x=[[None],[None]]
    assert dense.ddm_ilu_solve(x,combined,combined,[],rhs) is None
    assert values(x)==[[Fraction(-4)],[Fraction(9,2)]]


@pytest.mark.parametrize('case',['inconsistent','singular_first','singular_last','zero_columns','zero_rhs'])
def test_noninvertibility_and_completed_x_prefix_match_old_function(case):
    if case=='inconsistent':l=qq([[1,0,0],[2,1,0],[3,4,1]]);u=qq([[2,1],[0,3],[0,0]]);b=qq([[4],[14],[37]]);x=[[None],[None]]
    elif case=='zero_columns':l=qq([[1,0],[0,1]]);u=[[],[]];b=qq([[0],[1]]);x=[]
    elif case=='zero_rhs':l=qq([[1,0],[0,1]]);u=qq([[0,1],[0,0]]);b=[[],[]];x=[[],[]]
    else:l=qq([[1,0],[0,1]]);u=qq([[0,1],[0,2]] if case=='singular_first' else [[2,1],[0,0]]);b=qq([[1,2],[4,6]]);x=[[None,None],[None,None]]
    old_x=copy.deepcopy(x)
    try:expected=old_function()(old_x,l,u,[],b)
    except Exception as error:
        with pytest.raises(type(error)) as actual:dense.ddm_ilu_solve(x,l,u,[],b)
        assert str(actual.value)==str(error)
    else:assert dense.ddm_ilu_solve(x,l,u,[],b)==expected
    assert values(x)==values(old_x)


@pytest.mark.parametrize('kind',['mismatch','empty_b','underdetermined','ragged_x','ragged_l','ragged_u','ragged_b'])
def test_shape_failures_preserve_order_message_and_mutation_prefix(kind):
    l=qq([[1,0],[2,1]]);u=qq([[2,1],[0,3]]);b=qq([[1],[2]]);x=[[None],[None]]
    if kind=='mismatch':b=b[:1]
    elif kind=='empty_b':b=[]
    elif kind=='underdetermined':u=qq([[1,2,3],[0,1,2]])
    elif kind=='ragged_x':x=[[None],[]]
    elif kind=='ragged_l':l=qq([[1],[2]])
    elif kind=='ragged_u':u=qq([[2,1],[0]])
    elif kind=='ragged_b':b=qq([[1],[]])
    old_x=copy.deepcopy(x)
    try:expected=old_function()(old_x,l,u,[],b)
    except Exception as error:
        with pytest.raises(type(error)) as actual:dense.ddm_ilu_solve(x,l,u,[],b)
        assert str(actual.value)==str(error)
    else:assert dense.ddm_ilu_solve(x,l,u,[],b)==expected
    assert values(x)==values(old_x)


@pytest.mark.parametrize('alias',['xrows','x_l','x_u','x_b'])
def test_aliases_decline_without_losing_existing_feedback(alias,monkeypatch):
    def inputs():
        l=qq([[1,0],[2,1]]);u=qq([[2,1],[0,3]]);b=qq([[1,2],[3,4]]);x=qq([[0,0],[0,0]])
        if alias=='xrows':x[1]=x[0]
        elif alias=='x_l':x=l
        elif alias=='x_u':x=u
        else:x=b
        return x,l,u,b
    old_x,ol,ou,ob=inputs();x,l,u,b=inputs();calls=[];native=sx._core._ddm_ilu_solve
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_ilu_solve',record)
    try:expected=old_function()(old_x,ol,ou,[],ob)
    except Exception as error:
        with pytest.raises(type(error)):dense.ddm_ilu_solve(x,l,u,[],b)
    else:assert dense.ddm_ilu_solve(x,l,u,[],b)==expected
    assert values(x)==values(old_x) and all(value is False for value in calls)


@pytest.mark.parametrize('name',['len','range','reversed'])
def test_real_builtins_and_cloned_globals_keep_their_actual_callbacks(name):
    events=[];real=getattr(builtins,name)
    def hook(*args):events.append(name);return real(*args)
    namespace=dict(dense.__dict__);namespace[name]=hook
    current=types.FunctionType(dense.ddm_ilu_solve.__code__,namespace);old=old_function(namespace)
    l=qq([[1,0],[2,1]]);u=qq([[2,1],[0,3]]);b=qq([[1,2],[3,4]])
    a=[[None,None],[None,None]];old(a,l,u,[],b);before=events[:];events.clear()
    x=[[None,None],[None,None]];current(x,l,u,[],b)
    assert values(x)==values(a) and events==before and events


@pytest.mark.parametrize('name',['__bool__','__mul__','__sub__','__truediv__','parent','__isub__'])
def test_modified_arithmetic_callbacks_and_added_inplace_subtraction_are_not_skipped(name,monkeypatch):
    events=[];original=getattr(sx._core.DomainElement,name,None) or sx._core.DomainElement.__sub__
    def hook(*args):events.append(name);return original(*args)
    monkeypatch.setattr(sx._core.DomainElement,name,hook,raising=False)
    l=qq([[1,0],[2,1]]);u=qq([[2,1],[0,3]]);b=qq([[1],[3]])
    a=[[None],[None]];old_function()(a,l,u,[],b);before=events[:];events.clear()
    x=[[None],[None]];dense.ddm_ilu_solve(x,l,u,[],b);after=events[:]
    assert values(x)==values(a) and after==before and after


def test_actual_swap_iterator_and_rhs_copy_hooks_are_not_replayed():
    events=[]
    class Pair(tuple):
        def __iter__(self):events.append('pair');return super().__iter__()
    class Swaps(list):
        def __bool__(self):events.append('bool');return True
        def __iter__(self):events.append('iter');return super().__iter__()
    swaps=Swaps([Pair((0,1))]);l=qq([[1,0],[2,1]]);u=qq([[2,1],[0,3]]);b=qq([[1],[3]])
    a=[[None],[None]];old_function()(a,l,u,swaps,b);before=events[:];events.clear()
    x=[[None],[None]];dense.ddm_ilu_solve(x,l,u,swaps,b)
    assert values(x)==values(a) and events==before


def test_custom_exception_constructor_sees_original_completed_scratch_and_partial_x(monkeypatch):
    events=[]
    def construct(self):
        frame=inspect.currentframe().f_back
        while frame is not None and frame.f_code.co_name!='ddm_ilu_solve':frame=frame.f_back
        assert frame is not None
        state=frame.f_locals
        events.append((values(state['y']),values(state['x']),state['y'][0][0] is state['b'][0][0]))
    monkeypatch.setattr(DMNonInvertibleMatrixError,'__init__',construct)
    l=qq([[1,0],[0,1]]);u=qq([[0,1],[0,2]]);b=qq([[1,2],[4,6]])
    with pytest.raises(DMNonInvertibleMatrixError):old_function()([[None,None],[None,None]],l,u,[],b)
    before=events[:];events.clear()
    with pytest.raises(DMNonInvertibleMatrixError):dense.ddm_ilu_solve([[None,None],[None,None]],l,u,[],b)
    assert events==before and events[0][2] is True


def test_overwritten_placeholder_destructors_can_mutate_future_pivots(monkeypatch):
    def run(function):
        events=[];l=qq([[1,0],[0,1]]);u=qq([[1,0],[0,1]]);b=qq([[2],[3]])
        class Token:
            def __del__(self):events.append('delete');u[0][0]=sx.QQ(2)
        x=[[Token()],[Token()]];function(x,l,u,[],b)
        return values(x),events
    assert run(dense.ddm_ilu_solve)==run(old_function())==([[Fraction(1)],[Fraction(3)]],['delete','delete'])


@pytest.mark.parametrize('missing',['reversed','DMNonInvertibleMatrixError'])
def test_missing_late_names_do_not_precede_forward_substitution_or_singularity(missing):
    namespace=dict(dense.__dict__);lookups=dict(vars(builtins));lookups.pop(missing,None)
    namespace['__builtins__']=lookups;namespace.pop(missing,None)
    current=types.FunctionType(dense.ddm_ilu_solve.__code__,namespace);old=old_function(namespace)
    l=qq([[1,0],[0,1]]);u=qq([[0,1],[0,2]]);b=qq([[1,2],[4,6]])
    def failure(function):
        x=[[None,None],[None,None]]
        try:function(x,l,u,[],b)
        except Exception as error:
            trace=error.__traceback__
            while trace is not None and trace.tb_frame.f_code.co_name!='ddm_ilu_solve':trace=trace.tb_next
            return type(error),str(error),values(x),values(trace.tb_frame.f_locals['y'])
        raise AssertionError('expected failure')
    assert failure(current)==failure(old)


@pytest.mark.parametrize('domain',['integers','fractions','ZZ','GF','RR','QQ_I'])
def test_other_coefficient_domains_retain_original_types_and_complete_failures(domain):
    if domain=='integers':convert=int
    elif domain=='fractions':convert=Fraction
    else:convert={'ZZ':sx.ZZ,'GF':sx.GF(7),'RR':sx.RealField(139),'QQ_I':sx.QQ_I}[domain]
    l=[[convert(v) for v in row] for row in [[1,0],[2,1]]]
    u=[[convert(v) for v in row] for row in [[2,1],[0,3]]]
    b=[[convert(v) for v in row] for row in [[1],[3]]]
    a=[[None],[None]];x=[[None],[None]]
    try:expected=old_function()(a,l,u,[],b)
    except Exception as error:
        with pytest.raises(type(error)) as actual:dense.ddm_ilu_solve(x,l,u,[],b)
        assert str(actual.value)==str(error)
    else:assert dense.ddm_ilu_solve(x,l,u,[],b)==expected
    assert x==a and [[type(v) for v in row] for row in x]==[[type(v) for v in row] for row in a]


def test_readonly_tuple_sources_and_zero_unknowns_have_complete_exact_routes(monkeypatch):
    calls=[];native=sx._core._ddm_ilu_solve
    def record(*args):
        result=native(*args);calls.append(result);return result
    monkeypatch.setattr(sx._core,'_ddm_ilu_solve',record)
    l=tuple(tuple(row) for row in qq([[1,0],[2,1]]))
    u=tuple(tuple(row) for row in qq([[2,1],[0,3]]))
    b=tuple(tuple(row) for row in qq([[1],[3]]));x=[[None],[None]]
    dense.ddm_ilu_solve(x,l,u,[],b)
    assert values(x)==[[Fraction(1,3)],[Fraction(1,3)]] and calls==[True]
    calls.clear();dense.ddm_ilu_solve([],qq([[1,0],[0,1]]),[[],[]],[],qq([[0],[0]]))
    assert calls==[True]


@pytest.mark.parametrize('length_call',range(4))
def test_actual_length_prefix_is_not_replayed_when_dimensions_are_modified(length_call):
    events=[];namespace=dict(dense.__dict__)
    def length(value):
        index=len(events);events.append(index)
        return len(value)-1 if index==length_call else len(value)
    namespace['len']=length
    current=types.FunctionType(dense.ddm_ilu_solve.__code__,namespace);old=old_function(namespace)
    l=qq([[1,0],[2,1]]);u=qq([[2,1],[0,3]]);b=qq([[1,2],[3,4]])
    def run(function):
        events.clear();x=[[None,None],[None,None]]
        try:result=function(x,l,u,[],b)
        except Exception as error:return type(error),str(error),values(x),events[:]
        return result,values(x),events[:]
    assert run(current)==run(old)


def test_self_restoring_range_is_consumed_only_in_the_actual_scratch_prefix(monkeypatch):
    events=[];namespace=dict(dense.__dict__)
    def hook(*args):
        events.append(args);namespace.pop('range');return range(*args)
    namespace['range']=hook
    current=types.FunctionType(dense.ddm_ilu_solve.__code__,namespace)
    x=[[None],[None]]
    current(x,qq([[1,0],[2,1]]),qq([[2,1],[0,3]]),[],qq([[1],[3]]))
    assert events==[(2,)] and values(x)==[[Fraction(1,3)],[Fraction(1,3)]]


@pytest.mark.parametrize('digits',[100,500,1500])
def test_large_exact_coefficients_and_multiple_rhs_do_not_cross_binary64(digits):
    large=10**digits;tiny=Fraction(1,large)
    raw_l=[[1,0],[Fraction(large,3),1]];raw_u=[[Fraction(large,7),tiny],[0,Fraction(13,17)]]
    known=[[tiny,Fraction(large,19)],[Fraction(23,29),-tiny]]
    y=[[sum(Fraction(raw_u[i][j])*known[j][k] for j in range(2)) for k in range(2)] for i in range(2)]
    rhs=[[sum(Fraction(raw_l[i][j])*y[j][k] for j in range(2)) for k in range(2)] for i in range(2)]
    x=[[None,None],[None,None]]
    assert dense.ddm_ilu_solve(x,qq(raw_l),qq(raw_u),[],qq(rhs)) is None
    assert values(x)==known


def test_non_string_namespace_key_does_not_add_early_user_equality_calls():
    events=[]
    class Key:
        def __hash__(self):return hash('reversed')
        def __eq__(self,other):events.append('key');raise RuntimeError('late lookup')
    namespace=dict(dense.__dict__);namespace[Key()]=None
    current=types.FunctionType(dense.ddm_ilu_solve.__code__,namespace);old=old_function(namespace)
    l=qq([[1,0],[2,1]]);u=qq([[2,1],[0,3]]);b=qq([[1],[3]])
    def run(function):
        events.clear();x=[[None],[None]]
        try:function(x,l,u,[],b)
        except RuntimeError as error:return str(error),values(x),events[:]
        raise AssertionError('expected equality hook')
    assert run(current)==run(old)
