"""Exact bound transformations, native preparation and untouched host protocols."""
import importlib
from itertools import combinations
import random
import sys

import pytest
import sympy as sp
import symox as sx

lp=importlib.import_module('symox.solvers.simplex')


def native(c,A=None,b=None,A_eq=None,b_eq=None,bounds=None):
    return sx._core._native_linprog(c,A,b,A_eq,b_eq,bounds)


def foreign(value):
    if isinstance(value,(tuple,list)):
        return type(value)(foreign(v) for v in value)
    return sx.to_sympy(value)


@pytest.mark.parametrize('bounds',[
    None,(0,4),[(-3,2),(1,5)],{0:(-3,2),1:(1,5)},
    [(None,2),(-1,3)],[(-sx.oo,sx.oo),(1,3)],[(2,2),(-1,3)],
])
def test_bound_forms_keep_the_old_exact_assignment(bounds,monkeypatch):
    args=([-1,-2],[[1,1]],[6],None,None,bounds)
    result=native(*args)
    assert result is not None and result[0]=='ok'
    actual=lp.linprog(*args)
    monkeypatch.setattr(sx._core,'_native_linprog',lambda *args:None)
    expected=lp.linprog(*args)
    assert actual==expected


@pytest.mark.parametrize('seed',range(16))
def test_shifted_bound_program_has_independent_exact_vertex_certificate(seed,monkeypatch):
    rng=random.Random(88710+seed)
    c=[rng.randrange(-4,5),rng.randrange(-4,5)]
    bounds=[(rng.randrange(-4,0),rng.randrange(1,5)) for _ in c]
    rows=[[rng.randrange(-3,4) for _ in c] for _ in range(3)]
    rhs=[rng.randrange(0,5) for _ in rows]
    eq=[[1,1]] if seed%2 else None;eq_rhs=[0] if eq else None
    optimum,solution=foreign(lp.linprog(c,rows,rhs,eq,eq_rhs,bounds))
    inequalities=rows+[[1,0],[-1,0],[0,1],[0,-1]]
    limits=rhs+[bounds[0][1],-bounds[0][0],bounds[1][1],-bounds[1][0]]
    assert all(sum(a*x for a,x in zip(row,solution))<=b for row,b in zip(inequalities,limits))
    if eq:assert sum(solution)==0
    assert optimum==sum(a*x for a,x in zip(c,solution))
    values=[]
    equalities=eq or [];equal_limits=eq_rhs or []
    for indices in combinations(range(len(inequalities)),2-len(equalities)):
        matrix=sp.Matrix(equalities+[inequalities[i] for i in indices])
        if matrix.det()==0:continue
        vertex=matrix.inv()*sp.Matrix(equal_limits+[limits[i] for i in indices])
        if all(sum(a*x for a,x in zip(row,vertex))<=b for row,b in zip(inequalities,limits)):
            values.append(sum(a*x for a,x in zip(c,vertex)))
    assert optimum==min(values)
    # SymPy1.14 keeps an unintended nonnegative restriction for some negative
    # lower bounds. The complete vertex certificate above proves the actual
    # optimum; the existing general route must keep the same feasible point.
    monkeypatch.setattr(sx._core,'_native_linprog',lambda *args:None)
    assert (optimum,solution)==foreign(lp.linprog(c,rows,rhs,eq,eq_rhs,bounds))


def test_free_and_upper_only_variables_have_primal_dual_certificates():
    optimum,solution=lp.linprog([-1,-2],[[1,1]],[6],bounds=[(None,None),(1,3)])
    assert (optimum,solution)==(-9,[3,3])
    # x+y<=6 and y<=3 imply -x-2y>=-9, attained by the returned point.
    assert -solution[0]-2*solution[1]==-9
    assert lp.linprog([-3],bounds=(None,-2))==(6,[-2])
    with pytest.raises(lp.UnboundedLPError):lp.linprog([1],bounds=(None,-2))
    with pytest.raises(lp.InfeasibleLPError,match='lower bound exceeds upper bound'):
        lp.linprog([1],bounds=(2,1))


@pytest.mark.parametrize('args',[
    ([],None,None,None,None,None),([],[],[],[],[],(None,None)),
    ([],[[]],[0],None,None,None),([],None,None,[[]],[0],{}),
])
def test_zero_variable_programs_keep_their_mass_and_shape(args):
    assert lp.linprog(*args)==(0,[])


@pytest.mark.parametrize('A,b,A_eq,b_eq',[
    ([[]],[-1],None,None),(None,None,[[]],[1]),
])
def test_zero_variable_contradictions_are_infeasible(A,b,A_eq,b_eq):
    with pytest.raises(lp.InfeasibleLPError):lp.linprog([],A,b,A_eq,b_eq)


def test_exact_route_has_no_python_preparation_or_sympify_calls():
    codes={lp._matrix_rows.__code__,lp._bounds.__code__,lp._simplex.__code__,lp._optimize.__code__,sx.sympify.__code__}
    calls=[]
    def profile(frame,event,arg):
        if event=='call' and frame.f_code in codes:calls.append(frame.f_code.co_name)
    previous=sys.getprofile()
    try:
        sys.setprofile(profile)
        result=lp.linprog([1,2],[[-1,-1]],[-3],[[1,-1]],[1],bounds={0:(None,None),1:(0,4)})
    finally:sys.setprofile(previous)
    assert result==(4,[2,1])
    assert calls==[]


class Untouched:
    def __init__(self):self.calls=0
    def __iter__(self):self.calls+=1;raise AssertionError('custom iterator consumed')
    def _symox_(self):self.calls+=1;raise AssertionError('custom conversion consumed')


@pytest.mark.parametrize('position',range(6))
def test_generator_and_custom_protocols_decline_without_any_consumption(position):
    args=[[1],[[1]],[2],None,None,None]
    probe=Untouched();args[position]=probe
    assert native(*args) is None
    assert probe.calls==0
    source=(x for x in [1]);args[position]=source
    assert native(*args) is None
    assert next(source)==1


def test_nested_container_subclasses_are_not_inspected():
    class Row(list):
        def __iter__(self):raise AssertionError('subclass iteration')
    assert native([1],[Row([1])],[2]) is None
    class Coefficient(int):
        def __int__(self):raise AssertionError('int subclass conversion')
    assert native([Coefficient(1)],[[1]],[2]) is None
    class Rational(sx.Rational):
        def __new__(cls):return sx._core.Expr.__new__(cls,sx.Rational(1,2))
        def __mul__(self,other):raise KeyboardInterrupt
    assert native([Rational()],[[1]],[2]) is None
    with pytest.raises(KeyboardInterrupt):lp.linprog([Rational()],[[1]],[2])


def test_dictionary_key_protocols_remain_untouched_during_eligibility():
    class Key:
        calls=0
        def __hash__(self):return 0
        def __eq__(self,other):self.calls+=1;raise AssertionError('dictionary equality hook')
    key=Key();bounds={key:(0,1)}
    assert native([1],bounds=bounds) is None
    assert key.calls==0
    ignored=Untouched()
    assert lp.linprog([1],bounds={-1:ignored,3:ignored})==(0,[0])
    assert ignored.calls==0


@pytest.mark.parametrize('value',[sx.Float('0.5'),sx.sqrt(2),sx.Symbol('coefficient_host')])
def test_other_scalar_domains_keep_the_general_solver(value):
    assert native([value],[[1]],[2]) is None


@pytest.mark.parametrize('value',[sx.Float('0.5',60),sx.sqrt(2)])
def test_computable_nonrational_domains_keep_the_full_result(value,monkeypatch):
    args=([value],[[1]],[2],None,None,(1,2))
    actual=lp.linprog(*args)
    monkeypatch.setattr(sx._core,'_native_linprog',lambda *args:None)
    expected=lp.linprog(*args)
    assert actual==expected
    assert actual[0]==value and actual[1]==[1]
    assert getattr(actual[0],'_prec',None)==getattr(expected[0],'_prec',None)


def test_changed_preparation_and_scalar_hooks_still_run(monkeypatch):
    calls=[];original=lp._bounds
    def bounds(value,n):calls.append('bounds');return original(value,n)
    monkeypatch.setattr(lp,'_bounds',bounds)
    assert lp.linprog([1],bounds=(2,3))==(2,[2])
    assert calls==['bounds']
    def interrupted(*args,**kwargs):raise KeyboardInterrupt
    monkeypatch.setattr(sx._core,'_native_simplex',interrupted)
    with pytest.raises(KeyboardInterrupt):lp.linprog([1],[[1]],[2])


def test_registered_integer_converters_keep_the_general_sympification(monkeypatch):
    from symox.core.sympify import converter
    calls=[]
    def interrupted(value):calls.append(value);raise KeyboardInterrupt
    monkeypatch.setitem(converter,int,interrupted)
    with pytest.raises(KeyboardInterrupt):lp.linprog([1],[[1]],[2])
    assert calls==[1]


@pytest.mark.parametrize('name',['sympify','simplify'])
def test_changed_public_scalar_functions_keep_interruptions(name,monkeypatch):
    def interrupted(*args,**kwargs):raise KeyboardInterrupt
    monkeypatch.setattr(sx,name,interrupted)
    with pytest.raises(KeyboardInterrupt):lp.linprog([1],[[1]],[2])


@pytest.mark.parametrize('args',[
    ([1],[[1,2]],[2],None,None,None),([1],[[1]],[],None,None,None),
    ([1],None,None,[[1,2]],[1],None),([1],None,None,None,None,[]),
    ([1],None,None,None,None,[(0,1,2)]),([1],None,None,None,None,[(0,)]),
])
def test_malformed_exact_shapes_keep_original_error_type_and_message(args,monkeypatch):
    with pytest.raises(Exception) as actual:lp.linprog(*args)
    monkeypatch.setattr(sx._core,'_native_linprog',lambda *args:None)
    with pytest.raises(type(actual.value)) as expected:lp.linprog(*args)
    assert str(actual.value)==str(expected.value)


def test_shape_validation_precedes_later_rhs_protocols():
    probe=Untouched()
    with pytest.raises(ValueError,match='wrong number of columns'):
        lp.linprog([1],[[1,2]],probe)
    assert probe.calls==0


def test_arbitrary_size_integer_input_is_exact():
    n=2**200
    assert lp.linprog([n],A_eq=[[n]],b_eq=[3])==(3,[sx.Rational(3,n)])
