"""Exact LP pivots, full feasible solutions and bounded-polytope certificates."""
from itertools import combinations
import random

import pytest
import sympy as sp
import symox as sx
from symox.solvers.simplex import linprog, lpmin, lpmax, InfeasibleLPError, UnboundedLPError


def foreign(value):
    if isinstance(value, (tuple, list)):
        return type(value)(foreign(v) for v in value)
    return sx.to_sympy(value)


@pytest.mark.parametrize('seed',range(24))
def test_bounded_rational_program_has_an_exact_vertex_certificate(seed):
    rng=random.Random(261000+seed)
    n=2+seed%2
    c=[rng.randrange(-4,5) for _ in range(n)]
    rows=[[int(i==j) for j in range(n)] for i in range(n)]
    rhs=[rng.randrange(2,6) for _ in range(n)]
    for _ in range(3):
        row=[rng.randrange(-3,4) for _ in range(n)]
        rows.append(row);rhs.append(sum(row)+rng.randrange(4))
    optimum,solution=foreign(linprog(c,rows,rhs))
    assert len(solution)==n and all(v>=0 for v in solution)
    assert all(sum(a*x for a,x in zip(row,solution))<=b for row,b in zip(rows,rhs))
    assert optimum==sum(a*x for a,x in zip(c,solution))
    # The explicit box bounds certify boundedness. Enumerate every vertex
    # exactly, independently of either library's simplex implementation.
    inequalities=rows+[[-int(i==j) for j in range(n)] for i in range(n)]
    bounds=rhs+[0]*n
    candidates=[]
    for indices in combinations(range(len(inequalities)),n):
        matrix=sp.Matrix([inequalities[i] for i in indices])
        if matrix.det()==0:continue
        vertex=matrix.inv()*sp.Matrix([bounds[i] for i in indices])
        if all(sum(a*x for a,x in zip(row,vertex))<=b for row,b in zip(inequalities,bounds)):
            candidates.append(sum(a*x for a,x in zip(c,vertex)))
    assert optimum==min(candidates)
    from sympy.solvers.simplex import linprog as reference
    assert optimum==reference(c,rows,rhs)[0]


@pytest.mark.parametrize('cost,rows,rhs,equalities,equality_rhs,bounds',[
    ([1,2],[[-1,-1]],[-3],[[1,-1],[2,-2]],[1,2],None),
    ([-1,-2],[[1,1]],[6],None,None,[(None,None),(1,3)]),
    ([1],[[0]],[1],[[1]],[sp.Rational(1,3)],None),
    ([-1],[[1]],[2],None,None,(None,1)),
    ([2,1],[[1,1]],[4],[[1,-1]],[0],{0:(1,2),1:(None,None)}),
])
def test_equalities_redundancy_and_variable_bounds_match_reference(cost,rows,rhs,equalities,equality_rhs,bounds):
    from sympy.solvers.simplex import linprog as reference
    args=(cost,rows,rhs,equalities,equality_rhs,bounds)
    # Inputs are Python rationals; each library owns its normalization.
    def convert(v):
        if isinstance(v,list):return [convert(a) for a in v]
        if isinstance(v,tuple):return tuple(convert(a) for a in v)
        if isinstance(v,dict):return {k:convert(a) for k,a in v.items()}
        if isinstance(v,sp.Rational):return sx.Rational(int(v.p),int(v.q))
        return v
    args_native=convert(args)
    actual=foreign(linprog(*args_native))
    expected=reference(*args)
    assert actual[0]==expected[0]


def test_infeasible_and_unbounded_remain_distinct_errors():
    with pytest.raises(InfeasibleLPError):linprog([1],[[1]],[-1])
    with pytest.raises(UnboundedLPError):linprog([-1])


def test_rational_public_program_bypasses_python_pivots(monkeypatch):
    import importlib
    module=importlib.import_module('symox.solvers.simplex')
    def forbidden(*args,**kwargs):
        raise AssertionError('the exact rational program reached Python pivots')
    monkeypatch.setattr(module,'_optimize',forbidden)
    assert linprog([1],[[-1]],[-2])==(2,[2])


def test_symbolic_frontends_keep_the_complete_optimum_assignment():
    x,y=sx.symbols('x y')
    value,assignment=lpmin(x+2*y,[x+y>=3,x-y>=1,x>=0,y>=0])
    assert value==3 and assignment=={x:sx.Integer(3),y:sx.Integer(0)}
    value,assignment=lpmax(x+y,[x+y<=4,x>=1,y>=1])
    assert value==4 and assignment[x]+assignment[y]==4


def test_bland_pivots_terminate_on_the_classic_cycling_program():
    q=sx.Rational
    cost=[-10,57,9,24]
    rows=[[q(1,2),q(-11,2),q(-5,2),9],
          [q(1,2),q(-3,2),q(-1,2),1],[1,0,0,0]]
    rhs=[0,0,1]
    optimum,solution=linprog(cost,rows,rhs)
    assert optimum==-1
    assert all(v>=0 for v in solution)
    assert all(sum(a*x for a,x in zip(row,solution))<=b for row,b in zip(rows,rhs))
    assert sum(c*x for c,x in zip(cost,solution))==optimum
    # Exact primal/dual equality certifies the optimum independently.
    dual=[0,-18,-1]
    assert all(sum(rows[i][j]*dual[i] for i in range(3))<=cost[j] for j in range(4))
    assert sum(b*y for b,y in zip(rhs,dual))==optimum


@pytest.mark.parametrize('value',[sx.Float('0.5'),sx.sqrt(2),sx.Symbol('unknown')])
def test_other_scalar_domains_decline_only_the_optimization(value):
    assert sx._core._native_simplex([value],[]) is None



def test_free_variable_recession_ray_proves_unboundedness():
    # x -> x-t preserves x+y<=6 and 1<=y<=3, reducing x-2*y.
    with pytest.raises(UnboundedLPError):
        linprog([1,-2],[[1,1]],[6],bounds=[(None,None),(1,3)])


def test_equality_without_inequalities_keeps_its_unique_exact_solution():
    assert linprog([1],A_eq=[[1]],b_eq=[sx.Rational(1,3)])==(sx.Rational(1,3),[sx.Rational(1,3)])
    # SymPy1.14's linprog currently fails constructing its block matrix for
    # this accepted native input; the unique equality certifies the result.
    from sympy.solvers.simplex import linprog as reference
    with pytest.raises(ValueError,match='mismatched dimensions'):
        reference([1],A_eq=[[1]],b_eq=[sp.Rational(1,3)])
