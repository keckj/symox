"""Matrix extraction and residual certificates for native ODE systems."""
import pytest
import sympy as sp
import symox as sx
from symox.solvers.ode.systems import (
    linear_ode_to_matrix, canonical_odes, linodesolve, ODEOrderError, ODENonlinearError,
)


@pytest.mark.parametrize('sources,order', [
    (['diff(f,t)-f-g-1','diff(g,t)-f+g'],1),
    (['diff(f,t,2)+2*diff(g,t)-f','diff(g,t,2)+g'],2),
    (['diff(f,t)/t+f-g','2*diff(g,t)/(t+1)-f'],1),
    (['2*diff(f,t)+3*diff(g,t)-f-t','diff(f,t)-diff(g,t)-g'],1),
    (['f+2*g-1','2*f-g+3'],0),
])
def test_linear_system_matrices_match_reference_and_reconstruct(sources, order):
    from sympy.solvers.ode.systems import linear_ode_to_matrix as reference
    results = []
    for module, converter in [(sx,linear_ode_to_matrix),(sp,reference)]:
        t = module.Symbol('t')
        f,g = module.Function('f')(t),module.Function('g')(t)
        equations = [eval(source,vars(module)|{'t':t,'f':f,'g':g}) for source in sources]
        matrices,forcing = converter(equations,[f,g],t,order)
        coordinates = module.Matrix([f,g])
        reconstructed = matrices[0]*coordinates.diff(t,order)-forcing
        for i,matrix in enumerate(matrices[1:]):
            reconstructed -= matrix*coordinates.diff(t,order-i-1)
        assert all(module.simplify(value)==0 for value in reconstructed-module.Matrix(equations))
        results.append((matrices,forcing))
    (actual,b),(expected,B) = results
    assert all(a.to_sympy()==e for a,e in zip(actual,expected))
    assert b.to_sympy()==B


def test_matrix_extraction_rejects_nonlinearity_and_wrong_order():
    t = sx.Symbol('t')
    f,g = sx.Function('f')(t),sx.Function('g')(t)
    for equations in [[f.diff(t)-g*g,g.diff(t)-f], [1/f.diff(t)-g,g.diff(t)-f]]:
        with pytest.raises(ODENonlinearError):
            linear_ode_to_matrix(equations,[f,g],t,1)
    with pytest.raises(ODEOrderError):
        linear_ode_to_matrix([f.diff(t,2),g.diff(t)],[f,g],t,1)


@pytest.mark.parametrize('sources', [
    ['diff(f,t)-f-g-1','diff(g,t)-f+g'],
    ['2*diff(f,t)+diff(g,t)-f','diff(g,t)-g'],
    ['diff(f,t)**2-1','diff(g,t)-f'],
])
def test_canonical_highest_derivative_isolation(sources):
    t = sx.Symbol('t')
    f,g = sx.Function('f')(t),sx.Function('g')(t)
    equations = [eval(source,vars(sx)|{'t':t,'f':f,'g':g}) for source in sources]
    branches = canonical_odes(equations,[f,g],t)
    assert branches
    for branch in branches:
        substitutions = {item.lhs:item.rhs for item in branch}
        assert all(sx.simplify(item.xreplace(substitutions))==0 for item in equations)


@pytest.mark.parametrize('matrix,forcing', [
    ([[1,1],[0,1]],None),
    ([[0,1],[-1,0]],None),
    ([[1,0],[0,-2]],['t','exp(t)']),
    ([[0,1],[0,0]],['1','t*t']),
    ([['t',0],[0,'2*t']],None),
    ([['t',0],[0,'2*t']],['t','0']),
    ([[1,1,0],[0,1,1],[0,0,1]],None),
])
@pytest.mark.parametrize('doit',[False,True])
def test_linear_matrix_solutions_have_zero_residual(matrix, forcing, doit):
    t = sx.Symbol('t')
    def convert(value):
        return eval(value,vars(sx)|{'t':t}) if isinstance(value,str) else value
    A = sx.Matrix([[convert(value) for value in row] for row in matrix])
    b = sx.zeros(A.rows,1) if forcing is None else sx.Matrix([convert(value) for value in forcing])
    solution = sx.Matrix(linodesolve(A,t,b,doit=doit))
    residual = solution.diff(t)-A*solution-b
    assert all(sx.simplify(value)==0 for value in residual)
    # The homogeneous solution spans the full state space.
    constants = sorted(solution.free_symbols-A.free_symbols-b.free_symbols-{t},key=sx.default_sort_key)
    assert len(constants)==A.rows
    fundamental = solution.jacobian(constants)
    assert sx.simplify(fundamental.det()) != 0


def test_noncommuting_time_coefficients_and_invalid_antiderivatives_are_rejected():
    t = sx.Symbol('t')
    with pytest.raises(NotImplementedError):
        linodesolve(sx.Matrix([[0,t],[1,0]]),t)
    # SymPy does not check that a supplied B is an antiderivative of A.
    assert len(linodesolve(sx.eye(2),t,B=sx.eye(2)))==2
    with pytest.raises(ValueError):
        linodesolve(sx.eye(2),t,b=sx.Matrix([[1,2]]))


@pytest.mark.parametrize('forcing',[False,True])
def test_public_dsolve_accepts_linear_systems_and_initial_conditions(forcing):
    t,C1 = sx.symbols('t C1')
    f,g = sx.Function('f')(t),sx.Function('g')(t)
    equations = [f.diff(t)-f-g-(C1 if forcing else 0),g.diff(t)-g]
    solutions = sx.dsolve(equations,[f,g],ics={f.subs(t,0):1,g.subs(t,0):2})
    values = {item.lhs:item.rhs for item in solutions}
    for equation in equations:
        residual = equation.xreplace({f.diff(t):values[f].diff(t),g.diff(t):values[g].diff(t),**values})
        assert sx.simplify(residual)==0
    assert values[f].subs(t,0)==1
    assert values[g].subs(t,0)==2
    if forcing:
        assert C1 in values[f].free_symbols


def test_inferred_system_functions_match_explicit_function_order():
    t = sx.Symbol('t')
    f,g = sx.Function('f')(t),sx.Function('g')(t)
    equations = [f.diff(t)-f-g,g.diff(t)-g]
    assert sx.dsolve(equations) == sx.dsolve(equations,[f,g])


def test_initial_condition_symbols_cannot_become_generated_constants():
    t,C1,C2 = sx.symbols('t C1 C2')
    f,g = sx.Function('f')(t),sx.Function('g')(t)
    solutions = sx.dsolve([f.diff(t)-f,g.diff(t)-g],[f,g],
                         ics={f.subs(t,0):C1,g.subs(t,0):C2})
    assert solutions == [sx.Eq(f,C1*sx.exp(t)),sx.Eq(g,C2*sx.exp(t))]
