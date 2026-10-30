"""Equation helper results have denominator, determinant and root certificates."""
import itertools
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
import sympy.solvers.solvers as reference
from symox.solvers import solvers


@pytest.mark.parametrize('source',[
    'x/y','x/(y*z)','3/x+y/z','x/2+y/z','1/(x+y)+1/(x-y)',
    'sin(1/x)+cos(1/(y+1))','Eq(x/y,1/z)','1/sqrt(x)',
    '[1/x,1/y]','[Eq(x/y,1/z),1/(x+z)]',
])
@pytest.mark.parametrize('selected',[(),('x',),('y','z')])
def test_denominators_match_reference_recursively(source,selected):
    native=vars(sx)|dict(zip(['x','y','z'],sx.symbols('x y z')))
    foreign=vars(sp)|dict(zip(['x','y','z'],sp.symbols('x y z')))
    actual=solvers.denoms(eval(source,native),*[native[name] for name in selected])
    expected=reference.denoms(eval(source,foreign),*[foreign[name] for name in selected])
    assert {value.to_sympy() for value in actual}==expected


def test_recasting_preserves_order_overlaps_and_round_trips():
    x,y=sx.symbols('x y');f=sx.Function('f')(x)
    equations=[f+x*x+y,sx.Eq(f,x*x),x,y]
    unknowns=[f,y,x*x,f]
    changed,symbols,restore=solvers.recast_to_symbols(equations,unknowns)
    assert all(symbol.is_Symbol for symbol in symbols)
    assert symbols[0]==symbols[-1] and symbols[1]==y
    assert [value.xreplace(restore) for value in changed]==equations
    assert [restore.get(value,value) for value in symbols]==unknowns
    with pytest.raises(ValueError):solvers.recast_to_symbols(x,[x])


@pytest.mark.parametrize('size',range(6))
@pytest.mark.parametrize('method',[solvers.det_perm,solvers.det_minor,solvers.det_quick])
def test_determinants_agree_with_independent_bareiss_and_do_not_mutate(size,method):
    x=sx.Symbol('x')
    values=[(i+1)*(j+2)+(x if i==j else 0) for i in range(size) for j in range(size)]
    matrix=sx.Matrix(size,size,values)
    expected=matrix.to_sympy().det(method='bareiss')
    actual=method(matrix)
    assert sp.cancel(actual.to_sympy()-expected)==0
    assert list(matrix)==values


@pytest.mark.parametrize('values',[
    [[2]],[[1,2],[3,5]],[[1,2,0],[0,1,3],[2,0,1]],
    [['x',1],[1,'x']],[[1,'x',0],[0,1,'y'],[0,0,1]],
])
def test_quick_inverse_has_left_and_right_identity_certificates(values):
    matrix=sx.Matrix(values)
    actual=solvers.inv_quick(matrix)
    for product in [actual*matrix,matrix*actual]:
        assert all(sx.cancel(product[i,j]-(1 if i==j else 0))==0
                   for i in range(matrix.rows) for j in range(matrix.cols))
    expected=matrix.to_sympy().inv()
    assert all(sp.cancel(actual[i,j].to_sympy()-expected[i,j])==0
               for i in range(matrix.rows) for j in range(matrix.cols))
    assert sx.det_quick is solvers.det_quick and sx.inv_quick is solvers.inv_quick


def test_matrix_helper_validation_and_native_fallbacks():
    matrix=sx.Matrix([[1,2,3]])
    for method in [solvers.det_perm,solvers.det_minor,solvers.det_quick,solvers.inv_quick]:
        with pytest.raises(sx.NonSquareMatrixError):method(matrix)
    with pytest.raises(sx.NonInvertibleMatrixError):solvers.inv_quick(sx.Matrix([[1,2],[2,4]]))
    assert solvers.det_perm(sx.eye(9))==1
    assert solvers.det_minor(sx.eye(13))==1
    assert solvers.det_quick(sx.eye(3),method='bareiss')==1
    assert solvers.inv_quick(sx.eye(7))==sx.eye(7)


@pytest.mark.parametrize('rows',[
    [[1,1,1,2]],[[1,1,0,1],[0,1,1,1]],
    [[1,2,3,0]],[[1,-1,0,0],[0,1,-1,0]],
    [[0,0,0,1]],[[0,0,0,0]],[[1,0,0,0]],
    [[1,1,1,1],[1,2,3,2]],
])
@pytest.mark.parametrize('quick',[False,True])
def test_sparse_linear_solutions_have_valid_residuals_and_minimal_support(rows,quick):
    variables=sx.symbols('x y z');matrix=sx.Matrix(rows)
    result=solvers.minsolve_linear_system(matrix,*variables,quick=quick)
    standard=sx.solve_linear_system(matrix,*variables)
    if standard is None:
        assert result is None
        return
    if not standard or all(value==0 for value in standard.values()):
        assert result==standard
        return
    assert set(result)==set(variables)
    assert matrix[:,:-1]*sx.Matrix([result[value] for value in variables])==matrix[:,-1:]
    if not quick:
        support=sum(value!=0 for value in result.values())
        homogeneous=all(row[-1]==0 for row in rows)
        # Independently rule out every smaller support using matrix ranks.
        for size in range(1,support):
            for subset in itertools.combinations(range(3),size):
                coefficients=matrix[:,list(subset)]
                if homogeneous:
                    assert coefficients.rank()==size
                else:
                    assert coefficients.row_join(matrix[:,-1:]).rank()>coefficients.rank()


@pytest.mark.parametrize('source',[
    'sqrt(x)*x**Rational(1,3)+2','sqrt(x)+root(x+1,3)',
    'sqrt(x+1)+sqrt(x-1)-2','sqrt(1+sqrt(x))-2','1/sqrt(x)-2',
    'sqrt(x)+sqrt(x+1)+sqrt(x+2)-3',
    'Eq(sqrt(x+1),x-1)','sqrt(x)+y','sqrt(x+1/x)-2',
])
def test_radical_elimination_matches_reference_polynomial_root_equations(source):
    x,y=sx.symbols('x y');X,Y=sp.symbols('x y')
    expression=eval(source,vars(sx)|{'x':x,'y':y})
    foreign=eval(source,vars(sp)|{'x':X,'y':Y})
    actual,cov=solvers.unrad(expression)
    expected,reference_cov=reference.unrad(foreign)
    assert cov==reference_cov==[]
    # Elimination is defined up to a nonzero factor and repeated factors.
    a=sp.Poly(actual.to_sympy(),X,Y).sqf_part().monic()
    b=sp.Poly(expected,X,Y).sqf_part().monic()
    assert a==b
    assert not any(node.is_Pow and node.exp.is_Rational and node.exp.q!=1
                   and node.free_symbols for node in sp.preorder_traversal(actual.to_sympy()))


def test_shared_radical_coordinate_retains_its_covariant_relation():
    x=sx.Symbol('x')
    expression=sx.sqrt(x)+sx.root(x,3)-2
    polynomial,(p,relation)=solvers.unrad(expression)
    assert polynomial==p**3+p**2-2
    assert relation==p**6-x
    assert sx.cancel(expression.subs(x,1))==0
    assert polynomial.subs(p,1)==0 and relation.subs({p:1,x:1})==0
    assert solvers.unrad(x*x+1) is None
    assert solvers.unrad(sx.sqrt(x)+1,sx.Symbol('y')) is None
    # SymPy 1.14 values
    assert solvers.unrad(x**sx.Rational(1,65)+1)==(x+1,[])
    assert solvers.unrad(sx.sin(sx.sqrt(x))) is None


def test_four_distinct_radicals_eliminate_even_when_reference_search_stops():
    x=sx.Symbol('x')
    expression=sx.sqrt(x)+sx.sqrt(x+1)+sx.sqrt(x+4)+sx.sqrt(x+9)-6
    # The public unrad is SymPy's, which stops here; solve and solveset use
    # the established resultant elimination and keep the root.
    import pytest
    from symox.solvers._algebraic import _unrad_resultant
    with pytest.raises(NotImplementedError, match='Cannot remove all radicals'):
        solvers.unrad(expression)
    assert sx.solve(expression,x)==[0]
    assert sx.solveset(expression,x)==sx.FiniteSet(0)
    polynomial,cov=_unrad_resultant(expression)
    assert not cov and sx.Poly(polynomial,x).degree()>0
    assert expression.subs(x,0)==0 and polynomial.subs(x,0)==0
    assert all(not(node.is_Pow and node.exp.is_Rational and node.exp.q!=1)
               for node in sx.preorder_traversal(polynomial))


def test_solver_helpers_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
from symox.solvers.solvers import *
x,y=symbols('x y')
assert denoms(x/y)=={y}
assert det_quick(Matrix([[x,1],[1,x]]))==x*x-1
assert unrad(sqrt(1+sqrt(x))-2)==(x-9,[])
assert minsolve_linear_system(Matrix([[1,1,2]]),x,y)=={x:0,y:2}
'''
    subprocess.run([sys.executable,'-c',source],check=True)
