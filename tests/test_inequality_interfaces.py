"""Polynomial sign charts, rational poles and nested absolute-value branches."""
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx


RELATIONS=['==','!=','<','<=','>','>=']


@pytest.mark.parametrize('relation',RELATIONS)
@pytest.mark.parametrize('source',[
    '0','1','-2','x','x*x-1','(x*x-1)**2',
    '(x*x-2)*(x*x-1)','-(x-1)**2*(x+2)',
    'x*x+1','(x-1)*(x+2)*(x-3)',
])
def test_polynomial_interval_components_match_reference(source,relation):
    x=sx.Symbol('x');X=sp.Symbol('x')
    expression=sx.sympify(eval(source,{'x':x}));foreign=sp.sympify(eval(source,{'x':X}))
    result=sx.solve_poly_inequality(sx.Poly(expression,x),relation)
    expected=sp.solve_poly_inequality(sp.Poly(foreign,X),relation)
    assert [part.to_sympy() for part in result]==expected


@pytest.mark.parametrize('relation',RELATIONS)
@pytest.mark.parametrize('numerator,denominator',[
    ('x-1','x+1'),('x*x-1','x-1'),('x*x+1','(x-1)**2'),
    ('0','x'),('1','0'),('x*x-2','x*x-1'),
])
def test_rational_systems_preserve_denominator_poles(numerator,denominator,relation):
    x=sx.Symbol('x');X=sp.Symbol('x')
    def problem(module,variable):
        p=module.Poly(module.sympify(eval(numerator,{'x':variable})),variable)
        q=module.Poly(module.sympify(eval(denominator,{'x':variable})),variable)
        return [[((p,q),relation)]]
    result=sx.solve_rational_inequalities(problem(sx,x))
    expected=sp.solve_rational_inequalities(problem(sp,X))
    assert result.to_sympy()==expected


def test_rational_conjunctions_disjunctions_and_empty_groups():
    def problem(module):
        x=module.Symbol('x');P=lambda expr:module.Poly(expr,x)
        return [[],[((P(x*x-1),P(x-1)),'>='),((P(x-2),P(1)),'<=')],
                [((P(x+4),P(1)),'<')]]
    assert sx.solve_rational_inequalities(problem(sx)).to_sympy()==sp.solve_rational_inequalities(problem(sp))
    assert sx.solve_rational_inequalities([])==sx.S.EmptySet
    assert sx.solve_rational_inequalities([[]])==sx.S.EmptySet


@pytest.mark.parametrize('relation',RELATIONS)
@pytest.mark.parametrize('source',[
    'Abs(x-5)-3',
    'Abs(Abs(x)-2)-1',
    'Abs(x-4)+Abs(3*x-5)-7',
    'Abs(2*x+1)*Abs(x-2)-3',
    'Abs(x)**2-4',
    '1/Abs(x)-1',
    'Abs(Abs(2*x-1)-3)+Abs(x+2)-5',
])
def test_absolute_value_inequalities_have_the_same_real_solution_set(source,relation):
    from symox.solvers.inequalities import _abs_domain
    x=sx.Symbol('x',real=True);X=sp.Symbol('x',real=True)
    expression=eval(source,dict(vars(sx),x=x));foreign=eval(source,dict(vars(sp),x=X))
    result=sx.reduce_abs_inequality(expression,relation,x)
    domain=_abs_domain(expression,relation,x)
    expected=sp.reduce_abs_inequality(foreign,relation,X)
    assert domain.to_sympy()==expected.as_set()
    assert result.to_sympy().as_set()==expected.as_set()
    condition=sx.Relational(expression,0,relation)
    assert sx.solve_univariate_inequality(condition,x,relational=False).to_sympy()==expected.as_set()
    assert sx.reduce_inequalities([condition],x).to_sympy().as_set().intersect(sp.S.Reals)==expected.as_set()


def test_absolute_value_system_and_module_exports():
    from symox.solvers.inequalities import solve_poly_inequalities
    x=sx.Symbol('x',real=True);X=sp.Symbol('x',real=True)
    exprs=[(sx.Abs(3*x-5)-7,'<'),(sx.Abs(x+25)-13,'>')]
    foreign=[(sp.Abs(3*X-5)-7,'<'),(sp.Abs(X+25)-13,'>')]
    result=sx.reduce_abs_inequalities(iter(exprs),x)
    assert result.to_sympy().as_set()==sp.reduce_abs_inequalities(foreign,X).as_set()
    assert solve_poly_inequalities([(sx.Poly(x-1,x),'<'),(sx.Poly(x-2,x),'>')]).to_sympy()==sp.Union(sp.Interval.open(-sp.oo,1),sp.Interval.open(2,sp.oo))


@pytest.mark.parametrize('source',[
    '[[x*x<=0]]','[[x+2>0]]','[[(x+2,">")]]','[[x+2]]',
    '[[x*x>1,x<3],[x<-2]]','[[True]]','[[False]]','[[],[x>=1]]',
    '[[1/(x-1)>0,x<=3]]','[[Float(".5")*x+1>0]]',
])
@pytest.mark.parametrize('real',[False,True])
def test_rational_reduction_accepts_relations_expression_pairs_and_equations(source,real):
    from symox.solvers.inequalities import reduce_rational_inequalities
    from sympy.solvers.inequalities import reduce_rational_inequalities as reference
    x=sx.Symbol('x',**({'real':True} if real else {}))
    X=sp.Symbol('x',**({'real':True} if real else {}))
    exprs=eval(source,dict(vars(sx),x=x));foreign=eval(source,dict(vars(sp),x=X))
    foreign=[[sp.sympify(item) if isinstance(item,bool) else item for item in group] for group in foreign]
    actual=reduce_rational_inequalities(exprs,x,relational=False)
    expected=reference(foreign,X,relational=False)
    assert actual.to_sympy()==expected
    assert reduce_rational_inequalities(exprs,x).to_sympy().as_set().intersect(sp.S.Reals)==expected


def test_rational_reduction_known_sign_parameters_and_relational_metadata():
    from symox.solvers.inequalities import reduce_rational_inequalities
    x=sx.Symbol('x',real=True);a=sx.Symbol('a',positive=True)
    result=reduce_rational_inequalities([[a*x>1]],x,relational=False)
    assert result==sx.Interval.open(1/a,sx.oo)
    for coefficient in [sx.sqrt(2),sx.pi,sx.E]:
        assert reduce_rational_inequalities([[coefficient*x>1]],x,relational=False)==sx.Interval.open(1/coefficient,sx.oo)
    for relation,condition in zip(RELATIONS,[sx.Eq(x,0),sx.Ne(x,0),x<0,x<=0,x>0,x>=0]):
        assert condition.rel_op==relation
    with pytest.raises(AttributeError):x.rel_op
    with pytest.raises(TypeError):reduce_rational_inequalities([[x>0]],x+1)


def test_nested_absolute_values_respect_supplied_domains_and_system_constraints():
    x=sx.Symbol('x',real=True)
    condition=sx.Abs(sx.Abs(x)-2)<1
    assert sx.solve_univariate_inequality(condition,x,relational=False,domain=sx.Interval(0,2))==sx.Interval.Lopen(1,2)
    result=sx.reduce_inequalities(iter([condition,x>0]),iter([x]))
    assert result==sx.And(sx.Lt(1,x),sx.Lt(x,3))
    assert sx.reduce_inequalities([condition])==sx.reduce_abs_inequality(sx.Abs(sx.Abs(x)-2)-1,'<',x)


def test_invalid_variables_domains_relations_and_branch_limits():
    x,y=sx.symbols('x y')
    with pytest.raises(ValueError):sx.solve_poly_inequality(x,'>')
    with pytest.raises(ValueError):sx.solve_poly_inequality(sx.Poly(x,x),'bad')
    with pytest.raises(ValueError):sx.solve_poly_inequality(sx.Poly(x+y,x,y),'>')
    with pytest.raises(NotImplementedError):sx.solve_poly_inequality(sx.Poly(sx.sqrt(2)*x,x,extension=sx.sqrt(2)),'>')
    with pytest.raises(NotImplementedError):sx.solve_poly_inequality(sx.Poly(x,x,modulus=3),'>')
    with pytest.raises(ValueError):sx.solve_rational_inequalities([[((sx.Poly(x,x),sx.Poly(y,y)),'>')]])
    with pytest.raises(ValueError):sx.reduce_abs_inequality(sx.sqrt(x)-1,'<',x)
    with pytest.raises(TypeError):sx.reduce_abs_inequality(sx.Abs(x)-1,'<',x+1)
    z=sx.Symbol('z',imaginary=True)
    with pytest.raises(TypeError):sx.reduce_abs_inequality(sx.Abs(z)-1,'<',z)
    expression=sum(sx.Abs(x-k) for k in range(11))
    with pytest.raises(NotImplementedError):sx.reduce_abs_inequality(expression-10,'<',x)


def test_inequality_interfaces_do_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
x=Symbol('x',real=True)
assert solve_poly_inequality(Poly(x*x-1,x),'<')==[Interval.open(-1,1)]
assert solve_rational_inequalities([[((Poly(x-1,x),Poly(x+1,x)),'<')]])==Interval.open(-1,1)
assert reduce_abs_inequality(Abs(Abs(x)-2)-1,'<',x)==Or(And(Lt(-3,x),Lt(x,-1)),And(Lt(1,x),Lt(x,3)))
assert solve_univariate_inequality(Abs(Abs(x)-2)<1,x,relational=False)==Union(Interval.open(-3,-1),Interval.open(1,3))
assert reduce_inequalities([Abs(Abs(x)-2)<1,x>0],x)==And(Lt(1,x),Lt(x,3))
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=30)
