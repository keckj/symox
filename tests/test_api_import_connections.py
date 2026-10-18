"""Import contracts, initialization order and computations through public paths."""
import importlib
import json
from pathlib import Path
import subprocess
import sys
import types
import pytest
import symox as S
import sympy as R


@pytest.mark.parametrize('path',[
    'core', 'core.add', 'core.mod', 'core.power', 'core.relational',
    'functions', 'functions.special.elliptic_integrals', 'functions.special.zeta_functions',
    'integrals', 'integrals.transforms', 'calculus', 'sets', 'sets.powerset',
    'simplify', 'simplify.fu', 'simplify.combsimp', 'simplify.epathtools',
    'simplify.powsimp', 'simplify.ratsimp', 'ntheory.factor_', 'utilities',
    'utilities.enumerative', 'tensor', 'tensor.index_methods',
    'solvers.recurr', 'combinatorics.group_constructs', 'printing.pretty',
])
def test_complete_import_surfaces_against_reference(path):
    reference=importlib.import_module('sympy.'+path)
    native=importlib.import_module('symox.'+path)
    if hasattr(reference,'__all__'):names=reference.__all__
    else:
        names=[n for n in dir(reference) if not n.startswith('_')
               and callable(getattr(reference,n))
               and getattr(getattr(reference,n),'__module__','')==reference.__name__]
    assert names
    assert all(hasattr(native,name) for name in names)
    assert set(names)<=set(dir(native))


def test_function_exports_keep_late_installed_native_dispatch():
    functions=importlib.import_module('symox.functions')
    for name in importlib.import_module('sympy.functions').__all__:
        assert getattr(functions,name) is getattr(S,name)
    # Piecewise and the context-aware elementary functions are installed after
    # functions is first imported while symox initializes.
    x=S.Symbol('x')
    assert functions.Piecewise((1,x>0),(0,True))==S.Piecewise((1,x>0),(0,True))
    assert functions.conjugate(S.Matrix([[S.I]]))==S.Matrix([[-S.I]])
    with S.evaluate(False):assert functions.sqrt(4).is_Pow


@pytest.mark.parametrize('path,names',[
    ('matrices',['Matrix','MatrixSymbol','GramSchmidt','banded']),
    ('matrices.expressions',['MatrixSymbol','BlockMatrix','trace']),
    ('geometry.polygon',['Polygon','Triangle','RegularPolygon','rad','deg']),
    ('printing.pretty',['PrettyPrinter','pretty','pprint','pager_print']),
    ('core',['Symbol','Expr','sympify','common_assumptions']),
    ('polys',['Poly','groebner','QQ','ring']),
    ('tensor',['Array','tensorproduct','Indexed']),
])
def test_wildcard_import_preserves_existing_and_new_exports(path,names):
    namespace={}
    exec('from symox.'+path+' import *',namespace)
    assert all(name in namespace for name in names)
    assert all(not isinstance(namespace[n],types.ModuleType) for n in names)


def test_core_symbol_relational_and_assumption_imports():
    from symox.core.symbol import Symbol, symbols, Dummy
    from symox.core.relational import Eq, Ge
    from symox.core import assumptions,common_assumptions
    from symox.core.intfunc import igcd,integer_nthroot
    x=Symbol('x',positive=True)
    assert assumptions(x,['positive','negative'])=={'positive':True,'negative':False}
    assert assumptions(2,['integer','positive'])=={'integer':True,'positive':True}
    assert common_assumptions([0,1,2],['positive','integer'])=={'integer':True}
    assert common_assumptions([],['integer'])=={}
    assert Eq(x,x) is S.true or Eq(x,x)==S.true
    assert Ge(3,2)==S.true
    assert Dummy('x')!=Dummy('x')
    assert symbols('a b')==S.symbols('a b')
    assert igcd(36,24)==12 and integer_nthroot(81,4)==(3,True)
    for values in ([0,1,2],[-4,0,2],[S.Rational(1,2),2]):
        native=common_assumptions(values,['integer','positive','negative'])
        reference=R.core.common_assumptions([R.sympify(str(v)) for v in values],['integer','positive','negative'])
        assert native==reference


def test_transforms_and_integration_paths_compute_native_results():
    from symox.integrals.integrals import integrate, Integral
    from symox.integrals.transforms import laplace_transform
    from symox.integrals.laplace import inverse_laplace_transform
    from symox.integrals.singularityfunctions import singularityintegrate
    t=S.Symbol('t',positive=True);s=S.Symbol('s',positive=True)
    assert integrate(t**2,t)==t**3/3
    assert Integral(t,(t,0,1)).doit()==S.Rational(1,2)
    assert laplace_transform(S.exp(-2*t),t,s,noconds=True)==1/(s+2)
    assert inverse_laplace_transform(1/(s+2),s,t)==S.exp(-2*t)*S.Heaviside(t)
    x=S.Symbol('x')
    assert singularityintegrate(S.SingularityFunction(x,0,0),x)==S.SingularityFunction(x,0,1)


def test_simplification_rule_and_package_imports_preserve_callables():
    from symox.simplify import simplify,fu,cse
    from symox.simplify.fu import TR0,TR2,TR8
    from symox.simplify.radsimp import fraction
    x=S.Symbol('x')
    assert simplify is S.simplify
    assert TR0((x*x-1)/(x-1))==x+1
    assert TR2(S.tan(x))==S.sin(x)/S.cos(x)
    assert TR8(S.cos(x)**2)==S.cos(2*x)/2+S.Rational(1,2)
    assert fraction((x+1)/x)==(x+1,x)
    assert callable(fu) and callable(cse)
    assert S.simplify((x*x-1)/(x-1))==x+1


def test_domain_matrix_and_linear_solver_paths():
    from symox.polys import QQ,ring
    from symox.polys.matrices.domainmatrix import DM
    from symox.matrices import GramSchmidt,Matrix
    from symox.solvers.solveset import linear_eq_to_matrix
    from symox.solvers.recurr import rsolve
    x,y=S.symbols('x y')
    parent,rx=ring('x',QQ)
    assert (rx+1)**2==rx**2+2*rx+1
    assert DM([[1,2],[3,4]],QQ).to_Matrix()==Matrix([[1,2],[3,4]])
    a,b=linear_eq_to_matrix([x+y-3,x-y-1],[x,y]);assert a.LUsolve(b)==Matrix([2,1])
    vectors=GramSchmidt([Matrix([1,0]),Matrix([1,1])],orthonormal=True)
    assert vectors==[Matrix([1,0]),Matrix([0,1])]
    n=S.Symbol('n',integer=True);f=S.Function('f')
    solution=rsolve(f(n+1)-2*f(n),f(n),{f(0):3})
    assert S.simplify(solution-3*2**n)==0


def test_logic_and_indexed_paths_use_existing_algorithms():
    from symox.logic import And,Or,Not,satisfiable,to_cnf
    from symox.tensor.indexed import Idx,IndexedBase
    from symox.tensor.index_methods import get_indices
    p,q=S.symbols('p q')
    assert satisfiable(And(p,Not(p))) is False
    assert to_cnf(Or(p,And(p,q)))==And(p,Or(p,q))
    i=Idx('i',3);a=IndexedBase('a',shape=(3,))
    assert a[i]==S.IndexedBase('a',shape=(3,))[i]
    free,symmetry=get_indices(a[i]);assert i in free


def test_abc_symbols_and_enumerative_utilities():
    from symox.abc import x,E,I,beta,lamda
    from symox.utilities.enumerative import multiset_partitions_taocp,list_visitor
    assert x==S.Symbol('x') and E==S.Symbol('E') and I==S.Symbol('I')
    assert E!=S.E and I!=S.I and beta==S.Symbol('beta') and lamda==S.Symbol('lamda')
    assert len(list(multiset_partitions_taocp([2,1])))==4


_IMPORT_PATHS=[
    'functions','functions.special.elliptic_integrals','functions.special.error_functions',
    'functions.special.polynomials','functions.special.zeta_functions',
    'core','core.symbol','core.relational','core.assumptions','core.evalf','core.sympify',
    'integrals','integrals.integrals','integrals.transforms','integrals.laplace',
    'simplify','simplify.simplify','simplify.fu','simplify.hyperexpand','simplify.powsimp',
    'polys','polys.matrices.domainmatrix','matrices','matrices.immutable','matrices.sparse',
    'sets','sets.sets','sets.fancysets','solvers','solvers.diophantine.diophantine',
    'solvers.recurr','tensor.indexed','utilities','utilities.lambdify','utilities.enumerative',
    'interactive.printing','logic','logic.inference','printing.codeprinter','abc',
]
@pytest.mark.parametrize('reverse',[False,True])
def test_import_order_with_runtime_sympy_blocked(reverse):
    paths=_IMPORT_PATHS[::-1] if reverse else _IMPORT_PATHS
    code='''
import builtins,importlib
original=builtins.__import__
def blocked(name,*args,**kwargs):
    if name=='sympy' or name.startswith('sympy.'):
        raise ImportError('runtime SymPy delegation is forbidden')
    return original(name,*args,**kwargs)
builtins.__import__=blocked
for path in PATHS:importlib.import_module('symox.'+path)
import symox as S
assert S.core.sympify is S.sympify
assert S.utilities.lambdify is S.lambdify
assert S.solvers.diophantine is S.diophantine
assert importlib.import_module('symox.simplify').simplify is S.simplify
assert callable(S.series) and callable(S.simplify)
x=S.Symbol('x')
assert S.simplify((x*x-1)/(x-1))==x+1
assert S.functions.sqrt is S.sqrt
'''.replace('PATHS',repr(paths))
    result=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode==0,result.stderr
