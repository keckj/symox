"""Behavioral regressions for the final 88 root exports (no SymPy runtime)."""
import contextvars
import importlib
import sys
import types
from fractions import Fraction

import mpmath as mp
import pytest
import symox as s

x,y,z = s.symbols('x y z')


def equivalent(actual,expected):
    assert s.simplify(s.expand(actual-expected)) == 0


def test_root_surface_is_complete_in_fresh_interpreter():
    import subprocess
    result = subprocess.run([sys.executable,'-c',
        'import symox, sympy; assert not [n for n in sympy.__all__ if not hasattr(symox,n)]'],capture_output=True,text=True)
    assert result.returncode == 0,result.stderr


def test_evaluation_context_restores_state_and_overrides():
    with s.evaluate(False):
        assert s.Add(1,1).args == (s.Integer(1),s.Integer(1))
        assert (s.Integer(1)+1).is_Add
        assert s.sin(0).is_Function
        assert s.sqrt(4).is_Pow
        assert s.ln(1).class_name == "log"
        assert s.log(x,2).args == (x,s.Integer(2))
        assert s.Eq(x,x).is_Relational
        assert s.Mul(2,3,evaluate=True) == 6
        assert contextvars.copy_context().run(lambda:s.Add(1,1).is_Add)
        with s.evaluate(True):assert s.Add(1,1) == 2
        assert s.Add(1,1).is_Add
    assert s.Add(1,1) == 2
    with pytest.raises(RuntimeError):
        with s.evaluate(False):raise RuntimeError
    assert s.Add(1,1) == 2


def test_division_error_and_core_protocols():
    try:
        s.seterr(divide=True)
        with pytest.raises(ValueError):s.Integer(0)/0
    finally:s.seterr(divide=False)
    assert s.Integer(0)/0 == s.nan
    with pytest.raises(s.SympifyError):s.sympify(object())
    d = s.Dict((x,1),(y,x+1))
    assert dict(d) == {x:s.Integer(1),y:x+1}
    assert d.subs(x,2)[y] == 3
    assert d.free_symbols == {x,y}
    assert hash(d) == hash(s.Dict({y:x+1,x:1}))
    assert isinstance(s.Function('f'),s.FunctionClass)
    assert isinstance(s.sin,s.FunctionClass)
    assert isinstance(s.log,s.FunctionClass)
    F = s.WildFunction('F',nargs=1)
    assert s.sin(x).match(F) == {F:s.sin(x)}
    assert (x+1).match(F) is None
    assert s.Function('g')(x,y).match(F) is None
    assert s.RealNumber('1.5') == s.Float('1.5')
    result = s.evalf(s.pi,100,{})
    with mp.workprec(100):assert abs(mp.mpf(result[0])-mp.pi) < mp.mpf(2)**-95
    with pytest.raises(s.PrecisionExhausted):s.evalf(s.pi,100,{'maxprec':50})


@pytest.mark.parametrize('family',['FF_gmpy','FF_python','GMPYFiniteField','PythonFiniteField'])
def test_finite_domain_aliases(family):
    domain = getattr(s,family)(7)
    assert domain(8) == domain(1)


@pytest.mark.parametrize('family',['ZZ_gmpy','ZZ_python'])
def test_integer_domains(family):
    assert getattr(s,family)()(12) == s.ZZ(12)


@pytest.mark.parametrize('family',['QQ_gmpy','QQ_python'])
def test_rational_domains(family):
    assert getattr(s,family)()(1,3) == s.QQ(1,3)


def test_domain_constructors_and_options():
    assert s.GMPYIntegerRing() == s.ZZ
    assert s.PythonIntegerRing() == s.ZZ
    assert s.GMPYRationalField() == s.QQ
    assert s.PythonRational(1,3)+s.PythonRational(1,6) == Fraction(1,2)
    options = s.Options((x,),{'domain':'QQ','order':'lex'})
    assert options.domain == s.QQ and options.gens == (x,)
    assert options.clone({'expand':False}).expand is False
    with pytest.raises(ValueError):s.Options((x,x),{})
    with pytest.raises(ValueError):s.Options((),{'bad':True})


def test_root_order_precision_and_symmetric_sums():
    assert s.RootOf(x*x-2,0) == -s.sqrt(2)
    assert s.rootof(x*x-2,-1) == s.sqrt(2)
    root = s.CRootOf(x**5-x-1,0)
    assert root.free_symbols == set() and root.index == 0 and root.is_real is True
    assert root.is_number and root.is_algebraic is True
    assert root.poly.as_expr() == x**5-x-1
    with mp.workdps(65):
        actual = mp.mpf(str(root.evalf(60)))
        assert abs(actual**5-actual-1) < mp.mpf('1e-58')
    roots = [s.CRootOf(x*x+1,i).evalf(40) for i in range(2)]
    assert abs(complex(roots[0])+1j) < 1e-30
    assert abs(complex(roots[1])-1j) < 1e-30
    assert s.CRootOf((x-2)**3,1) == 2
    assert s.RootSum(x**3-x-1,s.Lambda(x,x*x)) == 2
    equivalent(s.RootSum(x*x-2,s.Lambda(x,1/(y-x))),2*y/(y*y-2))
    held = s.RootSum(x*x-2,s.Lambda(x,s.exp(x)),auto=False)
    equivalent(held.doit(),s.exp(-s.sqrt(2))+s.exp(s.sqrt(2)))
    with pytest.raises(IndexError):s.CRootOf(x*x-1,2)
    with pytest.raises(NotImplementedError):s.CRootOf(x*x-y,x,0)


def test_matrix_exports_and_deferred_lambdify():
    M = s.Matrix([[1,2],[3,4]])
    assert s.per(M) == 10
    assert s.Permanent(M).doit() == 10
    A = s.MatrixSymbol('A',4,4)
    assert s.Permanent(A).doit() == s.Permanent(A)
    blocks = s.blockcut(A,(1,3),(2,2))
    assert blocks.shape == (4,4)
    assert s.matrix_symbols(blocks) == {A}
    assert s.diagonalize_vector(s.Matrix([1,2])) == s.diag(1,2)
    X = s.DeferredVector('X')
    fn = s.lambdify(X,X[0]**2+X[2])
    assert fn([3,99,4]) == 13
    with pytest.raises(IndexError):X[-1]


def test_indexed_free_and_contracted_indices():
    i,j,k = [s.Idx(name,3) for name in ('i','j','k')]
    A,B = s.IndexedBase('A'),s.IndexedBase('B')
    assert s.get_indices(A[i,j]) == ({i,j},{})
    expr = A[i,j]*B[j,k]
    assert s.get_indices(expr) == ({i,k},{})
    assert s.get_contraction_structure(expr)[(j,)] == {expr}
    assert s.get_indices(A[i,i]) == (set(),{})
    with pytest.raises(ValueError):s.get_indices(A[i]+B[j])


@pytest.mark.parametrize('function,args,expected',[
    ('chebyshevt_root',(3,1),s.Integer(0)),
    ('chebyshevu_root',(3,1),s.Integer(0)),
    ('hermite_prob',(4,x),x**4-6*x*x+3),
    ('betainc',(2,3,0,x),x*x/2-2*x**3/3+x**4/4),
    ('erf2inv',(x,0),x),
    ('riemann_xi',(0,),s.Rational(1,2)),
    ('riemann_xi',(1,),s.Rational(1,2)),
])
def test_special_exact_identities(function,args,expected):
    equivalent(getattr(s,function)(*args),expected)


def test_special_derivatives_and_high_precision():
    a,b = s.symbols('a b',positive=True)
    f = s.betainc(a,b,0,x)
    equivalent(f.diff(x),x**(a-1)*(1-x)**(b-1))
    equivalent(f.fdiff(4),f.diff(x))
    equivalent(s.riemann_xi(x*x).diff(x),s.riemann_xi(y).diff(y).subs(y,x*x)*2*x)
    equivalent(s.sin(f).diff(x),s.cos(f)*x**(a-1)*(1-x)**(b-1))
    equivalent(s.hermite_prob(s.Symbol('n',integer=True),x).diff(x),s.Symbol('n',integer=True)*s.hermite_prob(s.Symbol('n',integer=True)-1,x))
    with mp.workdps(65):
        beta = s.betainc(s.Rational(1,3),s.Rational(2,3),0,s.Rational(1,2)).evalf(55)
        assert abs(mp.mpf(str(beta))-mp.betainc(mp.mpf(1)/3,mp.mpf(2)/3,0,mp.mpf('.5'))) < mp.mpf('1e-53')
        assert abs(mp.mpf(str(s.stieltjes(1).evalf(55)))-mp.stieltjes(1)) < mp.mpf('1e-53')
        xi = s.riemann_xi(2).evalf(55)
        assert abs(mp.mpf(str(xi))-mp.pi/6) < mp.mpf('1e-53')
    assert s.betainc_regularized(a,b,0,1) == 1
    assert s.sin(x).rewrite() == s.sin(x)
    assert s.sin(x).rewrite(s.sin,s.exp) == s.sin(x).rewrite(s.exp)
    equivalent(s.jacobi_normalized(0,0,0,x),1/s.sqrt(2))
    theta,phi = s.symbols('theta phi',real=True)
    equivalent(s.Znm(0,0,theta,phi),s.Ynm(0,0,theta,phi))
    equivalent(s.Ynm_c(1,1,theta,phi),s.conjugate(s.Ynm(1,1,theta,phi)))


def test_singularity_integrals_and_distribution_products():
    F = s.SingularityFunction
    assert F(3,1,2) == 4 and F(0,1,-1) == 0
    assert F(x,1,3).diff(x) == 3*F(x,1,2)
    assert s.singularityintegrate(F(x,1,-1),x) == F(x,1,0)
    equivalent(s.singularityintegrate(x*F(x,1,-2),x),F(x,1,-1)-F(x,1,0))
    result = s.singularityintegrate(x*F(x,1,1),x)
    for val in (2,3,4):
        equivalent(result.diff(x).subs(x,val), (x*F(x,1,1)).subs(x,val))
    assert s.singularityintegrate(s.sin(x),x) is None
    equivalent(F(x,1,2).rewrite(s.Heaviside),(x-1)**2*s.Heaviside(x-1,1))


def test_splines_partition_and_interpolation():
    basis = s.bspline_basis_set(2,(0,0,0,1,2,3,3,3),x)
    for point in (s.Rational(1,2),s.Rational(3,2),s.Rational(5,2)):
        assert sum(b.subs(x,point) for b in basis) == 1
    for degree in (1,2,3):
        spline = s.interpolating_spline(degree,x,[0,1,2,3],[0,1,8,27])
        for point,value in zip([0,1,2,3],[0,1,8,27]):assert spline.subs(x,point) == value
    with pytest.raises(ValueError):s.bspline_basis(2,[0,1,0,2],0,x)


def test_piecewise_folding_and_exclusive_branches():
    p = s.Piecewise((x,x<0),(x*x,True))
    folded = s.piecewise_fold(2*p+p**2)
    for v in (-2,0,3):assert folded.subs(x,v) == (2*p+p**2).subs(x,v)
    exclusive = s.piecewise_exclusive(s.Piecewise((1,x<1),(2,x<2)))
    for v,w in ((0,1),(1,2),(3,s.nan)):assert exclusive.subs(x,v) == w


def test_variational_calculus_and_line_integral():
    f = s.Function('f')(x)
    equations = s.euler_equations(f.diff(x)**2/2-f**2/2,f,x)
    equivalent(equations[0].lhs,-f-f.diff(x,2))
    higher = s.euler_equations(f.diff(x,2)**2/2,f,x)[0]
    equivalent(higher.lhs,f.diff(x,4))
    curve = s.Curve((x,0),(x,0,1))
    assert s.line_integrate(y,curve,(y,z)) == s.Rational(1,2)


def test_factor_system_triangulation_and_substitution():
    result = s.factor_system([x*x-1,y-x],(x,y))
    assert len(result) == 2
    assert s.factor_system([x,x-1],(x,)) == []
    assert s.factor_system([0],(x,)) == [[]]
    assert set(s.solve_triangulated([x+y-3,x-y-1],x,y)) == {(s.Integer(2),s.Integer(1))}
    assert s.substitution([x+y],[x],[{y:s.Integer(1)}],[y],[],[x,y]) == s.FiniteSet(s.Tuple(-1,1))
    assert s.substitution([x+y],[x],[{y:s.Integer(1)}],[y],[x+1],[x,y]) == s.EmptySet


def test_diophantine_families_are_solutions_and_complete_for_finite_cases():
    for equation,variables in [(2*x+3*y-5,(x,y)),(6*x+10*y+15*z-1,(x,y,z)),(x*y-12,(x,y)),(x*x+y*y-25,(x,y))]:
        rows = s.diophantine(equation,syms=variables)
        assert rows
        for row in rows:equivalent(equation.subs(dict(zip(variables,row)),simultaneous=True),0)
    actual = {(int(a),int(b)) for a,b in s.diophantine(x*x+y*y-25)}
    expected = {(a,b) for a in range(-5,6) for b in range(-5,6) if a*a+b*b == 25}
    assert actual == expected
    assert s.diophantine(2*x+4*y-1) == set()
    for row in s.diophantine(x*x-2*y*y-1):
        for parameter in row[0].free_symbols | row[1].free_symbols:
            for n in range(-3,4):
                values = [v.subs(parameter,n) for v in row]
                equivalent(values[0]**2-2*values[1]**2,1)
    with pytest.raises(NotImplementedError):s.diophantine(x**3+y**3-7)


def test_quaternion_hamilton_identities_rotations_and_calculus():
    Q = s.Quaternion;i,j,k = Q(0,1),Q(0,0,1),Q(0,0,0,1)
    assert i*j == k and j*i == -k and i*i == Q(-1)
    q = Q(1,2,3,4)
    assert (q*q.inverse()).simplify() == Q(1)
    assert q**0 == Q(1)
    assert Q(0,0,0,s.pi).exp().simplify() == Q(-1)
    r = Q.from_axis_angle(((0,0,1),s.pi/2))
    assert Q.rotate_point((1,0,0),r) == (0,1,0)
    matrix = r.to_rotation_matrix()
    assert matrix == s.Matrix([[0,-1,0],[1,0,0],[0,0,1]])
    assert Q.from_rotation_matrix(matrix).to_rotation_matrix() == matrix
    assert Q.from_rotation_matrix(s.diag(1,-1,-1)).to_rotation_matrix() == s.diag(1,-1,-1)
    assert Q(x,x*x).diff(x) == Q(1,2*x)
    with pytest.raises(ValueError):Q().inverse()


def test_utilities_assumptions_and_registered_handlers():
    assert s.capture(lambda:print('hello')) == 'hello\n'
    assert s.filldedent('    hello\n    world').strip() == 'hello world'
    assert s.comp(s.Float('1.0001'),1,tol=.001)
    assert s.arity(lambda x,y=0:None) == (1,2)
    assert s.arity(s.sin) == 1
    assert s.arity(s.atan2) == 2
    assert s.arity(s.Heaviside) == (1,2)
    assert s.arity(s.log) == (1,2)
    @s.vectorize(0,1)
    def plus(a,b):return a+b
    assert plus([1,2],[3,4]) == [[4,5],[5,6]]
    @s.threaded
    def square(e):return e**2
    assert square([x,y]) == [x*x,y*y]
    assert square(x+y) == x*x+y*y
    assert s.xthreaded(lambda e:e**2)(x+y) == (x+y)**2
    assert s.check_assumptions(s.Integer(2),positive=True) is True
    assert s.failing_assumptions(s.Integer(-2),positive=True) == {'positive':False}
    assert s.check_assumptions(x,positive=True) is None
    handler = lambda value,assumptions=True:True if value == 42 else None
    predicate = s.register_handler('answer42',handler)
    assert s.ask(predicate(42)) is True
    s.remove_handler('answer42',handler)
    assert s.ask(predicate(42)) is None
    a,b,c,d = s.symbols('a b c d')
    answer = s.bool_map(s.And(a,b),s.And(c,d))
    assert answer is not False
    assert s.satisfiable(s.Xor(answer[0].xreplace(answer[1]),s.And(c,d))) is False
    equivalent(s.expand_multinomial((x+y+z)**3),(x+y+z)**3)
    count,seconds,scaled,unit = s.timed(lambda:1+1,limit=1)
    assert count == 1 and seconds >= 0 and unit in ('s','ms','μs','ns')


def test_dispatch_strategies_external_and_interactive(monkeypatch):
    from symox.multipledispatch import dispatch
    namespace = {}
    @dispatch(object,namespace=namespace)
    def identify(value):return 'object'
    @dispatch(int,namespace=namespace)
    def identify(value):return 'integer'
    assert identify(2) == 'integer' and identify('x') == 'object'
    from symox.strategies import chain,exhaust,condition
    assert chain(lambda x:x+1,lambda x:x*2)(3) == 8
    assert exhaust(condition(lambda x:x>0,lambda x:x-1))(4) == 0
    assert s.external.import_module('math').sqrt(4) == 2
    assert s.external.import_module('symox_missing_dependency_123') is None
    assert s.concrete.summation(x,(x,1,3)) == 6
    assert s.release.__version__ == s.__version__
    previous = sys.displayhook
    s.init_printing(pretty_print=False)
    assert s.capture(lambda:sys.displayhook(x+1)).strip() in ('1 + x','x + 1')
    monkeypatch.setattr(sys,'displayhook',previous)
    assert 'matplotlib' in s.plot_backends


def test_test_runners_report_real_failures(tmp_path):
    good,bad = tmp_path/'good.txt',tmp_path/'bad.txt'
    good.write_text('>>> 1+1\n2\n')
    bad.write_text('>>> 1+1\n3\n')
    assert s.doctest(str(good)) is True
    assert s.doctest(str(bad)) is False


def test_remaining_surface_operates_without_runtime_sympy():
    import subprocess
    code = r"""
import sys,pickle
class BlockSymPy:
    def find_spec(self,name,*args):
        if name == 'sympy' or name.startswith('sympy.'):
            raise RuntimeError('unexpected runtime SymPy import: '+name)
sys.meta_path.insert(0,BlockSymPy())
import symox as s
from symox.polys.rootoftools import CRootOf,RootSum
from symox.calculus.euler import euler_equations
from symox.functions.special.beta_functions import betainc
from symox.core.parameters import evaluate
from symox.solvers.diophantine import diophantine
x,y = s.symbols('x y')
r = CRootOf(x**5-x-1,0)
assert str(r.evalf(35)).startswith('1.167303978261418684256045899854')
assert pickle.loads(pickle.dumps(r)) == r
assert RootSum(x*x-2,s.Lambda(x,x*x)) == 4
assert betainc(1,1,0,x) == x
assert s.Quaternion(0,1)*s.Quaternion(0,0,1) == s.Quaternion(0,0,0,1)
assert diophantine(x*x+y*y-1)
assert len(s.bspline_basis_set(1,[0,1,2],x)) == 1
assert 'sympy' not in sys.modules
"""
    result = subprocess.run([sys.executable,'-c',code],capture_output=True,text=True)
    assert result.returncode == 0,result.stderr


def test_runner_invokes_pytest_and_propagates_failure(tmp_path):
    import subprocess
    good = tmp_path/'test_good.py'
    bad = tmp_path/'test_bad.py'
    good.write_text('def test_good(): assert 2+2 == 4\n')
    bad.write_text('def test_bad(): assert 2+2 == 5\n')
    for path,expected in ((good,True),(bad,False)):
        result = subprocess.run([sys.executable,'-c',
            'import symox; assert symox.test(%r) is %r' % (str(path),expected)],
            capture_output=True,text=True)
        assert result.returncode == 0,result.stdout+result.stderr


def test_public_decorator_and_canonical_handler_path():
    module = types.ModuleType('symox_test_public_module')
    sys.modules[module.__name__] = module
    try:
        exec('from symox import public\n@public\ndef exposed(): return 42',module.__dict__)
        assert module.__all__ == ['exposed'] and module.exposed() == 42
    finally:sys.modules.pop(module.__name__)
    handler = lambda value,assumptions:True if value == 7 else None
    predicate = s.register_handler('custom_canonical',handler)
    try:
        assert importlib.import_module('symox.assumptions.ask').ask(predicate(7)) is True
    finally:s.remove_handler('custom_canonical',handler)


@pytest.mark.parametrize('degree',[1,2,3,4])
def test_spline_against_reference_at_noninteger_points(degree):
    import sympy as reference
    knots = [0,1,2,4,6,7]
    values = [2,-1,4,0,3,2]
    actual = s.interpolating_spline(degree,x,knots,values)
    t = reference.Symbol('t')
    expected = reference.interpolating_spline(degree,t,knots,values)
    for numerator in range(1,14):
        point = s.Rational(numerator,2)
        assert abs(float(actual.subs(x,point))-float(expected.subs(t,reference.Rational(numerator,2)))) < 1e-12


def test_euler_coupled_field_and_mixed_derivatives():
    f,g = s.Function('f')(x,y),s.Function('g')(x,y)
    equations = s.euler_equations((f.diff(x)**2+f.diff(y)**2)/2-f*g,(f,g),(x,y))
    equivalent(equations[0].lhs,-g-f.diff(x,2)-f.diff(y,2))
    equivalent(equations[1].lhs,-f)
    higher = s.euler_equations(f.diff(x,y)**2/2,f,(x,y))[0]
    equivalent(higher.lhs,f.diff(x,2,y,2))


def test_runner_default_is_source_tests(monkeypatch):
    from pathlib import Path
    observed = []
    monkeypatch.setattr(pytest,'main',lambda args:observed.extend(args) or 0)
    assert s.test() is True
    assert observed == [str(Path(__file__).resolve().parent)]


def test_algebraic_root_scale_and_partial_radical_list():
    assert s.CRootOf(x*(x**5-x-1),0) == 0
    tiny = s.CRootOf(x*x-s.Rational(2,10**60),1)
    with mp.workdps(65):
        actual = mp.mpf(str(tiny.evalf(50)))
        expected = mp.sqrt(2)*mp.mpf('1e-30')
        assert abs(actual/expected-1) < mp.mpf('1e-49')
    repeated = s.CRootOf((x*x-2)**3,3)
    assert abs(float(repeated.evalf(30))-2**.5) < 1e-15
