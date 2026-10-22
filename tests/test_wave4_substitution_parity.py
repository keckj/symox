"""Reference substitution and the nonlinsolve fallback keep every periodic family."""
import itertools

import symox as s

x,y,z = s.symbols('x y z')


def points(result,lo=-3,hi=3):
    """Sample a FiniteSet of tuples, expanding ImageSet components over n in lo..hi."""
    pts = set()
    for row in result.args:
        comps = []
        for e in row.args:
            if isinstance(e,s.ImageSet) or type(e).__name__ == 'ImageSet':
                lam = e.lamda
                comps.append([lam.expr.subs(dict(zip(lam.variables,ns)))
                              for ns in itertools.product(range(lo,hi+1),repeat=len(lam.variables))])
            else:
                comps.append([e])
        for combo in itertools.product(*comps):
            vals = [complex(s.N(v,20)) for v in combo]
            pts.add(tuple((round(v.real,8)+0.0,round(v.imag,8)+0.0) for v in vals))
    return pts


def contains_point(result,point):
    target = tuple((round(complex(v).real,8)+0.0,round(complex(v).imag,8)+0.0) for v in point)
    return target in points(result)


def test_reference_docstring_examples():
    a,b = s.symbols('a b',real=True)
    assert s.substitution([a+b],[a],[{b:1}],[b],set(),[a,b]) == s.FiniteSet(s.Tuple(-1,1))
    assert s.substitution([a+b],[a],[{b:1}],[b],{a+1},[b,a]) == s.EmptySet
    assert s.substitution([a+b],[a],[{b:1}],[b],{a-1},[b,a]) == s.FiniteSet(s.Tuple(1,-1))
    assert s.substitution([a+b-1,b-a**2+5],[a,b]) == s.FiniteSet(s.Tuple(-3,4),s.Tuple(2,-1))


def test_exp_sin_system_keeps_imageset_families():
    for solver in (s.substitution,s.nonlinsolve):
        r = solver([s.exp(x)-s.sin(y),y**2-4],[x,y])
        assert len(r.args) == 2
        assert all(type(row.args[0]).__name__ == 'ImageSet' for row in r.args)
        lsin2 = complex(s.N(s.log(s.sin(2))))
        for n in (-2,0,1):
            assert contains_point(r,(lsin2+2j*n*3.141592653589793,2))
            assert contains_point(r,(lsin2+1j*(2*n+1)*3.141592653589793,-2))


def test_docstring_system_with_free_parameter():
    r = s.substitution([z**2+s.exp(2*x)-s.sin(y),-3+s.exp(-y)],[y,z])
    assert len(r.args) == 4
    families = [row for row in r.args if type(row.args[0]).__name__ == 'ImageSet']
    assert len(families) == 2
    for row in families:
        lam = row.args[0].lamda
        assert s.simplify(s.exp(-lam.expr)-3) == 0


def test_imageset_inside_tuple_exposes_lamda():
    r = s.substitution([s.exp(x)-1,y-x],[x,y])
    assert s.Tuple(0,0) in r.args
    row = next(row for row in r.args if row != s.Tuple(0,0))
    for e in row.args:
        lam = e.lamda
        assert e.base_set == s.S.Integers
        assert s.simplify(s.exp(lam.expr.subs(lam.variables[0],5))-1) == 0


def test_nonlinsolve_on_transcendental_systems_does_not_raise():
    r = s.nonlinsolve([s.exp(x)-y,y**2-1],[x,y])
    assert contains_point(r,(0,1))
    assert contains_point(r,(1j*3.141592653589793,-1))
    assert contains_point(r,(2j*3.141592653589793,1))
    r = s.nonlinsolve([s.exp(x)-2,s.exp(y)-3],[x,y])
    assert contains_point(r,(complex(s.N(s.log(2))),complex(s.N(s.log(3)))+4j*3.141592653589793))


def test_polynomial_nonlinsolve_unchanged():
    r = s.nonlinsolve([x**2-2*y**2-2,x*y-2],[x,y])
    assert r == s.FiniteSet(s.Tuple(-2,-1),s.Tuple(2,1),
                            s.Tuple(-s.sqrt(2)*s.I,s.sqrt(2)*s.I),s.Tuple(s.sqrt(2)*s.I,-s.sqrt(2)*s.I))


def test_no_spurious_free_row_when_equation_constrains_it():
    # The reference returns (exp(y), y) here, which is not a solution for most y.
    r = s.nonlinsolve([s.exp(y)-x,x**2+1],[x,y])
    assert len(r.args) == 2
    for row in r.args:
        assert row.args[0] in (s.I,-s.I)
        lam = row.args[1].lamda
        assert s.simplify(s.exp(lam.expr)-row.args[0]) == 0
    # The reference returns (x, 1), which claims every x.
    r = s.substitution([s.sin(x)**2-1,y-1],[x,y])
    pi = 3.141592653589793
    assert contains_point(r,(pi/2,1))
    assert contains_point(r,(3*pi/2,1))
    assert not contains_point(r,(0,1))


def test_unsolved_equation_gives_conditionset():
    # Scalar solveset leaves sin(x) + cos(x) unsolved; a (x, 0) row would claim every x.
    r = s.substitution([s.sin(x)+s.cos(x),y],[x,y])
    assert isinstance(r,s.ConditionSet)
