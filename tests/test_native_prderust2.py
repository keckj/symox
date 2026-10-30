"""Second prde wave: limited_integrate_reduce (base case), constant_system in
every DE state, and the parallel_poly_from_expr field domain used by the
Risch tower. Expected values are copied from SymPy 1.14."""
import importlib

import pytest
import symox as sx

prde=importlib.import_module('symox.integrals.prde')
risch=importlib.import_module('symox.integrals.risch')
native=importlib.import_module('symox.integrals._prde_native')
polymatrix=importlib.import_module('symox.polys.polymatrix')
polytools=importlib.import_module('symox.polys.polytools')

Poly,Matrix,Rational=sx.Poly,polymatrix.PolyMatrix,sx.Rational
DifferentialExtension=risch.DifferentialExtension
x,t=sx.symbols('x t')


def sig(v):
    if isinstance(v,Poly):
        return ('Poly',tuple(str(g)for g in v.gens),str(v.domain),
                tuple(sorted((m,sx.srepr(c))for m,c in v.as_dict().items())))
    if isinstance(v,Matrix):
        return ('Matrix',v._rows,v._cols,tuple(str(g)for g in v.gens),str(v.domain),
                tuple(sx.srepr(e)for e in v._items))
    if isinstance(v,(list,tuple)):
        return (type(v).__name__,tuple(sig(e)for e in v))
    if isinstance(v,int):
        return ('int',v)
    return ('expr',sx.srepr(v))


@pytest.fixture
def kinds(monkeypatch):
    seen=[];call=native._call

    def spy(kind,args,fields):
        r=call(kind,args,fields)
        if r is not None:seen.append(kind)
        return r
    monkeypatch.setattr(native,'_call',spy)
    return seen


def python_only(monkeypatch,f,*args):
    with monkeypatch.context()as patch:
        patch.setattr(native,'_ROOT',None)
        return f(*args)


def base():
    return DifferentialExtension(extension={'D':[Poly(1,x)]})


def test_limited_integrate_reduce_base_values(kinds):
    DE=base()
    r=prde.limited_integrate_reduce(Poly(x,x),Poly(x**2,x),[(Poly(1,x),Poly(x,x))],DE)
    assert sig(r)==sig((Poly(x,x),Poly(-1,x),Poly(x,x),1,(Poly(x,x),Poly(1,x)),
                        [(Poly(-x,x),Poly(1,x))]))
    assert r[0] is r[2]
    G=[(Poly(x/3,x),Poly(2*x+1,x)),(Poly(5,x),Poly(x**3,x))]
    r=prde.limited_integrate_reduce(Poly(Rational(1,2)*x+3,x),Poly(3*x**2-x,x),G,DE)
    assert sig(r)==sig((Poly(x**2,x),Poly(-2*x,x),Poly(x**2,x),3,
                        (Poly(x**4+6*x**3,x,domain='QQ'),Poly(6*x-2,x,domain='QQ')),
                        [(Poly(-x**5,x,domain='QQ'),Poly(6*x+3,x,domain='QQ')),
                         (Poly(-5*x,x),Poly(1,x))]))
    assert kinds==['limited_integrate_reduce']*2


def test_sympy_limited_integrate_reduce_primitive_keeps_python(kinds):
    # SymPy 1.14 test_prde.py::test_limited_integrate_reduce.
    DE=DifferentialExtension(extension={'D':[Poly(1,x),Poly(1/x,t)]})
    r=prde.limited_integrate_reduce(Poly(x,t),Poly(t**2,t),[(Poly(x,t),Poly(t,t))],DE)
    assert r==(Poly(t,t),Poly(-1/x,t),Poly(t,t),1,(Poly(x,t),Poly(1,t,domain='ZZ[x]')),
               [(Poly(-x*t,t),Poly(1,t,domain='ZZ[x]'))])
    assert kinds==[]


@pytest.mark.parametrize('seed',range(16))
def test_generated_limited_integrate_reduce_matches_python(monkeypatch,seed):
    import random
    rng=random.Random(4401+seed);field=seed%2==1

    def poly(deg):
        c=[Rational(rng.randrange(-6,7),rng.randrange(1,4)if field else 1)for _ in range(deg+1)]
        return Poly(sum(ci*x**i for i,ci in enumerate(c)),x,domain='QQ'if field else'ZZ')
    fd=poly(rng.randrange(0,3))+Poly(x**3,x)
    G=[(poly(rng.randrange(0,3)),poly(rng.randrange(0,2))+Poly(x**2,x))for _ in range(1+seed%3)]
    args=(poly(rng.randrange(0,4)),fd,G,base())
    expected=python_only(monkeypatch,prde.limited_integrate_reduce,*args)
    assert sig(prde.limited_integrate_reduce(*args))==sig(expected)


def test_constant_system_primitive_case(monkeypatch,kinds):
    DE=DifferentialExtension(extension={'D':[Poly(1,x),Poly(1/x,t)]})
    A=Matrix([[1,2,Rational(1,2)],[3,4,5]],t);u=Matrix([[1],[Rational(2,3)]],t)
    r=prde.constant_system(A,u,DE)
    assert r==(Matrix([[1,0,4],[0,1,Rational(-7,4)]],t),Matrix([[Rational(-4,3)],[Rational(7,6)]],t))
    assert sig(r)==sig(python_only(monkeypatch,prde.constant_system,A,u,DE))
    assert kinds==['constant_system']


def test_constant_system_declines_symbolic_entries_and_bad_tower(monkeypatch,kinds):
    DE=DifferentialExtension(extension={'D':[Poly(1,x),Poly(1/x,t)]})
    A=Matrix([[x,1],[1,2]],t);u=Matrix([[1],[0]],t)
    assert sig(prde.constant_system(A,u,DE))==sig(python_only(monkeypatch,prde.constant_system,A,u,DE))
    A=Matrix([[1,2]],t);u=Matrix([[1],[2]],t)
    with pytest.raises(Exception)as info:
        python_only(monkeypatch,prde.constant_system,A,u,DE)
    with pytest.raises(type(info.value)):
        prde.constant_system(A,u,DE)
    assert kinds==[]


def test_parallel_poly_from_expr_field_domain():
    t0=sx.Symbol('_t0')
    for items,dom in (([0,x],'ZZ(x)'),([0,1/x],'ZZ(x)'),([Rational(1,2),x],'ZZ(x)'),
                      ([x,1/x+1],'ZZ(x)'),([Rational(1,3),2],'QQ')):
        ps,_=polytools.parallel_poly_from_expr(items,t0,field=True)
        assert [str(p.domain)for p in ps]==[dom,dom]


def test_log_log_integrand_matches_reference():
    # SymPy 1.14: x*log(log(x)) + NonElementaryIntegral(-1/log(x), x); the
    # base wave raised CoercionFailed (parallel_poly_from_expr field domain).
    r=risch.risch_integrate(sx.log(sx.log(x)),x)
    assert str(r)=='x*log(log(x)) + Integral(-1/log(x), x)'
    assert sorted(type(a).__name__ for a in r.args)==['Mul','NonElementaryIntegral']


# Tower integrands: every admitted constant_system call (base, exp, primitive and
# DE construction) against the Python body, then the whole result against the
# Python-only run.
TOWER=['exp(x)*x','exp(x**2)','x*log(x)','log(x)**2','exp(x)/x','1/log(x)','log(log(x))',
       'exp(exp(x))','x*exp(x)*log(x)','(x+1)*exp(x)/(x*exp(x)+1)','exp(x)/(exp(2*x)+1)',
       'log(x)/x','log(x)**3/x','exp(x)*log(x)','exp(2*x)+exp(x)','1/(exp(x)+1)','x**2*exp(x**2)',
       'log(x+1)/(x+1)','exp(x)*(1+1/x)*log(x)','log(x)*exp(log(x)**2)/x','exp(-x**2)*x',
       '1/(x*log(x))','exp(x)/(1+exp(x))**2','x*exp(x)**2','log(x**2+1)','exp(x)*x**3',
       'log(x)**2/x**2','exp(3*x)*x**2','1/(x*log(x)*log(log(x)))','exp(x)/x**2',
       'log(x)/(x*log(log(x)))','x*log(x)**2']


@pytest.mark.parametrize('text',TOWER)
def test_tower_integrands_match_python(monkeypatch,text):
    xs=sx.Symbol('x');f=sx.sympify(text,locals={'x':xs})
    records=[];depth=[0]
    for name in ('constant_system','limited_integrate_reduce'):
        entry=getattr(native,name);body=getattr(prde,name)

        def wrapper(*args,_entry=entry,_body=body,_name=name):
            if depth[0]:return None
            r=_entry(*args)
            if r is None:return None
            depth[0]+=1
            try:
                with monkeypatch.context()as patch:
                    patch.setattr(native,'_ROOT',None)
                    expected=_body(*args)
            finally:
                depth[0]-=1
            records.append(_name)
            assert sig(r)==sig(expected),_name
            return r
        monkeypatch.setattr(native,name,wrapper)
    try:
        actual=risch.risch_integrate(f,xs)
    except Exception as exc:
        actual=exc
    monkeypatch.undo()
    try:
        expected=python_only(monkeypatch,risch.risch_integrate,f,xs)
    except Exception as exc:
        assert type(actual)is type(exc)and str(actual)==str(exc)
    else:
        # str: RootSum results carry fresh Dummy indices on every run.
        assert type(actual)is type(expected) and str(actual)==str(expected)

