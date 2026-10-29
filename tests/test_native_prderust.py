"""Native base-case parametric Risch DE kernels; unchanged complete old fallback."""
import importlib

import pytest
import symox as sx

prde=importlib.import_module('symox.integrals.prde')
rde=importlib.import_module('symox.integrals.rde')
risch=importlib.import_module('symox.integrals.risch')
native=importlib.import_module('symox.integrals._prde_native')
polymatrix=importlib.import_module('symox.polys.polymatrix')

Poly,Matrix,Rational,S=sx.Poly,polymatrix.PolyMatrix,sx.Rational,sx.S
DifferentialExtension,derivation=risch.DifferentialExtension,risch.derivation
x,t,n=sx.symbols('x t n')
QQ=polymatrix.QQ

ENTRIES=('param_rischDE','param_poly_rischDE','prde_no_cancel_b_large','prde_no_cancel_b_small',
         'prde_spde','prde_normal_denom','prde_special_denom','prde_linear_constraints',
         'poly_linear_constraints','constant_system')


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
    """Record the native kernels that ran and returned a result."""
    seen=[];call=native._call

    def spy(kind,args,fields):
        r=call(kind,args,fields)
        if r is not None:seen.append(kind)
        return r
    monkeypatch.setattr(native,'_call',spy)
    original=native._ROOT[0]

    def spy_root(kind,*args):
        r=original(kind,*args)
        if r is not None and kind=='poly_linear_constraints':seen.append(kind)
        return r
    monkeypatch.setattr(native,'_ROOT',(spy_root,)+native._ROOT[1:])
    monkeypatch.setattr(native,'_verified',lambda:native._ROOT)
    return seen


def python_only(monkeypatch,f,*args):
    with monkeypatch.context()as patch:
        patch.setattr(native,'_ROOT',None)
        return f(*args)


def base():
    return DifferentialExtension(extension={'D':[Poly(1,x)]})


# SymPy 1.14 sympy/integrals/tests/test_prde.py, base-case examples with copied expected values.

def test_sympy_prde_linear_constraints(kinds):
    G=[(Poly(2*x**3+3*x+1,x),Poly(x**2-1,x)),(Poly(1,x),Poly(x-1,x)),(Poly(1,x),Poly(x+1,x))]
    assert prde.prde_linear_constraints(Poly(1,x),Poly(0,x),G,base())==\
        ((Poly(2*x,x,domain='QQ'),Poly(0,x,domain='QQ'),Poly(0,x,domain='QQ')),
         Matrix([[1,1,-1],[5,1,1]],x))
    assert kinds==['prde_linear_constraints']


def test_sympy_prde_no_cancel_b_large(kinds):
    DE=base()
    assert prde.prde_no_cancel_b_large(Poly(1,x),[Poly(x**2,x),Poly(1,x)],2,DE)==\
        ([Poly(x**2-2*x+2,x),Poly(1,x)],Matrix([[1,0,-1,0],[0,1,0,-1]],x))
    assert prde.prde_no_cancel_b_large(Poly(1,x),[Poly(x**3,x),Poly(1,x)],3,DE)==\
        ([Poly(x**3-3*x**2+6*x-6,x),Poly(1,x)],Matrix([[1,0,-1,0],[0,1,0,-1]],x))
    assert prde.prde_no_cancel_b_large(Poly(x,x),[Poly(x**2,x),Poly(1,x)],1,DE)==\
        ([Poly(x,x,domain='ZZ'),Poly(0,x,domain='ZZ')],Matrix([[1,-1,0,0],[1,0,-1,0],[0,1,0,-1]],x))
    assert kinds.count('prde_no_cancel_b_large')==3


def test_sympy_param_poly_rischDE(kinds):
    DE=base()
    a=Poly(x**2-x,x,field=True);b=Poly(1,x,field=True)
    q=[Poly(x,x,field=True),Poly(x**2,x,field=True)]
    h,A=prde.param_poly_rischDE(a,b,q,3,DE)
    assert A.nullspace()==[Matrix([0,1,1,1],DE.t)]
    assert h[0]+h[1]==Poly(x,x,domain='QQ')
    b=Poly(x**2-5*x+3,x,field=True)
    q=[Poly(1,x,field=True),Poly(x,x,field=True),Poly(x**2,x,field=True)]
    h,A=prde.param_poly_rischDE(a,b,q,3,DE)
    assert A.nullspace()==[Matrix([3,-5,1,-5,1,1],DE.t)]
    p=-Poly(5,DE.t)*h[0]+h[1]+h[2]
    assert a*derivation(p,DE)+b*p==Poly(x**2-5*x+3,x,domain='QQ')
    assert kinds.count('param_poly_rischDE')==2


def test_sympy_param_rischDE(kinds):
    DE=base()
    p1,px=Poly(1,x,field=True),Poly(x,x,field=True)
    G=[(p1,px),(p1,p1),(px,p1)]
    h,A=prde.param_rischDE(-p1,Poly(x**2,x,field=True),G,DE)
    assert len(h)==3
    p=[hi[0].as_expr()/hi[1].as_expr()for hi in h]
    V=A.nullspace()
    assert len(V)==2
    assert V[0]==Matrix([-1,1,0,-1,1,0],DE.t)
    y=-p[0]+p[1]+0*p[2]
    assert y.diff(x)-y/x**2==1-1/x
    assert kinds==['param_rischDE']


def test_sympy_limited_integrate(kinds):
    DE=base()
    G=[(Poly(x,x),Poly(x+1,x))]
    assert prde.limited_integrate(Poly(-(1+x+5*x**2-3*x**3),x),Poly(1-x-x**2+x**3,x),G,DE)==\
        ((Poly(x**2-x+2,x),Poly(x-1,x,domain='QQ')),[2])
    G=[(Poly(1,x),Poly(x,x))]
    assert prde.limited_integrate(Poly(5*x**2,x),Poly(3,x),G,DE)==\
        ((Poly(5*x**3/9,x),Poly(1,x,domain='QQ')),[0])
    assert kinds.count('param_rischDE')==2


def test_sympy_non_base_examples_keep_python_results(kinds):
    DE=DifferentialExtension(extension={'D':[Poly(1,x),Poly(t,t)]})
    G=[(Poly(t,t),Poly(1,t)),(Poly(t**2,t),Poly(1,t)),(Poly(t**3,t),Poly(1,t))]
    assert prde.prde_linear_constraints(Poly(t+1,t),Poly(t**2,t),G,DE)==\
        ((Poly(t,t,domain='QQ'),Poly(t**2,t,domain='QQ'),Poly(t**3,t,domain='QQ')),Matrix(0,3,[],t))
    DE=DifferentialExtension(extension={'D':[Poly(1,x),Poly(1/x,t)]})
    G=[(Poly(2*x,t),Poly(t,t)),(Poly(-x,t),Poly(t,t))]
    assert prde.prde_linear_constraints(Poly(1,t),Poly(0,t),G,DE)==\
        ((Poly(0,t,domain='QQ[x]'),Poly(0,t,domain='QQ[x]')),Matrix([[2*x,-x]],t))
    D=[Poly(x,t),Poly(-x*t,t)]
    assert prde.prde_spde(Poly(t,t),Poly(-1/x,t),D,n,DE)==\
        (Poly(t,t),Poly(0,t,domain='ZZ(x)'),[Poly(2*x,t,domain='ZZ(x)'),Poly(-x,t,domain='ZZ(x)')],
         [Poly(-x**2,t,domain='ZZ(x)'),Poly(0,t,domain='ZZ(x)')],n-1)
    DE=DifferentialExtension(extension={'D':[Poly(1,x),Poly(1+t**2,t)]})
    G=[(Poly(t,t),Poly(1+t**2,t)),(Poly(1,t),Poly(x+x*t**2,t))]
    assert prde.prde_normal_denom(Poly(1,t),Poly(x,t),G,DE)==\
        (Poly(x,t,domain='ZZ(x)'),(Poly(1,t,domain='ZZ(x)'),Poly(1,t,domain='ZZ(x)')),
         [(Poly(x*t,t,domain='ZZ(x)'),Poly(t**2+1,t,domain='ZZ(x)')),
          (Poly(1,t,domain='ZZ(x)'),Poly(t**2+1,t,domain='ZZ(x)'))],Poly(1,t,domain='ZZ(x)'))
    assert kinds==[]


# Every kernel call reached from risch_integrate, compared with the complete Python body.

INTEGRANDS=['log(x)**2','x*log(x)','log(x)**3/x','x**2*log(x)**2','x*log(x+1)','log(x**2+1)',
            'x*log(x)**3','log(x)**4','(x+1)*log(x)**2','log(x)**2/x**2','log(x)/x**2',
            'log(x)**2/x','x**3*log(x)','log(x)/(x+1)**2','log(2*x+3)**2','x*log(x**2-1)',
            'log(x)*log(x+1)','exp(x)*x','x**2*exp(x)','exp(x)*log(x)','log(x)/(x**2+1)',
            '(log(x)+1)/x**2','exp(2*x)*x**2','log(x)**2*x/3','log(x+1)**2/(x+1)',
            'exp(x)/x+exp(x)*log(x)','x*exp(x**2)','atan(x)','x*atan(x)','log(log(x))/x']


@pytest.mark.parametrize('text',INTEGRANDS)
def test_risch_integrate_kernel_calls_match_python(monkeypatch,text):
    xs=sx.Symbol('x');f=sx.sympify(text,locals={'x':xs})
    records=[];depth=[0]
    for name in ENTRIES:
        entry=getattr(native,name);body=getattr(prde,name)

        def wrapper(*args,_entry=entry,_body=body,_name=name):
            if depth[0]:return None
            fresh=tuple(list(a)if type(a)is list else a for a in args)
            r=_entry(*args)
            if r is None:return None
            depth[0]+=1
            try:
                with monkeypatch.context()as patch:
                    patch.setattr(native,'_ROOT',None)
                    expected=_body(*fresh)
            finally:
                depth[0]-=1
            records.append(_name)
            assert sig(r)==sig(expected),_name
            assert sig(list(args))==sig(list(fresh)),_name
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
        assert sx.srepr(actual)==sx.srepr(expected)
    if text in INTEGRANDS[:10]:
        assert records


def test_integrands_reach_every_first_wave_kernel(monkeypatch):
    seen=set();call=native._call

    def spy(kind,args,fields):
        r=call(kind,args,fields)
        if r is not None:seen.add(kind)
        return r
    monkeypatch.setattr(native,'_call',spy)
    xs=sx.Symbol('x')
    for text in INTEGRANDS[:27]:
        risch.risch_integrate(sx.sympify(text,locals={'x':xs}),xs)
    assert {'param_rischDE','param_poly_rischDE','prde_normal_denom','prde_special_denom',
            'prde_linear_constraints','constant_system'}<=seen


# Direct kernels against the Python bodies on generated base-case inputs.

def generated(seed):
    import random
    rng=random.Random(9127+seed)

    def poly(deg,field):
        c=[Rational(rng.randrange(-9,10),rng.randrange(1,5)if field else 1)for _ in range(deg+1)]
        return Poly(sum(ci*x**i for i,ci in enumerate(c)),x,domain='QQ'if field else'ZZ')
    field=seed%2==1
    return poly,field,rng


@pytest.mark.parametrize('seed',range(24))
def test_generated_kernels_match_python(monkeypatch,seed):
    poly,field,rng=generated(seed)
    DE=base()
    b=poly(rng.randrange(0,3),field)
    Q=[poly(rng.randrange(0,4),field)for _ in range(1+seed%3)]
    nn=rng.randrange(-1,5)
    cases=[('prde_no_cancel_b_large',(b,Q,nn,DE)),('prde_no_cancel_b_small',(b,Q,nn,DE)),
           ('param_poly_rischDE',(poly(rng.randrange(0,3),field),b,Q,nn,DE)),
           ('prde_linear_constraints',(b,b,[(q,poly(rng.randrange(0,3),field)+Poly(x**3,x))for q in Q],DE))]
    a=poly(rng.randrange(0,3),field)
    if not a.is_zero:
        cases.append(('prde_spde',(a,b,Q,nn,DE)))
    fd=poly(rng.randrange(0,3),field)
    if not fd.is_zero:
        G=[(q,fd*q+Poly(x**2+1,x))for q in Q]
        cases.append(('prde_normal_denom',(b,fd,G,DE)))
        cases.append(('param_rischDE',(b,fd,G,DE)))
    for name,args in cases:
        f=getattr(prde,name)
        fresh=lambda:tuple(list(v)if type(v)is list else v for v in args)
        try:
            expected=python_only(monkeypatch,f,*fresh())
        except Exception as exc:
            with pytest.raises(type(exc)):
                f(*fresh())
            continue
        q_new=fresh();actual=f(*q_new)
        q_old=fresh();python_only(monkeypatch,f,*q_old)
        assert sig(actual)==sig(expected),name
        assert sig(list(q_new))==sig(list(q_old)),name


def test_constant_system_and_poly_linear_constraints(monkeypatch,kinds):
    DE=base()
    A=Matrix([[1,2,3],[2,4,Rational(1,3)],[0,0,5]],x)
    u=Matrix([[1],[Rational(2,7)],[0]],x)
    expected=python_only(monkeypatch,prde.constant_system,A,u,DE)
    assert sig(prde.constant_system(A,u,DE))==sig(expected)
    p=[Poly(x**3+2,x),Poly(Rational(1,2)*x,x),Poly(0,x)]
    d=Poly(x**2-1,x)
    expected=python_only(monkeypatch,prde.poly_linear_constraints,p,d)
    assert sig(prde.poly_linear_constraints(p,d))==sig(expected)
    assert kinds==['constant_system','poly_linear_constraints']


# Protocol: any change to a consumed function, global or class keeps the Python body.

@pytest.mark.parametrize('mode',['rde','risch','polymatrix','global','builtin','code','cls','poly'])
def test_changed_protocol_declines(monkeypatch,mode):
    calls=[]
    if mode=='rde':
        original=rde.weak_normalizer
        monkeypatch.setattr(rde,'weak_normalizer',lambda *a:(calls.append(1),original(*a))[1])
        monkeypatch.setattr(prde,'weak_normalizer',rde.weak_normalizer)
    elif mode=='risch':
        original=risch.gcdex_diophantine
        monkeypatch.setattr(prde,'gcdex_diophantine',lambda *a:(calls.append(1),original(*a))[1])
    elif mode=='polymatrix':
        original=Matrix.nullspace
        monkeypatch.setattr(Matrix,'nullspace',lambda self:(calls.append(1),original(self))[1])
    elif mode=='global':
        monkeypatch.setattr(prde,'zeros',Matrix.zeros)
    elif mode=='builtin':
        monkeypatch.setitem(prde.__dict__,'len',lambda v:(calls.append(1),len(v))[1])
    elif mode=='code':
        f=prde.prde_spde
        monkeypatch.setattr(f,'__code__',f.__code__.replace())
    elif mode=='cls':
        monkeypatch.setattr(DifferentialExtension,'extra_attribute',1,raising=False)
    elif mode=='poly':
        monkeypatch.setattr(Poly,'extra_attribute',1,raising=False)
    assert native._verified()is None
    DE=base()
    h,A=prde.param_poly_rischDE(Poly(x**2-x,x,field=True),Poly(1,x,field=True),
                                [Poly(x,x,field=True),Poly(x**2,x,field=True)],3,DE)
    assert A.nullspace()==[Matrix([0,1,1,1],DE.t)]
    if mode in('risch','polymatrix','builtin'):
        assert calls


def test_stock_protocol_is_restored_after_mutation(monkeypatch):
    with monkeypatch.context()as patch:
        patch.setattr(prde,'zeros',Matrix.zeros)
        assert native._verified()is None
    assert native._verified()is not None


@pytest.mark.parametrize('mode',['subclass','property','case','gens','level'])
def test_extension_reads_decline(monkeypatch,mode):
    DE=base()
    if mode=='subclass':
        class Sub(DifferentialExtension):
            __slots__=()
        DE=object.__new__(Sub)
        for name in DifferentialExtension.__slots__:
            try:setattr(DE,name,getattr(base(),name))
            except AttributeError:pass
    elif mode=='property':
        monkeypatch.setattr(DifferentialExtension,'t',property(lambda self:x))
    elif mode=='case':
        DE=DifferentialExtension(extension={'D':[Poly(1,x),Poly(t,t)]})
    elif mode=='gens':
        DE.D=[Poly(1,x,t)]
    elif mode=='level':
        DE.level=-2
    assert native.param_poly_rischDE(Poly(x,x),Poly(1,x),[Poly(x,x)],2,DE)is None


@pytest.mark.parametrize('mode',['float','symbolic','extra_gen','tuple','bign'])
def test_nonadmitted_payloads_decline(mode):
    DE=base()
    b,Q,nn=Poly(1,x),[Poly(x**2,x),Poly(1,x)],2
    if mode=='float':b=Poly(1.5,x)
    elif mode=='symbolic':b=Poly(t,x)
    elif mode=='extra_gen':b=Poly(x,x,t)
    elif mode=='tuple':Q=tuple(Q)
    elif mode=='bign':nn=1<<40
    assert native.prde_no_cancel_b_large(b,Q,nn,DE)is None


def test_float_payload_keeps_python_error():
    DE=base()
    with pytest.raises(sx.CoercionFailed):
        prde.prde_no_cancel_b_large(Poly(1.5,x),[Poly(x**2,x),Poly(1,x)],2,DE)


def test_no_cancel_mutates_q_like_python(monkeypatch):
    DE=base()
    Q1=[Poly(x**2,x),Poly(1,x)];Q2=list(Q1)
    prde.prde_no_cancel_b_large(Poly(1,x),Q1,2,DE)
    python_only(monkeypatch,prde.prde_no_cancel_b_large,Poly(1,x),Q2,2,DE)
    assert sig(Q1)==sig(Q2)
