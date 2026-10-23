"""Global polynomial division honors construction and rendering options."""
import pickle
import subprocess
import sys

import pytest
import symox as sm

OPS = ['div', 'rem', 'quo', 'exquo']


def compare(actual, expected):
    import sympy as sy
    if isinstance(expected, tuple):
        assert isinstance(actual, tuple) and len(actual) == len(expected)
        for a, e in zip(actual, expected):
            compare(a, e)
    elif isinstance(expected, sy.Poly):
        assert isinstance(actual, sm.Poly)
        assert tuple(sm.to_sympy(g) for g in actual.gens) == expected.gens
        assert sy.expand(sm.to_sympy(actual.as_expr())-expected.as_expr()) == 0
        assert str(actual.domain) == str(expected.domain)
    else:
        assert not isinstance(actual, sm.Poly)
        assert sy.expand(sm.to_sympy(actual)-expected) == 0


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('kind', ['expr', 'poly', 'mixed'])
@pytest.mark.parametrize('auto', [False, True])
@pytest.mark.parametrize('polys', [None, False, True])
def test_global_auto_polys_and_exact_quotient_oracle(op, kind, auto, polys):
    import sympy as sy
    x = sm.Symbol('x')
    sx = sm.to_sympy(x)
    f, g = x+1, 2*x+2
    of, og = sx+1, 2*sx+2
    if kind != 'expr':
        f, of = sm.Poly(f,x), sy.Poly(of,sx)
    if kind == 'poly':
        g, og = sm.Poly(g,x), sy.Poly(og,sx)
    kwargs = {'auto':auto, 'polys':polys}
    try:
        expected = getattr(sy,op)(of,og,*((sx,) if kind == 'expr' else ()),**kwargs)
    except sy.polys.polyerrors.ExactQuotientFailed:
        with pytest.raises(sm.ExactQuotientFailed):
            getattr(sm,op)(f,g,*((x,) if kind == 'expr' else ()),**kwargs)
    else:
        compare(getattr(sm,op)(f,g,*((x,) if kind == 'expr' else ()),**kwargs), expected)


def selected(kind, order='lex'):
    import sympy as sy
    a = sm.Symbol('a', positive=True)
    b = sm.Dummy('unused', real=True)
    sa, sb = sm.to_sympy(a), sm.to_sympy(b)
    native = {'ZZ':sm.ZZ, 'QQ':sm.QQ, 'GF':sm.GF(5), 'AN':sm.QQ.algebraic_field(sm.sqrt(2)),
              'QQRF':sm.QQ.frac_field(a,b,order=order), 'ZZRF':sm.ZZ.frac_field(a,b,order=order),
              'QQring':sm.QQ.poly_ring(a,b,order=order), 'ZZring':sm.ZZ.poly_ring(a,b,order=order)}[kind]
    oracle = {'ZZ':sy.ZZ, 'QQ':sy.QQ, 'GF':sy.GF(5), 'AN':sy.QQ.algebraic_field(sy.sqrt(2)),
              'QQRF':sy.QQ.frac_field(sa,sb,order=order), 'ZZRF':sy.ZZ.frac_field(sa,sb,order=order),
              'QQring':sy.QQ.poly_ring(sa,sb,order=order), 'ZZring':sy.ZZ.poly_ring(sa,sb,order=order)}[kind]
    return native, oracle, a, b


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('kind', ['ZZ','QQ','GF','AN','QQRF','ZZRF','QQring','ZZring'])
@pytest.mark.parametrize('auto', [False, True])
def test_explicit_typed_domain_and_dense_recursive_division(op, kind, auto):
    import sympy as sy
    x,y = sm.symbols('x y')
    k, ok, a, unused = selected(kind, 'grevlex')
    c = a if kind.endswith('RF') or kind.endswith('ring') else sm.sqrt(2) if kind == 'AN' else sm.Integer(2)
    g = c*x+c
    f = (x+y)*g + (0 if op == 'exquo' else 1)
    kwargs = {'domain':k, 'auto':auto, 'polys':True}
    expected = getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(g),sm.to_sympy(x),sm.to_sympy(y),domain=ok,auto=auto,polys=True)
    actual = getattr(sm,op)(f,g,x,y,**kwargs)
    compare(actual, expected)
    results = actual if op == 'div' else (actual,)
    for value in results:
        assert value.domain == k
        if kind.endswith('RF') or kind.endswith('ring'):
            assert value.domain.symbols == (a,unused)
            assert value.domain.order == 'grevlex'
        assert all(coef in value.domain for coef in value.coeffs())
        restored = pickle.loads(pickle.dumps(value))
        assert restored.domain == value.domain and restored.as_expr() == value.as_expr()
    q,r = sm.div(f,g,x,y,domain=k,auto=auto,polys=True)
    residue = sm.expand(q.as_expr()*g+r.as_expr()-f)
    assert residue == 0 if kind != 'GF' else sm.Poly(residue,x,y,domain=k).is_zero


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('order', ['lex','grlex','grevlex'])
@pytest.mark.parametrize('auto', [False, True])
def test_ring_promotion_retraction_and_parameter_identity(op, order, auto):
    import sympy as sy
    x = sm.Symbol('x')
    k, ok, a, b = selected('ZZring',order)
    f,g = x+1, a*x+a
    try:
        expected=getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(g),sm.to_sympy(x),domain=ok,auto=auto,polys=True)
    except sy.polys.polyerrors.ExactQuotientFailed:
        with pytest.raises(sm.ExactQuotientFailed):
            getattr(sm,op)(f,g,x,domain=k,auto=auto,polys=True)
    else:
        actual=getattr(sm,op)(f,g,x,domain=k,auto=auto,polys=True)
        compare(actual,expected)
        for p in actual if op == 'div' else (actual,):
            assert p.domain.symbols == (a,b) and p.domain.order == order
            assert all(c in p.domain for c in p.coeffs())
    if auto:
        q,r=sm.div(f,g,x,domain=k,auto=True,polys=True)
        assert q.domain == r.domain == k.get_field()
        assert sm.rem(f,g,x,domain=k,polys=True).domain == k


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('keyword', [False, True])
def test_explicit_generators_reframe_poly_inputs(op, keyword):
    import sympy as sy
    x,y=sm.symbols('x y')
    sx,syy=sm.to_sympy(x),sm.to_sympy(y)
    f,g=sm.Poly(x+1,x),sm.Poly(2*x+2,x)
    of,og=sy.Poly(sx+1,sx),sy.Poly(2*sx+2,sx)
    actual=getattr(sm,op)(f,g,**{'gens':(y,),'polys':True}) if keyword else getattr(sm,op)(f,g,y,polys=True)
    expected=getattr(sy,op)(of,og,**{'gens':(syy,),'polys':True}) if keyword else getattr(sy,op)(of,og,syy,polys=True)
    compare(actual,expected)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('kind', ['ZZ','QQ','GF','AN','QQRF','QQring'])
def test_existing_typed_poly_inputs_keep_coefficient_descriptor(op, kind):
    x=sm.Symbol('x')
    k,_,a,b=selected(kind,'grlex')
    c=a if kind in ['QQRF','QQring'] else sm.sqrt(2) if kind == 'AN' else sm.Integer(2)
    f=sm.Poly((x+c)*(x+1),x,domain=k)
    g=sm.Poly(x+c,x,domain=k)
    result=getattr(sm,op)(f,g)
    for p in result if op == 'div' else (result,):
        assert isinstance(p,sm.Poly) and p.domain==k
        assert p.gens==(x,)
    expr=getattr(sm,op)(f,g,polys=False)
    for value in expr if op == 'div' else (expr,):
        assert not isinstance(value,sm.Poly)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('kind', ['ZZ','QQring','ZZring'])
def test_field_flag_promotes_before_auto_false(op, kind):
    import sympy as sy
    x=sm.Symbol('x')
    k,ok,a,b=selected(kind)
    f,g=sm.Poly(x+1,x,domain=k),sm.Poly(2*x+2,x,domain=k)
    of,og=sy.Poly(sm.to_sympy(x+1),sm.to_sympy(x),domain=ok),sy.Poly(sm.to_sympy(2*x+2),sm.to_sympy(x),domain=ok)
    compare(getattr(sm,op)(f,g,field=True,auto=False),getattr(sy,op)(of,og,field=True,auto=False))


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('option', ['bogus','auto','polys','field','order','extension'])
def test_option_errors_are_native_named_and_pickleable(op, option):
    from symox.polys.polyerrors import OptionError
    x=sm.Symbol('x')
    value={'bogus':None,'auto':2,'polys':'yes','field':2,'order':'potato','extension':False}[option]
    with pytest.raises(OptionError) as failure:
        getattr(sm,op)(x+1,x,x,**{option:value})
    assert sm.OptionError is OptionError
    restored=pickle.loads(pickle.dumps(failure.value))
    assert type(restored) is OptionError and str(restored)==str(failure.value)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('case', range(6))
def test_contradictory_options_are_not_silently_ignored(op, case):
    x=sm.Symbol('x')
    kwargs=[{'domain':sm.ZZ,'field':False},{'domain':sm.GF(5),'modulus':5},
            {'domain':sm.QQ,'extension':sm.sqrt(2)},{'extension':True,'modulus':5},
            {'gens':(x,)},{'gens':(x,x)}][case]
    args=(x,) if case != 5 else ()
    with pytest.raises(sm.OptionError):
        getattr(sm,op)(x+1,x,*args,**kwargs)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('case', range(4))
def test_incompatible_poly_domains_fail_explicitly(op, case):
    x=sm.Symbol('x');a=sm.Symbol('a',positive=True)
    domains=[(sm.GF(5),sm.GF(7)),(sm.GF(5),sm.QQ),
             (sm.QQ.frac_field(a),sm.QQ.algebraic_field(sm.sqrt(2))),
             (sm.QQ.poly_ring(a,order='lex'),sm.QQ.poly_ring(a,order='grlex'))][case]
    p,q=(sm.Poly(x+1,x,domain=k) for k in domains)
    with pytest.raises(sm.UnificationFailed):
        getattr(sm,op)(p,q)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('extension', [sm.sqrt(2),True,(sm.sqrt(2),sm.I)])
def test_explicit_algebraic_extension_selection_is_joint(op, extension):
    import sympy as sy
    x=sm.Symbol('x')
    c=sm.sqrt(2)+sm.I if isinstance(extension,tuple) else sm.sqrt(2)
    f,g=(x+c)*(x+1),x+c
    se=tuple(sm.to_sympy(v) for v in extension) if isinstance(extension,tuple) else extension if extension is True else sm.to_sympy(extension)
    compare(getattr(sm,op)(f,g,x,extension=extension,polys=True),getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(g),sm.to_sympy(x),extension=se,polys=True))


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('coefficient', ['function','radical','imaginary','inexact'])
def test_opaque_historical_coefficients_still_compute(op, coefficient):
    import sympy as sy
    x,a=sm.symbols('x a')
    c={'function':sm.sin(a),'radical':sm.sqrt(2),'imaginary':sm.I,'inexact':0.5*a}[coefficient]
    f,g=2*c*(x+1),c*(x+1)
    actual=getattr(sm,op)(f,g,x)
    expected=getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(g),sm.to_sympy(x))
    compare(actual,expected)


def test_no_sympy_runtime_for_typed_global_division():
    code='''
import importlib.abc,sys,pickle
class Block(importlib.abc.MetaPathFinder):
 def find_spec(self,name,path=None,target=None):
  if name=='sympy' or name.startswith('sympy.'):raise ImportError('blocked')
sys.meta_path.insert(0,Block())
import symox as s
x,y=s.symbols('x y');a=s.Dummy('a',positive=True)
for k in [s.ZZ,s.QQ,s.GF(5),s.QQ.algebraic_field(s.sqrt(2)),s.QQ.frac_field(a),s.QQ.poly_ring(a)]:
 p=s.Poly((x+1)*(x+y),x,y,domain=k);q=s.Poly(x+1,x,y,domain=k)
 u,r=s.div(p,q,polys=True);assert u.domain==r.domain==k
 assert s.exquo(p,q).domain==k
 assert s.rem(p,q)==r and s.quo(p,q)==u
 assert pickle.loads(pickle.dumps(u)).domain==k
for fn in [s.div,s.rem,s.quo,s.exquo]:
 try:fn(x+1,x,x,wrong=True)
 except s.OptionError as e:assert type(pickle.loads(pickle.dumps(e))) is s.OptionError
 else:raise AssertionError('unknown option ignored')
assert not any(n=='sympy' or n.startswith('sympy.') for n in sys.modules)
'''
    result=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert result.returncode==0,result.stderr


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('prime', [2,5,7,101])
def test_modulus_option_replaces_inferred_coefficient_domain(op, prime):
    import sympy as sy
    x,y=sm.symbols('x y')
    f,g=(x+1)*(x+y),x+1
    compare(getattr(sm,op)(f,g,x,y,modulus=prime,polys=True),getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(g),sm.to_sympy(x),sm.to_sympy(y),modulus=prime,polys=True))
    # A divisor which vanishes in the selected field is zero, even when its
    # input expression was not the zero polynomial over ZZ.
    with pytest.raises(ZeroDivisionError):
        getattr(sm,op)(f,prime*g,x,y,modulus=prime)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('domain', ['QQ','ZZ','GF(5)','QQ(a,b)','ZZ[a,b]'])
def test_legacy_domain_strings_are_resolved_once(op, domain):
    import sympy as sy
    x,a,b=sm.symbols('x a b')
    f,g=((x+a)*(x+1),x+1) if 'a' in domain else ((x+2)*(x+1),x+1)
    compare(getattr(sm,op)(f,g,x,domain=domain,polys=True),getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(g),sm.to_sympy(x),domain=domain,polys=True))


@pytest.mark.parametrize('op', OPS)
def test_dense_division_cannot_be_replaced_by_sparse_reduction(op):
    import sympy as sy
    x,y=sm.symbols('x y')
    f,g=x*y+x,y
    for k,ok in [(sm.QQ,sy.QQ),(sm.GF(5),sy.GF(5)),(sm.QQ.algebraic_field(sm.sqrt(2)),sy.QQ.algebraic_field(sy.sqrt(2)))]:
        if op == 'exquo':
            with pytest.raises(sm.ExactQuotientFailed):
                sm.exquo(f,g,x,y,domain=k)
        else:
            compare(getattr(sm,op)(f,g,x,y,domain=k,polys=True),getattr(sy,op)(sm.to_sympy(f),sm.to_sympy(g),sm.to_sympy(x),sm.to_sympy(y),domain=ok,polys=True))


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('kind', ['GF','AN','QQRF'])
def test_mixed_poly_expression_joint_inference_and_explicit_override(op, kind):
    import sympy as sy
    x=sm.Symbol('x')
    k,ok,a,b=selected(kind)
    c=sm.sqrt(2) if kind == 'AN' else a if kind == 'QQRF' else sm.Integer(2)
    f=sm.Poly((x+c)*(x+1),x,domain=k)
    g=x+c
    of=sy.Poly(sm.to_sympy(f.as_expr()),sm.to_sympy(x),domain=ok)
    # A mixed expression/Poly call re-polifies both expressions, as SymPy
    # does; a supplied domain preserves that exact coefficient descriptor.
    try:
        expected=getattr(sy,op)(of,sm.to_sympy(g))
    except sy.polys.polyerrors.ExactQuotientFailed:
        with pytest.raises(sm.ExactQuotientFailed):
            getattr(sm,op)(f,g)
    else:
        compare(getattr(sm,op)(f,g),expected)
    result=getattr(sm,op)(f,g,domain=k)
    for p in result if op == 'div' else (result,):
        assert p.domain==k and p.gens==(x,)


@pytest.mark.parametrize('op', OPS)
def test_unimplemented_valid_options_are_rejected_explicitly(op):
    x=sm.Symbol('x')
    for option in ['greedy','symmetric']:
        with pytest.raises(NotImplementedError):
            getattr(sm,op)(x+1,x,x,**{option:True})


@pytest.mark.parametrize('op', OPS)
def test_generator_order_options_follow_reference(op):
    # explicit generators ignore wrt/sort, whose values are still validated
    x=sm.Symbol('x')
    for option in ['sort','wrt']:
        with pytest.raises(sm.OptionError,match="invalid argument for '%s' option"%option):
            getattr(sm,op)(x+1,x,x,**{option:True})
    if op!='exquo':
        assert getattr(sm,op)(x+1,x,x,wrt=x)==getattr(sm,op)(x+1,x,x)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('target', ['ZZ','GFsame','GFother','QQ'])
def test_finite_poly_domain_conversion_uses_typed_residues(op, target):
    import sympy as sy
    x=sm.Symbol('x');sx=sm.to_sympy(x)
    p=sm.Poly((3*x+1)*(x+1),x,domain=sm.GF(5))
    q=sm.Poly(x+1,x,domain=sm.GF(5))
    sp=sy.Poly((3*sx+1)*(sx+1),sx,domain=sy.GF(5))
    sq=sy.Poly(sx+1,sx,domain=sy.GF(5))
    k={'ZZ':sm.ZZ,'QQ':sm.QQ,'GFsame':sm.GF(5),'GFother':sm.GF(7)}[target]
    ok={'ZZ':sy.ZZ,'QQ':sy.QQ,'GFsame':sy.GF(5),'GFother':sy.GF(7)}[target]
    if target=='QQ':
        with pytest.raises(sy.polys.polyerrors.CoercionFailed):getattr(sy,op)(sp,sq,domain=ok)
        with pytest.raises(sm.CoercionFailed):getattr(sm,op)(p,q,domain=k)
    else:
        compare(getattr(sm,op)(p,q,domain=k),getattr(sy,op)(sp,sq,domain=ok))
        if target.startswith('GF'):
            prime=5 if target=='GFsame' else 7
            compare(getattr(sm,op)(p,q,modulus=prime),getattr(sy,op)(sp,sq,modulus=prime))


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('kind', ['GF','AN','QQRF','QQring'])
@pytest.mark.parametrize('reorder', [False,True])
def test_explicit_existing_generators_preserve_typed_coefficient_fields(op, kind, reorder):
    import sympy as sy
    x,y=sm.symbols('x y')
    k,ok,a,b=selected(kind,'grevlex')
    c=a if kind.startswith('QQ') else sm.sqrt(2) if kind=='AN' else sm.Integer(3)
    f=sm.Poly((x+c)*(y+1),x,y,domain=k);g=sm.Poly(x+c,x,y,domain=k)
    sf=sy.Poly(sm.to_sympy(f.as_expr()),sm.to_sympy(x),sm.to_sympy(y),domain=ok)
    sg=sy.Poly(sm.to_sympy(g.as_expr()),sm.to_sympy(x),sm.to_sympy(y),domain=ok)
    generators=(y,x) if reorder else (x,y)
    actual=getattr(sm,op)(f,g,*generators)
    expected=getattr(sy,op)(sf,sg,*(sm.to_sympy(v) for v in generators))
    compare(actual,expected)
    for p in actual if op=='div' else (actual,):
        assert p.domain==k and p.gens==generators


@pytest.mark.parametrize('op', OPS)
def test_legacy_finite_polys_with_opaque_coefficients_cannot_claim_scalar_field(op):
    x,a=sm.symbols('x a')
    for coefficient in [sm.sqrt(2),sm.sin(a)]:
        # These legacy constructors predate typed coefficient domains. Global
        # division must reject them rather than compute in an undeclared
        # extension or label symbolic coefficients as elements of GF(5).
        p=sm.Poly(coefficient*x+1,x,modulus=5)
        with pytest.raises(sm.UnificationFailed,match='scalar coefficients'):
            getattr(sm,op)(p,p)
        with pytest.raises(ValueError,match='scalar coefficients'):
            sm.Poly(x+1,x,domain=sm.GF(5)).div(p)


@pytest.mark.parametrize('op', OPS)
@pytest.mark.parametrize('kind', ['ZZ','QQ','GF','AN','QQRF','ZZRF','QQring','ZZring'])
@pytest.mark.parametrize('generators', ['implicit','same','reordered','reframed'])
def test_extension_true_preserves_existing_domains_before_joint_inference(op, kind, generators):
    import sympy as sy
    x,y,z=sm.symbols('x y z')
    k,ok,a,b=selected(kind,'grevlex')
    f,g=sm.Poly(x+y+1,x,y,domain=k),sm.Poly(2*x+2*y+2,x,y,domain=k)
    sf=sy.Poly(sm.to_sympy(x+y+1),sm.to_sympy(x),sm.to_sympy(y),domain=ok)
    sg=sy.Poly(sm.to_sympy(2*x+2*y+2),sm.to_sympy(x),sm.to_sympy(y),domain=ok)
    args={'implicit':(),'same':(x,y),'reordered':(y,x),'reframed':(z,)}[generators]
    actual=getattr(sm,op)(f,g,*args,extension=True)
    expected=getattr(sy,op)(sf,sg,*(sm.to_sympy(v) for v in args),extension=True)
    compare(actual,expected)
    if generators!='reframed' and kind in ['GF','AN','QQRF','ZZRF','QQring']:
        for value in actual if op=='div' else (actual,):
            assert value.domain==k
            assert all(c in value.domain for c in value.coeffs())


@pytest.mark.parametrize('op', OPS)
def test_extension_true_preserves_algebraic_alias_metadata_and_float_arithmetic(op):
    import sympy as sy
    x=sm.Symbol('x')
    k=sm.QQ.algebraic_field(sm.sqrt(2),alias='alpha')
    p,q=sm.Poly(x+1,x,domain=k),sm.Poly(2*x+2,x,domain=k)
    actual=getattr(sm,op)(p,q,extension=True)
    for value in actual if op=='div' else (actual,):
        assert value.domain==k and value.domain.ext.alias==sm.Symbol('alpha')
    p,q=sm.Poly(0.5*x+0.5,x),sm.Poly(x+1,x,domain=sm.QQ)
    sx=sm.to_sympy(x)
    expected=getattr(sy,op)(sy.Poly(0.5*sx+0.5,sx),sy.Poly(sx+1,sx,domain=sy.QQ),extension=True)
    compare(getattr(sm,op)(p,q,extension=True),expected)
