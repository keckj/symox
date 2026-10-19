"""Explicit coefficient-domain choices retain the reference auto defaults."""
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('operation',['div','quo','rem','exquo'])
@pytest.mark.parametrize('choice',['ZZ','QQ','ring','field_false','field_true',
                                   'domain_none','extension_true','radical','modulus'])
@pytest.mark.parametrize('auto',['absent',None,False,True])
def test_explicit_domain_defaults_and_auto_override(operation,choice,auto):
    def inputs(lib):
        x,a=lib.symbols('x a')
        options={'ZZ':dict(domain=lib.ZZ),'QQ':dict(domain=lib.QQ),
                 'ring':dict(domain=lib.ZZ.poly_ring(a)),
                 'field_false':dict(field=False),'field_true':dict(field=True),
                 'domain_none':dict(domain=None),'extension_true':dict(extension=True),
                 'radical':dict(extension=lib.sqrt(2)),'modulus':dict(modulus=5)}[choice]
        if auto!='absent':options['auto']=auto
        g=(a+1)*x+a if choice=='ring' else 2*x-4
        return x,x*x+1,g,options
    x,f,g,options=inputs(sp)
    try:
        expected=getattr(sp,operation)(f,g,x,polys=True,**options)
    except sp.ExactQuotientFailed:
        x,f,g,options=inputs(sx)
        with pytest.raises(sx.ExactQuotientFailed):
            getattr(sx,operation)(f,g,x,polys=True,**options)
        return
    x,f,g,options=inputs(sx)
    actual=getattr(sx,operation)(f,g,x,polys=True,**options)
    pairs=zip(actual,expected) if isinstance(expected,tuple) else [(actual,expected)]
    for a,e in pairs:
        assert sp.expand(a.as_expr().to_sympy()-e.as_expr())==0
        assert str(a.domain)==str(e.domain)
