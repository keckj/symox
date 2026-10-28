"""CC quotient components certified from exact stored dyades, independently."""
from fractions import Fraction as F
import random
import pytest
import symox as sx
import sympy as sp

PRECISIONS=(1,53,139,140,211)


def float_dyad(man,exp,p):
    return sx.Float((int(man<0),abs(man),exp,abs(man).bit_length()),precision=p)


def dyad(value):
    if value.is_Float:
        sign,man,exp,_=value._mpf_
        return F((-1 if sign else 1)*man)*(F(2)**exp)
    assert value.is_Rational
    return F(int(value.p),int(value.q))


def parts(value):
    e=value.as_expr()
    return dyad(sx.re(e)),dyad(sx.im(e))


def rounded(q,p):
    if not q:return F(0)
    a,b=abs(q.numerator),q.denominator
    e=a.bit_length()-b.bit_length()
    if (a < b*(1<<e)) if e>=0 else (a*(1<<-e)<b):e-=1
    spacing=F(2)**(e-p+1);scaled=abs(q)/spacing
    n,r=divmod(scaled.numerator,scaled.denominator)
    if 2*r>scaled.denominator or 2*r==scaled.denominator and n%2:n+=1
    return (-n if q<0 else n)*spacing


def quotient(a,b,p):
    ar,ai=parts(a);br,bi=parts(b);d=br*br+bi*bi
    assert d
    return rounded((ar*br+ai*bi)/d,p),rounded((ai*br-ar*bi)/d,p)


def certify(k,a,b):
    expected=quotient(a,b,k.precision)
    result=a/b
    assert parts(result)==expected
    assert result.parent()==k
    assert all(f._prec==k.precision for f in result.as_expr().atoms(sx.Float))
    return result


@pytest.mark.parametrize('p',PRECISIONS)
@pytest.mark.parametrize('seed',range(18))
def test_full_complex_division_has_correctly_rounded_components(p,seed):
    rng=random.Random(195793+seed);k=sx.ComplexField(p)
    def value():
        return k.convert(sx.Rational(rng.randrange(-53,54),rng.randrange(1,18))+
                         sx.I*sx.Rational(rng.randrange(-53,54),rng.randrange(1,18)))
    a=value();b=value()
    if not b:b=k.one
    certify(k,a,b)
    certify(k,k.one,b)
    assert parts(b**-1)==quotient(k.one,b,p)


@pytest.mark.parametrize('p',PRECISIONS)
@pytest.mark.parametrize('kind',('real','imaginary','zero','identical','conjugate','cancellation','near_pole'))
def test_structural_zeros_cancellation_and_poles_preserve_bits(p,kind):
    k=sx.ComplexField(p)
    tiny=sx.Rational(1,2**(p+80));near=sx.Rational(2**p-1,2**p)
    if kind=='real':a,b=3,7
    elif kind=='imaginary':a,b=3*sx.I,7*sx.I
    elif kind=='zero':a,b=0,3+7*sx.I
    elif kind=='identical':a=b=1+near*sx.I
    elif kind=='conjugate':a,b=1+sx.I,1-sx.I
    elif kind=='cancellation':a,b=1+near*sx.I,1-near*sx.I
    else:a,b=1+sx.I,tiny+2*tiny*sx.I
    certify(k,k.convert(a),k.convert(b))


@pytest.mark.parametrize('p',(1,2,3,53,139,140,211))
def test_midpoint_ties_and_binade_transitions(p):
    k=sx.ComplexField(p)
    for a,b in [(1+2*sx.I,2+2*sx.I),(-1-2*sx.I,2+2*sx.I),
                (3+sx.I,2-2*sx.I),(1+sx.I,1-sx.I)]:
        certify(k,k.convert(a),k.convert(b))


@pytest.mark.parametrize('p',PRECISIONS)
@pytest.mark.parametrize('exp',(-10000,10000))
def test_extreme_magnitudes_with_exact_fraction_oracle(p,exp):
    k=sx.ComplexField(p)
    a=k.convert(sx.Float((0,3,exp,2),precision=p)+sx.I*sx.Float((1,5,exp-17,3),precision=p))
    b=k.convert(sx.Float((0,7,exp+11,3),precision=p)+sx.I*sx.Float((0,11,exp-13,4),precision=p))
    certify(k,a,b)


@pytest.mark.parametrize('p',PRECISIONS)
@pytest.mark.parametrize('exp',(10**12,-10**12,2**190,-2**190))
def test_unbounded_exponent_gaps_never_require_dense_shifts(p,exp):
    # An independent symbolic identity certifies these huge exponents; never
    # construct Fraction(2)**exp to materialize an enormous integer.
    k=sx.ComplexField(p);u=sx.Float((0,1,exp,1),precision=p)
    a=k.convert(u+sx.I*u);b=k.convert(u-sx.I*u)
    assert parts(a/b)==(F(0),F(1))
    assert parts(a/a)==(F(1),F(0))
    c=k.convert(1+sx.I*u)
    result=c/c
    assert parts(result)==(F(1),F(0)) and result.parent()==k


@pytest.mark.parametrize('p',PRECISIONS)
def test_same_selected_p_independent_of_presentation_precision(p):
    k=sx.ComplexField(p)
    a=k.convert(sx.Rational(2,7)+sx.I*sx.Rational(1,11))
    expected=quotient(k.one,a,p)
    assert parts(k.one/a)==expected
    # Test the reciprocal and negative-power protocol separately.
    assert parts(a**-1)==expected
    assert k.precision==p


@pytest.mark.parametrize('p',(140,211))
def test_first_incorrect_reciprocal_from_linear_system_regression(p):
    k=sx.ComplexField(p)
    if p==211:
        b=k.convert(sx.Rational(2,7)+sx.I*sx.Rational(1,11))
    else:
        re=float_dyad(338082524870514170850472812864291864830345,-139,p)
        im=float_dyad(453049261919023271490509285968143601185637,-138,p)
        b=k.convert(re+sx.I*im)
    result=certify(k,k.one,b)
    if p==140:
        assert parts(result)[1]==F(-940993504466598865261755653150685443041553,2**140)
    else:
        assert parts(result)[0]==F(1307450619184860710792925795146412643264894257215424150334213737,2**208)


@pytest.mark.parametrize('p',PRECISIONS)
def test_zero_division_error_and_nonfinite_conversion_are_not_hidden(p):
    k=sx.ComplexField(p)
    errors=[]
    for operation in (lambda:k.one/k.zero,lambda:k.zero**-1):
        with pytest.raises(Exception) as caught:operation()
        errors.append((type(caught.value),str(caught.value)))
    assert errors[0]==errors[1]
    assert 'division by zero' in errors[0][1]
    for value in (sx.nan,sx.oo,sx.zoo):
        with pytest.raises(sx.CoercionFailed):k.convert(value)


def test_custom_division_and_attribute_protocols_are_preserved(monkeypatch):
    cls=sx._core.DomainElement;native=cls.__truediv__;calls=[]
    def replacement(self,other):
        calls.append('division');return native(self,other)
    monkeypatch.setattr(cls,'__truediv__',replacement)
    k=sx.ComplexField(140);a=k.convert(2+3*sx.I);b=k.convert(7-11*sx.I)
    assert parts(a/b)==quotient(a,b,140) and calls==['division']
    def getter(self):
        calls.append('matrix');return False
    monkeypatch.setattr(cls,'is_Matrix',property(getter),raising=False)
    assert parts(a/b)==quotient(a,b,140) and calls[-2:]==['division','matrix']


@pytest.mark.parametrize('p',(53,140,211))
def test_parameter_domains_keep_selected_domain_and_complete_arithmetic(p):
    t=sx.Symbol('t');ground=sx.ComplexField(p)
    for k in (ground.poly_ring(t),ground.frac_field(t)):
        a=k.convert((1+sx.I)*t+2);one=k.one
        assert (a/one).parent()==k
        assert sx.simplify((a/one).as_expr()-a.as_expr())==0
        if k.is_FractionField:
            assert sx.simplify((a/a).as_expr()-1)==0 and (a/a).parent()==k


def test_reference_incorrect_multiplication_is_not_adopted():
    k=sx.ComplexField(140)
    a=k.convert(float_dyad(821110927807060936496507450470401678339335,-139,140)
                +sx.I*float_dyad(-704743296318480644540792222617112268510991,-139,140))
    b=k.convert(float_dyad(1393796574908163946345982392040522594123775,-140,140)
                +sx.I*float_dyad(-20115188238709136211484894206477219480047,-279,140))
    ar,ai=parts(a);br,bi=parts(b)
    assert parts(a*b)==(rounded(ar*br-ai*bi,140),rounded(ar*bi+ai*br,140))


@pytest.mark.parametrize('p',PRECISIONS)
@pytest.mark.parametrize('seed',range(12))
@pytest.mark.parametrize('op',('+','-','*'))
def test_other_primitives_use_exact_component_rounding(p,seed,op):
    rng=random.Random(991+seed);k=sx.ComplexField(p)
    def value():
        return k.convert(sx.Rational(rng.randrange(-15,16),rng.randrange(1,14))+
                         sx.I*sx.Rational(rng.randrange(-15,16),rng.randrange(1,14)))
    a=value();b=value();ar,ai=parts(a);br,bi=parts(b)
    if op=='+':exact=(ar+br,ai+bi);result=a+b
    elif op=='-':exact=(ar-br,ai-bi);result=a-b
    else:exact=(ar*br-ai*bi,ar*bi+ai*br);result=a*b
    assert parts(result)==tuple(rounded(q,p) for q in exact)
    assert result.parent()==k


@pytest.mark.parametrize('p',PRECISIONS)
def test_scalar_step_rounding_is_not_fused_across_operators(p):
    k=sx.ComplexField(p)
    a=k.convert(sx.Rational(3,7)+sx.I*sx.Rational(5,11))
    b=k.convert(sx.Rational(13,17)-sx.I*sx.Rational(19,23))
    c=k.convert(sx.Rational(29,31)+sx.I*sx.Rational(37,41))
    br,bi=parts(b);cr,ci=parts(c)
    product=(rounded(br*cr-bi*ci,p),rounded(br*ci+bi*cr,p))
    ar,ai=parts(a)
    expected=(rounded(ar-product[0],p),rounded(ai-product[1],p))
    assert parts(a-b*c)==expected
