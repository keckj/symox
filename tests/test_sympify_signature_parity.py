"""Sympify options are executed, not merely accepted by **kwargs."""
from decimal import Decimal
import inspect
import pytest
import sympy as sp
import symox as sx


def equivalent(a,b):
    if isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b)
        for x,y in zip(a,b):equivalent(x,y)
    else:
        assert (a.to_sympy() if hasattr(a,'to_sympy') else a)==b


@pytest.mark.parametrize('source',['x^2+1','x+x','0.1234567890123456789','0.125',
    'sin(x)**2+cos(x)**2','1/2','f(x)','2*(x+1)'])
@pytest.mark.parametrize('options',[{}, {'rational':True},{'evaluate':False},{'convert_xor':False}])
def test_text_options_match_reference(source,options):
    try:expected=sp.sympify(a=source,**options)
    except TypeError:
        with pytest.raises(TypeError):sx.sympify(a=source,**options)
    else:equivalent(sx.sympify(a=source,**options),expected)


def test_local_bindings_and_container_converter_rules():
    y=sx.Symbol('Y');Y=sp.Symbol('Y')
    for a in [[1,'x'],(1,'x'),{'x'}, {'x':'x'},range(0,5,2),Decimal('1.234567890123456789')]:
        equivalent(sx.sympify(a,locals={'x':y}),sp.sympify(a,locals={'x':Y}))
    equivalent(sx.S('x/2',locals={'x':y},rational=True),sp.S('x/2',locals={'x':Y},rational=True))
    assert str(inspect.signature(sx.sympify))==str(inspect.signature(sp.sympify))


@pytest.mark.parametrize('a',[None,[1,'x'],'x+1'])
def test_strict_rejects_parsing_and_nonconvertible_lists(a):
    with pytest.raises(sx.SympifyError):sx.sympify(a,strict=True)
    with pytest.raises(sp.SympifyError):sp.sympify(a,strict=True)
    if a is None:assert sx.sympify(a) is None


def test_user_hooks_and_registered_converters():
    class FloatLike:
        def __float__(self):return 1.25
    class IntLike:
        def __int__(self):return 3
    for strict in [False,True]:
        for value in [FloatLike(),IntLike()]:
            equivalent(sx.sympify(value,strict=strict),sp.sympify(value,strict=strict))
    class Hook:
        def _sympy_(self):return sx.Rational(1,3)
    assert sx.sympify(Hook())==sx.Rational(1,3)
    from symox.core.sympify import converter
    class Container:pass
    converter[Container]=lambda value:sx.Integer(7)
    try:assert sx.sympify(Container())==7
    finally:converter.pop(Container)


def test_native_structural_objects_never_use_interchange_hooks():
    from symox.geometry import Point
    for value in [Point(1,2),sx.Function('f')]:
        assert sx.sympify(value) is value
