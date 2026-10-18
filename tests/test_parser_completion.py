"""Token transforms and expression semantics against the SymPy reference."""
import ast
import builtins
import tokenize
from io import StringIO
import pytest
import symox as sx
import sympy as sp
from symox.parsing import sympy_parser as ours
from sympy.parsing import sympy_parser as reference


CASES=[('x+1',()),('x+y*z',()),('1/3',()),('sqrt(2)',()),('sin(x)**2+cos(x)**2',()),
       ('0.2[1]',()),('0.[3]',()),('2.1[6]',()),('x!',()),('x!!',()),
       ('lambda x: x**2',()),('f(x,y)',()),('x^2',('convert_xor',)),
       ('2x',('implicit_multiplication',)),('2(x+1)',('implicit_multiplication',)),
       ('sin x',('implicit_application',)),('sin 2x',('implicit_multiplication_application',)),
       ('sin^2 x',('convert_xor','implicit_multiplication_application')),
       ('xyz',('split_symbols','implicit_multiplication')),('theta',('split_symbols','implicit_multiplication')),
       ('x=1',('convert_equals_signs',)),('sin**2(x)',('function_exponentiation',)),
       ('0.125',('rationalize',)),('3j',()),('2**3**2',()),('α+β',())]


@pytest.mark.parametrize('source,extra',CASES)
def test_parse_reference(source,extra):
    a=ours.standard_transformations+tuple(getattr(ours,n) for n in extra)
    b=reference.standard_transformations+tuple(getattr(reference,n) for n in extra)
    got=ours.parse_expr(source,transformations=a)
    want=reference.parse_expr(source,transformations=b)
    assert got.to_sympy()==want


@pytest.mark.parametrize('source',['sin^2 x+3xyz','(x+1)(x-1)','.3x','10sin(x)','tan theta','x^3','x=1'])
def test_all_transformations(source):
    assert ours.parse_expr(source,transformations='all').to_sympy()==reference.parse_expr(source,transformations='all')


@pytest.mark.parametrize('source',['x+1','1+x','1+2','x*x','2**3','sin(0)','cos(0)','sqrt(4)','1==1','1<x<2','x/2','-x/2'])
def test_evaluate_false_tree(source):
    got=ours.parse_expr(source,evaluate=False)
    want=reference.parse_expr(source,evaluate=False)
    assert got.class_name==type(want).__name__
    assert tuple(str(a) for a in got.args)==tuple(str(a) for a in want.args)
    if got.is_Function:assert not got.is_UndefinedFunction


@pytest.mark.parametrize('name',['lambda_notation','auto_symbol','repeated_decimals','auto_number','factorial_notation',
                                 'convert_xor','implicit_application','implicit_multiplication','convert_equals_signs',
                                 'function_exponentiation','rationalize','implicit_multiplication_application','split_symbols'])
def test_transformed_code(name):
    source={'lambda_notation':'lambda x: x+1','factorial_notation':'x!!','repeated_decimals':'0.2[1]',
            'implicit_multiplication':'2x','implicit_application':'sin x','function_exponentiation':'sin**2(x)',
            'convert_equals_signs':'x=1','rationalize':'Float("0.125")'}.get(name,'x^2+2.5')
    def tokens():return [(t.type,t.string) for t in tokenize.generate_tokens(StringIO(source).readline)]
    assert getattr(ours,name)(tokens(),{},vars(sx).copy())==getattr(reference,name)(tokens(),{},vars(sp).copy())


def test_namespaces_custom_split_and_selectors():
    x=sx.Symbol('x',positive=True)
    assert ours.parse_expr('x+1',local_dict={'x':x})==x+1
    local={'x':''};assert ours.parse_expr('x+1',local_dict=local)==sx.Symbol('x')+1
    assert local=={'x':''}
    assert ours.parse_expr('sin(x)',local_dict={'sin':lambda a:a*a})==sx.Symbol('x')**2
    transform=ours.standard_transformations+(ours.split_symbols_custom(lambda n:n!='xy'),ours.implicit_multiplication)
    assert ours.parse_expr('xy+yz',transformations=transform)==sx.Symbol('xy')+sx.Symbol('y')*sx.Symbol('z')
    assert ours.T[:5]==ours.standard_transformations
    assert len(ours.T[0,2:4,-1])==4
    assert str(ours.T)==ours.transformations
    assert ours.eval_expr(ours.stringify_expr('x+1',{},vars(sx),ours.standard_transformations),{},vars(sx))==sx.Symbol('x')+1
    assert isinstance(ours.evaluateFalse('x+1'),ast.Expression)
    group=ours.ParenthesisGroup([(tokenize.NAME,'x')])
    applied=ours.AppliedFunction((tokenize.NAME,'f'),group)
    assert applied.expand()==[(tokenize.NAME,'f'),(tokenize.NAME,'x')]
    assert applied[0]==(tokenize.NAME,'f')


@pytest.mark.parametrize('kwargs',[{'local_dict':[]},{'global_dict':[]},{'local_dict':{'':1}},{'transformations':'unknown'}])
def test_invalid_options(kwargs):
    with pytest.raises((ValueError,TypeError)):ours.parse_expr('x+1',**kwargs)


def test_no_sympy_dependency(monkeypatch):
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name=='sympy' or name.startswith('sympy.'):raise AssertionError('imported SymPy')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    assert ours.parse_expr('2x+0.[3]',transformations='all')==2*sx.Symbol('x')+sx.Rational(1,3)
    assert ours.parse_expr('sin(0)',evaluate=False).args==(sx.Integer(0),)


def test_held_undefined_function_with_builtin_name():
    expression=ours.parse_expr('sin(0)',local_dict={'sin':sx.Function('sin')},evaluate=False)
    assert expression.is_UndefinedFunction
    assert expression.args==(sx.Integer(0),)
