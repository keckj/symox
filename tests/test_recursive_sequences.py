"""Recurrence sequences are evaluated from their defining finite history."""
import itertools
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.series.sequences import RecursiveSeq,SeqExpr,SeqExprOp,SeqBase
from sympy.series.sequences import RecursiveSeq as Reference


@pytest.mark.parametrize('formula,initial',[
    ('y(n-1)+y(n-2)',[0,1]),('2*y(n-1)+n',[1]),
    ('y(n-1)**2+1',[0]),('y(n-1)+a*y(n-3)',[1,2,3]),
    ('n*n',[]),
])
@pytest.mark.parametrize('start',[0,2,-3])
def test_recurrence_values_and_metadata_match_definition(formula,initial,start):
    n,a=sx.symbols('n a');y=sx.Function('y')
    N,A=sp.symbols('n a');Y=sp.Function('y')
    f=eval(formula,{'y':y,'n':n,'a':a})
    g=eval(formula,{'y':Y,'n':N,'a':A})
    actual=RecursiveSeq(f,y(n),n,initial,start=start)
    # SymPy 1.14 has an offset bug for nonzero starts. Compare its zero-origin
    # sequence after translating the explicit index in the defining recurrence.
    translated=g.subs(N,N+start).replace(lambda node:node.func==Y,
                                       lambda node:Y(node.args[0]-start))
    expected=Reference(translated,Y(N),N,initial,start=0)
    assert [value.to_sympy() for value in actual[:8]]==expected[:8]
    assert actual.degree==len(initial)
    assert actual.start==start and actual.stop==sx.oo
    assert list(itertools.islice(actual,8))==actual[:8]
    assert actual.coeff(start+3)==actual[3]
    assert actual.recurrence==sx.Eq(y(n),f)
    assert pickle.loads(pickle.dumps(actual))==actual
    assert actual.func(*actual.args)==actual
    assert actual.to_sympy().recurrence==sp.Eq(Y(N),g)


def test_symbolic_initial_values_parameters_scaling_and_operations():
    n,a,b=sx.symbols('n a b');y=sx.Function('y')
    value=RecursiveSeq(a*y(n-1)+n,y(n),n,[b])
    assert value.free_symbols=={a,b}
    assert value.subs({a:2,b:1})[:4]==[1,3,8,19]
    assert value.subs(n,99)==value
    changed=value.subs(a,n)
    assert changed.variables[0]!=n and n in changed.free_symbols
    assert all(sx.cancel(left-2*right)==0
               for left,right in zip(value.coeff_mul(2)[:5],value[:5]))
    assert value.coeff_mul(n)[:5]==[i*value[i] for i in range(5)]
    assert (value+sx.sequence(n))[:5]==[value[i]+i for i in range(5)]
    assert (value*sx.sequence(n))[:5]==[value[i]*i for i in range(5)]
    auto=RecursiveSeq(y(n-1)+y(n-2),y(n),n)
    assert len(auto.initial)==2 and auto.initial[0]!=auto.initial[1]
    assert auto[2]==sum(auto.initial)


def test_sequence_expression_bases_have_their_own_metadata_and_hierarchy():
    n,m,a=sx.symbols('n m a')
    expression=SeqExpr(sx.Tuple(1,2,3),(n,0,10))
    assert expression.gen==sx.Tuple(1,2,3)
    assert expression.interval==sx.Interval(0,10) and expression.length==11
    assert expression.variables==(n,)
    assert expression.to_sympy().gen==sp.Tuple(1,2,3)
    assert isinstance(sx.SeqFormula(n),SeqExpr) and isinstance(sx.SeqPer((1,2)),SeqExpr)
    left=sx.sequence(n*n,(n,0,10));right=sx.sequence((1,2),(m,5,9))
    op=SeqExprOp(left,right)
    assert isinstance(op,SeqBase) and op.start==5 and op.stop==9 and op.length==5
    assert set(op.gen)=={left.gen,right.gen}
    assert op.to_sympy().length==op.length.to_sympy()
    assert op.variables==(n,m)
    assert isinstance(sx.SeqAdd(left,right,evaluate=False),SeqExprOp)
    assert isinstance(sx.SeqMul(left,right,evaluate=False),SeqExprOp)
    with pytest.raises(NotImplementedError):expression.coeff(0)
    with pytest.raises(NotImplementedError):op.coeff(5)


def test_recurrence_argument_validation_and_evaluation_bound():
    n,a=sx.symbols('n a');y=sx.Function('y')
    for recurrence in [y(n),y(n+1),y(2*n),y(n-a),y(n-sx.Rational(1,2)),y(n-1,n-2)]:
        with pytest.raises(TypeError):RecursiveSeq(recurrence,y(n),n,[1])
    with pytest.raises(TypeError):RecursiveSeq(n,n,n)
    with pytest.raises(ValueError):RecursiveSeq(y(n-2),y(n),n,[1])
    with pytest.raises(ValueError):RecursiveSeq(y(n-65),y(n),n,[1]*65)
    value=RecursiveSeq(y(n-1)+1,y(n),n,[0])
    with pytest.raises(IndexError):value.coeff(-1)
    with pytest.raises(ValueError):value.coeff(100000)
    with pytest.raises(NotImplementedError):value.coeff(sx.Rational(1,2))


def test_recurrence_api_does_not_import_sympy():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import *
n=Symbol('n');y=Function('y')
f=RecursiveSeq(y(n-1)+y(n-2),y(n),n,[0,1])
assert f[:8]==[0,1,1,2,3,5,8,13]
assert isinstance(SeqFormula(n),SeqExpr)
'''
    subprocess.run([sys.executable,'-c',source],check=True)
