"""Sequence terms, interval arithmetic, scope and fitted recurrences."""
import itertools
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
from symox.series.sequences import (SeqBase,SeqFormula,SeqPer,SeqAdd,SeqMul,
                                    EmptySequence,sequence)


def test_public_sequence_exports():
    assert sx.S.EmptySequence is EmptySequence
    for name in ['SeqFormula','SeqPer','SeqAdd','SeqMul','EmptySequence','sequence']:
        assert getattr(sx,name) is globals()[name]


def equivalent(actual,expected):
    assert sp.simplify(actual.to_sympy()-expected)==0


@pytest.mark.parametrize('formula',['n**2','2**n','(-1)**n','n*(n+1)/2','factorial(n)'])
@pytest.mark.parametrize('bounds',[(0,8),(2,10),(-4,4),(0,'oo'),('-oo',0)])
def test_formula_terms_and_metadata_match_reference(formula,bounds):
    n=sx.Symbol('n');N=sp.Symbol('n')
    start,stop=map(sx.sympify,bounds);a=sx.sympify(formula)
    # Factorials at negative integer positions are poles; compare native
    # values directly through the common infinity representation.
    native=SeqFormula(a,(n,start,stop));reference=sp.SeqFormula(sp.sympify(formula),(N,start.to_sympy(),stop.to_sympy()))
    assert native.variables==(n,) and native.gen==a
    assert native.free_symbols==set()
    assert native.length.to_sympy()==reference.length
    assert native.interval.to_sympy()==reference.interval
    assert [v.to_sympy() for v in native[:7]]==reference[:7]
    assert list(itertools.islice(native,7))==native[:7]
    assert pickle.loads(pickle.dumps(native))==native
    assert native.to_sympy()==reference


@pytest.mark.parametrize('period',[(1,2,3),(0,-1),(1,),('n','n**2','n**3')])
@pytest.mark.parametrize('bounds',[(0,12),(2,14),(-4,8),(0,'oo'),('-oo',0)])
def test_periodic_terms_and_reverse_infinite_direction(period,bounds):
    n=sx.Symbol('n');N=sp.Symbol('n')
    values=tuple(sx.sympify(value) for value in period)
    start,stop=map(sx.sympify,bounds)
    native=SeqPer(values,(n,start,stop));reference=sp.SeqPer(tuple(v.to_sympy() for v in values),(N,start.to_sympy(),stop.to_sympy()))
    assert native.period==len(period)
    assert [v.to_sympy() for v in native[:10]]==reference[:10]
    assert native.gen==native.periodical
    assert pickle.loads(pickle.dumps(native))==native


def test_default_indices_and_sequence_factory():
    n,a=sx.symbols('n a')
    assert isinstance(sequence(n*n),SeqFormula)
    assert sequence(n*n).variables==(n,)
    assert isinstance(sequence((1,2,3)),SeqPer)
    assert SeqFormula(2)[:4]==[2,2,2,2]
    assert SeqPer([[1,2],[3,4]])[:6]==[1,2,3,4,1,2]
    assert SeqPer(sx.Tuple(1,2))[:4]==[1,2,1,2]
    seq=sequence(n)
    assert sequence(seq) is seq
    with pytest.raises(ValueError):sequence(n*a)
    assert sequence(n*a,(n,0,3)).free_symbols=={a}
    assert sequence(n,(n,3,0)) is EmptySequence
    assert sequence((),(n,0,3)) is EmptySequence


def test_sequence_sort_keys_do_not_iterate_infinite_series():
    n=sx.Symbol('n');N=sp.Symbol('n')
    native=SeqFormula(n*n,(n,0,sx.oo));reference=sp.SeqFormula(N*N,(N,0,sp.oo))
    def normalized(value):
        if isinstance(value,tuple):return tuple(normalized(item) for item in value)
        return value.to_sympy() if isinstance(value,sx.Expr) else value
    assert normalized(sx.default_sort_key(native))==sp.default_sort_key(reference)


@pytest.mark.parametrize('left_kind',['formula','periodic','mixed'])
@pytest.mark.parametrize('right_kind',['formula','periodic','mixed'])
@pytest.mark.parametrize('operation',['add','mul','sub'])
def test_sequence_arithmetic_uses_actual_coefficients_and_intervals(left_kind,right_kind,operation):
    n,m,a=sx.symbols('n m a')
    left={'formula':SeqFormula(n*n+a,(n,0,12)),
        'periodic':SeqPer((n,n+1),(n,0,12)),
        'mixed':SeqAdd(SeqFormula(n,(n,0,12)),SeqPer((1,2),(n,0,12)))}[left_kind]
    right={'formula':SeqFormula(m**3+a,(m,2,11)),
        'periodic':SeqPer((m,2*m,3*m),(m,2,11)),
        'mixed':SeqMul(SeqFormula(m,(m,2,11)),SeqPer((3,4),(m,2,11)))}[right_kind]
    if operation=='add':result=left+right;function=lambda a,b:a+b
    elif operation=='mul':result=left*right;function=lambda a,b:a*b
    else:result=left-right;function=lambda a,b:a-b
    assert result.start==2 and result.stop==11
    expected_symbols={a} if left_kind=='formula' or right_kind=='formula' else set()
    if operation=='sub' and left_kind==right_kind=='formula':expected_symbols=set()
    assert result.free_symbols==expected_symbols
    for point in range(2,12):
        assert sx.cancel(result.coeff(point)-function(left.coeff(point),right.coeff(point)))==0
    assert pickle.loads(pickle.dumps(result))==result


def test_periodic_arithmetic_aligns_different_origins():
    n=sx.Symbol('n')
    left=SeqPer((1,2),(n,0,5));right=SeqPer((3,4,5),(n,1,5))
    # SymPy 1.14 combines these periods without accounting for their
    # different origins. Verify the defining termwise identity directly.
    assert (left+right)[:5]==[5,5,7,4,6]
    assert (left*right)[:5]==[6,4,10,3,8]
    assert (left+right).coeff(1)==left.coeff(1)+right.coeff(1)


def test_periodic_fusion_handles_opposite_infinite_directions():
    n=sx.Symbol('n')
    left=SeqPer((n,n*n,1),(n,-sx.oo,10))
    right=SeqPer((2*n,3),(n,0,sx.oo))
    for operation in [SeqAdd,SeqMul]:
        actual=operation(left,right)
        assert actual.start==0 and actual.stop==10
        for point in range(11):
            expected=(left.coeff(point)+right.coeff(point) if operation is SeqAdd else
                      left.coeff(point)*right.coeff(point))
            assert actual.coeff(point)==expected


def test_bound_index_scope_and_capture_avoidance():
    n,m,a=sx.symbols('n m a')
    seq=SeqFormula(n+a,(n,0,m))
    assert seq.free_symbols=={a,m}
    assert seq.subs({a:3,m:4}).coeff(2)==5
    assert seq.subs(n,100).formula==seq.formula
    replaced=seq.subs(a,n)
    assert replaced.variables[0]!=n and replaced.free_symbols=={n,m}
    assert replaced.subs(m,4).coeff(2)==n+2
    left=SeqFormula(n+m,(n,0,4));right=SeqFormula(m+n,(m,0,4))
    total=left+right
    assert total.variables[0] not in (n,m)
    assert total.free_symbols=={n,m}
    assert total.coeff(2)==n+m+4
    assert seq.xreplace({a:2,m:4}).coeff(3)==5
    assert seq.xreplace({seq:left}) is left


@pytest.mark.parametrize('formula,expected',[
    ('n**2',[3,-3,1]),('2**n',[2]),('3*5**n+12',[6,-5]),
    ('n**4+3*n**2',[5,-10,10,-5,1]),('(-1)**n',[-1]),
    ('2**n+3**n',[5,-6]),('fibonacci(n)',[1,1]),('lucas(n)',[1,1]),
])
def test_linear_recurrences_and_generating_functions(formula,expected):
    n,t=sx.symbols('n t');N,T=sp.symbols('n t')
    seq=SeqFormula(sx.sympify(formula),(n,0,sx.oo))
    coefficients,generating=seq.find_linear_recurrence(14,gfvar=t)
    assert coefficients==expected
    reference=sp.SeqFormula(sp.sympify(formula),(N,0,sp.oo)).find_linear_recurrence(14,gfvar=T)
    equivalent(generating,reference[1])
    # Check unseen terms independently of the fitting sample.
    for k in range(14,24):
        assert seq.coeff(k)==sum(value*seq.coeff(k-j-1) for j,value in enumerate(coefficients))


def test_periodic_recurrence_bounds_and_failed_search():
    n=sx.Symbol('n')
    seq=SeqPer((1,2,3),(n,0,sx.oo))
    coefficients=seq.find_linear_recurrence(12)
    assert coefficients==[0,0,1]
    assert SeqFormula(n*n).find_linear_recurrence(10,d=2)==[]
    assert SeqFormula(sx.factorial(n)).find_linear_recurrence(10,d=3)==[]
    assert SeqFormula(0).find_linear_recurrence(10,gfvar=sx.Symbol('t'))==([],None)
    for size in [-1,1000]:
        with pytest.raises(ValueError):seq.find_linear_recurrence(size)


def test_empty_sequences_and_unevaluated_operations():
    n=sx.Symbol('n');seq=SeqFormula(n,(n,0,4))
    assert list(EmptySequence)==[] and EmptySequence[:]==[]
    assert seq+EmptySequence is seq
    assert seq*EmptySequence is EmptySequence
    assert SeqAdd() is EmptySequence and SeqMul() is EmptySequence
    assert SeqAdd(seq,SeqFormula(n,(n,5,10))) is EmptySequence
    assert SeqAdd(seq,seq,evaluate=False).coeff(2)==4
    assert pickle.loads(pickle.dumps(EmptySequence)) is EmptySequence
    assert seq.coeff_mul(2)[:5]==[0,2,4,6,8]


def test_invalid_inputs_and_bounded_slices():
    n=sx.Symbol('n');seq=SeqFormula(n,(n,0,4))
    for limits in [(n,-sx.oo,sx.oo),(n,),(n,0,1,2),(2,0,4),(n,0,sx.Rational(1,2))]:
        with pytest.raises(ValueError):SeqFormula(n,limits)
    for index in [-1,5]:
        with pytest.raises(IndexError):seq[index]
    with pytest.raises(ValueError):SeqFormula(n)[:]
    with pytest.raises(ValueError):SeqFormula(n)[:100001]
    with pytest.raises(ValueError):SeqPer(list(range(4097)))
    with pytest.raises(ValueError):seq[::0]
    with pytest.raises(TypeError):seq*2
    with pytest.raises(TypeError):seq+2


def test_no_sympy_import():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import symbols
from symox import sequence,SeqPer
n,t=symbols('n t')
seq=sequence(n*n)
assert seq[:4]==[0,1,4,9]
assert seq.find_linear_recurrence(10)==[3,-3,1]
assert (seq+SeqPer((1,2),(n,0,100)))[:4]==[1,3,5,11]
'''
    subprocess.run([sys.executable,'-c',source],check=True)
