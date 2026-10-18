"""Component evaluation, typed wildcard reconstruction and legacy interfaces."""
import itertools
import pickle
import warnings

import pytest
import sympy as sp
import symox as sx
from symox.tensor import tensor as tt
from symox.combinatorics import Permutation


def setup():
    L=tt.TensorIndexType('L',dim=2);i,j,k,l=tt.tensor_indices('i j k l',L)
    A,B=tt.tensor_heads('A B',[L]);M,N=tt.tensor_heads('M N',[L,L])
    S=tt.TensorHead('S',[L,L],tt.TensorSymmetry.fully_symmetric(2))
    F=tt.TensorHead('F',[L,L],tt.TensorSymmetry.fully_symmetric(-2))
    g=L.metric;delta=L.delta
    data={A(i):[1,2],B(-i):[3,4],M(i,j):[[1,2],[3,4]],
          N(-i,-j):[[2,0],[1,2]],S(i,j):[[1,2],[2,4]],F(i,j):[[0,2],[-2,0]],
          L:[[2,0],[0,3]]}
    return locals()


@pytest.mark.parametrize('source,axes',[
    ('A(i)',None),('A(-i)',None),('A(i)',('-i',)),
    ('M(i,j)',('j','i')),('M(i,-j)',None),('M(-i,-j)',None),
    ('A(i)*B(-i)',None),('M(i,-j)*A(j)',None),
    ('M(i,j)*N(-j,k)',None),('M(i,j)+M(j,i)',None),
    ('S(i,j)*B(-j)',None),('F(i,j)*B(-j)',None),
    ('F(j,i)',('i','j')),('g(i,j)*B(-j)',None),
    ('g(i,-i)',None),('delta(i,-i)',None),('M(i,-i)',None),
    ('A(i)*A(-i)',None),('(A(i)+2*A(i))*B(-i)',None),
])
def test_component_arrays_match_reference_including_variance_and_axis_order(source,axes):
    env=setup();expression=eval(source,env);data=env['data']
    wanted=None if axes is None else [eval(axis,env) for axis in axes]
    actual=expression.replace_with_arrays(data,wanted)
    foreign={key.to_sympy():sx.Array(value).to_sympy() for key,value in data.items()}
    # SymPy requires explicit metric and delta tensor entries in addition to
    # the type's covariant metric. Native evaluation derives those entries.
    foreign[env['g'](-env['i'],-env['j']).to_sympy()]=sp.Array([[2,0],[0,3]])
    foreign[env['delta'](env['i'],-env['j']).to_sympy()]=sp.Array([[1,0],[0,1]])
    expected=expression.to_sympy().replace_with_arrays(foreign,None if wanted is None else [i.to_sympy() for i in wanted])
    if getattr(expected,'shape',None)==():expected=expected[()]
    assert actual.to_sympy()==expected
    canonical=tt.canon_bp(expression)
    if isinstance(canonical,tt.TensExpr):
        # Align to the same output axes before comparing a canonical form.
        free=expression.get_free_indices() if wanted is None else wanted
        assert canonical.replace_with_arrays(data,free)==actual
    else:assert canonical==actual


@pytest.mark.parametrize('mapping',[{'i':0},{'i':1},{'j':0},{'i':1,'j':0}])
def test_selected_components_remove_only_the_selected_free_axes(mapping):
    env=setup();M,i,j=env['M'],env['i'],env['j']
    positions={env[name]:value for name,value in mapping.items()}
    value=tt.TensorElement(M(i,j),positions)
    actual=value.replace_with_arrays(env['data'])
    expected=M(i,j).to_sympy().replace_with_arrays({M(i,j).to_sympy():sp.Array([[1,2],[3,4]])})
    selectors=tuple(mapping.get(name,slice(None)) for name in ('i','j'))
    assert actual.to_sympy()==expected[selectors]
    assert value.get_free_indices()==[index for index in (i,j) if index not in positions]
    assert pickle.loads(pickle.dumps(value))==value
    assert value.to_sympy().get_free_indices()==[index.to_sympy() for index in value.get_free_indices()]


@pytest.mark.parametrize('dimension',[2,3])
def test_riemann_component_certificates_for_canonicalization_and_bianchi(dimension):
    L=tt.TensorIndexType('L',dim=dimension);i,j,k,l=tt.tensor_indices('i j k l',L)
    R=tt.TensorHead('R',[L]*4,tt.TensorSymmetry.riemann())
    values=[int(a==c and b==d)-int(a==d and b==c)
            for a,b,c,d in itertools.product(range(dimension),repeat=4)]
    data={R(i,j,k,l):sx.Array(values,(dimension,)*4),L:sx.eye(dimension)}
    contraction=R(i,j,k,l)*R(-i,-j,-k,-l)
    assert contraction.replace_with_arrays(data)==2*dimension*(dimension-1)
    assert tt.canon_bp(contraction).replace_with_arrays(data)==2*dimension*(dimension-1)
    projected=tt.riemann_cyclic(R(i,j,k,l))
    assert projected.replace_with_arrays(data,[i,j,k,l])==sx.Array(values,(dimension,)*4)


def test_component_shape_metric_symmetry_and_bound_validation():
    env=setup();L,A,M,S,F,i,j=(env[name] for name in ['L','A','M','S','F','i','j'])
    with pytest.raises(ValueError):A(i).replace_with_arrays({A(i):[1,2,3]})
    with pytest.raises(ValueError):A(i).replace_with_arrays({})
    with pytest.raises(ValueError):A(-i).replace_with_arrays({A(i):[1,2]})
    with pytest.raises(ValueError):A(-i).replace_with_arrays({A(i):[1,2],L:[[1,2],[0,1]]})
    with pytest.raises(ValueError):S(i,j).replace_with_arrays({S(i,j):[[1,2],[3,4]]})
    with pytest.raises(ValueError):F(i,j).replace_with_arrays({F(i,j):[[1,2],[-2,0]]})
    with pytest.raises(ValueError):M(i,j).replace_with_arrays(env['data'],[i,i])
    with pytest.raises(ValueError):tt.TensorElement(A(i),{i:2})
    with pytest.raises(ValueError):tt.TensorElement(A(i),{i:-1})
    n=sx.Symbol('n',integer=True)
    with pytest.raises(NotImplementedError):tt.TensorElement(A(i),{i:n}).replace_with_arrays(env['data'])
    odd=tt.TensorHead('odd',[L],comm=1)
    with pytest.raises(NotImplementedError):odd(i).replace_with_arrays({odd(i):[1,2]})
    nonsym=tt.TensorIndexType('U',dim=2,metric_symmetry=0);p,q=tt.tensor_indices('p q',nonsym)
    with pytest.raises(NotImplementedError):tt.contract_metric(nonsym.metric(p,-p),nonsym.metric)


def test_large_sparse_components_validate_stored_symmetries_without_dense_expansion():
    L=tt.TensorIndexType('L',dim=10**6);i,j=tt.tensor_indices('i j',L)
    S=tt.TensorHead('S',[L,L],tt.TensorSymmetry.fully_symmetric(2))
    sparse=sx.ImmutableSparseNDimArray({(0,1):2,(1,0):2},(10**6,10**6))
    result=S(j,i).replace_with_arrays({S(i,j):sparse},[i,j])
    assert result.shape==(10**6,10**6) and result[0,1]==result[1,0]==2
    bad=sx.ImmutableSparseNDimArray({(0,1):2},(10**6,10**6))
    with pytest.raises(ValueError):S(i,j).replace_with_arrays({S(i,j):bad})
    A,B=tt.tensor_heads('A B',[L])
    vector=sx.ImmutableSparseNDimArray({8123:3},(10**6,))
    other=sx.ImmutableSparseNDimArray({8123:5},(10**6,))
    assert (A(i)*B(-i)).replace_with_arrays({A(i):vector,B(-i):other})==15
    with pytest.raises(ValueError):A(-i).replace_with_arrays({A(i):vector,L:sparse})
    E=tt.TensorIndexType('E',dim=2);k=tt.TensorIndex('k',E)
    matrix=tt.TensorHead('M',[E,L])
    wide=sx.ImmutableSparseNDimArray({(0,1):1},(2,10**6))
    with pytest.raises(ValueError):matrix(-k,i).replace_with_arrays({matrix(k,i):wide,E:sx.eye(2)})


def test_tensor_index_type_special_heads_and_automatic_indices():
    L=tt.TensorIndexType('L',dim=3);i=tt.TensorIndex(True,L);j=tt.TensorIndex(True,L)
    assert i!=j and -(-i)==i
    assert L.delta is L.delta and L.metric is L.metric and L.epsilon is L.epsilon
    assert L.epsilon.rank==3 and L.epsilon.symmetry==tt.TensorSymmetry.fully_symmetric(-3)
    assert L.get_kronecker_delta()==L.delta and L.get_epsilon()==L.epsilon
    assert tt.TensorIndexType('N',metric_symmetry=None).metric is None
    assert tt.TensorIndexType('D').epsilon is None
    with pytest.raises(ValueError):tt.TensorIndexType('')
    with pytest.raises(ValueError):tt.TensorIndexType('N',dim=-1)


def test_typed_index_wildcards_enforce_variance_and_consistent_bindings():
    env=setup();L,i,j=env['L'],env['i'],env['j']
    p=tt.WildTensorIndex('p',L);q=tt.WildTensorIndex('q',L,ignore_updown=True)
    assert p.matches(i)=={p:i} and p.matches(-i) is None
    assert q.matches(-i)=={q:-i}
    assert p.matches(j,{p:i}) is None
    assert (-p).matches(-i,{p:i})=={p:i,-p:-i}
    assert (-p).matches(-j,{p:i}) is None
    assert pickle.loads(pickle.dumps(p))==p
    foreign=tt.TensorIndex('z',tt.TensorIndexType('Other',dim=2))
    assert p.matches(foreign) is None


def test_wild_tensor_heads_and_factor_patterns_reconstruct_the_matched_expression():
    env=setup();L,A,B,M,i,j=(env[name] for name in ['L','A','B','M','i','j'])
    W=tt.WildTensorHead('W');U=tt.WildTensorHead('U',unordered_indices=True)
    p=tt.WildTensorIndex('p',L)
    for pattern,expression in [(W(),A(i)),(W(),M(i,j)),(W(i),A(i)),
                               (A(p),A(i)),(A(i)*W(),A(i)*B(j)),
                               (A(i)*W(),2*A(i)*B(j)),
                               (M(p,j)*B(-j),M(i,j)*B(-j))]:
        result=pattern.matches(expression)
        assert result is not None
        assert pattern.xreplace(result)==expression
        assert expression.match(pattern)==result
    assert W(i).matches(A(j)) is None
    assert W(i,j).matches(M(j,i)) is None
    result=U(i,j).matches(M(j,i))
    assert result is not None and result[U.component if hasattr(U,'component') else tt.WildTensorHead('U',[L,L],unordered_indices=True)].expr==M(j,i)
    assert U(i,j).xreplace(result)==M(j,i)
    assert isinstance(W(i),tt.WildTensor)
    assert tt.WildTensor(A,[i])==A(i)
    P,Q=tt.tensor_heads('P Q',[L],comm=2)
    pattern=W()*P(i);expression=Q(j)*P(i)
    result=pattern.matches(expression)
    assert result is not None and pattern.xreplace(result)==expression


def test_legacy_factories_deprecation_signals_and_postprocessors():
    env=setup();L,A,i,j=env['L'],env['A'],env['i'],env['j']
    with pytest.warns(DeprecationWarning):sym=tt.tensorsymmetry([1,1])
    assert sym==tt.TensorSymmetry.fully_symmetric(2)
    with pytest.warns(DeprecationWarning):anti=tt.tensorsymmetry([2])
    assert anti==tt.TensorSymmetry.fully_symmetric(-2)
    with pytest.warns(DeprecationWarning):R=tt.tensorsymmetry([2,2])
    assert R==tt.TensorSymmetry.riemann()
    with pytest.warns(DeprecationWarning):factory=tt.TensorType([L],tt.TensorSymmetry.no_symmetry(1))
    assert factory('V')==tt.TensorHead('V',[L])
    with pytest.warns(DeprecationWarning):head=tt.tensorhead('S',[L,L],[[1,1]])
    assert head.symmetry==sym
    for function in [tt.deprecate_call,tt.deprecate_data,tt.deprecate_fun_eval]:
        with pytest.warns(DeprecationWarning):function()
    with pytest.warns(DeprecationWarning):assert A(i)(j)==A(j)
    with pytest.warns(DeprecationWarning):assert A(i).fun_eval((i,j))==A(j)
    assert tt.tensor_mul(2,A(i))==2*A(i)
    assert sx.Add(A(i),A(i))==2*A(i)
    assert sx.Mul(2,A(i))==2*A(i)
    assert tt.TensAdd.fromiter([A(i),A(i)])==2*A(i)
    class Add:
        args=(A(i),A(i))
    class Mul:
        args=(2,A(i))
    assert tt.get_postprocessor(Add)(Add())==2*A(i)
    assert tt.get_postprocessor(Mul)(Mul())==2*A(i)
    assert tt.get_postprocessor(Add)(sx.Integer(2))==2
