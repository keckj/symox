"""Einstein algebra is checked against independent tensor identities."""
import itertools
import pickle
import subprocess
import sys

import pytest
import sympy as sp
import symox as sx
import sympy.tensor.tensor as reference
from symox.tensor import tensor as tt


def environment():
    L=tt.TensorIndexType('L',dim=3)
    i,j,k,l,m,n=tt.tensor_indices('i j k l m n',L)
    A,B,C=tt.tensor_heads('A B C',[L])
    M,N=tt.tensor_heads('M N',[L,L])
    S=tt.TensorHead('S',[L,L],tt.TensorSymmetry.fully_symmetric(2))
    F=tt.TensorHead('F',[L,L],tt.TensorSymmetry.fully_symmetric(-2))
    R=tt.TensorHead('R',[L]*4,tt.TensorSymmetry.riemann())
    x=sx.Symbol('x')
    return locals()


@pytest.mark.parametrize('source',[
    'A(i)','M(j,i)','A(i)*B(-i)','M(i,j)*N(-j,-i)',
    'A(i)*B(-i)*C(j)','A(i)*B(-i)*A(j)*B(-j)',
    'M(i,-i)','M(-i,i)','M(i,-j)*N(j,-k)',
    '2*A(i)+3*A(i)','x*A(i)+A(i)*x','A(i)-A(i)',
    'M(i,j)+M(j,i)','(A(i)+B(i))*C(j)',
    'A(i)*(B(-i)+C(-i))','(A(i)+B(i))*(A(-i)+B(-i))',
])
def test_tensor_algebra_metadata_pickle_and_scalar_dispatch(source):
    env=environment();actual=eval(source,env)
    foreign={name:value.to_sympy() for name,value in env.items() if hasattr(value,'to_sympy')}
    expected=eval(source,foreign).expand()
    assert actual.to_sympy()==expected
    assert pickle.loads(pickle.dumps(actual))==actual
    if isinstance(actual,tt.TensExpr):
        assert actual.rank==expected.rank
        assert [index.to_sympy() for index in actual.get_indices()]==expected.get_indices()
        assert [index.to_sympy() for index in actual.get_free_indices()]==expected.get_free_indices()
        if not isinstance(actual,tt.TensAdd):
            assert {(index.to_sympy(),pos) for index,pos in actual.free}==set(expected.free)
            assert set(actual.dum)==set(expected.dum)
        assert actual.free_symbols=={sx.sympify(symbol) for symbol in expected.free_symbols}


@pytest.mark.parametrize('source',[
    'S(j,i)','F(j,i)','F(i,j)*F(-j,-i)',
    'S(i,j)*F(-i,-j)','M(i,j)*N(-j,-i)',
    'B(-j)*A(j)','M(i,-j)*N(j,-k)',
    'R(l,k,j,i)','R(j,i,k,l)','R(k,l,i,j)',
    'R(i,j,k,l)*R(-i,-j,-k,-l)',
    'S(i,j)*A(-j)*B(-i)','F(i,j)*A(-j)*B(-i)',
    'S(j,i)-S(i,j)','F(j,i)+F(i,j)',
])
def test_signed_orbit_canonicalization_matches_reference(source):
    env=environment();expression=eval(source,env)
    actual=tt.canon_bp(expression)
    expected=reference.canon_bp(expression.to_sympy())
    assert actual.to_sympy()==expected
    assert tt.canon_bp(actual)==actual
    assert sx.simplify(expression)==actual


@pytest.mark.parametrize('rank',[0,1,2,3,4,-1,-2,-3,-4])
def test_slot_symmetry_groups_have_the_same_signed_permutations(rank):
    from symox.combinatorics import PermutationGroup
    native=tt.TensorSymmetry.fully_symmetric(rank)
    foreign=reference.TensorSymmetry.fully_symmetric(rank)
    actual={tuple(p.array_form) for p in PermutationGroup(*native.generators).elements}
    from sympy.combinatorics import PermutationGroup as Group
    expected={tuple(p.array_form) for p in Group(*foreign.generators).generate_schreier_sims()}
    assert native.rank==abs(rank) and actual==expected
    assert pickle.loads(pickle.dumps(native))==native


@pytest.mark.parametrize('blocks',[(1,1),(2,-2),(3,1),(-2,-2),(2,0,1)])
def test_direct_product_symmetries_match_reference(blocks):
    actual=tt.TensorSymmetry.direct_product(*blocks)
    expected=reference.TensorSymmetry.direct_product(*blocks)
    assert actual.to_sympy()==expected


@pytest.mark.parametrize('symmetry',[1,-1])
@pytest.mark.parametrize('source',[
    'g(i,-i)','g(-i,i)','g(i,j)*A(-i)','g(i,j)*A(-j)',
    'g(-i,-j)*A(i)','g(-i,-j)*A(j)',
    'g(i,j)*A(-i)*B(-j)','A(-i)*g(i,j)*B(-j)',
    'g(i,j)*g(-j,k)*M(-i,-k)',
])
def test_metric_contractions_match_reference(symmetry,source):
    L=tt.TensorIndexType('L',dim=4,metric_symmetry=symmetry)
    i,j,k=tt.tensor_indices('i j k',L);g=L.metric
    A,B=tt.tensor_heads('A B',[L]);M=tt.TensorHead('M',[L,L])
    expression=eval(source,locals())
    result=tt.contract_metric(expression,g)
    expected=expression.to_sympy().contract_metric(g.to_sympy())
    if symmetry==-1 and source=='A(-i)*g(i,j)*B(-j)':
        # SymPy 1.14's antisymmetric contraction sign depends on the metric's
        # factor position. Commuting scalar tensors cannot change that sign.
        expected=(g(i,j)*A(-i)*B(-j)).to_sympy().contract_metric(g.to_sympy())
    assert reference.canon_bp(result.to_sympy()-expected)==0


def test_commutation_groups_and_antisymmetric_zero_detection():
    env=environment();L,i,j=env['L'],env['i'],env['j']
    A,B=tt.tensor_heads('A B',[L],comm=1)
    assert tt.canon_bp(A(i)*A(-i))==0
    assert tt.canon_bp(B(i)*A(j))==-tt.canon_bp(A(j)*B(i))
    P,Q=tt.tensor_heads('P Q',[L],comm=2)
    assert tt.canon_bp(P(i)*Q(j))!=tt.canon_bp(Q(j)*P(i))
    assert P.commutes_with(Q) is None and A.commutes_with(B)==1


@pytest.mark.parametrize('metric',[1,-1,None])
@pytest.mark.parametrize('heads',[('A','A'),('A','F'),('F','F')])
@pytest.mark.parametrize('up',[(True,True),(True,False),(False,True),(False,False)])
def test_canonicalization_respects_metric_variance_and_absent_metrics(metric,heads,up):
    L=tt.TensorIndexType('L',dim=4,metric_symmetry=metric);i,j=tt.tensor_indices('i j',L)
    A=tt.TensorHead('A',[L,L]);F=tt.TensorHead('F',[L,L],tt.TensorSymmetry.fully_symmetric(-2))
    first,second=[{'A':A,'F':F}[name] for name in heads]
    p=i if up[0] else -i;q=j if up[1] else -j
    expression=first(p,q)*second(-p,-q)
    result=tt.canon_bp(expression)
    assert reference.canon_bp(result.to_sympy()-expression.to_sympy())==0
    assert tt.canon_bp(result)==result


def test_dummy_renaming_cannot_capture_free_indices():
    env=environment();L,A,B,i,j=map(env.get,['L','A','B','i','j'])
    collision=tt.TensorIndex('L_0',L)
    expression=A(i)*B(-i)*A(collision)
    assert expression.get_free_indices()==[collision]
    assert all(index.name!='L_0' for index in tt.get_dummy_indices(expression))
    product=(A(i)*B(-i))*(A(i)*B(-i))
    assert product.rank==0 and len(product.dum)==2
    value=A(i)*B(-i)*A(j)
    changed=tt.substitute_indices(value,(j,collision))
    assert changed.get_free_indices()==[collision] and len(changed.dum)==1
    assert changed.equals(A(i)*B(-i)*A(collision))
    assert (sx.Symbol('x')*A(i)).subs(A,B)==sx.Symbol('x')*B(i)


def test_index_structure_helpers_signed_permutations_and_lines():
    env=environment();L,i,j,k,A,B,M,N=(env[name] for name in ['L','i','j','k','A','B','M','N'])
    expression=2*M(i,-j)*N(j,-k)
    structure=tt.get_index_structure(expression)
    permutation,dummies,metrics=structure.indices_canon_args()
    foreign=reference.get_index_structure(expression.to_sympy())
    ref_perm,ref_dummies,ref_metrics=foreign.indices_canon_args()
    assert permutation.array_form==ref_perm.array_form
    assert dummies==ref_dummies and metrics==ref_metrics
    assert tt.perm2tensor(expression,permutation)==expression
    from symox.combinatorics import Permutation
    negative=permutation.array_form;negative[-2:]=reversed(negative[-2:])
    assert tt.perm2tensor(expression,Permutation(negative))==-expression
    assert tt.get_coeff(expression)==2
    assert tt.get_indices(3)==tt.get_free_indices(3)==tt.get_dummy_indices(3)==()
    assert tt.get_index_structure(3).indices==[]
    assert tt.get_lines(expression,L)==reference.get_lines(expression.to_sympy(),L.to_sympy())
    trace=M(i,-j)*N(j,-i)
    assert tt.get_lines(trace,L)==reference.get_lines(trace.to_sympy(),L.to_sympy())
    args=tt.components_canon_args([M,M,N])
    assert [row[2:] for row in args]==[(2,0),(1,0)]
    with pytest.raises(ValueError):tt.get_coeff(A(i)+B(i))


def test_riemann_cyclic_identity_and_replacement_formula():
    env=environment();i,j,k,l,R=(env[name] for name in ['i','j','k','l','R'])
    expression=R(i,j,k,l)
    actual=tt.riemann_cyclic_replace(expression)
    expected=reference.riemann_cyclic_replace(expression.to_sympy())
    assert reference.canon_bp(actual.to_sympy()-expected)==0
    assert tt.riemann_cyclic(R(i,j,k,l)+R(i,k,l,j)+R(i,l,j,k))==0


def test_invalid_index_structures_and_explicit_canonicalization_bound(monkeypatch):
    env=environment();L,A,M,i,j=(env[name] for name in ['L','A','M','i','j'])
    other=tt.TensorIndexType('Other',dim=3);p=tt.TensorIndex('p',other)
    with pytest.raises(ValueError):A(p)
    with pytest.raises(ValueError):M(i,i)
    with pytest.raises(ValueError):A(i)*A(i)
    with pytest.raises(ValueError):tt.TensMul(A(i),A(-i),A(i))
    with pytest.raises(ValueError):A(i)+A(j)
    with pytest.raises(ValueError):A(i)+1
    with pytest.raises(ValueError):A(i)/A(j)
    with pytest.raises(ValueError):tt.substitute_indices(A(i),(i,p))
    with pytest.raises(ValueError):tt.TensorHead('Bad',[L,L],tt.TensorSymmetry.no_symmetry(1))
    with pytest.raises(ValueError):tt.TensorHead('Bad',[L,other],tt.TensorSymmetry.fully_symmetric(2))
    from symox.tensor import _canonical
    monkeypatch.setattr(_canonical,'MAX_ORBIT',1)
    with pytest.raises(NotImplementedError):tt.canon_bp(env['S'](j,i))


def test_indexed_tensor_algebra_runs_without_sympy_imports():
    source='''
import sys
class Block:
    def find_spec(self,name,*args):
        if name=='sympy' or name.startswith('sympy.'):
            raise AssertionError('SymPy import attempted')
sys.meta_path.insert(0,Block())
from symox import Symbol
from symox.tensor.tensor import *
L=TensorIndexType('L',dim=2);i,j=tensor_indices('i j',L)
A=TensorHead('A',[L]);S=TensorHead('S',[L,L],TensorSymmetry.fully_symmetric(2))
F=TensorHead('F',[L,L],TensorSymmetry.fully_symmetric(-2))
assert canon_bp(S(i,j)*F(-i,-j))==0
assert contract_metric(L.metric(i,j)*A(-j),L.metric)==A(i)
assert (A(i)*A(-i)).replace_with_arrays({A(i):[1,2],L:[[1,0],[0,1]]})==5
assert (Symbol('x')*A(i)).coeff==Symbol('x')
'''
    subprocess.run([sys.executable,'-c',source],check=True)
