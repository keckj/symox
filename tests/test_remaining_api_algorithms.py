"""Exact optimization and domain matrix calculations, compared with SymPy."""
import importlib
from fractions import Fraction
import random
import pytest
import symox as S
import sympy as R


@pytest.mark.parametrize('name',['bernoulli_poly','bernoulli_c_poly','euler_poly','genocchi_poly','andre_poly'])
@pytest.mark.parametrize('n',range(9))
def test_appell_sequences(name,n):
    ours=importlib.import_module('symox.polys.appellseqs');ref=importlib.import_module('sympy.polys.appellseqs')
    x=S.Symbol('x');rx=R.Symbol('x')
    result=getattr(ours,name)(n,x)
    assert R.expand(R.sympify(str(result))-getattr(ref,name)(n,rx))==0
    assert getattr(ours,name)(n,x,polys=True).as_expr()==result

@pytest.mark.parametrize('case',[
 ([1,2],[[-1,-1]],[-3],None,None,None),
 ([-1,-2],[[1,1]],[3],None,None,None),
 ([1,2],None,None,[[1,1]],[3],None),
 ([2,-1],[[-1,-1],[1,2]],[-2,8],None,None,[(None,None),(0,3)]),
 ([-1],[[-1]],[-2],None,None,[(0,5)]),
 ([1,1],[[1,0],[0,1]],[3,4],[[1,-1]],[1],[(-2,None),(-3,None)]),
])
def test_exact_two_phase_simplex(case):
    from symox.solvers.simplex import linprog
    from sympy.solvers.simplex import linprog as ref
    actual,point=linprog(*case)
    c,A,b,Aeq,beq,bounds=case;variables=R.symbols('x:'+str(len(c)));conditions=[]
    for matrix,values,equality in ((A,b,False),(Aeq,beq,True)):
        for row,value in zip(matrix or [],values or []):
            lhs=sum(a*x for a,x in zip(row,variables));conditions.append(R.Eq(lhs,value) if equality else lhs<=value)
    boundlist=[(0,None)]*len(c) if bounds is None else bounds
    for x,(lower,upper) in zip(variables,boundlist):
        if lower is not None:conditions.append(x>=lower)
        if upper is not None:conditions.append(x<=upper)
    from sympy.solvers.simplex import lpmin as reference_min
    expected,_=reference_min(sum(a*x for a,x in zip(c,variables)),conditions)
    assert R.sympify(str(actual))==expected
    c,A,b,Aeq,beq,bounds=case
    if A:
        for row,limit in zip(A,b):assert bool(sum(a*v for a,v in zip(row,point))<=limit)
    if Aeq:
        for row,limit in zip(Aeq,beq):assert sum(a*v for a,v in zip(row,point))==limit


def test_simplex_symbolic_free_variables_and_error_cases():
    from symox.solvers.simplex import lpmin,lpmax,linprog,InfeasibleLPError,UnboundedLPError
    x,y=S.symbols('x y')
    assert lpmin(x+y,[x>=-5,y>=2])==(-3,{x:-5,y:2})
    assert lpmax(x+2*y,[x+y<=3,x>=0,y>=0])==(6,{x:0,y:3})
    assert lpmin(2*x+4,[S.Eq(x,3)])==(10,{x:3})
    with pytest.raises(InfeasibleLPError):lpmin(x,[x>=2,x<=1])
    with pytest.raises(UnboundedLPError):lpmin(-x,[x>=0])
    with pytest.raises(InfeasibleLPError):linprog([1],bounds=[(2,1)])
    with pytest.raises(UnboundedLPError):linprog([-1])


def test_simplex_random_bounded_exact_problems():
    from symox.solvers.simplex import linprog
    from sympy.solvers.simplex import linprog as ref
    rng=random.Random(142)
    for _ in range(25):
        n=3;c=[rng.randint(-4,4) for _ in range(n)]
        A=[[rng.randint(-2,4) for _ in range(n)] for _ in range(3)];b=[rng.randint(3,12) for _ in A]
        ours,_=linprog(c,A,b,bounds=(0,5));expected,_=ref(c,A,b,bounds=(0,5))
        assert R.sympify(str(ours))==expected


@pytest.mark.parametrize('fmt',['ddm','sdm'])
def test_domain_storage_algebra(fmt):
    name=fmt.upper();C=getattr(importlib.import_module('symox.polys.matrices.'+fmt),name)
    rows=[[1,2],[3,5]] if fmt=='ddm' else {0:{0:1,1:2},1:{0:3,1:5}}
    a=C(rows,(2,2),S.QQ)
    matrix=lambda m:m._dm().to_Matrix()
    assert matrix(a.inv())==S.Matrix([[-5,2],[3,-1]])
    assert matrix(a.matmul(a.inv()))==S.eye(2)
    assert a.det()==S.QQ(-1)
    rref,pivots=a.rref();assert matrix(rref)==S.eye(2) and pivots==[0,1]
    lower,upper,swaps=a.lu();assert not swaps and matrix(lower.matmul(upper))==matrix(a)
    assert matrix(a.to_ddm())==matrix(a.to_sdm())==matrix(a.to_dfm())
    copied=a.copy();copied.setitem(0,0,7);assert a.getitem(0,0)==S.QQ(1)
    b=C.from_list([[1],[0]],(2,1),S.QQ)
    assert matrix(a.lu_solve(b))==S.Matrix([-5,3])
    singular=C.from_list([[1,2],[2,4]],(2,2),S.QQ)
    kernel,free=singular.nullspace();assert free==[1]
    assert matrix(singular) * matrix(kernel).T==S.zeros(2,1)


def test_low_level_dense_and_sparse_algorithms():
    from symox.polys.matrices.dense import ddm_irref,ddm_idet,ddm_iinv,ddm_berk,ddm_irref_den
    from symox.polys.matrices.sdm import sdm_irref,sdm_matmul,sdm_berk
    a=[[S.QQ(1),S.QQ(2)],[S.QQ(3),S.QQ(5)]]
    assert ddm_idet([row[:] for row in a],S.QQ)==S.QQ(-1)
    out=[[S.QQ.zero]*2 for _ in range(2)];ddm_iinv(out,a,S.QQ)
    assert out==[[S.QQ(-5),S.QQ(2)],[S.QQ(3),S.QQ(-1)]]
    reduced=[row[:] for row in a];assert ddm_irref(reduced)==[0,1]
    assert reduced==[[S.QQ.one,S.QQ.zero],[S.QQ.zero,S.QQ.one]]
    z=[[S.ZZ(2),S.ZZ(1)],[S.ZZ(0),S.ZZ(3)]];den,pivots=ddm_irref_den(z,S.ZZ)
    assert pivots==[0,1] and z[0][0]==z[1][1]==den
    sparse={0:{0:S.QQ(1),1:S.QQ(2)},1:{0:S.QQ(3),1:S.QQ(5)}}
    rr,pivots,nonzero=sdm_irref(sparse)
    assert pivots==[0,1] and rr=={0:{0:S.QQ.one},1:{1:S.QQ.one}} and nonzero=={}
    assert sdm_matmul(rr,sparse,S.QQ,2,2)==sparse
    assert [v[0] for v in ddm_berk(a,S.QQ)]==[S.QQ(1),S.QQ(-6),S.QQ(-1)]
    assert sdm_berk(sparse,2,S.QQ)=={0:S.QQ(1),1:S.QQ(-6),2:S.QQ(-1)}


def test_propositional_knowledge_base_and_three_valued_evaluation():
    from symox.logic import gateinputcount
    from symox.logic.inference import PropKB,pl_true,entails,valid,literal_symbol
    p,q=S.symbols('p q');kb=PropKB(S.Implies(p,q));kb.tell(p)
    assert kb.ask(q);kb.retract(p);assert not kb.ask(q)
    assert valid(S.Or(p,S.Not(p)))
    assert entails(q,[p,S.Implies(p,q)])
    assert pl_true(S.Or(p,q),{p:True}) is True
    assert pl_true(S.And(p,q),{p:False}) is False
    assert pl_true(S.And(p,q),{p:True}) is None
    assert literal_symbol(S.Not(p))==p
    assert gateinputcount(S.Or(p,S.And(p,q)))==4


def test_context_local_polynomial_configuration():
    from symox.polys.polyconfig import setup,query,using
    setup('GROEBNER','buchberger')
    with using(GROEBNER='f5b'):
        assert query('groebner')=='f5b'
        with using(GROEBNER='buchberger'):assert query('GROEBNER')=='buchberger'
        assert query('GROEBNER')=='f5b'
    assert query('GROEBNER')=='buchberger'


def test_integer_storage_nullspace_preserves_its_ring():
    from symox.polys.matrices.ddm import DDM
    from symox.polys.matrices.sdm import SDM
    for cls in (DDM,SDM):
        source=cls.from_list([[S.ZZ(2),S.ZZ(3),S.ZZ(5)]],(1,3),S.ZZ)
        kernel,free=source.nullspace()
        assert kernel.domain==S.ZZ and len(free)==2
        assert source.matmul(kernel.transpose()).is_zero_matrix


def test_module_morphism_integral_kernels_images_and_minpoly():
    # The reference API (SymPy 1.14): elements from ZZ columns, ``coeffs``,
    # ``matrix(modulus)`` and ``kernel(modulus)``.
    from symox.polys.numberfields.modules import PowerBasis,ModuleHomomorphism,InnerEndomorphism,find_min_poly
    from symox.polys.matrices import DomainMatrix
    def col(*c):return DomainMatrix([[S.ZZ(v)] for v in c],(len(c),1),S.ZZ)
    x=S.Symbol('module_x');basis=PowerBasis(S.Poly(x**3-2,x))
    f=ModuleHomomorphism(basis,basis,lambda a:basis(col(2*a.coeffs[0]+3*a.coeffs[1]+5*a.coeffs[2],0,0)))
    assert f.matrix().to_Matrix()==S.Matrix([[2,3,5],[0,0,0],[0,0,0]])
    # Over ZZ the reference kernel converts a rational nullspace to ZZ and
    # raises CoercionFailed here.
    with pytest.raises(S.CoercionFailed):f.kernel()
    alpha=basis(col(0,1,0));multiplication=InnerEndomorphism(basis,alpha)
    assert multiplication.matrix().to_Matrix()==S.Matrix([[0,0,2],[1,0,0],[0,1,0]])
    assert find_min_poly(alpha,S.QQ,x).as_expr()==x**3-2
    residues=f.kernel(modulus=2)
    assert residues.n==2 and residues.matrix.to_Matrix()==S.Matrix([[1,0],[0,1],[0,1]])
    assert residues.represent(basis(col(1,0,0))).to_Matrix()==S.Matrix([1,0])
