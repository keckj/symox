"""Polynomial boundary integration and numerical source parser regressions."""
import ctypes
import importlib
import itertools
import math
import random
import subprocess
import pytest
import symox as S
import sympy as R


@pytest.mark.parametrize('vertices',[
    [(0,0),(0,1),(1,1),(1,0)],[(0,0),(2,0),(3,2),(1,4),(-1,2)],
    [(0,0),(0,3),(1,3),(1,1),(3,1),(3,0)],[(0,3),(5,3),(1,1)],
])
def test_planar_polynomial_integrals_against_reference(vertices):
    from symox.integrals.intpoly import polytope_integrate
    from sympy.integrals.intpoly import polytope_integrate as reference
    x,y=S.symbols('x y');rx,ry=R.symbols('x y');native,ref=S.Polygon(*vertices),R.Polygon(*vertices)
    for i,j in itertools.product(range(4),repeat=2):
        expression=x**i*y**j
        assert S.simplify(polytope_integrate(native,expression)-S.sympify(str(reference(ref,rx**i*ry**j))))==0
    assert polytope_integrate(native,[1,x,y,x*y],max_degree=2)=={e:polytope_integrate(native,e) for e in [S.Integer(1),x,y,x*y]}


def test_spatial_polynomial_moments_and_orientation():
    from symox.integrals.intpoly import polytope_integrate
    cube=[[(0,0,0),(0,0,5),(0,5,0),(0,5,5),(5,0,0),(5,0,5),(5,5,0),(5,5,5)],[2,6,7,3],[3,7,5,1],[7,6,4,5],[1,5,4,0],[3,1,0,2],[0,4,6,2]]
    x,y,z=S.symbols('x y z')
    for i,j,k in itertools.product(range(3),repeat=3):
        expected=-S.Rational(5**(i+j+k+3),(i+1)*(j+1)*(k+1))
        assert polytope_integrate(cube,x**i*y**j*z**k)==expected
    assert polytope_integrate(cube,1+x+y+z)==-S.Rational(2125,2)
    reversed_cube=[cube[0],*[list(reversed(face)) for face in cube[1:]]]
    assert polytope_integrate(reversed_cube,1)==125


def test_polytope_helpers_and_hyperplane_representation():
    native=importlib.import_module('symox.integrals.intpoly');reference=importlib.import_module('sympy.integrals.intpoly')
    x,y=S.symbols('x y');triangle=S.Polygon((0,3),(5,3),(1,1))
    parameters=native.hyperplane_parameters(triangle)
    assert parameters==[((0,1),3),((1,-2),-1),((-2,-1),-3)]
    assert native.integration_reduction(triangle.sides,0,*parameters[0],1,(x,y),0)==5
    assert native.left_integral2D(3,0,triangle.sides,triangle.sides[0].p1,1,(x,y))==5
    assert native.distance_to_side((0,0,0),[(0,0,1),(0,1,0)],(1,0,0))==-S.sqrt(2)/2
    assert native.intersection(((-1,0),0),((0,1),1),'plane2D')==(0,1)
    assert native.decompose(x*x+x*y+x+y)=={1:x+y,2:x*x+x*y}
    assert native.gradient_terms(2)==[[1,0,0,0],[y,0,1,0],[y*y,0,2,0],[x,1,0,0],[x*y,1,1,0],[x*x,2,0,0]]
    hrep=[((-1,0),0),((0,1),1),((1,0),1),((0,-1),0)]
    assert native.polytope_integrate(hrep,x*y)==S.Rational(1,4)


def test_resolved_meijer_adapters_and_convergence():
    from symox.integrals.meijerint import meijerint_definite,meijerint_indefinite,meijerint_inversion
    x,t,a=S.symbols('x t a')
    assert meijerint_definite(S.exp(-x*x),x,-S.oo,S.oo)==(S.sqrt(S.pi),S.true)
    assert meijerint_indefinite(S.sin(x),x)==-S.cos(x)
    assert meijerint_inversion(1/x,x,t)==S.Heaviside(t)
    assert meijerint_definite(1/x,x,-1,1) is None
    assert meijerint_indefinite(S.Function('unsolved')(x),x) is None


def test_numeric_program_parsing_runs_python_and_compiled_c(tmp_path):
    from symox.parsing.sym_expr import SymPyExpression
    programs=[SymPyExpression('double f(double x){ double r=1; for(int i=0;i<3;i++){r+=x;} return r;}','c'),
              SymPyExpression('real function f(x)\nreal :: x\nreal :: r\nr = 1 + 3*x\nf = r\nend function','f')]
    for index,program in enumerate(programs):
        namespace={};exec('\n'.join(program.convert_to_python()),namespace)
        assert namespace['f'](2)==7
        source=tmp_path/('parsed%d.c'%index);binary=tmp_path/('parsed%d.so'%index)
        source.write_text('#include <math.h>\n'+'\n'.join(program.convert_to_c()))
        result=subprocess.run(['cc','-shared','-fPIC',str(source),'-o',str(binary),'-lm'],capture_output=True,text=True)
        assert result.returncode==0,result.stderr
        function=ctypes.CDLL(str(binary)).f;function.argtypes=[ctypes.c_double];function.restype=ctypes.c_double
        assert function(2)==7
    with pytest.raises(SyntaxError):SymPyExpression('int f(int x){if(x){return 1;}return 0;}','c')


def test_custom_binary_relations_and_fact_metadata():
    from symox.assumptions import BinaryRelation,AppliedBinaryRelation,Predicate
    from symox.assumptions.ask import get_known_facts_dict,get_all_known_facts
    class Successor(BinaryRelation):name='remaining_successor';is_reflexive=False
    relation=Successor()
    @relation.register(S.Number,S.Number)
    def numeric(left,right,assumptions):return S.ask(S.Q.zero(left-right-1),assumptions)
    assert isinstance(relation,Predicate)
    assert isinstance(relation(2,1),AppliedBinaryRelation)
    assert S.ask(relation(2,1)) is True and S.ask(relation(1,2)) is False
    x=S.Symbol('x');assert S.ask(relation(x,x)) is False
    facts=get_known_facts_dict()
    assert S.Q.rational in facts[S.Q.integer][0] and S.Q.irrational in facts[S.Q.integer][1]
    assert get_all_known_facts()
