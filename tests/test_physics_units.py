"""Dimensional arithmetic, exact conversions and the reference unit catalog."""
import importlib
import math
import subprocess
import sys
import pytest
import symox as S
from symox.physics import units as U
from symox.physics.units._core import _native,_analyze
from symox.physics.units._catalog import QUANTITIES,ALIASES
import sympy as R
import sympy.physics.units as ref


def same(a,b):assert S.simplify(a-S.sympify(b))==0


@pytest.mark.parametrize('name',list(QUANTITIES))
def test_catalog_conversion_against_reference(name):
    ours=U._quantities[name]
    theirs=next(getattr(ref,n) for n,(kind,key) in ALIASES.items() if kind=='quantity' and key==name)
    bases=[U.meter,U.gram,U.second,U.ampere,U.kelvin,U.mole,U.candela,U.bit,U.steradian,U.angular_mil,S.Integer(1)]
    rbases=[ref.meter,ref.gram,ref.second,ref.ampere,ref.kelvin,ref.mole,ref.candela]
    factor,dimensions=_analyze(S.sympify(ours),U.SI)
    value=U.convert_to(ours,bases)
    ratio=value/S.Mul(*(S.sympify(q)**dimensions.get(d,0) for q,d in zip(bases,['length','mass','time','current','temperature','amount_of_substance','luminous_intensity','information','steradian','angular_mil',None])))
    assert float(ratio.evalf())==pytest.approx(float(theirs.scale_factor.evalf()),rel=2e-12)
    assert U.SI.get_dimension_system().get_dimensional_dependencies(ours.dimension)=={U.Dimension(n):S.sympify(p) for n,p in dimensions.items()}


@pytest.mark.parametrize('expr,target,result',[
 (U.mile,U.meter,S.Rational(201168,125)*U.meter),
 (U.kilometer+U.mile,U.meter,S.Rational(326168,125)*U.meter),
 (U.joule,[U.kilogram,U.meter,U.second],U.kilogram*U.meter**2/U.second**2),
 (U.speed_of_light,[U.meter,U.second],299792458*U.meter/U.second),
 (S.pi*U.radian,U.degree,180*U.degree),
 (U.inch**2,U.meter,S.Rational(16129,25000000)*U.meter**2),
 (U.newton,U.second,S.sympify(U.newton)),
 (U.kilometer*(U.kilometer+U.mile),U.meter,2609344*U.meter**2),
 (2*U.hour,U.minute,120*U.minute),
 (3,U.meter,3),
])
def test_exact_conversions(expr,target,result):same(U.convert_to(expr,target),result)


def test_prefixes_and_custom_quantity():
    from symox.physics.units.prefixes import kilo,milli,mega,kibi,Prefix,prefix_unit
    assert kilo*milli==1 and kilo*kilo==mega
    assert kibi.scale_factor==1024
    same(kilo*U.meter,1000*U.meter)
    smoot=U.Quantity('smoot','smoot');smoot.set_global_relative_scale_factor(S.Rational(67,1),U.inch)
    same(U.convert_to(smoot,U.meter),S.Rational(8509,5000)*U.meter)
    assert len(prefix_unit(smoot,{'k':kilo,'m':milli}))==2
    same(U.quantity_simplify(U.centimeter/U.millimeter),10)
    same(U.quantity_simplify(5*U.kilometer/U.hour,across_dimensions=True,unit_system='SI'),S.Rational(25,18)*U.meter/U.second)
    same(U.quantity_simplify(U.ampere*U.ohm,across_dimensions=True,unit_system='SI'),U.volt)


def test_dimensions_and_rejection():
    dims=U.SI.get_dimension_system()
    assert dims.equivalent_dims(U.force,U.mass*U.length/U.time**2)
    assert dims.is_dimensionless(U.length/U.length)
    assert (U.length+U.length)==U.length
    assert U.length.has_integer_powers(dims)
    assert not (U.length**S.Rational(1,2)).has_integer_powers(dims)
    assert dims.get_dimensional_dependencies(S.sqrt(U.length))=={U.length:S.Rational(1,2)}
    assert U.SI.is_consistent
    U.check_dimensions(U.newton*U.meter+U.joule)
    with pytest.raises(ValueError):U.check_dimensions(U.meter+U.second)
    with pytest.raises(ValueError):U.check_dimensions(U.meter+1)
    with pytest.raises(ValueError):U.convert_to(U.meter,U.second,'unknown')
    x=S.Symbol('x');assert U.check_dimensions(x+U.meter)==x+U.meter


def test_aliases_modules_and_no_runtime_reference():
    assert set(U.__all__)==set(ref.__all__)
    assert all(hasattr(U,n) for n in U.__all__)
    assert U.m is U.meter and U.meters is U.meter
    from symox.physics.units.systems.si import SI,dimsys_SI
    assert SI is U.SI and dimsys_SI is SI.get_dimension_system()
    code="""
import sys
class Block:
 def find_spec(self,name,path=None,target=None):
  if name=='sympy' or name.startswith('sympy.'):raise AssertionError(name)
sys.meta_path.insert(0,Block())
from symox.physics.units import convert_to,kilometer,meter
assert convert_to(kilometer,meter)==1000*meter
"""
    p=subprocess.run([sys.executable,'-c',code],capture_output=True,text=True)
    assert p.returncode==0,p.stderr
