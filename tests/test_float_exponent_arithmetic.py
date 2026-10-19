"""Exponent-counter limits preserve exact expressions instead of wrapping values."""
import random
import subprocess
import sys

import symox as sx
import sympy as sp


def test_huge_products_preserve_exact_held_arithmetic_at_counter_boundary():
    literal='1e999999999999999999'
    value=sx.Float(literal,20);reference=sp.Float(literal,20)
    square=value*value;expected_square=reference*reference
    assert square.is_Float and square.to_sympy()==expected_square
    assert (value**2).to_sympy()==reference**2
    cube=square*value
    assert cube.is_Mul and not cube.is_Float
    assert cube.to_sympy().evalf(20)==expected_square*reference
    evaluated=sx.N(cube,30)
    assert evaluated.is_Mul
    held=sp.Mul(expected_square,reference,evaluate=False)
    assert evaluated.to_sympy().evalf(30)==held.evalf(30)
    matrix=sx.Matrix([[value]])*sx.Matrix([[value]])*sx.Matrix([[value]])
    assert matrix[0,0].is_Mul
    assert matrix[0,0].to_sympy().evalf(20)==expected_square*reference


def test_exponent_gaps_do_not_overflow_in_addition_or_division():
    value=sx.Float('1e999999999999999999',20)
    reference=sp.Float('1e999999999999999999',20)
    square=value*value;R=reference*reference
    assert (square+1/square).to_sympy()==R+1/R
    assert (square-square).is_zero
    assert (value/value).to_sympy()==reference/reference
    outside=square/sx.Float('1e-999999999999999999',20)
    assert not outside.is_Float
    expected=sp.Mul(R,1/sp.Float('1e-999999999999999999',20),evaluate=False)
    assert outside.to_sympy().evalf(20)==expected.evalf(20)
    complex_literal=value+sx.Float('.125',20)*sx.I
    assert sx.N(complex_literal,30).to_sympy()==sp.N(reference+sp.Float('.125',20)*sp.I,30)


def test_float_division_rounds_the_quotient_once():
    rng=random.Random(1811)
    for precision in [3,7,24,53,107,257]:
        for _ in range(80):
            numerator=str(rng.randrange(-10**10,10**10))
            denominator=str(rng.randrange(1,10**10))
            a=sx.Float(numerator,precision=precision);b=sx.Float(denominator,precision=precision)
            A=sp.Float(numerator,precision=precision);B=sp.Float(denominator,precision=precision)
            assert (a/b).to_sympy()==A/B


def test_inverse_of_large_integer_retains_sign_and_exactness():
    for sign in [-1,1]:
        denominator=sign*10**20000
        value=sx.Rational(1,denominator)
        assert value.to_sympy()==sp.Rational(1,denominator)
        assert (value*denominator)==1


def test_counter_boundary_fallback_never_imports_reference():
    source='''
import sys
class Block:
    def find_spec(self,fullname,*args):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise AssertionError('reference import')
sys.meta_path.insert(0,Block())
from symox import *
value=Float('1e999999999999999999',20)
square=value*value
assert square.is_Float
assert (square*value).is_Mul
assert N(square*value,30).is_Mul
assert (square+1/square)==square
assert (Matrix([[square]])*Matrix([[value]]))[0,0].is_Mul
'''
    subprocess.run([sys.executable,'-c',source],check=True,timeout=15)
