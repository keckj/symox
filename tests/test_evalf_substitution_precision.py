"""Numerical substitutions must not round arithmetic before evalf's refinement."""
import subprocess
import sys

import mpmath
import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('digits', [15, 50, 120])
@pytest.mark.parametrize('scale', ['1e16', '1e100', '1e1000'])
@pytest.mark.parametrize('operation', ['linear', 'square', 'sin'])
def test_substituted_cancellation_preserves_the_requested_precision(digits, scale, operation):
    def build(module):
        x,y,z=module.symbols('x y z')
        expression=x+y-z
        if operation=='square':expression=expression**2
        elif operation=='sin':expression=module.sin(expression)
        large=module.Float(scale,20)
        return expression,{x:large,y:1,z:large}
    expression,values=build(sx)
    reference,foreign_values=build(sp)
    # T + 1 - T = 1 for the identical finite binary literal T. The reference
    # needs an increased refinement budget for the 1000-digit cancellation.
    expected=(sp.sin(1) if operation=='sin' else sp.Integer(1)).evalf(digits)
    assert reference.evalf(digits,subs=foreign_values,maxn=1500)==expected
    actual=expression.evalf(digits,subs=values)
    assert actual.to_sympy()==expected
    assert expression.n(digits,subs=values)==actual
    assert sx.N(expression,digits,subs=values)==actual


def test_machine_float_substitutions_and_matrix_entries_are_not_eagerly_rounded():
    x,y,z=sx.symbols('x y z')
    expression=x+y-z;values={x:1e16,y:1,z:1e16}
    assert expression.evalf(subs=values)==sx.Float(1)
    assert sx.Matrix([expression]).evalf(60,subs=values)[0]==sx.Float(1,60)
    assert sx.N(sx.Matrix([expression]),60,subs=values)[0]==sx.Float(1,60)


def test_substitution_into_previous_values_remains_sequential():
    for module in [sx,sp]:
        x,y=module.symbols('x y')
        assert (x+y).evalf(subs={x:y,y:2})==module.Float(4)
    x,y=sx.symbols('x y')
    assert (x+y).evalf(subs={x:y+1,y:2})==sx.Float(5)


@pytest.mark.parametrize('method', ['evalf', 'n', 'N'])
def test_positional_substitution_mapping_has_the_same_precision(method):
    x,y,z=sx.symbols('x y z');expression=x+y-z
    values={x:1e16,y:1,z:1e16}
    actual=(sx.N(expression,80,subs=values) if method=='N'
            else getattr(expression,method)(80,values))
    assert actual==sx.Float(1,80)


def test_existing_float_terms_also_wait_for_numerical_refinement():
    x,y=sx.symbols('x y');large=sx.Float('1e100',20)
    expression=sx.Add(x,y,-large,evaluate=False)
    assert expression.evalf(80,subs={x:large,y:1})==sx.Float(1,80)


def test_both_complex_components_refine_after_substitution():
    def build(module):
        x,y,z=module.symbols('x y z')
        large=module.Float('1e100',20)+module.I*module.Float('1e200',20)
        value=1+module.I*module.Rational(1,10**50)
        return x+y-z,{x:large,y:value,z:large}
    expression,values=build(sx);reference,foreign_values=build(sp)
    expected=1+sp.I*sp.Rational(1,10**50)
    for actual in [reference.evalf(80,subs=foreign_values,maxn=800),
                   expression.evalf(80,subs=values).to_sympy()]:
        for a,b in zip(actual.as_real_imag(),expected.as_real_imag()):
            assert abs((a-b)/b)<sp.Rational(1,10**75)


def test_substituted_cancellation_keeps_a_tiny_component_at_a_branch_cut():
    def build(module):
        x,y,z=module.symbols('x y z')
        return module.sqrt(-1+module.I*(x+y-z)),{x:1e16,y:module.Rational(1,10**100),z:1e16}
    expression,values=build(sx)
    # The reference's ordinary evalf can return a low-accuracy scaled zero for
    # the real component here. An independent context resolves it first.
    context=mpmath.mp.clone();context.dps=400
    root=context.sqrt(context.mpc(-1,context.mpf('1e-100')))
    expected=(sp.Float(str(root.real),360),sp.Float(str(root.imag),360))
    actual=expression.evalf(80,subs=values).to_sympy()
    # Check components separately: a unit-sized imaginary part cannot mask
    # destruction of the positive real component of order 10**-100.
    for a,b in zip(actual.as_real_imag(),expected):
        assert b!=0
        assert abs((a-b)/b)<sp.Rational(1,10**75)


def test_unresolved_functions_and_algebraic_keys_keep_existing_substitution_semantics():
    x,y=sx.symbols('x y');f=sx.Function('f')
    assert f(x).evalf(subs={x:y})==f(y)
    expression=x**4-2
    values={x**2:3}
    assert expression.evalf(subs=values)==expression.subs(values).evalf()


def test_precision_substitutions_need_no_reference_imports():
    script='''
import sys
class NoReference:
    def find_spec(self,fullname,path=None,target=None):
        if fullname=='sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0,NoReference())
import symox as s
x,y,z=s.symbols('x y z')
large=s.Float('1e1000',20)
values={x:large,y:1,z:large}
assert (x+y-z).evalf(100,subs=values)==s.Float(1,100)
assert s.sin(x+y-z).evalf(100,subs=values)==s.sin(1).evalf(100)
assert s.Matrix([x+y-z]).evalf(100,subs=values)[0]==s.Float(1,100)
'''
    subprocess.run([sys.executable,'-c',script],check=True,timeout=20)
