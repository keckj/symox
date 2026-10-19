"""Public Piecewise arguments are reconstructible expression/condition pairs."""
import pickle
import subprocess
import sys

import symox as sx
import sympy as sp


def test_branch_pairs_expose_values_conditions_and_reconstruction():
    x = sx.Symbol('x')
    expression = sx.Piecewise((x, x>0), (sx.Rational(1, 3), True))
    assert len(expression.args) == 2
    first, second = expression.args
    assert first.expr == x and first.cond == (x>0)
    assert first.args == (x, x>0)
    assert first.func(*first.args) == first
    assert second.expr == sx.Rational(1, 3) and second.cond == sx.true
    assert expression.func(*expression.args) == expression
    assert pickle.loads(pickle.dumps(expression)) == expression
    assert first.to_sympy() == expression.to_sympy().args[0]


def test_pair_reconstruction_preserves_replacement_and_folding():
    x, y = sx.symbols('x y')
    expression = sx.Piecewise((x, x>0), (0, True))
    expected = sx.Piecewise((y, y>0), (0, True))
    assert expression.replace(x, y) == expected
    assert expression.xreplace({x: y}) == expected
    folded = sx.piecewise_fold(expression+1)
    reference = sp.piecewise_fold(expression.to_sympy()+1)
    assert sp.simplify(folded.to_sympy()-reference) == 0


def test_piecewise_constructor_accepts_native_and_python_sequences():
    x = sx.Symbol('x')
    expected = sx.Piecewise((x, x>0), (0, True))
    assert sx.Piecewise([x, x>0], [0, True]) == expected
    assert sx.Piecewise(sx.Tuple(x, x>0), sx.Tuple(0, True)) == expected


def test_pair_reflection_and_replacement_do_not_import_reference():
    script = '''
import sys
class NoReference:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'sympy' or fullname.startswith('sympy.'):
            raise ImportError('reference forbidden')
sys.meta_path.insert(0, NoReference())
import symox as sx
x,y=sx.symbols('x y')
p=sx.Piecewise((x,x>0),(0,True))
assert p.args[0].expr==x
assert p.func(*p.args)==p
assert p.replace(x,y)==sx.Piecewise((y,y>0),(0,True))
'''
    subprocess.run([sys.executable, '-c', script], check=True, timeout=20)
