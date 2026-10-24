"""Complete-root guards for exact distribution expansion hints."""
import subprocess
import sys

import pytest
import symox as sx
import sympy as sp


@pytest.mark.parametrize('formula', ['x', '2*x-2', '-3*x+2', 'x*x-1', 'x**3-x',
    '(x-1)*(x+2)*(x-3)', 'x*x', '(x-1)**2*(x+2)', 'sin(x)'])
@pytest.mark.parametrize('method', [False, True])
def test_dirac_expansion_matches_complete_reference_cases(formula, method):
    def run(module):
        x = module.Symbol('x', real=True)
        expression = module.DiracDelta(module.sympify(formula, locals={'x':x}))
        return expression.expand(diracdelta=True, wrt=x) if method else module.expand(expression, diracdelta=True, wrt=x)
    assert sp.simplify(run(sx).to_sympy()-run(sp)) == 0


def test_affine_parameters_have_proved_real_nonzero_scale():
    x = sx.Symbol('x', real=True)
    a = sx.Symbol('a', positive=True)
    b = sx.Symbol('b', real=True)
    assert sx.expand(sx.DiracDelta(a*(x-b)), diracdelta=True, wrt=x) == sx.DiracDelta(x-b)/a
    unknown = sx.Symbol('unknown', real=True)
    expression = sx.DiracDelta(unknown*x)
    # The singular zero-scale stratum is not silently divided away.
    assert sx.expand(expression, diracdelta=True, wrt=x) == expression


def test_variable_inference_errors_and_deep_hint():
    x, y = sx.symbols('x y', real=True)
    with pytest.raises(TypeError):
        sx.expand(sx.DiracDelta(x*y), diracdelta=True)
    with pytest.raises(TypeError):
        sx.expand(sx.DiracDelta(0), diracdelta=True)
    assert sx.expand(sx.DiracDelta(2*x), diracdelta=True) == sx.DiracDelta(x)/2
    expression = 1+sx.DiracDelta(2*x)
    assert sx.expand(expression, diracdelta=True, deep=False) == expression
    assert sx.expand(expression, diracdelta=True) == 1+sx.DiracDelta(x)/2
    assert sx.expand(expression, diracdelta=False) == expression
    undefined = sx.Function('DiracDelta')(2*x)
    assert sx.expand(undefined, diracdelta=True) == undefined


@pytest.mark.parametrize('order', [1, 2, 3])
def test_derivative_orders_remain_held_and_negative_argument_parity_is_preserved(order):
    x = sx.Symbol('x', real=True)
    expression = sx.DiracDelta(2*x, order)
    assert sx.expand(expression, diracdelta=True, wrt=x) == expression
    assert sx.DiracDelta(-x, order) == (-1)**order*sx.DiracDelta(x, order)


def test_explicit_matrix_entries_expand_and_keep_their_storage_contract():
    x = sx.Symbol('x', real=True)
    for constructor in [sx.Matrix, sx.ImmutableMatrix, sx.SparseMatrix, sx.ImmutableSparseMatrix]:
        matrix = constructor([[sx.DiracDelta(2*x), 0]])
        result = sx.expand(matrix, diracdelta=True, wrt=x)
        assert result[0, 0] == sx.DiracDelta(x)/2
        # SymPy 1.14: expand sympifies the matrix to its immutable class.
        assert result._immutable
        assert bool(getattr(result, '_symox_sparse', False)) == bool(getattr(matrix, '_symox_sparse', False))


def test_quintic_real_support_is_not_erased_by_incomplete_radical_roots():
    x = sx.Symbol('x', real=True)
    expression = sx.DiracDelta(x**5-x-1)
    expanded = sx.expand(expression, diracdelta=True, wrt=x)
    nodes = expanded.find(lambda node: getattr(node, 'class_name', None) == 'DiracDelta')
    assert len(nodes) == 1
    root = x-next(iter(nodes)).args[0]
    assert sx._core._exact_assumption(root, 'real') is True
    reference_x = sp.Symbol('x', real=True)
    polynomial = reference_x**5-reference_x-1
    assert sp.Poly(polynomial, reference_x).count_roots(-sp.oo, sp.oo) == 1
    assert sp.minimal_polynomial(root.to_sympy(), reference_x) == polynomial
    # SymPy 1.14's expansion uses its incomplete radical roots() dictionary
    # and returns zero here. This is an intentional mathematical difference.
    assert sp.expand(sp.DiracDelta(polynomial), diracdelta=True, wrt=reference_x) == 0


def test_distribution_expansion_does_not_import_reference_backend():
    script = '''
import sys
class BlockReference:
 def find_spec(self,fullname,path=None,target=None):
  if fullname=='sympy' or fullname.startswith('sympy.'):
   raise ImportError('reference forbidden')
sys.meta_path.insert(0,BlockReference())
import symox as s
x=s.Symbol('x',real=True)
assert s.expand(s.DiracDelta(x*x-1),diracdelta=True,wrt=x)==(s.DiracDelta(x-1)+s.DiracDelta(x+1))/2
assert 'sympy' not in sys.modules
'''
    subprocess.run([sys.executable, '-c', script], check=True, timeout=10)


def test_unqualified_root_domains_and_resource_limits_retain_their_contract():
    x = sx.Symbol('x', real=True)
    for argument in [sx.sin(x), 1/(x*x-1), sx.Float('1e300')*x-sx.Float('1e-300')]:
        expression = sx.DiracDelta(argument)
        assert sx.expand(expression, diracdelta=True, wrt=x) == expression
    expression = sx.DiracDelta(x**3-x)
    with pytest.raises(ValueError, match='max_terms'):
        sx.expand(expression, diracdelta=True, wrt=x, max_terms=2)
    # The complete square-free polynomial has no real support. This zero
    # follows from exact classification of all roots, not missing radicals.
    assert sx.expand(sx.DiracDelta(x*x+1), diracdelta=True, wrt=x) == 0
