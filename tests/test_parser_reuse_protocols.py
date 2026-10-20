"""Lexical reuse must not reuse mathematical values or namespace effects."""
import pytest
import symox as sx
from symox.parsing import sympy_parser as parser


def test_repeated_source_uses_current_local_functions_and_symbols():
    source = 'f(x) + 1'
    x, y = sx.symbols('x y')
    calls = []
    local = {'x': x, 'f': lambda value: calls.append(value) or value**2}
    assert sx.sympify(source, locals=local) == x**2+1
    local['x'] = y
    local['f'] = lambda value: calls.append(value) or value**3
    assert sx.sympify(source, locals=local) == y**3+1
    assert calls == [x, y]


def test_user_transformations_run_again_with_fresh_tokens():
    calls = []
    def mutate(tokens, local_dict, global_dict):
        calls.append(True)
        tokens[0] = (parser.NUMBER, str(len(calls)))
        return tokens
    transformations = (mutate,) + parser.standard_transformations
    assert parser.parse_expr('0+x', transformations=transformations) == 1+sx.Symbol('x')
    assert parser.parse_expr('0+x', transformations=transformations) == 2+sx.Symbol('x')
    assert parser.parse_expr('0+x') == sx.Symbol('x')


@pytest.mark.parametrize('space', ['\v', '\r', '\x1c', '\u00a0', '\u2003'])
def test_nonstandard_whitespace_still_uses_token_aware_normalization(space):
    x, y = sx.symbols('x y')
    assert sx.sympify('x'+space+'+y') == x+y


@pytest.mark.parametrize('space', [' ', '\t', '\f', '\u00a0', '\u2003'])
def test_quoted_whitespace_is_never_replaced(space):
    value = 'a'+space+'b'
    assert sx.sympify('"'+value+'"') == value


def test_lexical_reuse_preserves_current_numerical_and_held_options():
    source = 'x + 0.125 + x'
    sx.sympify(source)
    assert sx.sympify(source, rational=True) == 2*sx.Symbol('x')+sx.Rational(1,8)
    with sx.evaluate(False):
        held = sx.sympify(source)
    assert held.args == (sx.Symbol('x'), sx.Float('0.125'), sx.Symbol('x'))
    assert sx.sympify(source) == 2*sx.Symbol('x')+sx.Float('0.125')
