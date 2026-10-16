"""Unquoted UTF-8 names must never be sliced at a byte interior."""
import pytest
import sympy as sp
import symox as sx


@pytest.mark.parametrize('text', ['α', 'α+β', '变量+1', 'x+é', 'α_2**2+β3', 'sin(θ)+cos(φ)', 'α<=β', '2*变量/(é+1)'])
def test_unicode_names_match_sympy(text):
    assert sx.sympify(text).to_sympy() == sp.sympify(text)


def test_unicode_whitespace_and_invalid_characters():
    assert sx.sympify('α\u2003+\u00a0β') == sx.Symbol('α')+sx.Symbol('β')
    for text in ['x+🙂', '🙂', 'x+→', 'α@β']:
        with pytest.raises(ValueError):
            sx.sympify(text)
